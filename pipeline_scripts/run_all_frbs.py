"""Batch driver for master_run.py.

Iterates every FRB that has both ``<frb>_flux.fits`` and ``<frb>_invvar.fits``
in ``large_cutouts/``, looks up the localisation from
``master_frb_localization.csv`` (project root), and forwards to
``master_run.py`` with ``--outputs all``.

By default only FRBs with ``coord_semantics == "host"`` are run — those are
the entries where ra_deg / dec_deg refer to the optical host position. Pass
``--include-signal`` to also process ``signal``-semantics rows (RA/Dec is
the radio FRB position; AstroPath behaviour is less informative there).

Phase 3a host centre (``master_run`` / ``galfit_config.yaml``):
    ``--use-localization-host`` — CSV RA/Dec + nearest galaxy (SPREAD cut).
    ``--use-astropath-host`` — AstroPath ``sex_number`` from posteriors (default
    when ``cutouts.use_localization_host`` is false).

Localisation parameters passed to master_run:
    --ra            <- ra_deg
    --dec           <- dec_deg
    --err-a-arcsec  <- major_sigma_as if available, else ra_err_as
    --err-b-arcsec  <- minor_sigma_as if available, else dec_err_as
    --err-theta-deg <- pa_deg

If any of the error / PA columns is empty for an FRB, the corresponding
flag is **omitted** so master_run / the YAML defaults take over (1.0", 1.0",
0 deg). Coordinates are mandatory; FRBs without RA/Dec are skipped.

Usage (PowerShell):
    python pipeline_scripts/run_all_frbs.py                 # host-only
    python pipeline_scripts/run_all_frbs.py --include-signal
    python pipeline_scripts/run_all_frbs.py --skip-existing
    python pipeline_scripts/run_all_frbs.py --frb 20190608B 20180924B
    python pipeline_scripts/run_all_frbs.py --dry-run
    python pipeline_scripts/run_all_frbs.py --use-localization-host --frb 20190608B
    python pipeline_scripts/run_all_frbs.py --use-astropath-host --frb 20221101B
"""
import argparse
import base64
import concurrent.futures
import json
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import pandas as pd

try:
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    _HAS_ASTROPY = True
except ImportError:
    _HAS_ASTROPY = False


PIPELINE_DIR = Path(__file__).resolve().parent
REPO_ROOT = PIPELINE_DIR.parent
MASTER = PIPELINE_DIR / "master_run.py"
CUTOUT_DIR = REPO_ROOT / "large_cutouts"
LOC_CSV = REPO_ROOT / "master_frb_localization.csv"
OUTPUT_DIR = PIPELINE_DIR / "Output"


def parse_float(value):
    """float(value) or None if blank / NaN / unparseable."""
    if value is None:
        return None
    try:
        s = str(value).strip()
    except Exception:
        return None
    if s == "" or s.lower() == "nan":
        return None
    try:
        v = float(s)
        return None if v != v else v  # filter NaN
    except (ValueError, TypeError):
        return None


def get_localization(loc_df: pd.DataFrame, frb_name: str):
    """Return dict(ra, dec, err_a, err_b, err_theta, semantics) or None if RA/Dec missing.

    err_* are None when the CSV cell is blank — caller must omit the
    corresponding --err-* flag so the YAML default applies.
    """
    row = loc_df.loc[loc_df["frb"] == frb_name]
    if len(row) == 0:
        return None
    r = row.iloc[0]

    ra = parse_float(r.get("ra_deg"))
    dec = parse_float(r.get("dec_deg"))
    if ra is None or dec is None:
        return None

    # Prefer the explicit ellipse if present; fall back to RA/Dec sigmas.
    major = parse_float(r.get("major_sigma_as"))
    minor = parse_float(r.get("minor_sigma_as"))
    if major is not None and minor is not None:
        err_a, err_b = major, minor
    else:
        err_a = parse_float(r.get("ra_err_as"))
        err_b = parse_float(r.get("dec_err_as"))

    err_theta = parse_float(r.get("pa_deg"))
    semantics = (str(r.get("coord_semantics", "")).strip().lower() or "unknown")
    return {"ra": ra, "dec": dec, "err_a": err_a, "err_b": err_b,
            "err_theta": err_theta, "semantics": semantics}


def astropath_vs_csv_sep(posteriors_path: Path, csv_ra: float, csv_dec: float):
    """Separation (arcsec) between best AstroPath host and CSV host coords."""
    if not _HAS_ASTROPY or not posteriors_path.is_file():
        return None
    try:
        df = pd.read_csv(posteriors_path)
        if not len(df) or "posterior_O" not in df.columns:
            return None
        need = {"ra_deg", "dec_deg"}
        if not need.issubset(df.columns):
            return None
        best = df.sort_values("posterior_O", ascending=False).iloc[0]
        c_csv = SkyCoord(ra=csv_ra, dec=csv_dec, unit="deg", frame="icrs")
        c_ap = SkyCoord(ra=float(best["ra_deg"]), dec=float(best["dec_deg"]), unit="deg", frame="icrs")
        return float(c_csv.separation(c_ap).arcsec)
    except Exception:
        return None


def collect_summary(out_dir: Path, csv_ra=None, csv_dec=None):
    """Pull a few key science values from a completed run for the table."""
    s = {
        "posterior": None,
        "zp_40px": None,
        "n_stars": None,
        "ref_cat": None,
        "sep_astropath_vs_csv_arcsec": None,
    }
    posteriors = out_dir / "astropath_posteriors.csv"
    if posteriors.exists():
        try:
            df = pd.read_csv(posteriors)
            if len(df) and "posterior_O" in df.columns:
                s["posterior"] = float(df["posterior_O"].max())
            if csv_ra is not None and csv_dec is not None:
                s["sep_astropath_vs_csv_arcsec"] = astropath_vs_csv_sep(
                    posteriors, csv_ra, csv_dec
                )
        except Exception:
            pass
    zp = out_dir / "zero_points.json"
    if zp.exists():
        try:
            j = json.loads(zp.read_text())
            s["zp_40px"] = j.get("zp_aper") or j.get("zp_aper_40px")
            s["n_stars"] = j.get("n_calibration_stars")
            s["ref_cat"] = j.get("reference_catalog")
        except Exception:
            pass
    return s


def fmt(value, spec="", default="—"):
    if value is None:
        return default
    try:
        return format(value, spec) if spec else str(value)
    except Exception:
        return default


def format_eta(seconds):
    """Human-readable duration string."""
    s = int(round(seconds))
    if s < 60:
        return f"{s}s"
    m, s = divmod(s, 60)
    if m < 60:
        return f"{m}m {s}s" if s else f"{m}m"
    h, m = divmod(m, 60)
    return f"{h}h {m}m" if m else f"{h}h"


def _run_single_frb(cmd, out_subdir_str, log_path_str, loc_ra, loc_dec):
    """Execute one FRB pipeline run (top-level for ProcessPoolExecutor)."""
    out_subdir = Path(out_subdir_str)
    log_path = Path(log_path_str)
    out_subdir.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    with open(log_path, "w", encoding="utf-8") as lf:
        proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT,
                              stdin=subprocess.DEVNULL)
    elapsed = time.time() - t_start
    summary = collect_summary(out_subdir, loc_ra, loc_dec)
    return proc.returncode, elapsed, summary


def _img_to_base64(path: Path) -> str:
    """Read an image file and return a base64-encoded data URI."""
    data = path.read_bytes()
    return f"data:image/png;base64,{base64.b64encode(data).decode()}"


def generate_html_report(results, output_dir: Path, tag: str):
    """Write pipeline_scripts/Output/batch_report.html."""
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    png_names = ["galfit_results.png", "astropath_association.png", "qa_cutout_mask.png"]

    rows_html = []
    for frb, status, note, summary in results:
        rows_html.append(
            f"<tr><td>{frb}</td><td>{status}</td>"
            f"<td>{fmt(summary.get('posterior'), '.3f')}</td>"
            f"<td>{fmt(summary.get('sep_astropath_vs_csv_arcsec'), '.2f')}</td>"
            f"<td>{fmt(summary.get('zp_40px'), '.3f')}</td>"
            f"<td>{fmt(summary.get('n_stars'))}</td>"
            f"<td>{summary.get('ref_cat') or ''}</td>"
            f"<td>{note}</td></tr>"
        )

    frb_sections = []
    for frb, status, _note, _summary in results:
        out_dir = output_dir / f"{frb}_{tag}"
        if not out_dir.is_dir():
            continue
        images_html = []
        for png in png_names:
            p = out_dir / png
            if p.is_file():
                images_html.append(
                    f'<div class="thumb"><img src="{_img_to_base64(p)}" '
                    f'alt="{png}"><div class="caption">{png}</div></div>'
                )
        summary_json = out_dir / "pipeline_summary.json"
        link = ""
        if summary_json.is_file():
            rel = summary_json.resolve()
            link = f'<a href="file:///{rel}">pipeline_summary.json</a>'
        if images_html or link:
            frb_sections.append(
                f'<div class="frb-block"><h3>{frb} [{status}]</h3>'
                f'<div class="grid">{"".join(images_html)}</div>'
                f'{link}</div>'
            )

    html = f"""\
<!DOCTYPE html>
<html><head><meta charset="utf-8">
<title>Pipeline Batch Report &ndash; {now}</title>
<style>
  body {{ font-family: Arial, sans-serif; margin: 20px; background: #fafafa; }}
  h1 {{ color: #333; }}
  table {{ border-collapse: collapse; width: 100%; margin: 1em 0; }}
  th, td {{ border: 1px solid #ccc; padding: 6px 10px; text-align: left; }}
  th {{ background: #4a76a8; color: #fff; }}
  tr:nth-child(even) {{ background: #eef; }}
  .grid {{ display: flex; flex-wrap: wrap; gap: 12px; margin: 8px 0; }}
  .thumb img {{ max-width: 340px; border: 1px solid #aaa; }}
  .caption {{ font-size: 0.85em; color: #555; text-align: center; }}
  .frb-block {{ margin: 1.5em 0; padding: 10px; background: #fff;
               border: 1px solid #ddd; border-radius: 4px; }}
</style></head><body>
<h1>Pipeline Batch Report</h1>
<p>Generated: {now}</p>
<h2>Summary</h2>
<table>
<tr><th>FRB</th><th>Status</th><th>P(O)</th><th>d_host</th>
    <th>ZP_40</th><th>N*</th><th>Ref catalog</th><th>Note</th></tr>
{"".join(rows_html)}
</table>
<h2>FRB Details</h2>
{"".join(frb_sections) if frb_sections else "<p>No FRB output directories found.</p>"}
</body></html>"""

    report_path = output_dir / "batch_report.html"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(html, encoding="utf-8")
    return report_path


def main():
    parser = argparse.ArgumentParser(
        description="Run master_run.py against every FRB with cutouts in large_cutouts/.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--frb", nargs="*", default=None,
                        help="Restrict to specific FRB names (default: all FRBs with cutouts).")
    parser.add_argument("--include-signal", action="store_true",
                        help="Also process FRBs with coord_semantics='signal'. "
                             "Default: only coord_semantics='host' rows are run "
                             "(those are the ones where RA/Dec refers to the optical host).")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip FRBs whose Output/<frbname>_<tag>/ already exists.")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print the master_run command for each FRB without executing.")
    parser.add_argument("--outputs", nargs="+", default=["all"],
                        help="Forwarded verbatim to master_run --outputs (default: all).")
    parser.add_argument("--keep-workdir", action="store_true",
                        help="Forward --keep-workdir to master_run.")
    parser.add_argument("--condensed", action="store_true",
                        help="Forward --condensed to master_run (summary JSON + 3 PNGs only).")
    parser.add_argument(
        "--use-localization-host",
        action="store_true",
        help="Forward to master_run: Phase 3a uses CSV host RA/Dec, not AstroPath pick.",
    )
    parser.add_argument(
        "--use-astropath-host",
        action="store_true",
        help="Forward to master_run: Phase 3a uses AstroPath posteriors (overrides "
             "galfit_config.yaml cutouts.use_localization_host).",
    )
    parser.add_argument(
        "--list-file",
        type=Path,
        default=None,
        help="FRB list: .txt (one name per line) or .csv with a 'frb' column (e.g. new_hosts_master.csv).",
    )
    parser.add_argument(
        "--parallel", type=int, default=1, metavar="N",
        help="Number of parallel workers (default: 1 = sequential).",
    )
    parser.add_argument(
        "--no-auto-refresh", action="store_true",
        help="Skip automatic compare_pipeline_galfit_vs_master.py run after the batch.",
    )
    parser.add_argument(
        "--rerun-phase",
        choices=["1", "2", "3a", "statmorph", "3b"],
        default=None,
        help="Forwarded to master_run: skip earlier phases and re-execute from this phase.",
    )
    args = parser.parse_args()

    if not LOC_CSV.exists():
        raise SystemExit(f"[batch] localization CSV not found at {LOC_CSV}")
    if not CUTOUT_DIR.is_dir():
        raise SystemExit(f"[batch] large_cutouts/ not found at {CUTOUT_DIR}")
    loc_df = pd.read_csv(LOC_CSV, dtype=str).fillna("")

    # Validate --outputs upfront so a typo fails once here instead of once per FRB.
    _VALID_OUTPUTS = {"catalog", "psf", "photometry", "astropath", "statmorph", "galfit", "all"}
    bad_outputs = [o for o in args.outputs if o.lower() not in _VALID_OUTPUTS]
    if bad_outputs:
        raise SystemExit(
            f"[batch] Unknown --outputs entries {bad_outputs}. "
            f"Valid: {sorted(_VALID_OUTPUTS)}."
        )

    frb_names = list(args.frb or [])
    if args.list_file:
        if not args.list_file.is_file():
            raise SystemExit(f"[batch] --list-file not found: {args.list_file}")
        if args.list_file.suffix.lower() == ".csv":
            list_df = pd.read_csv(args.list_file)
            if "frb" not in list_df.columns:
                raise SystemExit(
                    f"[batch] --list-file {args.list_file} has no 'frb' column "
                    f"(columns: {list(list_df.columns)})."
                )
            frb_names.extend(list_df["frb"].astype(str).tolist())
        else:
            frb_names.extend(
                ln.strip()
                for ln in args.list_file.read_text().splitlines()
                if ln.strip() and not ln.startswith("#")
            )

    flux_files = sorted(CUTOUT_DIR.glob("*_flux.fits"))
    if frb_names:
        wanted = set(frb_names)
        flux_files = [f for f in flux_files if f.stem.removesuffix("_flux") in wanted]

    # Tag string mirrors master_run's resolve_outputs() (lowercase, deduped,
    # sorted) so we can pre-compute the per-FRB output folder name for
    # --skip-existing.
    outputs_norm = [o.lower() for o in args.outputs]
    is_all = ("all" in outputs_norm) or len(outputs_norm) == 0
    tag = "all" if is_all else "_".join(sorted(dict.fromkeys(outputs_norm)))

    print(f"[batch] {len(flux_files)} FRB(s) to process; outputs={args.outputs}")
    print(f"[batch] localization table: {LOC_CSV}")
    print(f"[batch] outputs root      : {OUTPUT_DIR}")
    print()

    results = []
    t0 = time.time()

    # Phase 1: filter FRBs and build commands
    tasks = []
    for i, flux in enumerate(flux_files, 1):
        frb = flux.stem.removesuffix("_flux")
        invvar = CUTOUT_DIR / f"{frb}_invvar.fits"
        prefix = f"[{i:2d}/{len(flux_files)}] {frb}"

        if not invvar.exists():
            print(f"{prefix} ... SKIP (no invvar)")
            results.append((frb, "SKIP", "no invvar", {}))
            continue

        loc = get_localization(loc_df, frb)
        if loc is None:
            print(f"{prefix} ... SKIP (no RA/Dec in CSV)")
            results.append((frb, "SKIP", "no localization", {}))
            continue
        if loc["semantics"] != "host" and not args.include_signal:
            print(f"{prefix} ... SKIP (coord_semantics={loc['semantics']!r}; pass --include-signal to run)")
            results.append((frb, "SKIP", f"semantics={loc['semantics']}", {}))
            continue

        out_subdir = OUTPUT_DIR / f"{frb}_{tag}"
        if args.skip_existing and out_subdir.exists():
            print(f"{prefix} ... SKIP (Output already exists)")
            results.append((frb, "SKIP", "exists", collect_summary(out_subdir, loc["ra"], loc["dec"])))
            continue

        cmd = [sys.executable, str(MASTER),
               "--image", str(flux), "--invvar", str(invvar),
               "--ra", str(loc["ra"]), "--dec", str(loc["dec"]),
               "--outputs", *args.outputs]
        if args.use_localization_host:
            cmd.append("--use-localization-host")
        if args.use_astropath_host:
            cmd.append("--use-astropath-host")
        if loc["err_a"] is not None:
            cmd += ["--err-a-arcsec", str(loc["err_a"])]
        if loc["err_b"] is not None:
            cmd += ["--err-b-arcsec", str(loc["err_b"])]
        if loc["err_theta"] is not None:
            cmd += ["--err-theta-deg", str(loc["err_theta"])]
        if args.keep_workdir:
            cmd.append("--keep-workdir")
        if args.condensed:
            cmd.append("--condensed")
        if args.rerun_phase:
            cmd += ["--rerun-phase", args.rerun_phase]

        if args.dry_run:
            print(f"{prefix} DRY: {' '.join(cmd)}")
            continue

        log_path = out_subdir / "master_run.log"
        tasks.append((frb, cmd, out_subdir, log_path, loc, prefix))

    # Phase 2: execute tasks (parallel or sequential)
    if not args.dry_run and tasks:
        n_tasks = len(tasks)
        if args.parallel > 1:
            print(f"\n[batch] Launching {n_tasks} FRB(s) with {args.parallel} workers...")
            with concurrent.futures.ProcessPoolExecutor(max_workers=args.parallel) as pool:
                future_map = {}
                for frb, cmd, out_subdir, log_path, loc, prefix in tasks:
                    fut = pool.submit(
                        _run_single_frb, cmd,
                        str(out_subdir), str(log_path),
                        loc["ra"], loc["dec"],
                    )
                    future_map[fut] = (frb, out_subdir, loc, prefix)
                for fut in concurrent.futures.as_completed(future_map):
                    frb, out_subdir, loc, prefix = future_map[fut]
                    rc, elapsed, summary = fut.result()
                    if rc == 0:
                        p_str = fmt(summary["posterior"], ".3f")
                        print(f"{prefix} ... OK ({elapsed:.0f}s, P={p_str})")
                        results.append((frb, "OK", f"{elapsed:.0f}s", summary))
                    else:
                        log_rel = (out_subdir / "master_run.log").relative_to(REPO_ROOT)
                        print(f"{prefix} ... FAIL rc={rc}  (log: {log_rel})")
                        results.append((frb, "FAIL", f"rc={rc}", summary))
        else:
            ema_elapsed = None
            for idx, (frb, cmd, out_subdir, log_path, loc, prefix) in enumerate(tasks):
                ell = (f"err_a={fmt(loc['err_a'], '.2f')} "
                       f"err_b={fmt(loc['err_b'], '.2f')} "
                       f"pa={fmt(loc['err_theta'], '.1f')}")
                print(f"{prefix} (ra={loc['ra']:.4f}, dec={loc['dec']:.4f}; {ell}) ... ",
                      end="", flush=True)

                out_subdir.mkdir(parents=True, exist_ok=True)
                t_start = time.time()
                with open(log_path, "w", encoding="utf-8") as lf:
                    proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT,
                                          stdin=subprocess.DEVNULL)
                elapsed = time.time() - t_start
                summary = collect_summary(out_subdir, loc["ra"], loc["dec"])

                ema_elapsed = elapsed if ema_elapsed is None else 0.3 * elapsed + 0.7 * ema_elapsed
                remaining = n_tasks - (idx + 1)
                eta_part = f", ETA {format_eta(ema_elapsed * remaining)} remaining" if remaining > 0 else ""

                if proc.returncode == 0:
                    p_str = fmt(summary["posterior"], ".3f")
                    print(f"OK ({elapsed:.0f}s{eta_part}, P={p_str}, ref={summary['ref_cat'] or '—'})")
                    results.append((frb, "OK", f"{elapsed:.0f}s", summary))
                else:
                    print(f"FAIL rc={proc.returncode}{eta_part}  (log: {log_path.relative_to(REPO_ROOT)})")
                    results.append((frb, "FAIL", f"rc={proc.returncode}", summary))

    # ---- Auto-refresh: compare pipeline vs master ----
    if not args.no_auto_refresh and not args.dry_run:
        compare_script = REPO_ROOT / "scripts" / "compare_pipeline_galfit_vs_master.py"
        if compare_script.is_file():
            print("\n[batch] Auto-running compare_pipeline_galfit_vs_master.py ...")
            try:
                subprocess.run([sys.executable, str(compare_script)], check=False)
                print("[batch] Compare script finished.")
            except Exception as e:
                print(f"[batch] Compare script failed: {e}")

    # ---- Summary table ----
    total_elapsed = time.time() - t0
    print()
    print(f"[batch] finished in {total_elapsed:.0f}s "
          f"({sum(1 for r in results if r[1] == 'OK')} ok, "
          f"{sum(1 for r in results if r[1] == 'FAIL')} fail, "
          f"{sum(1 for r in results if r[1] == 'SKIP')} skip)")
    print()
    print(f"{'FRB':<14} {'status':<6} {'P(O)':>7} {'d_host':>8} {'ZP_40':>8} {'N*':>4} "
          f"{'ref_catalog':<20} note")
    print("-" * 88)
    for frb, status, note, summary in results:
        p = fmt(summary.get("posterior"), ".3f")
        sep = fmt(summary.get("sep_astropath_vs_csv_arcsec"), ".2f")
        zp = fmt(summary.get("zp_40px"), ".3f")
        n = fmt(summary.get("n_stars"))
        ref = summary.get("ref_cat") or ""
        print(f"{frb:<14} {status:<6} {p:>7} {sep:>8} {zp:>8} {n:>4} {ref:<20} {note}")

    # ---- HTML batch report ----
    try:
        report_path = generate_html_report(results, OUTPUT_DIR, tag)
        print(f"\n[batch] HTML report: {report_path.relative_to(REPO_ROOT)}")
    except Exception as e:
        print(f"\n[batch] HTML report generation failed: {e}")


if __name__ == "__main__":
    main()
