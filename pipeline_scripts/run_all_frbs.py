"""Batch driver for master_run.py.

Iterates every FRB that has both ``<frb>_flux.fits`` and ``<frb>_invvar.fits``
in ``large_cutouts/``, looks up the localisation from
``master_frb_localization.csv`` (project root), and forwards to
``master_run.py`` with ``--outputs all``.

By default only FRBs with ``coord_semantics == "host"`` are run — those are
the entries where ra_deg / dec_deg refer to the optical host position. Pass
``--include-signal`` to also process ``signal``-semantics rows (RA/Dec is
the radio FRB position; AstroPath behaviour is less informative there).

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
"""
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd


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


def collect_summary(out_dir: Path):
    """Pull a few key science values from a completed run for the table."""
    s = {"posterior": None, "zp_40px": None, "n_stars": None, "ref_cat": None}
    posteriors = out_dir / "astropath_posteriors.csv"
    if posteriors.exists():
        try:
            df = pd.read_csv(posteriors)
            if len(df) and "posterior_O" in df.columns:
                s["posterior"] = float(df["posterior_O"].max())
        except Exception:
            pass
    zp = out_dir / "zero_points.json"
    if zp.exists():
        try:
            j = json.loads(zp.read_text())
            s["zp_40px"] = j.get("zp_aper_40px")
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
    args = parser.parse_args()

    if not LOC_CSV.exists():
        raise SystemExit(f"[batch] localization CSV not found at {LOC_CSV}")
    if not CUTOUT_DIR.is_dir():
        raise SystemExit(f"[batch] large_cutouts/ not found at {CUTOUT_DIR}")
    loc_df = pd.read_csv(LOC_CSV, dtype=str).fillna("")

    flux_files = sorted(CUTOUT_DIR.glob("*_flux.fits"))
    if args.frb:
        wanted = set(args.frb)
        flux_files = [f for f in flux_files if f.stem.removesuffix("_flux") in wanted]

    # Tag string mirrors master_run's logic so we can pre-compute the per-FRB
    # output folder name for --skip-existing.
    is_all = ("all" in [o.lower() for o in args.outputs]) or len(args.outputs) == 0
    tag = "all" if is_all else "_".join(sorted(args.outputs))

    print(f"[batch] {len(flux_files)} FRB(s) to process; outputs={args.outputs}")
    print(f"[batch] localization table: {LOC_CSV}")
    print(f"[batch] outputs root      : {OUTPUT_DIR}")
    print()

    results = []
    t0 = time.time()
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
            results.append((frb, "SKIP", "exists", collect_summary(out_subdir)))
            continue

        cmd = [sys.executable, str(MASTER),
               "--image", str(flux), "--invvar", str(invvar),
               "--ra", str(loc["ra"]), "--dec", str(loc["dec"]),
               "--outputs", *args.outputs]
        if loc["err_a"] is not None:
            cmd += ["--err-a-arcsec", str(loc["err_a"])]
        if loc["err_b"] is not None:
            cmd += ["--err-b-arcsec", str(loc["err_b"])]
        if loc["err_theta"] is not None:
            cmd += ["--err-theta-deg", str(loc["err_theta"])]
        if args.keep_workdir:
            cmd.append("--keep-workdir")

        if args.dry_run:
            print(f"{prefix} DRY: {' '.join(cmd)}")
            continue

        ell = (f"err_a={fmt(loc['err_a'], '.2f')} "
               f"err_b={fmt(loc['err_b'], '.2f')} "
               f"pa={fmt(loc['err_theta'], '.1f')}")
        print(f"{prefix} (ra={loc['ra']:.4f}, dec={loc['dec']:.4f}; {ell}) ... ",
              end="", flush=True)

        out_subdir.mkdir(parents=True, exist_ok=True)
        log_path = out_subdir / "master_run.log"
        t_start = time.time()
        with open(log_path, "w", encoding="utf-8") as lf:
            proc = subprocess.run(cmd, stdout=lf, stderr=subprocess.STDOUT)
        elapsed = time.time() - t_start
        summary = collect_summary(out_subdir)

        if proc.returncode == 0:
            p_str = fmt(summary["posterior"], ".3f")
            print(f"OK ({elapsed:.0f}s, P={p_str}, ref={summary['ref_cat'] or '—'})")
            results.append((frb, "OK", f"{elapsed:.0f}s", summary))
        else:
            print(f"FAIL rc={proc.returncode}  (log: {log_path.relative_to(REPO_ROOT)})")
            results.append((frb, "FAIL", f"rc={proc.returncode}", summary))

    # ---- Summary table ----
    total_elapsed = time.time() - t0
    print()
    print(f"[batch] finished in {total_elapsed:.0f}s "
          f"({sum(1 for r in results if r[1] == 'OK')} ok, "
          f"{sum(1 for r in results if r[1] == 'FAIL')} fail, "
          f"{sum(1 for r in results if r[1] == 'SKIP')} skip)")
    print()
    print(f"{'FRB':<14} {'status':<6} {'P(O)':>7} {'ZP_40':>8} {'N*':>4} "
          f"{'ref_catalog':<24} note")
    print("-" * 80)
    for frb, status, note, summary in results:
        p = fmt(summary.get("posterior"), ".3f")
        zp = fmt(summary.get("zp_40px"), ".3f")
        n = fmt(summary.get("n_stars"))
        ref = summary.get("ref_cat") or ""
        print(f"{frb:<14} {status:<6} {p:>7} {zp:>8} {n:>4} {ref:<24} {note}")


if __name__ == "__main__":
    main()
