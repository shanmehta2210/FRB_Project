"""
Parse GALFIT fit.log files produced by the new pipeline and benchmark them
against the published master_frb_galfit_from_logs.csv values.

What it does
------------
1. Walks `pipeline_scripts/Output/<FRB>_<tag>/fit.log` for every finished run
   (default tag = `all`; the FRB name is everything before the last `_<tag>`).
2. Extracts the with-PSF Sérsic block using `galfit_fitlog_parse.parse_fitlog_file`
   — the *same* parser `build_master_frb_galfit_from_logs.py` uses, so the
   block-selection policy (sane single-sérsic, free-n-before-fixed-n refine,
   etc.) matches the published numbers.
3. Writes `pipeline_galfit_results.csv` to the repo root with flat columns
   (includes `host_number`, `snr_win`, `snr_auto` from `sky_fit_audit.json` or
   `host_components.csv` row 0).
4. Joins on `frb` against `master_frb_galfit_from_logs.csv` (only the `_psf`
   columns — those are the only fits the new pipeline produces) and writes
   `pipeline_vs_master_galfit_diff.csv`.
5. Prints a per-parameter summary (count, median |Δ|, std |Δ|, top-N largest
   deviations) for every FRB with a parsed **host** row.

Comparison policy
-----------------
* **Host identification:** always GALFIT **component 1** — first ``sersic`` line in
  the selected ``fit.log`` block (``sersic_component_index=0``), matching row 0 in
  ``host_components.csv``. Neighbor Sérsics are ignored for results/deltas.
* **Magnitude / flux** is **not** compared (pipeline uses per-field ``zp_aper``;
  master uses mixed ``J)`` systems). ``mag`` is still written to
  ``pipeline_galfit_results.csv`` for reference only.
* ``n_sersic_components`` counts fitted galaxies (no sky). ``single_sersic`` is
  **informational only** (true when exactly one Sérsic was fit); multi-component
  stamps still export host structural parameters (especially ``b/a`` / inclination).
* Summary statistics include **all** matched FRBs with a parsed host.

Run from the repo root:
    python scripts/compare_pipeline_galfit_vs_master.py
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))
_PIPELINE_DIR = _SCRIPT_DIR.parent / "pipeline_scripts"
if str(_PIPELINE_DIR) not in sys.path:
    sys.path.insert(0, str(_PIPELINE_DIR))

from galfit_fitlog_parse import (  # noqa: E402  reuse the canonical parser
    count_fitted_sersic_components,
    inclination_err_from_b_a_err,
    inclination_from_b_a,
    parse_fitlog_file,
)
from reference_photometry import (  # noqa: E402
    final_host_mag,
    galfit_zp_used,
    load_photometry,
    pipeline_zp_ok,
)

REPO_ROOT = _SCRIPT_DIR.parent

MASTER_RUNS_ROOT = REPO_ROOT / "tools" / "galfit" / "runs"

_METRIC_KEYS = [
    "chi2nu",
    "mag", "mag_err",
    "re", "re_err",
    "n", "n_err",
    "b_a", "b_a_err",
    "pa", "pa_err",
    "x", "x_err",
    "y", "y_err",
]

# Parameters we actually want to compare scientifically. (x/y/x_err/y_err are
# cutout-frame pixel coordinates; comparing them apples-to-apples requires
# matching cutout origins, which is out of scope here.)
_COMPARE_KEYS = ["chi2nu", "re", "n", "b_a", "pa", "inc"]


def _cell(val):
    if val is None:
        return ""
    if isinstance(val, float) and (math.isnan(val) or not math.isfinite(val)):
        return ""
    return val


def _split_tag(folder_name: str, tag: str) -> str | None:
    """Return the FRB name for a folder named like '<FRB>_<tag>', else None."""
    suffix = f"_{tag}"
    if folder_name.endswith(suffix):
        return folder_name[: -len(suffix)]
    return None


def _host_snr_from_components(output_dir: Path) -> dict[str, float | int | None]:
    """Host SNR from host_components.csv row 0 (Phase 3a)."""
    comp_path = output_dir / "host_components.csv"
    out: dict[str, float | int | None] = {
        "snr_win": None,
        "snr_auto": None,
        "host_number": None,
    }
    if not comp_path.is_file():
        return out
    try:
        comp = pd.read_csv(comp_path, nrows=1)
    except (OSError, ValueError):
        return out
    if comp.empty:
        return out
    host = comp.iloc[0]
    if "NUMBER" in comp.columns and pd.notna(host["NUMBER"]):
        out["host_number"] = int(host["NUMBER"])
    if "SNR_WIN" in comp.columns and pd.notna(host["SNR_WIN"]):
        snr_win = float(host["SNR_WIN"])
        if math.isfinite(snr_win):
            out["snr_win"] = snr_win
    if "FLUX_AUTO" in comp.columns and "FLUXERR_AUTO" in comp.columns:
        flux = pd.to_numeric(host["FLUX_AUTO"], errors="coerce")
        ferr = pd.to_numeric(host["FLUXERR_AUTO"], errors="coerce")
        if pd.notna(flux) and pd.notna(ferr) and float(ferr) > 0:
            snr_auto = float(flux) / float(ferr)
            if math.isfinite(snr_auto):
                out["snr_auto"] = snr_auto
    return out


def _host_snr_from_audit(output_dir: Path) -> dict[str, float | int | None]:
    """Read host SNR fields from sky_fit_audit.json, else host_components.csv."""
    out: dict[str, float | int | None] = {
        "snr_win": None,
        "snr_auto": None,
        "host_number": None,
    }
    audit_path = output_dir / "sky_fit_audit.json"
    if audit_path.is_file():
        try:
            with audit_path.open(encoding="utf-8") as f:
                audit = json.load(f)
        except (json.JSONDecodeError, OSError):
            audit = {}
        for key in ("snr_win", "snr_auto"):
            val = audit.get(key)
            if val is not None:
                try:
                    fval = float(val)
                    if math.isfinite(fval):
                        out[key] = fval
                except (TypeError, ValueError):
                    pass
        hn = audit.get("host_number")
        if hn is not None:
            try:
                out["host_number"] = int(hn)
            except (TypeError, ValueError):
                pass
    fallback = _host_snr_from_components(output_dir)
    for key in out:
        if out[key] is None and fallback[key] is not None:
            out[key] = fallback[key]
    return out


def _reference_mag_fields(output_dir: Path, galfit_mag, galfit_mag_err) -> dict:
    """Reference-survey mag + final-mag substitution for one Output folder.

    ``mag_final`` is the GALFIT magnitude when the pipeline zero-point is
    trusted, otherwise the LS DR10 / PS1 reference magnitude stored by
    ``pipeline_scripts/reference_photometry.py`` (run automatically by
    master_run.py; backfill older folders with ``--backfill``).
    """
    ref = None
    ref_path = output_dir / "reference_photometry.json"
    if ref_path.is_file():
        try:
            ref = json.loads(ref_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            ref = None
    if ref is None:
        summary_path = output_dir / "pipeline_summary.json"
        if summary_path.is_file():
            try:
                summary = json.loads(summary_path.read_text(encoding="utf-8"))
                ref = summary.get("reference_photometry")
            except (json.JSONDecodeError, OSError):
                pass

    # pipeline_summary.json photometry section, or zero_points.json for
    # pre-summary Output folders.
    photometry = load_photometry(output_dir)

    ref_ok = isinstance(ref, dict) and ref.get("status") == "ok"
    mag_final, mag_final_err, source = final_host_mag(
        photometry, ref, galfit_mag, galfit_mag_err,
        zp_used=galfit_zp_used(output_dir),
    )
    return {
        "zp_ok": pipeline_zp_ok(photometry),
        "ref_survey": _cell(ref.get("survey")) if ref_ok else "",
        "ref_mag": _cell(ref.get("mag")) if ref_ok else "",
        "ref_mag_err": _cell(ref.get("mag_err")) if ref_ok else "",
        "ref_sep_arcsec": _cell(ref.get("sep_arcsec")) if ref_ok else "",
        "mag_final": _cell(mag_final),
        "mag_final_err": _cell(mag_final_err),
        "mag_final_source": source,
    }


def parse_pipeline_outputs(output_root: Path, tag: str) -> pd.DataFrame:
    """Walk Output/<FRB>_<tag>/ and parse every fit.log we find."""
    rows: list[dict] = []
    if not output_root.is_dir():
        print(f"[!] {output_root} does not exist — nothing to parse.")
        return pd.DataFrame()

    for sub in sorted(output_root.iterdir()):
        if not sub.is_dir():
            continue
        frb = _split_tag(sub.name, tag)
        if frb is None:
            continue
        log_path = sub / "fit.log"
        n_sersic = count_fitted_sersic_components(sub)
        host_snr = _host_snr_from_audit(sub)
        snr_fields = {
            "snr_win": _cell(host_snr["snr_win"]),
            "snr_auto": _cell(host_snr["snr_auto"]),
            "host_number": _cell(host_snr["host_number"]),
        }
        if not log_path.is_file():
            rows.append({
                "frb": frb,
                "n_sersic_components": n_sersic if n_sersic is not None else "",
                "single_sersic": False,
                "fit_log_path": "",
                "parse_strategy": "missing_fit_log",
                **snr_fields,
                **{k: "" for k in _METRIC_KEYS},
                "inc": "", "inc_err": "",
                **_reference_mag_fields(sub, None, None),
            })
            continue

        # sersic_component_index=0: FRB host is GALFIT component 1 (first sersic
        # line); see generate_galfit_cutouts.py row order + galfit_fitlog_parse docs.
        data, strategy = parse_fitlog_file(str(log_path), sersic_component_index=0)
        try:
            rel_path = str(log_path.resolve().relative_to(REPO_ROOT)).replace("\\", "/")
        except ValueError:
            rel_path = str(log_path).replace("\\", "/")

        single_sersic = n_sersic == 1
        row: dict = {
            "frb": frb,
            "n_sersic_components": n_sersic if n_sersic is not None else "",
            "single_sersic": single_sersic,
            **snr_fields,
            "fit_log_path": rel_path,
            "parse_strategy": strategy,
        }
        for k in _METRIC_KEYS:
            row[k] = _cell(data.get(k))

        ba = data.get("b_a")
        be = data.get("b_a_err")
        row["inc"] = _cell(inclination_from_b_a(ba)) if ba is not None else ""
        row["inc_err"] = (
            _cell(inclination_err_from_b_a_err(ba, be)) if ba is not None else ""
        )
        row.update(_reference_mag_fields(sub, data.get("mag"), data.get("mag_err")))
        rows.append(row)

    if not rows:
        return pd.DataFrame()

    cols = [
        "frb",
        "n_sersic_components",
        "single_sersic",
        "host_number",
        "snr_win",
        "snr_auto",
        *_METRIC_KEYS,
        "inc",
        "inc_err",
        "zp_ok",
        "ref_survey",
        "ref_mag",
        "ref_mag_err",
        "ref_sep_arcsec",
        "mag_final",
        "mag_final_err",
        "mag_final_source",
        "parse_strategy",
        "fit_log_path",
    ]
    return pd.DataFrame(rows)[cols]


def _to_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def build_diff_table(pipeline_df: pd.DataFrame, master_df: pd.DataFrame) -> pd.DataFrame:
    """Inner-join on `frb`, keep only the with-PSF master columns, compute deltas.

    Drops FRBs whose master row is present but has no numeric with-PSF data —
    those are radio-only localisations carried in master_frb_galfit_from_logs.csv
    as empty placeholders (e.g. coord_semantics != 'host'). Including them
    would silently inflate the "N matched" count without contributing to
    any per-parameter comparison.
    """
    keep_master_cols = ["frb"] + [
        c for k in _COMPARE_KEYS for c in (
            (f"{k}_psf", f"{k}_err_psf") if k != "inc"
            else ("inc_psf", "inc_err_psf")
        )
        if c in master_df.columns
    ]
    master_slim = master_df[keep_master_cols].copy()

    rename_master = {c: c.replace("_psf", "_master") for c in master_slim.columns if c != "frb"}
    master_slim = master_slim.rename(columns=rename_master)

    rename_pipeline = {}
    for k in _COMPARE_KEYS:
        rename_pipeline[k] = f"{k}_pipeline"
        if k == "inc":
            rename_pipeline["inc_err"] = "inc_err_pipeline"
        else:
            rename_pipeline[f"{k}_err"] = f"{k}_err_pipeline"
    meta_cols = ["frb", "n_sersic_components", "single_sersic"]
    keep_pipe_cols = [c for c in meta_cols + list(rename_pipeline.keys()) if c in pipeline_df.columns]
    pipe_slim = pipeline_df[keep_pipe_cols].copy()
    pipe_slim = pipe_slim.rename(columns=rename_pipeline)

    diff = pipe_slim.merge(master_slim, on="frb", how="inner")
    if diff.empty:
        return diff

    if "single_sersic" in diff.columns:
        diff["single_sersic"] = diff["single_sersic"].astype(str).str.lower().isin(
            ("true", "1", "yes")
        )

    ordered_cols = [c for c in meta_cols if c in diff.columns]
    for k in _COMPARE_KEYS:
        a = f"{k}_pipeline"
        b = f"{k}_master"
        if a in diff.columns and b in diff.columns:
            diff[a] = _to_numeric(diff[a])
            diff[b] = _to_numeric(diff[b])
            diff[f"{k}_delta"] = diff[a] - diff[b]
            ordered_cols += [a, b, f"{k}_delta"]
    diff = diff[ordered_cols]

    master_numeric_cols = [f"{k}_master" for k in _COMPARE_KEYS if f"{k}_master" in diff.columns]
    if master_numeric_cols:
        has_master_fit = diff[master_numeric_cols].notna().any(axis=1)
        n_dropped = int((~has_master_fit).sum())
        if n_dropped:
            dropped_frbs = diff.loc[~has_master_fit, "frb"].tolist()
            print(
                f"[*] Dropped {n_dropped} FRB(s) from diff: master row exists but "
                f"has no numeric with-PSF fit columns — {dropped_frbs}"
            )
        diff = diff[has_master_fit].reset_index(drop=True)
    return diff


def summarise_diff(diff: pd.DataFrame, top_n: int) -> None:
    if diff.empty:
        print("[!] No overlapping FRBs between pipeline outputs and master CSV.")
        return

    full_n = len(diff)
    if "chi2nu_pipeline" in diff.columns:
        comp = diff[_to_numeric(diff["chi2nu_pipeline"]).notna()].copy()
    else:
        comp = diff.copy()

    n_multi = 0
    if "n_sersic_components" in diff.columns:
        n_comp = _to_numeric(diff["n_sersic_components"])
        n_multi = int((n_comp > 1).sum()) if n_comp is not None else 0
    n_single = int(diff["single_sersic"].sum()) if "single_sersic" in diff.columns else 0
    print(
        f"[*] Host metrics: GALFIT component 1 (parser index 0) for all {len(comp)} FRBs "
        f"with a parsed fit ({full_n - len(comp)} skipped: no host chi2nu). "
        f"{n_multi} multi-Sérsic ({n_single} single-Sérsic)."
    )

    if comp.empty:
        print("[!] No FRBs with parsed host parameters to compare.")
        return

    print()
    print(f"Comparison summary (N = {len(comp)} FRBs, host component 1)")
    print("-" * 96)
    header = (
        f"{'param':>6}  {'N':>4}  {'median(d)':>12}  {'median|d|':>12}  "
        f"{'std(d)':>12}  {'max|d|':>12}  {'frb @ max|d|':>14}"
    )
    print(header)
    print("-" * len(header))
    for k in _COMPARE_KEYS:
        col = f"{k}_delta"
        if col not in diff.columns:
            continue
        vals = comp[col].dropna()
        if vals.empty:
            print(f"{k:>6}  {0:>4}  {'-':>12}  {'-':>12}  {'-':>12}  {'-':>12}  {'-':>14}")
            continue
        idx = vals.abs().idxmax()
        print(
            f"{k:>6}  "
            f"{len(vals):>4}  "
            f"{vals.median():>12.4f}  "
            f"{vals.abs().median():>12.4f}  "
            f"{vals.std():>12.4f}  "
            f"{vals.abs().max():>12.4f}  "
            f"{comp.loc[idx, 'frb']:>14}"
        )

    print()
    print(f"Top {top_n} largest |delta| per parameter")
    print("-" * 96)
    for k in _COMPARE_KEYS:
        col = f"{k}_delta"
        if col not in diff.columns:
            continue
        ranked = (
            comp.assign(_abs=comp[col].abs())
            .dropna(subset=[col])
            .sort_values("_abs", ascending=False)
        )
        if ranked.empty:
            continue
        print(f"  {k}:")
        for _, row in ranked.head(top_n).iterrows():
            print(
                f"    {row['frb']}: pipeline={row[f'{k}_pipeline']:.4f}, "
                f"master={row[f'{k}_master']:.4f}, delta={row[col]:+.4f}"
            )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-root",
        type=Path,
        default=REPO_ROOT / "pipeline_scripts" / "Output",
        help="Directory containing per-FRB pipeline output folders (default: pipeline_scripts/Output).",
    )
    parser.add_argument(
        "--tag",
        type=str,
        default="all",
        help="Folder-name suffix produced by master_run.py (default: 'all').",
    )
    parser.add_argument(
        "--master-csv",
        type=Path,
        default=REPO_ROOT / "master_frb_galfit_from_logs.csv",
        help="Path to master_frb_galfit_from_logs.csv (default: repo root).",
    )
    parser.add_argument(
        "--out-pipeline-csv",
        type=Path,
        default=REPO_ROOT / "pipeline_galfit_results.csv",
        help="Where to write the flat pipeline GALFIT results (default: repo root).",
    )
    parser.add_argument(
        "--out-diff-csv",
        type=Path,
        default=REPO_ROOT / "pipeline_vs_master_galfit_diff.csv",
        help="Where to write the per-FRB pipeline-vs-master delta table.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=5,
        help="How many of the largest |Δ| per parameter to print (default: 5).",
    )
    args = parser.parse_args()

    print(f"[*] Parsing pipeline outputs under {args.output_root} (tag='{args.tag}')")
    pipe_df = parse_pipeline_outputs(args.output_root, args.tag)
    if pipe_df.empty:
        print("[!] No FRB folders matched — nothing to do.")
        return

    pipe_df.to_csv(args.out_pipeline_csv, index=False)
    n_with_fit = (pipe_df["parse_strategy"] != "missing_fit_log").sum()
    print(f"[*] Wrote {args.out_pipeline_csv}  ({len(pipe_df)} FRBs, {n_with_fit} with fit.log)")

    if not args.master_csv.is_file():
        print(f"[!] {args.master_csv} not found — skipping comparison step.")
        return

    master_df = pd.read_csv(args.master_csv)
    diff = build_diff_table(pipe_df, master_df)
    if diff.empty:
        print("[!] No FRBs in common between pipeline and master CSV.")
        return

    diff.to_csv(args.out_diff_csv, index=False)
    n_single = int(diff["single_sersic"].sum()) if "single_sersic" in diff.columns else 0
    n_host = int(_to_numeric(diff.get("chi2nu_pipeline", pd.Series())).notna().sum())
    print(
        f"[*] Wrote {args.out_diff_csv}  ({len(diff)} matched FRBs; "
        f"{n_host} with host fit parsed; {n_single} single-Sérsic)"
    )
    print("[*] Deltas exclude mag/flux (per-field ZP vs legacy master conventions).")

    summarise_diff(diff, top_n=args.top_n)


if __name__ == "__main__":
    main()
