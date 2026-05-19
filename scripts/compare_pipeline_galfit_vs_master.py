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
3. Writes `pipeline_galfit_results.csv` to the repo root with flat columns.
4. Joins on `frb` against `master_frb_galfit_from_logs.csv` (only the `_psf`
   columns — those are the only fits the new pipeline produces) and writes
   `pipeline_vs_master_galfit_diff.csv`.
5. Prints a per-parameter summary (count, median |Δ|, std |Δ|, top-N largest
   deviations) for **single-Sérsic** fits only.

Comparison policy
-----------------
* **Magnitude / flux** is **not** compared (pipeline uses per-field ``zp_aper_40px``;
  master uses mixed ``J)`` systems). ``mag`` is still written to
  ``pipeline_galfit_results.csv`` for reference only.
* ``n_sersic_components`` counts fitted galaxies (no sky), from
  ``host_components.csv`` when present. ``compare_ok`` is true only when
  ``n_sersic_components == 1``; multi-component stamps are deblends where
  host-vs-master shape comparison is not meaningful.
* Summary statistics use ``compare_ok`` rows only.

By default omits 20171020A, 20220509G, 20240210A from the written CSVs (use
``--no-benchmark-exclusions`` to keep them). Parses the host as the first
``sersic`` line (``sersic_component_index=0``), matching GALFIT component 1.

Run from the repo root:
    python scripts/compare_pipeline_galfit_vs_master.py
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from galfit_fitlog_parse import (  # noqa: E402  reuse the canonical parser
    count_fitted_sersic_components,
    inclination_err_from_b_a_err,
    inclination_from_b_a,
    parse_fitlog_file,
)

REPO_ROOT = _SCRIPT_DIR.parent

# Divergent / non-comparable GALFIT logs dropped from the published benchmark table.
EXCLUDED_FROM_BENCHMARK = frozenset({"20171020A", "20220509G", "20240210A"})

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


def parse_pipeline_outputs(
    output_root: Path, tag: str, skip_frbs: frozenset[str]
) -> pd.DataFrame:
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
        if frb in skip_frbs:
            continue
        log_path = sub / "fit.log"
        n_sersic = count_fitted_sersic_components(sub)
        if not log_path.is_file():
            rows.append({
                "frb": frb,
                "n_sersic_components": n_sersic if n_sersic is not None else "",
                "compare_ok": False,
                "fit_log_path": "",
                "parse_strategy": "missing_fit_log",
                **{k: "" for k in _METRIC_KEYS},
                "inc": "", "inc_err": "",
            })
            continue

        # sersic_component_index=0: FRB host is GALFIT component 1 (first sersic
        # line); see generate_galfit_cutouts.py row order + galfit_fitlog_parse docs.
        data, strategy = parse_fitlog_file(str(log_path), sersic_component_index=0)
        try:
            rel_path = str(log_path.resolve().relative_to(REPO_ROOT)).replace("\\", "/")
        except ValueError:
            rel_path = str(log_path).replace("\\", "/")

        compare_ok = n_sersic == 1
        row: dict = {
            "frb": frb,
            "n_sersic_components": n_sersic if n_sersic is not None else "",
            "compare_ok": compare_ok,
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
        rows.append(row)

    if not rows:
        return pd.DataFrame()

    cols = [
        "frb", "n_sersic_components", "compare_ok", *_METRIC_KEYS,
        "inc", "inc_err", "parse_strategy", "fit_log_path",
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
    meta_cols = ["frb", "n_sersic_components", "compare_ok"]
    keep_pipe_cols = [c for c in meta_cols + list(rename_pipeline.keys()) if c in pipeline_df.columns]
    pipe_slim = pipeline_df[keep_pipe_cols].copy()
    pipe_slim = pipe_slim.rename(columns=rename_pipeline)

    diff = pipe_slim.merge(master_slim, on="frb", how="inner")
    if diff.empty:
        return diff

    if "compare_ok" in diff.columns:
        diff["compare_ok"] = diff["compare_ok"].astype(str).str.lower().isin(
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
    if "compare_ok" in diff.columns:
        comp = diff[diff["compare_ok"]].copy()
        n_comp = _to_numeric(diff["n_sersic_components"])
        n_multi = int((n_comp > 1).sum()) if n_comp is not None else 0
        print(
            f"[*] Comparison subset: {len(comp)} single-Sérsic FRBs "
            f"(excluded {full_n - len(comp)} multi/missing; {n_multi} with n_sersic>1)"
        )
    else:
        comp = diff

    if comp.empty:
        print("[!] No single-Sérsic FRBs to compare.")
        return

    print()
    print(f"Comparison summary (N = {len(comp)} comparable FRBs)")
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
    parser.add_argument(
        "--no-benchmark-exclusions",
        action="store_true",
        help=f"Keep all FRBs (default: omit {sorted(EXCLUDED_FROM_BENCHMARK)} from outputs).",
    )
    args = parser.parse_args()

    skip = frozenset() if args.no_benchmark_exclusions else EXCLUDED_FROM_BENCHMARK
    if skip:
        print(f"[*] Omitting from benchmark outputs: {sorted(skip)}")

    print(f"[*] Parsing pipeline outputs under {args.output_root} (tag='{args.tag}')")
    pipe_df = parse_pipeline_outputs(args.output_root, args.tag, skip)
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
    n_ok = int(diff["compare_ok"].sum()) if "compare_ok" in diff.columns else len(diff)
    print(
        f"[*] Wrote {args.out_diff_csv}  ({len(diff)} matched FRBs; "
        f"{n_ok} with compare_ok=True)"
    )
    print("[*] Deltas exclude mag/flux (per-field ZP vs legacy master conventions).")

    summarise_diff(diff, top_n=args.top_n)


if __name__ == "__main__":
    main()
