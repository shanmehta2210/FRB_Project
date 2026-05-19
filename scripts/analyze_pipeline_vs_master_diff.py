"""Basic statistics on pipeline_vs_master_galfit_diff.csv (run from repo root)."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# Keep in sync with scripts/compare_pipeline_galfit_vs_master.py
BENCHMARK_EXCLUDED = frozenset({"20171020A", "20220509G", "20240210A"})


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path("pipeline_vs_master_galfit_diff.csv"),
        help="Path to pipeline_vs_master_galfit_diff.csv",
    )
    parser.add_argument(
        "--no-benchmark-filter",
        action="store_true",
        help="Do not drop the default benchmark-excluded FRBs when summarising.",
    )
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    path = args.csv if args.csv.is_absolute() else root / args.csv

    df = pd.read_csv(path)
    if not args.no_benchmark_filter and "frb" in df.columns:
        present_excluded = sorted(set(df["frb"].astype(str)) & BENCHMARK_EXCLUDED)
        df = df[~df["frb"].isin(BENCHMARK_EXCLUDED)].reset_index(drop=True)
        if present_excluded:
            print(
                f"[*] Dropped {len(present_excluded)} benchmark-excluded row(s): "
                f"{present_excluded}"
            )
    delta_cols = [c for c in df.columns if c.endswith("_delta")]
    n = len(df)

    print("=" * 80)
    print(f"{path.name} - basic statistics (N = {n} FRBs)")
    print("=" * 80)

    rows = []
    for col in delta_cols:
        x = pd.to_numeric(df[col], errors="coerce")
        valid = x.dropna()
        if len(valid) == 0:
            continue
        rows.append(
            {
                "param": col.replace("_delta", ""),
                "n": len(valid),
                "mean_d": valid.mean(),
                "std_d": valid.std(ddof=1),
                "median_d": valid.median(),
                "mad_d": (valid - valid.median()).abs().median(),
                "q25": valid.quantile(0.25),
                "q75": valid.quantile(0.75),
                "iqr": valid.quantile(0.75) - valid.quantile(0.25),
                "min_d": valid.min(),
                "max_d": valid.max(),
                "mean_abs": valid.abs().mean(),
                "median_abs": valid.abs().median(),
                "rmse": float(np.sqrt((valid**2).mean())),
            }
        )
    sum_df = pd.DataFrame(rows).set_index("param")
    pd.set_option("display.width", 200)
    pd.set_option("display.max_columns", 20)
    pd.set_option("display.float_format", lambda v: f"{v:.6g}")

    print()
    print("--- Raw delta statistics (all N) ---")
    print(sum_df.to_string())

    chi = pd.to_numeric(df["chi2nu_delta"], errors="coerce")
    blow = df.loc[chi.abs() > 1e3, "frb"].tolist()
    print()
    print("--- chi2nu_delta: FRBs with |delta| > 1000 ---")
    print(blow if blow else "(none)")

    mask_ok = chi.abs() <= 1e3
    df_robust = df.loc[mask_ok]
    print()
    print(
        f"--- Subset excluding chi2nu blow-ups (N = {int(mask_ok.sum())}) - structural ---"
    )
    struct = ["re_delta", "n_delta", "b_a_delta", "pa_delta", "inc_delta"]
    for col in struct:
        if col not in df_robust.columns:
            continue
        x = pd.to_numeric(df_robust[col], errors="coerce").dropna()
        print(
            f"  {col.replace('_delta', '')}:  mean={x.mean():.4f}  "
            f"std={x.std(ddof=1):.4f}  median={x.median():.4f}  "
            f"median|d|={x.abs().median():.4f}  RMSE={np.sqrt((x**2).mean()):.4f}"
        )

    mag_d = pd.to_numeric(df["mag_delta"], errors="coerce")
    print()
    print(
        "--- mag_delta (pipeline J)=22.5; master J)=25 hosts corrected by -2.5 in "
        "compare_pipeline_galfit_vs_master.py) ---"
    )
    print(f"  mean = {mag_d.mean():.3f}  median = {mag_d.median():.3f}  std = {mag_d.std(ddof=1):.3f}")

    sub = df[delta_cols].apply(pd.to_numeric, errors="coerce")
    sub_nochi = sub.drop(columns=["chi2nu_delta"], errors="ignore")
    corr = sub_nochi.corr(method="spearman")
    print()
    print("--- Spearman correlation among deltas (excl. chi2nu_delta) ---")
    print(corr.round(3).to_string())

    print()
    print("--- Top 3 largest |delta| per parameter ---")
    for col in delta_cols:
        x = pd.to_numeric(df[col], errors="coerce")
        t = (
            df.assign(_abs=x.abs(), _d=x)
            .dropna(subset=["_d"])
            .nlargest(3, "_abs")
        )
        print(f"  {col}:")
        for _, r in t.iterrows():
            print(f"    {r['frb']}: delta={r['_d']:.6g}")


if __name__ == "__main__":
    main()
