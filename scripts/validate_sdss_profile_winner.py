#!/usr/bin/env python3
"""
Validate SDSS r-band profile winner: lnL vs mag-proxy vs fracDeV.

Writes diagnostics under plots/plots_null/.../sdss_profile_winner/.

Run from repo root (after patch_sdss_profile_winner.py):
    python scripts/validate_sdss_profile_winner.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (  # noqa: E402
    assign_sdss_profile_winner_columns,
    mag_proxy_sdss_exp_winner,
)
from pipeline_null_plot_utils import DEFAULT_SDSS, PLOTS_NULL  # noqa: E402

OUT_DIR = (
    PLOTS_NULL
    / "v1_null_cdf_inclination"
    / "diagnostics"
    / "sdss_profile_winner"
)
SAMPLE_N = 50_000
SEED = 42


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sdss-csv", type=Path, default=DEFAULT_SDSS)
    parser.add_argument("--sample-n", type=int, default=SAMPLE_N)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    usecols = [
        "modelMag_r",
        "deVMag_r",
        "expMag_r",
        "fracDeV_r",
        "lnLDeV_r",
        "lnLExp_r",
        "model_winner_is_exp",
    ]
    df = pd.read_csv(args.sdss_csv, usecols=lambda c: c in usecols or c == "rmag")
    if "lnLExp_r" not in df.columns:
        raise SystemExit(
            "Catalog missing lnL columns; run scripts/patch_sdss_profile_winner.py first."
        )

    if len(df) > args.sample_n:
        df = df.sample(n=args.sample_n, random_state=SEED)

    df = assign_sdss_profile_winner_columns(df)
    lnl_exp = df["model_winner_is_exp"].astype(bool)
    mag_exp = mag_proxy_sdss_exp_winner(df).astype(bool)
    frac = pd.to_numeric(df["fracDeV_r"], errors="coerce")
    frac_exp = frac < 0.5

    valid = lnl_exp.notna() & mag_exp.notna() & frac.notna()
    lnl_exp = lnl_exp[valid]
    mag_exp = mag_exp[valid]
    frac_exp = frac_exp[valid]

    agree_lnl_mag = (lnl_exp == mag_exp).mean()
    agree_lnl_frac = (lnl_exp == frac_exp).mean()
    agree_mag_frac = (mag_exp == frac_exp).mean()

    cross = pd.crosstab(
        pd.Series(lnl_exp, name="lnL_exp_wins"),
        pd.Series(mag_exp, name="mag_proxy_exp_wins"),
    )

    rows = [
        {"metric": "n_sample", "value": int(valid.sum())},
        {"metric": "frac_lnl_exp_winner", "value": float(lnl_exp.mean())},
        {"metric": "frac_mag_proxy_exp_winner", "value": float(mag_exp.mean())},
        {"metric": "frac_fracDeV_exp", "value": float(frac_exp.mean())},
        {"metric": "agreement_lnl_vs_mag_proxy", "value": float(agree_lnl_mag)},
        {"metric": "agreement_lnl_vs_fracDeV_lt_half", "value": float(agree_lnl_frac)},
        {"metric": "agreement_mag_proxy_vs_fracDeV", "value": float(agree_mag_frac)},
    ]
    summary = pd.DataFrame(rows)
    summary.to_csv(args.out_dir / "validation_metrics.csv", index=False)
    cross.to_csv(args.out_dir / "crosstab_lnl_vs_mag_proxy.csv")

    md_lines = [
        "# SDSS profile winner validation",
        "",
        f"Sample: {int(valid.sum())} rows from `{args.sdss_csv.name}`.",
        "",
        "## Primary rule",
        "",
        "`model_winner_is_exp` = (`lnLExp_r` > `lnLDeV_r`)",
        "",
        "## Agreement",
        "",
        f"- lnL vs mag-proxy (|expMag−modelMag| vs |deVMag−modelMag|): **{agree_lnl_mag:.4f}**",
        f"- lnL vs fracDeV < 0.5: **{agree_lnl_frac:.4f}**",
        f"- mag-proxy vs fracDeV < 0.5: **{agree_mag_frac:.4f}**",
        "",
        f"Fraction exp-winner (lnL): {lnl_exp.mean():.3f}",
        "",
        "Production CDF pools use lnL only; mag-proxy and fracDeV are audit cross-checks.",
    ]
    (args.out_dir / "validation_summary.md").write_text("\n".join(md_lines), encoding="utf-8")

    print(summary.to_string(index=False))
    print(f"Wrote {args.out_dir}")


if __name__ == "__main__":
    main()
