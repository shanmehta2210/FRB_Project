#!/usr/bin/env python3
"""Verify SDSS null pools: every galaxy has lnLExp_r > lnLDeV_r and uses expAB_r."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (
    SDSS_UR_MAX_CDF,
    cosi_array_from_df,
    prepare_null_strict_color_base,
    read_sdss_null_catalog,
    sdss_exp_wins_lnl_mask,
    slice_null_base_by_mag,
)

OUT_DIR = (
    Path(__file__).resolve().parents[1]
    / "plots"
    / "plots_null"
    / "v1_null_cdf_inclination"
    / "diagnostics"
    / "sdss_profile_winner"
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sdss-csv", default="SDSS_catalog_v1_allsky_modelmr.csv")
    parser.add_argument("--mag-limit", type=float, default=21.0)
    args = parser.parse_args()

    df = read_sdss_null_catalog(args.sdss_csv)
    n_cat = len(df)
    has_lnl = sdss_exp_wins_lnl_mask(df)
    print(f"Catalog N={n_cat}")
    print(f"  finite lnL: {has_lnl.notna().sum() if hasattr(has_lnl, 'sum') else 0}")  # noqa
    ln_exp = pd.to_numeric(df.get("lnLExp_r"), errors="coerce")
    ln_dev = pd.to_numeric(df.get("lnLDeV_r"), errors="coerce")
    finite = ln_exp.notna() & ln_dev.notna()
    print(f"  finite lnL pairs: {finite.sum()} ({finite.mean():.1%})")
    print(f"  exp wins (catalog): {(ln_exp > ln_dev).sum()} ({(ln_exp > ln_dev).sum() / max(finite.sum(),1):.1%} of finite)")

    base = prepare_null_strict_color_base(
        df,
        mag_column="modelMag_r",
        q_column="expAB_r",
        is_legacy=False,
        sdss_ur_max=SDSS_UR_MAX_CDF,
        sdss_exp_winner_only=True,
    )
    cut = slice_null_base_by_mag(base, mag_column="modelMag_r", mag_limit=args.mag_limit)
    le = pd.to_numeric(cut["lnLExp_r"], errors="coerce")
    ld = pd.to_numeric(cut["lnLDeV_r"], errors="coerce")

    violations = (le <= ld) | le.isna() | ld.isna()
    print(f"\nCDF pool m<{args.mag_limit}: N={len(cut)}")
    print(f"  lnL violations: {violations.sum()}")
    if violations.any():
        raise SystemExit("[FAIL] deV winners or missing lnL in CDF pool")

    cosi = cosi_array_from_df(cut, q_col="expAB_r")
    print(f"  mean cos(i): {np.mean(cosi):.4f}  median: {np.median(cosi):.4f}")
    print(f"  frac cos(i)>0.5: {(cosi > 0.5).mean():.3f}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {"metric": "catalog_n", "value": n_cat},
            {"metric": "finite_lnl_frac", "value": float(finite.mean())},
            {"metric": f"pool_n_mag_lt_{args.mag_limit}", "value": len(cut)},
            {"metric": "pool_lnl_violations", "value": int(violations.sum())},
            {"metric": "pool_mean_cosi", "value": float(np.mean(cosi))},
            {"metric": "pool_median_cosi", "value": float(np.median(cosi))},
        ]
    ).to_csv(OUT_DIR / "pool_audit.csv", index=False)
    print(f"Wrote {OUT_DIR / 'pool_audit.csv'}")
    print("[PASS] All SDSS CDF-pool galaxies have lnLExp_r > lnLDeV_r")


if __name__ == "__main__":
    main()
