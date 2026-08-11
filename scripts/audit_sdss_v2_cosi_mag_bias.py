#!/usr/bin/env python3
"""Audit cos(i) vs magnitude for SDSS v2 — cut-order and per-bin checks."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (
    Q0,
    SDSS_UR_MAX_CDF,
    apply_strict_q_cut,
    cosi_array_from_df,
    filter_sdss_drop_dev_winners,
    filter_sdss_ur,
    prepare_null_strict_color_base,
    read_sdss_null_catalog,
    slice_null_base_by_mag,
)
from pipeline_null_plot_utils import DEFAULT_SDSS_V2, REPO_ROOT

OUT = REPO_ROOT / "plots" / "plots_null" / "v2" / "sdss_audit" / "cosi_mag_bias_audit.csv"


def pool_stats(sub: pd.DataFrame, q_col: str = "expAB_r") -> dict:
    cosi = cosi_array_from_df(sub, q_col=q_col, q0=Q0)
    ba = pd.to_numeric(sub[q_col], errors="coerce")
    return {
        "n": len(sub),
        "mean_cosi": float(np.mean(cosi)),
        "median_cosi": float(np.median(cosi)),
        "median_expAB_r": float(ba.median()),
    }


def main() -> None:
    df = read_sdss_null_catalog(DEFAULT_SDSS_V2)
    rows: list[dict] = []

    # Order A: evolution diagnostic (strict -> ur -> lnL)
    s_a = apply_strict_q_cut(df, q_col="expAB_r", q0=Q0)
    s_a = filter_sdss_ur(s_a, SDSS_UR_MAX_CDF)
    s_a = filter_sdss_drop_dev_winners(s_a)
    for label, sub in [
        ("A_strict_ur_lnl", s_a),
        ("A_mag21", slice_null_base_by_mag(s_a, mag_column="modelMag_r", mag_limit=21)),
        ("A_mag22", slice_null_base_by_mag(s_a, mag_column="modelMag_r", mag_limit=22)),
    ]:
        rows.append({"stage": label, **pool_stats(sub)})

    # Order B: production (ur -> lnL -> strict)
    base_b = prepare_null_strict_color_base(
        df,
        mag_column="modelMag_r",
        q0=Q0,
        q_column="expAB_r",
        is_legacy=False,
        sdss_ur_max=SDSS_UR_MAX_CDF,
        sdss_exp_winner_only=True,
    )
    for label, sub in [
        ("B_ur_lnl_strict", base_b),
        ("B_mag21", slice_null_base_by_mag(base_b, mag_column="modelMag_r", mag_limit=21)),
        ("B_mag22", slice_null_base_by_mag(base_b, mag_column="modelMag_r", mag_limit=22)),
    ]:
        rows.append({"stage": label, **pool_stats(sub)})

    # Per 0.5 mag bin on production pool (not cumulative)
    mag = pd.to_numeric(base_b["modelMag_r"], errors="coerce")
    cosi = cosi_array_from_df(base_b, q_col="expAB_r", q0=Q0)
    for lo in np.arange(15, 24.5, 0.5):
        hi = lo + 0.5
        mask = (mag > lo) & (mag <= hi)
        n = int(mask.sum())
        if n < 50:
            continue
        rows.append(
            {
                "stage": f"bin_{lo:.1f}_{hi:.1f}",
                "n": n,
                "mean_cosi": float(cosi[mask].mean()),
                "median_cosi": float(np.median(cosi[mask])),
                "median_expAB_r": float(base_b.loc[mask, "expAB_r"].median()),
            }
        )

    # Cumulative mag limits (nested subsets)
    for mlim in [17, 18, 19, 20, 21, 22, 23]:
        sub = slice_null_base_by_mag(base_b, mag_column="modelMag_r", mag_limit=mlim)
        rows.append({"stage": f"cumulative_mag<={mlim}", **pool_stats(sub)})

    # Isotropic simulation: random cos(i) uniform, fixed b/a distribution from pool
    rng = np.random.default_rng(42)
    ba_samp = pd.to_numeric(base_b["expAB_r"], errors="coerce").to_numpy()
    ba_samp = ba_samp[np.isfinite(ba_samp)]
    sim_cosi = []
    for _ in range(50_000):
        inc = rng.uniform(0, 90)  # isotropic in solid angle -> uniform cos(i)? 
        # isotropic orientation: cos(i) uniform on [0,1]
    cosi_iso = rng.uniform(0, 1, size=50_000)
    # For each random cos(i), implied b/a from Hubble inverse
    q0 = Q0
    ba_implied = np.sqrt(cosi_iso**2 * (1 - q0**2) + q0**2)
    rows.append(
        {
            "stage": "sim_isotropic_uniform_cosi",
            "n": 50_000,
            "mean_cosi": float(np.mean(cosi_iso)),
            "median_cosi": float(np.median(cosi_iso)),
            "median_expAB_r": float(np.median(ba_implied)),
        }
    )

    # Wrong column sanity
    if "best_model_ba_r" in df.columns:
        wrong = prepare_null_strict_color_base(
            df,
            mag_column="modelMag_r",
            q0=Q0,
            q_column="best_model_ba_r",
            is_legacy=False,
            sdss_ur_max=SDSS_UR_MAX_CDF,
            sdss_exp_winner_only=True,
        )
        sub = slice_null_base_by_mag(wrong, mag_column="modelMag_r", mag_limit=21)
        rows.append(
            {
                "stage": "WRONG_best_model_ba_mag21",
                **pool_stats(sub, q_col="best_model_ba_r"),
            }
        )

    # Raw catalog: mag vs expAB_r before ur/lnL
    mag_raw = pd.to_numeric(df["modelMag_r"], errors="coerce")
    ba_raw = pd.to_numeric(df["expAB_r"], errors="coerce")
    ok = mag_raw.notna() & ba_raw.notna() & (ba_raw > Q0)
    corr = float(np.corrcoef(mag_raw[ok], ba_raw[ok])[0, 1])
    rows.append(
        {
            "stage": "raw_expAB_gt_q0_corr_mag_ba",
            "n": int(ok.sum()),
            "mean_cosi": corr,  # store correlation in mean_cosi column
            "median_cosi": float("nan"),
            "median_expAB_r": float(ba_raw[ok].median()),
        }
    )

    out_df = pd.DataFrame(rows)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(OUT, index=False)
    print(out_df.to_string(index=False))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
