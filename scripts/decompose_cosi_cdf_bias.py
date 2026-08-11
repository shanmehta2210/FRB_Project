#!/usr/bin/env python3
"""Decompose SDSS null cos(i) CDF bias vs uniform — per cut stage."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (  # noqa: E402
    Q0,
    SDSS_UR_MAX_CDF,
    apply_mag_cut,
    apply_strict_q_cut,
    cosi_array_from_df,
    filter_sdss_drop_dev_winners,
    filter_sdss_ur,
    hubble_cosi_from_ba,
    read_sdss_null_catalog,
    sdss_exp_wins_lnl_mask,
)

MAG = 20.0
OUT = Path(__file__).resolve().parents[1] / "cosi_cdf_audit"


def empirical_cdf(values: np.ndarray, x: np.ndarray) -> np.ndarray:
    v = np.sort(values[np.isfinite(values)])
    return np.searchsorted(v, x, side="right") / max(len(v), 1)


def median_cosi(values: np.ndarray) -> float:
    v = values[np.isfinite(values)]
    return float(np.median(v)) if len(v) else float("nan")


def cosi_from_isotropic(n: int, seed: int = 42) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.uniform(0.0, 1.0, size=n)


def cosi_isotropic_then_strict_ba(n: int, q0: float = Q0, seed: int = 42) -> np.ndarray:
    """Isotropic cos(i), then keep only Hubble-implied b/a > q0."""
    cosi = cosi_from_isotropic(n, seed)
    ba = np.array([math_ba(c, q0) for c in cosi])
    keep = ba > q0
    return cosi[keep]


def math_ba(cosi: float, q0: float) -> float:
    return float(np.sqrt(cosi**2 * (1 - q0**2) + q0**2))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    x = np.linspace(0, 1, 101)

    # --- Theory / Monte Carlo ---
    iso = cosi_from_isotropic(200_000)
    iso_ba = cosi_isotropic_then_strict_ba(200_000)
    rows_mc = [
        {
            "stage": "MC isotropic cos(i)",
            "n": len(iso),
            "median_cosi": median_cosi(iso),
            "cdf_at_cosi_0.5": float(empirical_cdf(iso, np.array([0.5]))[0]),
        },
        {
            "stage": "MC isotropic + strict b/a>0.2 via Hubble",
            "n": len(iso_ba),
            "median_cosi": median_cosi(iso_ba),
            "cdf_at_cosi_0.5": float(empirical_cdf(iso_ba, np.array([0.5]))[0]),
        },
    ]

    # --- SDSS catalog stages ---
    df = read_sdss_null_catalog(
        Path(__file__).resolve().parents[1] / "catalog/SDSS_catalog_v1_allsky_modelmr.csv"
    )
    stages: list[tuple[str, pd.DataFrame]] = [("catalog_all", df)]
    s = apply_mag_cut(df, mag_column="modelMag_r", limit=MAG)
    stages.append((f"mag_lt_{MAG}", s))
    s = filter_sdss_ur(s, SDSS_UR_MAX_CDF)
    stages.append((f"ur_lt_{SDSS_UR_MAX_CDF}", s))
    s = filter_sdss_drop_dev_winners(s)
    stages.append(("lnL_exp_wins", s))
    s = apply_strict_q_cut(s, q_col="expAB_r", q0=Q0)
    stages.append((f"strict_ba_gt_{Q0}", s))

    rows_sdss = []
    cdfs = {}
    for name, sub in stages:
        cosi = cosi_array_from_df(sub, q_col="expAB_r", q0=Q0)
        rows_sdss.append(
            {
                "stage": name,
                "n": len(sub),
                "median_cosi": median_cosi(cosi),
                "mean_cosi": float(np.mean(cosi)),
                "cdf_at_cosi_0.5": float(empirical_cdf(cosi, np.array([0.5]))[0]),
                "frac_cosi_gt_0.5": float(np.mean(cosi > 0.5)),
                "median_expAB_r": float(
                    pd.to_numeric(sub["expAB_r"], errors="coerce").median()
                ),
            }
        )
        cdfs[name] = empirical_cdf(cosi, x)

    tab = pd.DataFrame(rows_mc + rows_sdss)
    tab.to_csv(OUT / "stage_summary_mag20.csv", index=False)

    # CDF curves for plotting
    cdf_out = pd.DataFrame({"cos_i": x, "uniform": x})
    for name, curve in cdfs.items():
        cdf_out[name] = curve
    cdf_out.to_csv(OUT / "cdf_curves_mag20.csv", index=False)

    print(tab.to_string(index=False))
    print(f"\nWrote {OUT}")


if __name__ == "__main__":
    main()
