"""
Shared constants and helpers for COSMOS HST vs SDSS b/a audit.

See plots/plots_null/v2/sdss_audit/COSMOS/CATALOG_DECISIONS.md.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from null_catalog_utils import Q0, SDSS_UR_MAX_CDF

REPO_ROOT = Path(__file__).resolve().parent.parent
COSMOS_ROOT = REPO_ROOT / "plots" / "plots_null" / "v2" / "sdss_audit" / "COSMOS"
COSMOS_DATA = COSMOS_ROOT / "data"
COSMOS_PLOTS = COSMOS_ROOT / "plots"

# ACS/WFC F814W contiguous mosaic (Scoville et al. 2007)
COSMOS_RA_MIN = 149.43
COSMOS_RA_MAX = 150.80
COSMOS_DEC_MIN = 1.57
COSMOS_DEC_MAX = 2.89

COSMOS_UR_MAX = SDSS_UR_MAX_CDF  # 2.3 — symmetric u-r cut
COSMOS_BA_MIN = Q0  # 0.2 strict b/a

ZURICH_TBL = COSMOS_DATA / "cosmos_morph_zurich_1.0.tbl"
GIM2D_DISK_TYPE = 2  # Zurich TYPE: disk/spiral

HST_CSV = COSMOS_ROOT / "cosmos_hst_zurich_strict.csv"
HST_ENTIRE_CSV = COSMOS_ROOT / "cosmos_hst_zurich_entire.csv"
HST_DISK_CSV = COSMOS_ROOT / "cosmos_hst_zurich_disk_strict.csv"
HST_DISK_ENTIRE_CSV = COSMOS_ROOT / "cosmos_hst_zurich_disk_entire.csv"

SDSS_CSV = COSMOS_ROOT / "cosmos_sdss_dr17_nocolor_strict.csv"
SDSS_ENTIRE_CSV = COSMOS_ROOT / "cosmos_sdss_dr17_nocolor_entire.csv"
SDSS_DISK_CSV = COSMOS_ROOT / "cosmos_sdss_dr17_nocolor_disk_strict.csv"
SDSS_DISK_ENTIRE_CSV = COSMOS_ROOT / "cosmos_sdss_dr17_nocolor_disk_entire.csv"

COSMOS_CDFS = COSMOS_PLOTS / "cdfs"

# Zurich GIM2D reliable depth (Sargent+2007)
HST_MAG_MAX_RELIABLE = 22.5

MAG_STEP = 0.5
MAG_LO_START = 15.0
MIN_N_BIN = 25

HST_MAG_COL = "ACS_MAG_AUTO"
HST_BA_COL = "b_a"
HST_RE_COL = "Re_arcsec"

SDSS_MAG_COL = "modelMag_r"
SDSS_BA_COL = "expAB_r"
SDSS_RE_COL = "expRad_r"


def apply_ba_strict(df: pd.DataFrame, ba_col: str, q0: float = COSMOS_BA_MIN) -> pd.DataFrame:
    ba = pd.to_numeric(df[ba_col], errors="coerce")
    return df.loc[ba > q0].copy()


def in_footprint_mask(ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    return (
        np.isfinite(ra)
        & np.isfinite(dec)
        & (ra >= COSMOS_RA_MIN)
        & (ra <= COSMOS_RA_MAX)
        & (dec >= COSMOS_DEC_MIN)
        & (dec <= COSMOS_DEC_MAX)
    )


def mag_bin_table(
    mag: np.ndarray,
    ba: np.ndarray,
    *,
    pool_n: int,
    mag_col_label: str,
    step: float = MAG_STEP,
    min_n: int = MIN_N_BIN,
) -> pd.DataFrame:
    rows: list[dict] = []
    lo = MAG_LO_START
    while lo < 28.0:
        hi = lo + step
        if lo == MAG_LO_START:
            mask = np.isfinite(mag) & (mag <= hi)
        else:
            mask = np.isfinite(mag) & (mag > lo) & (mag <= hi)
        n = int(mask.sum())
        if n < min_n:
            lo = hi
            continue
        rows.append(
            {
                "mag_lo": lo,
                "mag_hi": hi,
                "n": n,
                "frac_pool_pct": 100.0 * n / max(1, pool_n),
                "median_b_a": float(np.median(ba[mask])),
                "mean_b_a": float(np.mean(ba[mask])),
            }
        )
        lo = hi
    out = pd.DataFrame(rows)
    out.attrs["mag_column"] = mag_col_label
    return out


def spearman_mag_ba(mag: np.ndarray, ba: np.ndarray) -> tuple[float, float]:
    from scipy import stats

    ok = np.isfinite(mag) & np.isfinite(ba)
    if ok.sum() < 10:
        return float("nan"), float("nan")
    rho, p = stats.spearmanr(mag[ok], ba[ok])
    return float(rho), float(p)


def cut_funnel_row(stage: str, n: int, note: str = "") -> dict:
    return {"stage": stage, "n_remaining": n, "note": note}
