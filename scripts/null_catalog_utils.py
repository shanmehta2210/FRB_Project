"""
Shared helpers for null-catalog construction and CDF plotting (v1).

Magnitude policy: compare SDSS ``rmag`` (model r) to Legacy ``tractor_mag_r``
(22.5 - 2.5*log10(flux_r) nanomaggies). Do not use Legacy ``petroMag_r`` (v0 misname).

Sample modes:
  strict    — require q > q0 before building cos(i) pools (mode A)
  inclusive — finite q in [0, 1]; q <= q0 maps to cos(i)=0 / i=90 deg (mode B)
"""

from __future__ import annotations

import math
from typing import Iterable

import numpy as np
import pandas as pd

Q0 = 0.2
MAG_LIMIT = 21.0

# Approximate Legacy DR10 ∩ SDSS DR16 imaging overlap (conservative Dec cut).
JOINT_DEC_MIN = -30.0
JOINT_DEC_MAX = 90.0


def hubble_cosi_from_ba(b_over_a: float, q0: float = Q0) -> float:
    """cos(i) from axis ratio using Hubble formula with intrinsic thickness q0."""
    val = (float(b_over_a) ** 2 - q0**2) / (1.0 - q0**2)
    if not math.isfinite(val) or val < 0:
        return 0.0
    if val > 1:
        return 1.0
    return math.sqrt(val)


def resolve_mag_column(df: pd.DataFrame, mag_column: str) -> str:
    """Pick magnitude column; Legacy may use tractor_mag_r or rmag."""
    if mag_column in df.columns:
        return mag_column
    if mag_column == "rmag" and "tractor_mag_r" in df.columns:
        return "tractor_mag_r"
    if mag_column == "tractor_mag_r" and "rmag" in df.columns:
        return "rmag"
    raise KeyError(f"Magnitude column {mag_column!r} not in dataframe: {list(df.columns)}")


def apply_mag_cut(
    df: pd.DataFrame,
    mag_column: str = "rmag",
    limit: float = MAG_LIMIT,
) -> pd.DataFrame:
    col = resolve_mag_column(df, mag_column)
    mag = pd.to_numeric(df[col], errors="coerce")
    return df.loc[mag <= limit].copy()


def apply_strict_q_cut(
    df: pd.DataFrame,
    q_col: str = "expAB_r",
    q0: float = Q0,
) -> pd.DataFrame:
    q = pd.to_numeric(df[q_col], errors="coerce")
    return df.loc[q > q0].copy()


def apply_inclusive_q(df: pd.DataFrame, q_col: str = "expAB_r") -> pd.DataFrame:
    q = pd.to_numeric(df[q_col], errors="coerce")
    return df.loc[np.isfinite(q) & (q >= 0.0) & (q <= 1.0)].copy()


def apply_legacy_type_cut(
    df: pd.DataFrame,
    exclude: str | Iterable[str] = "REX",
    type_col: str = "tractor_type",
) -> pd.DataFrame:
    if type_col not in df.columns:
        return df.copy()
    if isinstance(exclude, str):
        exclude_set = {t.strip().upper() for t in exclude.split(",") if t.strip()}
    else:
        exclude_set = {str(t).strip().upper() for t in exclude}
    if not exclude_set:
        return df.copy()
    types = df[type_col].astype(str).str.upper()
    return df.loc[~types.isin(exclude_set)].copy()


def resolve_q_column(df: pd.DataFrame, q_column: str) -> str:
    """Resolve axis-ratio column (SDSS may use best_model_ba_r vs expAB_r)."""
    if q_column in df.columns:
        return q_column
    if q_column == "best_model_ba" and "best_model_ba_r" in df.columns:
        return "best_model_ba_r"
    if q_column == "expAB_r":
        return "expAB_r"
    raise KeyError(f"Axis-ratio column {q_column!r} not in dataframe: {list(df.columns)}")


def prepare_null_sample(
    df: pd.DataFrame,
    *,
    sample_mode: str = "strict",
    mag_column: str = "rmag",
    mag_limit: float = MAG_LIMIT,
    q0: float = Q0,
    q_column: str = "expAB_r",
    exclude_legacy_types: str = "REX",
    is_legacy: bool = False,
) -> pd.DataFrame:
    """Apply standard null-catalog cuts for CDF construction."""
    q_col = resolve_q_column(df, q_column)
    out = apply_mag_cut(df, mag_column=mag_column, limit=mag_limit)
    if is_legacy:
        out = apply_legacy_type_cut(out, exclude=exclude_legacy_types)
    if sample_mode == "strict":
        out = apply_strict_q_cut(out, q_col=q_col, q0=q0)
    elif sample_mode == "inclusive":
        out = apply_inclusive_q(out, q_col=q_col)
    else:
        raise ValueError(f"Unknown sample_mode: {sample_mode!r}")
    return out.reset_index(drop=True)


def cosi_array_from_df(
    df: pd.DataFrame,
    q_col: str = "expAB_r",
    q0: float = Q0,
) -> np.ndarray:
    col = resolve_q_column(df, q_col)
    qvals = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
    return np.array([hubble_cosi_from_ba(v, q0=q0) for v in qvals], dtype=float)


def footprint_summary(df: pd.DataFrame, label: str = "") -> dict:
    ra = pd.to_numeric(df.get("RA_ICRS"), errors="coerce")
    dec = pd.to_numeric(df.get("DE_ICRS"), errors="coerce")
    prefix = f"{label}: " if label else ""
    summary = {
        "label": label,
        "n_rows": int(len(df)),
        "ra_min": float(ra.min()) if ra.notna().any() else float("nan"),
        "ra_max": float(ra.max()) if ra.notna().any() else float("nan"),
        "dec_min": float(dec.min()) if dec.notna().any() else float("nan"),
        "dec_max": float(dec.max()) if dec.notna().any() else float("nan"),
    }
    print(
        f"{prefix}n={summary['n_rows']}  "
        f"RA=[{summary['ra_min']:.3f}, {summary['ra_max']:.3f}]  "
        f"Dec=[{summary['dec_min']:.3f}, {summary['dec_max']:.3f}]"
    )
    return summary


def joint_footprint_sql_legacy() -> str:
    return f"dec >= {JOINT_DEC_MIN} AND dec <= {JOINT_DEC_MAX}"


def joint_footprint_sql_sdss() -> str:
    return f"dec >= {JOINT_DEC_MIN} AND dec <= {JOINT_DEC_MAX}"
