"""
Shared helpers for null-catalog construction and CDF plotting (v1).

Magnitude policy: SDSS null CDFs use ``modelMag_r`` and ``expAB_r`` (exponential
profile only; deV winners dropped via ``lnLExp_r > lnLDeV_r``). Legacy uses
``tractor_mag_r`` (alias ``rmag``) and ``expAB_r``
from Tractor shape_e1/e2. ``rmag`` on SDSS CSV is ``cmodelMag_r`` (composite), kept
for backward compatibility only. Do not use Legacy ``petroMag_r`` (v0 misname).

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
# SDSS PhotoObj stores expAB_r/deVAB_r with a catalog floor at 0.05 (~22% of v1 sample).
# Excluding b/a <= this value removes the quantization "rung" (see diagnostics/).
SDSS_BA_FLOOR_MIN = 0.05

# Approximate Legacy DR10 ∩ SDSS DR16 imaging overlap (conservative Dec cut).
JOINT_DEC_MIN = -30.0
JOINT_DEC_MAX = 90.0

# SDSS v2 full imaging footprint (no joint Dec clip); HTM random sampling.
SDSS_CATALOG_V2_DEFAULT = "SDSS_catalog_v2_fullsky_modelmr.csv"
SDSS_HTM_PRIME = 37
SDSS_HTM_MASK = 0x000000000000FFFF
SDSS_MIN_STRICT_MAG20_POOL = 50_000
SDSS_MAG20_LIMIT = 20.0
SDSS_LNL_COVERAGE_MIN = 0.95

# SDSS u-r color cuts for diagnostics (Strateva et al. 2001 uses u*-r* ~ 2.22).
UR_CUTS_DEFAULT = (3.5, 2.2, 1.5)
# Production strict null CDFs (late-type / spiral-like pools).
SDSS_UR_MAX_CDF = 2.3
LEGACY_GR_MAX_CDF = 0.75
LEGACY_MORPH_N_MIN = 0.75
LEGACY_MORPH_N_MAX = 2.0
LEGACY_CDF_TYPE_EXCLUDE = "REX,DEV"
SDSS_Q_COLUMN_CDF = "expAB_r"
COLOR_FIT_SUBSAMPLE = 50_000
COLOR_FIT_SEED = 42

# Minimal columns for RAM-safe null CDF work (one read per catalog).
LEGACY_NULL_USECOLS = (
    "tractor_mag_r",
    "rmag",
    "gmag",
    "expAB_r",
    "tractor_type",
    "rdVrad",
)
SDSS_NULL_USECOLS = (
    "modelMag_r",
    "modelMag_u",
    "modelMag_g",
    "expAB_r",
    "u_r",
    "g_r",
    "lnLDeV_r",
    "lnLExp_r",
    "model_winner_is_exp",
)
# Extra columns when pipeline diagnostics need Re/n CDFs.
LEGACY_NULL_USECOLS_EXTENDED = LEGACY_NULL_USECOLS + ("rPrad", "rdVrad")
SDSS_NULL_USECOLS_EXTENDED = SDSS_NULL_USECOLS + (
    "best_model_re_r",
    "fracDeV_r",
    "n_eff_r",
    "deVMag_r",
    "expMag_r",
    "deVAB_r",
    "best_model_ba_r",
)


def hubble_cosi_from_ba(b_over_a: float, q0: float = Q0) -> float:
    """cos(i) from axis ratio using Hubble formula with intrinsic thickness q0."""
    val = (float(b_over_a) ** 2 - q0**2) / (1.0 - q0**2)
    if not math.isfinite(val) or val < 0:
        return 0.0
    if val > 1:
        return 1.0
    return math.sqrt(val)


def resolve_mag_column(df: pd.DataFrame, mag_column: str) -> str:
    """Pick magnitude column; Legacy may use tractor_mag_r or rmag; SDSS modelMag_r."""
    if mag_column in df.columns:
        return mag_column
    if mag_column in ("best_model_mag_r", "best_model_mag") and "modelMag_r" in df.columns:
        return "modelMag_r"
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


def filter_frb_hosts_strict_ba(
    hosts: pd.DataFrame,
    *,
    q0: float = Q0,
    ba_col: str = "b_a",
) -> pd.DataFrame:
    """Keep FRB hosts with GALFIT axis ratio strictly above q0 (matches null strict pool)."""
    if ba_col not in hosts.columns:
        raise KeyError(f"FRB hosts missing {ba_col!r} for strict cut")
    ba = pd.to_numeric(hosts[ba_col], errors="coerce")
    out = hosts.loc[ba > q0].copy()
    return out.reset_index(drop=True)


def filter_frb_hosts_mag(
    hosts: pd.DataFrame,
    *,
    mag_limit: float,
    mag_column: str = "mag",
) -> pd.DataFrame:
    """Keep FRB hosts with GALFIT magnitude <= mag_limit (pipeline ``mag`` column)."""
    return apply_mag_cut(hosts, mag_column=mag_column, limit=mag_limit).reset_index(drop=True)


def apply_inclusive_q(df: pd.DataFrame, q_col: str = "expAB_r") -> pd.DataFrame:
    q = pd.to_numeric(df[q_col], errors="coerce")
    return df.loc[np.isfinite(q) & (q >= 0.0) & (q <= 1.0)].copy()


def apply_sdss_ba_floor_cut(
    df: pd.DataFrame,
    q_col: str = "best_model_ba_r",
    min_ba: float = SDSS_BA_FLOOR_MIN,
) -> pd.DataFrame:
    """Drop SDSS rows at the PhotoObj axis-ratio storage floor (default b/a > 0.05)."""
    q = pd.to_numeric(df[q_col], errors="coerce")
    return df.loc[q > min_ba].copy()


def sdss_exp_wins_lnl_mask(df: pd.DataFrame) -> pd.Series:
    """
    True where r-band exponential profile strictly beats deVaucouleurs in lnL.

    Requires finite ``lnLExp_r`` and ``lnLDeV_r``; rows without lnL are False.
    """
    if "lnLExp_r" not in df.columns or "lnLDeV_r" not in df.columns:
        raise KeyError(
            "SDSS catalog missing lnLExp_r / lnLDeV_r; run scripts/patch_sdss_profile_winner.py"
        )
    ln_exp = pd.to_numeric(df["lnLExp_r"], errors="coerce")
    ln_dev = pd.to_numeric(df["lnLDeV_r"], errors="coerce")
    return ln_exp.notna() & ln_dev.notna() & (ln_exp > ln_dev)


def assign_sdss_profile_winner_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Set ``model_winner_is_exp`` (1/0) from ``sdss_exp_wins_lnl_mask``."""
    out = df.copy()
    wins = sdss_exp_wins_lnl_mask(out)
    out["model_winner_is_exp"] = wins.astype(np.int8)
    return out


def mag_proxy_sdss_exp_winner(df: pd.DataFrame) -> pd.Series:
    """Audit helper: exp wins if |expMag - modelMag| < |deVMag - modelMag|."""
    model_mag = pd.to_numeric(df["modelMag_r"], errors="coerce")
    dev_mag = pd.to_numeric(df["deVMag_r"], errors="coerce")
    exp_mag = pd.to_numeric(df["expMag_r"], errors="coerce")
    d_dev = (dev_mag - model_mag).abs()
    d_exp = (exp_mag - model_mag).abs()
    return d_exp < d_dev


def filter_sdss_drop_dev_winners(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep only galaxies with finite lnL and ``lnLExp_r > lnLDeV_r`` (strict).

    Never uses ``best_model_ba_r`` or ``model_winner_is_exp`` without re-checking lnL.
    """
    keep = sdss_exp_wins_lnl_mask(df)
    out = df.loc[keep].copy()
    if len(out):
        ln_exp = pd.to_numeric(out["lnLExp_r"], errors="coerce")
        ln_dev = pd.to_numeric(out["lnLDeV_r"], errors="coerce")
        if not (ln_exp > ln_dev).all():
            raise RuntimeError("SDSS pool contains deV-winning or tied lnL profiles")
    return out


def legacy_spiral_morph_mask(df: pd.DataFrame) -> pd.Series:
    """True for EXP type or Sérsic n in [LEGACY_MORPH_N_MIN, LEGACY_MORPH_N_MAX]."""
    if "tractor_type" not in df.columns:
        raise KeyError("Legacy catalog missing tractor_type")
    types = df["tractor_type"].astype(str).str.upper()
    if "rdVrad" not in df.columns:
        raise KeyError("Legacy catalog missing rdVrad (Sérsic index)")
    n = pd.to_numeric(df["rdVrad"], errors="coerce")
    is_exp = types == "EXP"
    in_n = n.notna() & (n >= LEGACY_MORPH_N_MIN) & (n <= LEGACY_MORPH_N_MAX)
    return is_exp | in_n


def filter_legacy_spiral_morph(df: pd.DataFrame) -> pd.DataFrame:
    """Late-type morphology: Tractor EXP or disk-like Sérsic index."""
    return df.loc[legacy_spiral_morph_mask(df)].copy()


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
    exclude_sdss_ba_floor: bool = False,
    sdss_min_ba: float = SDSS_BA_FLOOR_MIN,
    sdss_ur_max: float | None = None,
    legacy_gr_max: float | None = None,
    sdss_exp_winner_only: bool = True,
    legacy_spiral_morph_only: bool = True,
) -> pd.DataFrame:
    """Apply standard null-catalog cuts for CDF construction."""
    q_col = resolve_q_column(df, q_column)
    out = apply_mag_cut(df, mag_column=mag_column, limit=mag_limit)
    if not is_legacy and exclude_sdss_ba_floor:
        out = apply_sdss_ba_floor_cut(out, q_col=q_col, min_ba=sdss_min_ba)
    if is_legacy:
        out = apply_legacy_type_cut(out, exclude=exclude_legacy_types)
    else:
        if sdss_ur_max is not None:
            out = filter_sdss_ur(out, sdss_ur_max)
        if sdss_exp_winner_only:
            out = filter_sdss_drop_dev_winners(out)
            if "expAB_r" not in out.columns:
                raise KeyError("SDSS CDF pools must use expAB_r after lnL exp-winner cut")
    if is_legacy and legacy_gr_max is not None:
        out = filter_legacy_gr(out, legacy_gr_max)
    if is_legacy and legacy_spiral_morph_only:
        out = filter_legacy_spiral_morph(out)
    if sample_mode == "strict":
        out = apply_strict_q_cut(out, q_col=q_col, q0=q0)
    elif sample_mode == "inclusive":
        out = apply_inclusive_q(out, q_col=q_col)
    else:
        raise ValueError(f"Unknown sample_mode: {sample_mode!r}")
    return out.reset_index(drop=True)


def prepare_null_strict_color_base(
    df: pd.DataFrame,
    *,
    mag_column: str,
    q0: float = Q0,
    q_column: str,
    exclude_legacy_types: str = "REX",
    is_legacy: bool = False,
    exclude_sdss_ba_floor: bool = False,
    sdss_min_ba: float = SDSS_BA_FLOOR_MIN,
    sdss_ur_max: float | None = None,
    legacy_gr_max: float | None = None,
    sdss_exp_winner_only: bool = True,
    legacy_spiral_morph_only: bool = True,
) -> pd.DataFrame:
    """
    Strict + color (+ morphology) cuts without magnitude limit (for mag-sliced CDF loops).

    Returns a trimmed frame: magnitude column, axis-ratio column, and any
    extra columns already present (e.g. rPrad for Re CDFs).
    """
    q_col = resolve_q_column(df, q_column)
    out = df
    if not is_legacy and exclude_sdss_ba_floor:
        out = apply_sdss_ba_floor_cut(out, q_col=q_col, min_ba=sdss_min_ba)
    if is_legacy:
        out = apply_legacy_type_cut(out, exclude=exclude_legacy_types)
    else:
        if sdss_ur_max is not None:
            out = filter_sdss_ur(out, sdss_ur_max)
        if sdss_exp_winner_only:
            out = filter_sdss_drop_dev_winners(out)
            if "expAB_r" not in out.columns:
                raise KeyError("SDSS CDF pools must use expAB_r after lnL exp-winner cut")
    if is_legacy and legacy_gr_max is not None:
        out = filter_legacy_gr(out, legacy_gr_max)
    if is_legacy and legacy_spiral_morph_only:
        out = filter_legacy_spiral_morph(out)
    out = apply_strict_q_cut(out, q_col=q_col, q0=q0)
    mag_col = resolve_mag_column(out, mag_column)
    keep = {mag_col, q_col}
    if is_legacy:
        keep.update(c for c in ("rPrad", "rdVrad", "tractor_type") if c in out.columns)
    else:
        keep.update(
            c
            for c in (
                "best_model_re_r",
                "n_eff_r",
                "fracDeV_r",
                "lnLDeV_r",
                "lnLExp_r",
                "model_winner_is_exp",
            )
            if c in out.columns
        )
    return out.loc[:, sorted(keep)].reset_index(drop=True)


def prepare_null_inclusive_color_base(
    df: pd.DataFrame,
    *,
    mag_column: str,
    q0: float = Q0,
    q_column: str,
    exclude_legacy_types: str = "REX",
    is_legacy: bool = False,
    exclude_sdss_ba_floor: bool = False,
    sdss_min_ba: float = SDSS_BA_FLOOR_MIN,
    sdss_ur_max: float | None = None,
    legacy_gr_max: float | None = None,
    sdss_exp_winner_only: bool = True,
    legacy_spiral_morph_only: bool = False,
) -> pd.DataFrame:
    """
    Inclusive b/a + color cuts without magnitude limit (for mag-sliced tests).

    Finite b/a in [0, 1]; no strict b/a > q0 requirement.
    """
    q_col = resolve_q_column(df, q_column)
    out = df
    if not is_legacy and exclude_sdss_ba_floor:
        out = apply_sdss_ba_floor_cut(out, q_col=q_col, min_ba=sdss_min_ba)
    if is_legacy:
        out = apply_legacy_type_cut(out, exclude=exclude_legacy_types)
    else:
        if sdss_ur_max is not None:
            out = filter_sdss_ur(out, sdss_ur_max)
        if sdss_exp_winner_only:
            out = filter_sdss_drop_dev_winners(out)
            if "expAB_r" not in out.columns:
                raise KeyError("SDSS CDF pools must use expAB_r after lnL exp-winner cut")
    if is_legacy and legacy_gr_max is not None:
        out = filter_legacy_gr(out, legacy_gr_max)
    if is_legacy and legacy_spiral_morph_only:
        out = filter_legacy_spiral_morph(out)
    out = apply_inclusive_q(out, q_col=q_col)
    mag_col = resolve_mag_column(out, mag_column)
    keep = {mag_col, q_col}
    if is_legacy:
        keep.update(c for c in ("rPrad", "rdVrad") if c in out.columns)
    else:
        keep.update(
            c
            for c in ("best_model_re_r", "n_eff_r", "fracDeV_r")
            if c in out.columns
        )
    return out.loc[:, sorted(keep)].reset_index(drop=True)


def slice_null_base_by_mag(
    base: pd.DataFrame,
    *,
    mag_column: str,
    mag_limit: float,
) -> pd.DataFrame:
    """Apply magnitude cut to a strict+color base pool (view-friendly copy)."""
    col = resolve_mag_column(base, mag_column)
    mag = pd.to_numeric(base[col], errors="coerce")
    return base.loc[mag <= mag_limit].copy()


def cosi_array_from_df(
    df: pd.DataFrame,
    q_col: str = "expAB_r",
    q0: float = Q0,
) -> np.ndarray:
    col = resolve_q_column(df, q_col)
    qvals = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
    return np.array([hubble_cosi_from_ba(v, q0=q0) for v in qvals], dtype=float)


def inc_deg_from_cosi(cosi: np.ndarray) -> np.ndarray:
    """Inclination in degrees from cos(i), clipped to [0, 1]."""
    c = np.asarray(cosi, dtype=float)
    return np.degrees(np.arccos(np.clip(c, 0.0, 1.0)))


def equal_count_quantile_edges(values: np.ndarray, n_bins: int = 8) -> np.ndarray:
    """
    Bin edges splitting ``values`` into ``n_bins`` pools of (nearly) equal size.

    Uses ``pandas.qcut``; duplicate edges are dropped if the distribution is discrete.
    """
    v = pd.Series(np.asarray(values, dtype=float))
    v = v[np.isfinite(v)]
    if len(v) < n_bins:
        raise ValueError(f"Need at least {n_bins} finite values for {n_bins} bins; got {len(v)}")
    _, edges = pd.qcut(v, q=n_bins, retbins=True, duplicates="drop")
    out = np.asarray(edges, dtype=float)
    if len(out) < 2:
        raise ValueError("qcut produced no usable bin edges")
    return out


def assign_values_to_bin_edges(values: np.ndarray, edges: np.ndarray) -> np.ndarray:
    """
    Integer bin index per value (0 .. n_bins-1) for half-open bins except the last (closed).

    Returns -1 for non-finite values.
    """
    v = np.asarray(values, dtype=float)
    idx = np.full(v.shape, -1, dtype=int)
    finite = np.isfinite(v)
    if not finite.any():
        return idx
    edges = np.asarray(edges, dtype=float)
    n_bins = len(edges) - 1
    for b in range(n_bins):
        lo, hi = edges[b], edges[b + 1]
        if b < n_bins - 1:
            mask = finite & (v >= lo) & (v < hi)
        else:
            mask = finite & (v >= lo) & (v <= hi)
        idx[mask] = b
    return idx


def _finite_positive_array(df: pd.DataFrame, col: str) -> np.ndarray:
    if col not in df.columns:
        raise KeyError(f"Column {col!r} not in dataframe")
    vals = pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)
    return vals[np.isfinite(vals) & (vals > 0)]


def re_arcsec_from_legacy_df(df: pd.DataFrame, col: str = "rPrad") -> np.ndarray:
    """Legacy Tractor effective radius (arcsec)."""
    return _finite_positive_array(df, col)


def n_from_legacy_df(df: pd.DataFrame, col: str = "rdVrad") -> np.ndarray:
    """Legacy Tractor Sérsic index."""
    return _finite_positive_array(df, col)


def re_arcsec_from_sdss_df(df: pd.DataFrame, col: str = "best_model_re_r") -> np.ndarray:
    """SDSS best-model effective radius (arcsec)."""
    return _finite_positive_array(df, col)


def n_from_sdss_df(df: pd.DataFrame, col: str = "n_eff_r") -> np.ndarray:
    """SDSS effective Sérsic index proxy from fracDeV."""
    return _finite_positive_array(df, col)


def ensure_sdss_n_eff(df: pd.DataFrame) -> pd.DataFrame:
    """Add ``n_eff_r`` from ``fracDeV_r`` if missing."""
    out = df.copy()
    if "n_eff_r" not in out.columns or out["n_eff_r"].isna().all():
        frac = pd.to_numeric(out.get("fracDeV_r"), errors="coerce")
        out["n_eff_r"] = 1.0 + 3.0 * frac
    return out


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


def sdss_footprint_sql_v2() -> str:
    """No Dec clip — SDSS-imaged area is defined by PhotoObj flags in the SQL query."""
    return ""


def sdss_htm_hash_sql(alias: str = "p", prime: int = SDSS_HTM_PRIME) -> str:
    """SDSS SkyServer HTM hash expression for unbiased random subsampling."""
    return f"({alias}.htmid * {int(prime)} & {SDSS_HTM_MASK:#018x})"


def sdss_htm_stratum_edges(n_strata: int) -> list[tuple[int, int]]:
    """
    Partition the 16-bit HTM hash into ``n_strata`` contiguous slices.

    Returns list of (lo, hi) inclusive bounds on the hash integer.
    """
    if n_strata < 1:
        raise ValueError(f"n_strata must be >= 1, got {n_strata}")
    n_bins = SDSS_HTM_MASK + 1  # 65536
    edges = np.linspace(0, n_bins, n_strata + 1, dtype=int)
    strata: list[tuple[int, int]] = []
    for i in range(n_strata):
        lo = int(edges[i])
        hi = int(edges[i + 1]) - 1
        if hi < lo:
            hi = lo
        strata.append((lo, hi))
    strata[-1] = (strata[-1][0], n_bins - 1)
    return strata


def sdss_htm_stratum_clause(
    stratum_lo: int,
    stratum_hi: int,
    *,
    alias: str = "p",
    prime: int = SDSS_HTM_PRIME,
) -> str:
    """SQL WHERE fragment restricting to one HTM hash stratum."""
    h = sdss_htm_hash_sql(alias=alias, prime=prime)
    return f"{h} >= {int(stratum_lo)} AND {h} <= {int(stratum_hi)}"


def count_strict_mag20_pool(
    df: pd.DataFrame,
    *,
    mag_limit: float = SDSS_MAG20_LIMIT,
) -> int:
    """Production strict null pool size at ``modelMag_r <= mag_limit``."""
    base = prepare_null_strict_color_base(
        df,
        mag_column="modelMag_r",
        q0=Q0,
        q_column=SDSS_Q_COLUMN_CDF,
        is_legacy=False,
        sdss_ur_max=SDSS_UR_MAX_CDF,
        sdss_exp_winner_only=True,
    )
    cut = slice_null_base_by_mag(base, mag_column="modelMag_r", mag_limit=mag_limit)
    return len(cut)


def cut_funnel_rows(
    df: pd.DataFrame,
    *,
    survey: str,
    mag_limit: float,
    sample_mode: str,
    mag_column: str,
    q_column: str,
    q0: float = Q0,
    exclude_legacy_types: str = "REX",
    is_legacy: bool = False,
    exclude_sdss_ba_floor: bool = False,
    sdss_min_ba: float = SDSS_BA_FLOOR_MIN,
    sdss_ur_max: float | None = None,
    legacy_gr_max: float | None = None,
    sdss_exp_winner_only: bool = True,
    legacy_spiral_morph_only: bool = True,
    legacy_cdf_type_exclude: str = LEGACY_CDF_TYPE_EXCLUDE,
) -> list[dict]:
    """Row counts after each null-catalog cut stage (for audit tables)."""
    rows: list[dict] = []
    n0 = len(df)

    def _row(stage: str, n: int, dropped: int) -> dict:
        return {
            "survey": survey,
            "mag_limit": mag_limit,
            "sample_mode": sample_mode,
            "stage": stage,
            "n_remaining": n,
            "n_dropped": dropped,
            "mag_column": mag_column,
            "q_column": q_column,
            "exclude_sdss_ba_floor": exclude_sdss_ba_floor,
        }

    rows.append(_row("catalog_raw", n0, 0))

    after_mag = apply_mag_cut(df, mag_column=mag_column, limit=mag_limit)
    rows.append(_row("after_mag_cut", len(after_mag), n0 - len(after_mag)))

    pool = after_mag
    if not is_legacy and exclude_sdss_ba_floor:
        q_col_pre = resolve_q_column(pool, q_column)
        pool = apply_sdss_ba_floor_cut(pool, q_col=q_col_pre, min_ba=sdss_min_ba)
        rows.append(
            _row(
                "after_sdss_ba_floor_exclude",
                len(pool),
                len(after_mag) - len(pool),
            )
        )

    if is_legacy:
        excl = legacy_cdf_type_exclude if legacy_spiral_morph_only else exclude_legacy_types
        after_type = apply_legacy_type_cut(pool, exclude=excl)
        rows.append(
            _row("after_type_exclude", len(after_type), len(pool) - len(after_type))
        )
        pool = after_type
    else:
        if sdss_ur_max is not None:
            before = len(pool)
            pool = filter_sdss_ur(pool, sdss_ur_max)
            rows.append(_row("after_sdss_ur", len(pool), before - len(pool)))
        if sdss_exp_winner_only:
            before = len(pool)
            pool = filter_sdss_drop_dev_winners(pool)
            rows.append(
                _row("after_sdss_exp_lnl_winner", len(pool), before - len(pool))
            )

    if is_legacy and legacy_gr_max is not None:
        before = len(pool)
        pool = filter_legacy_gr(pool, legacy_gr_max)
        rows.append(_row("after_legacy_gr", len(pool), before - len(pool)))
    if is_legacy and legacy_spiral_morph_only:
        before = len(pool)
        pool = filter_legacy_spiral_morph(pool)
        rows.append(_row("after_legacy_spiral_morph", len(pool), before - len(pool)))

    q_col = resolve_q_column(pool, q_column)
    before_q = len(pool)
    if sample_mode == "strict":
        pool = apply_strict_q_cut(pool, q_col=q_col, q0=q0)
        rows.append(_row("after_strict_ba", len(pool), before_q - len(pool)))
    elif sample_mode == "inclusive":
        pool = apply_inclusive_q(pool, q_col=q_col)
        rows.append(_row("after_inclusive_q", len(pool), before_q - len(pool)))
    else:
        raise ValueError(f"Unknown sample_mode: {sample_mode!r}")

    final = pool
    rows.append(_row("final_pool", len(final), 0))
    return rows


def _read_csv_usecols(path: Path, usecols: tuple[str, ...]) -> pd.DataFrame:
    """Read only columns present in the CSV (skip missing optional cols)."""
    header = pd.read_csv(path, nrows=0).columns.tolist()
    cols = [c for c in usecols if c in header]
    missing = set(usecols) - set(cols)
    if missing:
        print(f"[!] {path.name}: skipping missing columns {sorted(missing)}")
    return pd.read_csv(path, usecols=cols)


def read_legacy_null_catalog(
    path: Path,
    *,
    extended: bool = False,
) -> pd.DataFrame:
    usecols = LEGACY_NULL_USECOLS_EXTENDED if extended else LEGACY_NULL_USECOLS
    return _read_csv_usecols(path, usecols)


def read_sdss_null_catalog(
    path: Path,
    *,
    extended: bool = False,
) -> pd.DataFrame:
    usecols = SDSS_NULL_USECOLS_EXTENDED if extended else SDSS_NULL_USECOLS
    df = _read_csv_usecols(path, usecols)
    return ensure_sdss_colors(df)


def ensure_sdss_colors(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure u_r and g_r columns exist from modelMag_u/g and modelMag_r."""
    out = df.copy()
    mr = pd.to_numeric(out["modelMag_r"], errors="coerce")
    if "modelMag_u" in out.columns:
        mu = pd.to_numeric(out["modelMag_u"], errors="coerce")
        out["u_r"] = mu - mr
    if "modelMag_g" in out.columns:
        mg = pd.to_numeric(out["modelMag_g"], errors="coerce")
        out["g_r"] = mg - mr
    return out


def fit_ur_vs_gr(
    sdss: pd.DataFrame,
    *,
    subsample: int = COLOR_FIT_SUBSAMPLE,
    seed: int = COLOR_FIT_SEED,
) -> tuple[float, float]:
    """
    Robust linear fit u_r = alpha + beta * g_r on SDSS galaxies.

    WARNING: expects an in-memory DataFrame. Do not pass the full ~500k-row v1 CSV
    unless you have spare RAM. For production color work use plot_sdss_color_cuts.py.

    Returns (alpha, beta). Uses Theil-Sen if scipy is available, else polyfit.
    """
    sdss = ensure_sdss_colors(sdss)
    ur = pd.to_numeric(sdss["u_r"], errors="coerce")
    gr = pd.to_numeric(sdss["g_r"], errors="coerce")
    ok = ur.notna() & gr.notna() & np.isfinite(ur) & np.isfinite(gr)
    if ok.sum() < 1000:
        raise RuntimeError(
            f"Too few finite u_r/g_r for color fit ({ok.sum()}); run augment_sdss_v1_colors.py"
        )
    ur_v = ur.loc[ok].to_numpy(dtype=float)
    gr_v = gr.loc[ok].to_numpy(dtype=float)
    if len(ur_v) > subsample:
        rng = np.random.default_rng(seed)
        idx = rng.choice(len(ur_v), size=subsample, replace=False)
        ur_v = ur_v[idx]
        gr_v = gr_v[idx]
    try:
        from scipy.stats import theilslopes

        beta, alpha = theilslopes(ur_v, gr_v)[:2]
    except ImportError:
        beta, alpha = np.polyfit(gr_v, ur_v, 1)
    return float(alpha), float(beta)


def legacy_gr_limits_from_ur_cuts(
    ur_cuts: Iterable[float],
    *,
    alpha: float,
    beta: float,
) -> dict[float, float]:
    """Map SDSS u-r maxima to Legacy g-r proxy cuts: g_r_max = (u_r_max - alpha) / beta."""
    if abs(beta) < 1e-6:
        raise ValueError(f"Unstable color fit: beta={beta}")
    return {float(ur): (float(ur) - alpha) / beta for ur in ur_cuts}


def filter_sdss_ur(sdss: pd.DataFrame, ur_max: float) -> pd.DataFrame:
    sdss = ensure_sdss_colors(sdss)
    ur = pd.to_numeric(sdss["u_r"], errors="coerce")
    return sdss.loc[ur < ur_max].copy()


def filter_legacy_gr(legacy: pd.DataFrame, gr_max: float) -> pd.DataFrame:
    gmag = pd.to_numeric(legacy["gmag"], errors="coerce")
    rmag = pd.to_numeric(legacy["rmag"], errors="coerce")
    gr = gmag - rmag
    return legacy.loc[gr < gr_max].copy()


def ur_cut_folder_tag(ur_max: float) -> str:
    """e.g. 3.5 -> ur_lt_3p5"""
    s = f"{ur_max:g}".replace(".", "p")
    return f"ur_lt_{s}"


def gr_cut_folder_tag(gr_max: float) -> str:
    """e.g. 0.75 -> gr_lt_0p75"""
    s = f"{gr_max:g}".replace(".", "p")
    return f"gr_lt_{s}"
