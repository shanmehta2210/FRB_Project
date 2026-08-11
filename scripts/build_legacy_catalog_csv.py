"""
Build Legacy Survey DR10 Tractor null catalog (v1).

Output uses ``tractor_mag_r`` (22.5 - 2.5*log10(flux_r) nanomaggies) and ``expAB_r`` = b/a
from shape_e1/e2. Compare to SDSS ``rmag`` (model r), not Petrosian.

Default: random sample over joint Legacy∩SDSS Dec footprint via TAP + ORDER BY RANDOM().
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyvo

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (
    JOINT_DEC_MAX,
    JOINT_DEC_MIN,
    MAG_LIMIT,
    Q0,
    footprint_summary,
    prepare_null_sample,
)

TAP_URL = "https://datalab.noirlab.edu/tap"
TABLE = "ls_dr10.tractor"
DEFAULT_OUT = "catalog/LS_catalog_v1_allsky_modelmr.csv"
MIN_STRICT_POOL = 10_000


def flux_to_mag(flux_nmgy: np.ndarray) -> np.ndarray:
    flux = np.asarray(flux_nmgy, dtype=float)
    out = np.full(flux.shape, np.nan, dtype=float)
    good = np.isfinite(flux) & (flux > 0)
    out[good] = 22.5 - 2.5 * np.log10(flux[good])
    return out


def sigma_from_ivar(ivar: np.ndarray) -> np.ndarray:
    iv = np.asarray(ivar, dtype=float)
    out = np.full(iv.shape, np.nan, dtype=float)
    good = np.isfinite(iv) & (iv > 0)
    out[good] = 1.0 / np.sqrt(iv[good])
    return out


def q_from_e1e2(e1: np.ndarray, e2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    e1a = np.asarray(e1, dtype=float)
    e2a = np.asarray(e2, dtype=float)
    eabs = np.hypot(e1a, e2a)
    q = np.full(eabs.shape, np.nan, dtype=float)
    good = np.isfinite(eabs) & (eabs < 1.0)
    q[good] = (1.0 - eabs[good]) / (1.0 + eabs[good])
    return q, eabs


def q_err_from_e1e2_ivar(
    e1: np.ndarray,
    e2: np.ndarray,
    e1_ivar: np.ndarray,
    e2_ivar: np.ndarray,
) -> np.ndarray:
    e1a = np.asarray(e1, dtype=float)
    e2a = np.asarray(e2, dtype=float)
    s1 = sigma_from_ivar(e1_ivar)
    s2 = sigma_from_ivar(e2_ivar)

    eabs = np.hypot(e1a, e2a)
    out = np.full(eabs.shape, np.nan, dtype=float)

    valid = np.isfinite(eabs) & (eabs < 1.0) & np.isfinite(s1) & np.isfinite(s2)
    if not np.any(valid):
        return out

    idx = np.where(valid)[0]
    ev = eabs[idx]

    sigma_e = np.zeros_like(ev)
    nonzero = ev > 1e-12

    if np.any(nonzero):
        j = idx[nonzero]
        ev_nz = eabs[j]
        term1 = (e1a[j] / ev_nz) ** 2 * (s1[j] ** 2)
        term2 = (e2a[j] / ev_nz) ** 2 * (s2[j] ** 2)
        sigma_e[nonzero] = np.sqrt(term1 + term2)

    if np.any(~nonzero):
        j0 = idx[~nonzero]
        sigma_e[~nonzero] = np.sqrt(0.5 * (s1[j0] ** 2 + s2[j0] ** 2))

    dq_de = 2.0 / (1.0 + ev) ** 2
    q_err = dq_de * sigma_e

    out[idx] = q_err
    return out


def incl_from_q(q: np.ndarray, q0: float = Q0) -> np.ndarray:
    qa = np.asarray(q, dtype=float)
    out = np.full(qa.shape, np.nan, dtype=float)
    good = np.isfinite(qa)
    if not np.any(good):
        return out

    out[good & (qa <= q0)] = 90.0

    hi = good & (qa > q0)
    val = (qa[hi] ** 2 - q0**2) / (1.0 - q0**2)
    val = np.clip(val, 0.0, 1.0)
    out[hi] = np.degrees(np.arccos(np.sqrt(val)))
    return out


def _region_clause(region: str) -> str:
    if region == "joint":
        return f"AND dec >= {JOINT_DEC_MIN} AND dec <= {JOINT_DEC_MAX}"
    if region == "all":
        return ""
    if region == "legacy":
        return ""
    raise ValueError(f"Unknown region: {region!r}")


def run_query(
    service: pyvo.dal.TAPService,
    top_n: int,
    retries: int,
    region: str,
    random_order: bool,
) -> pd.DataFrame:
    region_sql = _region_clause(region)
    order_sql = "ORDER BY RANDOM()" if random_order else "ORDER BY objid"

    query = f"""
    SELECT TOP {int(top_n)}
        objid, ra, dec, type, brick_primary,
        flux_g, flux_r, flux_i, flux_z,
        sersic, shape_r,
        shape_e1, shape_e2,
        shape_e1_ivar, shape_e2_ivar
    FROM {TABLE}
    WHERE brick_primary = 1
      AND type <> 'PSF'
      AND flux_r > 0
      AND shape_e1 IS NOT NULL
      AND shape_e2 IS NOT NULL
      AND shape_e1_ivar IS NOT NULL
      AND shape_e2_ivar IS NOT NULL
      AND shape_e1 > -0.999999 AND shape_e1 < 0.999999
      AND shape_e2 > -0.999999 AND shape_e2 < 0.999999
      {region_sql}
    {order_sql}
    """

    last_err = None
    for attempt in range(1, retries + 1):
        try:
            tbl = service.search(query).to_table()
            return tbl.to_pandas()
        except Exception as exc:
            last_err = exc
            if random_order and attempt == 1:
                print(f"[!] TAP query with ORDER BY RANDOM() failed ({exc}); retrying without random order.")
                order_sql = "ORDER BY objid"
                query = query.replace("ORDER BY RANDOM()", "ORDER BY objid")
                random_order = False
            if attempt < retries:
                time.sleep(2.0 * attempt)
    raise RuntimeError(f"Legacy TAP query failed after {retries} attempts: {last_err}")


def build_catalog(df: pd.DataFrame) -> pd.DataFrame:
    e1 = pd.to_numeric(df["shape_e1"], errors="coerce").to_numpy(dtype=float)
    e2 = pd.to_numeric(df["shape_e2"], errors="coerce").to_numpy(dtype=float)
    e1_ivar = pd.to_numeric(df["shape_e1_ivar"], errors="coerce").to_numpy(dtype=float)
    e2_ivar = pd.to_numeric(df["shape_e2_ivar"], errors="coerce").to_numpy(dtype=float)

    q, eabs = q_from_e1e2(e1, e2)
    q_err = q_err_from_e1e2_ivar(e1, e2, e1_ivar, e2_ivar)
    inc = incl_from_q(q)

    gmag = flux_to_mag(pd.to_numeric(df["flux_g"], errors="coerce").to_numpy(dtype=float))
    tractor_mag_r = flux_to_mag(pd.to_numeric(df["flux_r"], errors="coerce").to_numpy(dtype=float))
    imag = flux_to_mag(pd.to_numeric(df["flux_i"], errors="coerce").to_numpy(dtype=float))
    zmag = flux_to_mag(pd.to_numeric(df["flux_z"], errors="coerce").to_numpy(dtype=float))

    out = pd.DataFrame(
        {
            "RA_ICRS": pd.to_numeric(df["ra"], errors="coerce"),
            "DE_ICRS": pd.to_numeric(df["dec"], errors="coerce"),
            "tractor_mag_r": tractor_mag_r,
            "rmag": tractor_mag_r,
            "expAB_r": q,
            "b_a": q,
            "b_a_err": q_err,
            "q_lt_q0": q <= Q0,
            "inclination_deg_q0_0p2": inc,
            "rPrad": pd.to_numeric(df["shape_r"], errors="coerce"),
            "rdVrad": pd.to_numeric(df["sersic"], errors="coerce"),
            "rdVell": 1.0 - q,
            "gmag": gmag,
            "imag": imag,
            "zmag": zmag,
            "tractor_objid": pd.to_numeric(df["objid"], errors="coerce"),
            "tractor_type": df["type"].astype(str),
            "brick_primary": pd.to_numeric(df["brick_primary"], errors="coerce"),
            "shape_e1": e1,
            "shape_e2": e2,
            "shape_e1_ivar": e1_ivar,
            "shape_e2_ivar": e2_ivar,
            "ellipticity_abs": eabs,
            "flux_r_nmgy": pd.to_numeric(df["flux_r"], errors="coerce"),
        }
    )

    out = out[np.isfinite(out["expAB_r"]) & (out["expAB_r"] >= 0.0) & (out["expAB_r"] <= 1.0)].copy()
    out.reset_index(drop=True, inplace=True)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build Legacy DR10 Tractor null catalog (v1, model mr, joint footprint)."
    )
    parser.add_argument("--top", type=int, default=500_000, help="TAP TOP N rows to fetch.")
    parser.add_argument("--retries", type=int, default=3, help="Retry count for TAP query.")
    parser.add_argument(
        "--region",
        choices=("joint", "all", "legacy"),
        default="joint",
        help="joint = Dec cut for Legacy∩SDSS overlap (default).",
    )
    parser.add_argument(
        "--no-random-order",
        action="store_true",
        help="Disable ORDER BY RANDOM(); shuffle after fetch with --seed instead.",
    )
    parser.add_argument("--seed", type=int, default=42, help="RNG seed for post-query shuffle.")
    parser.add_argument("--out", default=DEFAULT_OUT, help="Output CSV path.")
    parser.add_argument(
        "--min-strict-pool",
        type=int,
        default=MIN_STRICT_POOL,
        help="Minimum rows after strict null cuts (mag + q>q0 + no REX).",
    )
    args = parser.parse_args()

    if args.top <= 0:
        raise ValueError("--top must be > 0")

    svc = pyvo.dal.TAPService(TAP_URL)
    print(f"Querying {TABLE} region={args.region} TOP {args.top}...")
    raw = run_query(
        svc,
        top_n=args.top,
        retries=args.retries,
        region=args.region,
        random_order=not args.no_random_order,
    )
    print(f"Fetched rows: {len(raw)}")

    if args.no_random_order and len(raw) > 0:
        raw = raw.sample(frac=1.0, random_state=args.seed).reset_index(drop=True)
        print(f"Shuffled rows with seed={args.seed}")

    catalog = build_catalog(raw)
    footprint_summary(catalog, "Legacy v1 full")

    strict = prepare_null_sample(
        catalog,
        sample_mode="strict",
        mag_column="tractor_mag_r",
        is_legacy=True,
    )
    footprint_summary(strict, "Legacy v1 strict pool")

    if len(strict) < args.min_strict_pool:
        raise SystemExit(
            f"[!] Strict null pool has {len(strict)} rows < {args.min_strict_pool}. "
            "Increase --top or check footprint cuts."
        )

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    catalog.to_csv(args.out, index=False)
    print(f"Wrote: {args.out} ({len(catalog)} rows, strict pool={len(strict)})")


if __name__ == "__main__":
    main()
