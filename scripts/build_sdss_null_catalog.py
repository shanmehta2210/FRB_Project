"""
Build SDSS DR16 null catalog (v1) for inclination CDF comparisons.

Magnitude: ``rmag`` = ``cmodelMag_r`` (composite, kept for compatibility). Null CDF cuts
use ``modelMag_r`` and ``expAB_r`` after dropping deV profile winners
(``lnLExp_r > lnLDeV_r``).

Footprint: joint Legacy∩SDSS Dec range (dec -30 to +90 deg), primary detections,
type=GALAXY. Random sampling via chunked SQL + in-memory shuffle.

For full SDSS footprint and unbiased HTM random sampling, use
``build_sdss_null_catalog_v2.py`` (writes ``SDSS_catalog_v2_fullsky_modelmr.csv``).
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from astroquery.sdss import SDSS

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (
    JOINT_DEC_MAX,
    JOINT_DEC_MIN,
    MAG_LIMIT,
    Q0,
    assign_sdss_profile_winner_columns,
    footprint_summary,
    prepare_null_sample,
)

DEFAULT_OUT = "SDSS_catalog_v1_allsky_modelmr.csv"
MIN_STRICT_POOL = 10_000


def _sql_chunk(ra_min: float, ra_max: float, top_n: int) -> str:
    return f"""
    SELECT TOP {int(top_n)}
        p.ra AS ra,
        p.dec AS dec,
        p.cmodelMag_r AS cmodelMag_r,
        p.petroMag_r AS petroMag_r,
        p.modelMag_r AS modelMag_r,
        p.modelMag_u AS modelMag_u,
        p.modelMag_g AS modelMag_g,
        p.deVMag_r AS deVMag_r,
        p.expMag_r AS expMag_r,
        p.lnLDeV_r AS lnLDeV_r,
        p.lnLExp_r AS lnLExp_r,
        p.deVAB_r AS deVAB_r,
        p.expAB_r AS expAB_r,
        p.fracDeV_r AS fracDeV_r,
        p.expRad_r AS expRad_r,
        p.deVRad_r AS deVRad_r,
        p.type AS type
    FROM PhotoObj AS p
    WHERE p.ra >= {ra_min} AND p.ra < {ra_max}
      AND p.dec >= {JOINT_DEC_MIN} AND p.dec <= {JOINT_DEC_MAX}
      AND p.type = 3
      AND p.clean = 1
      AND p.cmodelMag_r > 0 AND p.cmodelMag_r < 90
      AND p.deVAB_r > 0 AND p.deVAB_r <= 1
      AND p.expAB_r > 0 AND p.expAB_r <= 1
    """


def query_sdss_chunked(
    total_target: int,
    chunk_size: int,
    n_ra_bins: int,
    retries: int,
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    ra_edges = np.linspace(0.0, 360.0, n_ra_bins + 1)
    per_bin = max(chunk_size, total_target // n_ra_bins + 1)

    for i in range(n_ra_bins):
        ra_min = float(ra_edges[i])
        ra_max = float(ra_edges[i + 1])
        sql = _sql_chunk(ra_min, ra_max, per_bin)
        last_err = None
        for attempt in range(1, retries + 1):
            try:
                tbl = SDSS.query_sql(sql, timeout=600)
                df = tbl.to_pandas()
                if len(df):
                    frames.append(df)
                    print(f"  RA bin [{ra_min:.1f}, {ra_max:.1f}): {len(df)} rows")
                break
            except Exception as exc:
                last_err = exc
                if attempt < retries:
                    time.sleep(2.0 * attempt)
        else:
            print(f"[!] RA bin [{ra_min:.1f}, {ra_max:.1f}) failed: {last_err}")

    if not frames:
        raise RuntimeError("SDSS query returned no rows from any RA chunk.")

    out = pd.concat(frames, ignore_index=True)
    out = out.drop_duplicates(subset=["ra", "dec", "cmodelMag_r"], keep="first")
    return out


def build_catalog(df: pd.DataFrame) -> pd.DataFrame:
    rmag = pd.to_numeric(df["cmodelMag_r"], errors="coerce")
    expab = pd.to_numeric(df["expAB_r"], errors="coerce")
    devab = pd.to_numeric(df["deVAB_r"], errors="coerce")
    model_mag = pd.to_numeric(df["modelMag_r"], errors="coerce")
    model_mag_u = pd.to_numeric(df.get("modelMag_u"), errors="coerce")
    model_mag_g = pd.to_numeric(df.get("modelMag_g"), errors="coerce")
    dev_mag = pd.to_numeric(df["deVMag_r"], errors="coerce")
    exp_mag = pd.to_numeric(df["expMag_r"], errors="coerce")
    frac_dev = pd.to_numeric(df["fracDeV_r"], errors="coerce")

    d_dev = np.abs(dev_mag - model_mag)
    d_exp = np.abs(exp_mag - model_mag)
    use_dev = d_dev <= d_exp
    best_model_ba = np.where(use_dev, devab, expab)
    exp_re = pd.to_numeric(df.get("expRad_r"), errors="coerce")
    dev_re = pd.to_numeric(df.get("deVRad_r"), errors="coerce")
    best_model_re = np.where(use_dev, dev_re, exp_re)
    n_eff = 1.0 + 3.0 * frac_dev

    ln_dev = pd.to_numeric(df.get("lnLDeV_r"), errors="coerce")
    ln_exp = pd.to_numeric(df.get("lnLExp_r"), errors="coerce")

    out = pd.DataFrame(
        {
            "RA_ICRS": pd.to_numeric(df["ra"], errors="coerce"),
            "DE_ICRS": pd.to_numeric(df["dec"], errors="coerce"),
            "rmag": rmag,
            "petroMag_r": pd.to_numeric(df.get("petroMag_r"), errors="coerce"),
            "modelMag_r": model_mag,
            "modelMag_u": model_mag_u,
            "modelMag_g": model_mag_g,
            "deVMag_r": dev_mag,
            "expMag_r": exp_mag,
            "lnLDeV_r": ln_dev,
            "lnLExp_r": ln_exp,
            "u_r": model_mag_u - model_mag,
            "g_r": model_mag_g - model_mag,
            "expAB_r": expab,
            "deVAB_r": devab,
            "best_model_ba_r": best_model_ba,
            "fracDeV_r": frac_dev,
            "expRad_r": exp_re,
            "deVRad_r": dev_re,
            "best_model_re_r": best_model_re,
            "n_eff_r": n_eff,
            "b_a": expab,
            "sdss_type": pd.to_numeric(df.get("type"), errors="coerce"),
        }
    )

    good = (
        np.isfinite(out["RA_ICRS"])
        & np.isfinite(out["DE_ICRS"])
        & np.isfinite(out["rmag"])
        & np.isfinite(out["expAB_r"])
        & np.isfinite(out["best_model_ba_r"])
        & (out["expAB_r"] >= 0.0)
        & (out["expAB_r"] <= 1.0)
        & (out["best_model_ba_r"] >= 0.0)
        & (out["best_model_ba_r"] <= 1.0)
    )
    out = out.loc[good].copy().reset_index(drop=True)
    if out["lnLExp_r"].notna().any():
        out = assign_sdss_profile_winner_columns(out)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Build SDSS DR16 null catalog (v1).")
    parser.add_argument("--top", type=int, default=500_000, help="Target rows to collect (before dedupe).")
    parser.add_argument("--chunk-size", type=int, default=50_000, help="TOP N per RA bin query.")
    parser.add_argument("--ra-bins", type=int, default=12, help="Number of RA chunks (0-360 deg).")
    parser.add_argument("--retries", type=int, default=3, help="Retries per chunk.")
    parser.add_argument("--seed", type=int, default=42, help="Shuffle seed after fetch.")
    parser.add_argument("--out", default=DEFAULT_OUT, help="Output CSV path.")
    parser.add_argument(
        "--min-strict-pool",
        type=int,
        default=MIN_STRICT_POOL,
        help="Minimum rows after strict null cuts.",
    )
    args = parser.parse_args()

    print(
        f"Querying SDSS PhotoObj (type=GALAXY, clean=1, "
        f"Dec [{JOINT_DEC_MIN}, {JOINT_DEC_MAX}])..."
    )
    raw = query_sdss_chunked(
        total_target=args.top,
        chunk_size=args.chunk_size,
        n_ra_bins=args.ra_bins,
        retries=args.retries,
    )
    print(f"Fetched {len(raw)} rows (deduped)")

    if len(raw) > args.top:
        raw = raw.sample(n=args.top, random_state=args.seed).reset_index(drop=True)
    else:
        raw = raw.sample(frac=1.0, random_state=args.seed).reset_index(drop=True)

    catalog = build_catalog(raw)
    footprint_summary(catalog, "SDSS v1 full")

    strict = prepare_null_sample(catalog, sample_mode="strict", mag_column="rmag", is_legacy=False)
    footprint_summary(strict, "SDSS v1 strict pool")

    if len(strict) < args.min_strict_pool:
        raise SystemExit(
            f"[!] Strict null pool has {len(strict)} rows < {args.min_strict_pool}. "
            "Increase --top or --ra-bins."
        )

    catalog.to_csv(args.out, index=False)
    print(f"Wrote: {args.out} ({len(catalog)} rows, strict pool={len(strict)})")


if __name__ == "__main__":
    main()
