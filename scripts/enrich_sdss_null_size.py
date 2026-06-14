"""
Enrich SDSS v1 null catalog with effective radius and n proxy columns.

Always (no network):
  - ``n_eff_r`` = 1 + 3 * ``fracDeV_r`` (SDSS bulge fraction proxy for Sérsic n)

With ``--query-radii`` (SDSS SQL):
  - ``expRad_r``, ``deVRad_r`` (arcsec effective radii)
  - ``best_model_re_r`` from profile whose magnitude is closer to ``modelMag_r``
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

from null_catalog_utils import JOINT_DEC_MAX, JOINT_DEC_MIN  # noqa: E402

DEFAULT_CSV = "SDSS_catalog_v1_allsky_modelmr.csv"


def _sql_chunk(ra_min: float, ra_max: float, top_n: int) -> str:
    return f"""
    SELECT TOP {int(top_n)}
        p.ra AS ra,
        p.dec AS dec,
        p.cmodelMag_r AS cmodelMag_r,
        p.modelMag_r AS modelMag_r,
        p.deVMag_r AS deVMag_r,
        p.expMag_r AS expMag_r,
        p.expRad_r AS expRad_r,
        p.deVRad_r AS deVRad_r
    FROM PhotoObj AS p
    WHERE p.ra >= {ra_min} AND p.ra < {ra_max}
      AND p.dec >= {JOINT_DEC_MIN} AND p.dec <= {JOINT_DEC_MAX}
      AND p.type = 3
      AND p.clean = 1
      AND p.cmodelMag_r > 0 AND p.cmodelMag_r < 90
      AND p.expRad_r > 0
      AND p.deVRad_r > 0
    """


def query_radii_chunked(n_ra_bins: int, chunk_size: int, retries: int) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    ra_edges = np.linspace(0.0, 360.0, n_ra_bins + 1)
    for i in range(n_ra_bins):
        ra_min = float(ra_edges[i])
        ra_max = float(ra_edges[i + 1])
        sql = _sql_chunk(ra_min, ra_max, chunk_size)
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
            print(f"[!] RA bin failed: {last_err}")
    if not frames:
        raise RuntimeError("No rows from radius enrichment query.")
    return pd.concat(frames, ignore_index=True).drop_duplicates(
        subset=["ra", "dec", "cmodelMag_r"], keep="first"
    )


def assign_best_model_re(df: pd.DataFrame) -> pd.DataFrame:
    exp_re = pd.to_numeric(df["expRad_r"], errors="coerce")
    dev_re = pd.to_numeric(df["deVRad_r"], errors="coerce")
    model_mag = pd.to_numeric(df["modelMag_r"], errors="coerce")
    dev_mag = pd.to_numeric(df["deVMag_r"], errors="coerce")
    exp_mag = pd.to_numeric(df["expMag_r"], errors="coerce")
    d_dev = np.abs(dev_mag - model_mag)
    d_exp = np.abs(exp_mag - model_mag)
    use_dev = d_dev <= d_exp
    out = df.copy()
    out["best_model_re_r"] = np.where(use_dev, dev_re, exp_re)
    return out


def add_n_eff_r(df: pd.DataFrame) -> pd.DataFrame:
    frac = pd.to_numeric(df["fracDeV_r"], errors="coerce")
    out = df.copy()
    out["n_eff_r"] = 1.0 + 3.0 * frac
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", default=DEFAULT_CSV)
    parser.add_argument("--query-radii", action="store_true", help="Fetch expRad/deVRad via SDSS SQL.")
    parser.add_argument("--ra-bins", type=int, default=12)
    parser.add_argument("--chunk-size", type=int, default=50_000)
    parser.add_argument("--retries", type=int, default=3)
    args = parser.parse_args()

    base = pd.read_csv(args.csv)
    base = add_n_eff_r(base)
    print(f"Added n_eff_r (median={base['n_eff_r'].median():.3f})")

    if args.query_radii:
        print("Querying SDSS for expRad_r / deVRad_r...")
        radii = query_radii_chunked(args.ra_bins, args.chunk_size, args.retries)
        radii = assign_best_model_re(radii)
        base["RA_ICRS"] = pd.to_numeric(base["RA_ICRS"], errors="coerce")
        base["DE_ICRS"] = pd.to_numeric(base["DE_ICRS"], errors="coerce")
        base["rmag"] = pd.to_numeric(base["rmag"], errors="coerce")
        radii["ra"] = pd.to_numeric(radii["ra"], errors="coerce")
        radii["dec"] = pd.to_numeric(radii["dec"], errors="coerce")
        radii["cmodelMag_r"] = pd.to_numeric(radii["cmodelMag_r"], errors="coerce")

        merged = base.drop(
            columns=[c for c in ("expRad_r", "deVRad_r", "best_model_re_r") if c in base.columns],
            errors="ignore",
        ).merge(
            radii[["ra", "dec", "cmodelMag_r", "expRad_r", "deVRad_r", "best_model_re_r"]],
            left_on=["RA_ICRS", "DE_ICRS", "rmag"],
            right_on=["ra", "dec", "cmodelMag_r"],
            how="left",
        )
        drop_cols = [c for c in ("ra", "dec", "cmodelMag_r") if c in merged.columns]
        merged = merged.drop(columns=drop_cols, errors="ignore")
        n_match = int(merged["best_model_re_r"].notna().sum())
        print(f"Matched radii for {n_match} / {len(base)} rows ({n_match / len(base):.1%})")
        base = merged
    elif "best_model_re_r" not in base.columns:
        print("[!] best_model_re_r missing; run with --query-radii for SDSS Re null plots.")

    base.to_csv(args.csv, index=False)
    print(f"Wrote {args.csv}")


if __name__ == "__main__":
    main()
