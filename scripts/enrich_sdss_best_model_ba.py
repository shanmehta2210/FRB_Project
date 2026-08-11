"""
Add SDSS best-model r-band axis ratio to the v1 null catalog.

PhotoObj fits pure deVaucouleurs and pure exponential profiles per band (PSF-convolved;
Stoughton et al. 2002; SDSS DR7 model-magnitude algorithms). ``modelMag_r`` is the
magnitude from the higher-likelihood single profile in r. We set ``best_model_ba_r`` to
``deVAB_r`` or ``expAB_r`` according to which profile's magnitude is closer to
``modelMag_r`` (same spirit as SDSS ``modelMag``).

Also stores ``fracDeV_r`` and ``best_model_ba_fracDeV`` (deV if fracDeV_r >= 0.5 else exp)
for diagnostics. Brightness cut still uses ``rmag`` = ``cmodelMag_r``.
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

from null_catalog_utils import JOINT_DEC_MAX, JOINT_DEC_MIN

DEFAULT_IN = "catalog/SDSS_catalog_v1_allsky_modelmr.csv"
DEFAULT_OUT = "catalog/SDSS_catalog_v1_allsky_modelmr.csv"


def _sql_chunk(ra_min: float, ra_max: float, top_n: int) -> str:
    return f"""
    SELECT TOP {int(top_n)}
        p.ra AS ra,
        p.dec AS dec,
        p.cmodelMag_r AS cmodelMag_r,
        p.modelMag_r AS modelMag_r,
        p.deVMag_r AS deVMag_r,
        p.expMag_r AS expMag_r,
        p.deVAB_r AS deVAB_r,
        p.expAB_r AS expAB_r,
        p.fracDeV_r AS fracDeV_r
    FROM PhotoObj AS p
    WHERE p.ra >= {ra_min} AND p.ra < {ra_max}
      AND p.dec >= {JOINT_DEC_MIN} AND p.dec <= {JOINT_DEC_MAX}
      AND p.type = 3
      AND p.clean = 1
      AND p.cmodelMag_r > 0 AND p.cmodelMag_r < 90
      AND p.deVAB_r > 0 AND p.deVAB_r <= 1
      AND p.expAB_r > 0 AND p.expAB_r <= 1
    """


def query_shape_chunked(n_ra_bins: int, chunk_size: int, retries: int) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    ra_edges = np.linspace(0.0, 360.0, n_ra_bins + 1)
    per_bin = chunk_size

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
            print(f"[!] RA bin failed: {last_err}")

    if not frames:
        raise RuntimeError("No rows returned from shape enrichment query.")
    out = pd.concat(frames, ignore_index=True)
    return out.drop_duplicates(subset=["ra", "dec", "cmodelMag_r"], keep="first")


def assign_best_model_ba(df: pd.DataFrame) -> pd.DataFrame:
    deVab = pd.to_numeric(df["deVAB_r"], errors="coerce")
    expab = pd.to_numeric(df["expAB_r"], errors="coerce")
    model_mag = pd.to_numeric(df["modelMag_r"], errors="coerce")
    dev_mag = pd.to_numeric(df["deVMag_r"], errors="coerce")
    exp_mag = pd.to_numeric(df["expMag_r"], errors="coerce")
    frac = pd.to_numeric(df["fracDeV_r"], errors="coerce")

    d_dev = np.abs(dev_mag - model_mag)
    d_exp = np.abs(exp_mag - model_mag)
    use_dev = d_dev <= d_exp

    best_mag = np.where(use_dev, deVab, expab)
    best_frac = np.where(frac >= 0.5, deVab, expab)

    out = df.copy()
    out["best_model_ba_r"] = best_mag
    out["best_model_ba_fracDeV"] = best_frac
    out["model_winner_is_deV"] = use_dev.astype(int)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Enrich SDSS v1 CSV with best-model b/a.")
    parser.add_argument("--in-csv", default=DEFAULT_IN)
    parser.add_argument("--out-csv", default=DEFAULT_OUT)
    parser.add_argument("--ra-bins", type=int, default=12)
    parser.add_argument("--chunk-size", type=int, default=50_000)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--match-tol-arcsec", type=float, default=1.0)
    args = parser.parse_args()

    print("Querying PhotoObj for modelMag / deVAB / expAB / fracDeV...")
    shape = query_shape_chunked(args.ra_bins, args.chunk_size, args.retries)
    shape = assign_best_model_ba(shape)
    print(f"Shape query rows: {len(shape)}")

    base = pd.read_csv(args.in_csv)
    base["RA_ICRS"] = pd.to_numeric(base["RA_ICRS"], errors="coerce")
    base["DE_ICRS"] = pd.to_numeric(base["DE_ICRS"], errors="coerce")
    base["rmag"] = pd.to_numeric(base["rmag"], errors="coerce")

    shape["ra"] = pd.to_numeric(shape["ra"], errors="coerce")
    shape["dec"] = pd.to_numeric(shape["dec"], errors="coerce")
    shape["cmodelMag_r"] = pd.to_numeric(shape["cmodelMag_r"], errors="coerce")

    merged = base.merge(
        shape[
            [
                "ra",
                "dec",
                "cmodelMag_r",
                "modelMag_r",
                "deVMag_r",
                "expMag_r",
                "deVAB_r",
                "fracDeV_r",
                "best_model_ba_r",
                "best_model_ba_fracDeV",
                "model_winner_is_deV",
            ]
        ],
        left_on=["RA_ICRS", "DE_ICRS", "rmag"],
        right_on=["ra", "dec", "cmodelMag_r"],
        how="left",
    )

    n_match = int(merged["best_model_ba_r"].notna().sum())
    print(f"Matched {n_match} / {len(base)} catalog rows ({n_match / len(base):.1%})")

    if n_match < 0.95 * len(base):
        print("[!] Low match fraction; check coordinates or re-run full SDSS build.")

    # Keep expAB_r as exponential-only; store best-model column separately.
    drop_cols = [c for c in ("ra", "dec", "cmodelMag_r") if c in merged.columns]
    merged = merged.drop(columns=drop_cols, errors="ignore")

    # Diagnostics
    expab = pd.to_numeric(merged["expAB_r"], errors="coerce")
    best = pd.to_numeric(merged["best_model_ba_r"], errors="coerce")
    valid = expab.notna() & best.notna()
    diff = (best - expab)[valid]
    print(f"median(best_model_ba - expAB_r) = {diff.median():.4f}")
    print(f"frac deV winner (modelMag): {merged['model_winner_is_deV'].mean():.3f}")

    merged.to_csv(args.out_csv, index=False)
    print(f"Wrote {args.out_csv}")


if __name__ == "__main__":
    main()
