#!/usr/bin/env python3

"""

Patch SDSS null catalog with r-band profile log-likelihoods and exp-winner flag.



Merges ``lnLDeV_r``, ``lnLExp_r`` from PhotoObj and sets ``model_winner_is_exp``

from ``lnLExp_r > lnLDeV_r`` (strict; ties count as deV).



Run from repo root::



    python scripts/patch_sdss_profile_winner.py

    python scripts/patch_sdss_profile_winner.py --footprint full --in-csv catalog/SDSS_catalog_v2_fullsky_modelmr.csv

"""



from __future__ import annotations



import argparse

import sys

import time

from pathlib import Path



import numpy as np

import pandas as pd



_SCRIPTS = Path(__file__).resolve().parent

if str(_SCRIPTS) not in sys.path:

    sys.path.insert(0, str(_SCRIPTS))



from merge_lnl_patch_into_sdss import merge_patch

from null_catalog_utils import (

    JOINT_DEC_MAX,

    JOINT_DEC_MIN,

    SDSS_CATALOG_V2_DEFAULT,

    assign_sdss_profile_winner_columns,

    sdss_htm_stratum_clause,

    sdss_htm_stratum_edges,

)



DEFAULT_IN = "catalog/SDSS_catalog_v1_allsky_modelmr.csv"

DEFAULT_OUT = "catalog/SDSS_catalog_v1_allsky_modelmr.csv"





def _dec_clause(footprint: str) -> str:

    if footprint == "joint":

        return f"AND p.dec >= {JOINT_DEC_MIN} AND p.dec <= {JOINT_DEC_MAX}"

    if footprint == "full":

        return ""

    raise ValueError(f"Unknown footprint: {footprint!r}")





def _sql_chunk_ra(

    ra_min: float,

    ra_max: float,

    top_n: int,

    *,

    footprint: str,

) -> str:

    dec_sql = _dec_clause(footprint)

    return f"""

    SELECT TOP {int(top_n)}

        p.ra AS ra,

        p.dec AS dec,

        p.cmodelMag_r AS cmodelMag_r,

        p.modelMag_r AS modelMag_r,

        p.deVMag_r AS deVMag_r,

        p.expMag_r AS expMag_r,

        p.lnLDeV_r AS lnLDeV_r,

        p.lnLExp_r AS lnLExp_r

    FROM PhotoObj AS p

    WHERE p.ra >= {ra_min} AND p.ra < {ra_max}

      {dec_sql}

      AND p.type = 3

      AND p.clean = 1

      AND p.mode = 1

      AND p.cmodelMag_r > 0 AND p.cmodelMag_r < 90

    """





def _sql_chunk_htm(

    stratum_lo: int,

    stratum_hi: int,

    top_n: int,

) -> str:

    htm_clause = sdss_htm_stratum_clause(stratum_lo, stratum_hi)

    return f"""

    SELECT TOP {int(top_n)}

        p.ra AS ra,

        p.dec AS dec,

        p.cmodelMag_r AS cmodelMag_r,

        p.modelMag_r AS modelMag_r,

        p.deVMag_r AS deVMag_r,

        p.expMag_r AS expMag_r,

        p.lnLDeV_r AS lnLDeV_r,

        p.lnLExp_r AS lnLExp_r

    FROM PhotoObj AS p

    WHERE p.type = 3

      AND p.clean = 1

      AND p.mode = 1

      AND p.cmodelMag_r > 0 AND p.cmodelMag_r < 90

      AND {htm_clause}

    """





def query_lnl_chunked_ra(

    n_ra_bins: int,

    chunk_size: int,

    retries: int,

    *,

    footprint: str,

) -> pd.DataFrame:

    from astroquery.sdss import SDSS



    frames: list[pd.DataFrame] = []

    ra_edges = np.linspace(0.0, 360.0, n_ra_bins + 1)



    for i in range(n_ra_bins):

        ra_min = float(ra_edges[i])

        ra_max = float(ra_edges[i + 1])

        sql = _sql_chunk_ra(ra_min, ra_max, chunk_size, footprint=footprint)

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

        raise RuntimeError("No rows returned from lnL patch query.")

    out = pd.concat(frames, ignore_index=True)

    return out.drop_duplicates(subset=["ra", "dec", "cmodelMag_r"], keep="first")





def query_lnl_chunked_htm(

    n_strata: int,

    chunk_size: int,

    retries: int,

) -> pd.DataFrame:

    from astroquery.sdss import SDSS



    frames: list[pd.DataFrame] = []

    for idx, (lo, hi) in enumerate(sdss_htm_stratum_edges(n_strata)):

        sql = _sql_chunk_htm(lo, hi, chunk_size)

        last_err = None

        for attempt in range(1, retries + 1):

            try:

                tbl = SDSS.query_sql(sql, timeout=900)

                df = tbl.to_pandas()

                if len(df):

                    frames.append(df)

                    print(f"  HTM stratum {idx} [{lo}, {hi}]: {len(df)} rows")

                break

            except Exception as exc:

                last_err = exc

                if attempt < retries:

                    time.sleep(2.0 * attempt)

        else:

            print(f"[!] HTM stratum {idx} failed: {last_err}")



    if not frames:

        raise RuntimeError("No rows returned from lnL HTM patch query.")

    out = pd.concat(frames, ignore_index=True)

    return out.drop_duplicates(subset=["ra", "dec", "cmodelMag_r"], keep="first")





def main() -> None:

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--in-csv", default=DEFAULT_IN)

    parser.add_argument("--out-csv", default=DEFAULT_OUT)

    parser.add_argument(

        "--footprint",

        choices=("joint", "full"),

        default="joint",

        help="joint: Legacy∩SDSS Dec clip; full: entire SDSS imaging footprint.",

    )

    parser.add_argument("--ra-bins", type=int, default=12)

    parser.add_argument("--n-strata", type=int, default=32)

    parser.add_argument("--chunk-size", type=int, default=80_000)

    parser.add_argument("--retries", type=int, default=3)

    parser.add_argument(

        "--cache-csv",

        default=None,

        help="Save/load SkyServer patch table (default: catalog/SDSS_lnl_patch_cache.csv or _v2).",

    )

    parser.add_argument(

        "--merge-only",

        action="store_true",

        help="Merge from --cache-csv without new SQL queries.",

    )

    args = parser.parse_args()



    if args.cache_csv is None:

        if args.footprint == "full" or args.in_csv == SDSS_CATALOG_V2_DEFAULT:

            args.cache_csv = "catalog/SDSS_lnl_patch_cache_v2.csv"

        else:

            args.cache_csv = "catalog/SDSS_lnl_patch_cache.csv"



    cache_path = Path(args.cache_csv)

    if args.merge_only:

        if not cache_path.is_file():

            raise SystemExit(f"--merge-only requires cache at {cache_path}")

        print(f"Loading patch cache {cache_path}...")

        patch = pd.read_csv(cache_path)

    else:

        if args.footprint == "full":

            print("Querying PhotoObj for lnL (full footprint, HTM strata)...")

            patch = query_lnl_chunked_htm(args.n_strata, args.chunk_size, args.retries)

        else:

            print("Querying PhotoObj for lnL (joint footprint, RA bins)...")

            patch = query_lnl_chunked_ra(

                args.ra_bins, args.chunk_size, args.retries, footprint=args.footprint

            )

        patch.to_csv(cache_path, index=False)

        print(f"Wrote patch cache {cache_path}")

    for c in ("lnLDeV_r", "lnLExp_r"):

        patch[c] = pd.to_numeric(patch[c], errors="coerce")



    base = pd.read_csv(args.in_csv)

    merged = merge_patch(base, patch)



    n_lnl = int(pd.to_numeric(merged["lnLExp_r"], errors="coerce").notna().sum())

    print(f"Matched lnL columns for {n_lnl} / {len(base)} rows ({n_lnl / len(base):.1%})")

    if n_lnl < 0.95 * len(base):

        print("[!] Low match fraction; increase --chunk-size or --n-strata and re-run.")



    merged = assign_sdss_profile_winner_columns(merged)

    finite = pd.to_numeric(merged["lnLExp_r"], errors="coerce").notna()

    frac_exp = float(merged.loc[finite, "model_winner_is_exp"].mean())

    print(f"Fraction lnL exp-winner (finite lnL only): {frac_exp:.3f}")



    merged.to_csv(args.out_csv, index=False)

    print(f"Wrote {args.out_csv}")





if __name__ == "__main__":

    main()


