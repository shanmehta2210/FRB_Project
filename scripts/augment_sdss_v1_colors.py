#!/usr/bin/env python3
"""
Add modelMag_u, modelMag_g and derived u_r, g_r to an existing SDSS v1 null CSV.

Matches rows by sky position to SDSS DR16 PhotoObj (type=GALAXY, clean=1)
in the same joint footprint as build_sdss_null_catalog.py.

Run from repo root:
    python scripts/augment_sdss_v1_colors.py
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

DEFAULT_CSV = "catalog/SDSS_catalog_v1_allsky_modelmr.csv"
MATCH_ARCSEC = 1.0


def _sql_colors(ra_min: float, ra_max: float, top_n: int) -> str:
    return f"""
    SELECT TOP {int(top_n)}
        p.ra AS ra,
        p.dec AS dec,
        p.modelMag_r AS modelMag_r,
        p.modelMag_u AS modelMag_u,
        p.modelMag_g AS modelMag_g
    FROM PhotoObj AS p
    WHERE p.ra >= {ra_min} AND p.ra < {ra_max}
      AND p.dec >= {JOINT_DEC_MIN} AND p.dec <= {JOINT_DEC_MAX}
      AND p.type = 3
      AND p.clean = 1
      AND p.modelMag_r > 0 AND p.modelMag_r < 90
      AND p.modelMag_u > 0 AND p.modelMag_u < 90
      AND p.modelMag_g > 0 AND p.modelMag_g < 90
    """


def query_colors_chunked(
    n_ra_bins: int,
    chunk_size: int,
    retries: int,
) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    ra_edges = np.linspace(0.0, 360.0, n_ra_bins + 1)
    per_bin = chunk_size

    for i in range(n_ra_bins):
        ra_min = float(ra_edges[i])
        ra_max = float(ra_edges[i + 1])
        sql = _sql_colors(ra_min, ra_max, per_bin)
        last_err = None
        for attempt in range(1, retries + 1):
            try:
                tbl = SDSS.query_sql(sql, timeout=600)
                df = tbl.to_pandas()
                if len(df):
                    frames.append(df)
                    print(f"  colors RA [{ra_min:.1f}, {ra_max:.1f}): {len(df)} rows")
                break
            except Exception as exc:
                last_err = exc
                if attempt < retries:
                    time.sleep(2.0 * attempt)
        else:
            print(f"[!] colors RA [{ra_min:.1f}, {ra_max:.1f}) failed: {last_err}")

    if not frames:
        raise RuntimeError("SDSS color query returned no rows.")
    out = pd.concat(frames, ignore_index=True)
    out = out.drop_duplicates(subset=["ra", "dec", "modelMag_r"], keep="first")
    return out


def match_colors_to_catalog(
    catalog: pd.DataFrame,
    photo: pd.DataFrame,
    *,
    max_sep_arcsec: float = MATCH_ARCSEC,
) -> pd.DataFrame:
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    ra_c = pd.to_numeric(catalog["RA_ICRS"], errors="coerce")
    dec_c = pd.to_numeric(catalog["DE_ICRS"], errors="coerce")
    ra_p = pd.to_numeric(photo["ra"], errors="coerce")
    dec_p = pd.to_numeric(photo["dec"], errors="coerce")

    ok_c = ra_c.notna() & dec_c.notna()
    ok_p = ra_p.notna() & dec_p.notna()
    cat = catalog.loc[ok_c].copy()
    ph = photo.loc[ok_p].copy()

    c_coord = SkyCoord(ra=ra_c.loc[ok_c].to_numpy(), dec=dec_c.loc[ok_c].to_numpy(), unit="deg")
    p_coord = SkyCoord(ra=ra_p.loc[ok_p].to_numpy(), dec=dec_p.loc[ok_p].to_numpy(), unit="deg")

    idx, sep2d, _ = c_coord.match_to_catalog_sky(p_coord)
    sep_arcsec = sep2d.to(u.arcsec).value

    model_u = np.full(len(cat), np.nan)
    model_g = np.full(len(cat), np.nan)
    good = sep_arcsec <= max_sep_arcsec
    mu = pd.to_numeric(ph["modelMag_u"], errors="coerce").to_numpy()
    mg = pd.to_numeric(ph["modelMag_g"], errors="coerce").to_numpy()
    model_u[good] = mu[idx[good]]
    model_g[good] = mg[idx[good]]

    out = catalog.copy()
    out["modelMag_u"] = np.nan
    out["modelMag_g"] = np.nan
    out.loc[ok_c, "modelMag_u"] = model_u
    out.loc[ok_c, "modelMag_g"] = model_g

    mr = pd.to_numeric(out["modelMag_r"], errors="coerce")
    out["u_r"] = out["modelMag_u"] - mr
    out["g_r"] = out["modelMag_g"] - mr
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=Path(DEFAULT_CSV))
    parser.add_argument("--out", type=Path, default=None, help="Default: overwrite --csv")
    parser.add_argument("--ra-bins", type=int, default=12)
    parser.add_argument("--chunk-size", type=int, default=80_000)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--match-arcsec", type=float, default=MATCH_ARCSEC)
    parser.add_argument("--force", action="store_true", help="Re-query even if columns exist.")
    args = parser.parse_args()

    header = pd.read_csv(args.csv, nrows=0)
    if not args.force and "modelMag_u" in header.columns:
        probe = pd.read_csv(args.csv, usecols=["modelMag_u"], nrows=5000)
        if probe["modelMag_u"].notna().mean() > 0.9:
            print(
                f"[*] {args.csv} already has modelMag_u "
                f"({probe['modelMag_u'].notna().mean():.1%} finite in first 5k rows); skip."
            )
            return

    print(f"[!] Loading full catalog into RAM: {args.csv} (required for sky match merge)")
    catalog = pd.read_csv(args.csv)

    print("Querying SDSS for modelMag_u, modelMag_g...")
    photo = query_colors_chunked(
        n_ra_bins=args.ra_bins,
        chunk_size=args.chunk_size,
        retries=args.retries,
    )
    print(f"PhotoObj color rows: {len(photo):,}")

    merged = match_colors_to_catalog(
        catalog,
        photo,
        max_sep_arcsec=args.match_arcsec,
    )
    n_u = merged["modelMag_u"].notna().sum()
    print(f"Matched modelMag_u: {n_u:,} / {len(merged):,} ({100 * n_u / len(merged):.1f}%)")

    out_path = args.out or args.csv
    merged.to_csv(out_path, index=False)
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
