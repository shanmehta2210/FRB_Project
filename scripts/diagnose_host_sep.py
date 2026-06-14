#!/usr/bin/env python3
"""Nearest SExtractor sources to a localization position."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

_SCRIPTS = Path(__file__).resolve().parent
_REPO = _SCRIPTS.parent
if str(_SCRIPTS / ".." / "pipeline_scripts" / "galfit_fitting") not in sys.path:
    pass

# LDAC reader lives in pipeline
sys.path.insert(0, str(_REPO / "pipeline_scripts" / "galfit_fitting"))
from generate_galfit_cutouts import get_table_from_ldac  # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("frb")
    p.add_argument("--ra", type=float, default=None)
    p.add_argument("--dec", type=float, default=None)
    p.add_argument("--radius-arcsec", type=float, default=30.0)
    args = p.parse_args()

    out = _REPO / "pipeline_scripts" / "Output" / f"{args.frb}_all"
    cat_path = out / "image.cat"
    flux_path = _REPO / "large_cutouts" / f"{args.frb}_flux.fits"
    if not cat_path.is_file():
        raise SystemExit(f"No catalog: {cat_path}")

    loc = __import__("pandas").read_csv(_REPO / "master_frb_localization.csv")
    row = loc.loc[loc["frb"] == args.frb].iloc[0]
    ra = args.ra if args.ra is not None else float(row["ra_deg"])
    dec = args.dec if args.dec is not None else float(row["dec_deg"])
    semantics = row.get("coord_semantics", "")

    cat = get_table_from_ldac(str(cat_path))
    host = SkyCoord(ra=ra, dec=dec, unit="deg")
    c = SkyCoord(ra=cat["ALPHAWIN_J2000"], dec=cat["DELTAWIN_J2000"], unit="deg")
    sep = host.separation(c).arcsec

    order = np.argsort(sep)
    print(f"FRB {args.frb}  host RA/Dec={ra:.6f} {dec:.6f}  semantics={semantics}")
    print(f"Catalog: {cat_path}  N={len(cat)}")
    print(f"\nNearest 10 SExtractor sources:")
    print(f"{'#':>6} {'sep\"':>7} {'MAG_AUTO':>9} {'CLASS_STAR':>11} {'FLUX_RADIUS':>12}")
    for i in order[:10]:
        print(
            f"{int(cat['NUMBER'][i]):6d} {sep[i]:7.2f} "
            f"{float(cat['MAG_AUTO'][i]):9.2f} {float(cat['CLASS_STAR'][i]):11.3f} "
            f"{float(cat['FLUX_RADIUS'][i]):12.2f}"
        )

    within = sep <= args.radius_arcsec
    print(f"\nWithin {args.radius_arcsec}\" : {int(within.sum())} sources")

    if flux_path.is_file():
        with fits.open(flux_path) as hdul:
            w = WCS(hdul[0].header).celestial
        x, y = w.world_to_pixel(host)
        print(f"Host pixel (flux WCS): x={float(x):.1f} y={float(y):.1f}  shape={np.squeeze(hdul[0].data).shape}")


if __name__ == "__main__":
    main()
