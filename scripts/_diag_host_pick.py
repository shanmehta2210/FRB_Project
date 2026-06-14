"""Diagnose SExtractor nearest-host pick vs CSV coords."""
import sys
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table

REPO = Path(__file__).resolve().parents[1]
frbs = [
    ("20180301A", 93.227, 4.671),
    ("20230718A", 128.16191666666666, -40.45194444444445),
    ("20180924B", 326.105235868384, -40.9002526146074),
]
for frb, ra, dec in frbs:
    cat_path = REPO / "pipeline_scripts" / "Output" / f"{frb}_all" / "image.cat"
    if not cat_path.is_file():
        print(frb, "no catalog")
        continue
    cat = Table.read(cat_path, hdu=2)
    host = SkyCoord(ra, dec, unit="deg")
    c = SkyCoord(ra=cat["ALPHAWIN_J2000"], dec=cat["DELTAWIN_J2000"], unit="deg")
    sep = host.separation(c).arcsec
    order = np.argsort(sep)
    print(f"\n=== {frb}  n_src={len(cat)} ===")
    for j in order[:8]:
        print(
            f"  #{int(cat['NUMBER'][j]):4d}  sep={sep[j]:6.2f}\"  "
            f"MAG_AUTO={cat['MAG_AUTO'][j]:6.2f}  FLUX_AUTO={cat['FLUX_AUTO'][j]:12.3g}  "
            f"CLASS_STAR={cat['CLASS_STAR'][j]:.2f}"
        )
    j0 = order[0]
    print(f"  -> Phase3a picks #{int(cat['NUMBER'][j0])} at sep={sep[j0]:.2f}\"")

    flux_p = REPO / "large_cutouts" / f"{frb}_flux.fits"
    if flux_p.is_file():
        from astropy.io import fits

        with fits.open(flux_p) as h:
            d = h[0].data
            iv = fits.getdata(REPO / "large_cutouts" / f"{frb}_invvar.fits")
        print(
            f"  cutout flux_med={np.nanmedian(d):.4g}  invvar>0={np.mean(iv > 0):.3f}  "
            f"inv_med={np.nanmedian(iv[iv > 0]) if np.any(iv > 0) else 0:.3g}"
        )
