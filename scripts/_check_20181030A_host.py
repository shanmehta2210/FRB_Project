from astropy.table import Table
from astropy.coordinates import SkyCoord

cat = Table.read(
    "CHIME/pipeline_scripts/Output/20181030A_all/.workdir/image.cat", hdu=2, format="fits"
)
psf = Table.read(
    "CHIME/pipeline_scripts/Output/20181030A_all/.workdir/image.psf.cat", hdu=2, format="fits"
)

loc = SkyCoord(158.594, 73.763767, unit="deg")
ap = SkyCoord(158.59505369837754, 73.76390142054908, unit="deg")

for num in [105]:
    for name, t in [("image.cat", cat), ("image.psf.cat", psf)]:
        sel = t["NUMBER"] == num
        if not sel.any():
            print(name, num, "NOT FOUND")
            continue
        c = SkyCoord(t["ALPHAWIN_J2000"][sel][0], t["DELTAWIN_J2000"][sel][0], unit="deg")
        print(
            f"{name} #{num}: RA={c.ra.deg:.6f} Dec={c.dec.deg:.6f}  "
            f"sep_loc={c.separation(loc).arcsec:.2f}\"  sep_ap={c.separation(ap).arcsec:.2f}\"  "
            f"X={float(t['X_IMAGE'][sel][0]):.1f} Y={float(t['Y_IMAGE'][sel][0]):.1f}"
        )

coords = SkyCoord(cat["ALPHAWIN_J2000"], cat["DELTAWIN_J2000"], unit="deg")
seps = ap.separation(coords).arcsec
idx = int(seps.argmin())
print(
    f"Nearest image.cat to AstroPath: #{int(cat['NUMBER'][idx])} sep={float(seps[idx]):.3f}\" "
    f"at ({float(cat['X_IMAGE'][idx]):.1f},{float(cat['Y_IMAGE'][idx]):.1f})"
)
