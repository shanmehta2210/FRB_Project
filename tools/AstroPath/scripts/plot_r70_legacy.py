import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.visualization import ZScaleInterval
import astropy.units as u

cdf = pd.read_csv("tools/AstroPath/results/r70_posterior.csv")
best = cdf.sort_values("posterior_O", ascending=False).iloc[0]
target_center = SkyCoord(ra=64.3996075, dec=7.931106, unit="deg")

fig, ax_img = plt.subplots(1, 1, figsize=(10, 10))
interval = ZScaleInterval()

fits_path = "tools/Photometry/coadded_astrometrically_corrected_rband_r70.fits"
with fits.open(fits_path) as hh:
    w = WCS(hh[0].header).celestial
    data_2d = np.squeeze(hh[0].data)
    while data_2d.ndim > 2: data_2d = data_2d[0]
    px_scale = 0.262

vmin, vmax = interval.get_limits(data_2d)
ax_img.imshow(data_2d, origin='lower', cmap='bone', vmin=vmin, vmax=vmax)

bx, by = w.world_to_pixel(SkyCoord(ra=best['ra'], dec=best['dec'], unit='deg'))
ax_img.scatter(bx, by, s=150, facecolors='none', edgecolors='red', lw=2, label=f"Host (P={best['posterior_O']:.2f})")
ax_img.text(bx, by - 15, f"{best['posterior_O']:.2f}", color='red', fontsize=10, ha='center', va='top', weight='bold')

cx, cy = w.world_to_pixel(target_center)
ax_img.scatter(cx, cy, s=150, marker='+', color='cyan', lw=2, label=f"FRB Center")

px_per_arcmin = 60.0 / px_scale
cx_min, cx_max = cx - 1.0*px_per_arcmin, cx + 1.0*px_per_arcmin
cy_min, cy_max = cy - 1.0*px_per_arcmin, cy + 1.0*px_per_arcmin

# Plot all other candidates
for _, row in cdf.iterrows():
    if row['objid'] != best['objid'] and row['posterior_O'] > 0:
        px, py = w.world_to_pixel(SkyCoord(ra=row['ra'], dec=row['dec'], unit='deg'))
        if cx_min <= px <= cx_max and cy_min <= py <= cy_max:
            ax_img.scatter(px, py, s=50, facecolors='none', edgecolors='orange', alpha=0.6)
            ax_img.text(px, py - 10, f"{row['posterior_O']:.2f}", color='orange', fontsize=8, ha='center', va='top')

ax_img.set_xlim(cx_min, cx_max)
ax_img.set_ylim(cy_min, cy_max)
ax_img.legend(loc='lower left')
ax_img.set_title("r-70 AstroPath Association (10 arcsec bounds)")

plt.savefig("tools/AstroPath/results/r70_10arcsec_plot.png", bbox_inches='tight', dpi=300)
print("Plot generated at tools/AstroPath/results/r70_10arcsec_plot.png")
