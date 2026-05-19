import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.visualization import ZScaleInterval, AsymmetricPercentileInterval
import matplotlib.patches as patches

# Define paths
IMAGE_PATH = 'tools/Photometry/coadded_astrometrically_corrected_rband_r70.fits'
PHOTOMETRY_PATH = 'tools/Photometry/r70_target_comparison_photometry.csv'
POSTERIOR_PATH = 'tools/AstroPath/results/r70_posterior.csv'
OUTPUT_DIR = 'plots'
OUTPUT_PATH = os.path.join(OUTPUT_DIR, 'r70_posterior_overlay.png')

# Create output directory if it doesn't exist
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 1. Load Data
with fits.open(IMAGE_PATH) as hdul:
    data = hdul[0].data
    header = hdul[0].header
    wcs = WCS(header)

phot_df = pd.read_csv(PHOTOMETRY_PATH)
post_df = pd.read_csv(POSTERIOR_PATH)

# Clean column names just in case
phot_df.columns = phot_df.columns.str.strip()
post_df.columns = post_df.columns.str.strip()

# 2. Merge Catalogs
# Merge on objid_ls (Photometry) == objid (AstroPath)
merged_df = pd.merge(phot_df, post_df, left_on='objid_ls', right_on='objid', how='left')

# 3. Generate Cutout
center_coord = SkyCoord('04h17m35.9058s', '+07d55m51.9812s', frame='icrs')
size = u.Quantity((1, 1), u.arcmin)

cutout = Cutout2D(data, center_coord, size, wcs=wcs)

# 4. Plot Structure
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, projection=cutout.wcs)

# Enhance contrast using percentile interval
interval = AsymmetricPercentileInterval(0.5, 99.5)
vmin, vmax = interval.get_limits(cutout.data)
norm = plt.Normalize(vmin=vmin, vmax=vmax)

ax.imshow(cutout.data, origin='lower', cmap='gray_r', norm=norm)

ax.set_xlabel('Right Ascension')
ax.set_ylabel('Declination')
ax.set_title('r70 Field - 1 Arcmin Cutout (40px Mag & Posterior)')

# Plot a crosshair for the FRB target center
target_x, target_y = cutout.wcs.world_to_pixel(center_coord)
ax.plot(target_x, target_y, 'g+', markersize=15, markeredgewidth=2, label='Target Center')

# 5. Overlay Sources and Annotations
for _, row in merged_df.iterrows():
    # check if the source has valid RA/Dec
    if pd.isna(row['RA']) or pd.isna(row['Dec']):
        continue
    
    coord = SkyCoord(row['RA'], row['Dec'], unit='deg', frame='icrs')
    # Check if the coordinate is within the cutout
    if cutout.wcs.footprint_contains(coord):
        # Convert to pixel coordinates of the cutout
        x, y = cutout.wcs.world_to_pixel(coord)
        
        flux_radius = row.get('FLUX_RADIUS', 5.0)
        if pd.isna(flux_radius):
            flux_radius = 5.0
            
        radius = flux_radius * 2

        # Draw circle
        circle = patches.Circle((x, y), radius=radius, edgecolor='red', facecolor='none', lw=1.5)
        ax.add_patch(circle)
        
        # Prepare annotation text
        mag_40 = row.get('MAG_CALIB_APER_40PX')
        posterior = row.get('posterior_O')
        
        texts = []
        if not pd.isna(mag_40):
            texts.append(f"{mag_40:.2f}")
        if not pd.isna(posterior):
            texts.append(f"P: {posterior:.3f}")
            
        if texts:
            annotation_text = "\n".join(texts)
            # Offset text slightly below the circle
            ax.text(x, y - radius - 2, annotation_text, color='blue', fontsize=8,
                    ha='center', va='top', bbox=dict(facecolor='white', alpha=0.6, edgecolor='none', pad=1))

# Save output
plt.savefig(OUTPUT_PATH, dpi=300, bbox_inches='tight')
print(f"Plot saved successfully to {OUTPUT_PATH}")

