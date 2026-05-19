import os
import sys
import numpy as np
import pyvo
import warnings
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.table import Table
from astropy.nddata import NDData
from photutils.psf import extract_stars, EPSFBuilder, EPSFStars

warnings.simplefilter('ignore', category=FITSFixedWarning)
warnings.simplefilter('ignore', category=RuntimeWarning)

# Hardcoded exclusions from manual user inspection (1-indexed based on _initial plot)
MANUAL_EXCLUSIONS_CATALOG = {
    "20190608B": [4],
    "20211203C": [10],
    "20220310F": [4, 8, 9],
    "20220509G": [1],
    "20221012A": [3, 8],
}

def _save_star_grid(frb_name, extracted_stars, fluxes, tag, title, original_indices=None):
    if len(extracted_stars) == 0:
        return
    import matplotlib.pyplot as plt
    from astropy.visualization import simple_norm
    n_stars = len(extracted_stars)
    n_cols = min(5, n_stars)
    n_rows = (n_stars + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 3 * n_rows))
    fig.subplots_adjust(hspace=0.5, wspace=0.3)
    if n_rows == 1 and n_cols == 1:
        axes = [axes]
    axes = list(np.array(axes).flatten())
    for i, ax in enumerate(axes):
        if i < n_stars:
            star = extracted_stars[i]
            norm = simple_norm(star.data, 'log', percent=99.)
            ax.imshow(star.data, origin='lower', cmap='viridis', norm=norm)
            flux = fluxes[i]
            if original_indices is not None:
                orig = original_indices[i] if i < len(original_indices) else '?'
                label = f"#{i+1} (was #{orig})"
            else:
                label = f"Star {i+1}"
            ax.set_title(f"{label}\nFlux: {flux:.1f}", fontsize=8)
            ax.axis('off')
        else:
            ax.axis('off')
    fig.suptitle(title, fontsize=14)
    diag_dir = os.path.join("psfs", "diagnostic")
    os.makedirs(diag_dir, exist_ok=True)
    outpath = os.path.join(diag_dir, f"{frb_name}_catalog_stars_{tag}.png")
    plt.savefig(outpath, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"Saved {tag} diagnostic to {outpath}")

def get_catalog_stars(ra_center, dec_center, radius_deg=0.1):
    """Query NOIRLab's TAP service for DESI Legacy Survey DR10 Tractor objects."""
    print(f"Querying Legacy Survey DR10 Tractor Catalog near RA={ra_center:.4f}, Dec={dec_center:.4f}...")
    
    # Expand RA range to account for cos(dec) compression at high latitudes
    ra_range = radius_deg / np.cos(np.radians(np.clip(dec_center, -85, 85)))
    ra_min = ra_center - ra_range
    ra_max = ra_center + ra_range
    dec_min = dec_center - radius_deg
    dec_max = dec_center + radius_deg

    # Build RA clause, handling wraparound at 0/360
    if ra_min < 0:
        ra_clause = f"(ra > {ra_min + 360} OR ra < {ra_max})"
    elif ra_max > 360:
        ra_clause = f"(ra > {ra_min} OR ra < {ra_max - 360})"
    else:
        ra_clause = f"ra > {ra_min} AND ra < {ra_max}"
    
    query = f"""
    SELECT TOP 30 ra, dec, type, flux_r, fracflux_r, anymask_r
    FROM ls_dr10.tractor
    WHERE {ra_clause}
      AND dec > {dec_min} AND dec < {dec_max}
      AND type = 'PSF'
      AND fracflux_r < 0.05
      AND anymask_r = 0
      AND flux_r > 5
    ORDER BY flux_r DESC
    """
    
    service = pyvo.dal.TAPService("https://datalab.noirlab.edu/tap")
    try:
        result = service.search(query)
        table = result.to_table()
        print(f"Found {len(table)} pristine stars passing quality cuts.")
        return table
    except Exception as e:
        print(f"Error querying catalog: {e}")
        return None

def build_psf(file_path, output_suffix="_flux_psf.fits", oversampling=4, size=25, max_stars=25):
    frb_name = os.path.basename(file_path).split('_')[0]
    print(f"\n=== Building PSF for {frb_name} using Survey Catalog ===")
    
    try:
        with fits.open(file_path) as hdu:
            data = hdu[0].data
            header = hdu[0].header
            
            # Replace NaNs
            if np.isnan(data).any():
                data = np.nan_to_num(data, nan=np.nanmedian(data))
                
            # Perform Local Sky Subtraction to strip background pedestal
            from astropy.stats import sigma_clipped_stats
            _, median, _ = sigma_clipped_stats(data, sigma=3.0)
            data_sub = data - median
                
            wcs = WCS(header)
            
            # Get center coords in RA/Dec based on image size
            ny, nx = data.shape
            center_x, center_y = nx / 2, ny / 2
            ra_center, dec_center = wcs.all_pix2world(center_x, center_y, 0)
            
            # Try to infer pixel scale to roughly set search radius
            pixel_scale_deg = np.abs(wcs.pixel_scale_matrix[0,0]) if hasattr(wcs, 'pixel_scale_matrix') else 0.262/3600
            radius_deg = (max(nx, ny) / 2) * pixel_scale_deg * 1.5 # 1.5x padding
            
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        return

    # 1. Fetch catalog
    cat_stars = get_catalog_stars(float(ra_center), float(dec_center), radius_deg=radius_deg)
    if cat_stars is None or len(cat_stars) == 0:
        print("No suitable PSF objects found in the catalog!")
        return

    # Limit to top N stars
    cat_stars = cat_stars[:max_stars]

    # 2. Convert RA/Dec to Pixel Coordinates
    x_pix, y_pix = wcs.all_world2pix(cat_stars['ra'], cat_stars['dec'], 0)
    
    # Filter points outside image or too close to edge
    valid = (x_pix > size) & (x_pix < nx - size) & (y_pix > size) & (y_pix < ny - size)
    x_pix = x_pix[valid]
    y_pix = y_pix[valid]
    cat_stars = cat_stars[valid]
    
    print(f"Kept {len(valid.nonzero()[0])} stars within the image boundaries.")
    if len(x_pix) == 0:
        print("No stars inside image area!")
        return
        
    # Create Table for EPSF extraction
    stars_tbl = Table()
    stars_tbl['x'] = x_pix
    stars_tbl['y'] = y_pix
    
    # 3. Extract Stars initially
    nddata = NDData(data=data_sub)
    extracted_stars = extract_stars(nddata, catalogs=stars_tbl, size=size)
    valid_fluxes = cat_stars['flux_r']

    # --- SAVE INITIAL DIAGNOSTIC ---
    _save_star_grid(frb_name, extracted_stars, valid_fluxes, tag="initial",
                    title=f"{frb_name} - Initial Catalog Stars (use these numbers for exclusions)")

    # --- MANUAL EXCLUSIONS ---
    original_indices = list(range(1, len(extracted_stars) + 1))
    if frb_name in MANUAL_EXCLUSIONS_CATALOG:
        exclude_1_based = MANUAL_EXCLUSIONS_CATALOG[frb_name]
        print(f"Applying manual exclusions for {frb_name}: Removing stars {exclude_1_based} (1-indexed)")
        
        exclude_0_based = sorted([i - 1 for i in exclude_1_based], reverse=True)
        
        # Convert to lists to pop
        star_list = list(extracted_stars)
        flux_list = list(valid_fluxes)
        
        for idx in exclude_0_based:
            if 0 <= idx < len(star_list):
                print(f"  -> Popped catalog star {idx+1}")
                star_list.pop(idx)
                flux_list.pop(idx)
                original_indices.pop(idx)
                
        extracted_stars = EPSFStars(star_list)
        valid_fluxes = flux_list

    print(f"Final Selection: {len(extracted_stars)} catalog stars for PSF building.")
    if len(extracted_stars) == 0:
        print("No suitable stars left after exclusions.")
        return

    # --- SAVE FINAL DIAGNOSTIC ---
    _save_star_grid(frb_name, extracted_stars, valid_fluxes, tag="final",
                    title=f"{frb_name} - Final Catalog Stars Used for EPSF",
                    original_indices=original_indices)
        
    # 4. Build EPSF
    print("Building EPSF...")
    epsf_builder = EPSFBuilder(oversampling=oversampling, maxiters=10, progress_bar=False)
    epsf, fitted_stars = epsf_builder(extracted_stars)
    
    # 5. Save Result
    output_filename = os.path.join("psfs", f"{frb_name}_flux_psf.fits")
    # Ensure strictly normalized
    epsf.data = epsf.data / np.sum(epsf.data)
    fits.writeto(output_filename, epsf.data, overwrite=True)
    print(f"Saved PSF to {output_filename}")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python build_catalog_psf.py <fits_file>")
        sys.exit(1)
        
    os.makedirs("psfs", exist_ok=True)
    build_psf(sys.argv[1])


