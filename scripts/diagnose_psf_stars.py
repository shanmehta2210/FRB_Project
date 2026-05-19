import os
import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.nddata import NDData
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from photutils.psf import extract_stars
import warnings
from astropy.utils.exceptions import AstropyWarning

warnings.simplefilter('ignore', category=AstropyWarning)

import glob

def get_all_frbs():
    # Only the 8 PanSTARRS-sourced FRBs that use DAOStarFinder fallback
    return [
        "20171020A", "20210807D", "20211127I", "20220207C",
        "20220307B", "20220319D", "20220825A", "20220912A"
    ]

TARGET_FRBS = get_all_frbs()

def visualize_selected_stars(frb_name, input_dir="large_cutouts", output_dir="psfs/diagnostic"):
    # Replicate exact logic from build_psf.py to get the specific stars it chose
    file_path = os.path.join(input_dir, f"{frb_name}_flux.fits")
    
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return
        
    print(f"\n--- Diagnosing {frb_name} ---")
    
    try:
        with fits.open(file_path) as hdul:
            data = hdul[0].data
            if data is None and len(hdul) > 1:
                data = hdul[1].data
    except Exception as e:
        print(f"Error loading: {e}")
        return

    # Clean data
    if np.any(np.isnan(data)):
        med_val = np.nanmedian(data)
        data = np.nan_to_num(data, nan=med_val)

    # Background
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    data_sub = data - median

    # Find stars (current parameters in build_psf.py)
    threshold = 35.0 * std
    daofind = DAOStarFinder(fwhm=4.0, threshold=threshold) 
    sources = daofind(data_sub)
    
    if sources is None or len(sources) == 0:
        threshold = 15.0 * std
        daofind = DAOStarFinder(fwhm=4.0, threshold=threshold)
        sources = daofind(data_sub)
        if sources is None or len(sources) == 0:
            print("No stars found.")
            return

    if 'xcentroid' in sources.colnames:
        sources.rename_column('xcentroid', 'x')
    if 'ycentroid' in sources.colnames:
        sources.rename_column('ycentroid', 'y')
        
    # Standard spatial bounding
    extract_size = 25
    margin = extract_size
    mask = (
        (sources['x'] > margin) & (sources['x'] < data.shape[1] - margin) &
        (sources['y'] > margin) & (sources['y'] < data.shape[0] - margin)
    )
    sources = sources[mask]
    
    # Tighter Morphology filters (Relaxed from 0.20/0.45)
    if 'roundness1' in sources.colnames:
        mask_round = (np.abs(sources['roundness1']) < 0.30) & (np.abs(sources['roundness2']) < 0.30)
        sources = sources[mask_round]

    if 'sharpness' in sources.colnames:
        mask_sharp = (sources['sharpness'] > 0.40) & (sources['sharpness'] < 0.85)
        sources = sources[mask_sharp]

    # Existing Isolation Check
    if len(sources) > 1:
        from scipy.spatial.distance import cdist
        coords = np.array([sources['x'], sources['y']]).T
        dist_matrix = cdist(coords, coords)
        np.fill_diagonal(dist_matrix, np.inf)
        min_dists = np.min(dist_matrix, axis=1)
        mask_isolated = min_dists > extract_size
        sources = sources[mask_isolated]

    # Existing Saturation filter
    if len(sources) > 0:
        cleaned_sources = sources.copy()
        cleaned_sources.sort('peak')
        if len(cleaned_sources) > 5:
             cleaned_sources = cleaned_sources[:-1] 
        sources = cleaned_sources

    # --- NEW: Strict Pixel-Level Edge Isolation Check ---
    # Eagerly reject stars that have bright neighbors encroaching on their 25x25 box
    valid_sources = []
    # Use 10 sigma above background as the limit for the box edge
    edge_limit = 10.0 * std 
    
    for row in sources:
        x_c = int(np.round(row['x']))
        y_c = int(np.round(row['y']))
        half_box = extract_size // 2
        
        # 25x25 Cutout
        cutout = data_sub[y_c-half_box : y_c+half_box+1, x_c-half_box : x_c+half_box+1]
        
        if cutout.shape != (extract_size, extract_size):
            continue
            
        # Check the 2-pixel wide outer boundary of this cutout
        top_edge = cutout[-2:, :]
        bot_edge = cutout[:2, :]
        left_edge = cutout[:, :2]
        right_edge = cutout[:, -2:]
        
        # If the max pixel value on any edge is greater than edge_limit, reject as not isolated
        edge_max = max(np.max(top_edge), np.max(bot_edge), np.max(left_edge), np.max(right_edge))
        
        if edge_max < edge_limit:
            valid_sources.append(row)
            
    if len(valid_sources) > 0:
        from astropy.table import vstack
        sources = vstack(valid_sources)
    else:
        sources = sources[:0] # Empty table
        
    sources.sort('peak')
    sources.reverse()
    
    # Show ALL stars that passed the criteria — no cap.
    selection = sources
    
    print(f"Selected {len(selection)} strictly filtered stars to put into EPSFBuilder.")
        
    # Extract
    stars = extract_stars(NDData(data=data_sub), selection, size=extract_size)
    
    # --- PLOT THEM ---
    n_stars = len(stars)
    n_cols = 5
    n_rows = (n_stars + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(15, 3 * n_rows))
    fig.subplots_adjust(hspace=0.5, wspace=0.3)
    if n_rows == 1 and n_cols == 1:
        axes = np.array([axes])
    axes = axes.flatten()
    
    for i in range(len(axes)):
        if i < n_stars:
            star = stars[i]
            # Plot the star cutout
            from astropy.visualization import simple_norm
            norm = simple_norm(star.data, 'log', percent=99.)
            
            im = axes[i].imshow(star.data, origin='lower', cmap='viridis', norm=norm)
            
            # Subplot title info
            x_real = selection['x'][i]
            y_real = selection['y'][i]
            peak = selection['peak'][i]
            axes[i].set_title(f"Star {i+1}\n(x={x_real:.0f}, y={y_real:.0f})\nPeak: {peak:.1f}")
            axes[i].axis('off')
        else:
            axes[i].axis('off')
            
    fig.suptitle(f"{frb_name} - Selected Stars for EPSF", fontsize=16)
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    outpath = os.path.join(output_dir, f"{frb_name}_star_selection_initial.png")
    plt.savefig(outpath, bbox_inches='tight', dpi=150)
    plt.close()
    print(f"Saved visualization to {outpath}")

def main():
    for frb in TARGET_FRBS:
        visualize_selected_stars(frb)

if __name__ == "__main__":
    main()
