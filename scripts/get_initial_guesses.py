import os
import glob
import numpy as np
import warnings
import csv
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.segmentation import detect_sources, SourceCatalog

def calculate_robust_properties(data, segment, label_id, center):
    """
    Calculate r_e (half-light radius) and n (Sersic index) using 
    Curve of Growth on the segmented pixels.
    """
    # 1. Mask data to only this segment
    mask = (segment.data == label_id)
    # Get pixel values and coordinates
    y_indices, x_indices = np.where(mask)
    flux_values = data[mask]
    
    # 2. Sort pixels by distance from centroid
    xc, yc = center
    distances = np.sqrt((x_indices - xc)**2 + (y_indices - yc)**2)
    
    # Sort by distance
    sorted_indices = np.argsort(distances)
    sorted_flux = flux_values[sorted_indices]
    sorted_dist = distances[sorted_indices]
    
    # 3. Cumulative Flux
    cum_flux = np.cumsum(sorted_flux)
    total_flux = cum_flux[-1]
    
    if total_flux <= 0:
        return 3.0, 1.0 # Fallback
        
    # 4. Find radii
    # r_e is where cumulative flux is 50%
    
    def get_radius_at_fraction(fraction):
        target_flux = total_flux * fraction
        # Find index where cum_flux >= target_flux
        idx = np.searchsorted(cum_flux, target_flux)
        if idx < len(sorted_dist):
            return sorted_dist[idx]
        return sorted_dist[-1]
        
    r20 = get_radius_at_fraction(0.2)
    r50 = get_radius_at_fraction(0.5)
    r80 = get_radius_at_fraction(0.8)
    
    # Safety
    if r20 == 0: r20 = 0.5
    if r50 < 0.5: r50 = 0.5
    if r80 <= r20: r80 = r20 + 0.1
    
    # 5. Estimate n from Concentration
    # C = 5 * log10(r80 / r20)
    C = 5 * np.log10(r80 / r20)
    
    if C < 2.5: n_guess = 0.5
    elif C < 3.3: n_guess = 1.0 # Spiral / Disk
    elif C < 4.5: n_guess = 2.5 # Intermediate
    else: n_guess = 4.0       # Elliptical
    
    return float(r50), float(n_guess)

def get_initial_guesses(fits_filename, threshold_sigma=3.0, plot=False):
    """
    Automatically detects the central object in a FITS image and 
    returns initial guesses for GALFIT.
    """
    
    # 1. Load Data
    with fits.open(fits_filename) as hdu:
        data = hdu[0].data
        header = hdu[0].header
        zp = header.get('MAGZERO', 25.0) 
        if zp is None: zp = 25.0

    # 2. Estimate Background & Threshold
    mean, median, std = sigma_clipped_stats(data, sigma=3.0)
    threshold = median + (threshold_sigma * std)
    
    # 3. Detect Sources (Segmentation)
    segm = detect_sources(data, threshold, npixels=5)
    
    if segm is None:
        warnings.warn(f"No sources detected in {fits_filename}! Defaulting to center.")
        h, w = data.shape
        return {
            'x': w/2, 'y': h/2, 'mag': 20.0, 'r_e': 5.0, 
            'n': 1.0, 'axis_ratio': 0.5, 'pa': 0.0
        }

    # 4. Measure Properties
    cat = SourceCatalog(data, segm)
    
    # 5. Identify the "Target" (Closest to center)
    h, w = data.shape
    center = (w/2, h/2)
    
    distances = np.sqrt((cat.xcentroid - center[0])**2 + (cat.ycentroid - center[1])**2)
    target_idx = np.argmin(distances)
    target = cat[target_idx]
    
    # 6. Extract Parameters (With Robust Estimation)
    robust_re, robust_n = calculate_robust_properties(data, segm, target.label, (target.xcentroid, target.ycentroid))
    
    r_e_guess = robust_re
    n_guess = robust_n
    
    # Parametric fallback (just for b/a and PA)
    try:
        axis_ratio = 1.0 - target.ellipticity
        if axis_ratio <= 0.05: axis_ratio = 0.1
        pa_deg = target.orientation.value 
        if np.isnan(pa_deg): pa_deg = 0.0
    except:
        axis_ratio = 0.5
        pa_deg = 0.0

    # B. Magnitude
    flux = target.segment_flux
    if flux <= 0: 
        mag_guess = 25.0
    else:
        mag_guess = zp - 2.5 * np.log10(flux)

    # 7. Visualization
    if plot:
        import matplotlib.pyplot as plt
        from matplotlib.patches import Ellipse
        
        plt.figure(figsize=(6,6))
        plt.imshow(data, origin='lower', cmap='gray', vmin=median, vmax=threshold + 5*std)
        plt.plot(target.xcentroid, target.ycentroid, 'rx', ms=10, label='Center')
        
        ell = Ellipse(xy=(target.xcentroid, target.ycentroid),
                      width=r_e_guess*2, height=r_e_guess*2*axis_ratio,
                      angle=pa_deg, edgecolor='red', facecolor='none', lw=2)
        plt.gca().add_patch(ell)
        plt.title(f"Mag: {mag_guess:.1f}, Re: {r_e_guess:.1f}, n: {n_guess}, C: {5 * np.log10(robust_re*2 / robust_re):.2f}") # C is Approx
        plt.legend()
        plt.show()

    # 8. Return Dictionary
    return {
        'x': target.xcentroid,
        'y': target.ycentroid,
        'mag': mag_guess,
        'r_e': r_e_guess,
        'n': n_guess,
        'axis_ratio': axis_ratio,
        'pa': pa_deg
    }

def process_directory(input_dir, output_csv="initial_guesses.csv"):
    """
    Process all galaxy files in the directory and save guesses to CSV.
    """
    search_pattern = os.path.join(input_dir, "*_flux.fits")
    files = sorted(glob.glob(search_pattern))
    print(f"Found {len(files)} files in {input_dir}")
    
    exclude_frbs = []
    
    results = []
    
    for f in files:
        basename = os.path.basename(f)
        frb_name = basename.replace("_flux.fits", "")
        
        if frb_name in exclude_frbs:
            print(f"Skipping {basename} (Excluded)")
            continue
            
        try:
            print(f"Processing {basename}...")
            guesses = get_initial_guesses(f, plot=False)
            guesses['filename'] = basename
            results.append(guesses)
        except Exception as e:
            print(f"Error processing {basename}: {e}")
            import traceback
            traceback.print_exc()
            
    if results:
        keys = ['filename', 'x', 'y', 'mag', 'r_e', 'n', 'axis_ratio', 'pa']
        with open(output_csv, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=keys)
            writer.writeheader()
            writer.writerows(results)
        print(f"Saved initial guesses to {output_csv}")
    else:
        print("No results to save.")

if __name__ == "__main__":
    possible_dirs = [
        "cropped_host_galaxies",
        "../cropped_host_galaxies",
        "/cropped_host_galaxies"
    ]
    
    target_dir = None
    for d in possible_dirs:
        if os.path.isdir(d):
            target_dir = d
            break
            
    if target_dir:
        process_directory(target_dir, output_csv="initial_guesses.csv")
    else:
        print("Could not find 'cropped_host_galaxies' directory.")
