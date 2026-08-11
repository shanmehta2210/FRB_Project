import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.nddata import block_reduce
import csv

def compute_fwhm(profile, r_vals):
    """
    Compute FWHM given a 1D profile and corresponding radius values.
    Assumes the peak is at r=0.
    """
    peak = np.max(profile)
    half_max = peak / 2.0
    
    # Find points where profile crosses half_max
    # We take the first crossing looking outwards from the peak
    # Peak is at index 0
    idx = np.where(profile < half_max)[0]
    if len(idx) == 0:
        return np.nan # Couldn't find half max
    
    cross_idx = idx[0]
    
    # Linear interpolation for better precision
    if cross_idx > 0:
        r1, r2 = r_vals[cross_idx-1], r_vals[cross_idx]
        p1, p2 = profile[cross_idx-1], profile[cross_idx]
        
        # p = m*r + c => r = (p - c)/m
        # m = (p2-p1)/(r2-r1)
        r_hm = r1 + (half_max - p1) * (r2 - r1) / (p2 - p1)
    else:
        r_hm = r_vals[0]
        
    return r_hm * 2 # FWHM is 2 * HWHM



def extract_radial_profile(image, center, angle_deg, max_radius):
    """
    Extract a radial profile from an image at a specific angle.
    """
    cx, cy = center
    r_vals = np.arange(0, max_radius, 1.0) # step by 1 pixel
    
    angle_rad = np.deg2rad(angle_deg)
    
    x_vals = cx + r_vals * np.cos(angle_rad)
    y_vals = cy + r_vals * np.sin(angle_rad)
    
    # Use closest pixel integer coordinates
    x_idx = np.round(x_vals).astype(int)
    y_idx = np.round(y_vals).astype(int)
    
    # Valid indices check
    valid = (x_idx >= 0) & (x_idx < image.shape[1]) & (y_idx >= 0) & (y_idx < image.shape[0])
    
    extracted_profile = image[y_idx[valid], x_idx[valid]]
    
    return r_vals[valid], extracted_profile



def main():
    psf_dir = "psfs"
    plot_dir = os.path.join(psf_dir, "fwhm")
    
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    
    psf_files = sorted(glob.glob(os.path.join(psf_dir, "*_psf.fits")))
    print(f"Found {len(psf_files)} PSF files.")
    
    results = []
    angles = [0, 45, 90, 135]
    
    for psf_file in psf_files:
        frb_name = os.path.basename(psf_file).replace("_flux_psf.fits", "")
        print(f"Analyzing {frb_name}...")
        
        with fits.open(psf_file) as hdu:
            data = hdu[0].data
            
        if data is None:
            print(f"Skipping {frb_name}, no data found.")
            continue
            
        # 1. Handle oversampling by 4x. We use block_reduce to downsample.
        # This averages blocks of 4x4 into 1 pixel.
        try:
           downsampled_data = block_reduce(data, 4, func=np.mean)
        except Exception as e:
           print(f"Failed to downsample {frb_name}: {e}")
           continue
           
        # Find the center (assuming center pixel is the brightest or physically center)
        # Using peak pixel for center
        cy, cx = np.unravel_index(np.argmax(downsampled_data), downsampled_data.shape)
        center = (cx, cy)
        
        max_r = min(cx, cy, downsampled_data.shape[1]-cx, downsampled_data.shape[0]-cy)
        if max_r < 2:
            print(f"Warning: {frb_name} has very a small bounding box.")
            max_r = min(downsampled_data.shape[0], downsampled_data.shape[1]) // 2 
            
        fig, ax = plt.subplots(figsize=(8,6))
        
        fwhms = []
        for angle in angles:
            r_vals, profile = extract_radial_profile(downsampled_data, center, angle, max_r)
            
            # normalize profile to 1 for easier FWHM plotting visually
            if len(profile) > 0:
                norm_profile = profile / np.max(profile)
            else:
                norm_profile = profile
                
            ax.plot(r_vals, norm_profile, label=f"{angle} deg")
            
            fwhm_val = compute_fwhm(profile, r_vals)
            fwhms.append(fwhm_val)
            
        
        ax.axhline(0.5, color='gray', linestyle='--', label="Half-Max")
        ax.set_title(f"PSF Radial Profiles - {frb_name}")
        ax.set_xlabel("Radius (pixels, original scale)")
        ax.set_ylabel("Normalized Flux")
        ax.legend()
        ax.grid(True)
        
        plot_path = os.path.join(plot_dir, f"{frb_name}_psf_profile.png")
        plt.savefig(plot_path)
        plt.close(fig)
        
        avg_fwhm = np.nanmean(fwhms)
        print(f"  FWHM: {avg_fwhm:.2f} pixels")
        
        results.append({
            "FRB": frb_name,
            "FWHM_0": fwhms[0],
            "FWHM_45": fwhms[1],
            "FWHM_90": fwhms[2],
            "FWHM_135": fwhms[3],
            "Avg_FWHM": avg_fwhm
        })
        
    csv_file = "Archive/csv/psf/psf_fwhm_summary.csv"
    keys = ["FRB", "FWHM_0", "FWHM_45", "FWHM_90", "FWHM_135", "Avg_FWHM"]
    with open(csv_file, "w", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(results)
        
    print(f"Done. Saved summary to {csv_file}")

if __name__ == "__main__":
    main()
