import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.modeling import models, fitting
import matplotlib

matplotlib.use('Agg')

def evaluate_epsf_with_moffat(epsf_data, frb_name, out_dir):
    y_shape, x_shape = epsf_data.shape
    y, x = np.mgrid[:y_shape, :x_shape]
    
    center_y, center_x = y_shape // 2, x_shape // 2
    max_amp = np.max(epsf_data)
    
    moffat_init = models.Moffat2D(amplitude=max_amp, 
                                  x_0=center_x, y_0=center_y, 
                                  gamma=2.0, alpha=3.0)
    
    fitter = fitting.LevMarLSQFitter()
    moffat_fit = fitter(moffat_init, x, y, epsf_data)
    
    model_data = moffat_fit(x, y)
    residual_data = epsf_data - model_data
    max_fractional_residual = np.max(np.abs(residual_data)) / max_amp
    
    fit_alpha = moffat_fit.alpha.value
    fit_gamma = moffat_fit.gamma.value
    
    try:
        analytical_fwhm = 2 * fit_gamma * np.sqrt(2**(1/fit_alpha) - 1)
    except:
        analytical_fwhm = np.nan
    
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
    
    im1 = ax1.imshow(epsf_data, origin='lower', cmap='viridis')
    ax1.set_title("Input ePSF")
    plt.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04)
    
    im2 = ax2.imshow(model_data, origin='lower', cmap='viridis')
    ax2.set_title(f"Moffat Model\n(alpha={fit_alpha:.2f})")
    plt.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)
    
    vmax = np.max(np.abs(residual_data))
    im3 = ax3.imshow(residual_data, origin='lower', cmap='RdBu', vmin=-vmax, vmax=vmax)
    ax3.set_title("Residual (Data - Model)")
    plt.colorbar(im3, ax=ax3, fraction=0.046, pad=0.04)
    
    plt.suptitle(f"FRB {frb_name} PSF Moffat Fit")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{frb_name}_moffat_fit.png"), bbox_inches='tight', dpi=150)
    plt.close()

    return {
        'frb_name': frb_name,
        'moffat_alpha': fit_alpha,
        'moffat_gamma': fit_gamma,
        'moffat_fwhm': analytical_fwhm,
        'max_frac_resid_pct': max_fractional_residual * 100
    }

def main():
    base_dir = "Galfit/runs"
    out_dir = "psfs/moffat_diagnostics"
    os.makedirs(out_dir, exist_ok=True)
    
    results = []
    
    frb_dirs = glob.glob(os.path.join(base_dir, "*"))
    for frb_dir in frb_dirs:
        frb_name = os.path.basename(frb_dir)
        psf_dir = os.path.join(frb_dir, "with_psf")
        
        if not os.path.isdir(psf_dir):
            continue
            
        psf_files = glob.glob(os.path.join(psf_dir, "*down*.fits"))
        if not psf_files:
            # Fallback to general psf if down isn't found
            psf_files = glob.glob(os.path.join(psf_dir, "*psf*.fits"))
            if not psf_files:
                continue
            
        target_psf = psf_files[0]

                
        try:
            with fits.open(target_psf) as hdul:
                # Find first valid data extension
                psf_data = None
                for hdu in hdul:
                    if hdu.data is not None and isinstance(hdu.data, np.ndarray) and len(hdu.data.shape) == 2:
                        psf_data = hdu.data
                        break
                        
                if psf_data is not None:
                    print(f"Evaluating Moffat profile for {frb_name} using {os.path.basename(target_psf)}...")
                    stats = evaluate_epsf_with_moffat(psf_data, frb_name, out_dir)
                    results.append(stats)
        except Exception as e:
            print(f"Error processing {frb_name}: {e}")
            
    if results:
        df_new = pd.DataFrame(results)
        
        csv_path = "psf_fwhm_summary.csv"
        if os.path.exists(csv_path):
            df_old = pd.read_csv(csv_path)
            # rename frb column if it exists under a different name
            merge_col = None
            if 'FRB' in df_old.columns: merge_col = 'FRB'
            elif 'frb_name' in df_old.columns: merge_col = 'frb_name'
            
            if merge_col:
                # To handle potential mismatches, rename the new df key temporarily
                df_new = df_new.rename(columns={'frb_name': merge_col})
                # Merge
                merged = pd.merge(df_old, df_new, on=merge_col, how='left')
                merged.to_csv(csv_path, index=False)
                print(f"Appended Moffat statistics to {csv_path}")
            else:
                df_new.to_csv("moffat_fwhm_stats.csv", index=False)
                print("Could not merge with psf_fwhm_summary.csv (no matching ID column). Saved to moffat_fwhm_stats.csv")
        else:
            df_new.to_csv("moffat_fwhm_stats.csv", index=False)
            print("psf_fwhm_summary.csv not found. Saved standalone to moffat_fwhm_stats.csv")

if __name__ == "__main__":
    main()
