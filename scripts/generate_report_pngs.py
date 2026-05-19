import os
import glob
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.visualization import MinMaxInterval, LinearStretch, ImageNormalize
import matplotlib

matplotlib.use('Agg') # Off-screen rendering

def save_scaled_png(data, out_path, vmin=None, vmax=None):
    if vmin is not None and vmax is not None:
        norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=LinearStretch())
    else:
        norm = ImageNormalize(data, interval=MinMaxInterval(), stretch=LinearStretch())
    plt.figure(figsize=(5, 5))
    plt.imshow(data, cmap='gray', origin='lower', norm=norm)
    plt.axis('off')
    plt.savefig(out_path, bbox_inches='tight', pad_inches=0, dpi=150)
    plt.close()

def main():
    base_dir = "Bias_Report"
    os.makedirs(base_dir, exist_ok=True)
    
    runs_dir = os.path.join("Galfit", "runs")
    frb_dirs = glob.glob(os.path.join(runs_dir, "*"))
    
    for frb_dir in frb_dirs:
        frb_name = os.path.basename(frb_dir)
        if not os.path.isdir(frb_dir):
            continue
            
        for fit_type in ['no_psf', 'with_psf']:
            fit_dir = os.path.join(frb_dir, fit_type)
            if not os.path.isdir(fit_dir):
                continue
                
            out_fits_path = os.path.join(fit_dir, "out.fits")
            if os.path.exists(out_fits_path):
                try:
                    with fits.open(out_fits_path) as hdul:
                        # Assuming standard GALFIT format: [1]=Image, [2]=Model, [3]=Residual
                        if len(hdul) >= 4:
                            image_data = hdul[1].data
                            model_data = hdul[2].data
                            resid_data = hdul[3].data
                            
                            prefix = f"{frb_name}_{fit_type}"
                            
                            # Use image min/max to scale model and residual identically
                            img_vmin = float(image_data.min())
                            img_vmax = float(image_data.max())
                            
                            save_scaled_png(image_data, os.path.join(base_dir, f"{prefix}_image.png"))
                            save_scaled_png(model_data, os.path.join(base_dir, f"{prefix}_model.png"), vmin=img_vmin, vmax=img_vmax)
                            save_scaled_png(resid_data, os.path.join(base_dir, f"{prefix}_residual.png"), vmin=img_vmin, vmax=img_vmax)
                            
                            print(f"Processed out.fits for {frb_name} ({fit_type})")
                except Exception as e:
                    print(f"Error processing {out_fits_path}: {e}")
            
            # Look for downsampled PSF
            psf_files = glob.glob(os.path.join(fit_dir, "*psf*.fits"))
            for psf_file in psf_files:
                try:
                    with fits.open(psf_file) as hdul:
                        # Find first extension with data
                        psf_data = None
                        for hdu in hdul:
                            if hdu.data is not None:
                                psf_data = hdu.data
                                break
                        
                        if psf_data is not None:
                            psf_name = os.path.splitext(os.path.basename(psf_file))[0]
                            out_path = os.path.join(base_dir, f"{frb_name}_{fit_type}_{psf_name}.png")
                            if not os.path.exists(out_path): # Don't overwrite if multiple matched
                                save_scaled_png(psf_data, out_path)
                                print(f"Processed PSF {psf_file} for {frb_name} ({fit_type})")
                except Exception as e:
                    print(f"Error processing {psf_file}: {e}")

if __name__ == "__main__":
    main()
