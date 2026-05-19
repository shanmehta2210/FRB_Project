import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.nddata import block_reduce
from astropy.visualization import simple_norm

def downsample_psf(input_array, factor=4, output_filename="psf_native_vis.fits"):
    """
    Downsamples an oversampled PSF to native resolution and saves it.
    
    Args:
        input_array (numpy.ndarray): The 2D oversampled PSF data.
        factor (int): The oversampling factor to reverse (default 4).
        output_filename (str): Name of the output FITS file.
    """
    # 1. Downsample using sum to conserve flux
    # This takes every factor x factor grid of pixels and sums them into 1 pixel
    native_psf = block_reduce(input_array, block_size=factor, func=np.sum)
    
    # 2. Sanity Check: Ensure flux is still ~1.0
    original_flux = np.sum(input_array)
    
    # 2.5: Strip the 11% background pedestal baked into the 25x25 array
    # Measure the outer 2-pixel border
    mask = np.ones(native_psf.shape, dtype=bool)
    mask[2:-2, 2:-2] = False
    pedestal = np.median(native_psf[mask])
    
    native_psf = native_psf - pedestal
    native_psf[native_psf < 0.0] = 0.0
    
    new_flux = np.sum(native_psf)
    print(f"  Original Flux: {original_flux:.5f}")
    print(f"  New Flux:      {new_flux:.5f}")
    
    # 3. Restrict flux explicitly to 1.0 to prevent GALFIT bloating
    native_psf = native_psf / new_flux
    
    # 4. Write to FITS
    fits.writeto(output_filename, native_psf, overwrite=True)
    return native_psf

def main():
    base_dir = "psfs"
    output_dir = os.path.join(base_dir, "downsampled_psfs")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")
        
    psf_files = sorted(glob.glob(os.path.join(base_dir, "*_flux_psf.fits")))
    print(f"Found {len(psf_files)} PSF FITS files to process.")
    
    for psf_file in psf_files:
        filename = os.path.basename(psf_file)
        frb_name = filename.split("_flux_psf.fits")[0]
        print(f"\nProcessing {frb_name}...")
        
        try:
            with fits.open(psf_file) as hdu:
                data = hdu[0].data
                
            out_fits = os.path.join(output_dir, f"{frb_name}_1x_psf.fits")
            native_psf = downsample_psf(data, factor=4, output_filename=out_fits)
            
            # Create PNG visualization
            out_png = os.path.join(output_dir, f"{frb_name}_1x_psf.png")
            fig, ax = plt.subplots(figsize=(6, 6))
            
            # Use Log or Asinh stretching to enhance faint wings
            norm = simple_norm(native_psf, 'log', percent=99.9)
            
            im = ax.imshow(native_psf, origin='lower', cmap='viridis', norm=norm)
            fig.colorbar(im, ax=ax, label='Flux')
            ax.set_title(f"Native PSF (1x) - {frb_name}")
            
            plt.savefig(out_png, dpi=150, bbox_inches='tight')
            plt.close(fig)
            print(f"  Saved PNG to {out_png}")
            
        except Exception as e:
            print(f"  Error processing {filename}: {e}")

if __name__ == "__main__":
    main()
