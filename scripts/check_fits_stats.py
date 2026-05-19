import os
from astropy.io import fits
import numpy as np

def check_stats(fits_path):
    if not os.path.exists(fits_path):
        print(f"File not found: {fits_path}")
        return

    print(f"Stats for {fits_path}:")
    try:
        with fits.open(fits_path) as hdul:
            for i, hdu in enumerate(hdul):
                print(f"  Extension {i} ({hdu.name}):")
                if hdu.data is None:
                    print("    No data")
                    continue
                
                data = hdu.data
                print(f"    Shape: {data.shape}")
                print(f"    Min: {np.nanmin(data):.4e}")
                print(f"    Max: {np.nanmax(data):.4e}")
                print(f"    Mean: {np.nanmean(data):.4e}")
                print(f"    Std: {np.nanstd(data):.4e}")
                
                # Check for zero count
                zeros = np.sum(data == 0)
                print(f"    Zero count: {zeros} / {data.size} ({zeros/data.size:.1%})")
                
    except Exception as e:
        print(f"Error reading FITS: {e}")

if __name__ == "__main__":
    check_stats("Galfit/galfit_output/20220914A/out.fits")
