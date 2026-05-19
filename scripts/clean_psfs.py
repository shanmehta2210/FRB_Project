import os
import glob
import numpy as np
from astropy.io import fits

def clean_psf(psf_path):
    with fits.open(psf_path, mode='update') as hdul:
        data = hdul[0].data
        
        # 1. Measure and subtract the local background pedestal from the outer 5-pixel edge
        mask = np.ones(data.shape, dtype=bool)
        mask[5:-5, 5:-5] = False
        pedestal = np.median(data[mask])
        data = data - pedestal
        
        # 2. Hard clip ANY negative pixels to exactly 0.0
        # GALFIT convolution diverges violently if PSFs contain negative values
        data[data < 0.0] = 0.0
        
        # 3. Strictly normalize the strictly positive array to 1.0
        total_flux = np.sum(data)
        data = data / total_flux
        
        hdul[0].data = data
        hdul.flush()
        print(f"Cleaned {os.path.basename(psf_path)}: Pedestal={pedestal:.6f}, Negative pixels removed, Sum=1.0")

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    psf_dir = os.path.join(base_dir, "..", "psfs")
    
    psf_files = glob.glob(os.path.join(psf_dir, "*_flux_psf.fits"))
    for file in psf_files:
        clean_psf(file)
        
    print(f"Successfully cleaned {len(psf_files)} 4x PSFs.")

if __name__ == "__main__":
    main()
