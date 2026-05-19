import os
import glob
import numpy as np
from astropy.io import fits

def normalize_psfs(file_list):
    count = 0
    for psf_path in file_list:
        with fits.open(psf_path, mode='update') as hdul:
            data = hdul[0].data
            total_flux = np.sum(data)
            
            if np.isclose(total_flux, 1.0, atol=1e-4):
                continue
                
            print(f"Normalizing {os.path.basename(psf_path)} (Original Sum: {total_flux:.4f})")
            hdul[0].data = data / total_flux
            hdul.flush()
            count += 1
    return count

def main():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    # Target 4x PSFs
    psfs_4x_dir = os.path.join(base_dir, "psfs")
    psfs_4x_files = glob.glob(os.path.join(psfs_4x_dir, "*_psf.fits"))
    # Exclude anything in subdirs or anything that explicitly says 1x (just in case)
    psfs_4x_files = [f for f in psfs_4x_files if os.path.isfile(f) and "_1x_" not in f]
    
    # Target 1x PSFs
    psfs_1x_dir = os.path.join(base_dir, "psfs", "downsampled_psfs")
    psfs_1x_files = glob.glob(os.path.join(psfs_1x_dir, "*_1x_psf.fits"))
    
    count_4x = normalize_psfs(psfs_4x_files)
    count_1x = normalize_psfs(psfs_1x_files)
            
    print(f"Successfully normalized {count_4x} 4x PSFs to exactly 1.0")
    print(f"Successfully normalized {count_1x} 1x PSFs to exactly 1.0")

if __name__ == "__main__":
    main()
