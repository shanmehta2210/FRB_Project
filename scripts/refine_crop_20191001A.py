
import os
from astropy.io import fits
from astropy.wcs import WCS

INPUT_FILE = "cropped_host_galaxies/20191001A_flux.fits"

def main():
    if not os.path.exists(INPUT_FILE):
        print(f"Error: {INPUT_FILE} not found.")
        return

    print(f"Refining crop for {INPUT_FILE}...")
    
    # Read data first, then close handles to release file lock on Windows
    with fits.open(INPUT_FILE) as hdul:
        data = hdul[0].data.copy()  # Copy data to keep it after closing
        header = hdul[0].header
        wcs = WCS(header)
        print(f"Original shape: {data.shape}")
        
    # Crop x(0, 40) and y(0, 42)
    # Assuming user means 0-based indices or range 0 to 40/42
    # Python slice: y[0:42], x[0:40]
    x_start, x_end = 0, 40
    y_start, y_end = 0, 42
    
    cropped_data = data[y_start:y_end, x_start:x_end]
    cropped_wcs = wcs[y_start:y_end, x_start:x_end]
    
    print(f"New shape: {cropped_data.shape}")
    
    new_header = cropped_wcs.to_header()
    
    # Save
    new_hdu = fits.PrimaryHDU(data=cropped_data, header=new_header)
    new_hdu.writeto(INPUT_FILE, overwrite=True)
    print(f"Overwrote {INPUT_FILE} with refined crop.")

if __name__ == "__main__":
    main()
