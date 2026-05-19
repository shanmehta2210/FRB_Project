
import os
import csv
import shutil
import numpy as np
from astropy.io import fits
import glob

# Constants
INPUT_DIR = "large_cutouts"
OUTPUT_DIR = "cropped_host_galaxies"
CSV_FILE = "master_frb_summary.csv"

def main():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
        print(f"Created output directory: {OUTPUT_DIR}")

    # Read CSV
    crops = {}
    with open(CSV_FILE, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            frb = row['FRB']
            # Check if coordinates are present
            if row['xmin'] and row['xmax'] and row['ymin'] and row['ymax']:
                try:
                    crops[frb] = {
                        'xmin': int(float(row['xmin'])),
                        'xmax': int(float(row['xmax'])),
                        'ymin': int(float(row['ymin'])),
                        'ymax': int(float(row['ymax']))
                    }
                except ValueError:
                    print(f"Skipping {frb}: Invalid coordinates.")
            else:
                print(f"Skipping {frb}: Missing coordinates (Notes: {row.get('notes', '')})")

    print(f"Found {len(crops)} images to crop.")

    # Process files
    for frb, crop in crops.items():
        flux_file = os.path.join(INPUT_DIR, f"{frb}_flux.fits")
        psf_file = os.path.join(INPUT_DIR, f"{frb}_flux_psf.fits")
        
        if not os.path.exists(flux_file):
            print(f"Warning: File not found {flux_file}")
            continue

        print(f"Processing {frb}...")
        
        # 1. Crop Flux Image
        try:
            with fits.open(flux_file) as hdu:
                data = hdu[0].data
                header = hdu[0].header
                
                # DS9 coordinates are 1-based, inclusive.
                # Numpy is 0-based, exclusive at the end.
                # xmin, xmax (DS9) -> xmin-1 : xmax (Numpy)
                
                x1 = crop['xmin'] - 1
                x2 = crop['xmax']
                y1 = crop['ymin'] - 1
                y2 = crop['ymax']
                
                # Bounds check
                if x1 < 0 or y1 < 0 or x2 > data.shape[1] or y2 > data.shape[0]:
                     print(f"  Error: Crop bounds [{x1}:{x2}, {y1}:{y2}] out of image size {data.shape[::-1]}")
                     continue

                cropped_data = data[y1:y2, x1:x2]
                
                # Update header WCS CRPIX if present (Shift the reference pixel)
                if 'CRPIX1' in header:
                    header['CRPIX1'] -= x1
                if 'CRPIX2' in header:
                    header['CRPIX2'] -= y1
                
                output_flux = os.path.join(OUTPUT_DIR, f"{frb}_flux.fits")
                fits.writeto(output_flux, cropped_data, header, overwrite=True)
                print(f"  Saved crop to {output_flux}")
                
        except Exception as e:
            print(f"  Error processing flux file: {e}")
            continue

        # 2. Copy PSF File
        if os.path.exists(psf_file):
            output_psf = os.path.join("psfs", f"{frb}_flux_psf.fits")
            shutil.copy2(psf_file, output_psf)
            print(f"  Copied PSF to {output_psf}")
        else:
            print(f"  Warning: PSF file not found {psf_file}")

    print("Done.")

if __name__ == "__main__":
    main()
