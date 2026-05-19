
import os
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

# Configuration
LARGE_CUTOUTS_DIR = "large_cutouts"
OUTPUT_DIR = "cropped_host_galaxies"
CROPPING_CSV = "cropping_regions.csv"

# List of FRBs to update
TARGET_FRBS = ["20191001A", "20220319D"]

def crop_image(frb_name, xmin, xmax, ymin, ymax):
    # Construct file paths
    # Note: large cutouts usually have format {FRB}_flux.fits
    input_filename = f"{frb_name}_flux.fits"
    input_path = os.path.join(LARGE_CUTOUTS_DIR, input_filename)
    
    if not os.path.exists(input_path):
        print(f"Error: Input file at {input_path} not found.")
        return

    # Output path
    output_path = os.path.join(OUTPUT_DIR, input_filename)
    
    print(f"Processing {frb_name}...")
    print(f"  Input: {input_path}")
    print(f"  Crop: x[{xmin}:{xmax}], y[{ymin}:{ymax}]")

    try:
        with fits.open(input_path) as hdul:
            data = hdul[0].data
            header = hdul[0].header
            wcs = WCS(header)
            
            # Arrays are y, x in numpy (index 0 is y, index 1 is x)
            # Fits coordinates are 1-based, numpy is 0-based
            # The CSV likely uses 0-based or 1-based pixel coordinates from DS9?
            # DS9 uses 1-based coordinates usually.
            # Let's assume the CSV contains standard DS9 image coordinates (1-based), 
            # or if they are from the previous script logic, we need to match that.
            # Generally: slice = data[ymin-1 : ymax, xmin-1 : xmax]
            
            # Let's verify bounds
            if ymax > data.shape[0] or xmax > data.shape[1]:
                 print(f"  Warning: Crop bounds exceed image dimensions {data.shape}")
            
            # Perform crop
            # Ensure indices are integers
            x1, x2 = int(xmin), int(xmax)
            y1, y2 = int(ymin), int(ymax)
            
            # Assuming 1-based inclusive coordinates (standard DS9)
            # Numpy slice: [start:stop] where stop is exclusive.
            # So [y1-1 : y2, x1-1 : x2]
            cropped_data = data[y1-1:y2, x1-1:x2]
            
            # Update WCS
            cropped_wcs = wcs[y1-1:y2, x1-1:x2]
            new_header = cropped_wcs.to_header()
            
            # Update header with necessary keys from old header that might be lost?
            # usually wcs.to_header() keeps the WCS keywords. We might want to keep others.
            # But usually for just flux maps, WCS is the most important.
            
            # Create new HDU
            new_hdu = fits.PrimaryHDU(data=cropped_data, header=new_header)
            
            # Save
            new_hdu.writeto(output_path, overwrite=True)
            print(f"  Saved to {output_path}")

    except Exception as e:
        print(f"  Failed: {e}")

def main():
    # Read CSV
    df = pd.read_csv(CROPPING_CSV)
    
    # Filter for targets
    # Clean up column names just in case (strip whitespace)
    df.columns = [c.strip() for c in df.columns]
    
    for frb in TARGET_FRBS:
        row = df[df['FRB'] == frb]
        if row.empty:
            print(f"Warning: {frb} not found in {CROPPING_CSV}")
            continue
        
        # Get coordinates
        # Handle potential string/int types
        xmin = row.iloc[0]['xmin']
        xmax = row.iloc[0]['xmax']
        ymin = row.iloc[0]['ymin']
        ymax = row.iloc[0]['ymax']
        
        crop_image(frb, xmin, xmax, ymin, ymax)

    print("Done.")

if __name__ == "__main__":
    main()
