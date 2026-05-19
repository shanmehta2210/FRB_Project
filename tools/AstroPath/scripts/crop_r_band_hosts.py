import pandas as pd
import numpy as np
import os
from astropy.io import fits
from astropy.wcs import WCS

# Configuration
base_dir = r"C:\Users\lenovo\Desktop\Bhardwajetal_2024_nature_inclination_angle-main"
master_csv = os.path.join(base_dir, "master_frb_summary.csv")
large_dir = os.path.join(base_dir, "large_cutouts")
out_dir = os.path.join(base_dir, "cropped_host_galaxies")

# User provided bounds: frb -> (xmin, xmax), (ymin, ymax)
bounds = {
    "20190611B": {"xmin": 1135, "xmax": 1175, "ymin": 1135, "ymax": 1175},
    "20190711A": {"xmin": 1103, "xmax": 1154, "ymin": 1108, "ymax": 1156},
    "20200430A": {"xmin": 1130, "xmax": 1165, "ymin": 1130, "ymax": 1165},
    "20220105A": {"xmin": 1123, "xmax": 1155, "ymin": 1135, "ymax": 1162},
    "20220725A": {"xmin": 1108, "xmax": 1177, "ymin": 1111, "ymax": 1179},
    "20221106A": {"xmin": 1115, "xmax": 1175, "ymin": 1121, "ymax": 1173},
    "20230526A": {"xmin": 1125, "xmax": 1175, "ymin": 1114, "ymax": 1163},
    "20230708A": {"xmin": 1133, "xmax": 1155, "ymin": 1138, "ymax": 1156},
    "20230902A": {"xmin": 1128, "xmax": 1158, "ymin": 1133, "ymax": 1163},
    "20231226A": {"xmin": 1113, "max": 1171, "ymin": 1116, "ymax": 1164},
    "20240201A": {"xmin": 1073, "xmax": 1189, "ymin": 1126, "ymax": 1199},
    "20240208A": {"xmin": 1130, "xmax": 1169, "ymin": 1120, "ymax": 1159},
    "20240210A": {"xmin": 1064, "xmax": 1292, "ymin": 1027, "ymax": 1228},
    "20240304A": {"xmin": 1133, "xmax": 1165, "ymin": 1135, "ymax": 1167},
    "20240310A": {"xmin": 1114, "xmax": 1160, "ymin": 1128, "ymax": 1170},
    "20240318A": {"xmin": 1114, "xmax": 1167, "ymin": 1116, "ymax": 1180},
}

# Fix typo in dictionary (xmax for 20231226A was missing 'x')
bounds["20231226A"] = {"xmin": 1113, "xmax": 1171, "ymin": 1116, "ymax": 1164}

# Step 1: Update master CSV
print("Updating master CSV...")
df = pd.read_csv(master_csv)
for frb, b in bounds.items():
    mask = df['FRB'] == frb
    if mask.any():
        df.loc[mask, 'xmin'] = b['xmin']
        df.loc[mask, 'xmax'] = b['xmax']
        df.loc[mask, 'ymin'] = b['ymin']
        df.loc[mask, 'ymax'] = b['ymax']
df.to_csv(master_csv, index=False)
print("CSV update complete.")

# Step 2: Create output dir (already exists as per user, but let's be safe)
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# Step 3: Process FITS files
print("Cropping FITS files...")
for frb, b in bounds.items():
    print(f"  Processing {frb}...")
    
    # 1-based to 0-based indexing
    # Fits indexing [ymin:ymax, xmin:xmax]
    y_start, y_end = b['ymin'] - 1, b['ymax']
    x_start, x_end = b['xmin'] - 1, b['xmax']
    
    for suffix in ['_flux', '_invvar']:
        in_path = os.path.join(large_dir, f"{frb}{suffix}.fits")
        if not os.path.exists(in_path):
            print(f"    Warning: Missing {in_path}")
            continue
            
        with fits.open(in_path) as hdul:
            hdr = hdul[0].header
            data = hdul[0].data
            
            # Data is (4, 2290, 2290) -> g, r, i, z
            # Extract Slice 2 (r-band) -> index 1
            if data.ndim == 3:
                r_band_data = data[1, y_start:y_end, x_start:x_end]
            else:
                r_band_data = data[y_start:y_end, x_start:x_end]
                
            # Handle invvar to sigma conversion
            out_suffix = "_flux" if suffix == "_flux" else "_sigma"
            if out_suffix == "_sigma":
                # sigma = 1 / sqrt(invvar)
                # prevent division by zero
                r_band_data = np.where(r_band_data > 0, 1.0 / np.sqrt(r_band_data), 0)
            
            # Create new HDU and update WCS
            new_hdr = hdr.copy()
            # Remove 3D parts of WCS if present
            for key in ['NAXIS3', 'CTYPE3', 'CRPIX3', 'CRVAL3', 'CDELT3', 'CUNIT3']:
                if key in new_hdr: del new_hdr[key]
            
            new_hdr['NAXIS'] = 2
            
            # WCS update: CRPIX shift
            # CRPIX_new = CRPIX_old - x_start
            if 'CRPIX1' in new_hdr: new_hdr['CRPIX1'] -= x_start
            if 'CRPIX2' in new_hdr: new_hdr['CRPIX2'] -= y_start
            
            out_path = os.path.join(out_dir, f"{frb}{out_suffix}.fits")
            fits.writeto(out_path, r_band_data, new_hdr, overwrite=True)

print("Batch processing complete.")
