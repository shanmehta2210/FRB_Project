"""
Fetch Pan-STARRS r-band cutouts for the 3 FRBs that have no Legacy Survey coverage.
Saves directly to large_cutouts/ as {FRB}_flux.fits, matching the project convention.
"""

import os
import requests
import numpy as np
from io import BytesIO
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

OUTPUT_DIR = "large_cutouts"

# The 3 FRBs that failed on Legacy Survey
FAILED_FRBS = {
    "20171020A": (333.75, -19.6667),
    "20210807D": (299.2214, -0.7624),
    "20211127I": (199.8082, -18.8378),
}

# Pan-STARRS cutout size in pixels (240 pixels at 0.25"/pixel gives ~1' FOV)
# But we want to match Legacy Survey's 1000 pixels at 0.262"/pixel ~ 4.37'
# PS1 pixel scale is 0.25"/pixel. 1000 * 0.262 / 0.25 = 1048 pixels.
PS1_SIZE = 1048


def get_ps1_filename(ra, dec, filters="r"):
    """Query PS1 for the image filename covering this position."""
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = f"{service}?ra={ra}&dec={dec}&filters={filters}&sep=comma"
    try:
        r = requests.get(url, timeout=30)
        if r.status_code == 200:
            lines = r.text.strip().split('\n')
            if len(lines) > 1:
                keys = lines[0].split(',')
                values = lines[1].split(',')
                data = dict(zip(keys, values))
                return data.get('filename')
    except Exception as e:
        print(f"  Error querying filename: {e}")
    return None


def download_ps1_cutout(filename, ra, dec, size=PS1_SIZE):
    """Download a FITS cutout from PS1."""
    service = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"
    url = f"{service}?ra={ra}&dec={dec}&size={size}&format=fits&red={filename}"
    print(f"  Downloading from: {url}")
    try:
        r = requests.get(url, timeout=120)
        if r.status_code == 200:
            if b"<html>" in r.content[:100]:
                print(f"  Received HTML error instead of FITS")
                return None
            return fits.open(BytesIO(r.content))
        else:
            print(f"  Failed: HTTP {r.status_code}")
    except Exception as e:
        print(f"  Error: {e}")
    return None


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    success = 0
    failed = []
    
    for frb_name, (ra, dec) in FAILED_FRBS.items():
        print(f"\n{'='*60}")
        print(f"Fetching {frb_name} from Pan-STARRS (RA={ra}, Dec={dec})...")
        
        # Check declination: PS1 covers Dec > -30
        if dec < -30:
            print(f"  WARNING: Dec={dec} is below PS1 coverage (Dec > -30)")
        
        # Get filename
        ps1_file = get_ps1_filename(ra, dec, filters="r")
        if not ps1_file:
            print(f"  FAILED: No PS1 filename found")
            failed.append(frb_name)
            continue
        
        print(f"  Found PS1 file: {ps1_file}")
        
        # Download cutout
        hdul = download_ps1_cutout(ps1_file, ra, dec)
        if hdul is None:
            print(f"  FAILED: Could not download cutout")
            failed.append(frb_name)
            continue
        
        data = hdul[0].data
        header = hdul[0].header
        
        if data is None:
            print(f"  FAILED: No data in FITS")
            failed.append(frb_name)
            continue
        
        # Add metadata
        header['OBJECT'] = frb_name
        header['SURVEY'] = 'Pan-STARRS DR1'
        
        # Save
        output_file = os.path.join(OUTPUT_DIR, f"{frb_name}_flux.fits")
        fits.writeto(output_file, data, header, overwrite=True)
        
        mean, median, std = sigma_clipped_stats(data[~np.isnan(data)], sigma=3.0)
        print(f"  Saved: {output_file}")
        print(f"  Shape: {data.shape}")
        print(f"  Stats: mean={mean:.4f}, median={median:.4f}, std={std:.4f}")
        
        hdul.close()
        success += 1
    
    print(f"\n{'='*60}")
    print(f"Completed: {success}/{len(FAILED_FRBS)} succeeded")
    if failed:
        print(f"Failed: {', '.join(failed)}")


if __name__ == "__main__":
    main()
