
import os
import requests
import numpy as np
from astropy.io import fits
from io import BytesIO

# Configuration
OUTPUT_DIR = "host_galaxies_data"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Missing FRBs List
MISSING_FRBS = [
    {"FRB": "20220207C", "RA": 310.1995416666667, "DEC": 72.88232777777776},
    {"FRB": "20220825A", "RA": 311.98145833333336, "DEC": 72.58496944444444},
    {"FRB": "20220912A", "RA": 347.2704166666667, "DEC": 48.706944444444446},
    {"FRB": "20220307B", "RA": 350.8744999999999, "DEC": 72.19238611111112},
    {"FRB": "20220319D", "RA": 32.17791666666667, "DEC": 71.03526111111111}
]

def get_ps1_filename(ra, dec, filters="r"):
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = f"{service}?ra={ra}&dec={dec}&filters={filters}&sep=comma"
    try:
        r = requests.get(url, timeout=30)
        if r.status_code == 200:
            lines = r.text.strip().split('\n')
            if len(lines) > 1:
                # The first line is header: "filename,shortname,..."
                keys = lines[0].split(',')
                values = lines[1].split(',')
                data = dict(zip(keys, values))
                return data.get('filename')
    except Exception as e:
        print(f"  Error querying filename: {e}")
    return None

def download_fits(filename, ra, dec, size=256):
    service = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"
    # Note: red={filename} is the parameter to specify the file
    url = f"{service}?ra={ra}&dec={dec}&size={size}&format=fits&red={filename}"
    print(f"  Downloading from: {url}")
    try:
        r = requests.get(url, timeout=60)
        if r.status_code == 200:
            # Check if content is HTML error
            if b"<html>" in r.content[:100]:
                print(f"  Received HTML instead of FITS: {r.content[:200]}")
                return None
            return fits.open(BytesIO(r.content))
        else:
            print(f"  Failed: Status {r.status_code}")
    except Exception as e:
        print(f"  Error: {e}")
    return None

def process_frb(frb):
    name = frb['FRB']
    ra = frb['RA']
    dec = frb['DEC']
    print(f"Processing {name}...")

    # 1. Get Flux Filename
    # Query for 'r' band
    flux_filename = get_ps1_filename(ra, dec, filters="r")
    if not flux_filename:
        print(f"  Failed to find filename for {name}")
        return

    print(f"  Found flux filename: {flux_filename}")
    
    # 2. Derive Weight Filename
    # PS1 pattern: ...stk.r.unconv.fits -> ...stk.r.unconv.wt.fits
    if ".fits" in flux_filename and ".wt.fits" not in flux_filename:
         wt_filename = flux_filename.replace(".fits", ".wt.fits")
    else:
         # Fallback or weird naming
         wt_filename = flux_filename + ".wt" 
         
    print(f"  Derived weight filename: {wt_filename}")

    # 3. Download Flux
    flux_hdul = download_fits(flux_filename, ra, dec)
    
    # 4. Download Weight
    wt_hdul = download_fits(wt_filename, ra, dec)

    if flux_hdul and wt_hdul:
        # Save Flux
        flux_path = os.path.join(OUTPUT_DIR, f"{name}_flux.fits")
        flux_header = flux_hdul[0].header
        flux_header['OBJECT'] = name
        flux_header['SURVEY'] = 'Pan-STARRS'
        flux_hdul.writeto(flux_path, overwrite=True)
        print(f"  Saved flux: {flux_path}")

        # Save Invvar (Weight)
        wt_path = os.path.join(OUTPUT_DIR, f"{name}_invvar.fits")
        wt_data = wt_hdul[0].data
        wt_header = wt_hdul[0].header
        wt_header['OBJECT'] = name
        wt_header['SURVEY'] = 'Pan-STARRS'
        wt_hdul.writeto(wt_path, overwrite=True)
        print(f"  Saved invvar: {wt_path}")

        # Calculate and Save Sigma
        # Sigma = 1 / sqrt(Weight)
        with np.errstate(divide='ignore', invalid='ignore'):
            sigma_data = 1.0 / np.sqrt(wt_data)
        
        # Mask invalids
        sigma_data[wt_data <= 0] = np.nan 
        
        sigma_hdu = fits.PrimaryHDU(data=sigma_data, header=flux_header)
        sigma_path = os.path.join(OUTPUT_DIR, f"{name}_sigma.fits")
        sigma_hdu.writeto(sigma_path, overwrite=True)
        print(f"  Saved sigma: {sigma_path}")

        flux_hdul.close()
        wt_hdul.close()
    else:
        print(f"  Failed to retrieve data for {name}")

if __name__ == "__main__":
    for frb in MISSING_FRBS:
        process_frb(frb)
    print("Done.")
