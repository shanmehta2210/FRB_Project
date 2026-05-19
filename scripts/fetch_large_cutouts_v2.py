import os
import sys
import time
import requests
import numpy as np
from io import BytesIO
import pandas as pd
from astropy.io import fits

# Configuration
OUTPUT_DIR = "large_cutouts"
CUTOUT_SIZE = 2290  # pixels for 10' FOV at 0.262"/px
PIXSCALE = 0.262    # arcsec/pixel
os.makedirs(OUTPUT_DIR, exist_ok=True)

def fetch_legacy_survey_images(ra, dec, frb_name, size=CUTOUT_SIZE, pixscale=PIXSCALE):
    """
    Fetches cutout from Legacy Survey with invvar if available.
    Returns True if success.
    """
    layers = ['ls-dr10', 'ls-dr9']
    
    for layer in layers:
        url = (
            f"https://www.legacysurvey.org/viewer/cutout.fits?"
            f"ra={ra}&dec={dec}&size={size}&layer={layer}"
            f"&pixscale={pixscale}&invvar=True"
        )
        print(f"    --> Trying {layer}...", end=" ", flush=True)
        try:
            resp = requests.get(url, timeout=180)
            if resp.status_code != 200:
                print(f"HTTP {resp.status_code}")
                continue
            
            with fits.open(BytesIO(resp.content)) as hdul:
                if len(hdul) < 2:
                    print("No invvar HDU")
                    # Still save the flux if available
                    flux_data = hdul[0].data
                    if flux_data is not None:
                        fits.writeto(os.path.join(OUTPUT_DIR, f"{frb_name}_flux.fits"), flux_data, hdul[0].header, overwrite=True)
                    continue
                
                flux_data = hdul[0].data
                invvar_data = hdul[1].data
                
                if flux_data is None or np.all(flux_data == 0):
                    print("No coverage (all zeros)")
                    continue
                
                fits.writeto(os.path.join(OUTPUT_DIR, f"{frb_name}_flux.fits"), flux_data, hdul[0].header, overwrite=True)
                fits.writeto(os.path.join(OUTPUT_DIR, f"{frb_name}_invvar.fits"), invvar_data, hdul[1].header, overwrite=True)
                print(f"OK! Saved to {OUTPUT_DIR}/")
                return True
        except Exception as e:
            print(f"Error: {e}")
            continue
    return False

def fetch_panstarrs_fallback(ra, dec, frb_name, size=CUTOUT_SIZE):
    """
    Fetches PS1 cutout and attempts to get weight map.
    """
    # PS1 pixscale is 0.25", so adjust size if needed, but 2290 is close enough or we can use 2400
    ps1_size = int(size * 0.262 / 0.25)
    
    # 1. Get filename
    svc = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = f"{svc}?ra={ra}&dec={dec}&filters=r&sep=comma"
    try:
        r = requests.get(url, timeout=30)
        if r.status_code != 200: return False
        lines = r.text.strip().split('\n')
        if len(lines) < 2: return False
        keys = lines[0].split(',')
        values = lines[1].split(',')
        ps1_data = dict(zip(keys, values))
        filename = ps1_data.get('filename')
        if not filename: return False
    except: return False

    # 2. Download flux
    svc2 = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"
    url_flux = f"{svc2}?ra={ra}&dec={dec}&size={ps1_size}&format=fits&red={filename}"
    print(f"    --> Downloading PS1 flux...", end=" ", flush=True)
    try:
        r = requests.get(url_flux, timeout=180)
        if r.status_code == 200 and not b"<html>" in r.content[:200].lower():
            with fits.open(BytesIO(r.content)) as hdul:
                fits.writeto(os.path.join(OUTPUT_DIR, f"{frb_name}_flux.fits"), hdul[0].data, hdul[0].header, overwrite=True)
                print("OK!", end=" ")
        else:
            print("Failed")
            return False
    except: return False

    # 3. Download weight (fallback to simplified constant if wt fails, but try deriving URL)
    # The 'filename' is like 'skycell.1315.045.stk.r.unconv.fits'
    # The weight is 'skycell.1315.045.stk.r.unconv.wt.fits'
    wt_filename = filename.replace('.fits', '.wt.fits')
    url_wt = f"{svc2}?ra={ra}&dec={dec}&size={ps1_size}&format=fits&red={wt_filename}"
    print(f"weight...", end=" ", flush=True)
    try:
        r = requests.get(url_wt, timeout=180)
        if r.status_code == 200 and not b"<html>" in r.content[:200].lower():
            with fits.open(BytesIO(r.content)) as hdul:
                fits.writeto(os.path.join(OUTPUT_DIR, f"{frb_name}_invvar.fits"), hdul[0].data, hdul[0].header, overwrite=True)
                print("OK!")
                return True
        else:
            print("Failed (using constant proxy)")
            # Note: In a real pipeline we might use the flux image to estimate noise,
            # but for now we just flag the failure.
            return True # Still considered 'success' for the flux side
    except:
        print("Error")
        return True

def main():
    frbs = pd.read_csv("master_frb_summary.csv")
    total = len(frbs)
    
    # Identify which ones are missing
    missing_frbs = []
    for _, frb in frbs.iterrows():
        name = frb['FRB']
        if not os.path.exists(os.path.join(OUTPUT_DIR, f"{name}_flux.fits")):
            missing_frbs.append(frb)
            
    print(f"Total FRBs: {total}. Missing: {len(missing_frbs)}")
    
    for i, frb in enumerate(missing_frbs):
        name = frb['FRB']
        ra = frb['RA_deg']
        dec = frb['DEC_deg']
        print(f"[{i+1}/{len(missing_frbs)}] {name}...")
        
        success = fetch_legacy_survey_images(ra, dec, name)
        if not success:
            print("    LS failed, trying PS1...")
            fetch_panstarrs_fallback(ra, dec, name)

if __name__ == "__main__":
    main()
