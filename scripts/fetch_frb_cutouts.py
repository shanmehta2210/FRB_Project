
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
CUTOUT_SIZE = 1000  # pixels
PIXSCALE = 0.262    # arcsec/pixel
os.makedirs(OUTPUT_DIR, exist_ok=True)


def progress_bar(current, total, frb_name, status="", width=40):
    """Print a progress bar to terminal."""
    pct = current / total
    filled = int(width * pct)
    bar = "█" * filled + "░" * (width - filled)
    print(f"\r  [{bar}] {current}/{total} ({pct:.0%}) {frb_name} {status}", end="", flush=True)


def fetch_legacy_survey_images(ra, dec, size=CUTOUT_SIZE, pixscale=PIXSCALE):
    """
    Fetches cutout from Legacy Survey.
    Tries bands in order: r to handle missing bands.
    Returns (flux_data, header, layer) or (None, None, None).
    """
    layers = ['ls-dr10', 'ls-dr9']
    bands_to_try = ['r']
    
    for layer in layers:
        for band in bands_to_try:
            url = (
                f"https://www.legacysurvey.org/viewer/cutout.fits?"
                f"ra={ra}&dec={dec}&size={size}&layer={layer}"
                f"&pixscale={pixscale}&bands={band}"
            )
            print(f"    --> Trying {layer} {band}-band...", end=" ", flush=True)
            try:
                resp = requests.get(url, timeout=120)
                if resp.status_code != 200:
                    print(f"HTTP {resp.status_code}")
                    continue
                
                with fits.open(BytesIO(resp.content)) as hdul:
                    data = hdul[0].data
                    header = hdul[0].header.copy()
                    if data is None or np.all(data == 0):
                        print("No coverage (all zeros)")
                        continue
                    print(f"OK! shape={data.shape}")
                    header['BAND'] = band
                    return data, header, layer
            except requests.exceptions.RequestException as e:
                print(f"Error: {e}")
                continue
            
    return None, None, None


def fetch_panstarrs_images(ra, dec, size=1048):
    """
    Fetches cutout from Pan-STARRS DR1.
    PS1 pixel scale = 0.25"/pixel.
    size=1048 gives similar FOV to 1000px at 0.262"/pixel.
    Coverage: Dec > -30.
    Returns (flux_data, header) or (None, None).
    """
    if dec < -30:
        print(f"    --> PS1: Dec={dec:.4f} below coverage (Dec > -30)")
        return None, None

    # Step 1: Get filename
    svc = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = f"{svc}?ra={ra}&dec={dec}&filters=r&sep=comma"
    print(f"    --> Querying PS1 filename...", end=" ", flush=True)
    try:
        r = requests.get(url, timeout=30)
        if r.status_code != 200:
            print(f"HTTP {r.status_code}")
            return None, None
        lines = r.text.strip().split('\n')
        if len(lines) < 2:
            print("No image found")
            return None, None
        keys = lines[0].split(',')
        values = lines[1].split(',')
        ps1_data = dict(zip(keys, values))
        filename = ps1_data.get('filename')
        if not filename:
            print("No filename")
            return None, None
        print(f"found: {filename}")
    except Exception as e:
        print(f"Error: {e}")
        return None, None

    # Step 2: Download cutout
    svc2 = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"
    url2 = f"{svc2}?ra={ra}&dec={dec}&size={size}&format=fits&red={filename}"
    print(f"    --> Downloading PS1 cutout ({size}px)...", end=" ", flush=True)
    try:
        r = requests.get(url2, timeout=180)
        if r.status_code != 200:
            print(f"HTTP {r.status_code}")
            return None, None
        if b"<html>" in r.content[:200].lower():
            print("HTML error response")
            return None, None
        with fits.open(BytesIO(r.content)) as hdul:
            data = hdul[0].data
            header = hdul[0].header.copy()
            if data is None:
                print("No data")
                return None, None
            nan_frac = np.sum(np.isnan(data)) / data.size if data.size > 0 else 0
            print(f"OK! shape={data.shape}, NaN={nan_frac:.1%}")
            return data, header
    except Exception as e:
        print(f"Error: {e}")
        return None, None


def main():
    """
    Main: reads master_frb_summary.csv, fetches cutouts for all FRBs.
    Priority: Legacy Survey DR10 > DR9 > Pan-STARRS DR1
    """
    frb_coord_file = "master_frb_summary.csv"
    if not os.path.exists(frb_coord_file):
        print(f"Error: {frb_coord_file} not found")
        return

    frbs = pd.read_csv(frb_coord_file)
    total = len(frbs)
    
    print(f"=" * 60)
    print(f"  Fetching cutouts for {total} FRBs")
    print(f"  Output: {OUTPUT_DIR}/")
    print(f"  Cutout size: {CUTOUT_SIZE}px at {PIXSCALE}\"/px")
    print(f"  Priority: Legacy Survey DR10 > DR9 > Pan-STARRS DR1")
    print(f"=" * 60)
    
    results = {"success": [], "failed": []}
    start_time = time.time()
    
    for index, frb in frbs.iterrows():
        frb_name = frb['FRB']
        ra = frb['RA_deg']
        dec = frb['DEC_deg']
        
        i = index + 1
        elapsed = time.time() - start_time
        avg_time = elapsed / i if i > 1 else 0
        eta = avg_time * (total - i) if i > 1 else 0
        
        print(f"\n[{i}/{total}] {frb_name}  RA={ra:.6f}  Dec={dec:.6f}  (ETA: {eta:.0f}s)")
        
        # Try Legacy Survey first
        flux_data, header, layer = fetch_legacy_survey_images(ra, dec)
        source = f"Legacy Survey ({layer})" if layer else None
        
        # Fallback to Pan-STARRS
        if flux_data is None:
            print(f"    Legacy Survey failed, trying Pan-STARRS...")
            flux_data, header = fetch_panstarrs_images(ra, dec)
            source = "Pan-STARRS DR1" if flux_data is not None else None
        
        if flux_data is not None and header is not None:
            header['OBJECT'] = frb_name
            header['SURVEY'] = source
            
            output_file = os.path.join(OUTPUT_DIR, f"{frb_name}_flux.fits")
            fits.writeto(output_file, flux_data, header, overwrite=True)
            
            print(f"    [OK] Saved: {output_file} ({source})")
            results["success"].append(frb_name)
        else:
            print(f"    [FAIL] No data from any source")
            results["failed"].append(frb_name)
    
    # Summary
    elapsed = time.time() - start_time
    print(f"\n{'=' * 60}")
    print(f"  DONE in {elapsed:.0f}s")
    print(f"  Success: {len(results['success'])}/{total}")
    if results["failed"]:
        print(f"  Failed:  {', '.join(results['failed'])}")
    print(f"{'=' * 60}")

if __name__ == "__main__":
    main()
