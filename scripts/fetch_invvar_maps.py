"""
Fetch inverse variance (invvar) maps for all FRB large cutouts.

For Legacy Survey sources: re-fetches with invvar=True and extracts HDU[1].
For Pan-STARRS sources: fetches the .wt.fits weight map cutout.

Saves as {FRB}_invvar.fits in large_cutouts/, matching the existing flux files.
"""

import os
import time
import requests
import numpy as np
from io import BytesIO
from astropy.io import fits
import pandas as pd

MAX_RETRIES = 3
RETRY_DELAY = 10  # seconds

OUTPUT_DIR = "large_cutouts"

# Legacy Survey parameters (must match original fetch)
LS_CUTOUT_SIZE = 1000
LS_PIXSCALE = 0.262

# Pan-STARRS parameters
PS1_SIZE = 1048


def fetch_legacy_invvar(ra, dec, size=LS_CUTOUT_SIZE, pixscale=LS_PIXSCALE):
    """Fetch invvar map from Legacy Survey by requesting with invvar=True."""
    layers = ['ls-dr10', 'ls-dr9']
    for layer in layers:
        url = (
            f"https://www.legacysurvey.org/viewer/cutout.fits?"
            f"ra={ra}&dec={dec}&size={size}&layer={layer}"
            f"&pixscale={pixscale}&bands=r&invvar"
        )
        for attempt in range(1, MAX_RETRIES + 1):
            print(f"    LS {layer}: fetching with invvar (attempt {attempt})...", end=" ", flush=True)
            try:
                resp = requests.get(url, timeout=300)
                if resp.status_code != 200:
                    print(f"HTTP {resp.status_code}")
                    break  # try next layer
                with fits.open(BytesIO(resp.content)) as hdul:
                    if len(hdul) < 2:
                        print(f"Only {len(hdul)} HDU(s), expected >=2")
                        break
                    invvar_data = hdul[1].data
                    header = hdul[1].header.copy()
                    if invvar_data is None or np.all(invvar_data == 0):
                        print("No invvar coverage (all zeros)")
                        break
                    print(f"OK! shape={invvar_data.shape}")
                    return invvar_data, header, layer
            except requests.exceptions.RequestException as e:
                print(f"Error: {e}")
                if attempt < MAX_RETRIES:
                    print(f"    Retrying in {RETRY_DELAY}s...")
                    time.sleep(RETRY_DELAY)
                continue
    return None, None, None


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
        print(f"    Error querying PS1 filename: {e}")
    return None


def fetch_ps1_invvar(ra, dec, size=PS1_SIZE):
    """Fetch weight map from Pan-STARRS (equivalent to invvar)."""
    if dec < -30:
        print(f"    PS1: Dec={dec:.4f} below coverage (Dec > -30)")
        return None, None

    flux_filename = get_ps1_filename(ra, dec, filters="r")
    if not flux_filename:
        print("    PS1: No filename found")
        return None, None

    # Derive weight filename: .fits -> .wt.fits
    if ".fits" in flux_filename and ".wt.fits" not in flux_filename:
        wt_filename = flux_filename.replace(".fits", ".wt.fits")
    else:
        wt_filename = flux_filename + ".wt"

    print(f"    PS1 weight file: {wt_filename}")

    service = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"
    url = f"{service}?ra={ra}&dec={dec}&size={size}&format=fits&red={wt_filename}"
    print(f"    Downloading PS1 weight cutout...", end=" ", flush=True)
    try:
        r = requests.get(url, timeout=180)
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
            print(f"OK! shape={data.shape}")
            return data, header
    except Exception as e:
        print(f"Error: {e}")
        return None, None


def main():
    # Read coordinates
    frb_coord_file = "master_frb_summary.csv"
    if not os.path.exists(frb_coord_file):
        print(f"Error: {frb_coord_file} not found")
        return

    frbs = pd.read_csv(frb_coord_file)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    total = len(frbs)
    success = 0
    failed = []

    print(f"{'='*60}")
    print(f"  Fetching invvar maps for {total} FRBs")
    print(f"  Output: {OUTPUT_DIR}/")
    print(f"{'='*60}")

    for _, frb in frbs.iterrows():
        frb_name = frb['FRB']
        ra = frb['RA_deg']
        dec = frb['DEC_deg']

        # Check existing flux file to determine source survey
        flux_path = os.path.join(OUTPUT_DIR, f"{frb_name}_flux.fits")
        if not os.path.exists(flux_path):
            print(f"\n[{frb_name}] SKIP - no flux file found")
            failed.append(frb_name)
            continue

        flux_header = fits.getheader(flux_path)
        survey = flux_header.get('SURVEY', '')

        output_path = os.path.join(OUTPUT_DIR, f"{frb_name}_invvar.fits")

        # Skip if already downloaded
        if os.path.exists(output_path):
            print(f"\n[{frb_name}] SKIP - invvar already exists")
            success += 1
            continue

        print(f"\n[{frb_name}] RA={ra:.6f} Dec={dec:.6f} source={survey}")

        if 'Legacy' in survey:
            invvar_data, header, layer = fetch_legacy_invvar(ra, dec)
            if invvar_data is not None:
                header['OBJECT'] = frb_name
                header['SURVEY'] = f'Legacy Survey ({layer})'
                header['HISTORY'] = 'Inverse Variance Map'
                fits.writeto(output_path, invvar_data, header, overwrite=True)
                print(f"    SAVED: {output_path}")
                success += 1
            else:
                print(f"    FAILED: Could not fetch invvar from Legacy Survey")
                failed.append(frb_name)

        elif 'Pan-STARRS' in survey:
            wt_data, header = fetch_ps1_invvar(ra, dec)
            if wt_data is not None:
                header['OBJECT'] = frb_name
                header['SURVEY'] = 'Pan-STARRS DR1'
                header['HISTORY'] = 'Weight (Inverse Variance) Map from Pan-STARRS'
                fits.writeto(output_path, wt_data, header, overwrite=True)
                print(f"    SAVED: {output_path}")
                success += 1
            else:
                print(f"    FAILED: Could not fetch weight map from Pan-STARRS")
                failed.append(frb_name)
        else:
            print(f"    SKIP: Unknown survey '{survey}'")
            failed.append(frb_name)

    print(f"\n{'='*60}")
    print(f"  Done! {success}/{total} invvar maps fetched successfully")
    if failed:
        print(f"  Failed: {', '.join(failed)}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
