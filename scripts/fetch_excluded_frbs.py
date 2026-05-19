"""
Fetch r-band cutouts for the 6 previously excluded FRBs.
Uses coordinates from the FRB Inclination angle estimate sheet.
Tries Legacy Survey DR10 first, falls back to Pan-STARRS DR1.

Output: large_cutouts/{FRB}_flux.fits (1000x1000px for LS, ~1048px for PS1)
"""

import os
import sys
import requests
import numpy as np
from io import BytesIO
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

OUTPUT_DIR = "large_cutouts"

# Coordinates from the Excel sheet
FRBS = {
    "20171020A": (333.8531, -19.5853),
    "20210320C": (204.458854, -16.122552),
    "20210807D": (310.19767916805, 72.8826580623927),
    "20211127I": (199.808176085586, -18.8380065058089),
    "20211203C": (204.565260359076, -31.3791066801824),
    "20211212A": (157.350941401671, 1.36088942102641),
}


def fetch_legacy_survey(ra, dec, size=1000, pixscale=0.262):
    """Try Legacy Survey layers. Returns (data, header) or (None, None)."""
    layers = ['ls-dr10', 'ls-dr9']
    for layer in layers:
        url = (
            f"https://www.legacysurvey.org/viewer/cutout.fits?"
            f"ra={ra}&dec={dec}&size={size}&layer={layer}"
            f"&pixscale={pixscale}&bands=r"
        )
        print(f"    LS {layer}: {url}")
        try:
            resp = requests.get(url, timeout=120)
            if resp.status_code != 200:
                print(f"    HTTP {resp.status_code}")
                continue
            with fits.open(BytesIO(resp.content)) as hdul:
                data = hdul[0].data
                header = hdul[0].header.copy()
                if data is None or np.all(data == 0):
                    print(f"    No coverage (all zeros)")
                    continue
                print(f"    Success! Shape={data.shape}")
                header['SURVEY'] = f'Legacy Survey ({layer})'
                return data, header
        except Exception as e:
            print(f"    Error: {e}")
            continue
    return None, None


def fetch_panstarrs(ra, dec, size=1048):
    """
    Fetch from Pan-STARRS DR1. Coverage: Dec > -30.
    PS1 pixel scale = 0.25"/pixel.
    size=1048 pixels ~ same FOV as 1000px at 0.262"/pixel.
    Returns (data, header) or (None, None).
    """
    if dec < -30:
        print(f"    PS1: Dec={dec:.4f} is below coverage limit (Dec > -30)")
        return None, None

    # Step 1: Get the image filename
    svc = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = f"{svc}?ra={ra}&dec={dec}&filters=r&sep=comma"
    print(f"    PS1 filename query: {url}")
    try:
        r = requests.get(url, timeout=30)
        if r.status_code != 200:
            print(f"    HTTP {r.status_code}")
            return None, None
        lines = r.text.strip().split('\n')
        if len(lines) < 2:
            print(f"    No PS1 image found at this position")
            return None, None
        keys = lines[0].split(',')
        values = lines[1].split(',')
        ps1_data = dict(zip(keys, values))
        filename = ps1_data.get('filename')
        if not filename:
            print(f"    No filename in response")
            return None, None
        print(f"    PS1 file: {filename}")
    except Exception as e:
        print(f"    Filename query error: {e}")
        return None, None

    # Step 2: Download the cutout
    svc2 = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"
    url2 = f"{svc2}?ra={ra}&dec={dec}&size={size}&format=fits&red={filename}"
    print(f"    PS1 cutout: {url2}")
    try:
        r = requests.get(url2, timeout=180)
        if r.status_code != 200:
            print(f"    HTTP {r.status_code}")
            return None, None
        # Check for HTML error response
        if b"<html>" in r.content[:200].lower():
            print(f"    Received HTML error instead of FITS")
            return None, None
        with fits.open(BytesIO(r.content)) as hdul:
            data = hdul[0].data
            header = hdul[0].header.copy()
            if data is None:
                print(f"    No data in FITS")
                return None, None
            # Handle NaN values
            nan_frac = np.sum(np.isnan(data)) / data.size
            print(f"    Success! Shape={data.shape}, NaN fraction={nan_frac:.2%}")
            header['SURVEY'] = 'Pan-STARRS DR1'
            return data, header
    except Exception as e:
        print(f"    Cutout download error: {e}")
        return None, None


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    results = {}
    for frb_name, (ra, dec) in FRBS.items():
        print(f"\n{'='*60}")
        print(f"[{frb_name}] RA={ra:.6f}, Dec={dec:.6f}")

        # Try Legacy Survey first
        print(f"  Trying Legacy Survey...")
        data, header = fetch_legacy_survey(ra, dec)

        # Fallback to Pan-STARRS
        if data is None:
            print(f"  Legacy Survey failed. Trying Pan-STARRS...")
            data, header = fetch_panstarrs(ra, dec)

        if data is not None:
            header['OBJECT'] = frb_name
            output_file = os.path.join(OUTPUT_DIR, f"{frb_name}_flux.fits")
            fits.writeto(output_file, data, header, overwrite=True)

            valid = data[~np.isnan(data)] if np.any(np.isnan(data)) else data
            mean, median, std = sigma_clipped_stats(valid, sigma=3.0)
            print(f"  SAVED: {output_file}")
            print(f"  Shape: {data.shape}, Mean={mean:.6f}, Median={median:.6f}, Std={std:.6f}")
            results[frb_name] = "SUCCESS"
        else:
            print(f"  FAILED: No data from any source")
            results[frb_name] = "FAILED"

    print(f"\n{'='*60}")
    print("SUMMARY:")
    for frb, status in results.items():
        print(f"  {frb}: {status}")
    print(f"Succeeded: {sum(1 for s in results.values() if s == 'SUCCESS')}/{len(FRBS)}")


if __name__ == "__main__":
    main()
