import os
import requests
import numpy as np
from io import BytesIO
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

OUTPUT_DIR = "large_cutouts"
CUTOUT_SIZE = 1000  # pixels
PIXSCALE = 0.262    

# 20210807D: 299.2214 -0.7624
# 20220319D: 02:08:42.7, +71:02:06.9 -> compute deg
# 20220825A: 311.9815, 72.585

c_20220319D = SkyCoord('02h08m42.7s +71d02m06.9s', frame='icrs')

FRBS_TO_FETCH = {
    "20210807D": (299.2214, -0.7624),
    "20220319D": (c_20220319D.ra.deg, c_20220319D.dec.deg),
    "20220825A": (311.9815, 72.5850)
}

def fetch_legacy_survey_images(ra, dec, size=CUTOUT_SIZE, pixscale=PIXSCALE):
    layers = ['ls-dr10', 'ls-dr9']
    for layer in layers:
        url = (f"https://www.legacysurvey.org/viewer/cutout.fits?"
               f"ra={ra}&dec={dec}&size={size}&layer={layer}&pixscale={pixscale}&bands=r")
        print(f"    --> Trying {layer}...", end=" ", flush=True)
        try:
            resp = requests.get(url, timeout=120)
            if resp.status_code == 200:
                with fits.open(BytesIO(resp.content)) as hdul:
                    data = hdul[0].data
                    header = hdul[0].header.copy()
                    if data is not None and not np.all(data == 0):
                        print(f"OK! shape={data.shape}")
                        return data, header, layer
            print(f"HTTP {resp.status_code} or empty")
        except Exception as e:
            print(f"Error: {e}")
    return None, None, None

def fetch_panstarrs_images(ra, dec, size=1048):
    if dec < -30:
        return None, None
    svc = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = f"{svc}?ra={ra}&dec={dec}&filters=r&sep=comma"
    print(f"    --> Querying PS1 filename...", end=" ", flush=True)
    try:
        r = requests.get(url, timeout=30)
        if r.status_code != 200:
            return None, None
        lines = r.text.strip().split('\n')
        if len(lines) < 2: return None, None
        ps1_data = dict(zip(lines[0].split(','), lines[1].split(',')))
        filename = ps1_data.get('filename')
        if not filename: return None, None
        print(f"found: {filename}")
        
        url2 = f"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?ra={ra}&dec={dec}&size={size}&format=fits&red={filename}"
        print(f"    --> Downloading PS1 cutout...", end=" ", flush=True)
        r2 = requests.get(url2, timeout=180)
        if r2.status_code == 200 and b"<html>" not in r2.content[:200].lower():
            with fits.open(BytesIO(r2.content)) as hdul:
                print(f"OK! shape={hdul[0].data.shape}")
                return hdul[0].data, hdul[0].header.copy()
    except Exception as e:
        print(f"Error: {e}")
    return None, None

def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for frb, (ra, dec) in FRBS_TO_FETCH.items():
        print(f"\n[{frb}]   RA={ra:.6f}  Dec={dec:.6f}")
        
        data, header, layer = fetch_legacy_survey_images(ra, dec)
        source = f"Legacy Survey ({layer})" if layer else None
        
        if data is None:
            print("    Legacy Survey failed, trying Pan-STARRS...")
            data, header = fetch_panstarrs_images(ra, dec)
            source = "Pan-STARRS DR1" if data is not None else None
            
        if data is not None and header is not None:
            header['OBJECT'] = frb
            header['SURVEY'] = source
            out_path = os.path.join(OUTPUT_DIR, f"{frb}_flux.fits")
            fits.writeto(out_path, data, header, overwrite=True)
            print(f"    [OK] Saved: {out_path} ({source})")
        else:
            print(f"    [FAIL] No data found.")

if __name__ == "__main__":
    main()
