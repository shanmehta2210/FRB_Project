import os
import subprocess
import shutil
import pandas as pd
import numpy as np
from astropy.io import fits
import sys
import gc

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from tools.photutils.scripts.run_photutils_ellipse_nopsf import choose_half_light_isophote
from photutils.isophote import Ellipse

def basic_photutils_fit(frb, flux_path):
    try:
        img = fits.getdata(flux_path).squeeze()
        img = np.nan_to_num(img)
        ell = Ellipse(img)
        iso = ell.fit_image(step=0.5)
        res = choose_half_light_isophote(iso)
        if res is not None:
             return {"FRB": frb, "status": "ok", "sma_rep_pix": res['sma_pix'], "b_over_a": 1.0 - res['eps'], "pa": res['pa_rad']}
        return {"FRB": frb, "status": "failed"}
    except Exception as e:
        return {"FRB": frb, "status": f"failed: {e}"}

def to_wsl(path):
    if path is None or path == 'none': return 'none'
    if ':' in path:
        drive, rest = os.path.splitdrive(path)
        drive_letter = drive[0].lower()
        wsl_path = f"/mnt/{drive_letter}" + rest.replace('\\', '/')
        return wsl_path
    return path.replace('\\', '/')

def parse_fitlog(log_path):
    with open(log_path, 'r') as f: content = f.read()
    blocks = content.split('--------------------------------------------')
    target_block = None
    for block in reversed(blocks):
        if "sersic" in block and "Chi^2/nu =" in block:
            target_block = block
            break
            
    result = {'re': '', 'n': '', 'b_a': '', 'pa': ''}
    if not target_block: return result
    
    lines = target_block.strip().split('\n')
    for line in lines:
        if "sersic" in line and ":" in line:
            clean = line.replace('(', ' ').replace(')', ' ').replace(',', ' ').replace('[', ' ').replace(']', ' ')
            parts = clean.split()
            if len(parts) >= 9 and parts[0].strip() == 'sersic':
                try:
                    result['re'] = float(parts[5].replace('*', ''))
                    result['n'] = float(parts[6].replace('*', ''))
                    result['b_a'] = float(parts[7].replace('*', ''))
                    result['pa'] = float(parts[8].replace('*', ''))
                except:
                    pass
    return result

def run_simulation(base_path="tools/simulation"):
    catalog_path = os.path.join(base_path, "mock_catalog.csv")
    if not os.path.exists(catalog_path):
        print("Mock catalog not found!")
        return

    df = pd.read_csv(catalog_path)
    galfit_results = []
    photutils_results = []
    
    outroot = os.path.join(base_path, "mock_galfit_runs")
    os.makedirs(outroot, exist_ok=True)
    
    for i, row in df.iterrows():
        frb = row['FRB']
        flux_path = os.path.abspath(row['flux_path'])
        sigma_path = os.path.abspath(row['sigma_path'])
        psf_path = os.path.abspath(os.path.join(base_path, "mock_data", "mock_psf.fits"))
        
        # --- PHOTUTILS EXECUTION ---
        print(f"[{i+1}/{len(df)}] Running Photutils on {frb}...")
        phot_res = basic_photutils_fit(frb, flux_path)
        photutils_results.append(phot_res)
            
        # --- GALFIT SETUP AND RECOVERY ---
        print(f"[{i+1}/{len(df)}] Running GALFIT on {frb}...")
        frb_dir = os.path.join(outroot, frb)
        os.makedirs(frb_dir, exist_ok=True)
        
        # Set guesses very close to truth to observe ideal convergence behavior
        init_re = row['true_re_pix'] * 1.01
        init_n = row['true_n'] * 1.01
        init_q = np.clip(row['true_q'] * 1.01, 0.1, 1.0)
        
        # Copy files to local dir to avoid GALFIT's Fortran path length limits
        shutil.copy(flux_path, os.path.join(frb_dir, "flux.fits"))
        shutil.copy(sigma_path, os.path.join(frb_dir, "sigma.fits"))
        shutil.copy(psf_path, os.path.join(frb_dir, "psf.fits"))
        
        feedme = os.path.join(frb_dir, 'galfit.feedme')
        with open(feedme, 'w') as f:
            f.write(f'''===============================================================================
# IMAGE and GALFIT CONTROL PARAMETERS
A) flux.fits
B) out.fits
C) none
D) psf.fits
E) 1
F) none
G) none
H) 1 100 1 100
I) 120 120
J) 25.0
K) 0.262 0.262
O) regular
P) 0

# INITIAL FITTING PARAMETERS
# Component number: 1
0) sersic
1) 50.0 50.0 1 1
3) 12.5 1
4) {init_re:.2f} 1
5) {init_n:.2f} 1
6) 0.0000 0
7) 0.0000 0
8) 0.0000 0
9) {init_q:.2f} 1
10) 45.0 1
Z) 0

# Component number: 2
0) sky
1) 0.0000 1
2) 0 0
3) 0 0
Z) 0
================================================================================
''')
        
        # Run GALFIT
        saved_dir = os.getcwd()
        os.chdir(frb_dir)
        subprocess.run(["wsl", "galfit", "galfit.feedme"], capture_output=True, timeout=120)
        os.chdir(saved_dir)
        
        # Parse logs
        fitlog = os.path.join(frb_dir, "fit.log")
        if os.path.exists(fitlog):
            gres = parse_fitlog(fitlog)
            gres['FRB'] = frb
            gres['status'] = 'ok'
            galfit_results.append(gres)
        else:
            galfit_results.append({'FRB': frb, 'status': 'failed'})

        gc.collect()

    pd.DataFrame(photutils_results).to_csv(os.path.join(base_path, "mock_photutils_results.csv"), index=False)
    pd.DataFrame(galfit_results).to_csv(os.path.join(base_path, "mock_galfit_results.csv"), index=False)
    print("Done executing pipelines on simulated data.")

if __name__ == "__main__":
    run_simulation()
