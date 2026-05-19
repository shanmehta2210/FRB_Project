import os
import sys
import csv
import shutil
import subprocess
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

# Final Orchestrator including both no-PSF and with-PSF runs

def to_wsl(path):
    if not path or path == 'none': return 'none'
    if ':' in path:
        drive, rest = os.path.splitdrive(path)
        return f"/mnt/{drive[0].lower()}" + rest.replace('\\', '/')
    return path.replace('\\', '/')

def get_platescale(wcs):
    if not wcs.is_celestial: wcs = wcs.celestial
    return np.mean(np.sum(wcs.pixel_scale_matrix**2, axis=0)**0.5) * 3600

def parse_fitlog(log_path):
    if not os.path.exists(log_path): return None
    with open(log_path, 'r') as f:
        content = f.read()
    blocks = content.split('--------------------------------------------')
    target_block = None
    for block in reversed(blocks):
        if "sersic" in block and "Chi^2/nu =" in block:
            target_block = block
            break
    if not target_block: return None
    
    res = {'chi2nu': ''}
    lines = target_block.strip().split('\n')
    for i, line in enumerate(lines):
        if "Chi^2/nu =" in line:
            try: res['chi2nu'] = line.split('=')[1].split(',')[0].strip()
            except: pass
        if "sersic" in line and ":" in line:
            clean = line.replace('(', ' ').replace(')', ' ').replace(',', ' ').replace('[', ' ').replace(']', ' ')
            p = clean.split()
            if len(p) >= 9 and p[0].strip() == 'sersic':
                res.update({'x': p[2], 'y': p[3], 'mag': p[4], 're': p[5], 'n': p[6], 'b_a': p[7], 'pa': p[8]})
                if i + 1 < len(lines):
                    err_p = lines[i+1].replace('(', ' ').replace(')', ' ').replace(',', ' ').split()
                    if len(err_p) >= 7:
                        res.update({'x_err': err_p[0], 'y_err': err_p[1], 'mag_err': err_p[2], 're_err': err_p[3], 'n_err': err_p[4], 'b_a_err': err_p[5], 'pa_err': err_p[6]})
    return res

# --- Main Configuration ---
base_dir = r"C:\Users\lenovo\Desktop\Bhardwajetal_2024_nature_inclination_angle-main"
crop_dir = os.path.join(base_dir, "cropped_host_galaxies")
psf_base_dir = os.path.join(base_dir, "psfs", "PSFEX + SExtractor", "final_center_psfs")
runs_dir = os.path.join(base_dir, "tools", "galfit", "runs")

targets = [
    "20190611B", "20190711A", "20200430A", "20220105A", "20220725A", "20221106A", 
    "20230526A", "20230708A", "20230902A", "20231226A", "20240201A", "20240208A", 
    "20240210A", "20240304A", "20240310A", "20240318A"
]

# Import get_initial_guesses logic
sys.path.append(os.path.join(base_dir, "scripts"))
try:
    from get_initial_guesses import get_initial_guesses
except ImportError:
    print("Could not import get_initial_guesses. Ensure scripts/get_initial_guesses.py exists.")
    sys.exit(1)

all_results = []
magzero = 22.5

for frb in targets:
    print(f"\n--- Processing {frb} ---")
    img_path = os.path.join(crop_dir, f"{frb}_flux.fits")
    sigma_path = os.path.join(crop_dir, f"{frb}_sigma.fits")
    psf_src = os.path.join(psf_base_dir, frb, f"{frb}_center_psf_25.fits")
    
    if not all(os.path.exists(p) for p in [img_path, sigma_path, psf_src]):
        print(f"Skipping {frb}: Missing files.")
        continue

    # 1. Get initial guesses
    try:
        guesses = get_initial_guesses(img_path)
        mag_val = guesses['mag'] - 2.5 # Shift to 22.5 from 25.0 system
    except Exception as e:
        print(f"Guess error for {frb}: {e}")
        continue

    # 2. Setup folder and copy files
    run_modes = ["no_psf_sigma", "with_psf_sigma"]
    for mode in run_modes:
        d = os.path.join(runs_dir, frb, mode)
        os.makedirs(d, exist_ok=True)
        shutil.copy(sigma_path, os.path.join(d, "sigma.fits"))
        if "with_psf" in mode:
            shutil.copy(psf_src, os.path.join(d, "psf.fits"))

        # Create constraints (using provided 1.5 to 100 fixed range)
        with open(os.path.join(d, "constraints.txt"), "w") as f:
            f.write("1 n 0.5 to 6.0\n1 re 1.5 to 100.0\n")

        # Generate feedme
        img_data, hdr = fits.getdata(img_path, header=True)
        ny, nx = img_data.shape
        pscale = get_platescale(WCS(hdr))
        psf_val = "psf.fits" if "with_psf" in mode else "none"
        
        feedme_txt = f"""===============================================================================
# IMAGE and GALFIT CONTROL PARAMETERS
A) {to_wsl(img_path)}  # Input data image (FITS file)
B) out.fits  # Output data image block
C) sigma.fits  # Sigma image name
D) {psf_val}  # Input PSF file
E) 1  # PSF fine sampling factor
F) none  # Bad pixel mask
G) constraints.txt  # File with parameter constraints
H) 1 {nx} 1 {ny}  # Image region to fit (1-based full image)
I) {nx+25} {ny+25}  # Size of convolution box
J) {magzero:.4f}  # Photometric zeropoint
K) {pscale:.4f} {pscale:.4f}  # Plate scale [arcsec/pixel]
O) regular  # Display type
P) 0  # Choose: 0=optimize

# INITIAL FITTING PARAMETERS
# Component number: 1
0) sersic  # Component type
1) {guesses['x']:.4f} {guesses['y']:.4f} 1 1  # position x y
3) {mag_val:.4f} 1  # Integrated magnitude
4) {guesses['r_e']:.4f} 1  # effective radius (pix)
5) {guesses['n']:.4f} 1  # sersic index
6) 0.0000 0  # ----
7) 0.0000 0  # ----
8) 0.0000 0  # ----
9) {guesses['axis_ratio']:.4f} 1  # Axis ratio (b/a)
10) {guesses['pa']:.4f} 1  # Position angle
Z) 0  # Skip this model

# Component number: 2
0) sky  # component type
1) 0.0000 1  # Sky background
2) 0 0  # dsky/dx
3) 0 0  # dsky/dy
Z) 0  # Skip this model
================================================================================
"""
        with open(os.path.join(d, "galfit.feedme"), "w") as f: f.write(feedme_txt)

    # 3. Execute GALFIT for BOTH modes
    for mode in run_modes:
        print(f"  Running GALFIT for {frb} ({mode})...")
        d = os.path.join(runs_dir, frb, mode)
        prev_cwd = os.getcwd()
        os.chdir(d)
        subprocess.run(["wsl", "galfit", "galfit.feedme"], capture_output=True, timeout=180)
        os.chdir(prev_cwd)

    # 4. Parse results from BOTH modes
    res_no = parse_fitlog(os.path.join(runs_dir, frb, "no_psf_sigma", "fit.log"))
    res_psf = parse_fitlog(os.path.join(runs_dir, frb, "with_psf_sigma", "fit.log"))
    
    if res_no or res_psf:
        combined = {'FRB': frb}
        if res_no:
            for k, v in res_no.items(): combined[f"{k}_nopsf"] = v
        if res_psf:
            for k, v in res_psf.items(): combined[f"{k}_psf"] = v
        all_results.append(combined)

# Save final CSV with dual-run results
if all_results:
    keys = ['FRB']
    params = ['chi2nu', 'mag', 'mag_err', 're', 're_err', 'n', 'n_err', 'b_a', 'b_a_err', 'pa', 'pa_err', 'x', 'x_err', 'y', 'y_err']
    f_keys = keys + [f"{p}_nopsf" for p in params] + [f"{p}_psf" for p in params]
    
    out_csv = os.path.join(base_dir, "new_16_frbs_galfit_results.csv")
    with open(out_csv, 'w', newline='') as f:
        # Filter keys to only those that exist in the dictionary
        writer = csv.DictWriter(f, fieldnames=f_keys, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(all_results)
    print(f"\nFinalized results for {len(all_results)} targets in {out_csv}")
