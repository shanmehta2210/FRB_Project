"""
Full GALFIT pipeline with proper sigma maps from inverse variance images.

Steps:
1. Convert large_cutouts invvar -> sigma, crop to match cropped_host_galaxies
2. Setup & run GALFIT Stage 1 (no-PSF) with sigma in C) parameter
3. Setup & run GALFIT Stage 2 (with-PSF) using Stage 1 results as initial guesses
4. Compile results to galfit_sigma_metrics_summary.csv

Constraints: n in [0.5, 6.0], Re in [0.5*FWHM, 500] — parameters are FREE for fitting.
"""

import os
import sys
import csv
import shutil
import subprocess
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


# ---------- path helpers ----------
def project_root():
    return os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))

def to_wsl(path):
    if path is None or path == 'none':
        return 'none'
    if ':' in path:
        drive, rest = os.path.splitdrive(path)
        drive_letter = drive[0].lower()
        return f"/mnt/{drive_letter}" + rest.replace('\\', '/')
    return path.replace('\\', '/')

def get_platescale(wcs):
    if not wcs.is_celestial:
        wcs = wcs.celestial
    return np.mean(np.sum(wcs.pixel_scale_matrix**2, axis=0)**0.5) * 3600  # arcsec


# ================================================================
# STEP 1: Convert invvar -> sigma, crop, and save
# ================================================================
def step1_create_cropped_sigma(root):
    print("\n" + "="*70)
    print("STEP 1: Creating cropped sigma maps from inverse variance images")
    print("="*70)

    large_dir = os.path.join(root, "large_cutouts")
    crop_dir = os.path.join(root, "cropped_host_galaxies")
    csv_path = os.path.join(root, "master_frb_summary.csv")

    # Read crop coordinates
    crops = {}
    with open(csv_path, 'r') as f:
        for row in csv.DictReader(f):
            frb = row['FRB']
            if row['xmin'] and row['xmax'] and row['ymin'] and row['ymax']:
                crops[frb] = {
                    'xmin': int(float(row['xmin'])),
                    'xmax': int(float(row['xmax'])),
                    'ymin': int(float(row['ymin'])),
                    'ymax': int(float(row['ymax']))
                }

    created = 0
    for frb, c in crops.items():
        invvar_path = os.path.join(large_dir, f"{frb}_invvar.fits")
        sigma_out = os.path.join(crop_dir, f"{frb}_sigma.fits")

        if not os.path.exists(invvar_path):
            print(f"  [SKIP] {frb}: invvar not found")
            continue

        with fits.open(invvar_path) as hdu:
            invvar = hdu[0].data.astype(np.float64)
            header = hdu[0].header

        # Open flux image to check units (no scaling needed for nanomaggies)
        flux_path = os.path.join(crop_dir, f"{frb}_flux.fits")
        
        scalar = 1.0

        # sigma = scalar / sqrt(invvar), handle zeros and negatives
        valid = invvar > 0
        sigma = np.zeros_like(invvar)
        sigma[valid] = scalar / np.sqrt(invvar[valid])
        # For invalid pixels (invvar <= 0 or NaN), use a large sigma to downweight
        sigma[~valid] = 1e10
        sigma[np.isnan(sigma)] = 1e10

        # Crop with same DS9 convention as crop_images.py
        x1 = c['xmin'] - 1
        x2 = c['xmax']
        y1 = c['ymin'] - 1
        y2 = c['ymax']

        if x1 < 0 or y1 < 0 or x2 > sigma.shape[1] or y2 > sigma.shape[0]:
            print(f"  [ERROR] {frb}: crop bounds out of range")
            continue

        cropped_sigma = sigma[y1:y2, x1:x2]

        # Update WCS CRPIX
        if 'CRPIX1' in header:
            header['CRPIX1'] -= x1
        if 'CRPIX2' in header:
            header['CRPIX2'] -= y1

        fits.writeto(sigma_out, cropped_sigma.astype(np.float32), header, overwrite=True)
        created += 1
        print(f"  [OK] {frb}: sigma map {cropped_sigma.shape} saved")

    print(f"\nCreated {created} cropped sigma maps in {crop_dir}")
    return created


# ================================================================
# STEP 2: Setup GALFIT Stage 1 (no-PSF) with sigma
# ================================================================
def step2_setup_stage1(root):
    print("\n" + "="*70)
    print("STEP 2: Setting up GALFIT Stage 1 (no-PSF) with sigma maps")
    print("="*70)

    img_dir = os.path.join(root, "cropped_host_galaxies")
    runs_dir = os.path.join(root, "Galfit", "runs")
    guesses_csv = os.path.join(root, "csv_archive/initial_guesses.csv")
    fwhm_csv = os.path.join(root, "psf_fwhm_summary.csv")

    # Read FWHMs
    fwhm_dict = {}
    if os.path.exists(fwhm_csv):
        with open(fwhm_csv, 'r') as f:
            for row in csv.DictReader(f):
                try:
                    fwhm_dict[row['FRB']] = float(row['Avg_FWHM'])
                except (KeyError, ValueError):
                    pass

    # Read initial guesses
    with open(guesses_csv, 'r') as f:
        rows = list(csv.DictReader(f))

    setup_count = 0
    for row in rows:
        filename = row['filename']
        frb = filename.replace("_flux.fits", "")
        imgfile = os.path.join(img_dir, filename)
        sigma_file = os.path.join(img_dir, f"{frb}_sigma.fits")

        if not os.path.exists(imgfile):
            print(f"  [SKIP] {frb}: flux image not found")
            continue
        if not os.path.exists(sigma_file):
            print(f"  [SKIP] {frb}: sigma map not found")
            continue

        mag = float(row['mag'])
        r_e = float(row['r_e'])
        n_idx = float(row['n'])
        axis_ratio = float(row['axis_ratio'])
        pa = float(row['pa'])
        
        # Enforce n=1 lock for all except 20171020A
        if frb != "20171020A":
            n_idx = 1.0
            n_fit_flag = 0  # 0 means parameter is locked
        else:
            n_fit_flag = 1  # 1 means free to fit
        cx = float(row['x'])
        cy = float(row['y'])
        fwhm_val = fwhm_dict.get(frb, 3.0)

        # Create run directory
        no_psf_dir = os.path.join(runs_dir, frb, "no_psf_sigma")
        os.makedirs(no_psf_dir, exist_ok=True)

        # Copy sigma map locally to avoid WSL long-path issues
        local_sigma = os.path.join(no_psf_dir, "sigma.fits")
        shutil.copy(sigma_file, local_sigma)

        # Write constraints
        with open(os.path.join(no_psf_dir, 'constraints.txt'), 'w') as f:
            f.write("1 n 0.5 to 6.0\n")
            f.write(f"1 re {0.5 * fwhm_val:.2f} to 500.0\n")

        # Get image shape and platescale
        img_data, hdr = fits.getdata(imgfile, header=True)
        ny, nx = img_data.shape
        region = (0, nx - 1, 0, ny - 1)
        conv_x = nx + 15
        conv_y = ny + 15
        platescale = get_platescale(WCS(hdr))
        zeropoint = 25.0

        img_wsl = to_wsl(imgfile)

        feedme_path = os.path.join(no_psf_dir, 'galfit.feedme')
        with open(feedme_path, 'w') as f:
            f.write("===============================================================================\n")
            f.write("# IMAGE and GALFIT CONTROL PARAMETERS\n")
            f.write(f"A) {img_wsl}  # Input data image (FITS file)\n")
            f.write(f"B) out.fits  # Output data image block\n")
            f.write(f"C) sigma.fits  # Sigma image name\n")
            f.write(f"D) none  # Input PSF file\n")
            f.write(f"E) 1  # PSF fine sampling factor\n")
            f.write(f"F) none  # Bad pixel mask\n")
            f.write(f"G) constraints.txt  # File with parameter constraints\n")
            f.write(f"H) {region[0]} {region[1]} {region[2]} {region[3]}  # Image region to fit\n")
            f.write(f"I) {conv_x} {conv_y}  # Size of convolution box\n")
            f.write(f"J) {zeropoint:.4f}  # Photometric zeropoint\n")
            f.write(f"K) {platescale:.4f} {platescale:.4f}  # Plate scale [arcsec/pixel]\n")
            f.write(f"O) regular  # Display type\n")
            f.write(f"P) 0  # Choose: 0=optimize\n")
            f.write("\n# INITIAL FITTING PARAMETERS\n")
            f.write("# Component number: 1\n")
            f.write("0) sersic  # Component type\n")
            f.write(f"1) {cx:.4f} {cy:.4f} 1 1  # position x y\n")
            f.write(f"3) {mag:.4f} 1  # Integrated magnitude\n")
            f.write(f"4) {r_e:.4f} 1  # effective radius (pix)\n")
            f.write(f"5) {n_idx:.4f} {n_fit_flag}  # sersic index\n")
            f.write("6) 0.0000 0  # ----\n")
            f.write("7) 0.0000 0  # ----\n")
            f.write("8) 0.0000 0  # ----\n")
            f.write(f"9) {axis_ratio:.4f} 1  # Axis ratio (b/a)\n")
            f.write(f"10) {pa:.4f} 1  # Position angle\n")
            f.write("Z) 0  # Skip this model\n\n")
            f.write("# Component number: 2\n")
            f.write("0) sky  # component type\n")
            f.write("1) 0.0000 1  # Sky background\n")
            f.write("2) 0 0  # dsky/dx\n")
            f.write("3) 0 0  # dsky/dy\n")
            f.write("Z) 0  # Skip this model\n")
            f.write("================================================================================\n")

        setup_count += 1
        print(f"  [OK] {frb}: feedme written to {no_psf_dir}")

    print(f"\nSet up {setup_count} Stage 1 (no-PSF + sigma) runs")
    return setup_count


# ================================================================
# STEP 3: Run GALFIT Stage 1
# ================================================================
def step3_run_stage1(root):
    print("\n" + "="*70)
    print("STEP 3: Running GALFIT Stage 1 (no-PSF + sigma)")
    print("="*70)

    runs_dir = os.path.join(root, "Galfit", "runs")
    frbs = sorted([d for d in os.listdir(runs_dir) if os.path.isdir(os.path.join(runs_dir, d))])

    success = 0
    failed = []
    for frb in frbs:
        no_psf_dir = os.path.join(runs_dir, frb, "no_psf_sigma")
        feedme = os.path.join(no_psf_dir, "galfit.feedme")

        if not os.path.exists(feedme):
            continue

        print(f"  Running GALFIT for {frb} (no_psf_sigma)...", end=" ", flush=True)
        saved_dir = os.getcwd()
        os.chdir(no_psf_dir)

        result = subprocess.run(["wsl", "galfit", "galfit.feedme"],
                                capture_output=True, text=True, timeout=120)
        os.chdir(saved_dir)

        fitlog = os.path.join(no_psf_dir, "fit.log")
        if os.path.exists(fitlog):
            print("OK")
            success += 1
        else:
            print("FAILED")
            failed.append(frb)
            if result.stderr:
                print(f"    stderr: {result.stderr[:200]}")

    print(f"\nStage 1 complete: {success} succeeded, {len(failed)} failed")
    if failed:
        print(f"  Failed: {', '.join(failed)}")
    return success


# ================================================================
# STEP 4: Parse fit.log from Stage 1
# ================================================================
def parse_fitlog(log_path):
    with open(log_path, 'r') as f:
        content = f.read()

    blocks = content.split('--------------------------------------------')

    target_block = None
    for block in reversed(blocks):
        if "sersic" in block and "Chi^2/nu =" in block:
            target_block = block
            break

    result = {'chi2nu': '', 'x': '', 'y': '', 'mag': '', 're': '', 'n': '', 'b_a': '', 'pa': '',
              'x_err': '', 'y_err': '', 'mag_err': '', 're_err': '', 'n_err': '', 'b_a_err': '', 'pa_err': ''}

    if not target_block:
        return result

    lines = target_block.strip().split('\n')
    for i, line in enumerate(lines):
        if "Chi^2/nu =" in line:
            try:
                parts = line.split('=')[1].split(',')
                result['chi2nu'] = float(parts[0].strip())
            except (IndexError, ValueError):
                pass

        if "sersic" in line and ":" in line:
            clean = line.replace('(', ' ').replace(')', ' ').replace(',', ' ').replace('[', ' ').replace(']', ' ')
            parts = clean.split()
            if len(parts) >= 9 and parts[0].strip() == 'sersic':
                try:
                    result['x'] = parts[2].replace('*', '')
                    result['y'] = parts[3].replace('*', '')
                    result['mag'] = parts[4].replace('*', '')
                    result['re'] = parts[5].replace('*', '')
                    result['n'] = parts[6].replace('*', '')
                    result['b_a'] = parts[7].replace('*', '')
                    result['pa'] = parts[8].replace('*', '')
                except (IndexError, ValueError):
                    pass

                # Error line follows
                if i + 1 < len(lines):
                    err_line = lines[i + 1].replace('(', ' ').replace(')', ' ').replace(',', ' ')
                    err_parts = err_line.split()
                    if len(err_parts) >= 7:
                        try:
                            result['x_err'] = err_parts[0]
                            result['y_err'] = err_parts[1]
                            result['mag_err'] = err_parts[2]
                            result['re_err'] = err_parts[3]
                            result['n_err'] = err_parts[4]
                            result['b_a_err'] = err_parts[5]
                            result['pa_err'] = err_parts[6]
                        except IndexError:
                            pass
    return result


# ================================================================
# STEP 5: Setup GALFIT Stage 2 (with-PSF) using Stage 1 results
# ================================================================
def step5_setup_stage2(root):
    print("\n" + "="*70)
    print("STEP 5: Setting up GALFIT Stage 2 (with-PSF + sigma)")
    print("="*70)

    runs_dir = os.path.join(root, "Galfit", "runs")
    psf_dir = os.path.join(root, "psfs", "downsampled_psfs")
    crop_dir = os.path.join(root, "cropped_host_galaxies")

    frbs = sorted([d for d in os.listdir(runs_dir) if os.path.isdir(os.path.join(runs_dir, d))])

    setup_count = 0
    for frb in frbs:
        no_psf_dir = os.path.join(runs_dir, frb, "no_psf_sigma")
        with_psf_dir = os.path.join(runs_dir, frb, "with_psf_sigma")
        fitlog_path = os.path.join(no_psf_dir, "fit.log")
        feedme_src = os.path.join(no_psf_dir, "galfit.feedme")

        if not os.path.exists(fitlog_path) or not os.path.exists(feedme_src):
            continue

        # Parse Stage 1 results
        data = parse_fitlog(fitlog_path)
        if not data['x'] or not data['y']:
            print(f"  [SKIP] {frb}: could not parse Stage 1 fit.log")
            continue

        pos_x = float(data['x'])
        pos_y = float(data['y'])
        mag_val = float(data['mag'])
        re_val = float(data['re'])
        n_val = float(data['n'])
        b_a_val = float(data['b_a'])
        pa_val = float(data['pa'])

        # Enforce n=1 lock for all except 20171020A
        if frb != "20171020A":
            n_val = 1.0
            n_fit_flag = 0
        else:
            n_fit_flag = 1

        # Check PSF exists
        psf_file = os.path.join(psf_dir, f"{frb}_1x_psf.fits")
        if not os.path.exists(psf_file):
            print(f"  [SKIP] {frb}: PSF not found")
            continue

        # Check sigma exists
        sigma_src = os.path.join(crop_dir, f"{frb}_sigma.fits")
        if not os.path.exists(sigma_src):
            print(f"  [SKIP] {frb}: sigma map not found")
            continue

        os.makedirs(with_psf_dir, exist_ok=True)

        # Copy PSF and sigma locally to avoid WSL long-path buffer overflow
        shutil.copy(psf_file, os.path.join(with_psf_dir, "psf.fits"))
        shutil.copy(sigma_src, os.path.join(with_psf_dir, "sigma.fits"))

        # Copy constraints from Stage 1
        constraints_src = os.path.join(no_psf_dir, 'constraints.txt')
        if os.path.exists(constraints_src):
            shutil.copy(constraints_src, os.path.join(with_psf_dir, 'constraints.txt'))

        # Read Stage 1 feedme and modify for Stage 2
        with open(feedme_src, 'r') as f:
            lines = f.readlines()

        feedme_out = os.path.join(with_psf_dir, 'galfit.feedme')
        with open(feedme_out, 'w') as f:
            current_comp = None
            for line in lines:
                if line.startswith("0) sersic"):
                    current_comp = "sersic"
                elif line.startswith("0) sky"):
                    current_comp = "sky"

                if line.startswith("D)"):
                    f.write("D) psf.fits  # Input PSF file\n")
                elif current_comp == "sersic":
                    if line.startswith("1) ") and "position x y" in line:
                        f.write(f"1) {pos_x:.4f} {pos_y:.4f} 1 1  # position x y\n")
                    elif line.startswith("3) "):
                        f.write(f"3) {mag_val:.4f} 1  # Integrated magnitude\n")
                    elif line.startswith("4) "):
                        f.write(f"4) {re_val:.4f} 1  # effective radius (pix)\n")
                    elif line.startswith("5) "):
                        f.write(f"5) {n_val:.4f} {n_fit_flag}  # sersic index\n")
                    elif line.startswith("9) "):
                        f.write(f"9) {b_a_val:.4f} 1  # Axis ratio (b/a)\n")
                    elif line.startswith("10) "):
                        f.write(f"10) {pa_val:.4f} 1  # Position angle\n")
                    else:
                        f.write(line)
                else:
                    f.write(line)

        setup_count += 1
        print(f"  [OK] {frb}: Stage 2 feedme written")

    print(f"\nSet up {setup_count} Stage 2 (with-PSF + sigma) runs")
    return setup_count


# ================================================================
# STEP 6: Run GALFIT Stage 2
# ================================================================
def step6_run_stage2(root):
    print("\n" + "="*70)
    print("STEP 6: Running GALFIT Stage 2 (with-PSF + sigma)")
    print("="*70)

    runs_dir = os.path.join(root, "Galfit", "runs")
    frbs = sorted([d for d in os.listdir(runs_dir) if os.path.isdir(os.path.join(runs_dir, d))])

    success = 0
    failed = []
    for frb in frbs:
        with_psf_dir = os.path.join(runs_dir, frb, "with_psf_sigma")
        feedme = os.path.join(with_psf_dir, "galfit.feedme")

        if not os.path.exists(feedme):
            continue

        print(f"  Running GALFIT for {frb} (with_psf_sigma)...", end=" ", flush=True)
        saved_dir = os.getcwd()
        os.chdir(with_psf_dir)

        result = subprocess.run(["wsl", "galfit", "galfit.feedme"],
                                capture_output=True, text=True, timeout=120)
        os.chdir(saved_dir)

        fitlog = os.path.join(with_psf_dir, "fit.log")
        if os.path.exists(fitlog):
            print("OK")
            success += 1
        else:
            print("FAILED")
            failed.append(frb)
            if result.stderr:
                print(f"    stderr: {result.stderr[:200]}")

    print(f"\nStage 2 complete: {success} succeeded, {len(failed)} failed")
    if failed:
        print(f"  Failed: {', '.join(failed)}")
    return success


# ================================================================
# STEP 7: Compile all results to CSV
# ================================================================
def step7_compile_results(root):
    print("\n" + "="*70)
    print("STEP 7: Compiling results to galfit_sigma_metrics_summary.csv")
    print("="*70)

    runs_dir = os.path.join(root, "Galfit", "runs")
    output_csv = os.path.join(root, "galfit_sigma_metrics_summary.csv")

    frbs = sorted([d for d in os.listdir(runs_dir) if os.path.isdir(os.path.join(runs_dir, d))])

    results = []
    for frb in frbs:
        no_psf_log = os.path.join(runs_dir, frb, "no_psf_sigma", "fit.log")
        with_psf_log = os.path.join(runs_dir, frb, "with_psf_sigma", "fit.log")

        # Only include FRBs that have at least one sigma run
        if not os.path.exists(no_psf_log) and not os.path.exists(with_psf_log):
            continue

        row = {'FRB': frb}

        if os.path.exists(no_psf_log):
            data = parse_fitlog(no_psf_log)
            for k, v in data.items():
                row[f"{k}_nopsf"] = v
        else:
            for k in ['x', 'y', 'mag', 're', 'n', 'b_a', 'pa', 'chi2nu',
                       'x_err', 'y_err', 'mag_err', 're_err', 'n_err', 'b_a_err', 'pa_err']:
                row[f"{k}_nopsf"] = ''

        if os.path.exists(with_psf_log):
            data = parse_fitlog(with_psf_log)
            for k, v in data.items():
                row[f"{k}_psf"] = v
        else:
            for k in ['x', 'y', 'mag', 're', 'n', 'b_a', 'pa', 'chi2nu',
                       'x_err', 'y_err', 'mag_err', 're_err', 'n_err', 'b_a_err', 'pa_err']:
                row[f"{k}_psf"] = ''

        results.append(row)

    # Write CSV in same format as galfit_metrics_summary.csv
    fieldnames = ['FRB']
    keys = ['chi2nu', 'mag', 'mag_err', 're', 're_err', 'n', 'n_err',
            'b_a', 'b_a_err', 'pa', 'pa_err', 'x', 'x_err', 'y', 'y_err']
    for k in keys:
        fieldnames.append(f"{k}_nopsf")
    for k in keys:
        fieldnames.append(f"{k}_psf")

    with open(output_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"\nSaved {len(results)} rows to {output_csv}")
    return len(results)


# ================================================================
# MAIN
# ================================================================
def main():
    root = project_root()
    print(f"Project root: {root}")

    # Step 1: Create cropped sigma maps
    step1_create_cropped_sigma(root)

    # Step 2: Setup Stage 1 feedme files
    step2_setup_stage1(root)

    # Step 3: Run Stage 1
    step3_run_stage1(root)

    # Step 4 is parse_fitlog (used internally)

    # Step 5: Setup Stage 2 from Stage 1 results
    step5_setup_stage2(root)

    # Step 6: Run Stage 2
    step6_run_stage2(root)

    # Step 7: Compile results
    step7_compile_results(root)

    print("\n" + "="*70)
    print("PIPELINE COMPLETE")
    print("="*70)


if __name__ == "__main__":
    main()
