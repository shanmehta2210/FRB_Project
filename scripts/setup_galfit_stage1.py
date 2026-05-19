import os
import csv
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

def to_wsl(path):
    if path is None or path == 'none': return 'none'
    if ':' in path:
        drive, rest = os.path.splitdrive(path)
        drive_letter = drive[0].lower()
        wsl_path = f"/mnt/{drive_letter}" + rest.replace('\\', '/')
        return wsl_path
    return path.replace('\\', '/')

def get_platescale(wcs):
    if not wcs.is_celestial:
        wcs = wcs.celestial
    platescale = np.mean(np.sum(wcs.pixel_scale_matrix**2, axis=0)**0.5)*3600 #arcsec
    return platescale

def generate_feedme(outdir, imgfile, psffile, fwhm, position, mag, r_e, n, axis_ratio, pa):
    # Ensure paths are WSL compatible
    img_wsl = to_wsl(imgfile)
    psf_wsl = to_wsl(psffile) if psffile else 'none'
    
    # Write constraints.txt
    constraints_path = os.path.join(outdir, 'constraints.txt')
    with open(constraints_path, 'w') as f:
        f.write("1 n 0.5 to 6.0\n")
        f.write(f"1 re {0.5 * fwhm:.2f} to 500.0\n")

    # Read image to get shape and wcs
    img, hdr = fits.getdata(imgfile, header=True)
    region = (0, img.shape[1]-1, 0, img.shape[0]-1)
    conv_x = int(region[1] - region[0] + 1) + 15
    conv_y = int(region[3] - region[2] + 1) + 15
    platescale = get_platescale(WCS(hdr))
    zeropoint = 25.0 # default from previous wrappers

    feedme_path = os.path.join(outdir, 'galfit.feedme')
    with open(feedme_path, 'w') as f:
        f.write("===============================================================================\n")
        f.write("# IMAGE and GALFIT CONTROL PARAMETERS\n")
        f.write(f"A) {img_wsl}  # Input data image (FITS file)\n")
        f.write(f"B) out.fits  # Output data image block\n")
        f.write(f"C) none  # Sigma image name (made from data if blank or 'none')\n")
        f.write(f"D) {psf_wsl}  # Input PSF file\n")
        f.write(f"E) 1  # PSF fine sampling factor\n")
        f.write(f"F) none  #Bad pixel mask\n")
        f.write(f"G) constraints.txt  # File with parameter constraints (ASCII file)\n")
        f.write(f"H) {region[0]} {region[1]} {region[2]} {region[3]}  # Image region to fit (xmin xmax ymin ymax)\n")
        f.write(f"I) {conv_x} {conv_y} # Size of convolution box (x y)\n")
        f.write(f"J) {zeropoint:.4f}  # Photometric zeropoint (mag)\n")
        f.write(f"K) {platescale:.4f} {platescale:.4f} # Plate scale (dx dy) [arcsec/pixel]\n")
        f.write(f"O) regular   # Display type (regular, curses, both)\n")
        f.write(f"P) 0 # Choose: 0=optimize, 1=model, 2=imgblock, 3=subcomps\n")
        f.write("\n# INITIAL FITTING PARAMETERS\n")
        f.write("# Component number: 1\n")
        f.write("0) sersic # Component type\n")
        f.write(f"1) {position[0]:.4f} {position[1]:.4f} 1 1 # position x y\n")
        f.write(f"3) {mag:.4f} 1 # Integrated magnitude\n")
        f.write(f"4) {r_e:.4f} 1 # effective radius (pix)\n")
        f.write(f"5) {n:.4f} 1 # sersic index\n")
        f.write("6) 0.0000 0 # ----\n")
        f.write("7) 0.0000 0 # ----\n")
        f.write("8) 0.0000 0 # ----\n")
        f.write(f"9) {axis_ratio:.4f} 1 # Axis ratio (b/a)\n")
        f.write(f"10) {pa:.4f} 1 # Position angle (PA) [deg: Up=0, left =90]\n")
        f.write("Z) 0 # Skip this model? (yes=1,no=0)\n\n")
        f.write("# Component number: 2\n")
        f.write("0) sky # component type\n")
        f.write("1) 0.0000 1 # Sky background\n")
        f.write("2) 0 0 # dsky/dx\n")
        f.write("3) 0 0 # dsky/dy\n")
        f.write("Z) 0 # Skip this model\n")
        f.write("================================================================================\n")


def main():
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    
    img_dir_path = os.path.join(project_root, "cropped_host_galaxies")
    old_output_root = os.path.join(project_root, "Galfit")
    guesses_csv = os.path.join(project_root, "csv_archive/initial_guesses.csv")
    fwhm_csv = os.path.join(project_root, "psf_fwhm_summary.csv")

    # Read FWHMs
    fwhm_dict = {}
    if os.path.exists(fwhm_csv):
        with open(fwhm_csv, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                try:
                    fwhm_dict[row['FRB']] = float(row['Avg_FWHM'])
                except:
                    pass

    # Read Guesses
    if not os.path.exists(guesses_csv):
        print(f"Guesses CSV not found: {guesses_csv}")
        return

    with open(guesses_csv, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    for row in rows:
        filename = row['filename']
        frb_name = filename.replace("_flux.fits", "")
        imgfile = os.path.join(img_dir_path, filename)

        if not os.path.exists(imgfile):
            print(f"Skipping {frb_name}: Image not found")
            continue

        mag = float(row['mag'])
        r_e = float(row['r_e'])
        n_idx = float(row['n'])
        axis_ratio = float(row['axis_ratio'])
        pa = float(row['pa'])
        cx = float(row['x'])
        cy = float(row['y'])
        
        fwhm_val = fwhm_dict.get(frb_name, 3.0)

        # Create localized architecture
        base_dir = os.path.join(old_output_root, "runs", frb_name)
        no_psf_dir = os.path.join(base_dir, "no_psf")
        with_psf_dir = os.path.join(base_dir, "with_psf")
        
        os.makedirs(no_psf_dir, exist_ok=True)
        os.makedirs(with_psf_dir, exist_ok=True)

        generate_feedme(
            outdir=no_psf_dir,
            imgfile=imgfile,
            psffile=None,
            fwhm=fwhm_val,
            position=(cx, cy),
            mag=mag,
            r_e=r_e,
            n=n_idx,
            axis_ratio=axis_ratio,
            pa=pa
        )
        print(f"Set up Stage 1 for {frb_name} in {no_psf_dir}")

if __name__ == "__main__":
    main()
