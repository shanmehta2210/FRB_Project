import os
import sys
import csv

# Add the project root to the python path so we can import Galfit.galfit_wrap
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from Galfit.galfit_wrap import run

def test_1x_psf():
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    imgfile = os.path.join(project_root, "cropped_host_galaxies", "20171020A_flux.fits")
    psffile = os.path.join(project_root, "psfs", "downsampled_psfs", "20171020A_1x_psf.fits")
    outdir = os.path.join(project_root, "Galfit", "test_1x_psf_20171020A")
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        
    print(f"Testing GALFIT on 20171020A using the 1x Downsampled PSF...")
    print(f"Image: {imgfile}")
    print(f"PSF: {psffile}")
    
    # Needs initial guesses to run
    # Read them directly from initial_guesses.csv
    csv_path = os.path.join(project_root, "initial_guesses.csv")
    mag, r_e, n_idx, axis_ratio, pa, cx, cy = None, None, None, None, None, None, None
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['filename'] == "20171020A_flux.fits":
                mag = float(row['mag'])
                r_e = float(row['r_e'])
                n_idx = float(row['n'])
                axis_ratio = float(row['axis_ratio'])
                pa = float(row['pa'])
                cx = float(row['x'])
                cy = float(row['y'])
                break
                
    if mag is None:
        print("Could not find guesses.")
        return
        
    position = (cx, cy)
    
    ret = run(
        imgfile=imgfile,
        psffile=psffile,
        outdir=outdir,
        outfile="out.fits",
        configfile="galfit.feedme",
        finesample=1, 
        noisefile="none", 
        position=position,
        int_mag=mag,
        r_e=r_e,
        n=n_idx,
        axis_ratio=axis_ratio,
        pa=pa,
        skip_sky=False 
    )
    
    print(f"Exit code: {ret}")

if __name__ == "__main__":
    test_1x_psf()
