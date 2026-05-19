import os
import glob
import subprocess

LARGE_CUTOUTS_DIR = "large_cutouts"

def main():
    flux_files = sorted(glob.glob(os.path.join(LARGE_CUTOUTS_DIR, "*_flux.fits")))
    print(f"Found {len(flux_files)} files in {LARGE_CUTOUTS_DIR}")
    
    os.makedirs("psfs", exist_ok=True)
    
    for f in flux_files:
        basename = os.path.basename(f)
        frb_name = basename.replace("_flux.fits", "")
        
        cmd = ["python", os.path.join("scripts", "build_catalog_psf.py"), f]
        try:
            subprocess.run(cmd, check=True)
            print(f"    Finished -> {frb_name}")
        except subprocess.CalledProcessError as e:
            print(f"    ERROR running catalog script on {f}: {e}")

if __name__ == "__main__":
    main()
