
import os
import glob
import subprocess
import re
import csv

# Constants
INPUT_DIR = "cropped_host_galaxies"
OUTPUT_DIR = "sersic_fits_output"
SCRIPT = "scripts/calc_inclination_sersic.py"
OVERSAMPLING = 4

def main():
    # Find all PSF files in the psfs directory
    psf_files = sorted(glob.glob(os.path.join("psfs", "*_flux_psf.fits")))
    print(f"Found {len(psf_files)} PSF files to process.")
    
    results = []
    
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    for psf_file in psf_files:
        # Deduce galaxy file: psfs/20220207C_flux_psf.fits -> cropped_host_galaxies/20220207C_flux.fits
        galaxy_filename = os.path.basename(psf_file).replace("_psf.fits", ".fits")
        galaxy_file = os.path.join(INPUT_DIR, galaxy_filename)
        
        if not os.path.exists(galaxy_file):
            print(f"  [SKIP] Galaxy file not found for {psf_file}")
            continue
            
        frb_name = os.path.basename(galaxy_file).replace("_flux.fits", "")
        print(f"Processing {frb_name}...")
        
        output_plot = os.path.join(OUTPUT_DIR, f"{frb_name}_fit.png")
        output_log = os.path.join(OUTPUT_DIR, f"{frb_name}_fit.log")
        
        # Build command - NO CROP, NO WCS inputs (Image is already manually cropped)
        cmd = [
            "python", SCRIPT,
            galaxy_file,
            psf_file,
            "--oversampling", str(OVERSAMPLING),
            "--save_plot", output_plot
        ]
        
        try:
            # Run script
            # Use 'utf-8' for encoding to avoid previous issues
            result = subprocess.run(cmd, capture_output=True, text=True, encoding='utf-8', errors='replace')
            
            # Save log
            with open(output_log, "w", encoding='utf-8') as f:
                f.write(result.stdout)
                if result.stderr:
                    f.write("\n--- STDERR ---\n")
                    f.write(result.stderr)
            
            # Parse output for inclination
            inclination = "N/A"
            axis_ratio = "N/A"
            n_index = "N/A"
            r_eff = "N/A"
            
            stdout = result.stdout
            
            # Regex extraction
            inc_match = re.search(r"Inclination:\s+([\d\.]+)", stdout)
            if inc_match: inclination = inc_match.group(1)
            
            q_match = re.search(r"Axis Ratio \(q\):\s+([\d\.]+)", stdout)
            if q_match: axis_ratio = q_match.group(1)
            
            n_match = re.search(r"Sersic Index \(n\):\s+([\d\.]+)", stdout)
            if n_match: n_index = n_match.group(1)
            
            r_match = re.search(r"Effective Radius:\s+([\d\.]+)", stdout)
            if r_match: r_eff = r_match.group(1)
            
            print(f"  Result: i={inclination} deg, q={axis_ratio}, n={n_index}")
            
            results.append({
                "FRB": frb_name,
                "Inclination": inclination,
                "AxisRatio": axis_ratio,
                "SersicIndex": n_index,
                "Reff": r_eff,
                "Status": "Success" if result.returncode == 0 else "Failed"
            })
                
        except Exception as e:
            print(f"  [ERROR] {e}")
            results.append({
                "FRB": frb_name,
                "Inclination": "Error",
                "Status": "Error"
            })

    # Save summary CSV
    csv_file = os.path.join(OUTPUT_DIR, "summary.csv")
    keys = ["FRB", "Inclination", "AxisRatio", "SersicIndex", "Reff", "Status"]
    
    with open(csv_file, "w", newline="", encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(results)
        
    print(f"Processing complete. Summary saved to {csv_file}")

if __name__ == "__main__":
    main()
