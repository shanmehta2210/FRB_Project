import os
import glob
import re
import csv
import numpy as np

def calculate_inclination(b_a, q0=0.2):
    """
    Calculates inclination angle from axis ratio b/a.
    cos^2(i) = (q^2 - q0^2) / (1 - q0^2)
    """
    q = b_a
    if q < q0:
        return 90.0
    
    # cos^2(i)
    cos2_i = (q**2 - q0**2) / (1 - q0**2)
    
    # Avoid numerical errors slightly creating cos2_i > 1 or < 0
    if cos2_i > 1.0: cos2_i = 1.0
    if cos2_i < 0.0: cos2_i = 0.0
    
    i_rad = np.arccos(np.sqrt(cos2_i))
    return np.degrees(i_rad)

def parse_galfit_log(log_path):
    """
    Parses the LAST iteration block from fit.log, and grabs the 
    FIRST sersic component within that block (Component 1 = galaxy).
    Returns n, b_a, re, chi2nu.
    """
    with open(log_path, 'r') as f:
        content = f.read()
        
    # Split by horizontal rules to get blocks
    blocks = content.split("-----------------------------------------------------------------------------")
    
    # Filter out empty blocks
    blocks = [b for b in blocks if "Chi^2/nu" in b]
    
    if not blocks:
        return None, None, None, None
        
    last_block = blocks[-1]
    
    n_val = None
    b_a_val = None
    re_val = None
    chi2nu_val = None
    
    lines = last_block.split('\n')
    for line in lines:
        if "Chi^2/nu =" in line:
            try:
                parts = line.split('=')[1].split(',')
                chi2nu_val = float(parts[0].strip())
            except ValueError:
                pass
                
        if "sersic" in line and ":" in line and n_val is None:
            clean_line = line.replace('(', ' ').replace(')', ' ').replace(',', ' ')
            parts = clean_line.split()
            
            if len(parts) >= 8 and parts[0].strip() == 'sersic':
                try:
                    re_str = parts[5].replace('*', '').replace('[', '').replace(']', '')
                    n_str = parts[6].replace('*', '').replace('[', '').replace(']', '')
                    ba_str = parts[7].replace('*', '').replace('[', '').replace(']', '')
                    
                    re_val = float(re_str)
                    n_val = float(n_str)
                    b_a_val = float(ba_str)
                except ValueError:
                    continue
            
    return n_val, b_a_val, re_val, chi2nu_val

def main():
    root_dir = os.getcwd()
    runs_dir = os.path.join(root_dir, "Galfit", "runs")
    output_csv = os.path.join(root_dir, "galfit_inc_angle.csv")
    
    print(f"Scanning GALFIT output directories in {runs_dir}...")
    results = []
    
    if not os.path.exists(runs_dir):
        print(f"Error: Directory {runs_dir} does not exist.")
        return

    # Use no_psf as the basis for the subdirectories
    subdirs = [d for d in os.listdir(runs_dir) if os.path.isdir(os.path.join(runs_dir, d))]
    print(f"Found {len(subdirs)} subdirectories: {subdirs}")
    
    q0_val = 0.2
    
    for frb_name in subdirs:
        dir_no_psf = os.path.join(runs_dir, frb_name, "no_psf")
        dir_with_psf = os.path.join(runs_dir, frb_name, "with_psf")
        
        log_no_psf = os.path.join(dir_no_psf, "fit.log")
        log_with_psf = os.path.join(dir_with_psf, "fit.log")
        
        row = {
            'frb_name': frb_name,
            'q0': q0_val,
            'chi2nu_nopsf': '',
            'n_nopsf': '',
            're_nopsf': '',
            'b_a_nopsf': '',
            'inc_nopsf': '',
            'chi2nu_psf': '',
            'n_psf': '',
            're_psf': '',
            'b_a_psf': '',
            'inc_psf': '',
            'notes': ''
        }
        
        notes = []
        
        # Parse No PSF
        if os.path.exists(log_no_psf):
            n_1, ba_1, re_1, chi2nu_1 = parse_galfit_log(log_no_psf)
            if n_1 is not None and ba_1 is not None:
                row['n_nopsf'] = n_1
                row['re_nopsf'] = re_1
                row['b_a_nopsf'] = ba_1
                row['chi2nu_nopsf'] = chi2nu_1 if chi2nu_1 is not None else ''
                inc1 = calculate_inclination(ba_1, q0=q0_val)
                row['inc_nopsf'] = f"{inc1:.2f}"
            else:
                notes.append("NoPSF: Parse Error")
        else:
            notes.append("NoPSF: Missing log")
            
        # Parse With PSF
        if os.path.exists(log_with_psf):
            n_2, ba_2, re_2, chi2nu_2 = parse_galfit_log(log_with_psf)
            if n_2 is not None and ba_2 is not None:
                row['n_psf'] = n_2
                row['re_psf'] = re_2
                row['b_a_psf'] = ba_2
                row['chi2nu_psf'] = chi2nu_2 if chi2nu_2 is not None else ''
                inc2 = calculate_inclination(ba_2, q0=q0_val)
                row['inc_psf'] = f"{inc2:.2f}"
            else:
                notes.append("WithPSF: Parse Error")
        else:
            notes.append("WithPSF: Missing log")
            
        row['notes'] = " | ".join(notes)
        results.append(row)
        print(f"Processed {frb_name}")
            
    # Write CSV
    with open(output_csv, 'w', newline='') as f:
        fieldnames = [
            'frb_name', 'q0', 
            'chi2nu_nopsf', 'n_nopsf', 're_nopsf', 'b_a_nopsf', 'inc_nopsf', 
            'chi2nu_psf', 'n_psf', 're_psf', 'b_a_psf', 'inc_psf', 
            'notes'
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        
        writer.writeheader()
        for r in results:
            writer.writerow(r)
            
    print(f"Done. Results saved to {output_csv}")

if __name__ == "__main__":
    main()
