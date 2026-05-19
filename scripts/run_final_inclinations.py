import os
import subprocess
import glob
import csv
import numpy as np

def calculate_inclination(b_a, q0=0.2):
    q = b_a
    if q < q0:
        return 90.0
    cos2_i = (q**2 - q0**2) / (1 - q0**2)
    if cos2_i > 1.0: cos2_i = 1.0
    if cos2_i < 0.0: cos2_i = 0.0
    i_rad = np.arccos(np.sqrt(cos2_i))
    return np.degrees(i_rad)

def extract_ba_from_log(log_path):
    """
    STRICTLY extracts the final b/a from the LAST instance of a Sersic
    component block in the fit.log file. For multi-component fits,
    we extract the b/a of the first Sersic component in that block.
    """
    if not os.path.exists(log_path):
        return None, None, None
        
    with open(log_path, 'r') as f:
        content = f.read()
        
    blocks = content.split("-----------------------------------------------------------------------------")
    blocks = [b for b in blocks if "Chi^2/nu" in b]
    
    if not blocks:
        return None, None, None
        
    last_block = blocks[-1]
    
    n_val = None
    b_a_val = None
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
                    n_str = parts[6].replace('*', '').replace('[', '').replace(']', '')
                    ba_str = parts[7].replace('*', '').replace('[', '').replace(']', '')
                    n_val = float(n_str)
                    b_a_val = float(ba_str)
                except ValueError:
                    continue
                    
    return b_a_val, n_val, chi2nu_val


def force_n1_and_run():
    runs_dir = "Galfit/runs"
    frb_list = sorted(os.listdir(runs_dir))
    
    print("Step 1: Forcing n=1 and running GALFIT for all cases (except 20171020A no_psf)...")
    
    for frb in frb_list:
        frb_dir = os.path.join(runs_dir, frb)
        if not os.path.isdir(frb_dir): continue
        
        for stage in ['no_psf', 'with_psf']:
            # Skip 20171020A no_psf as requested by user
            if frb == "20171020A" and stage == "no_psf":
                print(f"Skipping {frb} {stage} (retaining n=0.55)...")
                continue
                
            stage_dir = os.path.join(frb_dir, stage)
            feedme = os.path.join(stage_dir, 'galfit.feedme')
            
            if not os.path.exists(feedme): continue
            
            with open(feedme, 'r') as f:
                content = f.read()
                
            lines = content.split('\n')
            in_comp1 = False
            modified = False
            for i, line in enumerate(lines):
                if line.startswith('# Component number: 1'):
                    in_comp1 = True
                elif line.startswith('# Component number: 2'):
                    in_comp1 = False
                    
                if in_comp1 and line.strip().startswith('5)'):
                    # Force n=1.0 locked
                    new_line = " 5) 1.0000 0 # sersic index"
                    if line.strip() != new_line.strip():
                        lines[i] = new_line
                        modified = True
                        
            if modified:
                with open(feedme, 'w') as f:
                    f.write('\n'.join(lines))
                    
            print(f"Running GALFIT for {frb} [{stage}]...")
            cwd = os.getcwd()
            try:
                os.chdir(stage_dir)
                subprocess.run(['wsl', 'galfit', 'galfit.feedme'], capture_output=True)
            except Exception as e:
                print(f"Failed to run GALFIT in {stage_dir}: {e}")
            finally:
                os.chdir(cwd)

def calculate_all():
    runs_dir = "Galfit/runs"
    frb_list = sorted(os.listdir(runs_dir))
    output_csv = "final_inclination_angles.csv"
    q0_val = 0.2
    
    print("\nStep 2: Parsing Final b/a and Calculating Inclination Angles...")
    results = []
    
    for frb in frb_list:
        frb_dir = os.path.join(runs_dir, frb)
        if not os.path.isdir(frb_dir): continue
        
        row = {
            'frb_name': frb,
            'q0': q0_val,
            'n_nopsf': '', 'b_a_nopsf': '', 'chi2nu_nopsf': '', 'inc_nopsf': '',
            'n_psf': '', 'b_a_psf': '', 'chi2nu_psf': '', 'inc_psf': ''
        }
        
        # No PSF
        log_nopsf = os.path.join(frb_dir, 'no_psf', 'fit.log')
        ba1, n1, chi1 = extract_ba_from_log(log_nopsf)
        if ba1 is not None:
            row['b_a_nopsf'] = ba1
            row['n_nopsf'] = n1
            row['chi2nu_nopsf'] = chi1
            row['inc_nopsf'] = f"{calculate_inclination(ba1, q0_val):.2f}"
            
        # With PSF
        log_psf = os.path.join(frb_dir, 'with_psf', 'fit.log')
        ba2, n2, chi2 = extract_ba_from_log(log_psf)
        if ba2 is not None:
            row['b_a_psf'] = ba2
            row['n_psf'] = n2
            row['chi2nu_psf'] = chi2
            row['inc_psf'] = f"{calculate_inclination(ba2, q0_val):.2f}"
            
        results.append(row)
        
    with open(output_csv, 'w', newline='') as f:
        fieldnames = ['frb_name', 'q0', 'chi2nu_nopsf', 'n_nopsf', 'b_a_nopsf', 'inc_nopsf', 
                      'chi2nu_psf', 'n_psf', 'b_a_psf', 'inc_psf']
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow(r)
            
    print(f"Done! Results written to {output_csv}")

if __name__ == "__main__":
    force_n1_and_run()
    calculate_all()
