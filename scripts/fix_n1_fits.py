import os
import glob
import subprocess
import re

runs_dir = "Galfit/runs"

def get_sersic_n(fit_log_path):
    if not os.path.exists(fit_log_path):
        return None
    with open(fit_log_path, 'r') as f:
        # Read the last few lines to find the final converged sersic index
        lines = f.readlines()
        
    # parse from the bottom up to find the last iteration
    for line in reversed(lines):
        if "sersic" in line and "Chi" not in line and "Init" not in line and "Restart" not in line and "Output" not in line:
            # typical line:
            # sersic    : (  104.50,   105.20)   18.05     24.29    0.51    0.22    43.81
            parts = line.split(')')
            if len(parts) > 1:
                params = parts[1].split()
                if len(params) >= 4:
                    try:
                        n_val_str = params[3].replace('[', '').replace(']', '').replace('*', '')
                        return float(n_val_str)
                    except:
                        pass
    return None

def fix_n_in_feedme(feedme_path):
    if not os.path.exists(feedme_path):
        return False
    with open(feedme_path, 'r') as f:
        content = f.read()
    
    # replace '5) <value> 1 # sersic index' with '5) 1.0000 0 # sersic index'
    # we only want to do this for the FIRST sersic component (component 1)
    
    lines = content.split('\n')
    in_comp1 = False
    modified = False
    for i, line in enumerate(lines):
        if line.startswith('# Component number: 1'):
            in_comp1 = True
        elif line.startswith('# Component number: 2'):
            in_comp1 = False
            
        if in_comp1 and line.startswith(' 5)'):
            if ' 0 0' not in line and ' 0 #' not in line:
                lines[i] = " 5) 1.0000 0 # sersic index"
                modified = True
            elif ' 1.0000 0' not in line:
                lines[i] = " 5) 1.0000 0 # sersic index"
                modified = True
    
    if modified:
        with open(feedme_path, 'w') as f:
            f.write('\n'.join(lines))
    return modified

for frb in os.listdir(runs_dir):
    frb_dir = os.path.join(runs_dir, frb)
    if not os.path.isdir(frb_dir): continue
    
    for stage in ['no_psf', 'with_psf']:
        stage_dir = os.path.join(frb_dir, stage)
        fit_log = os.path.join(stage_dir, 'fit.log')
        feedme = os.path.join(stage_dir, 'galfit.feedme')
        
        n_val = get_sersic_n(fit_log)
        
        if n_val is not None and n_val < 2.0:
            print(f"FRB {frb} {stage} has n={n_val:.2f}. Fixing n=1 and running GALFIT.")
            modified = fix_n_in_feedme(feedme)
            if modified:
                # Run galfit
                cwd = os.getcwd()
                os.chdir(stage_dir)
                subprocess.run(['wsl', 'galfit', 'galfit.feedme'], capture_output=True)
                os.chdir(cwd)
            else:
                print(f"  -> n=1 already locked in {stage} for {frb}.")

print("Completed checking and running n=1 fits.")
