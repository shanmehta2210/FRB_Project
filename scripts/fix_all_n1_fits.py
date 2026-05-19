import os
import subprocess

runs_dir = "Galfit/runs"

def fix_n_in_feedme(feedme_path):
    if not os.path.exists(feedme_path):
        return False
    with open(feedme_path, 'r') as f:
        content = f.read()
    
    # replace '5) <value> 1 # sersic index' with '5) 1.0000 0 # sersic index'
    # we only want to do this for the FIRST sersic component
    
    lines = content.split('\n')
    in_comp1 = False
    modified = False
    for i, line in enumerate(lines):
        if line.startswith('# Component number: 1'):
            in_comp1 = True
        elif line.startswith('# Component number: 2'):
            in_comp1 = False
            
        if in_comp1 and line.strip().startswith('5)'):
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

frb_list = sorted(os.listdir(runs_dir))

for frb in frb_list:
    frb_dir = os.path.join(runs_dir, frb)
    if not os.path.isdir(frb_dir): continue
    
    print(f"Processing FRB {frb}...")
    for stage in ['no_psf', 'with_psf']:
        stage_dir = os.path.join(frb_dir, stage)
        feedme = os.path.join(stage_dir, 'galfit.feedme')
        
        modified = fix_n_in_feedme(feedme)
        if modified:
            print(f"  [{stage}] Locking n=1 and running GALFIT...")
            cwd = os.getcwd()
            try:
                os.chdir(stage_dir)
                subprocess.run(['wsl', 'galfit', 'galfit.feedme'], capture_output=True)
            finally:
                os.chdir(cwd)
        else:
            print(f"  [{stage}] n=1 already locked. Running GALFIT just to be safe...")
            cwd = os.getcwd()
            try:
                os.chdir(stage_dir)
                subprocess.run(['wsl', 'galfit', 'galfit.feedme'], capture_output=True)
            finally:
                os.chdir(cwd)

print("Completed running n=1 fits for ALL galaxies.")
