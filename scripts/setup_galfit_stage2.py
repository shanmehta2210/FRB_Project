import os
import shutil

def parse_fitlog(log_path):
    with open(log_path, 'r') as f:
        content = f.read()
    
    # GALFIT separates fits with a long dashed line
    blocks = content.split('--------------------------------------------')
    
    # Work backwards to find the last complete block containing a sersic component
    target_block = None
    for block in reversed(blocks):
        if "sersic" in block and "Chi^2/nu =" in block:
            target_block = block
            break
            
    n_val = None
    b_a_val = None
    re_val = None
    mag_val = None
    pa_val = None
    pos_x = None
    pos_y = None
    
    if not target_block:
        return pos_x, pos_y, mag_val, re_val, n_val, b_a_val, pa_val

    lines = target_block.strip().split('\n')
    for line in lines:
        if "sersic" in line and ":" in line:
            clean_line = line.replace('(', ' ').replace(')', ' ').replace(',', ' ')
            parts = clean_line.split()
            if len(parts) >= 8 and parts[0].strip() == 'sersic':
                try:
                    pos_x = float(parts[2].replace('*', ''))
                    pos_y = float(parts[3].replace('*', ''))
                    mag_val = float(parts[4].replace('*', ''))
                    re_val = float(parts[5].replace('*', ''))
                    n_val = float(parts[6].replace('*', ''))
                    b_a_val = float(parts[7].replace('*', ''))
                    pa_val = float(parts[8].replace('*', ''))
                except ValueError:
                    continue
    return pos_x, pos_y, mag_val, re_val, n_val, b_a_val, pa_val

def generate_feedme_stage2(outdir, no_psf_feedme_path, local_psf_name, pos_x, pos_y, mag, r_e, n, axis_ratio, pa):
    with open(no_psf_feedme_path, 'r') as f:
        lines = f.readlines()

    feedme_path = os.path.join(outdir, 'galfit.feedme')
    with open(feedme_path, 'w') as f:
        current_comp = None
        for line in lines:
            if line.startswith("0) sersic"):
                current_comp = "sersic"
            elif line.startswith("0) sky"):
                current_comp = "sky"
                
            if line.startswith("D)"):
                f.write(f"D) {local_psf_name}  # Input PSF file\n")
            elif current_comp == "sersic":
                if line.startswith("1) ") and "position x y" in line:
                    f.write(f"1) {pos_x:.4f} {pos_y:.4f} 1 1 # position x y\n")
                elif line.startswith("3) "):
                    f.write(f"3) {mag:.4f} 1 # Integrated magnitude\n")
                elif line.startswith("4) "):
                    f.write(f"4) {r_e:.4f} 1 # effective radius (pix)\n")
                elif line.startswith("5) "):
                    f.write(f"5) {n:.4f} 1 # sersic index\n")
                elif line.startswith("9) "):
                    f.write(f"9) {axis_ratio:.4f} 1 # Axis ratio (b/a)\n")
                elif line.startswith("10) "):
                    f.write(f"10) {pa:.4f} 1 # Position angle (PA) [deg: Up=0, left =90]\n")
                else:
                    f.write(line)
            else:
                f.write(line)
                
    # Also copy constraints.txt
    constraints_src = os.path.join(os.path.dirname(no_psf_feedme_path), 'constraints.txt')
    constraints_dst = os.path.join(outdir, 'constraints.txt')
    if os.path.exists(constraints_src):
        shutil.copy(constraints_src, constraints_dst)

def main():
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    runs_dir = os.path.join(project_root, "Galfit", "runs")
    psf_dir = os.path.join(project_root, "psfs", "downsampled_psfs")
    
    if not os.path.exists(runs_dir):
        print("Runs directory not found.")
        return

    frbs = [d for d in os.listdir(runs_dir) if os.path.isdir(os.path.join(runs_dir, d))]
    
    for frb in frbs:
        no_psf_dir = os.path.join(runs_dir, frb, "no_psf")
        with_psf_dir = os.path.join(runs_dir, frb, "with_psf")
        fitlog_path = os.path.join(no_psf_dir, "fit.log")
        no_psf_feedme_path = os.path.join(no_psf_dir, "galfit.feedme")
        
        psf_file = os.path.join(psf_dir, f"{frb}_1x_psf.fits")
        if not os.path.exists(psf_file):
            print(f"PSF not found for {frb}: {psf_file}")
            continue
            
        # VERY IMPORTANT: To avoid GALFIT buffer overflow with WSL long path, copy PSF locally
        local_psf_name = "psf.fits"
        local_psf_path = os.path.join(with_psf_dir, local_psf_name)
        shutil.copy(psf_file, local_psf_path)
            
        if os.path.exists(fitlog_path) and os.path.exists(no_psf_feedme_path):
            vals = parse_fitlog(fitlog_path)
            if None in vals:
                print(f"Failed to parse {fitlog_path}")
                continue
            
            pos_x, pos_y, mag_val, re_val, n_val, b_a_val, pa_val = vals
            generate_feedme_stage2(with_psf_dir, no_psf_feedme_path, local_psf_name, pos_x, pos_y, mag_val, re_val, n_val, b_a_val, pa_val)
            print(f"Set up Stage 2 (with local localized PSF) for {frb} in {with_psf_dir}")
        else:
            print(f"fit.log or galfit.feedme missing for {frb} in {no_psf_dir}")

if __name__ == "__main__":
    main()
