import os
import csv

runs_dir = "Galfit/runs"

def extract_fit_blocks(fit_log_path):
    if not os.path.exists(fit_log_path):
        return []
    with open(fit_log_path, 'r') as f:
        content = f.read()
    
    blocks = content.split('-----------------------------------------------------------------------------')
    
    parsed_blocks = []
    for b in blocks:
        if 'Chi^2/nu' not in b: continue
        
        lines = b.strip().split('\n')
        n_val = None
        ba_val = None
        chi2nu = None
        for line in lines:
            if line.startswith(' sersic '):
                parts = line.split(')')
                if len(parts) > 1:
                    params = parts[1].split()
                    if len(params) >= 4:
                        try:
                            # Index 2 is Sersic n, Index 3 is b/a
                            n_val_str = params[2].replace('[', '').replace(']', '').replace('*', '')
                            n_val = float(n_val_str)
                            ba_val_str = params[3].replace('[', '').replace(']', '').replace('*', '')
                            ba_val = float(ba_val_str)
                        except:
                            pass
            elif 'Chi^2/nu' in line:
                try:
                    chi2_str = line.split('=')[1].strip()
                    chi2nu = float(chi2_str)
                except:
                    pass
        
        if chi2nu is not None and n_val is not None and ba_val is not None:
             parsed_blocks.append({'n': n_val, 'ba': ba_val, 'chi2nu': chi2nu})
             
    return parsed_blocks

csv_data = []

for frb in sorted(os.listdir(runs_dir)):
    frb_dir = os.path.join(runs_dir, frb)
    if not os.path.isdir(frb_dir): continue
    
    for stage in ['no_psf', 'with_psf']:
        log_path = os.path.join(frb_dir, stage, 'fit.log')
        blocks = extract_fit_blocks(log_path)
        
        if not blocks: continue
        
        if len(blocks) == 1:
            orig_block = blocks[0]
            final_block = blocks[0]
        else:
            final_block = blocks[-1]
            orig_block = blocks[-2]
            
        orig_n = orig_block['n']
        orig_ba = orig_block['ba']
        orig_chi2 = orig_block['chi2nu']
        
        new_n = final_block['n']
        new_ba = final_block['ba']
        new_chi2 = final_block['chi2nu']
        
        diff_pct = 0.0
        if orig_chi2 > 0:
            diff_pct = 100 * (new_chi2 - orig_chi2) / orig_chi2
            
        csv_data.append({
            'Galaxy': frb,
            'Stage': stage,
            'Original_n': orig_n,
            'New_n': new_n,
            'Original_ba': orig_ba,
            'New_ba': new_ba,
            'Original_Chi2nu': orig_chi2,
            'New_Chi2nu': new_chi2,
            'Chi2_Change_Pct': diff_pct
        })

output_csv = "n1_comparison_results.csv"
with open(output_csv, 'w', newline='') as csvfile:
    fieldnames = ['Galaxy', 'Stage', 'Original_n', 'New_n', 'Original_ba', 'New_ba', 'Original_Chi2nu', 'New_Chi2nu', 'Chi2_Change_Pct']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in csv_data:
        writer.writerow(row)

print(f"Results written to {output_csv}")
