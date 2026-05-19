import os
import glob
import pandas as pd
import numpy as np
import math

def incl_from_q(q, q0=0.2):
    if not np.isfinite(q) or q < 0: return np.nan
    if q <= q0: return 90.0
    val = (q**2 - q0**2) / (1.0 - q0**2)
    val = min(1.0, max(0.0, val))
    return math.degrees(math.acos(math.sqrt(val)))

def parse_unrestrained_fitlog(log_path):
    if not os.path.exists(log_path): return None
    with open(log_path, 'r') as f: content = f.read()
    
    blocks = content.split('--------------------------------------------\n--------------------------------------------')
    if len(blocks) <= 1:
        blocks = content.split('--------------------------------------------')
        
    for block in blocks:
        if " sersic " not in block or "Chi^2/nu" not in block:
            continue
            
        # Check component count
        if block.count(" sersic ") > 1:
            continue # Skip multi-component attempts
            
        lines = block.strip().split('\n')
        is_unrestrained = False
        data = None
        
        for i, line in enumerate(lines):
            if "sersic" in line and ":" in line:
                if "[" in line or "]" in line:
                    is_unrestrained = False
                    break # Restrained parameter, skip this block
                
                is_unrestrained = True
                
                parts = line.replace('(', ' ').replace(')', ' ').replace(',', ' ').replace('*', ' ').split()
                if len(parts) >= 9 and parts[0] == 'sersic':
                    try:
                        data = {
                            'n': float(parts[6]),
                            'b_a': float(parts[7]),
                            'b_a_err': np.nan # Default
                        }
                    except:
                        is_unrestrained = False
                        break
                        
                # Look at next line for errors
                if i+1 < len(lines) and "sky" not in lines[i+1]:
                    err_line = lines[i+1]
                    err_parts = err_line.replace('(', ' ').replace(')', ' ').replace(',', ' ').replace('*', ' ').split()
                    if len(err_parts) >= 6:
                        try:
                            data['b_a_err'] = float(err_parts[5])
                        except:
                            pass
        
        if is_unrestrained and data is not None:
             return data
             
    return None

def main():
    root = "."
    runs_dir = os.path.join(root, "tools/galfit/runs")
    master_path = os.path.join(root, "galfit_vs_legacy_master.csv")
    
    if not os.path.exists(master_path):
        print("Error: galfit_vs_legacy_master.csv not found.")
        return

    master = pd.read_csv(master_path)
    
    rows = []
    for frb in master['FRB'].unique():
        fitlog_path = os.path.join(runs_dir, frb, "with_psf_sigma", "fit.log")
        data = parse_unrestrained_fitlog(fitlog_path)
        
        if not data: 
            continue
            
        m_row = master[master['FRB'] == frb].iloc[0]
        inc_freen = incl_from_q(data['b_a'])
        
        # Propagate inclination error using numerical derivative approximation
        inc_err = ""
        if pd.notna(data['b_a_err']) and np.isfinite(inc_freen):
            try:
                # Delta inclination ~ |d(inc)/d(q)| * b_a_err
                q = data['b_a']
                dq = data['b_a_err']
                q0 = 0.2
                
                q_up = min(1.0, q + dq)
                q_down = max(0.0, q - dq)
                
                i_up = incl_from_q(q_up, q0)
                i_down = incl_from_q(q_down, q0)
                
                if np.isfinite(i_up) and np.isfinite(i_down):
                    inc_err = round(abs(i_up - i_down) / 2.0, 2)
            except:
                pass
        
        rows.append({
            'FRB': frb,
            'galfit_n': round(data['n'], 2),
            'galfit_inc': round(inc_freen, 2),
            'galfit_err': inc_err,
            'legacy_type': m_row['type_ls'],
            'legacy_n': round(m_row['sersic_ls_fit'], 2) if pd.notna(m_row['sersic_ls_fit']) else np.nan,
            'legacy_inc': round(m_row['ls_inc_deg'], 2),
            'legacy_err': round(m_row['ls_inc_err_deg'], 2) if pd.notna(m_row['ls_inc_err_deg']) else np.nan,
        })

    out = pd.DataFrame(rows)
    out_path = "galfit_vs_legacy_quick_read.csv"
    out.to_csv(out_path, index=False)
    print(f"Saved {len(rows)} rows to {out_path}")

    # Update master
    for row in rows:
        idx = master[master['FRB'] == row['FRB']].index
        if len(idx) > 0:
            master.loc[idx, 'galfit_inc_psf_deg'] = row['galfit_inc']
            master.loc[idx, 'delta_deg_ls_minus_galfit'] = row['legacy_inc'] - row['galfit_inc']
            master.loc[idx, 'sersic_n_fit'] = row['galfit_n']

    master.to_csv(master_path, index=False)
    print(f"Updated {master_path}")
    
    if len(out) > 0:
        out['delta'] = out['legacy_inc'] - out['galfit_inc']
        valid = out[out['delta'].notna()]
        print(f"Mean Delta vs Legacy: {valid['delta'].mean():.2f}")
        print(f"RMSE vs Legacy: {np.sqrt(np.mean(valid['delta']**2)):.2f}")

if __name__ == "__main__":
    main()
