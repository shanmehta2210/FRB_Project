import math, os
import numpy as np
import pandas as pd

Q0 = 0.2

def strip_ast(val):
    try:
        if pd.isna(val): return None
        return float(str(val).replace('*', '').strip())
    except:
        return None

def incl_from_q(q, q0=Q0):
    try:
        q = float(q)
        if not np.isfinite(q): return 0.0
        if q <= q0: return 90.0
        val = (q*q - q0*q0) / (1.0 - q0*q0)
        result = math.degrees(math.acos(math.sqrt(min(1.0, max(0.0, val)))))
        return result if np.isfinite(result) else 0.0
    except:
        return 0.0

def ba_err_to_inc_err(ba, ba_err, q0=Q0):
    try:
        ba, ba_err = float(ba), float(ba_err)
        if ba_err <= 0: return 0.0
        q_up   = min(1.0, ba + ba_err)
        q_down = max(1e-9, ba - ba_err)
        i_up   = incl_from_q(q_up)
        i_down = incl_from_q(q_down)
        err = abs(i_up - i_down) / 2.0
        if err == 0.0 and ba <= q0:
            i_boundary = incl_from_q(min(1.0, q0 + ba_err))
            err = abs(90.0 - i_boundary)
        return round(err, 4)
    except:
        return 0.0

def parse_second_to_last_block(log_path):
    if not os.path.exists(log_path): return {}
    with open(log_path, 'r') as f: content = f.read()
    
    blocks = content.split('--------------------------------------------\n--------------------------------------------')
    if len(blocks) <= 1:
        blocks = content.split('--------------------------------------------')
    
    sersic_blocks = [b for b in blocks if " sersic " in b and "Chi^2/nu" in b and b.count(" sersic ") == 1]
    
    if len(sersic_blocks) >= 2:
        block = sersic_blocks[-2]
    elif len(sersic_blocks) == 1:
        block = sersic_blocks[0]
    else:
        return {}
        
    lines = block.strip().split('\n')
    data = {}
    
    for i, line in enumerate(lines):
        if "Chi^2/nu" in line:
            try:
                data['chi2nu'] = float(line.split('=')[1].split(',')[0].strip())
            except: pass
        if "sersic" in line and ":" in line:
            parts = line.replace('(', ' ').replace(')', ' ').replace(',', ' ').replace('*', ' ').split()
            if len(parts) >= 9 and parts[0] == 'sersic':
                try:
                    data['x'] = float(parts[2])
                    data['y'] = float(parts[3])
                    data['mag'] = float(parts[4])
                    data['re'] = float(parts[5])
                    data['n'] = float(parts[6])
                    data['b_a'] = float(parts[7])
                    data['pa'] = float(parts[8])
                except Exception as e: pass
                
            if i+1 < len(lines) and "sky" not in lines[i+1]:
                err_parts = lines[i+1].replace('(', ' ').replace(')', ' ').replace(',', ' ').replace('*', ' ').split()
                if len(err_parts) >= 6:
                    try:
                        data['x_err'] = float(err_parts[0])
                        data['y_err'] = float(err_parts[1])
                        data['mag_err'] = float(err_parts[2])
                        data['re_err'] = float(err_parts[3])
                        data['n_err'] = float(err_parts[4])
                        data['b_a_err'] = float(err_parts[5])
                        data['pa_err'] = float(err_parts[6])
                    except: pass
                    
    return data

def main():
    root = "c:/Users/lenovo/Desktop/Bhardwajetal_2024_nature_inclination_angle-main"
    summary_path = os.path.join(root, "master_frb_summary.csv")
    new16_path = os.path.join(root, "Archive", "csv", "galfit", "new_16_frbs_galfit_results.csv")
    runs_dir = os.path.join(root, "tools/galfit/runs")
    
    m_df = pd.read_csv(summary_path)
    new_df = pd.read_csv(new16_path)
    
    # 1. Fix existing asterisks in new16
    for col in new_df.columns:
        if col != 'FRB':
            new_df[col] = new_df[col].apply(strip_ast)
    
    new_16_list = new_df['FRB'].tolist()
    
    # 2. Get old 23 FRBs (those with inc_psf in master summary but NOT in new16 list)
    old_frbs = m_df[m_df['inc_psf'].notna() & ~m_df['FRB'].isin(new_16_list)]['FRB'].tolist()
    
    # Exclude any FRBs missing from runs_dir
    old_frbs = [f for f in old_frbs if os.path.isdir(os.path.join(runs_dir, f))]
    print(f"Found {len(old_frbs)} old FRBs to process.")
    
    rows = []
    for frb in old_frbs:
        nopsf_log = os.path.join(runs_dir, frb, "no_psf_sigma", "fit.log")
        psf_log = os.path.join(runs_dir, frb, "with_psf_sigma", "fit.log")
        
        nopsf_data = parse_second_to_last_block(nopsf_log)
        psf_data = parse_second_to_last_block(psf_log)
        
        row = {'FRB': frb}
        for k, v in nopsf_data.items():
            row[f'{k}_nopsf'] = v
        for k, v in psf_data.items():
            row[f'{k}_psf'] = v
            
        rows.append(row)
        
    old_df = pd.DataFrame(rows)
    
    # Append
    combined_df = pd.concat([new_df, old_df], ignore_index=True)
    
    # Recompute inclinations for all
    inc_nopsf, inc_err_nopsf, inc_psf, inc_err_psf = [], [], [], []
    for _, row in combined_df.iterrows():
        ba_n = strip_ast(row.get('b_a_nopsf'))
        bae_n = strip_ast(row.get('b_a_err_nopsf')) or 0.0
        ba_p = strip_ast(row.get('b_a_psf'))
        bae_p = strip_ast(row.get('b_a_err_psf')) or 0.0
        
        inc_nopsf.append(round(incl_from_q(ba_n), 4) if ba_n is not None else 0.0)
        inc_err_nopsf.append(round(ba_err_to_inc_err(ba_n, bae_n), 4) if ba_n is not None else 0.0)
        inc_psf.append(round(incl_from_q(ba_p), 4) if ba_p is not None else 0.0)
        inc_err_psf.append(round(ba_err_to_inc_err(ba_p, bae_p), 4) if ba_p is not None else 0.0)
        
    combined_df['inc_nopsf'] = inc_nopsf
    combined_df['inc_err_nopsf'] = inc_err_nopsf
    combined_df['inc_psf'] = inc_psf
    combined_df['inc_err_psf'] = inc_err_psf
    
    # Save back
    combined_df.to_csv(new16_path, index=False)
    print(f"Saved merged results {len(combined_df)} rows directly back to {new16_path}")

if __name__ == "__main__":
    main()
