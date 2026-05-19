import pandas as pd
import os

master_path = 'c:/Users/lenovo/Desktop/Bhardwajetal_2024_nature_inclination_angle-main/master_frb_summary.csv'
success_path = 'c:/Users/lenovo/Desktop/Bhardwajetal_2024_nature_inclination_angle-main/tools/AstroPath/results/successful_associations.csv'

if os.path.exists(master_path) and os.path.exists(success_path):
    master = pd.read_csv(master_path)
    success = pd.read_csv(success_path)
    
    # Create notes column if missing
    if 'Notes' not in master.columns:
        master['Notes'] = ''
        
    # Mark exactly the 16 high-confidence FRBs
    high_conf_list = success['FRB'].tolist()
    mask = master['FRB'].isin(high_conf_list)
    
    master.loc[mask, 'Notes'] = 'Host association done (Posterior > 0.8)'
    
    master.to_csv(master_path, index=False)
    print(f"Successfully marked 'Notes' for {len(high_conf_list)} high-confidence FRBs in {master_path}.")
else:
    print("Required CSV files not found.")
