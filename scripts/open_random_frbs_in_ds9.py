
import pandas as pd
import random

def open_random_frbs_in_ds9():
    """
    Selects two random processed FRBs and prints the ds9 command to open their flux and sigma images.
    """
    try:
        frbs = pd.read_csv("master_frb_summary.csv")
        processed_frbs = frbs[frbs['status'] == 'processed']
        
        if len(processed_frbs) < 2:
            print("Not enough processed FRBs to select two.")
            return

        random_frbs = processed_frbs.sample(n=2)
        
        file_paths = []
        for index, frb in random_frbs.iterrows():
            frb_name = frb['FRB']
            file_paths.append(f"host_galaxies_data/{frb_name}_flux.fits")
            file_paths.append(f"host_galaxies_data/{frb_name}_sigma.fits")
            
        ds9_command = f"ds9 {' '.join(file_paths)}"
        print("Copy and paste the following command into your terminal to open the files in ds9:")
        print(ds9_command)

    except FileNotFoundError:
        print("Error: master_frb_summary.csv not found.")

if __name__ == '__main__':
    open_random_frbs_in_ds9()
