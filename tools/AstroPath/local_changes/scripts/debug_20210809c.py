import os
import sys
import numpy as np
import pandas as pd
import math
import pyvo
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy import units as u

# Ensure astropath is importable
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from astropath import path

def query_ls_tractor(ra_deg, dec_deg, radius_arcsec=60.0):
    url = 'https://datalab.noirlab.edu/tap'
    service = pyvo.dal.TAPService(url)
    d = radius_arcsec / 3600.0
    query = f"SELECT objid, ra, dec, type, flux_r, shape_r FROM ls_dr10.tractor WHERE ra > {ra_deg-d} AND ra < {ra_deg+d} AND dec > {dec_deg-d} AND dec < {dec_deg+d} AND type != 'PSF' AND flux_r > 0"
    try:
        result = service.search(query).to_table().to_pandas()
        return result
    except Exception as e:
        print(f"LS Query Error: {e}")
        return pd.DataFrame()

def debug_20210809c():
    ra, dec = 271.1570833333333, 1.3288888888888888
    frb_id = "20210809C"
    
    print(f"DEBUGGING {frb_id} at {ra}, {dec}...")
    
    # 1. Query LS
    cand_df = query_ls_tractor(ra, dec)
    if cand_df.empty:
        print("  LS returned no candidates.")
        # Try PS1 fallback here too
    else:
        print(f"  LS found {len(cand_df)} candidates.")
        
    if cand_df.empty:
        print("  Checking PS1 Fallback...")
        v = Vizier(columns=['RAJ2000', 'DEJ2000', 'objID', 'rmag', 'rKmag'], catalog='II/349/ps1')
        coord = SkyCoord(ra=ra, dec=dec, unit='deg', frame='icrs')
        try:
            res = v.query_region(coord, radius=60*u.arcsec)
            if res:
                df = res[0].to_pandas()
                df = df.rename(columns={'RAJ2000': 'ra', 'DEJ2000': 'dec', 'objID': 'objid'})
                df = df.dropna(subset=['rmag', 'rKmag'])
                df = df[(df['rmag'] - df['rKmag']) > 0.05]
                print(f"  PS1 found {len(df)} candidates after filter.")
                cand_df = df
                source_cat = 'PS1'
            else:
                print("  PS1 returned no results in region.")
        except Exception as e:
            print(f"  PS1 Query Error: {e}")

    if cand_df.empty:
        print("  No candidates found in either catalog.")
        return

    # Try running PATH logic
    try:
        cand_df['shape_r'] = cand_df.get('shape_r', 0.5).fillna(0.5)
        if 'source_catalog' not in cand_df:
            if 'rmag' in cand_df: cand_df['mag_r'] = cand_df['rKmag']
            else: cand_df['mag_r'] = -2.5 * np.log10(cand_df['flux_r']) + 22.5
            
        print("  Initializing PATH...")
        mypath = path.PATH()
        coord = SkyCoord(ra=ra, dec=dec, unit='deg')
        mypath.init_localization('eellipse', center_coord=coord, 
                                 eellipse={'a': 0.5, 'b': 0.5, 'theta': 0.0})
                                 
        mypath.init_candidates(ra=cand_df['ra'].values, 
                               dec=cand_df['dec'].values, 
                               ang_size=cand_df['shape_r'].values, 
                               mag=cand_df['mag_r'].values)
                               
        mypath.init_cand_prior(P_O_method='inverse', P_U=0.1)
        mypath.init_theta_prior(PDF='exp', max=30.0, scale=2.0)
        
        print("  Calculating Posteriors...")
        mypath.calc_priors()
        P_Oix, P_Ux = mypath.calc_posteriors(method='local', step_size=0.1, max_radius=30.0)
        print(f"  SUCCESS! Max Posterior: {np.max(P_Oix):.4f}")
        
    except Exception as e:
        print(f"  PATH EXECUTION ERROR: {e}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    debug_20210809c()
