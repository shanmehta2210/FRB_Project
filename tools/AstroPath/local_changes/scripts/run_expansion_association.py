import os
import sys
import numpy as np
import pandas as pd
import math
import pyvo
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table

# Ensure astropath is importable
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from astropath import path

def query_ls_tractor(ra_deg, dec_deg, radius_arcsec=60.0):
    url = 'https://datalab.noirlab.edu/tap'
    service = pyvo.dal.TAPService(url)
    
    dec_clip = max(-85.0, min(85.0, dec_deg))
    dra = (radius_arcsec / 3600.0) / math.cos(math.radians(dec_clip))
    ddec = radius_arcsec / 3600.0
    
    ra_min = ra_deg - dra
    ra_max = ra_deg + dra
    
    if ra_min < 0:
        ra_clause = f"(ra > {ra_min + 360:.8f} OR ra < {ra_max:.8f})"
    elif ra_max > 360:
        ra_clause = f"(ra > {ra_min:.8f} OR ra < {ra_max - 360:.8f})"
    else:
        ra_clause = f"ra > {ra_min:.8f} AND ra < {ra_max:.8f}"
        
    query = f"""
    SELECT objid, ra, dec, type, flux_r, shape_r
    FROM ls_dr10.tractor
    WHERE {ra_clause}
      AND dec > {dec_deg - ddec:.8f} AND dec < {dec_deg + ddec:.8f}
      AND type != 'PSF' 
      AND flux_r > 0
    """
    try:
        result = service.search(query).to_table().to_pandas()
        if not result.empty:
            result['source_catalog'] = 'LS_DR10'
        return result
    except Exception as e:
        print(f"Error querying LS for RA={ra_deg}, DEC={dec_deg}: {e}")
        return pd.DataFrame()

def query_ps1_candidates(ra_deg, dec_deg, radius_arcsec=60.0):
    """Fallback to Pan-STARRS1 if Legacy Surveys is missing."""
    v = Vizier(columns=['RAJ2000', 'DEJ2000', 'objID', 'rmag', 'rKmag'], catalog='II/349/ps1')
    coord = SkyCoord(ra=ra_deg, dec=dec_deg, unit='deg', frame='icrs')
    try:
        res = v.query_region(coord, radius=radius_arcsec*u.arcsec)
        if not res:
            return pd.DataFrame()
        
        df = res[0].to_pandas()
        # Rename for compatibility
        df = df.rename(columns={'RAJ2000': 'ra', 'DEJ2000': 'dec', 'objID': 'objid'})
        
        # Extendedness check: rmag - rKmag > 0.05
        # Also ensure we have a valid kron mag
        df = df.dropna(subset=['rmag', 'rKmag'])
        df = df[(df['rmag'] - df['rKmag']) > 0.05]
        
        # Map shape_r to a constant 0.5 (or use a better proxy if available)
        df['shape_r'] = 0.5
        df['flux_r'] = 10**((22.5 - df['rKmag'])/2.5) # Fake flux for internal PATH math if needed
        df['source_catalog'] = 'PS1'
        return df
    except Exception as e:
        print(f"Error querying PS1 for RA={ra_deg}, DEC={dec_deg}: {e}")
        return pd.DataFrame()

def process_frbs():
    master_df = pd.read_csv('master_frb_summary.csv')
    
    # Exclude confirmed gaps
    exclude_list = ['20230718A', '20230731A']
    
    expansion_df = master_df.iloc[23:].copy() 
    expansion_df = expansion_df.dropna(subset=['RA_deg', 'DEC_deg'])
    expansion_df = expansion_df[~expansion_df['FRB'].isin(exclude_list)]
    
    out_dir = os.path.join('tools', 'astropath', 'results')
    os.makedirs(out_dir, exist_ok=True)
    
    summary_results = []
    
    for idx, row in expansion_df.iterrows():
        frb_id = row['FRB']
        ra = row['RA_deg']
        dec = row['DEC_deg']
        
        print(f"Processing {frb_id} at {ra}, {dec}...")
        
        # 1. Parse Localization
        if not pd.isna(row.get('major_sigma_as')) and not pd.isna(row.get('minor_sigma_as')):
            a_arcsec = float(row['major_sigma_as'])
            b_arcsec = float(row['minor_sigma_as'])
            pa_deg = float(row.get('pa_deg', 0.0))
        else:
            a_arcsec = float(row.get('ra_err_as', 0.5))
            b_arcsec = float(row.get('dec_err_as', 0.5))
            pa_deg = float(row.get('pa_deg', 0.0))
            if pd.isna(a_arcsec) or a_arcsec == 0: a_arcsec = 0.5
            if pd.isna(b_arcsec) or b_arcsec == 0: b_arcsec = 0.5
                
        # 2. Query Candidates (Successive)
        cand_df = query_ls_tractor(ra, dec, radius_arcsec=60.0)
        
        if cand_df.empty:
            print(f"  No extended sources in LS DR10. Retrying with Pan-STARRS1...")
            cand_df = query_ps1_candidates(ra, dec, radius_arcsec=60.0)
            
        if cand_df.empty:
            print(f"  No extended sources found in LS or PS1 for {frb_id}.")
            continue
            
        # 3. Clean Candidates
        source_cat = cand_df['source_catalog'].iloc[0]
        cand_df['shape_r'] = cand_df['shape_r'].fillna(0.5)
        cand_df.loc[cand_df['shape_r'] == 0, 'shape_r'] = 0.5
        
        if source_cat == 'LS_DR10':
            cand_df['mag_r'] = -2.5 * np.log10(cand_df['flux_r']) + 22.5
        else:
            # PS1 already has rKmag
            cand_df['mag_r'] = cand_df['rKmag']
        
        # 4. Run PATH
        mypath = path.PATH()
        coord = SkyCoord(ra=ra, dec=dec, unit='deg')
        mypath.init_localization('eellipse', center_coord=coord, 
                                 eellipse={'a': a_arcsec, 'b': b_arcsec, 'theta': pa_deg})
                                 
        mypath.init_candidates(ra=cand_df['ra'].values, 
                               dec=cand_df['dec'].values, 
                               ang_size=cand_df['shape_r'].values, 
                               mag=cand_df['mag_r'].values)
                               
        mypath.init_cand_prior(P_O_method='inverse', P_U=0.1)
        mypath.init_theta_prior(PDF='exp', max=30.0, scale=2.0)
        
        mypath.calc_priors()
        P_Oix, P_Ux = mypath.calc_posteriors(method='local', step_size=0.1, max_radius=30.0)
        
        cand_df['prior_O'] = mypath.prior_Oi
        cand_df['posterior_O'] = P_Oix
        
        out_csv = os.path.join(out_dir, f'{frb_id}_posterior.csv')
        cand_df.to_csv(out_csv, index=False)
        
        best_idx = np.argmax(P_Oix)
        best_cand = cand_df.iloc[best_idx]
        
        summary_results.append({
            'FRB': frb_id,
            'Source_Catalog': source_cat,
            'Best_Host_objid': best_cand['objid'],
            'Best_Host_RA': best_cand['ra'],
            'Best_Host_DEC': best_cand['dec'],
            'Best_Host_Mag': best_cand['mag_r'],
            'Best_Host_AngSize': best_cand['shape_r'],
            'Posterior_P_Ox': best_cand['posterior_O'],
            'Unseen_P_U': P_Ux
        })
        
        print(f"  Success ({source_cat}). Best Match Posterior: {best_cand['posterior_O']:.4f}")
        
    if summary_results:
        sum_df = pd.DataFrame(summary_results)
        sum_df.to_csv(os.path.join(out_dir, 'astropath_expansion_summary.csv'), index=False)
        print("Completed processing. Global expansion summary updated.")
    else:
        print("No results to save.")

if __name__ == '__main__':
    process_frbs()
