import os
import sys
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u

# Add astropath to sys.path
sys.path.append('tools/AstroPath/astropath_pkg')
from astropath import path

def process_frb(row, data_dir, results_dir):
    frb_id = row['FRB']
    print(f"Analyzing {frb_id}...")
    
    # 1. Load Candidate Data
    cand_path = os.path.join(data_dir, f"{frb_id}_candidates.csv")
    if not os.path.exists(cand_path):
        print(f"  Warning: No candidate CSV found for {frb_id}. Skipping.")
        return None
        
    cand_df = pd.read_csv(cand_path)
    if cand_df.empty:
        print(f"  Warning: {frb_id} candidate list is empty. Skipping.")
        return None

    # 2. Clean Candidate Data
    survey = cand_df['source_catalog'].iloc[0]
    if survey == 'LS_DR10':
        # Filter out stars (PSF type)
        cand_df = cand_df[cand_df['type'] != 'PSF'].copy()
        # Convert flux to mag
        cand_df = cand_df[cand_df['flux_r'] > 0].copy()
        cand_df['mag'] = -2.5 * np.log10(cand_df['flux_r']) + 22.5
        cand_df['ang_size'] = cand_df['shape_r'].fillna(0.5)
        cand_df.loc[cand_df['ang_size'] <= 0, 'ang_size'] = 0.5
    else: # PS1
        # Use rmag - rKmag > 0.05 for extendedness if possible
        if 'rKmag' in cand_df.columns:
            cand_df = cand_df[(cand_df['rmag'] - cand_df['rKmag']) > 0.05].copy()
        # Use rmag directly
        cand_df = cand_df.rename(columns={'rmag': 'mag'})
        cand_df = cand_df.dropna(subset=['mag']).copy()
        # Fallback ang_size for PS1 Vizier
        cand_df['ang_size'] = 0.5
        
    if cand_df.empty:
        print(f"  Warning: {frb_id} has no valid magnitude candidates. Skipping.")
        return None

    # 3. Setup Localization Ellipse with Fallbacks
    ra_frb, dec_frb = row['RA_deg'], row['DEC_deg']
    frb_coord = SkyCoord(ra=ra_frb, dec=dec_frb, unit='deg', frame='icrs')
    
    # Cascade: major/minor -> ra_err/dec_err -> 0.5
    a = row.get('major_sigma_as', np.nan)
    b = row.get('minor_sigma_as', np.nan)
    pa = row.get('pa_deg', 0.0) # Explicit default to avoid NaN rotation matrix
    
    if pd.isna(a) or pd.isna(b):
        a = row.get('ra_err_as', 0.5)
        b = row.get('dec_err_as', 0.5)
    
    # Clamp to reasonable minimums
    a = max(a, 0.1) if not pd.isna(a) else 0.5
    b = max(b, 0.1) if not pd.isna(b) else 0.5
    if pd.isna(pa): pa = 0.0

    print(f"  Localization: a={a:.2f}\", b={b:.2f}\", PA={pa:.1f}deg")

    # 4. PATH Execution
    try:
        mypath = path.PATH()
        mypath.init_localization('eellipse', center_coord=frb_coord, 
                                 eellipse={'a': a, 'b': b, 'theta': pa})
        
        mypath.init_candidates(ra=cand_df['ra'].values, 
                               dec=cand_df['dec'].values, 
                               ang_size=cand_df['ang_size'].values, 
                               mag=cand_df['mag'].values)
                               
        mypath.init_cand_prior(P_O_method='inverse', P_U=0.1)
        mypath.init_theta_prior(PDF='exp', max=60.0, scale=2.0)
        mypath.calc_priors()
        
        P_Oix, P_Ux = mypath.calc_posteriors(method='local', step_size=0.1, max_radius=60.0)
        
        cand_df['prior_O'] = mypath.prior_Oi
        cand_df['posterior_O'] = P_Oix
        
        # Save per-FRB results
        cand_df.to_csv(os.path.join(results_dir, f"{frb_id}_posterior.csv"), index=False)
        
        # Extract best candidate info
        best_idx = cand_df['posterior_O'].idxmax()
        best = cand_df.loc[best_idx]
        
        # Calc distance to best candidate
        best_coord = SkyCoord(ra=best['ra'], dec=best['dec'], unit='deg')
        sep = frb_coord.separation(best_coord).arcsec
        
        return {
            'FRB': frb_id,
            'Source_Catalog': survey,
            'loc_a': a, 'loc_b': b, 'loc_pa': pa,
            'P_U': P_Ux,
            'Best_ObjID': best['objid'],
            'Best_Mag': best['mag'],
            'Separation_as': sep,
            'Best_Posterior': best['posterior_O']
        }
        
    except Exception as e:
        print(f"  Error processing {frb_id}: {e}")
        return None

def main():
    summary_path = 'master_frb_summary.csv'
    data_dir = 'tools/AstroPath/data'
    results_dir = 'tools/AstroPath/results'
    os.makedirs(results_dir, exist_ok=True)
    
    df_full = pd.read_csv(summary_path)
    # Expansion targets starts index 23 (20181112A)
    expansion_df = df_full.iloc[23:].copy()
    
    summary_data = []
    for idx, row in expansion_df.iterrows():
        res = process_frb(row, data_dir, results_dir)
        if res:
            summary_data.append(res)
            
    if summary_data:
        res_df = pd.DataFrame(summary_data)
        res_df.to_csv(os.path.join(results_dir, 'astropath_expansion_summary.csv'), index=False)
        print(f"\nAnalysis complete. Results saved to {results_dir}")

if __name__ == '__main__':
    main()
