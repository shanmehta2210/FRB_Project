import os
import sys
import numpy as np
import pandas as pd
import time
import pyvo
from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
from astropy import units as u

def query_ls_tractor(ra, dec, radius_arcsec=60.0):
    """Query LS DR10 Tractor for all sources in radius."""
    url = 'https://datalab.noirlab.edu/tap'
    service = pyvo.dal.TAPService(url)
    d = radius_arcsec / 3600.0
    # Refined query to remove band-specific shape columns that may not exist in all bricks.
    # shape_r is typically available as the reference band.
    query = f"""
    SELECT objid, ra, dec, type, flux_g, flux_r, flux_z, shape_r
    FROM ls_dr10.tractor 
    WHERE ra > {ra-d} AND ra < {ra+d} 
      AND dec > {dec-d} AND dec < {dec+d}
    """
    try:
        result = service.search(query).to_table().to_pandas()
        return result
    except Exception as e:
        print(f"  LS Query Error: {e}")
        return pd.DataFrame()

def query_ps1_vizier(ra, dec, radius_arcsec=60.0):
    """Query PS1 via Vizier for all sources in radius."""
    v = Vizier(columns=['RAJ2000', 'DEJ2000', 'objID', 'gmag', 'rmag', 'imag', 'zmag', 'gKmag', 'rKmag'], 
               catalog='II/349/ps1')
    coord = SkyCoord(ra=ra, dec=dec, unit='deg', frame='icrs')
    try:
        results = v.query_region(coord, radius=radius_arcsec*u.arcsec)
        if results and len(results) > 0:
            df = results[0].to_pandas()
            df = df.rename(columns={'RAJ2000': 'ra', 'DEJ2000': 'dec', 'objID': 'objid'})
            return df
        return pd.DataFrame()
    except Exception as e:
        print(f"  PS1 Query Error: {e}")
        return pd.DataFrame()

def main():
    summary_path = 'master_frb_summary.csv'
    if not os.path.exists(summary_path):
        print(f"Error: {summary_path} not found.")
        return
        
    df_full = pd.read_csv(summary_path)
    # Start from index 24 (Row 25), the first expansion FRB 20181112A
    expansion_df = df_full.iloc[24:].copy()
    
    data_dir = 'tools/AstroPath/data'
    os.makedirs(data_dir, exist_ok=True)
    
    extraction_log = []
    
    for idx, row in expansion_df.iterrows():
        frb_id = row['FRB']
        ra = row['RA_deg']
        dec = row['DEC_deg']
        
        print(f"Processing {frb_id} at {ra:.4f}, {dec:.4f}...")
        
        # 1. Try LS DR10
        cand_df = query_ls_tractor(ra, dec)
        source = "LS_DR10"
        
        # 2. If empty, try PS1
        if cand_df.empty:
            print(f"  LS empty/failed. Trying Pan-STARRS1 fallback...")
            cand_df = query_ps1_vizier(ra, dec)
            source = "PS1"
            
        if not cand_df.empty:
            out_path = os.path.join(data_dir, f"{frb_id}_candidates.csv")
            cand_df['source_catalog'] = source
            cand_df.to_csv(out_path, index=False)
            print(f"  Success: Found {len(cand_df)} sources in {source}. Saved to {out_path}")
            extraction_log.append({'FRB': frb_id, 'Status': 'Success', 'Count': len(cand_df), 'Catalog': source})
        else:
            print(f"  Warning: No sources found in either LS or PS1 for {frb_id}.")
            extraction_log.append({'FRB': frb_id, 'Status': 'Not Found', 'Count': 0, 'Catalog': 'None'})
            
    log_df = pd.DataFrame(extraction_log)
    log_df.to_csv('tools/AstroPath/data/extraction_summary.csv', index=False)
    print("\nExtraction complete. Summary saved to tools/AstroPath/data/extraction_summary.csv")

if __name__ == '__main__':
    main()

if __name__ == '__main__':
    main()
