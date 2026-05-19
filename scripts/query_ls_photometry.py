import pandas as pd
import numpy as np
import pyvo
import math
import os
from typing import Optional, Dict

def _ra_clause(ra: float, dra: float) -> str:
    ra_min = ra - dra
    ra_max = ra + dra
    if ra_min < 0:
        return f"(ra > {ra_min + 360:.8f} OR ra < {ra_max:.8f})"
    if ra_max > 360:
        return f"(ra > {ra_min:.8f} OR ra < {ra_max - 360:.8f})"
    return f"ra > {ra_min:.8f} AND ra < {ra_max:.8f}"

def query_ls_source(
    svc: pyvo.dal.TAPService,
    ra: float,
    dec: float,
    radius_arcsec: float = 2.0,
) -> Optional[Dict]:
    """Queries for the nearest LS DR10 Tractor source."""
    dec_clip = max(-85.0, min(85.0, dec))
    dra = (radius_arcsec / 3600.0) / math.cos(math.radians(dec_clip))
    ddec = radius_arcsec / 3600.0

    query = f"""
    SELECT TOP 10 objid, ra, dec, type, flux_r, flux_ivar_r
    FROM ls_dr10.tractor
    WHERE {_ra_clause(ra, dra)}
      AND dec > {dec - ddec:.8f} AND dec < {dec + ddec:.8f}
      AND flux_r > 0
    """
    
    try:
        tab = svc.search(query).to_table()
    except Exception as e:
        print(f"Error querying {ra, dec}: {e}")
        return None

    if len(tab) == 0:
        return None

    ra_arr = np.array(tab["ra"], dtype=float)
    dec_arr = np.array(tab["dec"], dtype=float)
    dra_as = (ra_arr - ra) * np.cos(np.radians(dec)) * 3600.0
    ddec_as = (dec_arr - dec) * 3600.0
    sep = np.hypot(dra_as, ddec_as)

    idx = int(np.argmin(sep))
    r = tab[idx]
    
    flux_r = float(r["flux_r"])
    flux_ivar_r = float(r["flux_ivar_r"])
    
    # Calculate magnitude and error
    # flux 1 nanomaggy = mag 22.5
    mag_ls = 22.5 - 2.5 * np.log10(flux_r)
    # sigma_flux = 1 / sqrt(ivar)
    # sigma_mag = (2.5 / ln(10)) * (sigma_flux / flux) = 1.0857 / (flux * sqrt(ivar))
    magerr_ls = 1.0857 / (flux_r * np.sqrt(flux_ivar_r)) if flux_ivar_r > 0 else np.nan

    return {
        "objid_ls": int(r["objid"]),
        "ra_ls": float(r["ra"]),
        "dec_ls": float(r["dec"]),
        "type_ls": str(r["type"]),
        "flux_r_ls": flux_r,
        "mag_ls": mag_ls,
        "magerr_ls": magerr_ls,
        "sep_arcsec": float(sep[idx]),
    }

def main():
    input_file = "Photometry/target_region_photometry.csv"
    output_file = "Photometry/photometry_ls_comparison.csv"
    
    if not os.path.exists(input_file):
        print(f"Input file {input_file} not found.")
        return

    df = pd.read_csv(input_file)
    print(f"Found {len(df)} sources in {input_file}")

    svc = pyvo.dal.TAPService("https://datalab.noirlab.edu/tap")
    
    comparison_data = []
    
    for i, row in df.iterrows():
        print(f"Processing source {i+1}/{len(df)} at RA={row.RA:.5f}, Dec={row.Dec:.5f}...")
        ls_match = query_ls_source(svc, row.RA, row.Dec)
        
        entry = row.to_dict()
        if ls_match:
            entry.update(ls_match)
            # Calculate sigma distance: (Our_Mag - LS_Mag) / LS_MagErr
            # Using LS sigmas as requested
            if not np.isnan(ls_match["magerr_ls"]) and ls_match["magerr_ls"] > 0:
                entry["ls_sigma_distance"] = (row.MAG_CALIB - ls_match["mag_ls"]) / ls_match["magerr_ls"]
            else:
                entry["ls_sigma_distance"] = np.nan
        else:
            entry["objid_ls"] = -1
            entry["ls_sigma_distance"] = np.nan
            
        comparison_data.append(entry)

    comp_df = pd.DataFrame(comparison_data)
    comp_df.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")
    
    # Statistical Summary
    valid_matches = comp_df[comp_df["objid_ls"] != -1]
    print(f"\nMatched {len(valid_matches)} out of {len(df)} sources.")
    if len(valid_matches) > 0:
        mean_sig = valid_matches["ls_sigma_distance"].mean()
        std_sig = valid_matches["ls_sigma_distance"].std()
        print(f"Mean LS Sigma Distance: {mean_sig:.3f}")
        print(f"Std LS Sigma Distance: {std_sig:.3f}")

if __name__ == "__main__":
    main()
