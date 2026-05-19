import pandas as pd
import os

frb_df = pd.read_csv('frb_coordinates.csv')
crop_df = pd.read_csv('cropping_regions.csv')

updates = {
    '20210807D': {'RA_deg': 299.2214, 'DEC_deg': -0.7624},
    '20220319D': {'RA_deg': 32.17791666666667, 'DEC_deg': 71.03525},
    '20220825A': {'RA_deg': 311.9815, 'DEC_deg': 72.585}
}

def deg_to_hms(ra_deg):
    ra_h = ra_deg / 15.0
    h = int(ra_h)
    m = int((ra_h - h) * 60)
    s = ((ra_h - h) * 60 - m) * 60
    return f"{h}:{m:02d}:{s:05.2f}"

def deg_to_dms(dec_deg):
    sign = '+' if dec_deg >= 0 else '-'
    dec_abs = abs(dec_deg)
    d = int(dec_abs)
    m = int((dec_abs - d) * 60)
    s = ((dec_abs - d) * 60 - m) * 60
    return f"{sign if sign == '-' else ''}{d}:{m:02d}:{s:04.1f}"

for frb, coords in updates.items():
    idx = frb_df[frb_df['FRB'] == frb].index
    if not idx.empty:
        frb_df.loc[idx, 'RA_deg'] = coords['RA_deg']
        frb_df.loc[idx, 'DEC_deg'] = coords['DEC_deg']
        frb_df.loc[idx, 'RA_hms'] = deg_to_hms(coords['RA_deg'])
        frb_df.loc[idx, 'DEC_dms'] = deg_to_dms(coords['DEC_deg'])

merged_df = pd.merge(frb_df, crop_df, on='FRB', how='left')
merged_df.to_csv('frb_coordinates.csv', index=False)
print('Merged frb_coordinates.csv successfully.')

if os.path.exists('cropping_regions.csv'):
    os.remove('cropping_regions.csv')
    print('Removed cropping_regions.csv.')
