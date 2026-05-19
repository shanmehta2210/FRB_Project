import pandas as pd
import os
import subprocess

df = pd.read_csv('tools/AstroPath/results/successful_associations.csv')

# Use absolute Windows path provided by the user
ds9_exe = r'C:\SAOImageDS9\ds9.exe'
ds9_cmd = [ds9_exe]
frame = 1

fits_found = []
missing = []

for idx, row in df.iterrows():
    frb = row['FRB']
    fname = os.path.abspath(f'large_cutouts/{frb}_flux.fits')
    
    if not os.path.exists(fname):
        missing.append(frb)
        continue
        
    ra = row['Host_RA']
    dec = row['Host_Dec']
    fits_found.append(frb)
    
    # Configure DS9 frame - Center on host ra/dec, no markers
    ds9_cmd.extend([
        '-frame', str(frame),
        '-fits', fname,
        '-zscale', 
        '-cmap', 'b',
        '-zoom', '2',
        '-pan', 'to', str(ra), str(dec), 'wcs', 'fk5'
    ])
    frame += 1

if missing:
    print(f"Missing large cutouts for {len(missing)} FRBs: {', '.join(missing)}")

if fits_found:
    print(f"Launching DS9 with {len(fits_found)} frames: {', '.join(fits_found)}")
    # Try setting a common DISPLAY for WSL if not set
    env = os.environ.copy()
    if 'DISPLAY' not in env:
        # Standard localhost or WSLg address
        env['DISPLAY'] = ':0'
    subprocess.Popen(ds9_cmd, env=env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
else:
    print("No large cutouts found.")
