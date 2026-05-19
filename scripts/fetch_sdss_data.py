
import pandas as pd
from astroquery.sdss import SDSS
from astropy.coordinates import SkyCoord
from astropy.nddata.utils import Cutout2D
from astropy.io import fits
import astropy.units as u
import os
from astropy.wcs import WCS

# Create a directory to store the new cutouts
if not os.path.exists('sdss_cutouts'):
    os.makedirs('sdss_cutouts')

# Read the FRB coordinates CSV
try:
    df = pd.read_csv('master_frb_summary.csv')
except FileNotFoundError:
    print("Error: master_frb_summary.csv not found. Please run the script to create it first.")
    exit()

# Loop through all FRBs
for index, row in df.iterrows():
    frb_name = row['FRB']
    ra = row['RA_deg']
    dec = row['DEC_deg']

    print(f"\nProcessing {frb_name} at (RA, Dec) = ({ra}, {dec})")

    # Create a SkyCoord object
    coords = SkyCoord(ra=ra*u.degree, dec=dec*u.degree, frame='icrs')

    hdu_list = None  # Initialize hdu_list to None
    try:
        # Get the image plates
        images = SDSS.get_images(coordinates=coords, band='r')

        if images:
            # Assuming the first image in the list is the one we want
            hdu_list = images[0]
            
            # Print HDU info to inspect for the error map
            print("\n--- HDU Info ---")
            hdu_list.info()
            print("----------------\n")

            # The primary image is usually in the first HDU
            primary_hdu = hdu_list[0]
            wcs = WCS(primary_hdu.header)

            # Let's assume for now that the error map is in another HDU.
            # A common practice is to have it in an HDU named 'ERR' or similar,
            # or just be the next HDU. We will need to confirm this by inspecting the hdu_list.info() output.
            
            error_hdu = None
            # Look for an HDU named 'ERR' or 'ERROR'
            for hdu in hdu_list:
                if hdu.name in ('ERR', 'ERROR'):
                    error_hdu = hdu
                    print(f"Found error map in HDU named '{hdu.name}'")
                    break
            
            # If no named error HDU, look for one with similar dimensions
            if not error_hdu:
                for hdu in hdu_list:
                    if hdu.data is not None and hdu.data.shape == primary_hdu.data.shape and hdu is not primary_hdu:
                        error_hdu = hdu
                        print(f"Found potential error map in HDU index {hdu_list.index(hdu)}")
                        break

            # Define the cutout size
            cutout_size = u.Quantity(1, u.arcmin)

            # Create the cutout from the primary image
            cutout_image = Cutout2D(primary_hdu.data, coords, cutout_size, wcs=wcs, mode='trim')
            
            # Save the cutout image
            cutout_image_filename = f"sdss_cutouts/{frb_name}_cutout_r.fits"
            # Update the WCS in the header for the cutout
            cutout_header = cutout_image.wcs.to_header()
            fits.writeto(cutout_image_filename, cutout_image.data, header=cutout_header, overwrite=True)
            print(f"Saved image cutout to {cutout_image_filename}")

            # If we found a potential error map, create a cutout from it too
            if error_hdu:
                cutout_error = Cutout2D(error_hdu.data, coords, cutout_size, wcs=wcs, mode='trim')
                cutout_error_filename = f"sdss_cutouts/{frb_name}_cutout_r_sigma.fits"
                # Update the WCS in the header for the cutout
                cutout_error_header = cutout_error.wcs.to_header()
                fits.writeto(cutout_error_filename, cutout_error.data, header=cutout_error_header, overwrite=True)
                print(f"Saved error map cutout to {cutout_error_filename}")
            else:
                print("Could not find a separate HDU for the error map.")

        else:
            print(f"No SDSS images found for {frb_name}")

    except Exception as e:
        print(f"An error occurred while processing {frb_name}: {e}")

    finally:
        if hdu_list:
            hdu_list.close()

print("\nFinished processing all FRBs.")
