
import os
import requests
import numpy as np
from astropy.io import fits
from io import BytesIO
import argparse

def fetch_and_process_legacy_survey_data(ra, dec, frb_name, output_dir='host_galaxies_data'):
    """
    Fetches, processes, and saves data from the Legacy Survey for a given RA and Dec.

    Args:
        ra (float): Right Ascension in degrees.
        dec (float): Declination in degrees.
        frb_name (str): The name of the FRB.
        output_dir (str): Directory to save the output FITS files.
    """
    print(f"Processing {frb_name} (RA={ra}, Dec={dec})")

    # 1. Construct URL
    url = (
        f"https://www.legacysurvey.org/viewer/cutout.fits?"
        f"ra={ra}&dec={dec}&layer=ls-dr10&pixscale=0.262&invvar=True"
    )
    print(f"Fetching data from: {url}")

    try:
        # 2. Fetch Data
        response = requests.get(url, timeout=60)
        response.raise_for_status()  # Raise an exception for bad status codes
        
        with fits.open(BytesIO(response.content)) as hdul:
            # 3. Process HDUs
            if len(hdul) < 2:
                print("Error: FITS file does not contain the expected number of HDUs (at least 2).")
                return

            flux_image = hdul[0].data
            invvar_map = hdul[1].data
            
            print("Successfully extracted Flux Image (HDU 0) and Inverse Variance Map (HDU 1).")

            # 4. Compute Sigma Map
            sigma_map = np.full_like(invvar_map, np.inf, dtype=np.float32)
            valid_invvar = invvar_map > 0
            sigma_map[valid_invvar] = 1.0 / np.sqrt(invvar_map[valid_invvar])
            
            print("Computed Sigma Map. Pixels with invvar==0 are set to infinity.")

            # 5. Output
            os.makedirs(output_dir, exist_ok=True)
            
            base_filename = frb_name
            
            # Save Flux Image
            flux_hdu = fits.PrimaryHDU(flux_image, header=hdul[0].header)
            flux_hdu.header['HISTORY'] = 'Flux Image (nanomaggies)'
            flux_filename = os.path.join(output_dir, f"{base_filename}_flux.fits")
            flux_hdu.writeto(flux_filename, overwrite=True)
            print(f"Saved Flux Image to: {flux_filename}")

            # Save Inverse Variance Map
            invvar_hdu = fits.PrimaryHDU(invvar_map, header=hdul[1].header)
            invvar_hdu.header['HISTORY'] = 'Inverse Variance Map'
            invvar_filename = os.path.join(output_dir, f"{base_filename}_invvar.fits")
            invvar_hdu.writeto(invvar_filename, overwrite=True)
            print(f"Saved Inverse Variance Map to: {invvar_filename}")

            # Save Sigma Map
            sigma_hdu = fits.PrimaryHDU(sigma_map)
            sigma_hdu.header['HISTORY'] = 'Sigma (RMS Noise) Map'
            sigma_filename = os.path.join(output_dir, f"{base_filename}_sigma.fits")
            sigma_hdu.writeto(sigma_filename, overwrite=True)
            print(f"Saved Sigma Map to: {sigma_filename}")

    except requests.exceptions.RequestException as e:
        print(f"Error fetching data: {e}")
        raise e
    except (IOError, KeyError, IndexError) as e:
        print(f"Error processing FITS file: {e}")
        raise e

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Fetch and process data from the Legacy Survey.")
    parser.add_argument("ra", type=float, help="Right Ascension in degrees.")
    parser.add_argument("dec", type=float, help="Declination in degrees.")
    parser.add_argument("frb_name", type=str, help="The name of the FRB.")
    parser.add_argument("--output_dir", type=str, default="host_galaxies_data", help="Directory to save the output FITS files.")
    
    args = parser.parse_args()
    
    fetch_and_process_legacy_survey_data(args.ra, args.dec, args.frb_name, args.output_dir)
