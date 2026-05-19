"""
One-shot photometry calibration for the r70 field (SExtractor, PSFEx, PS1 + Legacy comparison).

Run from ``tools/Photometry`` (or pass ``--workdir``):

    python calibrate_photometry.py

The companion ``photometry.ipynb`` explains the steps; it does **not** invoke this script.
"""
import os
import sys
import subprocess
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from astroquery.vizier import Vizier
import pyvo
import math
import time

# Directory containing this file (``tools/Photometry``) — configs live next to the script.
PACKAGE_DIR = Path(__file__).resolve().parent


def resolve_photometry_workdir() -> str:
    """
    Return the Photometry working directory: the folder that contains ``photomCat.sex``.

    Resolution order:
    1. Directory of this script (so the pipeline works no matter what the CWD is).
    2. Current working directory, if ``photomCat.sex`` is there.
    3. Walk parents of CWD for ``tools/Photometry`` or legacy ``Photometry`` with that file.
    4. Fall back to ``PACKAGE_DIR``.
    """
    if (PACKAGE_DIR / "photomCat.sex").is_file():
        return str(PACKAGE_DIR)
    cwd = Path.cwd().resolve()
    if (cwd / "photomCat.sex").is_file():
        return str(cwd)
    for base in [cwd, *cwd.parents]:
        for leaf in (base / "tools" / "Photometry", base / "Photometry"):
            if (leaf / "photomCat.sex").is_file():
                return str(leaf.resolve())
    return str(PACKAGE_DIR)


def test_wsl_dependency(cmd):
    try:
        subprocess.run(f"wsl {cmd} --help", shell=True, capture_output=True, check=False)
        return True
    except Exception:
        return False

def get_table_from_ldac(filename, frame=1):
    """Load an astropy table from a fits_ldac by frame."""
    if frame > 0:
        frame = frame * 2
    tbl = Table.read(filename, hdu=frame, format='fits')
    return tbl

def _ra_clause(ra, dra):
    ra_min = ra - dra
    ra_max = ra + dra
    if ra_min < 0:
        return f"(ra > {ra_min + 360:.8f} OR ra < {ra_max:.8f})"
    if ra_max > 360:
        return f"(ra > {ra_min:.8f} OR ra < {ra_max - 360:.8f})"
    return f"ra > {ra_min:.8f} AND ra < {ra_max:.8f}"

def query_ls_source(svc, ra, dec, radius_arcsec=2.0):
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
    
    if flux_r <= 0: return None
    
    mag_ls = 22.5 - 2.5 * np.log10(flux_r)
    magerr_ls = 1.0857 / (flux_r * np.sqrt(flux_ivar_r)) if flux_ivar_r > 0 else np.nan

    return {
        "objid_ls": int(r["objid"]),
        "ra_ls": float(r["ra"]),
        "dec_ls": float(r["dec"]),
        "type_ls": str(r["type"]),
        "mag_ls": mag_ls,
        "magerr_ls": magerr_ls,
        "sep_arcsec_ls": float(sep[idx]),
    }

def run_calibration(
    workdir=None,
    image="coadded_astrometrically_corrected_rband_r70.fits",
    center="04h17m35.9058s +07d55m51.9812s",
    output="r70_target_comparison_photometry.csv",
):
    """
    Full SExtractor + PSFEx + zero-point calibration pipeline; writes comparison CSV.

    If ``workdir`` is None, uses :func:`resolve_photometry_workdir` (same folder as this
    script by default). Returns the absolute path to the written CSV.
    """
    workDir = os.path.abspath(workdir if workdir is not None else resolve_photometry_workdir())
    imageName = image
    imagePath = os.path.join(workDir, imageName)

    if not os.path.exists(imagePath):
        raise FileNotFoundError(f"Science image not found: {imagePath}")

    print(f"Working Directory: {workDir}")
    print(f"Processing Image: {imageName}")

    # Check dependencies
    if not test_wsl_dependency("source-extractor"):
        raise RuntimeError("source-extractor not found in WSL.")
    if not test_wsl_dependency("psfex"):
        raise RuntimeError("psfex not found in WSL.")

    # Load image and WCS
    with fits.open(imagePath) as hdul:
        data = hdul[0].data
        header = hdul[0].header
    w = WCS(header)

    target_center = SkyCoord(center, frame='icrs')
    ra_center, dec_center = target_center.ra.deg, target_center.dec.deg

    # Step 1: Query PS1 via Vizier for calibration
    print("Querying Pan-STARRS1 for calibration stars...")
    width_arcmin = 5.47
    height_arcmin = 13.68
    catNum = 'II/349'
    v = Vizier(columns=['*'], column_filters={"rmag": "< 20", "Nd": "> 6", "e_rmag": f"< {1.086/30:.3f}"}, row_limit=-1)
    try:
        Q = v.query_region(target_center, width=width_arcmin * u.arcmin, height=height_arcmin * u.arcmin, catalog=catNum, cache=False)
        if not Q:
            raise RuntimeError("No PS1 stars found.")
        ps1_stars = Q[0]
        print(f"Found {len(ps1_stars)} PS1 stars.")
    except Exception as e:
        print(f"Error querying PS1: {e}")
        sys.exit(1)

    ps1CatCoords = SkyCoord(ra=ps1_stars['RAJ2000'], dec=ps1_stars['DEJ2000'], frame='icrs', unit='degree')

    # Step 2: Aperture Photometry (Pass 1)
    print("Running SExtractor for aperture photometry...")
    catalogName = imageName + ".cat"
    configFile = "photomCat.sex"
    paramFile = "photomCat.param"
    
    cmd = f'wsl --cd "{workDir}" source-extractor -c {configFile} {imageName} -CATALOG_NAME {catalogName} -PARAMETERS_NAME {paramFile}'
    subprocess.run(cmd, shell=True, check=True, capture_output=True)

    sourceTable = get_table_from_ldac(os.path.join(workDir, catalogName))
    cleanSources = sourceTable[(sourceTable['FLAGS'] == 0) & (sourceTable['FWHM_WORLD'] < 2.0/3600.0)]
    print(f"Detected {len(cleanSources)} clean sources.")

    # Step 3: Deriving 13-px Aperture ZP
    sourceCatCoords = SkyCoord(ra=cleanSources['ALPHAWIN_J2000'], dec=cleanSources['DELTAWIN_J2000'], frame='icrs', unit='degree')
    idx_image, idx_ps1, d2d, _ = ps1CatCoords.search_around_sky(sourceCatCoords, 0.6 * u.arcsec)
    
    aperture_diameter = np.array([4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 15.0, 20.0, 25.0, 30.0, 40.0])
    largest_aper_idx = 9 # 13px reference
    aper_40px_idx = 14 # 40px
    
    offsets = ps1_stars['rmag'][idx_ps1] - cleanSources['MAG_APER'][:, largest_aper_idx][idx_image]
    zp_mean, zp_med, zp_std = sigma_clipped_stats(offsets, sigma=3, maxiters=10)
    print(f"13px Aperture Zero Point: {zp_med:.4f} +/- {zp_std:.4f}")

    # Step 4: Run PSFEx
    print("Running PSFEx...")
    psfConfigFile = 'psfex_conf.psfex'
    cmd = f'wsl --cd "{workDir}" psfex -c {psfConfigFile} {catalogName}'
    subprocess.run(cmd, shell=True, check=True, capture_output=True)

    # Step 5: PSF Photometry (Pass 2)
    print("Running SExtractor with PSF model...")
    psfName = imageName + ".psf"
    psfCatalogName = imageName + ".psf.cat"
    psfParamFile = "photomPSF.param"
    psfSexFile = "PSFPhotom.sex"

    cmd = (f'wsl --cd "{workDir}" source-extractor -c {psfSexFile} {imageName} '
           f'-CATALOG_NAME {psfCatalogName} -PSF_NAME {psfName} '
           f'-PARAMETERS_NAME {psfParamFile} -MAG_ZEROPOINT 0')
    subprocess.run(cmd, shell=True, check=True, capture_output=True)

    psfSourceTable = get_table_from_ldac(os.path.join(workDir, psfCatalogName))
    cleanPSFSources = psfSourceTable[(psfSourceTable['FLAGS_MODEL'] == 0) & (psfSourceTable['FWHM_WORLD'] < 2.0/3600.0)]
    
    psfSourceCatCoords = SkyCoord(ra=cleanPSFSources['ALPHAWIN_J2000'], dec=cleanPSFSources['DELTAWIN_J2000'], frame='icrs', unit='degree')
    idx_psfimage, idx_psfps1, d2d_psf, _ = ps1CatCoords.search_around_sky(psfSourceCatCoords, 0.6 * u.arcsec)

    psf_offsets = ps1_stars['rmag'][idx_psfps1] - cleanPSFSources['MAG_POINTSOURCE'][idx_psfimage]
    zp_psf_mean, zp_psf_med, zp_psf_std = sigma_clipped_stats(psf_offsets, sigma=3, maxiters=10)
    print(f"PSF Zero Point: {zp_psf_med:.4f} +/- {zp_psf_std:.4f}")

    # Step 6: Query Legacy Surveys TAP for target region
    print("Querying Legacy Surveys TAP for target region...")
    svc = pyvo.dal.TAPService("https://datalab.noirlab.edu/tap")
    
    # Filter regional sources from PSF pass (more accurate)
    # Applying SNR > 10 filter as per user request
    separations = target_center.separation(psfSourceCatCoords)
    in_region = separations < (0.5 * u.arcmin)
    region_sources = cleanPSFSources[in_region & (cleanPSFSources['SNR_WIN'] > 10)]
    
    print(f"Found {len(region_sources)} sources with SNR > 10 in target region.")
    
    # We also need aperture mags for these same sources, let's match them from the aperture table
    region_coords = SkyCoord(ra=region_sources['ALPHAWIN_J2000'], dec=region_sources['DELTAWIN_J2000'], unit=u.deg)
    idx_aper, _, _ = region_coords.match_to_catalog_sky(sourceCatCoords)
    
    comparison_data = []
    for i, src in enumerate(region_sources):
        ra, dec = src['ALPHAWIN_J2000'], src['DELTAWIN_J2000']
        print(f"Processing target source {i+1}/{len(region_sources)} at {ra:.5f}, {dec:.5f}...")
        ls_match = query_ls_source(svc, ra, dec)
        
        # Calibrated mags and errors (Quadrature sum)
        calib_mag_psf = src['MAG_POINTSOURCE'] + zp_psf_med
        mag_err_psf_instr = src['MAGERR_POINTSOURCE']
        mag_err_psf_calib = np.sqrt(mag_err_psf_instr**2 + zp_psf_std**2)

        # Get matching aperture mags and errors
        matched_aper = cleanSources[idx_aper[i]]
        
        # 13px Reference
        calib_mag_aper_13 = matched_aper['MAG_APER'][largest_aper_idx] + zp_med
        mag_err_aper_13_instr = matched_aper['MAGERR_APER'][largest_aper_idx]
        mag_err_aper_13_calib = np.sqrt(mag_err_aper_13_instr**2 + zp_std**2)

        # 40px Large Aperture
        calib_mag_aper_40 = matched_aper['MAG_APER'][aper_40px_idx] + zp_med
        mag_err_aper_40_instr = matched_aper['MAGERR_APER'][aper_40px_idx]
        mag_err_aper_40_calib = np.sqrt(mag_err_aper_40_instr**2 + zp_std**2)

        entry = {
            'RA': round(float(ra), 6),
            'Dec': round(float(dec), 6),
            'X_IMAGE': round(float(src['XMODEL_IMAGE']), 2),
            'Y_IMAGE': round(float(src['YMODEL_IMAGE']), 2),
            'ZERO_PT_APER_13PX': round(float(zp_med), 4),
            'ZERO_PT_PSF': round(float(zp_psf_med), 4),
            'MAG_CALIB_APER_13PX': round(float(calib_mag_aper_13), 3),
            'MAGERR_CALIB_APER_13PX': round(float(mag_err_aper_13_calib), 3),
            'MAG_CALIB_APER_40PX': round(float(calib_mag_aper_40), 3),
            'MAGERR_CALIB_APER_40PX': round(float(mag_err_aper_40_calib), 3),
            'MAG_CALIB_PSF': round(float(calib_mag_psf), 3),
            'MAGERR_CALIB_PSF': round(float(mag_err_psf_calib), 3),
            'FLUX_RADIUS': round(float(src['FLUX_RADIUS']), 2),
            'SNR': round(float(src['SNR_WIN']), 2),
            'FLAGS': int(src['FLAGS'])
        }
        
        if ls_match:
            # Rounding LS values for clean output
            entry.update({
                "objid_ls": ls_match["objid_ls"],
                "ra_ls": round(float(ls_match["ra_ls"]), 6),
                "dec_ls": round(float(ls_match["dec_ls"]), 6),
                "type_ls": ls_match["type_ls"],
                "mag_ls": round(float(ls_match["mag_ls"]), 3),
                "magerr_ls": round(float(ls_match["magerr_ls"]), 3) if not np.isnan(ls_match["magerr_ls"]) else np.nan,
                "sep_arcsec_ls": round(float(ls_match["sep_arcsec_ls"]), 3)
            })
            if not np.isnan(ls_match["magerr_ls"]) and ls_match["magerr_ls"] > 0:
                entry["ls_sigma_dist_aper_13"] = round(float((calib_mag_aper_13 - ls_match["mag_ls"]) / ls_match["magerr_ls"]), 2)
                entry["ls_sigma_dist_aper_40"] = round(float((calib_mag_aper_40 - ls_match["mag_ls"]) / ls_match["magerr_ls"]), 2)
                entry["ls_sigma_dist_psf"] = round(float((calib_mag_psf - ls_match["mag_ls"]) / ls_match["magerr_ls"]), 2)
            else:
                entry["ls_sigma_dist_psf"] = np.nan
        else:
            entry["objid_ls"] = -1
            entry["ls_sigma_dist_psf"] = np.nan
            
        comparison_data.append(entry)

    df = pd.DataFrame(comparison_data)
    outputPath = os.path.join(workDir, output)
    df.to_csv(outputPath, index=False)
    print(f"Results saved to {outputPath}")
    return outputPath


def main():
    parser = argparse.ArgumentParser(description="Automated PSF and Aperture Photometry Calibration")
    parser.add_argument("--image", default="coadded_astrometrically_corrected_rband_r70.fits", help="Science image FITS file")
    parser.add_argument("--center", default="04h17m35.9058s +07d55m51.9812s", help="Target center (RA Dec string)")
    parser.add_argument(
        "--workdir",
        default=None,
        help="Working directory (default: directory containing this script, or auto-discovered Photometry folder)",
    )
    parser.add_argument("--output", default="r70_target_comparison_photometry.csv", help="Output CSV filename")
    args = parser.parse_args()
    try:
        run_calibration(
            workdir=args.workdir,
            image=args.image,
            center=args.center,
            output=args.output,
        )
    except (FileNotFoundError, RuntimeError) as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
