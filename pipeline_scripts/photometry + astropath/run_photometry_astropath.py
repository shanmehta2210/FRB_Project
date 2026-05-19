import os
import sys
import yaml
import subprocess
import argparse

TEMPLATE_CONV = """CONV NORM
# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1
"""

TEMPLATE_NNW = """NNW
# Neural Network Weights for the SExtractor star/galaxy classifier (V1.3)
# inputs:	9 for profile parameters + 1 for seeing.
# outputs:	``Stellarity index'' (0.0 to 1.0)
# Seeing FWHM range: from 0.025 to 5.5'' (images must have 1.5 < FWHM < 5 pixels)
# Optimized for Moffat profiles with 2<= beta <= 4.

 3 10 10  1

-1.56604e+00 -2.48265e+00 -1.44564e+00 -1.24675e+00 -9.44913e-01 -5.22453e-01  4.61342e-02  8.31957e-01  2.15505e+00  2.64769e-01
 3.03477e+00  2.69561e+00  3.16188e+00  3.34497e+00  3.51885e+00  3.65570e+00  3.74856e+00  3.84541e+00  4.22811e+00  3.27734e+00

-3.22480e-01 -2.12804e+00  6.50750e-01 -1.11242e+00 -1.40683e+00 -1.55944e+00 -1.84558e+00 -1.18946e-01  5.52395e-01 -4.36564e-01 -5.30052e+00
 4.62594e-01 -3.29127e+00  1.10950e+00 -6.01857e-01  1.29492e-01  1.42290e+00  2.90741e+00  2.44058e+00 -9.19118e-01  8.42851e-01 -4.69824e+00
-2.57424e+00  8.96469e-01  8.34775e-01  2.18845e+00  2.46526e+00  8.60878e-02 -6.88080e-01 -1.33623e-02  9.30403e-02  1.64942e+00 -1.01231e+00
 4.81041e+00  1.53747e+00 -1.12216e+00 -3.16008e+00 -1.67404e+00 -1.75767e+00 -1.29310e+00  5.59549e-01  8.08468e-01 -1.01592e-02 -7.54052e+00
 1.01933e+01 -2.09484e+01 -1.07426e+00  9.87912e-01  6.05210e-01 -6.04535e-02 -5.87826e-01 -7.94117e-01 -4.89190e-01 -8.12710e-02 -2.07067e+01
-5.31793e+00  7.94240e+00 -4.64165e+00 -4.37436e+00 -1.55417e+00  7.54368e-01  1.09608e+00  1.45967e+00  1.62946e+00 -1.01301e+00  1.13514e-01
 2.20336e-01  1.70056e+00 -5.20105e-01 -4.28330e-01  1.57258e-03 -3.36502e-01 -8.18568e-02 -7.16163e+00  8.23195e+00 -1.71561e-02 -1.13749e+01
 3.75075e+00  7.25399e+00 -1.75325e+00 -2.68814e+00 -3.71128e+00 -4.62933e+00 -2.13747e+00 -1.89186e-01  1.29122e+00 -7.49380e-01  6.71712e-01
-8.41923e-01  4.64997e+00  5.65808e-01 -3.08277e-01 -1.01687e+00  1.73127e-01 -8.92130e-01  1.89044e+00 -2.75543e-01 -7.72828e-01  5.36745e-01
-3.65598e+00  7.56997e+00 -3.76373e+00 -1.74542e+00 -1.37540e-01 -5.55400e-01 -1.59195e-01  1.27910e-01  1.91906e+00  1.42119e+00 -4.35502e+00

-1.70059e+00 -3.65695e+00  1.22367e+00 -5.74367e-01 -3.29571e+00  2.46316e+00  5.22353e+00  2.42038e+00  1.22919e+00 -9.22250e-01 -2.32028e+00

 0.00000e+00 
 1.00000e+00 
"""

TEMPLATE_PARAM_PSF = """NUMBER
VIGNET(27,27)
X_IMAGE
Y_IMAGE
XWIN_IMAGE
YWIN_IMAGE
ERRAWIN_IMAGE
ERRBWIN_IMAGE
ALPHAWIN_J2000
DELTAWIN_J2000
FLUX_RADIUS
FWHM_WORLD
FLUX_AUTO
FLUXERR_AUTO
SNR_WIN
ELONGATION
FLUX_MAX
MAG_AUTO
MAGERR_AUTO
FLUX_APER(15)
FLUXERR_APER(15)
MAG_APER(15)
MAGERR_APER(15)
FLAGS
BACKGROUND
CLASS_STAR
FLAGS_MODEL
NITER_MODEL
FLUX_MODEL
FLUXERR_MODEL
MAG_MODEL
MAGERR_MODEL
XMODEL_IMAGE
YMODEL_IMAGE
FLUX_POINTSOURCE
FLUXERR_POINTSOURCE
MAG_POINTSOURCE
MAGERR_POINTSOURCE
FLUXRATIO_POINTSOURCE
FLUXRATIOERR_POINTSOURCE
AWIN_IMAGE
BWIN_IMAGE
THETAWIN_IMAGE
SPREAD_MODEL
SPREADERR_MODEL
"""

TEMPLATE_SEX_PSF = """# PSF Photometry configuration for SExtractor
CATALOG_NAME     {CATALOG_NAME}
CATALOG_TYPE     FITS_LDAC
PARAMETERS_NAME  photomPSF.param

DETECT_TYPE      CCD
DETECT_MINAREA   5
THRESH_TYPE      RELATIVE
DETECT_THRESH    {DETECT_THRESH}
ANALYSIS_THRESH  {ANALYSIS_THRESH}

FILTER           Y
FILTER_NAME      default.conv

DEBLEND_NTHRESH  32
DEBLEND_MINCONT  {DEBLEND_MINCONT}

CLEAN            Y
CLEAN_PARAM      1.0

MASK_TYPE        CORRECT

PSF_NAME         {PSF_NAME}
PSF_NMAX         1

WEIGHT_TYPE      {WEIGHT_TYPE}
WEIGHT_IMAGE     {WEIGHT_IMAGE}
WEIGHT_GAIN      Y

PHOT_APERTURES   4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 15.0, 20.0, 25.0, 30.0, {PHOT_APERTURES_PX}
PHOT_FLUXFRAC    0.5
PHOT_AUTOPARAMS  2.5, 3.5
PHOT_PETROPARAMS 2.0, 3.5

SATUR_KEY        SATURATE
SATUR_LEVEL      90000

MAG_ZEROPOINT    {MAG_ZEROPOINT}
MAG_GAMMA        4.0
GAIN             {GAIN}
PIXEL_SCALE      {PIXEL_SCALE}

SEEING_FWHM      {SEEING_FWHM}
STARNNW_NAME     default.nnw

BACK_SIZE        128
BACK_FILTERSIZE  3

BACK_TYPE        AUTO
BACKPHOTO_TYPE   LOCAL

CHECKIMAGE_TYPE  MODELS,-MODELS
CHECKIMAGE_NAME  psf_models.fits, psf_resi.fits

NTHREADS         1
MEMORY_OBJSTACK  10000
MEMORY_PIXSTACK  5000000
MEMORY_BUFSIZE   1024

VERBOSE_TYPE     NORMAL
WRITE_XML        N
"""


TEMPLATE_ASTROPHYSICS = """# -*- coding: utf-8 -*-
import os
import sys
import json
import yaml
import argparse
import numpy as np
import pandas as pd
import astropy.units as u
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from astroquery.vizier import Vizier
import matplotlib.pyplot as plt

# AstroPath dependency. Path is injected at template render time as an absolute
# (Linux/WSL) path so the WSL script does not depend on its working directory.
_ASTROPATH_PKG = "__ASTROPATH_PKG__"
if _ASTROPATH_PKG and _ASTROPATH_PKG not in sys.path:
    sys.path.insert(0, _ASTROPATH_PKG)
try:
    from astropath import path
except ImportError as exc:
    print(f"Error: AstroPath package not importable from '{_ASTROPATH_PKG}': {exc}")
    sys.exit(1)

def get_table_from_ldac(filename, frame=1):
    if frame > 0: frame = frame * 2
    return Table.read(filename, hdu=frame, format='fits')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--image", required=True)
    parser.add_argument("--psfcat", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--ra", type=float, default=None)
    parser.add_argument("--dec", type=float, default=None)
    args = parser.parse_args()

    # 1. READ CONFIG
    with open(args.config, "r") as f: config = yaml.safe_load(f)
    astio = config.get("astropath", {})
    err_a = astio.get("err_a_arcsec", 1.0)
    err_b = astio.get("err_b_arcsec", 1.0)
    err_pa = astio.get("err_theta_deg", 0.0)
    px_scale = config.get("sextractor_psf", {}).get("pixel_scale", 0.262)
    
    # 2. READ WCS TARGET
    with fits.open(args.image) as hdul:
        w = WCS(hdul[0].header).celestial
        x_shape = hdul[0].header['NAXIS1']
        y_shape = hdul[0].header['NAXIS2']
        if args.ra is None or args.dec is None:
            # Parse target from WCS center
            target_center = w.pixel_to_world(x_shape / 2, y_shape / 2)
        else:
            target_center = SkyCoord(args.ra, args.dec, unit='deg', frame='icrs')

    print(f"[*] Target FRB set to RA={target_center.ra.deg:.5f}, DEC={target_center.dec.deg:.5f}")

    # 3. ZERO-POINT CALIBRATION (PS1 with SkyMapper southern fallback)
    # All photometry is sourced from the with-PSF catalog (image.psf.cat). The 40-px
    # aperture is PSF-model-agnostic but is intentionally taken from the same catalog
    # as MAG_AUTO and MAG_POINTSOURCE so every magnitude has identical detection,
    # deblending, and flag context.
    #
    # Reference catalog selection — PS1 (II/349) is preferred when available:
    #   * full North + footprint to Dec ~ -30 deg
    #   * deeper than SkyMapper, more r-band detections per star (Nd > 6 filter)
    # SkyMapper DR1.1 (II/358) is the southern fallback:
    #   * covers the full southern hemisphere; overlaps PS1 down to ~ -30 deg
    #   * AB PSF magnitudes (rPSF) directly comparable to PS1 r
    #   * Vizier-native cone search, ~1 s round-trip (same order as PS1 query)
    # Both queries return an Astropy Table normalised to the schema
    #   ['RAJ2000', 'DEJ2000', 'rmag']
    # so the rest of the calibration code is catalog-agnostic.
    def _query_ps1(center):
        v = Vizier(columns=['*'], column_filters={"rmag": "< 20", "Nd": "> 6"}, row_limit=-1)
        Q = v.query_region(center, width=5.5*u.arcmin, height=13.6*u.arcmin, catalog='II/349')
        if not Q or len(Q[0]) == 0:
            return None, None
        return Q[0], "II/349 (PS1 DR1)"

    def _query_skymapper(center):
        # SkyMapper DR1.1: r-band PSF magnitudes; e_rPSF cut keeps only well-measured stars.
        v = Vizier(columns=['RAICRS', 'DEICRS', 'rPSF', 'e_rPSF'],
                   column_filters={"rPSF": "< 20", "e_rPSF": "< 0.05"}, row_limit=-1)
        Q = v.query_region(center, width=5.5*u.arcmin, height=13.6*u.arcmin, catalog='II/358')
        if not Q or len(Q[0]) == 0:
            return None, None
        t = Q[0]
        # Drop masked rows so sigma_clipped_stats sees a clean float column.
        if hasattr(t['rPSF'], 'mask'):
            t = t[~t['rPSF'].mask]
            if len(t) == 0:
                return None, None
        # Normalise columns to the PS1-style schema used downstream.
        t.rename_columns(['RAICRS', 'DEICRS', 'rPSF'], ['RAJ2000', 'DEJ2000', 'rmag'])
        return t, "II/358 (SkyMapper DR1.1)"

    print("[*] Querying PS1 for Zero-Points...")
    ps1_stars, cal_catalog_id = _query_ps1(target_center)
    if ps1_stars is None:
        print("[*] PS1 returned no stars (likely southern field). Falling back to SkyMapper DR1.1...")
        ps1_stars, cal_catalog_id = _query_skymapper(target_center)
    if ps1_stars is None:
        raise RuntimeError("No PS1 or SkyMapper calibration stars found in this field")
    print(f"[*] Calibration reference: {cal_catalog_id} ({len(ps1_stars)} candidate stars)")
    ps1_coords = SkyCoord(ra=ps1_stars['RAJ2000'], dec=ps1_stars['DEJ2000'], frame='icrs', unit='degree')

    cat_psf = get_table_from_ldac(args.psfcat)
    # Calibration stars: clean PSF flag, sub-2'' FWHM (exclude extended sources)
    cat_psf_cln = cat_psf[(cat_psf['FLAGS_MODEL'] == 0) & (cat_psf['FWHM_WORLD'] < 2.0/3600.0)]
    psfCoords = SkyCoord(ra=cat_psf_cln['ALPHAWIN_J2000'], dec=cat_psf_cln['DELTAWIN_J2000'], unit='deg')
    idx_pimg, idx_pps1, _, _ = ps1_coords.search_around_sky(psfCoords, 0.6*u.arcsec)
    if len(idx_pimg) < 3:
        raise RuntimeError(f"Too few PS1 calibration matches ({len(idx_pimg)}); cannot anchor ZP.")

    # 40-px aperture ZP — this is the ZP that goes into MAG_CALIB_APER_40PX (the magnitude
    # passed to AstroPath when mag_mode='mag_40px') and is the recommended production value.
    zp_med, _, zp_std = sigma_clipped_stats(
        ps1_stars['rmag'][idx_pps1] - cat_psf_cln['MAG_APER'][:, 14][idx_pimg], sigma=3)

    # PSF model ZP — for MAG_POINTSOURCE flux estimator (PSF-fit photometry)
    zp_p_med, _, zp_p_std = sigma_clipped_stats(
        ps1_stars['rmag'][idx_pps1] - cat_psf_cln['MAG_POINTSOURCE'][idx_pimg], sigma=3)

    # Kron / MAG_AUTO ZP — independent because Kron is a different flux estimator
    zp_auto_med, _, zp_auto_std = sigma_clipped_stats(
        ps1_stars['rmag'][idx_pps1] - cat_psf_cln['MAG_AUTO'][idx_pimg], sigma=3)

    print(f"[*] Calibration Complete (N_stars={len(idx_pimg)}): "
          f"40-px ZP = {zp_med:.3f}, PSF ZP = {zp_p_med:.3f}, Auto ZP = {zp_auto_med:.3f}")

    # Persist ZPs as a machine-readable artefact so master_run.py can hand them
    # off as part of the 'photometry' deliverable without parsing stdout.
    with open("zero_points.json", "w") as zf:
        json.dump({
            "n_calibration_stars": int(len(idx_pimg)),
            "zp_aper_40px": float(zp_med),
            "zp_aper_40px_std": float(zp_std),
            "zp_psf": float(zp_p_med),
            "zp_psf_std": float(zp_p_std),
            "zp_auto": float(zp_auto_med),
            "zp_auto_std": float(zp_auto_std),
            "filter_band": astio.get("filter_band", "r"),
            "reference_catalog": cal_catalog_id,
            "match_radius_arcsec": 0.6,
        }, zf, indent=2)

    # 4. FIND & MATCH CANDIDATES (single catalog, no cross-cat alignment needed)
    psfCoords_all = SkyCoord(ra=cat_psf['ALPHAWIN_J2000'], dec=cat_psf['DELTAWIN_J2000'], unit='deg')
    seps_all = target_center.separation(psfCoords_all).to(u.arcsec).value
    if len(seps_all) == 0:
        print("[!!] WARNING: cat_psf is empty — Phase 2 SExtractor produced no detections.")
    else:
        nearest_idx = int(np.argmin(seps_all))
        nearest_sep = float(seps_all[nearest_idx])
        if nearest_sep > 5.0:
            # Mirrors Phase 3a's diagnostic: a real host should sit within a
            # few arcsec of any sub-arcsec-localised FRB. Larger gaps usually
            # mean SExtractor missed a faint host or detect_thresh is too high.
            print(f"[!!] WARNING: closest cat_psf source is {nearest_sep:.2f} arcsec "
                  f"from FRB (cat_psf NUMBER={int(cat_psf['NUMBER'][nearest_idx])}, "
                  f"SNR_WIN={float(cat_psf['SNR_WIN'][nearest_idx]):.2f}). "
                  f"Real host may be below detect_thresh — double-check the field.")

    in_region = seps_all < 60.0   # 1 arcmin
    snr_min = astio.get("target_snr_min", 1.0)
    snr_mask = cat_psf['SNR_WIN'] > snr_min
    # FLAGS_MODEL is intentionally NOT used to gate candidates. Extended / deblended /
    # neighbor-affected galaxies commonly have FLAGS_MODEL > 0 and would otherwise be
    # silently dropped before AstroPath ever sees them. FLAGS_MODEL stays strict only
    # for the ZP star sample above.
    reg_srcs = cat_psf[in_region & snr_mask]
    n_in_region = int(np.sum(in_region))
    n_dropped_snr = int(np.sum(in_region & ~snr_mask))
    print(f"[*] Cross-matching {len(reg_srcs)} valid candidates in a 1 arcminute radius bounds "
          f"(within 1': {n_in_region}; dropped by SNR_WIN <= {snr_min}: {n_dropped_snr}).")

    records = []
    candidates_for_astropath = []

    for i, src in enumerate(reg_srcs):
        ra, dec = float(src['ALPHAWIN_J2000']), float(src['DELTAWIN_J2000'])

        # All three calibrated magnitudes from the same row in cat_psf
        mag_40 = float(src['MAG_APER'][14]) + float(zp_med)
        mag_psf = float(src['MAG_POINTSOURCE']) + float(zp_p_med)
        mag_auto = float(src['MAG_AUTO']) + float(zp_auto_med)
        ang_size = float(src['FLUX_RADIUS']) * float(px_scale)
        spread_model = float(src['SPREAD_MODEL'])
        spread_model_err = float(src['SPREADERR_MODEL'])

        mag_mode = config.get("sextractor_psf", {}).get("mag_mode", "mag_40px")
        if mag_mode == "mag_psf":
            mag_for_path = mag_psf
        elif mag_mode == "mag_auto":
            mag_for_path = mag_auto
        else:
            mag_for_path = mag_40

        # Star / point-source exclusion — SPREAD_MODEL with 3-sigma uncertainty inflation.
        #
        # SPREAD_MODEL ~ 0 for unresolved point sources; positive values indicate
        # extended emission. The 3*SPREADERR inflation guards against faint /
        # low-SNR galaxies whose nominal SPREAD_MODEL barely scrapes above the
        # cut despite having very large uncertainties — those were getting
        # classified as stars and filtered out before AstroPath ever saw them.
        # The cut is now: only label as a point source if SPREAD_MODEL is
        # confidently below 0.005 even after pushing it 3-sigma upward.
        # Negative SPREAD_MODEL values (cosmic rays etc.) still satisfy the
        # condition and are excluded.
        is_star = (spread_model + 3.0 * spread_model_err) < 0.005

        # Sanity-check the magnitude before letting AstroPath see it.
        # AstroPath's "inverse" prior is P_O ~ 1/Sigma_m where Sigma_m comes from
        # the Driver+2016 r-band number counts.  The power-law extrapolates
        # blindly, so any catalog source with mag ~99 / 131 / NaN / negative
        # flux gets a 1/Sigma_m of order 1e+12, completely swamping the real
        # hosts after renormalisation.  Restrict mag to the physically
        # plausible range Driver+2016 was fit on.
        bad_mag = (not np.isfinite(mag_for_path)) or (mag_for_path < 12.0) or (mag_for_path > 28.0)

        if bad_mag and not is_star:
            print(f"    src {i:3d}: SPREAD_MODEL={spread_model:.4f}+/-{spread_model_err:.4f}  mag={mag_for_path:.2f}  -> EXCLUDED (mag outside [12, 28] — corrupt photometry)")
        else:
            print(f"    src {i:3d}: SPREAD_MODEL={spread_model:.4f}+/-{spread_model_err:.4f}  mag={mag_for_path:.2f}  -> {'EXCLUDED (point src)' if is_star else 'INCLUDED'}")

        if (not is_star) and (not bad_mag):
            candidates_for_astropath.append({
                "objid": i, "ra": ra, "dec": dec, "mag": mag_for_path, "ang_size": ang_size, "source": mag_mode
            })

        records.append({
            "objid": i,
            "sex_number": int(src['NUMBER']),
            "RA": ra, "Dec": dec,
            "MAG_CALIB_APER_40PX": round(mag_40, 3),
            "MAG_CALIB_PSF": round(mag_psf, 3),
            "MAG_CALIB_AUTO": round(mag_auto, 3),
            "FLUX_RADIUS": float(src['FLUX_RADIUS']),
            "SPREAD_MODEL": round(spread_model, 5),
            "SPREADERR_MODEL": round(spread_model_err, 5),
            "included_in_astropath": (not is_star) and (not bad_mag)
        })
        
    df = pd.DataFrame(records)
    df.to_csv("calibrated_photometry_results.csv", index=False)

    if len(candidates_for_astropath) == 0:
        print("[!] No valid AstroPath candidates found.")
        sys.exit(0)

    # 5. ASTROPATH ASSOCIATION
    print(f"[*] Spinning up AstroPath Evaluation Engine for {len(candidates_for_astropath)} candidates...")
    mypath = path.PATH()
    mypath.init_localization("eellipse", center_coord=target_center, eellipse={"a": err_a, "b": err_b, "theta": err_pa})
    
    cdf = pd.DataFrame(candidates_for_astropath)
    mypath.init_candidates(ra=cdf["ra"].to_numpy(), dec=cdf["dec"].to_numpy(),
                           ang_size=cdf["ang_size"].to_numpy(), mag=cdf["mag"].to_numpy())
    
    # ========================================================================
    # ASTROPATH PRIOR CONFIGURATION
    # ========================================================================
    # All knobs that control which Bayesian prior set AstroPath evaluates live
    # in this single block.  Defaults below reproduce the "adopted" prior set
    # from Aggarwal et al. 2021 (PATH paper, ApJ 911 95), also defined in
    # astropath/priors.py::load_std_priors().  Edit the values to experiment;
    # the rest of the script does not need to be touched.
    #
    # ---- 1. CANDIDATE PRIOR  P(O_i)  --------------------------------------
    # Per-galaxy weight before the FRB position is taken into account.
    #   "inverse"   : raw P(O_i) ~ 1 / ang_size   (Aggarwal+2021 "adopted")
    #   "identical" : every candidate equally likely  (Aggarwal+2021 "conservative")
    #   "linear"    : raw P(O_i) ~ ang_size
    #   "user"      : custom function, set astropath.priors.USR_raw_prior_Oi
    # See astropath/priors.py::raw_prior_Oi.
    P_O_METHOD = "inverse"
    # P(U): prior probability that the true host is *unseen* (below detection
    # limit or outside the candidate list).  0.0 disables the unseen channel.
    P_U = astio.get("p_u", 0.1)

    # ---- 2. OFFSET PRIOR  P(theta | O_i)  ---------------------------------
    # Radial profile of expected FRB-to-host offsets, evaluated in units of
    # each candidate's angular half-light radius phi (= ang_size, arcsec).
    # See astropath/bayesian.py::pw_Oi for the exact functional forms.
    #
    #   "exp"     : P(theta) ~ exp(-theta / (THETA_SCALE * phi))   (default)
    #   "uniform" : P(theta) = const inside theta < THETA_MAX * phi, else 0
    #   "core"    : P(theta) ~ phi / (theta + phi)                 (cored)
    THETA_PDF = "exp"
    # THETA_MAX: truncation / integration radius, in units of phi (NOT arcsec).
    # Aggarwal+2021 default = 6.0 for both adopted and conservative sets.
    THETA_MAX = 6.0
    # THETA_SCALE: only used by "exp"; multiplies phi inside the exponential.
    # 1.0 -> e-folding length = phi (Aggarwal default).  Increase for a
    # broader, more permissive offset prior.
    THETA_SCALE = 1.0

    # ---- 3. POSTERIOR INTEGRATION (numerical, not statistical) ------------
    # AstroPath's px_Oi_local sets the integration grid step in arcsec to
    # `phi * POSTERIOR_STEP`.  When phi is large compared to the localization
    # sigma (e.g. phi = 8" host with 0.5" ASKAP localization), the default
    # step (0.1 * phi = 0.8") under-resolves the localization peak — the
    # entire FRB error ellipse fits inside a single grid cell and the host
    # likelihood collapses to ~zero.  Symptom: real, obvious hosts get
    # posterior_O ~ 1e-5 while P(U) sits at ~1.
    #
    # We therefore pick POSTERIOR_STEP adaptively so that the absolute grid
    # step satisfies  step_arcsec <= sigma_loc / 5  for every candidate, where
    # sigma_loc is the smaller localization semi-axis.  Floor at 0.005 to
    # keep the integration tractable for tiny localizations.
    POSTERIOR_METHOD = "local"   # "local" | "fixed"; "local" handles wide localisations
    POSTERIOR_RMAX   = 60.0      # max radius (arcsec) used for P(U) normalisation

    sigma_loc_arcsec = float(min(err_a, err_b))
    phi_max_arcsec = float(max(2.0, np.max(cdf["ang_size"].to_numpy())))
    POSTERIOR_STEP = max(0.005, min(0.1, sigma_loc_arcsec / (5.0 * phi_max_arcsec)))
    print(f"[*] AstroPath integration step_size = {POSTERIOR_STEP:.4f} (phi units); "
          f"~{POSTERIOR_STEP * phi_max_arcsec:.3f} arcsec at largest candidate "
          f"(phi_max={phi_max_arcsec:.2f} arcsec, sigma_loc_min={sigma_loc_arcsec:.2f} arcsec).")
    # ========================================================================

    mypath.init_cand_prior(P_O_method=P_O_METHOD, P_U=P_U)
    mypath.init_theta_prior(PDF=THETA_PDF, max=THETA_MAX, scale=THETA_SCALE)
    mypath.calc_priors()
    p_oix, p_ux = mypath.calc_posteriors(method=POSTERIOR_METHOD,
                                         step_size=POSTERIOR_STEP,
                                         max_radius=POSTERIOR_RMAX)

    cdf["posterior_O"] = p_oix
    cdf["posterior_U"] = p_ux
    best = cdf.sort_values("posterior_O", ascending=False).iloc[0]
    sep = target_center.separation(SkyCoord(ra=best['ra'], dec=best['dec'], unit='deg')).to(u.arcsec).value
    
    print(f"[*] ASTROPATH SUCCESS!")
    print(f"    Most Probable Host ObjID: {best['objid']}")
    print(f"    Posterior P(O): {best['posterior_O']:.4f}")
    print(f"    Unseen P(U): {p_ux:.4f}")
    print(f"    Separation: {sep:.2f} arcsec")

    cdf.to_csv("astropath_posteriors.csv", index=False)

    # 6. PLOT
    try:
        from astropy.visualization import ZScaleInterval
        fig, (ax_img, ax_scatter) = plt.subplots(1, 2, figsize=(16, 8))
        
        # Image Panel
        data_2d = np.squeeze(fits.getdata(args.image))
        while data_2d.ndim > 2: data_2d = data_2d[0]
        data_2d = np.asarray(data_2d, dtype=float)

        bx, by = w.world_to_pixel(SkyCoord(ra=best['ra'], dec=best['dec'], unit='deg'))
        cx, cy = w.world_to_pixel(target_center)
        px_per_arcmin = 60.0 / px_scale
        cx_min, cx_max = cx - px_per_arcmin, cx + px_per_arcmin
        cy_min, cy_max = cy - px_per_arcmin, cy + px_per_arcmin

        # Display stretch must be based on the *zoomed* field, not the full
        # cutout.  ZScaleInterval on a 10' image with a bright source far from
        # the FRB drives vmin/vmax so high that the central arcminute maps to a
        # flat black panel (seen e.g. on 20210320C) even though the data are fine.
        ny, nx = data_2d.shape
        ix0 = int(np.clip(np.floor(min(cx_min, cx_max)), 0, max(0, nx - 1)))
        ix1 = int(np.clip(np.ceil(max(cx_min, cx_max)), 0, nx))
        iy0 = int(np.clip(np.floor(min(cy_min, cy_max)), 0, max(0, ny - 1)))
        iy1 = int(np.clip(np.ceil(max(cy_min, cy_max)), 0, ny))
        if ix1 <= ix0:
            ix1 = min(ix0 + 1, nx)
        if iy1 <= iy0:
            iy1 = min(iy0 + 1, ny)
        crop = data_2d[iy0:iy1, ix0:ix1]
        crop_z = np.nan_to_num(crop, nan=0.0, posinf=0.0, neginf=0.0)
        interval = ZScaleInterval()
        try:
            vmin, vmax = interval.get_limits(crop_z)
        except Exception:
            vmin, vmax = np.nanpercentile(crop_z, [1.0, 99.0])
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin >= vmax:
            vmin = float(np.nanmin(crop_z))
            vmax = float(np.nanmax(crop_z))
            if not np.isfinite(vmin) or not np.isfinite(vmax) or vmin >= vmax:
                vmin, vmax = 0.0, 1.0

        ax_img.imshow(data_2d, origin='lower', cmap='bone', vmin=vmin, vmax=vmax)
        ax_img.scatter(bx, by, s=150, facecolors='none', edgecolors='red', lw=2, label=f"Host (P={best['posterior_O']:.2f})")
        ax_img.text(bx, by - 15, f"{best['posterior_O']:.2f}", color='red', fontsize=10, ha='center', va='top', weight='bold')
        
        ax_img.scatter(cx, cy, s=100, marker='+', color='cyan', label=f"FRB Center")
        
        # Plot all other candidates — skip any whose displayed value would be "0.00"
        for _, row in cdf.iterrows():
            if row['objid'] != best['objid'] and row['posterior_O'] >= 0.005:
                px, py = w.world_to_pixel(SkyCoord(ra=row['ra'], dec=row['dec'], unit='deg'))
                if cx_min <= px <= cx_max and cy_min <= py <= cy_max:
                    ax_img.scatter(px, py, s=50, facecolors='none', edgecolors='orange', alpha=0.6)
                    ax_img.text(px, py - 10, f"{row['posterior_O']:.2f}", color='orange', fontsize=8, ha='center', va='top')
        
        # Limit axis to strictly 1 arcminute radius around center
        ax_img.set_xlim(cx_min, cx_max)
        ax_img.set_ylim(cy_min, cy_max)
        
        ax_img.legend(loc='lower left')
        ax_img.set_title("AstroPath Highest Probability Host Association")
        
        # Scatter Panel
        sc = ax_scatter.scatter(cdf['mag'], cdf['posterior_O'], c=cdf['posterior_O'], cmap='viridis', s=80, edgecolor='k')
        ax_scatter.plot(best['mag'], best['posterior_O'], 'ro', markersize=15, markerfacecolor='none', markeredgewidth=2, label='Selected Host')
        plt.colorbar(sc, ax=ax_scatter, label='Posterior P(O|x)')
        ax_scatter.set_xlabel('Calibrated 40px Aperture Magnitude')
        ax_scatter.set_ylabel('Posterior Probability P(O|x)')
        ax_scatter.set_title("AstroPath Posterior vs Magnitude")
        ax_scatter.grid(True, alpha=0.4)
        ax_scatter.legend()
        
        plt.tight_layout()
        plt.savefig("astropath_association.png", dpi=300)
    except Exception as e:
        print(f"[!] Plotting error: {e}")

if __name__ == "__main__":
    main()
"""


def measure_psf_fwhm_arcsec(proto_path, pixel_scale):
    """Measure PSF FWHM in arcseconds by sampling directional profiles through the centroid.

    Profiles are sampled along 8 evenly-spaced angles (0 to 157.5 deg). For each angle the
    half-maximum crossing distance is found in both the + and - directions via linear
    interpolation, giving a per-direction FWHM. All valid directional FWHMs are averaged,
    making the estimate robust to PSF asymmetry.
    """
    try:
        import numpy as np
        from astropy.io import fits
        from scipy.ndimage import map_coordinates

        with fits.open(proto_path) as hdul:
            data = hdul[0].data
        data = np.squeeze(data)
        while data.ndim > 2:
            data = data[0]
        data = np.clip(data.astype(float), 0, None)
        total = data.sum()
        if total <= 0:
            return None

        ny, nx = data.shape
        y_idx, x_idx = np.mgrid[0:ny, 0:nx]
        cx = (x_idx * data).sum() / total
        cy = (y_idx * data).sum() / total

        peak = data.max()
        half_max = peak / 2.0

        # Sample radial distance up to the nearest image edge from the centroid
        max_r = min(cx, cy, nx - 1 - cx, ny - 1 - cy)
        r_samples = np.linspace(0, max_r, 400)

        N_dirs = 8  # sample every 22.5 degrees
        angles = np.linspace(0, np.pi, N_dirs, endpoint=False)

        fwhm_values = []
        for angle in angles:
            cos_a, sin_a = np.cos(angle), np.sin(angle)
            half_widths = []
            for sign in (+1, -1):
                xs = cx + sign * cos_a * r_samples
                ys = cy + sign * sin_a * r_samples
                valid = (xs >= 0) & (xs < nx) & (ys >= 0) & (ys < ny)
                if valid.sum() < 2:
                    continue
                profile = map_coordinates(data, [ys[valid], xs[valid]], order=1)
                below = profile < half_max
                if not below.any():
                    continue
                k = int(np.argmax(below))
                if k == 0:
                    continue
                # Linear interpolation for sub-pixel crossing distance
                r_valid = r_samples[valid]
                slope = profile[k] - profile[k - 1]
                if slope == 0:
                    continue
                r_cross = r_valid[k - 1] + (half_max - profile[k - 1]) / slope * (r_valid[k] - r_valid[k - 1])
                half_widths.append(r_cross)
            if len(half_widths) == 2:
                fwhm_values.append(half_widths[0] + half_widths[1])
            elif len(half_widths) == 1:
                fwhm_values.append(2.0 * half_widths[0])

        if not fwhm_values:
            return None

        fwhm_px = float(np.mean(fwhm_values))
        return fwhm_px * pixel_scale
    except Exception:
        return None


def _to_wsl_path(p):
    """Convert a Windows absolute path to a /mnt/<drive>/... WSL path; leave POSIX paths unchanged."""
    p_abs = os.path.abspath(p)
    if len(p_abs) > 2 and p_abs[1] == ':':
        return '/mnt/' + p_abs[0].lower() + '/' + p_abs[3:].replace('\\', '/')
    return p_abs.replace('\\', '/')


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_config = os.path.join(script_dir, "photometry_astropath_config.yaml")

    parser = argparse.ArgumentParser(description="Phase 2 Orchestrator: Photometry + AstroPath")
    parser.add_argument("--image", default="image.fits", help="Base FITS footprint (canonical name: image.fits)")
    parser.add_argument("--psf", default="image.psf", help="PSF model generated by Phase 1")
    parser.add_argument("--ra", type=float, default=None, help="Specific FRB Target RA. Leaves as WCS center if omitted.")
    parser.add_argument("--dec", type=float, default=None, help="Specific FRB Target DEC.")
    parser.add_argument("--config", default=default_config,
                        help="Path to photometry_astropath_config.yaml (default: alongside this script). "
                             "master_run.py writes a per-run config into the workdir and points here.")
    args = parser.parse_args()

    image_path = os.path.abspath(args.image)
    image_dir = os.path.dirname(image_path)
    image_name = os.path.basename(image_path)
    psf_catalog_name = f"{os.path.splitext(image_name)[0]}.psf.cat"

    config_path = os.path.abspath(args.config)
    if not os.path.exists(config_path):
        print(f"Error: Config file {config_path} not found.")
        sys.exit(1)
    print(f"[*] Phase 2 config: {config_path}")

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    # Build Template
    s_conf = config.get("sextractor_psf", {})
    
    # Check invvar for map weighting — honours use_weight_map from YAML
    use_weight_map = s_conf.get("use_weight_map", True)
    invvar_name = "invvar.fits"
    if use_weight_map and os.path.exists(os.path.join(image_dir, invvar_name)):
        weight_type = "MAP_WEIGHT"
        weight_image = invvar_name
        print("[*] Weight map: invvar.fits found and enabled.")
    else:
        weight_type = "NONE"
        weight_image = "NONE"
        if use_weight_map and not os.path.exists(os.path.join(image_dir, invvar_name)):
            print("[!] Warning: use_weight_map=true but invvar.fits not found. Running without weight map.")
        else:
            print("[*] Weight map disabled by config (use_weight_map=false).")

    # Measure actual SEEING_FWHM from proto_image.fits for reliable CLASS_STAR classification.
    # The PSFEx proto_image is the 25x25 PSF stamp — moment-based FWHM is more accurate than
    # the config default and gets injected directly into SExtractor's neural-network classifier.
    pixel_scale_cfg = s_conf.get("pixel_scale", 0.262)
    proto_path = os.path.join(image_dir, "proto_image.fits")
    seeing_fwhm = s_conf.get("seeing_fwhm", 2.0)  # fallback
    if os.path.exists(proto_path):
        measured = measure_psf_fwhm_arcsec(proto_path, pixel_scale_cfg)
        if measured:
            seeing_fwhm = measured
            print(f"[*] Measured SEEING_FWHM from proto_image.fits: {seeing_fwhm:.3f} arcsec")
        else:
            print(f"[!] FWHM measurement returned None, falling back to config: {seeing_fwhm} arcsec")
    else:
        print(f"[!] proto_image.fits not found in {image_dir}, using config SEEING_FWHM={seeing_fwhm}")

    deblend_mincont = float(s_conf.get("deblend_mincont", 0.005))

    sex_content = TEMPLATE_SEX_PSF.format(
        CATALOG_NAME=psf_catalog_name,
        DETECT_THRESH=s_conf.get("detect_thresh", 10),
        ANALYSIS_THRESH=s_conf.get("analysis_thresh", 10),
        DEBLEND_MINCONT=deblend_mincont,
        PSF_NAME=args.psf,
        WEIGHT_TYPE=weight_type,
        WEIGHT_IMAGE=weight_image,
        PHOT_APERTURES_PX=s_conf.get("phot_apertures_px", 40.0),
        MAG_ZEROPOINT=s_conf.get("mag_zeropoint", 0.0),
        GAIN=s_conf.get("gain", 1.6),
        PIXEL_SCALE=s_conf.get("pixel_scale", 0.0),
        SEEING_FWHM=seeing_fwhm,
    )

    path_sex = os.path.join(image_dir, "default_psf.sex")
    path_param = os.path.join(image_dir, "photomPSF.param")
    with open(path_sex, "w", newline="\n") as f: f.write(sex_content)
    with open(path_param, "w", newline="\n") as f: f.write(TEMPLATE_PARAM_PSF)
    
    # Ensure standards
    with open(os.path.join(image_dir, "default.conv"), "w", newline="\n") as f: f.write(TEMPLATE_CONV)
    with open(os.path.join(image_dir, "default.nnw"), "w", newline="\n") as f: f.write(TEMPLATE_NNW)

    print("[Phase 2] Executing Subprocess SExtractor with PSF Model")
    try:
        cmd_sex = ["wsl", "source-extractor", image_name, "-c", "default_psf.sex"]
        subprocess.run(cmd_sex, cwd=image_dir, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing PSF Photometry: {e}")
        sys.exit(1)

    print(f"[*] Phase 2 SExtractor DEBLEND_MINCONT = {deblend_mincont}")

    print("[Phase 2] Generating WSL Native Conda Script")
    # Resolve the AstroPath package path absolutely (relative to THIS orchestrator file)
    # so the WSL-side script no longer depends on cwd. The repo layout is:
    #   pipeline_scripts/photometry + astropath/run_photometry_astropath.py
    #   tools/AstroPath/astropath_pkg/
    script_dir = os.path.dirname(os.path.abspath(__file__))
    astropath_pkg_dir = os.path.abspath(os.path.join(script_dir, "..", "..", "tools", "AstroPath", "astropath_pkg"))
    astropath_pkg_wsl = _to_wsl_path(astropath_pkg_dir)
    if not os.path.isdir(astropath_pkg_dir):
        print(f"[!] Warning: AstroPath package directory not found at {astropath_pkg_dir} — WSL import will fail.")

    astro_script = os.path.join(image_dir, "_run_astrophysics_wsl.py")
    rendered = TEMPLATE_ASTROPHYSICS.replace("__ASTROPATH_PKG__", astropath_pkg_wsl)
    with open(astro_script, "w", encoding="utf-8", newline="\n") as f: f.write(rendered)

    print("[Phase 2] Triggering Conda `frb_project` Environment OS Bridge")
    coord_args = f"--ra {args.ra} --dec {args.dec}" if args.ra is not None else ""
    wsl_config_path = _to_wsl_path(config_path)
    wsl_bash_cmd = (
        f"conda activate frb_project && python _run_astrophysics_wsl.py "
        f"--image {image_name} --psfcat {psf_catalog_name} "
        f"--config '{wsl_config_path}' {coord_args}"
    )
    
    wsl_error = None
    try:
        subprocess.run(["wsl", "-e", "bash", "-ic", wsl_bash_cmd], cwd=image_dir, check=True)
    except subprocess.CalledProcessError as e:
        wsl_error = e
        print(f"Error inside the Astropath WSL bridge: {e}")
    finally:
        print("[Phase 2] Cleaning up templates")
        for temp_file in [path_sex, path_param, os.path.join(image_dir, "default.conv"), os.path.join(image_dir, "default.nnw"), astro_script]:
            if os.path.exists(temp_file):
                os.remove(temp_file)

    # Propagate the WSL failure so callers (e.g. master_run.py) see a non-zero
    # exit code instead of a false-positive "Phase 2 OK".
    if wsl_error is not None:
        sys.exit(1)

if __name__ == "__main__":
    main()
