import os
import sys
import yaml
import subprocess
import argparse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline_shared import (  # noqa: E402
    TEMPLATE_CONV,
    TEMPLATE_NNW,
    format_phot_apertures,
    get_logger,
    render_param_template,
    resolve_apertures,
)

log = get_logger("phase2")

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
FLUX_APER({NAPER})
FLUXERR_APER({NAPER})
MAG_APER({NAPER})
MAGERR_APER({NAPER})
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

PHOT_APERTURES   {PHOT_APERTURES}
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
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.stats import sigma_clipped_stats
from astroquery.vizier import Vizier
import matplotlib.pyplot as plt

# AstroPath dependency. Path is injected at template render time as an absolute
# (Linux/WSL) path so the WSL script does not depend on its working directory.
_ASTROPATH_PKG = "__ASTROPATH_PKG__"
if _ASTROPATH_PKG and _ASTROPATH_PKG not in sys.path:
    sys.path.insert(0, _ASTROPATH_PKG)
_DIAG_PKG = "__DIAG_PKG__"
if _DIAG_PKG and _DIAG_PKG not in sys.path:
    sys.path.insert(0, _DIAG_PKG)
try:
    from astropath import path
except ImportError as exc:
    print(f"Error: AstroPath package not importable from '{_ASTROPATH_PKG}': {exc}")
    sys.exit(1)
try:
    from field_depth import measure_field_depth
    from pipeline_diagnostics import add_sep_arcsec, plot_candidate_geometry
except ImportError as exc:
    print(f"[!] Diagnostic modules not importable from '{_DIAG_PKG}': {exc}")
    measure_field_depth = None
    add_sep_arcsec = None
    plot_candidate_geometry = None

# Production aperture, injected by the Phase 2 orchestrator from the resolved
# YAML aperture ladder (default: the largest aperture). Index is the column in
# MAG_APER / FLUX_APER; diameter (px) is used for the field-depth calculation.
_PROD_APER_INDEX = int("__PROD_APER_INDEX__")
_PROD_APER_DIAM = float("__PROD_APER_DIAM__")

def _query_with_retry(query_fn, max_retries=3, base_delay=2.0, label="query"):
    # Retry a query function with exponential backoff.
    import time
    for attempt in range(max_retries):
        try:
            result = query_fn()
            return result
        except Exception as e:
            delay = base_delay * (2 ** attempt)
            if attempt < max_retries - 1:
                print(f"[!] {label} attempt {attempt+1} failed: {e}. Retrying in {delay:.0f}s...")
                time.sleep(delay)
            else:
                print(f"[!] {label} failed after {max_retries} attempts: {e}")
                raise

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
    diag_cfg = config.get("diagnostics", {})
    do_field_depth = bool(diag_cfg.get("field_depth", True))
    do_geometry_plots = bool(diag_cfg.get("geometry_plots", True))
    err_a = astio.get("err_a_arcsec", 1.0)
    err_b = astio.get("err_b_arcsec", 1.0)
    err_pa = astio.get("err_theta_deg", 0.0)
    s_conf = config.get("sextractor_psf", {})
    phot_aperture_px = _PROD_APER_DIAM

    # Candidate search radius (arcsec) — configurable; legacy default 60" (1 arcmin).
    search_radius_arcsec = float(astio.get("search_radius_arcsec", 60.0))

    # Candidate prior and offset prior from config (fall back to Aggarwal+2021 adopted).
    cfg_p_o_method = astio.get("p_o_method", "inverse")
    cfg_theta_pdf = astio.get("theta_pdf", "exp")
    cfg_theta_max = float(astio.get("theta_max", 6.0))
    cfg_theta_scale = float(astio.get("theta_scale", 1.0))

    # 2. READ WCS TARGET + PIXEL SCALE
    with fits.open(args.image) as hdul:
        w = WCS(hdul[0].header).celestial
        x_shape = hdul[0].header['NAXIS1']
        y_shape = hdul[0].header['NAXIS2']
        if args.ra is None or args.dec is None:
            target_center = w.pixel_to_world(x_shape / 2, y_shape / 2)
        else:
            target_center = SkyCoord(args.ra, args.dec, unit='deg', frame='icrs')

    # Pixel scale from the WCS (arcsec/pixel). Always computed, never hardcoded.
    wcs_scales_deg = proj_plane_pixel_scales(w)
    px_scale = float(np.mean(wcs_scales_deg) * 3600.0)
    print(f"[*] Target FRB set to RA={target_center.ra.deg:.5f}, DEC={target_center.dec.deg:.5f}")
    print(f"[*] Pixel scale from WCS: {px_scale:.6f} arcsec/px")

    # Calibration catalogs cover the full image footprint (not a legacy 5.5' x 13.6' box).
    img_center = w.pixel_to_world(x_shape / 2.0, y_shape / 2.0)
    scales_deg = proj_plane_pixel_scales(w)
    cal_width_arcmin = float(x_shape * scales_deg[0] * 60.0)
    cal_height_arcmin = float(y_shape * scales_deg[1] * 60.0)
    print(f"[*] Calibration query footprint: {cal_width_arcmin:.2f}' x {cal_height_arcmin:.2f}' "
          f"({x_shape} x {y_shape} px, image center RA={img_center.ra.deg:.5f}, Dec={img_center.dec.deg:.5f})")

    # 3. ZERO-POINT CALIBRATION (PS1 + Legacy Survey; pick survey with more matches)
    # All photometry is sourced from the with-PSF catalog (image.psf.cat). The 40-px
    # aperture is PSF-model-agnostic but is intentionally taken from the same catalog
    # as MAG_AUTO and MAG_POINTSOURCE so every magnitude has identical detection,
    # deblending, and flag context.
    #
    # PS1 (Vizier) and Legacy DR10 Tractor (NOIRLab TAP) are queried every time.
    # After matching to clean PSF stars, we use whichever survey yields more matches
    # (>= MIN_CAL_STARS). Ties favour PS1.
    #
    # TODO: Add more survey options as backup reference catalogs (e.g., SDSS,
    # 2MASS, Gaia) for fields where neither PS1 nor Legacy has coverage.
    MIN_CAL_STARS = 3
    MATCH_RADIUS = 0.6 * u.arcsec

    def _query_ps1(center, width_arcmin, height_arcmin):
        v = Vizier(columns=['*'], column_filters={"rmag": "< 20", "Nd": "> 6"}, row_limit=-1)
        Q = _query_with_retry(
            lambda: v.query_region(
                center,
                width=width_arcmin * u.arcmin,
                height=height_arcmin * u.arcmin,
                catalog='II/349',
            ),
            label="PS1 Vizier",
        )
        if not Q or len(Q[0]) == 0:
            return None, None
        return Q[0], "II/349 (PS1 DR1)"

    def _query_legacy(center, width_arcmin, height_arcmin):
        try:
            import pyvo
        except ImportError:
            print("[!] pyvo not available; skipping Legacy Survey calibration query")
            return None, None
        ra_c = center.ra.deg
        dec_c = center.dec.deg
        width_deg = width_arcmin / 60.0
        height_deg = height_arcmin / 60.0
        ra_half = width_deg / 2.0
        dec_half = height_deg / 2.0
        cos_dec = np.cos(np.radians(np.clip(dec_c, -85, 85)))
        ra_range = ra_half / cos_dec
        ra_min = ra_c - ra_range
        ra_max = ra_c + ra_range
        dec_min = dec_c - dec_half
        dec_max = dec_c + dec_half
        if ra_min < 0:
            ra_clause = "(ra > {:.6f} OR ra < {:.6f})".format(ra_min + 360, ra_max)
        elif ra_max > 360:
            ra_clause = "(ra > {:.6f} OR ra < {:.6f})".format(ra_min, ra_max - 360)
        else:
            ra_clause = "ra > {:.6f} AND ra < {:.6f}".format(ra_min, ra_max)
        query = (
            "SELECT ra, dec, flux_r FROM ls_dr10.tractor WHERE "
            + ra_clause
            + " AND dec > {:.6f} AND dec < {:.6f}".format(dec_min, dec_max)
            + " AND type = 'PSF' AND fracflux_r < 0.05 AND anymask_r = 0 AND flux_r > 10"
        )
        service = pyvo.dal.TAPService("https://datalab.noirlab.edu/tap")
        raw = _query_with_retry(
            lambda: service.search(query).to_table(),
            label="Legacy TAP",
        )
        if raw is None or len(raw) == 0:
            return None, None
        flux = np.asarray(raw['flux_r'], dtype=float)
        ok = np.isfinite(flux) & (flux > 0)
        if not np.any(ok):
            return None, None
        rmag = 22.5 - 2.5 * np.log10(flux[ok])
        bright = rmag < 20.0
        if not np.any(bright):
            return None, None
        t = Table()
        t['RAJ2000'] = np.asarray(raw['ra'][ok][bright])
        t['DEJ2000'] = np.asarray(raw['dec'][ok][bright])
        t['rmag'] = rmag[bright]
        return t, "LS DR10 Tractor (NOIRLab TAP)"

    cat_psf = get_table_from_ldac(args.psfcat)
    cat_psf_cln = cat_psf[(cat_psf['FLAGS_MODEL'] == 0) & (cat_psf['FWHM_WORLD'] < 2.0/3600.0)]
    psfCoords = SkyCoord(ra=cat_psf_cln['ALPHAWIN_J2000'], dec=cat_psf_cln['DELTAWIN_J2000'], unit='deg')

    def _match_count(ref_table):
        if ref_table is None or len(ref_table) == 0:
            return 0, None, None
        ref_coords = SkyCoord(ra=ref_table['RAJ2000'], dec=ref_table['DEJ2000'], frame='icrs', unit='degree')
        idx_pimg, idx_pref, _, _ = ref_coords.search_around_sky(psfCoords, MATCH_RADIUS)
        return len(idx_pimg), idx_pimg, idx_pref

    print("[*] Querying PS1 for calibration stars...")
    try:
        ps1_tab, ps1_cat_id = _query_ps1(img_center, cal_width_arcmin, cal_height_arcmin)
    except Exception as exc:
        print(f"[!] PS1 Vizier query failed after all retries: {exc}. Continuing with n_ps1=0.")
        ps1_tab, ps1_cat_id = None, None
    n_ps1_cand = len(ps1_tab) if ps1_tab is not None else 0
    print(f"    PS1 candidates in field: {n_ps1_cand}")

    print("[*] Querying Legacy Survey DR10 (Tractor) for calibration stars...")
    try:
        leg_tab, leg_cat_id = _query_legacy(img_center, cal_width_arcmin, cal_height_arcmin)
    except Exception as exc:
        print(f"[!] Legacy TAP query failed after all retries: {exc}. Continuing with n_legacy=0.")
        leg_tab, leg_cat_id = None, None
    n_leg_cand = len(leg_tab) if leg_tab is not None else 0
    print(f"    Legacy candidates in field: {n_leg_cand}")

    n_ps1, idx_ps1_img, idx_ps1_ref = _match_count(ps1_tab)
    n_leg, idx_leg_img, idx_leg_ref = _match_count(leg_tab)
    match_arcsec = float(MATCH_RADIUS.to(u.arcsec).value)
    print(f"[*] Matched PSF stars within {match_arcsec:.1f} arcsec: "
          f"PS1={n_ps1}, Legacy={n_leg} (need >={MIN_CAL_STARS})")

    if n_ps1 >= MIN_CAL_STARS and n_ps1 >= n_leg:
        ref_stars, cal_catalog_id = ps1_tab, ps1_cat_id
        idx_pimg, idx_pref = idx_ps1_img, idx_ps1_ref
    elif n_leg >= MIN_CAL_STARS:
        ref_stars, cal_catalog_id = leg_tab, leg_cat_id
        idx_pimg, idx_pref = idx_leg_img, idx_leg_ref
    else:
        raise RuntimeError(
            "Too few calibration matches: PS1={}, Legacy={} "
            "(need >={} within {:.1f} arcsec).".format(
                n_ps1, n_leg, MIN_CAL_STARS, match_arcsec)
        )
    print(f"[*] Using {cal_catalog_id} for ZP ({len(idx_pimg)} matched stars)")

    # Production-aperture ZP — this is the ZP that goes into MAG_CALIB_APER (the
    # magnitude passed to AstroPath when mag_mode='mag_aper') and is the recommended
    # production value. sigma_clipped_stats returns (mean, median, std): we adopt the
    # sigma-clipped MEDIAN, which is robust against asymmetric reference-star outliers.
    _, zp_med, zp_std = sigma_clipped_stats(
        ref_stars['rmag'][idx_pref] - cat_psf_cln['MAG_APER'][:, _PROD_APER_INDEX][idx_pimg], sigma=3)

    # PSF model ZP — for MAG_POINTSOURCE flux estimator (PSF-fit photometry)
    _, zp_p_med, zp_p_std = sigma_clipped_stats(
        ref_stars['rmag'][idx_pref] - cat_psf_cln['MAG_POINTSOURCE'][idx_pimg], sigma=3)

    # Kron / MAG_AUTO ZP — independent because Kron is a different flux estimator
    _, zp_auto_med, zp_auto_std = sigma_clipped_stats(
        ref_stars['rmag'][idx_pref] - cat_psf_cln['MAG_AUTO'][idx_pimg], sigma=3)

    print(f"[*] Calibration Complete (N_stars={len(idx_pimg)}): "
          f"Aper ZP = {zp_med:.3f} ({_PROD_APER_DIAM:g}px), PSF ZP = {zp_p_med:.3f}, Auto ZP = {zp_auto_med:.3f}")

    zp_payload = {
        "n_calibration_stars": int(len(idx_pimg)),
        "n_ps1_matches": int(n_ps1),
        "n_legacy_matches": int(n_leg),
        "production_aperture_px": float(_PROD_APER_DIAM),
        "pixel_scale_arcsec_px": float(px_scale),
        "zp_aper": float(zp_med),
        "zp_aper_std": float(zp_std),
        "zp_psf": float(zp_p_med),
        "zp_psf_std": float(zp_p_std),
        "zp_auto": float(zp_auto_med),
        "zp_auto_std": float(zp_auto_std),
        "filter_band": astio.get("filter_band", "r"),
        "reference_catalog": cal_catalog_id,
        "match_radius_arcsec": match_arcsec,
    }
    if do_field_depth and measure_field_depth is not None:
        try:
            from pathlib import Path as _Path
            seg_p = _Path("segmentation_map.fits")
            proto_p = _Path("proto_image.fits")
            sex_p = _Path("default_psf.sex")
            fd = measure_field_depth(
                _Path(args.image),
                float(zp_med),
                seg_path=seg_p if seg_p.is_file() else None,
                proto_path=proto_p if proto_p.is_file() else None,
                sex_path=sex_p if sex_p.is_file() else None,
                aperture_diameter_px=float(phot_aperture_px),
            )
            zp_payload["field_depth"] = fd
            print(f"[*] Field depth (5σ @ {fd['aperture_diameter_px']} px): "
                  f"m_lim = {fd['m_lim_5sigma_ab']} AB")
        except Exception as exc:
            print(f"[!] Field depth measurement failed: {exc}")
    with open("zero_points.json", "w") as zf:
        json.dump(zp_payload, zf, indent=2)

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

    in_region = seps_all < search_radius_arcsec
    snr_min = astio.get("target_snr_min", 0.0)
    snr_mask = cat_psf['SNR_WIN'] > snr_min
    reg_srcs = cat_psf[in_region & snr_mask]
    n_in_region = int(np.sum(in_region))
    n_dropped_snr = int(np.sum(in_region & ~snr_mask))
    print(f"[*] Cross-matching {len(reg_srcs)} valid candidates within {search_radius_arcsec:.0f} arcsec "
          f"(in region: {n_in_region}; dropped by SNR_WIN <= {snr_min}: {n_dropped_snr}).")

    records = []
    candidates_for_astropath = []

    for i, src in enumerate(reg_srcs):
        ra, dec = float(src['ALPHAWIN_J2000']), float(src['DELTAWIN_J2000'])

        # All three calibrated magnitudes from the same row in cat_psf
        mag_40 = float(src['MAG_APER'][_PROD_APER_INDEX]) + float(zp_med)
        mag_psf = float(src['MAG_POINTSOURCE']) + float(zp_p_med)
        mag_auto = float(src['MAG_AUTO']) + float(zp_auto_med)
        ang_size = float(src['FLUX_RADIUS']) * float(px_scale)
        spread_model = float(src['SPREAD_MODEL'])
        spread_model_err = float(src['SPREADERR_MODEL'])

        mag_mode = config.get("sextractor_psf", {}).get("mag_mode", "mag_aper")
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
                "objid": i,
                "sex_number": int(src["NUMBER"]),
                "ra": ra,
                "dec": dec,
                "mag": mag_for_path,
                "ang_size": ang_size,
                "source": mag_mode,
            })

        records.append({
            "objid": i,
            "sex_number": int(src['NUMBER']),
            "RA": ra, "Dec": dec,
            "MAG_CALIB_APER": round(mag_40, 3),
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
    # All knobs are now YAML-configurable under the "astropath:" block in
    # photometry_astropath_config.yaml.  Defaults reproduce the "adopted"
    # prior set from Aggarwal et al. 2021 (PATH paper, ApJ 911 95).
    #
    # ---- 1. CANDIDATE PRIOR  P(O_i)  --------------------------------------
    # Per-galaxy weight before the FRB position is taken into account.
    #   "inverse"      : P(O) ~ 1/Sigma_m              (Aggarwal+2021 "adopted")
    #   "inverse_ang"  : P(O) ~ 1/(Sigma_m * R_eff)    (R_eff-weighted)
    #   "inverse_ang2" : P(O) ~ 1/(Sigma_m * R_eff^2)  (R_eff^2-weighted)
    #   "identical"    : all equal                      (Aggarwal+2021 "conservative")
    # See astropath/priors.py::raw_prior_Oi.
    P_O_METHOD = cfg_p_o_method
    P_U = astio.get("p_u", 0.1)

    # ---- 2. OFFSET PRIOR  P(theta | O_i)  ---------------------------------
    # Radial profile of expected FRB-to-host offsets, evaluated in units of
    # each candidate's angular half-light radius phi (= ang_size, arcsec).
    # See astropath/bayesian.py::pw_Oi for the exact functional forms.
    #
    #   "exp"     : P(theta) ~ exp(-theta / (scale * phi))   (default)
    #   "uniform" : P(theta) = const inside theta < max * phi, else 0
    #   "core"    : P(theta) ~ phi / (theta + phi)           (cored)
    #   "flat"    : P(theta) = const in arcsec within max * phi
    THETA_PDF = cfg_theta_pdf
    THETA_MAX = cfg_theta_max
    THETA_SCALE = cfg_theta_scale

    # ---- 3. POSTERIOR INTEGRATION (numerical, not statistical) ------------
    # px_Oi_local sets grid step in arcsec = phi * POSTERIOR_STEP. Adaptive
    # step ensures step_arcsec <= sigma_loc / 5 for all candidates, with a
    # floor at 0.005 and ceiling at 0.1 to keep grids tractable.
    POSTERIOR_METHOD = "local"
    POSTERIOR_RMAX = search_radius_arcsec

    sigma_loc_arcsec = float(min(err_a, err_b))
    phi_max_arcsec = float(max(2.0, np.max(cdf["ang_size"].to_numpy())))
    POSTERIOR_STEP = max(0.005, min(0.1, sigma_loc_arcsec / (5.0 * phi_max_arcsec)))
    print(f"[*] AstroPath config: P_O={P_O_METHOD}, theta_pdf={THETA_PDF}, "
          f"theta_max={THETA_MAX}, theta_scale={THETA_SCALE}")
    print(f"[*] AstroPath integration step_size = {POSTERIOR_STEP:.4f} (phi units); "
          f"~{POSTERIOR_STEP * phi_max_arcsec:.3f} arcsec at largest candidate "
          f"(phi_max={phi_max_arcsec:.2f} arcsec, sigma_loc_min={sigma_loc_arcsec:.2f} arcsec).")
    print(f"[*] Search radius = {search_radius_arcsec:.0f} arcsec (P(U) normalisation radius)")
    # ========================================================================

    mypath.init_cand_prior(P_O_method=P_O_METHOD, P_U=P_U)
    mypath.init_theta_prior(PDF=THETA_PDF, max=THETA_MAX, scale=THETA_SCALE)
    mypath.calc_priors()
    p_oix, p_ux = mypath.calc_posteriors(method=POSTERIOR_METHOD,
                                         step_size=POSTERIOR_STEP,
                                         max_radius=POSTERIOR_RMAX)

    cdf["posterior_O"] = p_oix
    cdf["posterior_U"] = p_ux
    if add_sep_arcsec is not None:
        cdf = add_sep_arcsec(cdf, target_center)
    best = cdf.sort_values("posterior_O", ascending=False).iloc[0]
    sep = float(best["sep_arcsec"]) if "sep_arcsec" in cdf.columns else target_center.separation(
        SkyCoord(ra=best['ra'], dec=best['dec'], unit='deg')).to(u.arcsec).value
    
    print(f"[*] ASTROPATH SUCCESS!")
    print(f"    Most Probable Host ObjID: {best['objid']} (SExtractor NUMBER={int(best['sex_number'])})")
    print(f"    Posterior P(O): {best['posterior_O']:.4f}")
    print(f"    Unseen P(U): {p_ux:.4f}")
    print(f"    Separation: {sep:.2f} arcsec")

    cdf.to_csv("astropath_posteriors.csv", index=False)

    if do_geometry_plots and plot_candidate_geometry is not None:
        try:
            plot_candidate_geometry(
                cdf, target_center, THETA_MAX, ".",
                best_objid=int(best["objid"]),
            )
            print("[*] Wrote sep_vs_shape_r.png and sep_vs_x_max_reff.png")
        except Exception as exc:
            print(f"[!] Geometry diagnostic plots failed: {exc}")

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
        px_per_search = search_radius_arcsec / px_scale
        cx_min, cx_max = cx - px_per_search, cx + px_per_search
        cy_min, cy_max = cy - px_per_search, cy + px_per_search

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
        ax_scatter.set_xlabel('Calibrated Production Aperture Magnitude')
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
        log.error(f"Config file {config_path} not found.")
        sys.exit(1)
    log.info(f"Phase 2 config: {config_path}")

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    # Build Template
    s_conf = config.get("sextractor_psf", {})

    # Aperture ladder + production aperture (YAML-configurable; default = largest).
    apertures, prod_index, prod_diam = resolve_apertures(s_conf)
    log.info(
        f"Apertures (px): {format_phot_apertures(apertures)} "
        f"| production = {prod_diam:g} px (index {prod_index})"
    )

    # Check invvar for map weighting — honours use_weight_map from YAML
    use_weight_map = s_conf.get("use_weight_map", True)
    invvar_name = "invvar.fits"
    if use_weight_map and os.path.exists(os.path.join(image_dir, invvar_name)):
        weight_type = "MAP_WEIGHT"
        weight_image = invvar_name
        log.info("Weight map: invvar.fits found and enabled.")
    else:
        weight_type = "NONE"
        weight_image = "NONE"
        if use_weight_map and not os.path.exists(os.path.join(image_dir, invvar_name)):
            log.warning("use_weight_map=true but invvar.fits not found. Running without weight map.")
        else:
            log.info("Weight map disabled by config (use_weight_map=false).")

    # Measure actual SEEING_FWHM from proto_image.fits for reliable CLASS_STAR classification.
    # The PSFEx proto_image is the 25x25 PSF stamp — moment-based FWHM is more accurate than
    # the config default and gets injected directly into SExtractor's neural-network classifier.
    # Compute pixel scale from WCS for the orchestrator side as well.
    try:
        from astropy.io import fits as _afits
        from astropy.wcs import WCS as _WCS
        from astropy.wcs.utils import proj_plane_pixel_scales as _ppps
        with _afits.open(image_path) as _hdul:
            _w = _WCS(_hdul[0].header).celestial
            pixel_scale_cfg = float(np.mean(_ppps(_w)) * 3600.0)
    except Exception:
        pixel_scale_cfg = s_conf.get("pixel_scale", 0.262)
    proto_path = os.path.join(image_dir, "proto_image.fits")
    seeing_fwhm = s_conf.get("seeing_fwhm", 2.0)  # fallback
    if os.path.exists(proto_path):
        measured = measure_psf_fwhm_arcsec(proto_path, pixel_scale_cfg)
        if measured:
            seeing_fwhm = measured
            log.info(f"Measured SEEING_FWHM from proto_image.fits: {seeing_fwhm:.3f} arcsec")
        else:
            log.warning(f"FWHM measurement returned None, falling back to config: {seeing_fwhm} arcsec")
    else:
        log.warning(f"proto_image.fits not found in {image_dir}, using config SEEING_FWHM={seeing_fwhm}")

    deblend_mincont = float(s_conf.get("deblend_mincont", 0.005))

    sex_content = TEMPLATE_SEX_PSF.format(
        CATALOG_NAME=psf_catalog_name,
        DETECT_THRESH=s_conf.get("detect_thresh", 3),
        ANALYSIS_THRESH=s_conf.get("analysis_thresh", 3),
        DEBLEND_MINCONT=deblend_mincont,
        PSF_NAME=args.psf,
        WEIGHT_TYPE=weight_type,
        WEIGHT_IMAGE=weight_image,
        PHOT_APERTURES=format_phot_apertures(apertures),
        MAG_ZEROPOINT=s_conf.get("mag_zeropoint", 0.0),
        GAIN=s_conf.get("gain", 1.6),
        PIXEL_SCALE=pixel_scale_cfg,
        SEEING_FWHM=seeing_fwhm,
    )

    path_sex = os.path.join(image_dir, "default_psf.sex")
    path_param = os.path.join(image_dir, "photomPSF.param")
    with open(path_sex, "w", newline="\n") as f: f.write(sex_content)
    with open(path_param, "w", newline="\n") as f:
        f.write(render_param_template(TEMPLATE_PARAM_PSF, len(apertures)))

    # Ensure standards
    with open(os.path.join(image_dir, "default.conv"), "w", newline="\n") as f: f.write(TEMPLATE_CONV)
    with open(os.path.join(image_dir, "default.nnw"), "w", newline="\n") as f: f.write(TEMPLATE_NNW)

    log.info("Executing Subprocess SExtractor with PSF Model")
    try:
        cmd_sex = ["wsl", "source-extractor", image_name, "-c", "default_psf.sex"]
        subprocess.run(cmd_sex, cwd=image_dir, check=True)
    except subprocess.CalledProcessError as e:
        log.error(f"Error executing PSF Photometry: {e}")
        sys.exit(1)

    log.info(f"Phase 2 SExtractor DEBLEND_MINCONT = {deblend_mincont}")

    log.info("Generating WSL Native Conda Script")
    # Resolve the AstroPath package path absolutely (relative to THIS orchestrator file)
    # so the WSL-side script no longer depends on cwd. The repo layout is:
    #   pipeline_scripts/photometry + astropath/run_photometry_astropath.py
    #   tools/AstroPath/astropath_pkg/
    script_dir = os.path.dirname(os.path.abspath(__file__))
    astropath_pkg_dir = os.path.abspath(os.path.join(script_dir, "..", "..", "tools", "AstroPath", "astropath_pkg"))
    astropath_pkg_wsl = _to_wsl_path(astropath_pkg_dir)
    diag_pkg_wsl = _to_wsl_path(script_dir)
    if not os.path.isdir(astropath_pkg_dir):
        log.warning(f"AstroPath package directory not found at {astropath_pkg_dir} — WSL import will fail.")

    astro_script = os.path.join(image_dir, "_run_astrophysics_wsl.py")
    rendered = (
        TEMPLATE_ASTROPHYSICS.replace("__ASTROPATH_PKG__", astropath_pkg_wsl)
        .replace("__DIAG_PKG__", diag_pkg_wsl)
        .replace("__PROD_APER_INDEX__", str(int(prod_index)))
        .replace("__PROD_APER_DIAM__", str(float(prod_diam)))
    )
    with open(astro_script, "w", encoding="utf-8", newline="\n") as f: f.write(rendered)

    log.info("Triggering Conda `frb_project` Environment OS Bridge")
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
        log.error(f"Error inside the Astropath WSL bridge: {e}")
    finally:
        log.info("Cleaning up templates")
        for temp_file in [path_sex, path_param, os.path.join(image_dir, "default.conv"), os.path.join(image_dir, "default.nnw"), astro_script]:
            if os.path.exists(temp_file):
                os.remove(temp_file)

    # Propagate the WSL failure so callers (e.g. master_run.py) see a non-zero
    # exit code instead of a false-positive "Phase 2 OK".
    if wsl_error is not None:
        sys.exit(1)

if __name__ == "__main__":
    main()
