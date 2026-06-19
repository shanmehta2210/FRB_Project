import os
import sys
import yaml
import subprocess
import argparse

# Shared pipeline helpers (CONV/NNW templates, aperture resolution, logging).
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline_shared import (  # noqa: E402
    TEMPLATE_CONV,
    TEMPLATE_NNW,
    format_phot_apertures,
    get_logger,
    render_param_template,
    resolve_apertures,
)

log = get_logger("phase1")

TEMPLATE_PARAM = """NUMBER
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
FLAGS
BACKGROUND
CLASS_STAR
FLUX_APER({NAPER})
FLUXERR_APER({NAPER})
MAG_APER({NAPER})
MAGERR_APER({NAPER})
A_IMAGE
B_IMAGE
THETA_IMAGE
AWIN_IMAGE
BWIN_IMAGE
THETAWIN_IMAGE
ERRTHETAWIN_IMAGE
"""

TEMPLATE_SEX = """# Default configuration file for SExtractor
CATALOG_NAME     {CATALOG_NAME}
CATALOG_TYPE     FITS_LDAC
PARAMETERS_NAME  default.param

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

WEIGHT_TYPE      {WEIGHT_TYPE}
WEIGHT_IMAGE     {WEIGHT_IMAGE}

PHOT_APERTURES   {PHOT_APERTURES}
PHOT_FLUXFRAC    0.5
PHOT_AUTOPARAMS  1.0,2.0
PHOT_PETROPARAMS 1.0,2.0

SATUR_KEY        SATURATE
SATUR_LEVEL      90000

MAG_ZEROPOINT    {MAG_ZEROPOINT}
MAG_GAMMA        4.0
GAIN             {GAIN}
PIXEL_SCALE      {PIXEL_SCALE}

SEEING_FWHM      {SEEING_FWHM}
STARNNW_NAME     default.nnw

BACK_SIZE        1024
BACK_FILTERSIZE  3

BACK_TYPE        AUTO
BACKPHOTO_TYPE   LOCAL

CHECKIMAGE_TYPE  SEGMENTATION
CHECKIMAGE_NAME  segmentation_map.fits

NTHREADS         1
MEMORY_OBJSTACK  10000
MEMORY_PIXSTACK  5000000
MEMORY_BUFSIZE   1024

VERBOSE_TYPE     NORMAL
WRITE_XML        N
"""

TEMPLATE_PSFEX = """# Default configuration file for PSFEx
BASIS_TYPE      PIXEL
BASIS_NUMBER    11
PSF_SAMPLING    {PSF_SAMPLING}
PSF_ACCURACY    0.001
PSF_SIZE        25,25
CENTER_KEYS     X_IMAGE,Y_IMAGE
PHOTFLUX_KEY    FLUX_APER(1)
PHOTFLUXERR_KEY FLUXERR_APER(1)

PSF_RECENTER    Y

PSFVAR_KEYS     X_IMAGE,Y_IMAGE
PSFVAR_GROUPS   1,1
PSFVAR_DEGREES  {PSFVAR_DEGREES}
PSFVAR_NSNAP    1

SAMPLE_AUTOSELECT  Y
SAMPLEVAR_TYPE     SEEING
SAMPLE_FWHMRANGE   {SAMPLE_FWHMRANGE}
SAMPLE_VARIABILITY {SAMPLE_VARIABILITY}
SAMPLE_MINSN       {SAMPLE_MINSN}
SAMPLE_MAXELLIP    {SAMPLE_MAXELLIP}

HOMOBASIS_TYPE     GAUSS-LAGUERRE
HOMOBASIS_NUMBER   10
HOMOBASIS_SCALE    1.0
HOMOPSF_PARAMS     2.0, 3.0
HOMOKERNEL_DIR     
HOMOKERNEL_SUFFIX  .homo.fits

OUTCAT_TYPE        ASCII_HEAD
OUTCAT_NAME        psfex_out.cat

CHECKPLOT_DEV       NULL
CHECKPLOT_TYPE      NONE
CHECKPLOT_NAME      NONE

CHECKIMAGE_TYPE PROTOTYPES
CHECKIMAGE_NAME  proto.fits
CHECKIMAGE_CUBE N

PSF_DIR         
PSF_SUFFIX      .psf
VERBOSE_TYPE    QUIET
WRITE_XML       Y
XML_NAME        psfex.xml
NTHREADS        0
"""


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_config = os.path.join(script_dir, "pipeline_config.yaml")

    parser = argparse.ArgumentParser(description="Phase 1 / SExtractor + PSFEx orchestrator")
    parser.add_argument("image", help="Path to flux FITS image (canonical name: image.fits)")
    parser.add_argument("--config", default=default_config,
                        help="Path to pipeline_config.yaml (default: alongside this script). "
                             "master_run.py writes a per-run config into the workdir and points here.")
    parser.add_argument("--keep-templates", action="store_true",
                        help="Keep the generated default.{sex,param,conv,nnw,psfex} files after the "
                             "run (default: remove them; master_run deletes the whole workdir anyway).")
    args = parser.parse_args()

    image_path = os.path.abspath(args.image)
    if not os.path.exists(image_path):
        log.error(f"Image {image_path} not found.")
        sys.exit(1)

    image_dir = os.path.dirname(image_path)
    image_name = os.path.basename(image_path)
    # Catalog name is fixed so PSFEx writes deterministic image.psf and proto_image.fits
    # regardless of the user-supplied input image filename. Phase 2 / Phase 3 contracts
    # rely on those exact filenames.
    catalog_name = "image.cat"

    config_path = os.path.abspath(args.config)
    if not os.path.exists(config_path):
        log.error(f"Config file {config_path} not found.")
        sys.exit(1)
    log.info(f"Phase 1 config: {config_path}")

    with open(config_path, "r") as f:
        config = yaml.safe_load(f)

    # Read Sextractor Configs
    s_conf = config.get("sextractor", {})
    detect_thresh = s_conf.get("detect_thresh", 3)
    analysis_thresh = s_conf.get("analysis_thresh", 3)
    gain = s_conf.get("gain", 1.6)
    mag_zeropoint = s_conf.get("mag_zeropoint", 0.0)
    pixel_scale = s_conf.get("pixel_scale")
    if pixel_scale is None:
        try:
            from astropy.io import fits as _afits
            from astropy.wcs import WCS as _WCS
            from astropy.wcs.utils import proj_plane_pixel_scales as _ppps
            import numpy as _np
            with _afits.open(image_path) as _hdul:
                _w = _WCS(_hdul[0].header).celestial
                pixel_scale = float(_np.mean(_ppps(_w)) * 3600.0)
            log.info(f"Pixel scale from WCS: {pixel_scale:.6f} arcsec/px")
        except Exception:
            pixel_scale = 0.0
            log.warning("Could not compute pixel scale from WCS; using SExtractor default 0.0")
    seeing_fwhm = s_conf.get("seeing_fwhm", 2.0)
    deblend_mincont = float(s_conf.get("deblend_mincont", 0.005))

    # Aperture ladder + production aperture (YAML-configurable; default = largest).
    apertures, prod_index, prod_diam = resolve_apertures(s_conf)
    log.info(
        f"Apertures (px): {format_phot_apertures(apertures)} "
        f"| production = {prod_diam:g} px (index {prod_index})"
    )

    # Read PSFEx Configs
    p_conf = config.get("psfex", {})
    psf_sampling = p_conf.get("psf_sampling", 1.0)
    psfvar_degrees = p_conf.get("psfvar_degrees", 0)
    sample_fwhm_range = p_conf.get("sample_fwhm_range", "2.0,10")
    sample_variability = p_conf.get("sample_variability", 0.4)
    sample_minsn = p_conf.get("sample_minsn", 30)
    sample_max_ellp = p_conf.get("sample_max_ellp", 0.2)
    min_accepted_stars = int(p_conf.get("min_accepted_stars", 25))

    # Render Templates
    param_content = render_param_template(TEMPLATE_PARAM, len(apertures))

    # Check for weight map — honours use_weight_map from YAML first
    use_weight_map = s_conf.get("use_weight_map", True)
    invvar_name = "invvar.fits"
    invvar_path = os.path.join(image_dir, invvar_name)
    if use_weight_map and os.path.exists(invvar_path):
        weight_type = "MAP_WEIGHT"
        weight_image = invvar_name
        log.info("Weight map: invvar.fits found and enabled.")
    else:
        weight_type = "NONE"
        weight_image = "NONE"
        if use_weight_map and not os.path.exists(invvar_path):
            log.warning("use_weight_map=true but invvar.fits not found. Running without weight map.")
        else:
            log.info("Weight map disabled by config (use_weight_map=false).")

    sex_content = TEMPLATE_SEX.format(
        CATALOG_NAME=catalog_name,
        DETECT_THRESH=detect_thresh,
        ANALYSIS_THRESH=analysis_thresh,
        PHOT_APERTURES=format_phot_apertures(apertures),
        MAG_ZEROPOINT=mag_zeropoint,
        GAIN=gain,
        PIXEL_SCALE=pixel_scale,
        SEEING_FWHM=seeing_fwhm,
        DEBLEND_MINCONT=deblend_mincont,
        WEIGHT_TYPE=weight_type,
        WEIGHT_IMAGE=weight_image
    )

    psfex_content = TEMPLATE_PSFEX.format(
        PSF_SAMPLING=psf_sampling,
        PSFVAR_DEGREES=psfvar_degrees,
        SAMPLE_FWHMRANGE=sample_fwhm_range,
        SAMPLE_VARIABILITY=sample_variability,
        SAMPLE_MINSN=sample_minsn,
        SAMPLE_MAXELLIP=sample_max_ellp
    )

    # Write to target directory
    path_conv = os.path.join(image_dir, "default.conv")
    path_nnw = os.path.join(image_dir, "default.nnw")
    path_param = os.path.join(image_dir, "default.param")
    path_sex = os.path.join(image_dir, "default.sex")
    path_psfex = os.path.join(image_dir, "default.psfex")

    with open(path_conv, "w", newline="\n") as f: f.write(TEMPLATE_CONV)
    with open(path_nnw, "w", newline="\n") as f: f.write(TEMPLATE_NNW)
    with open(path_param, "w", newline="\n") as f: f.write(param_content)
    with open(path_sex, "w", newline="\n") as f: f.write(sex_content)
    with open(path_psfex, "w", newline="\n") as f: f.write(psfex_content)

    try:
        log.info("Running SExtractor...")
        # Note: image_name is used, we run inside image_dir
        cmd_sex = ["wsl", "source-extractor", image_name, "-c", "default.sex"]
        subprocess.run(cmd_sex, cwd=image_dir, check=True)

        # Iterative PSFEx Strategy
        snrs_to_try = [30, 20, 10, 5, 3, 2.5, 2]
        import xml.etree.ElementTree as ET

        for snr in snrs_to_try:
            log.info(f"Executing PSFEx with SAMPLE_MINSN = {snr}")
            # Render dynamic PSFEX for this loop iteration
            psfex_content = TEMPLATE_PSFEX.format(
                PSF_SAMPLING=psf_sampling,
                PSFVAR_DEGREES=psfvar_degrees,
                SAMPLE_FWHMRANGE=sample_fwhm_range,
                SAMPLE_VARIABILITY=sample_variability,
                SAMPLE_MINSN=snr,
                SAMPLE_MAXELLIP=sample_max_ellp
            )
            with open(path_psfex, "w", newline="\n") as f: f.write(psfex_content)

            cmd_psfex = ["wsl", "psfex", catalog_name, "-c", "default.psfex"]
            subprocess.run(cmd_psfex, cwd=image_dir, check=True)

            # Parse XML to find NStars_Accepted_Total
            xml_path = os.path.join(image_dir, "psfex.xml")
            accepted_stars = 0
            if os.path.exists(xml_path):
                try:
                    tree = ET.parse(xml_path)
                    root = tree.getroot()

                    # Account for XML namespaces commonly used by VO tools
                    idx = -1
                    fields = root.findall('.//FIELD') + root.findall('.//{*}FIELD')
                    for i, field in enumerate(fields):
                        if field.attrib.get('name') == 'NStars_Accepted_Total':
                            idx = i
                            break

                    if idx >= 0:
                        tds = root.findall('.//TR/TD') + root.findall('.//{*}TR/{*}TD')
                        if len(tds) > idx:
                            accepted_stars = int(tds[idx].text)
                    else:
                        log.warning("Could not locate NStars_Accepted_Total node in XML.")
                except Exception as e:
                    log.warning(f"XML parsing error: {e}")

            log.info(f"  -> Stars Accepted: {accepted_stars}")
            if accepted_stars >= min_accepted_stars:
                log.info(f"  -> Found >= {min_accepted_stars} optimal stars! Halting iteration.")
                break
            elif snr != snrs_to_try[-1]:
                log.info("  -> Insufficient stars. Attempting lower SNR...")

        log.info(f"Success! Catalog and PSF generated in {image_dir}")
        log.info(f"SExtractor DEBLEND_MINCONT = {deblend_mincont}")

    except subprocess.CalledProcessError as e:
        log.error(f"Subprocess Error: {e}")
        # Propagate the failure so master_run.py / callers see a non-zero exit
        # code instead of a false-positive "Phase 1 OK".
        sys.exit(1)
    finally:
        if args.keep_templates:
            log.info("Kept intermediate template files (--keep-templates).")
        else:
            for p in [path_conv, path_nnw, path_param, path_sex, path_psfex]:
                if os.path.exists(p):
                    os.remove(p)

if __name__ == "__main__":
    main()
