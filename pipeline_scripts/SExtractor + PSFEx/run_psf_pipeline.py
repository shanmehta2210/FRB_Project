import os
import sys
import yaml
import subprocess
import shutil
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
FLUX_APER(15)
FLUXERR_APER(15)
MAG_APER(15)
MAGERR_APER(15)
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

PHOT_APERTURES   4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 15.0, 20.0, 25.0, 30.0, {PHOT_APERTURES_PX}
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
    args = parser.parse_args()

    image_path = os.path.abspath(args.image)
    if not os.path.exists(image_path):
        print(f"Error: Image {image_path} not found.")
        sys.exit(1)

    image_dir = os.path.dirname(image_path)
    image_name = os.path.basename(image_path)
    # Catalog name is fixed so PSFEx writes deterministic image.psf and proto_image.fits
    # regardless of the user-supplied input image filename. Phase 2 / Phase 3 contracts
    # rely on those exact filenames.
    catalog_name = "image.cat"

    config_path = os.path.abspath(args.config)
    if not os.path.exists(config_path):
        print(f"Error: Config file {config_path} not found.")
        sys.exit(1)
    print(f"[*] Phase 1 config: {config_path}")
        
    with open(config_path, "r") as f:
        config = yaml.safe_load(f)
        
    # Read Sextractor Configs
    s_conf = config.get("sextractor", {})
    detect_thresh = s_conf.get("detect_thresh", 10)
    analysis_thresh = s_conf.get("analysis_thresh", 10)
    phot_apertures_px = s_conf.get("phot_apertures_px", 40.0)
    gain = s_conf.get("gain", 1.6)
    mag_zeropoint = s_conf.get("mag_zeropoint", 0.0)
    pixel_scale = s_conf.get("pixel_scale", 0.0)
    seeing_fwhm = s_conf.get("seeing_fwhm", 2.0)
    deblend_mincont = float(s_conf.get("deblend_mincont", 0.005))
    
    # Read PSFEx Configs
    p_conf = config.get("psfex", {})
    psf_sampling = p_conf.get("psf_sampling", 1.0)
    psfvar_degrees = p_conf.get("psfvar_degrees", 0)
    sample_fwhm_range = p_conf.get("sample_fwhm_range", "2.0,10")
    sample_variability = p_conf.get("sample_variability", 0.4)
    sample_minsn = p_conf.get("sample_minsn", 30)
    sample_max_ellp = p_conf.get("sample_max_ellp", 0.2)
    
    # Render Templates
    # 1. PARAM
    param_content = TEMPLATE_PARAM
    
    # Check for weight map — honours use_weight_map from YAML first
    use_weight_map = s_conf.get("use_weight_map", True)
    invvar_name = "invvar.fits"
    invvar_path = os.path.join(image_dir, invvar_name)
    if use_weight_map and os.path.exists(invvar_path):
        weight_type = "MAP_WEIGHT"
        weight_image = invvar_name
        print("[*] Weight map: invvar.fits found and enabled.")
    else:
        weight_type = "NONE"
        weight_image = "NONE"
        if use_weight_map and not os.path.exists(invvar_path):
            print("[!] Warning: use_weight_map=true but invvar.fits not found. Running without weight map.")
        else:
            print("[*] Weight map disabled by config (use_weight_map=false).")

    # 2. SEX
    sex_content = TEMPLATE_SEX.format(
        CATALOG_NAME=catalog_name,
        DETECT_THRESH=detect_thresh,
        ANALYSIS_THRESH=analysis_thresh,
        PHOT_APERTURES_PX=phot_apertures_px,
        MAG_ZEROPOINT=mag_zeropoint,
        GAIN=gain,
        PIXEL_SCALE=pixel_scale,
        SEEING_FWHM=seeing_fwhm,
        DEBLEND_MINCONT=deblend_mincont,
        WEIGHT_TYPE=weight_type,
        WEIGHT_IMAGE=weight_image
    )
    
    # 3. PSFEX
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
        print("[+] Running SExtractor...")
        # Note: image_name is used, we run inside image_dir
        cmd_sex = ["wsl", "source-extractor", image_name, "-c", "default.sex"]
        subprocess.run(cmd_sex, cwd=image_dir, check=True)
        
        # Iterative PSFEx Strategy
        snrs_to_try = [30, 20, 10, 5, 3, 2.5, 2]
        import xml.etree.ElementTree as ET
        
        for snr in snrs_to_try:
            print(f"[*] Executing PSFEx with SAMPLE_MINSN = {snr}")
            # 3. Render dynamic PSFEX for this loop iteration
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
                        print("  [!] Could not locate NStars_Accepted_Total node in XML.")
                except Exception as e:
                    print(f"  [!] XML parsing error: {e}")
            
            print(f"  -> Stars Accepted: {accepted_stars}")
            if accepted_stars >= 25:
                print("  -> Found > 25 optimal stars! Halting iteration.")
                break
            elif snr != snrs_to_try[-1]:
                print("  -> Insufficient stars. Attempting lower SNR...")
        
        print(f"[+] Success! Catalog and PSF generated in {image_dir}")
        print(f"[*] SExtractor DEBLEND_MINCONT = {deblend_mincont}")

    except subprocess.CalledProcessError as e:
        print(f"[!] Subprocess Error: {e}")
    finally:
        # Cleanup Disabled for Inspection phase!
        print("[*] Skipped final cleanup so intermediate templates can be manually inspected.")
        # for p in [path_conv, path_nnw, path_param, path_sex, path_psfex]:
        #     if os.path.exists(p):
        #         os.remove(p)

if __name__ == "__main__":
    main()
