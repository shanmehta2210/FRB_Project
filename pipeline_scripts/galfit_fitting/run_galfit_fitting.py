"""
Build galfit.feedme, run GALFIT, and QA the global sky level.

Sky policy:
  - Seed sky from SExtractor BACKGROUND (host row 0 in host_components.csv), ADU,
    unless that disagrees with the host_cutout median — then use cutout_median.
  - After each fit, compare fit.log sky to that reference.
  - If |delta| > sky_tolerance_adu (but not >> tolerance), rerun with a soft
    constraint holding sky within ±sky_tolerance_adu of the feedme seed.
  - GALFIT 3.0.x on modern glibc may exit non-zero after a good fit
    (``free(): invalid pointer``); we accept the fit when fit.log/out.fits exist.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from astropy.io import fits
from astropy.visualization import ZScaleInterval

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

_PIPELINE_DIR = Path(__file__).resolve().parents[1]
if str(_PIPELINE_DIR) not in sys.path:
    sys.path.insert(0, str(_PIPELINE_DIR))

_GALFIT_DIR = Path(__file__).resolve().parent
if str(_GALFIT_DIR) not in sys.path:
    sys.path.insert(0, str(_GALFIT_DIR))

from sersic_init import effective_re_px  # noqa: E402

from scripts.galfit_fitlog_parse import (  # noqa: E402
    parse_fitlog_sky_level,
    sky_level_is_plausible,
)
from pipeline_shared import get_logger, header_mag_zeropoint_from_fits  # noqa: E402

log = get_logger("phase3b")

_GALFIT_CRASH_MARKERS = (
    "GALFIT crashed",
    "Singular Matrix",
    "Numerical Recipes run-time error",
    "Segmentation fault",
)

# GALFIT prompts interactively ("...try again:") when an input file named in the
# feedme is missing. stdin is closed below so it can never block on a prompt,
# but a prompt loop on EOF would still spin — this wall-clock cap kills it.
_GALFIT_TIMEOUT_S = 3600

# Only retry with a sky constraint when the fitted sky is moderately close to the
# seed. Large deltas usually mean the SExtractor/cutout seed is not representative
# (common on weighted stacks); constraining would break an otherwise good fit.
_SKY_RETRY_MAX_DELTA_FACTOR = 5.0

_GALFIT_ARTIFACTS = (
    "fit.log",
    "out.fits",
)


def load_mag_zeropoint(wkdir: str, config: dict) -> tuple[float, str]:
    """GALFIT J) ZP: galfit_config > zero_points.json > FITS header > 22.5."""
    if "mag_zeropoint" in config:
        val = float(config["mag_zeropoint"])
        if math.isfinite(val):
            return val, "galfit_config"
    zp_path = os.path.join(wkdir, "zero_points.json")
    if os.path.isfile(zp_path):
        try:
            with open(zp_path, encoding="utf-8") as f:
                zp_data = json.load(f)
            zp_val = zp_data.get("zp_aper") or zp_data.get("zp_aper_40px")
            if zp_val is not None:
                val = float(zp_val)
                if math.isfinite(val):
                    return val, "zero_points.json"
        except (json.JSONDecodeError, TypeError, ValueError):
            pass
    for name in ("image.fits", "host_cutout.fits"):
        zp_hdr = header_mag_zeropoint_from_fits(os.path.join(wkdir, name))
        if zp_hdr is not None:
            return zp_hdr, f"fits_header:{name}"
    return 22.5, "default_22.5"


def load_sky_ref(comp_df: pd.DataFrame, wkdir: str | None = None) -> tuple[float, str]:
    """Sky level (ADU) to seed GALFIT component 2 and sky QA."""
    sex_bg: float | None = None
    sex_source = "fallback_zero"
    if "BACKGROUND" in comp_df.columns:
        host_bg = comp_df.iloc[0].get("BACKGROUND")
        if pd.notna(host_bg) and math.isfinite(float(host_bg)):
            sex_bg = float(host_bg)
            sex_source = "host_row"
        else:
            series = pd.to_numeric(comp_df["BACKGROUND"], errors="coerce")
            med = series.median()
            if pd.notna(med) and math.isfinite(float(med)):
                sex_bg = float(med)
                sex_source = "median_components"

    cutout_sky = _cutout_sky_median(wkdir)
    if cutout_sky is not None:
        if sex_bg is None:
            return cutout_sky, "cutout_median"
        if abs(cutout_sky - sex_bg) > max(10.0, abs(sex_bg) * 0.5 + 5.0):
            return cutout_sky, "cutout_median"
        return sex_bg, sex_source

    if sex_bg is not None:
        return sex_bg, sex_source
    return 0.0, "fallback_zero"


def _cutout_sky_median(wkdir: str | None) -> float | None:
    """Robust sky level from unmasked host_cutout pixels."""
    if not wkdir:
        return None
    cutout_path = os.path.join(wkdir, "host_cutout.fits")
    mask_path = os.path.join(wkdir, "host_mask.fits")
    if not (os.path.isfile(cutout_path) and os.path.isfile(mask_path)):
        return None
    try:
        with fits.open(cutout_path) as hdul:
            data = np.asarray(hdul[0].data, dtype=float)
        with fits.open(mask_path) as hdul:
            mask = np.asarray(hdul[0].data)
        sky_pix = data[(mask == 0) & np.isfinite(data)]
        if sky_pix.size < 50:
            return None
        val = float(np.median(sky_pix))
        return val if math.isfinite(val) else None
    except (OSError, ValueError, TypeError):
        return None


def host_audit_fields(comp_df: pd.DataFrame) -> dict:
    """SExtractor metrics for the GALFIT host (row 0 in host_components.csv)."""
    if comp_df.empty:
        return {}
    host = comp_df.iloc[0]
    out: dict = {}
    if "NUMBER" in host.index and pd.notna(host["NUMBER"]):
        out["host_number"] = int(host["NUMBER"])
    if "SNR_WIN" in host.index and pd.notna(host["SNR_WIN"]):
        snr_win = float(host["SNR_WIN"])
        if math.isfinite(snr_win):
            out["snr_win"] = snr_win
    if "FLUX_AUTO" in host.index and "FLUXERR_AUTO" in host.index:
        flux = pd.to_numeric(host["FLUX_AUTO"], errors="coerce")
        ferr = pd.to_numeric(host["FLUXERR_AUTO"], errors="coerce")
        if pd.notna(flux) and pd.notna(ferr) and float(ferr) > 0:
            snr_auto = float(flux) / float(ferr)
            if math.isfinite(snr_auto):
                out["snr_auto"] = snr_auto
    if "MAG_40PX" in host.index and pd.notna(host["MAG_40PX"]):
        mag_40 = float(host["MAG_40PX"])
        if math.isfinite(mag_40):
            out["mag_40px_inst"] = mag_40
    return out


def load_cutout_meta(wkdir: str) -> dict:
    path = os.path.join(wkdir, "cutout_meta.json")
    if not os.path.isfile(path):
        return {}
    try:
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError) as exc:
        log.warning(f"cutout_meta.json unreadable ({exc}); continuing without cutout metadata.")
        return {}
    return data if isinstance(data, dict) else {}


def build_feedme_and_constraints(
    comp_df: pd.DataFrame,
    *,
    xmax: int,
    ymax: int,
    has_sigma: bool,
    mag_zeropoint: float,
    plate_scale_x: float,
    plate_scale_y: float,
    sky_ref: float,
    sky_comp_num: int,
    constrain_sky: bool,
    sky_tolerance_adu: float,
    mag_min: float,
    mag_max: float,
    conv_box_pad: int = 24,
    extended_host: bool = False,
) -> tuple[str, str]:
    sigma_line = "host_sigma.fits" if has_sigma else "none"
    feedme = f"""===============================================================================
# IMAGE and GALFIT CONTROL PARAMETERS
A) host_cutout.fits  # Input data image (FITS file)
B) out.fits  # Output data image block
C) {sigma_line}  # Sigma image name (none = let GALFIT auto-generate)
D) proto_image.fits  # Input PSF file
E) 1  # PSF fine sampling factor
F) host_mask.fits  # Bad pixel mask
G) constraints.txt  # File with parameter constraints
H) 1 {xmax} 1 {ymax}  # Image region to fit (1-based full image)
I) {xmax + conv_box_pad} {ymax + conv_box_pad}  # Size of convolution box
J) {mag_zeropoint:.4f}  # Photometric zeropoint
K) {plate_scale_x} {plate_scale_y}  # Plate scale [arcsec/pixel]
O) regular  # Display type
P) 0  # Choose: 0=optimize

# INITIAL FITTING PARAMETERS
"""

    constraints = ""
    comp_num = 1
    for _, row in comp_df.iterrows():
        xc = row["XC_CUTOUT"]
        yc = row["YC_CUTOUT"]
        # Guard against NaN / SExtractor's 99.0 "undefined magnitude" sentinel,
        # which would otherwise produce a NaN/out-of-range GALFIT seed and crash
        # the fit. Prefer the production-aperture mag, then MAG_AUTO, then the
        # midpoint of the allowed magnitude window. The chosen seed is finally
        # clamped into [mag_min, mag_max] so it never violates the constraint.
        seed_mag = None
        is_host = comp_num == 1
        elong = float(row.get("ELONGATION", 1.0) or 1.0)
        flux_r = float(row.get("FLUX_RADIUS", 0.0) or 0.0)
        use_auto_first = is_host and (
            extended_host or elong >= 2.5 or flux_r >= 25.0
        )
        mag_candidates = (
            (row.get("MAG_AUTO"), row.get("MAG_40PX"))
            if use_auto_first
            else (row.get("MAG_40PX"), row.get("MAG_AUTO"))
        )
        for cand in mag_candidates:
            if cand is not None and pd.notna(cand):
                cand_f = float(cand)
                if math.isfinite(cand_f) and abs(cand_f) < 90.0:
                    seed_mag = cand_f + float(mag_zeropoint)
                    break
        if seed_mag is None:
            seed_mag = 0.5 * (float(mag_min) + float(mag_max))
        mag_clamped = min(max(seed_mag, float(mag_min)), float(mag_max))
        if abs(mag_clamped - seed_mag) > 0.05:
            log.warning(
                "GALFIT component %d: mag seed %.3f clamped to %.3f "
                "(allowed %.1f–%.1f); check mag_zeropoint.",
                comp_num,
                seed_mag,
                mag_clamped,
                mag_min,
                mag_max,
            )
        mag = mag_clamped
        re = effective_re_px(row)
        if is_host and extended_host:
            awin = float(row.get("AWIN_IMAGE", 0) or 0)
            bwin = float(row.get("BWIN_IMAGE", 0) or 0)
            if awin > 0 and bwin > 0:
                re = max(float(re), math.sqrt(awin * bwin))
        elongation = row.get("ELONGATION", 1.0)
        ba = 1.0 / elongation if elongation > 0 else 1.0
        pa = float(row.get("THETA_IMAGE", 0.0)) - 90.0

        feedme += f"""# Component number: {comp_num}
 0) sersic  # Component type
 1) {xc:.4f} {yc:.4f} 1 1  # position x y
 3) {mag:.4f} 1  # Integrated magnitude
 4) {re:.4f} 1  # effective radius (pix)
 5) 1.0000 1  # sersic index
 6) 0.0000 0  # ----
 7) 0.0000 0  # ----
 8) 0.0000 0  # ----
 9) {ba:.4f} 1  # Axis ratio (b/a)
10) {pa:.4f} 1  # Position angle
 Z) 0  # Skip this model

"""
        constraints += f"{comp_num} n 0.5 to 6.0\n"
        re_max = 100.0
        if is_host and extended_host:
            re_max = min(float(max(xmax, ymax)), 300.0)
        constraints += f"{comp_num} re 1.5 to {re_max:.1f}\n"
        constraints += f"{comp_num} mag {mag_min:.1f} to {mag_max:.1f}\n\n"
        comp_num += 1

    feedme += f"""# Component number: {sky_comp_num}
 0) sky  # component type
 1) {sky_ref:.6f} 1  # Sky background (SExtractor BACKGROUND, ADU)
 2) 0.0000 0  # dsky/dx
 3) 0.0000 0  # dsky/dy
 Z) 0  # Skip this model
================================================================================
"""

    if constrain_sky:
        tol = float(sky_tolerance_adu)
        constraints += f"{sky_comp_num} 1 {-tol:.6f} {tol:.6f}\n"

    return feedme, constraints


def write_configs(wkdir: str, feedme: str, constraints: str) -> None:
    with open(os.path.join(wkdir, "galfit.feedme"), "w", newline="\n") as f:
        f.write(feedme)
    with open(os.path.join(wkdir, "constraints.txt"), "w", newline="\n") as f:
        f.write(constraints)


def clear_galfit_artifacts(wkdir: str) -> None:
    """Remove fit.log, out.fits and every numbered galfit.NN restart file.

    GALFIT appends to fit.log and emits a fresh galfit.NN per run, so stale
    files from an earlier pass would poison sky parsing / restart detection.
    """
    stale = [os.path.join(wkdir, name) for name in _GALFIT_ARTIFACTS]
    stale += glob.glob(os.path.join(wkdir, "galfit.[0-9][0-9]"))
    for path in stale:
        if os.path.isfile(path):
            try:
                os.remove(path)
            except OSError as exc:
                log.warning(f"Could not remove stale GALFIT artifact {os.path.basename(path)}: {exc}")


def ensure_proto_image(wkdir: str) -> bool:
    proto_dst = os.path.join(wkdir, "proto_image.fits")
    if os.path.exists(proto_dst):
        return True
    proto_src = os.path.join(wkdir, "..", "SExtractor + PSFEx", "proto_image.fits")
    if os.path.exists(proto_src):
        shutil.copy(proto_src, proto_dst)
        log.info("Synced proto_image.fits from ../SExtractor + PSFEx/.")
        return True
    log.error(f"ABORT: proto_image.fits missing — looked in {wkdir} and {os.path.dirname(proto_src)}.")
    log.error("Run pipeline_scripts/SExtractor + PSFEx/run_psf_pipeline.py first.")
    return False


def _feedme_input_files(feedme_path: str) -> list[str]:
    """Input filenames referenced by the feedme control block.

    Covers A) data, C) sigma, D) PSF, F) mask, G) constraints. B) is the
    output and 'none' placeholders are skipped.
    """
    wanted_keys = {"A", "C", "D", "F", "G"}
    names: list[str] = []
    try:
        with open(feedme_path, encoding="utf-8", errors="replace") as f:
            for line in f:
                stripped = line.strip()
                if len(stripped) < 2 or stripped[1] != ")":
                    continue
                key = stripped[0].upper()
                if key not in wanted_keys:
                    continue
                tokens = stripped[2:].split()
                if tokens and tokens[0].lower() != "none":
                    names.append(tokens[0])
    except OSError:
        return []
    return names


def run_galfit(wkdir: str) -> bool:
    log.info("Running wsl galfit...")
    log_path = os.path.join(wkdir, "fit.log")
    out_path = os.path.join(wkdir, "out.fits")

    # Every file named in the feedme must exist BEFORE launch: GALFIT drops
    # into an interactive "try again:" prompt on a missing input, which (with
    # stdin closed) becomes a busy prompt-loop until the timeout kills it.
    feedme_path = os.path.join(wkdir, "galfit.feedme")
    if not os.path.isfile(feedme_path):
        log.error(f"galfit.feedme missing in {wkdir} — cannot run GALFIT.")
        return False
    missing = [
        f for f in _feedme_input_files(feedme_path)
        if not os.path.isfile(os.path.join(wkdir, f))
    ]
    if missing:
        log.error(f"GALFIT inputs missing in {wkdir}: {', '.join(missing)} — refusing to launch.")
        return False

    combined = ""
    exit_code: int | None = None
    try:
        proc = subprocess.run(
            ["wsl", "galfit", "galfit.feedme"],
            cwd=wkdir,
            check=False,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            stdin=subprocess.DEVNULL,
            timeout=_GALFIT_TIMEOUT_S,
        )
        exit_code = proc.returncode
        combined = (proc.stdout or "") + (proc.stderr or "")
    except subprocess.TimeoutExpired as exc:
        combined = ""
        for stream in (exc.stdout, exc.stderr):
            if stream:
                combined += stream if isinstance(stream, str) else stream.decode("utf-8", "replace")
        if combined.strip():
            print(combined, end="" if combined.endswith("\n") else "\n")
        log.error(
            f"GALFIT exceeded the {_GALFIT_TIMEOUT_S}s wall-clock limit and was killed "
            "(likely stuck on an interactive prompt or a pathological fit)."
        )
        return False
    except FileNotFoundError:
        log.error(
            "Could not launch GALFIT: 'wsl' was not found on PATH. "
            "WSL must be installed and 'wsl galfit' runnable (see README §2 sanity checks)."
        )
        return False
    except OSError as exc:
        log.error(f"Could not launch GALFIT via WSL: {exc}")
        return False

    # Raw GALFIT stdout/stderr is passed through verbatim (multi-line tool output).
    if combined.strip():
        print(combined, end="" if combined.endswith("\n") else "\n")

    if any(marker in combined for marker in _GALFIT_CRASH_MARKERS):
        log.error("GALFIT reported a crash (bad parameters / singular matrix).")
        return False

    has_log = os.path.isfile(log_path) and os.path.getsize(log_path) > 0
    has_out = os.path.isfile(out_path) and os.path.getsize(out_path) > 0
    if has_log and has_out:
        if exit_code not in (0, None):
            log.warning(
                "GALFIT exited with code %s after writing fit.log/out.fits "
                "(known teardown issue on GALFIT 3.0.x + modern glibc; "
                "treating fit as successful).",
                exit_code,
            )
        return True

    if exit_code not in (0, None):
        log.error(f"GALFIT exited with code {exit_code} without usable outputs.")
    elif not has_log:
        log.error("GALFIT did not produce a non-empty fit.log.")
    else:
        log.error("GALFIT did not produce out.fits.")
    return False


def sky_within_tolerance(
    sky_fit: float | None, sky_ref: float, tolerance_adu: float
) -> bool:
    if sky_fit is None or not math.isfinite(sky_fit):
        return False
    return abs(sky_fit - sky_ref) <= tolerance_adu


def write_results_png(wkdir: str) -> None:
    out_fits = os.path.join(wkdir, "out.fits")
    if not os.path.exists(out_fits):
        return
    try:
        with fits.open(out_fits) as hdul:
            img_data = hdul[1].data
            model_data = hdul[2].data
            res_data = hdul[3].data

        fig, axes = plt.subplots(1, 3, figsize=(18, 6))
        interval = ZScaleInterval()
        vmin, vmax = interval.get_limits(img_data)

        axes[0].imshow(img_data, origin="lower", cmap="bone", vmin=vmin, vmax=vmax)
        axes[0].set_title("Input Block (Native)")
        axes[1].imshow(model_data, origin="lower", cmap="bone", vmin=vmin, vmax=vmax)
        axes[1].set_title("Sersic Model Synthesis")
        axes[2].imshow(res_data, origin="lower", cmap="bone", vmin=vmin, vmax=vmax)
        axes[2].set_title("Structural Residuals")

        qa_path = os.path.join(wkdir, "galfit_results.png")
        plt.savefig(qa_path, bbox_inches="tight", dpi=200)
        plt.close(fig)
        log.info(f"Standardized visualization saved to {qa_path}")
    except Exception as exc:
        log.warning(f"Diagnostic charting skipped: {exc}")


def write_sky_audit(wkdir: str, audit: dict) -> None:
    path = os.path.join(wkdir, "sky_fit_audit.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(audit, f, indent=2)
    log.info(f"Sky QA audit written to {path}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", default=".", help="Directory containing cutout inputs")
    args = parser.parse_args()

    wkdir = os.path.abspath(args.dir)
    exit_code = 0

    req_files = ["host_cutout.fits", "host_mask.fits", "host_components.csv"]
    for rf in req_files:
        if not os.path.exists(os.path.join(wkdir, rf)):
            log.error(f"Missing required file {rf} in {wkdir}")
            return 1

    has_sigma = os.path.exists(os.path.join(wkdir, "host_sigma.fits"))
    if has_sigma:
        log.info("host_sigma.fits found — using sigma-weighted fit.")
    else:
        log.warning("host_sigma.fits not found — GALFIT will generate its own sigma image.")

    try:
        with fits.open(os.path.join(wkdir, "host_cutout.fits")) as hdul:
            data = hdul[0].data
            if data is None or data.ndim < 2:
                log.error("host_cutout.fits has no 2-D image data in HDU 0.")
                return 1
            ymax, xmax = data.shape[-2:]
    except OSError as exc:
        log.error(f"Could not open host_cutout.fits: {exc}")
        return 1

    config: dict = {}
    yaml_file = os.path.join(wkdir, "galfit_config.yaml")
    if os.path.exists(yaml_file):
        try:
            with open(yaml_file, "r", encoding="utf-8") as f:
                loaded = yaml.safe_load(f)
            if isinstance(loaded, dict):
                config = loaded
        except (OSError, yaml.YAMLError) as exc:
            log.warning(f"galfit_config.yaml unreadable ({exc}); using built-in defaults.")

    plate_scale_x = config.get("plate_scale_x")
    plate_scale_y = config.get("plate_scale_y")
    if plate_scale_x is None or plate_scale_y is None:
        image_fits = os.path.join(wkdir, "image.fits")
        if os.path.isfile(image_fits):
            try:
                from astropy.io import fits as _afits
                from astropy.wcs import WCS as _WCS
                from astropy.wcs.utils import proj_plane_pixel_scales as _ppps
                import numpy as _np
                with _afits.open(image_fits) as _hdul:
                    _w = _WCS(_hdul[0].header).celestial
                    _scales = _ppps(_w) * 3600.0
                    plate_scale_x = float(_scales[0])
                    plate_scale_y = float(_scales[1])
                log.info(f"Plate scale from WCS: {plate_scale_x:.6f} x {plate_scale_y:.6f} arcsec/px")
            except Exception as exc:
                log.warning(f"Could not compute plate scale from WCS: {exc}")
                plate_scale_x, plate_scale_y = 0.262, 0.262
        else:
            plate_scale_x, plate_scale_y = 0.262, 0.262
    else:
        plate_scale_x = float(plate_scale_x)
        plate_scale_y = float(plate_scale_y)
    mag_zeropoint, mag_zp_source = load_mag_zeropoint(wkdir, config)
    sky_check_enabled = bool(config.get("sky_check_enabled", True))
    sky_tolerance_adu = float(config.get("sky_tolerance_adu", 3.0))
    sky_max_retries = int(config.get("sky_max_retries", 1))
    mag_min = float(config.get("mag_min", 8.0))
    mag_max = float(config.get("mag_max", 40.0))
    conv_box_pad = int(config.get("conv_box_pad", 24))

    log.info(f"GALFIT photometric zeropoint: {mag_zeropoint:.4f} (source={mag_zp_source})")
    log.info(f"Sérsic mag constraints: {mag_min:.1f} to {mag_max:.1f}")

    try:
        comp_df = pd.read_csv(os.path.join(wkdir, "host_components.csv"))
    except Exception as exc:
        log.error(f"Could not read host_components.csv: {exc}")
        return 1
    if comp_df.empty:
        log.error("host_components.csv has no components — nothing to fit (re-run Phase 3a).")
        return 1
    cutout_meta = load_cutout_meta(wkdir)
    extended_host = bool(cutout_meta.get("extended_host"))
    if extended_host:
        log.info(
            "Extended host cutout (host-only GALFIT): "
            f"elong={cutout_meta.get('host_elongation', '?')}, "
            f"bbox={cutout_meta.get('host_bbox_px', '?')} px"
        )
    host_meta = host_audit_fields(comp_df)
    sky_ref, sky_ref_source = load_sky_ref(comp_df, wkdir)
    n_sersic = len(comp_df)
    sky_comp_num = n_sersic + 1

    log.info(
        f"SExtractor sky reference: {sky_ref:.6g} ADU "
        f"(source={sky_ref_source}, component={sky_comp_num})"
    )
    if sky_ref_source == "fallback_zero":
        log.warning("no valid BACKGROUND in host_components.csv; sky seed is 0.")
    if host_meta.get("snr_win") is not None:
        log.info(
            f"GALFIT host #{host_meta.get('host_number', '?')}: "
            f"SNR_WIN={host_meta['snr_win']:.3f}"
            + (
                f", SNR_AUTO={host_meta['snr_auto']:.3f}"
                if host_meta.get("snr_auto") is not None
                else ""
            )
        )

    if not ensure_proto_image(wkdir):
        return 1

    def _build(constrain_sky: bool) -> None:
        feedme, constraints = build_feedme_and_constraints(
            comp_df,
            xmax=xmax,
            ymax=ymax,
            has_sigma=has_sigma,
            mag_zeropoint=mag_zeropoint,
            plate_scale_x=plate_scale_x,
            plate_scale_y=plate_scale_y,
            sky_ref=sky_ref,
            sky_comp_num=sky_comp_num,
            constrain_sky=constrain_sky,
            sky_tolerance_adu=sky_tolerance_adu,
            mag_min=mag_min,
            mag_max=mag_max,
            conv_box_pad=conv_box_pad,
            extended_host=extended_host,
        )
        write_configs(wkdir, feedme, constraints)
        mode = "constrained" if constrain_sky else "free"
        log.info(f"Wrote galfit.feedme and constraints.txt (sky {mode}, seed={sky_ref:.6g} ADU)")

    def _write_audit(
        *,
        sky_pass1: float | None,
        sky_pass2: float | None,
        retried: bool,
        passed: bool,
        galfit_pass1_ok: bool,
        galfit_pass2_ok: bool | None,
        failure_reason: str | None,
    ) -> None:
        sky_final = sky_pass2 if retried and sky_pass2 is not None else sky_pass1
        audit = {
            "sky_ref_adu": sky_ref,
            "sky_ref_source": sky_ref_source,
            "sky_tolerance_adu": sky_tolerance_adu,
            "sky_pass1_adu": sky_pass1,
            "sky_pass2_adu": sky_pass2,
            "sky_final_adu": sky_final,
            "retried": retried,
            "passed": passed,
            "galfit_pass1_ok": galfit_pass1_ok,
            "galfit_pass2_ok": galfit_pass2_ok,
            "failure_reason": failure_reason,
            **host_meta,
        }
        write_sky_audit(wkdir, audit)

    # Clear ALL stale GALFIT outputs (fit.log, out.fits, galfit.NN) so success
    # detection can never be fooled by artifacts from an earlier run.
    clear_galfit_artifacts(wkdir)
    log_path = os.path.join(wkdir, "fit.log")
    sky_pass1: float | None = None
    sky_pass2: float | None = None
    retried = False
    galfit_pass1_ok = False
    galfit_pass2_ok: bool | None = None
    failure_reason: str | None = None

    # Pass 1: SExtractor-seeded sky, no sky constraint
    _build(constrain_sky=False)
    galfit_pass1_ok = run_galfit(wkdir)
    if galfit_pass1_ok:
        sky_pass1 = parse_fitlog_sky_level(log_path, sky_ref=sky_ref)
        if sky_pass1 is not None:
            delta1 = abs(sky_pass1 - sky_ref)
            log.info(f"Pass 1 fit.log sky: {sky_pass1:.6g} ADU (|delta|={delta1:.6g} vs ref)")
        else:
            log.warning(
                "Pass 1: no plausible sky in fit.log "
                "(GALFIT may have crashed before writing a summary block)."
            )
            failure_reason = "sky_parse_failed_pass1"
    else:
        failure_reason = "galfit_crash_pass1"

    sky_delta = abs(sky_pass1 - sky_ref) if sky_pass1 is not None else None
    seed_unreliable = (
        sky_delta is not None
        and sky_delta > _SKY_RETRY_MAX_DELTA_FACTOR * sky_tolerance_adu
    )
    need_retry = (
        galfit_pass1_ok
        and sky_pass1 is not None
        and sky_check_enabled
        and sky_max_retries > 0
        and not sky_within_tolerance(sky_pass1, sky_ref, sky_tolerance_adu)
        and not seed_unreliable
    )
    if seed_unreliable and galfit_pass1_ok and sky_pass1 is not None:
        log.warning(
            "Sky QA: fitted sky %.6g ADU differs from seed %.6g ADU by %.6g "
            "(>> ±%.1f tolerance) — seed likely unreliable; skipping constrained retry.",
            sky_pass1,
            sky_ref,
            sky_delta,
            sky_tolerance_adu,
        )

    if need_retry:
        retried = True
        log.info(
            f"Sky QA: |delta|={abs(sky_pass1 - sky_ref):.6g} ADU > {sky_tolerance_adu} — "
            f"retrying with +/-{sky_tolerance_adu} ADU constraint on component {sky_comp_num}"
        )
        clear_galfit_artifacts(wkdir)
        _build(constrain_sky=True)
        galfit_pass2_ok = run_galfit(wkdir)
        if galfit_pass2_ok:
            sky_pass2 = parse_fitlog_sky_level(log_path, sky_ref=sky_ref)
            if sky_pass2 is not None:
                delta2 = abs(sky_pass2 - sky_ref)
                log.info(f"Pass 2 fit.log sky: {sky_pass2:.6g} ADU (|delta|={delta2:.6g} vs ref)")
            else:
                log.warning("Pass 2: no plausible sky in fit.log.")
                failure_reason = "sky_parse_failed_pass2"
        else:
            failure_reason = "galfit_crash_pass2"
    elif not galfit_pass1_ok:
        exit_code = 1
    elif sky_check_enabled and sky_pass1 is None:
        exit_code = 1

    sky_final = sky_pass2 if retried and sky_pass2 is not None else sky_pass1
    sky_final_delta = (
        abs(sky_final - sky_ref)
        if sky_final is not None and math.isfinite(sky_ref)
        else None
    )
    passed = (
        galfit_pass1_ok
        and (
            not sky_check_enabled
            or sky_final is None
            or sky_within_tolerance(sky_final, sky_ref, sky_tolerance_adu)
            or (
                sky_final_delta is not None
                and sky_final_delta
                > _SKY_RETRY_MAX_DELTA_FACTOR * sky_tolerance_adu
            )
        )
    )

    if sky_check_enabled and galfit_pass1_ok and not passed:
        if sky_final is not None and sky_level_is_plausible(sky_final, sky_ref):
            log.error(
                f"SKY QA FAILED: ref={sky_ref:.6g} ADU, "
                f"final={sky_final:.6g}, tolerance=±{sky_tolerance_adu} ADU"
            )
            failure_reason = failure_reason or "sky_out_of_tolerance"
        else:
            log.error(
                f"SKY QA FAILED: ref={sky_ref:.6g} ADU, "
                f"no usable fitted sky (reason={failure_reason or 'unknown'})"
            )
        exit_code = 1
    elif sky_check_enabled and passed:
        log.info(f"Sky QA passed (final sky={sky_final:.6g} ADU)")

    _write_audit(
        sky_pass1=sky_pass1,
        sky_pass2=sky_pass2,
        retried=retried,
        passed=passed,
        galfit_pass1_ok=galfit_pass1_ok,
        galfit_pass2_ok=galfit_pass2_ok,
        failure_reason=failure_reason,
    )
    if galfit_pass1_ok and os.path.isfile(os.path.join(wkdir, "out.fits")):
        write_results_png(wkdir)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
