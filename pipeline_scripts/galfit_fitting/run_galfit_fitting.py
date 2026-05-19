"""
Build galfit.feedme, run GALFIT, and QA the global sky level.

Sky policy:
  - Seed sky from SExtractor BACKGROUND (host row 0 in host_components.csv), ADU.
  - After each fit, compare fit.log sky to that reference.
  - If |delta| > sky_tolerance_adu, rerun with a soft constraint holding sky
    within ±sky_tolerance_adu of the feedme seed (GALFIT constraints.txt syntax).
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import yaml
from astropy.io import fits
from astropy.visualization import ZScaleInterval

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from scripts.galfit_fitlog_parse import (  # noqa: E402
    parse_fitlog_sky_level,
    sky_level_is_plausible,
)

_GALFIT_CRASH_MARKERS = (
    "GALFIT crashed",
    "Singular Matrix",
    "Numerical Recipes run-time error",
)

_GALFIT_ARTIFACTS = (
    "fit.log",
    "out.fits",
    "galfit.01",
    "galfit.02",
    "galfit.03",
    "galfit.04",
)


def load_mag_zeropoint(wkdir: str, config: dict) -> tuple[float, str]:
    """GALFIT J) ZP: galfit_config > zero_points.json zp_aper_40px > 22.5."""
    if "mag_zeropoint" in config:
        val = float(config["mag_zeropoint"])
        if math.isfinite(val):
            return val, "galfit_config"
    zp_path = os.path.join(wkdir, "zero_points.json")
    if os.path.isfile(zp_path):
        try:
            with open(zp_path, encoding="utf-8") as f:
                zp_data = json.load(f)
            zp_val = zp_data.get("zp_aper_40px")
            if zp_val is not None:
                val = float(zp_val)
                if math.isfinite(val):
                    return val, "zero_points.json"
        except (json.JSONDecodeError, TypeError, ValueError):
            pass
    return 22.5, "default_22.5"


def load_sky_ref(comp_df: pd.DataFrame) -> tuple[float, str]:
    """SExtractor BACKGROUND in ADU for the global sky seed."""
    if "BACKGROUND" not in comp_df.columns:
        return 0.0, "missing_column"
    host_bg = comp_df.iloc[0].get("BACKGROUND")
    if pd.notna(host_bg) and math.isfinite(float(host_bg)):
        return float(host_bg), "host_row"
    series = pd.to_numeric(comp_df["BACKGROUND"], errors="coerce")
    med = series.median()
    if pd.notna(med) and math.isfinite(float(med)):
        return float(med), "median_components"
    return 0.0, "fallback_zero"


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
I) {xmax + 24} {ymax + 24}  # Size of convolution box
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
        raw_mag = row.get("MAG_40PX", row["MAG_AUTO"])
        mag = float(raw_mag) + float(mag_zeropoint)
        re = row["FLUX_RADIUS"]
        if pd.isna(re) or re <= 0:
            re = 1.0
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
        constraints += f"{comp_num} re 1.5 to 100.0\n"
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
    for name in _GALFIT_ARTIFACTS:
        path = os.path.join(wkdir, name)
        if os.path.isfile(path):
            os.remove(path)


def ensure_proto_image(wkdir: str) -> bool:
    proto_dst = os.path.join(wkdir, "proto_image.fits")
    if os.path.exists(proto_dst):
        return True
    proto_src = os.path.join(wkdir, "..", "SExtractor + PSFEx", "proto_image.fits")
    if os.path.exists(proto_src):
        shutil.copy(proto_src, proto_dst)
        print("[*] Synced proto_image.fits from ../SExtractor + PSFEx/.")
        return True
    print(f"[!] ABORT: proto_image.fits missing — looked in {wkdir} and {os.path.dirname(proto_src)}.")
    print("[!] Run pipeline_scripts/SExtractor + PSFEx/run_psf_pipeline.py first.")
    return False


def run_galfit(wkdir: str) -> bool:
    print("\n[*] Running wsl galfit...")
    log_path = os.path.join(wkdir, "fit.log")
    try:
        proc = subprocess.run(
            ["wsl", "galfit", "galfit.feedme"],
            cwd=wkdir,
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        combined = (proc.stdout or "") + (proc.stderr or "")
    except (subprocess.CalledProcessError, FileNotFoundError, OSError) as exc:
        print(f"[!] GALFIT execution failed: {exc}")
        return False

    if combined.strip():
        print(combined, end="" if combined.endswith("\n") else "\n")

    if any(marker in combined for marker in _GALFIT_CRASH_MARKERS):
        print("[!] GALFIT reported a crash (bad parameters / singular matrix).")
        return False
    if not os.path.isfile(log_path) or os.path.getsize(log_path) == 0:
        print("[!] GALFIT did not produce a non-empty fit.log.")
        return False
    return True


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
        print(f"[*] Standardized visualization saved to {qa_path}")
    except Exception as exc:
        print(f"[!] Diagnostic charting skipped: {exc}")


def write_sky_audit(wkdir: str, audit: dict) -> None:
    path = os.path.join(wkdir, "sky_fit_audit.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump(audit, f, indent=2)
    print(f"[*] Sky QA audit written to {path}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir", default=".", help="Directory containing cutout inputs")
    args = parser.parse_args()

    wkdir = os.path.abspath(args.dir)
    exit_code = 0

    req_files = ["host_cutout.fits", "host_mask.fits", "host_components.csv"]
    for rf in req_files:
        if not os.path.exists(os.path.join(wkdir, rf)):
            print(f"[!] Error: Missing required file {rf} in {wkdir}")
            return 1

    has_sigma = os.path.exists(os.path.join(wkdir, "host_sigma.fits"))
    if has_sigma:
        print("[*] host_sigma.fits found — using sigma-weighted fit.")
    else:
        print("[!] host_sigma.fits not found — GALFIT will generate its own sigma image.")

    with fits.open(os.path.join(wkdir, "host_cutout.fits")) as hdul:
        ymax, xmax = hdul[0].data.shape

    config: dict = {}
    yaml_file = os.path.join(wkdir, "galfit_config.yaml")
    if os.path.exists(yaml_file):
        with open(yaml_file, "r", encoding="utf-8") as f:
            loaded = yaml.safe_load(f)
            if isinstance(loaded, dict):
                config = loaded

    plate_scale_x = float(config.get("plate_scale_x", 0.262))
    plate_scale_y = float(config.get("plate_scale_y", 0.262))
    mag_zeropoint, mag_zp_source = load_mag_zeropoint(wkdir, config)
    sky_check_enabled = bool(config.get("sky_check_enabled", True))
    sky_tolerance_adu = float(config.get("sky_tolerance_adu", 3.0))
    sky_max_retries = int(config.get("sky_max_retries", 1))
    mag_min = float(config.get("mag_min", 8.0))
    mag_max = float(config.get("mag_max", 40.0))

    print(f"[*] GALFIT photometric zeropoint: {mag_zeropoint:.4f} (source={mag_zp_source})")
    print(f"[*] Sérsic mag constraints: {mag_min:.1f} to {mag_max:.1f}")

    comp_df = pd.read_csv(os.path.join(wkdir, "host_components.csv"))
    sky_ref, sky_ref_source = load_sky_ref(comp_df)
    n_sersic = len(comp_df)
    sky_comp_num = n_sersic + 1

    print(
        f"[*] SExtractor sky reference: {sky_ref:.6g} ADU "
        f"(source={sky_ref_source}, component={sky_comp_num})"
    )
    if sky_ref_source == "fallback_zero":
        print("[!] Warning: no valid BACKGROUND in host_components.csv; sky seed is 0.")

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
        )
        write_configs(wkdir, feedme, constraints)
        mode = "constrained" if constrain_sky else "free"
        print(f"[*] Wrote galfit.feedme and constraints.txt (sky {mode}, seed={sky_ref:.6g} ADU)")

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
        }
        write_sky_audit(wkdir, audit)

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
            print(f"[*] Pass 1 fit.log sky: {sky_pass1:.6g} ADU (|delta|={delta1:.6g} vs ref)")
        else:
            print(
                "[!] Pass 1: no plausible sky in fit.log "
                "(GALFIT may have crashed before writing a summary block)."
            )
            failure_reason = "sky_parse_failed_pass1"
    else:
        failure_reason = "galfit_crash_pass1"

    need_retry = (
        galfit_pass1_ok
        and sky_pass1 is not None
        and sky_check_enabled
        and sky_max_retries > 0
        and not sky_within_tolerance(sky_pass1, sky_ref, sky_tolerance_adu)
    )

    if need_retry:
        retried = True
        print(
            f"[*] Sky QA: |delta|={abs(sky_pass1 - sky_ref):.6g} ADU > {sky_tolerance_adu} — "
            f"retrying with +/-{sky_tolerance_adu} ADU constraint on component {sky_comp_num}"
        )
        clear_galfit_artifacts(wkdir)
        _build(constrain_sky=True)
        galfit_pass2_ok = run_galfit(wkdir)
        if galfit_pass2_ok:
            sky_pass2 = parse_fitlog_sky_level(log_path, sky_ref=sky_ref)
            if sky_pass2 is not None:
                delta2 = abs(sky_pass2 - sky_ref)
                print(f"[*] Pass 2 fit.log sky: {sky_pass2:.6g} ADU (|delta|={delta2:.6g} vs ref)")
            else:
                print("[!] Pass 2: no plausible sky in fit.log.")
                failure_reason = "sky_parse_failed_pass2"
        else:
            failure_reason = "galfit_crash_pass2"
    elif not galfit_pass1_ok:
        exit_code = 1
    elif sky_check_enabled and sky_pass1 is None:
        exit_code = 1

    sky_final = sky_pass2 if retried and sky_pass2 is not None else sky_pass1
    passed = (
        galfit_pass1_ok
        and (
            not sky_check_enabled
            or (
                sky_final is not None
                and sky_within_tolerance(sky_final, sky_ref, sky_tolerance_adu)
            )
        )
    )

    if sky_check_enabled and galfit_pass1_ok and not passed:
        if sky_final is not None and sky_level_is_plausible(sky_final, sky_ref):
            print(
                f"[!] SKY QA FAILED: ref={sky_ref:.6g} ADU, "
                f"final={sky_final:.6g}, tolerance=±{sky_tolerance_adu} ADU"
            )
            failure_reason = failure_reason or "sky_out_of_tolerance"
        else:
            print(
                f"[!] SKY QA FAILED: ref={sky_ref:.6g} ADU, "
                f"no usable fitted sky (reason={failure_reason or 'unknown'})"
            )
        exit_code = 1
    elif sky_check_enabled and passed:
        print(f"[*] Sky QA passed (final sky={sky_final:.6g} ADU)")

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
