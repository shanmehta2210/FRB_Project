"""
Master orchestrator for the FRB host pipeline.

Runs the pipeline phases (SExtractor + PSFEx, Photometry + AstroPath,
optional Statmorph, GALFIT cutouts + fit) end to end against a staged copy of
the user-supplied flux / inverse-variance FITS pair.

Execution modes:
    Full run (default)   — all phases execute sequentially (--outputs all).
    --outputs <tools>    — run ONLY the phases the requested tools depend on
                           and expose only their files. Dependency chain:
                             catalog / psf        -> Phase 1
                             photometry/astropath -> Phases 1 + 2
                             statmorph            -> Phases 1 + 2 + 3a + statmorph
                             galfit               -> Phases 1 + 2 + 3a + 3b
    --rerun-phase <N>    — skip phases before <N>; requires an existing workdir
                           with prior-phase outputs (useful for tweaking GALFIT
                           without re-running SExtractor + photometry).
    --dry-run            — print the commands that would execute, then exit.

Usage (PowerShell / cmd):
    python pipeline_scripts/master_run.py ^
        --image  large_cutouts/20240114A_flux.fits ^
        --invvar large_cutouts/20240114A_invvar.fits ^
        --ra 64.39632 --dec 7.93212 ^
        --outputs all

    # Re-run only Phase 3b on an existing workdir:
    python pipeline_scripts/master_run.py ^
        --image  large_cutouts/20240114A_flux.fits ^
        --invvar large_cutouts/20240114A_invvar.fits ^
        --ra 64.39632 --dec 7.93212 ^
        --rerun-phase 3b --keep-workdir

    # Override the FRB localisation ellipse + AstroPath unseen prior:
    python pipeline_scripts/master_run.py ^
        --image  large_cutouts/20180924B_flux.fits ^
        --invvar large_cutouts/20180924B_invvar.fits ^
        --ra 326.105235868384 --dec -40.9002526146074 ^
        --err-a-arcsec 0.16 --err-b-arcsec 0.16 --err-theta-deg 0.0 ^
        --p-u 0.05 ^
        --outputs astropath photometry

Layout produced:
    pipeline_scripts/Output/<frbname>_<tag>/
        <files for the chosen tools>
        .workdir/                 (only when --keep-workdir is given)

Tool keywords -> exposed files:
    catalog    : image.cat
    psf        : proto_image.fits, image.psf
    photometry : calibrated_photometry_results.csv, zero_points.json (field_depth)
    astropath  : astropath_association.png, astropath_posteriors.csv,
                 sep_vs_shape_r.png, sep_vs_x_max_reff.png
    galfit     : fit.log, out.fits, galfit_results.png, qa_cutout_mask.png
    all        : every artefact in the workdir except staged inputs and
                 SExtractor/PSFEx text templates.

YAML overrides:
    The default per-phase YAML configurations live alongside each phase script.
    master_run.py loads those defaults, applies any CLI overrides, and writes
    per-run YAML files into the workdir before each phase runs (the phase
    scripts themselves are pointed at the workdir copies via --config).
    The repo defaults are not modified.

The auto-derived <frbname> is the input flux filename's stem with a trailing
"_flux" or "_image" stripped (override with --frb-name). Tag is "all" when
'all' is requested, otherwise the sorted, underscore-joined tool keywords.
"""

import argparse
import datetime
import json
import math
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from pipeline_shared import get_logger, resolve_apertures  # noqa: E402

log = get_logger("master")


PIPELINE_DIR = Path(__file__).resolve().parent
PHASE1 = PIPELINE_DIR / "SExtractor + PSFEx" / "run_psf_pipeline.py"
PHASE2 = PIPELINE_DIR / "photometry + astropath" / "run_photometry_astropath.py"
PHASE3A = PIPELINE_DIR / "galfit_fitting" / "generate_galfit_cutouts.py"
PHASE_STATMORPH = PIPELINE_DIR / "galfit_fitting" / "run_statmorph_pipeline.py"
PHASE3B = PIPELINE_DIR / "galfit_fitting" / "run_galfit_fitting.py"
OUTPUT_DIR = PIPELINE_DIR / "Output"

VALID_RERUN_PHASES = {"1", "2", "3a", "statmorph", "3b"}

# Default YAML locations for each phase.
PHASE1_DEFAULT_YAML = PIPELINE_DIR / "SExtractor + PSFEx" / "pipeline_config.yaml"
PHASE2_DEFAULT_YAML = PIPELINE_DIR / "photometry + astropath" / "photometry_astropath_config.yaml"
PHASE3_DEFAULT_YAML = PIPELINE_DIR / "galfit_fitting" / "galfit_config.yaml"


# Files exposed for each --outputs keyword.
TOOL_FILES = {
    "catalog":    ["image.cat"],
    "psf":        ["proto_image.fits", "image.psf"],
    "photometry": ["calibrated_photometry_results.csv", "zero_points.json",
                   "reference_photometry.json"],
    "astropath":  [
        "astropath_association.png",
        "astropath_posteriors.csv",
        "sep_vs_shape_r.png",
        "sep_vs_x_max_reff.png",
    ],
    "galfit":     ["fit.log", "out.fits", "galfit_results.png", "qa_cutout_mask.png",
                   "sky_fit_audit.json", "reference_photometry.json"],
    "statmorph":  ["statmorph_results.json"],
}

# Phases each --outputs keyword depends on. Only the union of the phases
# required by the requested outputs is executed; everything else is skipped.
# Notes on the dependency chain:
#   * every tool needs Phase 1 (catalog / PSF model / segmentation map);
#   * astropath and photometry are both produced by the Phase 2 bridge;
#   * Phase 3a needs Phase 2's image.psf.cat for the SPREAD galaxy cut and
#     (in AstroPath mode) astropath_posteriors.csv for the host pick;
#   * statmorph and 3b both consume Phase 3a cutouts.
TOOL_PHASES = {
    "catalog":    {"1"},
    "psf":        {"1"},
    "photometry": {"1", "2"},
    "astropath":  {"1", "2"},
    "statmorph":  {"1", "2", "3a", "statmorph"},
    "galfit":     {"1", "2", "3a", "3b"},
}

# Condensed mode emits only these files (fast re-inspection without the full
# catalog / segmap / PSF / CSV weight). The consolidated JSON replaces the
# scattered zero_points.json + sky_fit_audit.json.
CONDENSED_FILES = [
    "pipeline_summary.json",
    "galfit_results.png",
    "astropath_association.png",
    "qa_cutout_mask.png",
]

# Files that must never be exposed in 'all' mode.
ALL_DENYLIST = {
    "default.conv", "default.nnw", "default.param", "default.sex", "default.psfex",
    "default_psf.sex", "photomPSF.param", "_run_astrophysics_wsl.py",
    "image.fits", "invvar.fits",
    # Per-run config snapshots written by master_run.py — kept inside the
    # workdir as provenance only; not part of the science deliverables.
    "pipeline_config.yaml", "photometry_astropath_config.yaml", "galfit_config.yaml",
}


def compute_pixel_scale(image_path: Path) -> float | None:
    """Derive the plate scale (arcsec/pixel) from the FITS WCS header.

    Returns the mean of the two pixel scale axes, or None on failure.
    """
    try:
        from astropy.io import fits as afits
        from astropy.wcs import WCS
        from astropy.wcs.utils import proj_plane_pixel_scales
        with afits.open(image_path) as hdul:
            w = WCS(hdul[0].header).celestial
            if w.naxis == 0:
                return None
            scales_deg = proj_plane_pixel_scales(w)
            scale_arcsec = float(np.mean(scales_deg) * 3600.0)
            return round(scale_arcsec, 6)
    except Exception as exc:
        log.warning(f"Could not compute pixel scale from WCS: {exc}")
        return None


def validate_coordinates(ra: float, dec: float) -> None:
    """Validate RA/Dec are within physical bounds."""
    if not math.isfinite(ra) or not (0.0 <= ra < 360.0):
        raise SystemExit(f"[master] Invalid RA={ra}: must be in [0, 360).")
    if not math.isfinite(dec) or not (-90.0 <= dec <= 90.0):
        raise SystemExit(f"[master] Invalid Dec={dec}: must be in [-90, 90].")


def validate_wcs_containment(image_path: Path, ra: float, dec: float) -> None:
    """Warn if RA/Dec falls outside the FITS image WCS footprint."""
    try:
        from astropy.io import fits as afits
        from astropy.wcs import WCS
        with afits.open(image_path) as hdul:
            w = WCS(hdul[0].header).celestial
            if w.naxis == 0:
                return
            ny, nx = hdul[0].data.shape[-2:]
            px, py = w.world_to_pixel_values(ra, dec)
            if not (0 <= px < nx and 0 <= py < ny):
                log.warning(
                    f"Target RA={ra:.6f}, Dec={dec:.6f} maps to pixel "
                    f"({px:.1f}, {py:.1f}) which is OUTSIDE the image "
                    f"footprint ({nx}x{ny}). Pipeline may fail or fit noise."
                )
            else:
                log.info(f"WCS containment OK: pixel ({px:.1f}, {py:.1f}) in {nx}x{ny} image.")
    except Exception as exc:
        log.warning(f"WCS containment check skipped: {exc}")


def _parse_psfex_xml(work_dir: Path) -> dict | None:
    """Extract PSF quality metrics from psfex.xml (Phase 1 output)."""
    xml_path = work_dir / "psfex.xml"
    if not xml_path.is_file():
        return None
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        fields = root.findall(".//{*}FIELD") or root.findall(".//FIELD")
        tds = root.findall(".//{*}TR/{*}TD") or root.findall(".//TR/TD")
        if not fields or not tds:
            return None
        field_map = {}
        for i, f in enumerate(fields):
            name = f.attrib.get("name", "")
            if name and i < len(tds):
                field_map[name] = tds[i].text
        metrics: dict = {}
        for key in ("FWHM_FromFluxRadius_Mean", "FWHM_FromFluxRadius_Min",
                     "FWHM_FromFluxRadius_Max", "FWHM_FromFluxRadius_StDev",
                     "Ellipticity_Mean", "Ellipticity_StDev",
                     "Residuals_Mean", "Residuals_StDev",
                     "NStars_Accepted_Total", "NStars_Rejected_Total",
                     "Chi2_Mean", "Chi2_StDev"):
            val = field_map.get(key)
            if val is not None:
                try:
                    metrics[key] = float(val) if "." in val or "e" in val.lower() else int(val)
                except (ValueError, TypeError):
                    pass
        return metrics if metrics else None
    except Exception:
        return None


def _mc_inclination(
    b_a: float, b_a_err: float, q0: float = 0.2, n_samples: int = 10000
) -> dict | None:
    """Monte Carlo inclination with asymmetric confidence intervals.

    Draws N samples from N(b/a, b/a_err), converts each to inclination via
    cos²i = (q² - q0²)/(1 - q0²), and reports median + 16th/84th percentile
    errors. This correctly captures the asymmetric error distribution near
    q → q0 (edge-on) where analytic symmetric propagation breaks down.

    References: Holmberg (1946), Hubble (1926), Padilla & Strauss (2008).
    """
    if not math.isfinite(b_a) or not math.isfinite(b_a_err) or b_a_err <= 0:
        return None
    rng = np.random.default_rng(42)
    q_samples = rng.normal(b_a, b_a_err, n_samples)
    q_samples = np.clip(q_samples, 1e-6, 1.0)
    val = (q_samples ** 2 - q0 ** 2) / (1.0 - q0 ** 2)
    val = np.clip(val, 0.0, 1.0)
    inc_samples = np.degrees(np.arccos(np.sqrt(val)))
    p16, p50, p84 = np.percentile(inc_samples, [16, 50, 84])
    return {
        "median_deg": round(float(p50), 2),
        "err_lo_deg": round(float(p50 - p16), 2),
        "err_hi_deg": round(float(p84 - p50), 2),
        "p16_deg": round(float(p16), 2),
        "p84_deg": round(float(p84), 2),
    }


def derive_frb_name(image_path: Path) -> str:
    """`20240114A_flux.fits` -> `20240114A`. Falls back to the bare stem."""
    name = image_path.stem
    lower = name.lower()
    for suf in ("_flux", "_image"):
        if lower.endswith(suf):
            return name[: -len(suf)]
    return name


def resolve_outputs(values):
    """Return (sorted_canonical_list, is_all_flag)."""
    if not values:
        return list(TOOL_FILES.keys()), True
    raw = [v.lower() for v in values]
    if "all" in raw:
        return list(TOOL_FILES.keys()), True
    chosen = []
    for v in raw:
        if v not in TOOL_FILES:
            raise SystemExit(
                f"[master] Unknown --outputs entry '{v}'. "
                f"Valid keywords: {sorted(TOOL_FILES)} or 'all'."
            )
        if v not in chosen:
            chosen.append(v)
    return sorted(chosen), False


def stage_inputs(work_dir: Path, image_path: Path, invvar_path):
    """Reset the workdir and copy inputs in under the canonical names."""
    if work_dir.exists():
        shutil.rmtree(work_dir, ignore_errors=True)
    work_dir.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(image_path, work_dir / "image.fits")
    if invvar_path is not None:
        shutil.copyfile(invvar_path, work_dir / "invvar.fits")


def _load_yaml(path: Path) -> dict:
    with open(path, "r") as f:
        return yaml.safe_load(f) or {}


def stage_configs(work_dir: Path, args, has_invvar: bool):
    """Build per-run YAML configs in the workdir from defaults + CLI overrides.

    Returns (phase1_cfg_path, phase2_cfg_path).  galfit_config.yaml is also
    written when --galfit-zp was given (Phase 3b reads it from the workdir).

    Logic:
      * Defaults are sourced from the in-repo YAMLs alongside each phase script.
      * Any CLI override that is not None replaces the matching field.
      * use_weight_map is set explicitly to has_invvar so the phase scripts
        do not need to second-guess (the per-script fallback warning still
        fires if invvar.fits is somehow absent at runtime).
    """
    cfg1 = _load_yaml(PHASE1_DEFAULT_YAML)
    cfg2 = _load_yaml(PHASE2_DEFAULT_YAML)
    cfg1.setdefault("sextractor", {})
    cfg1.setdefault("psfex", {})
    cfg2.setdefault("sextractor_psf", {})
    cfg2.setdefault("astropath", {})

    # ---- common SExtractor knobs (apply to both phases) ----
    if args.detect_thresh is not None:
        cfg1["sextractor"]["detect_thresh"] = args.detect_thresh
        cfg1["sextractor"]["analysis_thresh"] = args.detect_thresh
        cfg2["sextractor_psf"]["detect_thresh"] = args.detect_thresh
        cfg2["sextractor_psf"]["analysis_thresh"] = args.detect_thresh
    if args.seeing_fwhm is not None:
        cfg1["sextractor"]["seeing_fwhm"] = args.seeing_fwhm
        cfg2["sextractor_psf"]["seeing_fwhm"] = args.seeing_fwhm
    if args.gain is not None:
        cfg1["sextractor"]["gain"] = args.gain
        cfg2["sextractor_psf"]["gain"] = args.gain
    if args.deblend_mincont is not None:
        cfg1["sextractor"]["deblend_mincont"] = args.deblend_mincont
        cfg2["sextractor_psf"]["deblend_mincont"] = args.deblend_mincont

    # ---- weight-map: default true, auto-disable if no invvar provided ----
    cfg1["sextractor"]["use_weight_map"] = bool(has_invvar)
    cfg2["sextractor_psf"]["use_weight_map"] = bool(has_invvar)

    # ---- Phase 2 / AstroPath knobs ----
    if args.mag_mode is not None:
        cfg2["sextractor_psf"]["mag_mode"] = args.mag_mode
    if args.target_snr_min is not None:
        cfg2["astropath"]["target_snr_min"] = args.target_snr_min
    if args.err_a_arcsec is not None:
        cfg2["astropath"]["err_a_arcsec"] = args.err_a_arcsec
    if args.err_b_arcsec is not None:
        cfg2["astropath"]["err_b_arcsec"] = args.err_b_arcsec
    if args.err_theta_deg is not None:
        cfg2["astropath"]["err_theta_deg"] = args.err_theta_deg
    if args.p_u is not None:
        cfg2["astropath"]["p_u"] = args.p_u

    p1_path = work_dir / "pipeline_config.yaml"
    p2_path = work_dir / "photometry_astropath_config.yaml"
    with open(p1_path, "w") as f:
        yaml.safe_dump(cfg1, f, sort_keys=False)
    with open(p2_path, "w") as f:
        yaml.safe_dump(cfg2, f, sort_keys=False)

    return p1_path, p2_path


def write_galfit_config(work_dir: Path, args, plate_scale: float | None = None) -> Path | None:
    """Write Phase-3b config after Phase 2 (needs zero_points.json for ZP)."""
    cfg = _load_yaml(PHASE3_DEFAULT_YAML) if PHASE3_DEFAULT_YAML.is_file() else {}

    if plate_scale is not None:
        cfg["plate_scale_x"] = plate_scale
        cfg["plate_scale_y"] = plate_scale

    zp_source = "cli" if args.galfit_zp is not None else None
    if args.galfit_zp is not None:
        cfg["mag_zeropoint"] = float(args.galfit_zp)
    else:
        zp_path = work_dir / "zero_points.json"
        if zp_path.is_file():
            try:
                zp_data = json.loads(zp_path.read_text(encoding="utf-8"))
                zp_val = zp_data.get("zp_aper") or zp_data.get("zp_aper_40px")
                if zp_val is not None and math.isfinite(float(zp_val)):
                    cfg["mag_zeropoint"] = float(zp_val)
                    zp_source = "zero_points.json"
            except (json.JSONDecodeError, TypeError, ValueError):
                pass
        if "mag_zeropoint" not in cfg:
            from pipeline_shared import header_mag_zeropoint_from_fits

            zp_hdr = header_mag_zeropoint_from_fits(work_dir / "image.fits")
            if zp_hdr is not None:
                cfg["mag_zeropoint"] = float(zp_hdr)
                zp_source = "image.fits header"

    path = work_dir / "galfit_config.yaml"
    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)

    mag_zp = cfg.get("mag_zeropoint")
    if mag_zp is not None:
        log.info(f"GALFIT mag_zeropoint={mag_zp:.4f} (source={zp_source or 'galfit_config.yaml'})")
    else:
        log.info("GALFIT mag_zeropoint: not set (Phase 3b will use zero_points.json, FITS header, or 22.5)")
    return path


def run_phase(label: str, cmd, cwd: Path | None = None) -> int:
    """Invoke a phase as a subprocess and return its return code (-1 on launch error)."""
    log.info(f"========== {label} ==========")
    cmd_str = " ".join(repr(str(c)) for c in cmd)
    cwd_str = f"  (cwd={cwd})" if cwd else ""
    log.info(f"$ {cmd_str}{cwd_str}")
    try:
        # stdin is closed for the whole phase subtree so no child process
        # (including WSL tools like GALFIT) can ever hang on console input.
        proc = subprocess.run(
            [str(c) for c in cmd],
            cwd=str(cwd) if cwd else None,
            check=False,
            stdin=subprocess.DEVNULL,
        )
        rc = proc.returncode
    except Exception as e:  # launch-side errors (file missing, permissions, ...)
        log.error(f"{label}: subprocess could not start: {type(e).__name__}: {e}")
        return -1
    if rc != 0:
        log.warning(f"{label}: exit code {rc}")
    else:
        log.info(f"{label}: OK")
    return rc


def collect(work_dir: Path, out_dir: Path, chosen, is_all):
    collected = []
    missing = []
    if is_all:
        for entry in sorted(work_dir.iterdir()):
            if entry.is_file() and entry.name not in ALL_DENYLIST:
                shutil.copy2(entry, out_dir / entry.name)
                collected.append(entry.name)
        return collected, missing

    seen = set()
    for tool in chosen:
        for fname in TOOL_FILES[tool]:
            if fname in seen:
                continue
            seen.add(fname)
            src = work_dir / fname
            if src.exists():
                shutil.copy2(src, out_dir / fname)
                collected.append(fname)
            else:
                missing.append(fname)
    return collected, missing


def _safe_json(path: Path) -> dict | None:
    """Load a JSON file, returning None on any error."""
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        return None


def _parse_galfit_host(work_dir: Path, plate_scale: float | None = None) -> dict | None:
    """Extract host Sérsic parameters + inclination from fit.log.

    When ``plate_scale`` (arcsec/px) is known, re_arcsec and re_arcsec_err are
    added. MC inclination with asymmetric CIs supplements the analytic error.
    """
    _REPO_ROOT = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(_REPO_ROOT))
    try:
        from scripts.galfit_fitlog_parse import (
            parse_fitlog_file,
            inclination_from_b_a,
            inclination_err_from_b_a_err,
            count_fitted_sersic_components,
        )
    except ImportError:
        return None
    log_path = work_dir / "fit.log"
    data, strategy = parse_fitlog_file(str(log_path))
    if not data:
        return None
    inc = inclination_from_b_a(data.get("b_a"))
    inc_err = inclination_err_from_b_a_err(data.get("b_a"), data.get("b_a_err"))
    n_comp = count_fitted_sersic_components(str(work_dir))

    re_px = data.get("re")
    re_px_err = data.get("re_err")
    re_arcsec = None
    re_arcsec_err = None
    if plate_scale is not None and re_px is not None and math.isfinite(re_px):
        re_arcsec = round(float(re_px) * plate_scale, 4)
        if re_px_err is not None and math.isfinite(re_px_err):
            re_arcsec_err = round(float(re_px_err) * plate_scale, 4)

    mc_inc = None
    b_a = data.get("b_a")
    b_a_err = data.get("b_a_err")
    if b_a is not None and b_a_err is not None:
        mc_inc = _mc_inclination(b_a, b_a_err)

    return {
        "selection_strategy": strategy,
        "n_sersic_components": n_comp,
        "chi2nu": data.get("chi2nu"),
        "x": data.get("x"),
        "y": data.get("y"),
        "mag": data.get("mag"),
        "mag_err": data.get("mag_err"),
        "re_px": re_px,
        "re_px_err": re_px_err,
        "re_arcsec": re_arcsec,
        "re_arcsec_err": re_arcsec_err,
        "n": data.get("n"),
        "n_err": data.get("n_err"),
        "b_a": b_a,
        "b_a_err": b_a_err,
        "pa_deg": data.get("pa"),
        "pa_deg_err": data.get("pa_err"),
        "inclination_deg": inc,
        "inclination_deg_err": inc_err,
        "inclination_mc": mc_inc,
    }


def _best_posterior(work_dir: Path) -> dict | None:
    """Best AstroPath candidate from posteriors CSV."""
    try:
        import pandas as pd
    except ImportError:
        return None
    post_path = work_dir / "astropath_posteriors.csv"
    if not post_path.is_file():
        return None
    try:
        df = pd.read_csv(post_path)
    except Exception:
        return None
    if df.empty or "posterior_O" not in df.columns:
        return None
    best = df.sort_values("posterior_O", ascending=False).iloc[0]
    out: dict = {}
    for key in ("ra_deg", "dec_deg", "posterior_O", "posterior_U",
                "sep_arcsec", "mag", "ang_size", "sex_number"):
        if key in best.index:
            val = best[key]
            try:
                import numpy as np
                if np.isfinite(float(val)):
                    out[key] = float(val) if key != "sex_number" else int(val)
            except (TypeError, ValueError):
                pass
    return out if out else None


def build_pipeline_summary(
    work_dir: Path,
    frb_name: str,
    args,
    rc1: int, rc2: int, rc3a: int, rc3b: int,
    rc_statmorph: int = -1,
    plate_scale: float | None = None,
) -> dict:
    """Consolidate all pipeline JSON outputs + parsed fit.log into one dict."""
    zp = _safe_json(work_dir / "zero_points.json") or {}
    sky = _safe_json(work_dir / "sky_fit_audit.json") or {}
    galfit = _parse_galfit_host(work_dir, plate_scale=plate_scale)
    astropath = _best_posterior(work_dir)
    psf_metrics = _parse_psfex_xml(work_dir)
    statmorph_results = _safe_json(work_dir / "statmorph_results.json")

    summary: dict = {
        "frb": frb_name,
        "ra_deg": args.ra,
        "dec_deg": args.dec,
        "plate_scale_arcsec_px": plate_scale,
        "timestamp_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "phase_exit_codes": {
            "phase1_sextractor_psfex": rc1,
            "phase2_photometry_astropath": rc2,
            "phase3a_galfit_cutouts": rc3a,
            "phase_statmorph": rc_statmorph,
            "phase3b_galfit_fit": rc3b,
        },
        "psf_quality": psf_metrics,
        "photometry": {
            "reference_catalog": zp.get("reference_catalog"),
            "n_calibration_stars": zp.get("n_calibration_stars"),
            "n_ps1_matches": zp.get("n_ps1_matches"),
            "n_legacy_matches": zp.get("n_legacy_matches"),
            "filter_band": zp.get("filter_band"),
            "production_aperture_px": zp.get("production_aperture_px"),
            "match_radius_arcsec": zp.get("match_radius_arcsec"),
            "zp_aper": zp.get("zp_aper") or zp.get("zp_aper_40px"),
            "zp_aper_std": zp.get("zp_aper_std") or zp.get("zp_aper_40px_std"),
            "zp_psf": zp.get("zp_psf"),
            "zp_psf_std": zp.get("zp_psf_std"),
            "zp_auto": zp.get("zp_auto"),
            "zp_auto_std": zp.get("zp_auto_std"),
        },
        "field_depth": zp.get("field_depth"),
        "astropath": astropath,
        "sky_qa": {
            "sky_ref_adu": sky.get("sky_ref_adu"),
            "sky_ref_source": sky.get("sky_ref_source"),
            "sky_tolerance_adu": sky.get("sky_tolerance_adu"),
            "sky_pass1_adu": sky.get("sky_pass1_adu"),
            "sky_pass2_adu": sky.get("sky_pass2_adu"),
            "sky_final_adu": sky.get("sky_final_adu"),
            "retried": sky.get("retried"),
            "passed": sky.get("passed"),
            "failure_reason": sky.get("failure_reason"),
        } if sky else None,
        "statmorph": statmorph_results,
        "galfit_host": galfit,
    }

    # Trusted external r-band magnitude for the host (LS DR10 / PS1). Embeds
    # reference_photometry.json (written earlier in main()) and sets
    # galfit_host.mag_final{,_err,_source}: GALFIT mag when the pipeline ZP is
    # trusted, reference-survey mag when Phase 2 calibration failed.
    from reference_photometry import attach_reference_photometry
    attach_reference_photometry(summary, work_dir)
    return summary


def write_pipeline_summary(work_dir: Path, summary: dict) -> Path:
    """Write consolidated JSON to the workdir."""
    path = work_dir / "pipeline_summary.json"
    with open(path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, default=str)
    log.info(f"Consolidated summary written to {path.name}")
    return path


def main():
    parser = argparse.ArgumentParser(
        description="Master FRB pipeline orchestrator (Phases 1, 2, 3 always run).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--image", required=True, type=Path,
                        help="Path to flux FITS (e.g. large_cutouts/20240114A_flux.fits).")
    parser.add_argument("--invvar", type=Path, default=None,
                        help="Path to inverse-variance FITS (optional).")
    parser.add_argument("--ra", type=float, required=True, help="Host RA [deg].")
    parser.add_argument("--dec", type=float, required=True, help="Host Dec [deg].")
    parser.add_argument(
        "--outputs", nargs="+", default=["all"],
        help="One or more of: " + ", ".join(sorted(TOOL_FILES)) + ", all. Default: all. "
             "Selects both which files are exposed AND which phases execute: only "
             "the phases the requested outputs depend on are run (e.g. "
             "'--outputs astropath' runs Phases 1 + 2 and skips cutouts/statmorph/GALFIT).",
    )
    parser.add_argument("--frb-name", default=None,
                        help="Override the auto-derived FRB tag.")
    parser.add_argument("--keep-workdir", action="store_true",
                        help="Retain Output/<frbname>_<tag>/.workdir/ for inspection.")
    parser.add_argument("--condensed", action="store_true",
                        help="Condensed output: only pipeline_summary.json + the three "
                             "diagnostic PNGs (galfit_results, astropath_association, "
                             "qa_cutout_mask). Useful for quick visual re-inspection "
                             "without the full catalog/FITS weight. Default is the full "
                             "deliverable set.")
    parser.add_argument("--rerun-phase", default=None,
                        choices=list(VALID_RERUN_PHASES),
                        help="Skip phases before this one and re-execute from here. "
                             "Requires --keep-workdir from a prior run so the workdir "
                             "still contains earlier-phase outputs. Valid: "
                             + ", ".join(sorted(VALID_RERUN_PHASES)) + ".")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print the commands that would execute, then exit.")

    # ---- Optional YAML overrides (None = keep YAML default) ----
    sx = parser.add_argument_group("SExtractor / shared overrides (apply to both phases)")
    sx.add_argument("--detect-thresh", type=float, default=None,
                    help="Detection + analysis threshold (sigma).")
    sx.add_argument("--deblend-mincont", type=float, default=None,
                    help="SExtractor DEBLEND_MINCONT (default 0.005; SExtractor/DES norm).")
    sx.add_argument("--seeing-fwhm", type=float, default=None,
                    help="Approximate seeing FWHM in arcsec (improves CLASS_STAR + PSF FWHM init).")
    sx.add_argument("--gain", type=float, default=None,
                    help="CCD gain (e-/ADU).")

    ap = parser.add_argument_group("Phase 2 / AstroPath overrides")
    ap.add_argument("--mag-mode", choices=["mag_aper", "mag_psf", "mag_auto"], default=None,
                    help="Which calibrated magnitude is passed to AstroPath as the candidate mag.")
    ap.add_argument("--target-snr-min", type=float, default=None,
                    help="Minimum SNR_WIN for a source to be considered as an AstroPath candidate.")
    ap.add_argument("--err-a-arcsec", type=float, default=None,
                    help="FRB localization semi-major axis (arcsec).")
    ap.add_argument("--err-b-arcsec", type=float, default=None,
                    help="FRB localization semi-minor axis (arcsec).")
    ap.add_argument("--err-theta-deg", type=float, default=None,
                    help="FRB localization PA in degrees (E of N).")
    ap.add_argument("--p-u", type=float, default=None,
                    help="Prior probability that the true host is unseen.")

    gf = parser.add_argument_group("Phase 3 / GALFIT overrides")
    gf.add_argument("--galfit-zp", type=float, default=None,
                    help="Override GALFIT J) ZP (default: zp_aper from zero_points.json).")
    gf.add_argument(
        "--use-localization-host",
        action="store_true",
        default=None,
        help="Phase 3a: always centre on --ra/--dec (CSV host), ignore AstroPath host pick.",
    )
    gf.add_argument(
        "--use-astropath-host",
        action="store_true",
        default=None,
        help="Phase 3a: use AstroPath posteriors when available (overrides galfit_config cutouts).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Override output root (default: pipeline_scripts/Output). "
             "Writes <output-dir>/<frb>_<tag>/.",
    )

    args = parser.parse_args()

    # ---- Input validation ----
    validate_coordinates(args.ra, args.dec)

    image_path = args.image.resolve()
    if not image_path.is_file():
        raise SystemExit(f"[master] --image not found: {image_path}")
    invvar_path = args.invvar.resolve() if args.invvar else None
    if invvar_path is not None and not invvar_path.is_file():
        raise SystemExit(f"[master] --invvar not found: {invvar_path}")

    chosen, is_all = resolve_outputs(args.outputs)
    tag = "all" if is_all else "_".join(chosen)
    frb_name = args.frb_name or derive_frb_name(image_path)

    # Phases required by the requested outputs. Execution (not just file
    # collection) is selective: e.g. `--outputs astropath` runs Phase 1
    # (SExtractor + PSFEx) and Phase 2 (PSF-aware SExtractor + calibration +
    # AstroPath) because AstroPath needs them, but skips 3a/statmorph/3b.
    required_phases: set[str] = set()
    for tool in chosen:
        required_phases |= TOOL_PHASES[tool]

    out_root = args.output_dir.resolve() if args.output_dir else OUTPUT_DIR
    out_dir = out_root / f"{frb_name}_{tag}"
    work_dir = out_dir / ".workdir"
    # NOTE: out_dir is created only after the --dry-run early exit below, so
    # dry runs never leave an Output/<FRB>_<tag>/ shell behind (stray folders
    # pollute pipeline_galfit_results.csv).

    rerun = args.rerun_phase
    rerun_order = ["1", "2", "3a", "statmorph", "3b"]

    overrides_summary = {
        k: v for k, v in vars(args).items()
        if k in {"detect_thresh", "deblend_mincont", "seeing_fwhm", "gain",
                 "mag_mode", "target_snr_min", "err_a_arcsec", "err_b_arcsec",
                 "err_theta_deg", "p_u", "galfit_zp",
                 "use_localization_host", "use_astropath_host"} and v is not None
    }
    log.info("-- run configuration --")
    log.info(f"  FRB name      : {frb_name}")
    log.info(f"  Outputs (tag) : {tag} (is_all={is_all})")
    log.info(f"  Phases to run : {[p for p in rerun_order if p in required_phases]}")
    log.info(f"  Output folder : {out_dir}")
    log.info(f"  Workdir       : {work_dir}")
    log.info(f"  Image         : {image_path}")
    log.info(f"  Invvar        : {invvar_path}")
    log.info(f"  RA, Dec [deg] : {args.ra}, {args.dec}")
    log.info(f"  Weight map    : {'enabled (invvar provided)' if invvar_path else 'disabled (no invvar)'}")
    if rerun:
        log.info(f"  --rerun-phase : {rerun}")
    if args.dry_run:
        log.info("  --dry-run     : enabled")
    if overrides_summary:
        log.info(f"  YAML overrides: {overrides_summary}")
    else:
        log.info("  YAML overrides: <none, defaults used>")

    # WCS containment check (non-fatal warning)
    validate_wcs_containment(image_path, args.ra, args.dec)

    def _required(phase_id: str) -> bool:
        """True if the requested --outputs need this phase at all."""
        return phase_id in required_phases

    def _should_run(phase_id: str) -> bool:
        """True if this phase should execute given --outputs and --rerun-phase."""
        if not _required(phase_id):
            return False
        if rerun is None:
            return True
        return rerun_order.index(phase_id) >= rerun_order.index(rerun)

    def _skip_reason(phase_id: str) -> str:
        """Log label for a phase that will not execute."""
        return "not required by --outputs" if not _required(phase_id) else "--rerun-phase"

    # ---- Dry-run mode: print commands and exit BEFORE any staging so no
    # Output/<FRB>_<tag>/ shell is ever created by a dry run ----
    if args.dry_run:
        log.info("DRY RUN — commands that would execute:")
        for phase_id, label, script in (
            ("1", "Phase 1", PHASE1),
            ("2", "Phase 2", PHASE2),
            ("3a", "Phase 3a", PHASE3A),
            ("statmorph", "Statmorph", PHASE_STATMORPH),
            ("3b", "Phase 3b", PHASE3B),
        ):
            if _should_run(phase_id):
                log.info(f"  {label}: {sys.executable} {script} ...")
            else:
                log.info(f"  {label}: SKIPPED ({_skip_reason(phase_id)})")
        log.info("Exiting (--dry-run).")
        return

    out_dir.mkdir(parents=True, exist_ok=True)

    # Stage inputs (skip for rerun to preserve existing workdir)
    if rerun is None:
        stage_inputs(work_dir, image_path, invvar_path)
    else:
        if not work_dir.exists():
            raise SystemExit(
                f"[master] --rerun-phase={rerun} requires an existing workdir at "
                f"{work_dir}. Run the full pipeline first with --keep-workdir."
            )
        log.info(f"Re-run mode: reusing existing workdir {work_dir}")

    p1_cfg, p2_cfg = stage_configs(work_dir, args, has_invvar=invvar_path is not None)
    log.info(f"Wrote per-run configs: {p1_cfg.name}, {p2_cfg.name}")

    _, prod_aper_index, prod_aper_diam = resolve_apertures(_load_yaml(p1_cfg).get("sextractor", {}))
    log.info(f"Production aperture: {prod_aper_diam:g} px (MAG_APER index {prod_aper_index})")

    # Derive plate scale from the FITS WCS header (never hardcoded).
    plate_scale = compute_pixel_scale(image_path)
    if plate_scale is not None:
        log.info(f"Pixel scale from WCS: {plate_scale:.6f} arcsec/px")
        # Inject into workdir YAMLs so phase scripts have the computed value.
        for cfg_path in (p1_cfg, p2_cfg):
            cfg_tmp = _load_yaml(cfg_path)
            top_key = "sextractor" if "sextractor" in cfg_tmp else "sextractor_psf"
            if top_key in cfg_tmp:
                cfg_tmp[top_key]["pixel_scale"] = plate_scale
            with open(cfg_path, "w") as f_cfg:
                yaml.safe_dump(cfg_tmp, f_cfg, sort_keys=False)
    else:
        log.warning("Could not compute pixel scale from WCS. SExtractor will use its internal default.")
    galfit_defaults = _load_yaml(PHASE3_DEFAULT_YAML) if PHASE3_DEFAULT_YAML.is_file() else {}

    # ---------- Phase 1: SExtractor + PSFEx (hard dependency) ----------
    # Phase 1 is required by every --outputs keyword, so it is only ever
    # skipped via --rerun-phase (in which case the workdir must already
    # contain its outputs — verified below).
    rc1 = 0
    if _should_run("1"):
        rc1 = run_phase("Phase 1 / SExtractor + PSFEx",
                        [sys.executable, PHASE1, work_dir / "image.fits",
                         "--config", p1_cfg])
        if rc1 != 0:
            raise SystemExit("[master] Phase 1 failed; pipeline cannot continue.")
    else:
        log.info(f"Phase 1: SKIPPED ({_skip_reason('1')})")
    for f in ("image.cat", "image.psf", "proto_image.fits", "segmentation_map.fits"):
        if not (work_dir / f).exists():
            raise SystemExit(f"[master] Phase 1 did not produce expected file: {f}")

    # ---------- Phase 2: Photometry + AstroPath ----------
    rc2 = 0
    if _should_run("2"):
        rc2 = run_phase(
            "Phase 2 / Photometry + AstroPath",
            [sys.executable, PHASE2,
             "--image", "image.fits",
             "--ra", args.ra, "--dec", args.dec,
             "--config", p2_cfg.name],
            cwd=work_dir,
        )
    else:
        rc2 = -1 if not _required("2") else 0
        log.info(f"Phase 2: SKIPPED ({_skip_reason('2')})")

    # ---------- Phase 3a: GALFIT cutouts (best-effort) ----------
    cutouts_cfg = galfit_defaults.get("cutouts") or {}
    use_loc_host = bool(cutouts_cfg.get("use_localization_host", False))
    min_post = float(cutouts_cfg.get("min_astropath_posterior", 0.05))
    if args.use_localization_host is True:
        use_loc_host = True
    elif args.use_astropath_host is True:
        use_loc_host = False

    phase3a_cmd = [
        sys.executable, PHASE3A,
        "--image", work_dir / "image.fits",
        "--segmap", work_dir / "segmentation_map.fits",
        "--catalog", work_dir / "image.cat",
        "--ra", str(args.ra), "--dec", str(args.dec),
        "--outdir", str(work_dir),
        "--mag-aper-index", str(int(prod_aper_index)),
    ]
    for cli_flag, cfg_key in (
        ("--no-data-sigma", "no_data_sigma"),
        ("--sigma-rescale-min", "sigma_rescale_min"),
        ("--sigma-rescale-max", "sigma_rescale_max"),
        ("--psf-match-arcsec", "psf_match_arcsec"),
        ("--host-pad", "host_pad"),
        ("--re-sep-factor", "re_sep_factor"),
        ("--neighbor-class-star-max", "neighbor_class_star_max"),
        ("--max-roi-iterations", "max_roi_iterations"),
        ("--max-fit-components", "max_fit_components"),
        ("--max-cutout-side", "max_cutout_side"),
        ("--host-only-min-bbox-side", "host_only_min_bbox_side"),
        ("--host-only-min-elongation", "host_only_min_elongation"),
    ):
        if cfg_key in cutouts_cfg and cutouts_cfg[cfg_key] is not None:
            phase3a_cmd += [cli_flag, str(cutouts_cfg[cfg_key])]
    if use_loc_host:
        phase3a_cmd.append("--no-astropath-override")
    else:
        phase3a_cmd += ["--min-astropath-posterior", str(min_post)]

    rc3a = 0
    if _should_run("3a"):
        if use_loc_host:
            log.info("Phase 3a: use localization host (--ra/--dec); AstroPath posteriors ignored for centre")
        rc3a = run_phase("Phase 3a / GALFIT cutouts", phase3a_cmd)
    else:
        rc3a = -1 if not _required("3a") else 0
        log.info(f"Phase 3a: SKIPPED ({_skip_reason('3a')})")

    # ---------- Phase Statmorph: non-parametric morphology (best-effort) ----------
    # Runs AFTER cutouts are generated but BEFORE GALFIT fitting, so its
    # metrics (CAS, Gini-M20) can motivate future GALFIT configuration.
    rc_statmorph = -1
    if _should_run("statmorph") and (work_dir / "host_cutout.fits").exists():
        if PHASE_STATMORPH.is_file():
            statmorph_cmd = [
                sys.executable, str(PHASE_STATMORPH),
                "--cutout", str(work_dir / "host_cutout.fits"),
                "--sigma", str(work_dir / "host_sigma.fits"),
                "--outdir", str(work_dir),
            ]
            psf_path = work_dir / "proto_image.fits"
            if psf_path.is_file():
                statmorph_cmd += ["--psf", str(psf_path)]
            rc_statmorph = run_phase("Phase Statmorph / non-parametric morphology", statmorph_cmd)
        else:
            log.info("Statmorph: script not found, skipping.")
    elif not _should_run("statmorph"):
        log.info(f"Statmorph: SKIPPED ({_skip_reason('statmorph')})")
    else:
        log.warning("Statmorph: skipped (no host_cutout.fits from Phase 3a).")

    # ---------- Phase 3b: GALFIT fit (only if 3a produced cutouts) ----------
    rc3b = -1
    if _should_run("3b"):
        if (rc3a == 0 or not _should_run("3a")) and (work_dir / "host_cutout.fits").exists():
            write_galfit_config(work_dir, args, plate_scale=plate_scale)
            rc3b = run_phase(
                "Phase 3b / GALFIT fit",
                [sys.executable, PHASE3B, "--dir", work_dir],
            )
        else:
            log.warning("Skipping Phase 3b: Phase 3a outputs not available.")
    else:
        log.info(f"Phase 3b: SKIPPED ({_skip_reason('3b')})")

    # ---------- Reference photometry (always queried, best-effort) ----------
    # External r-band magnitude for the host from LS DR10 Tractor (PS1 DR1
    # fallback outside the LS footprint). Queried regardless of Phase 2
    # success so a failed zero-point calibration never loses the host mag:
    # downstream, galfit_host.mag_final substitutes this value whenever the
    # pipeline ZP is missing or unreliable. Network failure is non-fatal.
    try:
        from reference_photometry import fetch_reference_photometry

        ref = fetch_reference_photometry(
            work_dir, fallback_ra=args.ra, fallback_dec=args.dec
        )
        if ref.get("status") == "ok":
            log.info(
                f"Reference photometry: {ref['survey']} r={ref['mag']} "
                f"(err={ref.get('mag_err')}, sep={ref.get('sep_arcsec')}\", "
                f"coords from {ref.get('coord_source')})"
            )
        else:
            log.warning(f"Reference photometry unavailable: {ref.get('status')}")
    except Exception as e:
        log.warning(f"Reference photometry query failed: {type(e).__name__}: {e}")

    # ---------- Consolidated summary ----------
    summary = build_pipeline_summary(
        work_dir, frb_name, args, rc1, rc2, rc3a, rc3b,
        rc_statmorph=rc_statmorph, plate_scale=plate_scale,
    )
    write_pipeline_summary(work_dir, summary)

    # ---------- Collect ----------
    log.info("========== Collecting outputs ==========")
    if args.condensed:
        collected = []
        missing = []
        for fname in CONDENSED_FILES:
            src = work_dir / fname
            if src.exists():
                shutil.copy2(src, out_dir / fname)
                collected.append(fname)
            else:
                missing.append(fname)
    else:
        collected, missing = collect(work_dir, out_dir, chosen, is_all)
        # Always include the consolidated summary in non-condensed mode too
        summary_src = work_dir / "pipeline_summary.json"
        if summary_src.exists() and "pipeline_summary.json" not in collected:
            shutil.copy2(summary_src, out_dir / "pipeline_summary.json")
            collected.append("pipeline_summary.json")

    log.info(
        f"Phase exit codes: 1={rc1}, 2={rc2}, 3a={rc3a}, "
        f"statmorph={rc_statmorph}, 3b={rc3b} (-1 = not run)"
    )
    log.info(f"Collected ({len(collected)}): {collected}")
    if missing:
        log.warning(f"Missing  ({len(missing)}): {missing}")
    log.info(f"Output folder: {out_dir}")

    # ---------- Cleanup workdir ----------
    if not args.keep_workdir:
        shutil.rmtree(work_dir, ignore_errors=True)
        log.info(f"Removed workdir: {work_dir}")
    else:
        log.info(f"Retained workdir at {work_dir}")

    if not collected:
        raise SystemExit("[master] No requested outputs were produced (all phases failed?).")

    # Completeness gate: every explicitly requested deliverable must come from
    # a phase that succeeded. Full runs must finish GALFIT cutouts + fit;
    # partial Phase 1–2 success is not OK. (Statmorph stays best-effort even
    # when requested — it is skipped gracefully when the package is missing.)
    if is_all and (rc3a != 0 or rc3b != 0):
        raise SystemExit(
            f"[master] Incomplete pipeline: phase3a={rc3a}, phase3b={rc3b} "
            "(GALFIT deliverables missing)."
        )
    if not is_all:
        if "galfit" in chosen and (rc3a != 0 or rc3b != 0):
            raise SystemExit(
                f"[master] --outputs galfit requested but phase3a={rc3a}, "
                f"phase3b={rc3b} (GALFIT deliverables missing)."
            )
        if ("astropath" in chosen or "photometry" in chosen) and rc2 != 0:
            raise SystemExit(
                f"[master] --outputs {'/'.join(t for t in ('photometry', 'astropath') if t in chosen)} "
                f"requested but Phase 2 exited with code {rc2}."
            )


if __name__ == "__main__":
    main()
