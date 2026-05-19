"""
Master orchestrator for the FRB host pipeline.

Runs all three pipeline phases (SExtractor + PSFEx, Photometry + AstroPath,
GALFIT cutouts + fit) end to end against a staged copy of the user-supplied
flux / inverse-variance FITS pair. The full pipeline always executes for
provenance consistency; the user only chooses *which* outputs to expose in the
final deliverable folder.

Usage (PowerShell / cmd):
    python pipeline_scripts/master_run.py ^
        --image  large_cutouts/20240114A_flux.fits ^
        --invvar large_cutouts/20240114A_invvar.fits ^
        --ra 64.39632 --dec 7.93212 ^
        --outputs all

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
    photometry : calibrated_photometry_results.csv, zero_points.json
    astropath  : astropath_association.png, astropath_posteriors.csv
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
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import yaml


PIPELINE_DIR = Path(__file__).resolve().parent
PHASE1 = PIPELINE_DIR / "SExtractor + PSFEx" / "run_psf_pipeline.py"
PHASE2 = PIPELINE_DIR / "photometry + astropath" / "run_photometry_astropath.py"
PHASE3A = PIPELINE_DIR / "galfit_fitting" / "generate_galfit_cutouts.py"
PHASE3B = PIPELINE_DIR / "galfit_fitting" / "run_galfit_fitting.py"
OUTPUT_DIR = PIPELINE_DIR / "Output"

# Default YAML locations for each phase.
PHASE1_DEFAULT_YAML = PIPELINE_DIR / "SExtractor + PSFEx" / "pipeline_config.yaml"
PHASE2_DEFAULT_YAML = PIPELINE_DIR / "photometry + astropath" / "photometry_astropath_config.yaml"
PHASE3_DEFAULT_YAML = PIPELINE_DIR / "galfit_fitting" / "galfit_config.yaml"


# Files exposed for each --outputs keyword.
TOOL_FILES = {
    "catalog":    ["image.cat"],
    "psf":        ["proto_image.fits", "image.psf"],
    "photometry": ["calibrated_photometry_results.csv", "zero_points.json"],
    "astropath":  ["astropath_association.png", "astropath_posteriors.csv"],
    "galfit":     ["fit.log", "out.fits", "galfit_results.png", "qa_cutout_mask.png", "sky_fit_audit.json"],
}

# Files that must never be exposed in 'all' mode.
ALL_DENYLIST = {
    "default.conv", "default.nnw", "default.param", "default.sex", "default.psfex",
    "default_psf.sex", "photomPSF.param", "_run_astrophysics_wsl.py",
    "image.fits", "invvar.fits",
    # Per-run config snapshots written by master_run.py — kept inside the
    # workdir as provenance only; not part of the science deliverables.
    "pipeline_config.yaml", "photometry_astropath_config.yaml", "galfit_config.yaml",
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
    if args.pixel_scale is not None:
        cfg1["sextractor"]["pixel_scale"] = args.pixel_scale
        cfg2["sextractor_psf"]["pixel_scale"] = args.pixel_scale
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


def write_galfit_config(work_dir: Path, args) -> Path | None:
    """Write Phase-3b config after Phase 2 (needs zero_points.json for ZP)."""
    cfg = _load_yaml(PHASE3_DEFAULT_YAML) if PHASE3_DEFAULT_YAML.is_file() else {}

    if args.pixel_scale is not None:
        cfg["plate_scale_x"] = args.pixel_scale
        cfg["plate_scale_y"] = args.pixel_scale

    zp_source = "cli" if args.galfit_zp is not None else None
    if args.galfit_zp is not None:
        cfg["mag_zeropoint"] = float(args.galfit_zp)
    else:
        zp_path = work_dir / "zero_points.json"
        if zp_path.is_file():
            try:
                zp_data = json.loads(zp_path.read_text(encoding="utf-8"))
                zp_val = zp_data.get("zp_aper_40px")
                if zp_val is not None and math.isfinite(float(zp_val)):
                    cfg["mag_zeropoint"] = float(zp_val)
                    zp_source = "zero_points.json"
            except (json.JSONDecodeError, TypeError, ValueError):
                pass

    path = work_dir / "galfit_config.yaml"
    with open(path, "w", encoding="utf-8") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)

    mag_zp = cfg.get("mag_zeropoint")
    if mag_zp is not None:
        print(f"[master] GALFIT mag_zeropoint={mag_zp:.4f} (source={zp_source or 'galfit_config.yaml'})")
    else:
        print("[master] GALFIT mag_zeropoint: not set (Phase 3b will use zero_points.json or 22.5)")
    return path


def run_phase(label: str, cmd, cwd: Path = None) -> int:
    """Invoke a phase as a subprocess and return its return code (-1 on launch error)."""
    print(f"\n========== {label} ==========")
    cmd_str = " ".join(repr(str(c)) for c in cmd)
    cwd_str = f"  (cwd={cwd})" if cwd else ""
    print(f"$ {cmd_str}{cwd_str}")
    try:
        proc = subprocess.run(
            [str(c) for c in cmd],
            cwd=str(cwd) if cwd else None,
            check=False,
        )
        rc = proc.returncode
    except Exception as e:  # launch-side errors (file missing, permissions, ...)
        print(f"[master] {label}: subprocess could not start: {type(e).__name__}: {e}")
        return -1
    if rc != 0:
        print(f"[master] {label}: exit code {rc}")
    else:
        print(f"[master] {label}: OK")
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
        help="One or more of: " + ", ".join(sorted(TOOL_FILES)) + ", all. Default: all.",
    )
    parser.add_argument("--frb-name", default=None,
                        help="Override the auto-derived FRB tag.")
    parser.add_argument("--keep-workdir", action="store_true",
                        help="Retain Output/<frbname>_<tag>/.workdir/ for inspection.")

    # ---- Optional YAML overrides (None = keep YAML default) ----
    sx = parser.add_argument_group("SExtractor / shared overrides (apply to both phases)")
    sx.add_argument("--detect-thresh", type=float, default=None,
                    help="Detection + analysis threshold (sigma).")
    sx.add_argument("--deblend-mincont", type=float, default=None,
                    help="SExtractor DEBLEND_MINCONT (default 0.005; SExtractor/DES norm).")
    sx.add_argument("--pixel-scale", type=float, default=None,
                    help="Plate scale in arcsec/pixel.")
    sx.add_argument("--seeing-fwhm", type=float, default=None,
                    help="Approximate seeing FWHM in arcsec (improves CLASS_STAR + PSF FWHM init).")
    sx.add_argument("--gain", type=float, default=None,
                    help="CCD gain (e-/ADU).")

    ap = parser.add_argument_group("Phase 2 / AstroPath overrides")
    ap.add_argument("--mag-mode", choices=["mag_40px", "mag_psf", "mag_auto"], default=None,
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

    gf = parser.add_argument_group("Phase 3b / GALFIT overrides")
    gf.add_argument("--galfit-zp", type=float, default=None,
                    help="Override GALFIT J) ZP (default: zp_aper_40px from zero_points.json).")

    args = parser.parse_args()

    image_path = args.image.resolve()
    if not image_path.is_file():
        raise SystemExit(f"[master] --image not found: {image_path}")
    invvar_path = args.invvar.resolve() if args.invvar else None
    if invvar_path is not None and not invvar_path.is_file():
        raise SystemExit(f"[master] --invvar not found: {invvar_path}")

    chosen, is_all = resolve_outputs(args.outputs)
    tag = "all" if is_all else "_".join(chosen)
    frb_name = args.frb_name or derive_frb_name(image_path)

    out_dir = OUTPUT_DIR / f"{frb_name}_{tag}"
    work_dir = out_dir / ".workdir"
    out_dir.mkdir(parents=True, exist_ok=True)

    overrides_summary = {
        k: v for k, v in vars(args).items()
        if k in {"detect_thresh", "deblend_mincont", "pixel_scale", "seeing_fwhm", "gain",
                 "mag_mode", "target_snr_min", "err_a_arcsec", "err_b_arcsec",
                 "err_theta_deg", "p_u", "galfit_zp"} and v is not None
    }
    print("[master] -- run configuration --")
    print(f"[master]   FRB name      : {frb_name}")
    print(f"[master]   Outputs (tag) : {tag} (is_all={is_all})")
    print(f"[master]   Output folder : {out_dir}")
    print(f"[master]   Workdir       : {work_dir}")
    print(f"[master]   Image         : {image_path}")
    print(f"[master]   Invvar        : {invvar_path}")
    print(f"[master]   RA, Dec [deg] : {args.ra}, {args.dec}")
    print(f"[master]   Weight map    : {'enabled (invvar provided)' if invvar_path else 'disabled (no invvar)'}")
    if overrides_summary:
        print(f"[master]   YAML overrides: {overrides_summary}")
    else:
        print("[master]   YAML overrides: <none, defaults used>")

    stage_inputs(work_dir, image_path, invvar_path)
    p1_cfg, p2_cfg = stage_configs(work_dir, args, has_invvar=invvar_path is not None)
    print(f"[master] Wrote per-run configs: {p1_cfg.name}, {p2_cfg.name}")

    # ---------- Phase 1: SExtractor + PSFEx (hard dependency) ----------
    rc1 = run_phase("Phase 1 / SExtractor + PSFEx",
                    [sys.executable, PHASE1, work_dir / "image.fits",
                     "--config", p1_cfg])
    if rc1 != 0:
        raise SystemExit("[master] Phase 1 failed; pipeline cannot continue.")
    for f in ("image.cat", "image.psf", "proto_image.fits", "segmentation_map.fits"):
        if not (work_dir / f).exists():
            raise SystemExit(f"[master] Phase 1 did not produce expected file: {f}")

    # ---------- Phase 2: Photometry + AstroPath (best-effort) ----------
    rc2 = run_phase(
        "Phase 2 / Photometry + AstroPath",
        [sys.executable, PHASE2,
         "--image", "image.fits",
         "--ra", args.ra, "--dec", args.dec,
         "--config", p2_cfg.name],  # path is inside cwd=work_dir
        cwd=work_dir,
    )

    # ---------- Phase 3a: GALFIT cutouts (best-effort) ----------
    rc3a = run_phase(
        "Phase 3a / GALFIT cutouts",
        [sys.executable, PHASE3A,
         "--image", work_dir / "image.fits",
         "--segmap", work_dir / "segmentation_map.fits",
         "--catalog", work_dir / "image.cat",
         "--ra", args.ra, "--dec", args.dec,
         "--outdir", work_dir],
    )

    # ---------- Phase 3b: GALFIT fit (only if 3a produced cutouts) ----------
    rc3b = -1
    if rc3a == 0 and (work_dir / "host_cutout.fits").exists():
        write_galfit_config(work_dir, args)
        rc3b = run_phase(
            "Phase 3b / GALFIT fit",
            [sys.executable, PHASE3B, "--dir", work_dir],
        )
    else:
        print("[master] Skipping Phase 3b: Phase 3a outputs not available.")

    # ---------- Collect ----------
    print("\n========== Collecting outputs ==========")
    collected, missing = collect(work_dir, out_dir, chosen, is_all)

    print(f"[master] Phase exit codes: 1={rc1}, 2={rc2}, 3a={rc3a}, 3b={rc3b}")
    print(f"[master] Collected ({len(collected)}): {collected}")
    if missing:
        print(f"[master] Missing  ({len(missing)}): {missing}")
    print(f"[master] Output folder: {out_dir}")

    # ---------- Cleanup workdir ----------
    if not args.keep_workdir:
        shutil.rmtree(work_dir, ignore_errors=True)
        print(f"[master] Removed workdir: {work_dir}")
    else:
        print(f"[master] Retained workdir at {work_dir}")

    if not collected:
        raise SystemExit("[master] No requested outputs were produced (all phases failed?).")


if __name__ == "__main__":
    main()
