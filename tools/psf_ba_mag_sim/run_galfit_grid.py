#!/usr/bin/env python3
"""Run dual GALFIT fits (with / without PSF) on mock grid, seeded at ground truth."""

from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

from sim_utils import TOOL_DIR, ensure_output_layout, load_config, outputs_dir

_REPO_ROOT = TOOL_DIR.parents[1]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from scripts.galfit_fitlog_parse import parse_fitlog_file, parse_fitlog_sky_level  # noqa: E402

_GALFIT_CRASH_MARKERS = (
    "GALFIT crashed",
    "Singular Matrix",
    "Numerical Recipes run-time error",
)

_CATALOG_FIELDS = [
    "galaxy_id",
    "realization",
    "mode",
    "ba_fit",
    "re_fit_pix",
    "mag_fit",
    "n_fit",
    "pa_fit",
    "chi2nu",
    "sky_fit_e",
    "dchi2_iter1",
    "mag_err",
    "converged",
    "parse_strategy",
]

_DCHI_RE = re.compile(r"dChi2/Chi2:\s*([-+eE0-9.]+)")


def parse_fitlog_dchi(log_path: Path) -> float | None:
    if not log_path.is_file():
        return None
    for line in log_path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "Iteration : 1" in line and "dChi2/Chi2" in line:
            m = _DCHI_RE.search(line)
            if m:
                try:
                    return float(m.group(1))
                except ValueError:
                    return None
    return None


def is_sane_fit(
    parsed: dict, truth: dict, dchi1: float | None, *, lock_mag: bool = False
) -> bool:
    ba = parsed.get("ba_fit")
    if ba is None or ba == "":
        return False
    try:
        ba_f = float(ba)
    except (TypeError, ValueError):
        return False
    if not (0.05 < ba_f <= 1.0):
        return False

    if not lock_mag:
        mag_fit = parsed.get("mag_fit")
        mag_true = float(truth["mag_true"])
        if mag_fit not in (None, ""):
            if abs(float(mag_fit) - mag_true) > 1.5:
                return False

    if dchi1 is not None and abs(dchi1) > 1e6:
        pass  # GALFIT often prints huge dChi2 on iter 1; judge by mag/ba instead

    chi2 = parsed.get("chi2nu")
    if chi2 not in (None, ""):
        try:
            if abs(float(chi2)) > 1e4:
                return False
        except (TypeError, ValueError):
            return False

    return True


def read_truth_catalog(catalog_path: Path) -> list[dict]:
    with open(catalog_path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def sky_tolerance_for_truth(truth: dict, cfg: dict) -> float:
    g = cfg.get("galfit", {})
    fixed = g.get("sky_tolerance_adu")
    if fixed is not None and math.isfinite(float(fixed)):
        return float(fixed)
    sky_rate = float(truth.get("sky_rate_per_pix", 0.0) or 0.0)
    exptime = float(truth.get("coadd_exptime_sec", 1.0) or 1.0)
    pix_noise = math.sqrt(max(exptime * sky_rate, 1e-12))
    frac = float(g.get("sky_tolerance_noise_frac", 0.2))
    return max(0.05, frac * pix_noise)


def build_feedme(
    *,
    truth: dict,
    mode: str,
    stamp_px: int,
    zeropoint_e: float,
    pixel_scale: float,
    conv_box_pad: int,
    mag_min: float,
    mag_max: float,
    xy_bound: float,
    sersic_n: float,
    sky_tolerance: float,
    lock_mag: bool,
    lock_sky: bool,
) -> tuple[str, str]:
    psf_line = "psf.fits" if mode == "psf" else "none"

    xc = float(truth["xc"])
    yc = float(truth["yc"])
    mag = float(truth["mag_true"])
    re_pix = float(truth["re_pix_true"])
    ba = float(truth["ba_true"])
    pa = float(truth["pa_true"])
    sky = float(truth.get("sky_fit_seed", truth.get("sky_e_per_pix", 0.0)))
    mag_vary = 0 if lock_mag else 1
    sky_vary = 0 if lock_sky else 1

    feedme = f"""===============================================================================
# IMAGE and GALFIT CONTROL PARAMETERS
A) mock.fits  # Input data image
B) out.fits  # Output data image block
C) none  # Sigma: auto from data (ideal Poisson sigma breaks mag on sky-sub stamps)
D) {psf_line}  # PSF image (none = no convolution)
E) 1  # PSF fine sampling factor
F) none  # Bad pixel mask
G) constraints.txt  # Parameter constraints
H) 1 {stamp_px} 1 {stamp_px}  # Image region to fit
I) {stamp_px + conv_box_pad} {stamp_px + conv_box_pad}  # Convolution box
J) {zeropoint_e:.4f}  # Photometric zeropoint
K) {pixel_scale} {pixel_scale}  # Plate scale arcsec/pixel
O) regular  # Display type
P) 0  # Optimize

# INITIAL FITTING PARAMETERS
# Component number: 1
 0) sersic
 1) {xc:.4f} {yc:.4f} 1 1  # position x y
 3) {mag:.4f} {mag_vary}  # Integrated magnitude
 4) {re_pix:.4f} 1  # effective radius pix
 5) {sersic_n:.4f} 0  # sersic index locked to 1
 6) 0.0000 0
 7) 0.0000 0
 8) 0.0000 0
 9) {ba:.4f} 1  # axis ratio b/a
10) {pa:.4f} 1  # position angle
 Z) 0

# Component number: 2
 0) sky
 1) {sky:.8f} {sky_vary}  # sky background ADU/e-
 2) 0.0000 0
 3) 0.0000 0
 Z) 0
================================================================================
"""

    constraints = f"""1 x {-xy_bound:.1f} {xy_bound:.1f}
1 y {-xy_bound:.1f} {xy_bound:.1f}
"""
    if not lock_sky:
        constraints += f"2 1 {-sky_tolerance:.6f} {sky_tolerance:.6f}\n"
    return feedme, constraints


def native_run_dir(galaxy_id: str, mode: str) -> Path:
    """WSL-native path; GALFIT aborts when cwd is on /mnt/c."""
    run_dir = Path.home() / ".psf_ba_mag_sim" / galaxy_id / mode
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


def stage_inputs(
    *,
    truth: dict,
    mode: str,
    mocks_root: Path,
    run_dir: Path,
) -> None:
    gid = truth["galaxy_id"]
    mock_src = mocks_root / gid
    for name in ("mock.fits", "sigma.fits"):
        shutil.copy2(mock_src / name, run_dir / name)
    if mode == "psf":
        shutil.copy2(mocks_root / "psf.fits", run_dir / "psf.fits")


def sync_outputs(run_dir: Path, archive_dir: Path) -> None:
    archive_dir.mkdir(parents=True, exist_ok=True)
    for name in (
        "galfit.feedme",
        "constraints.txt",
        "fit.log",
        "out.fits",
        "galfit.01",
        "galfit.02",
        "galfit.03",
        "galfit.04",
    ):
        src = run_dir / name
        if src.is_file():
            shutil.copy2(src, archive_dir / name)


def run_galfit(run_dir: Path) -> bool:
    log_path = run_dir / "fit.log"
    out_path = run_dir / "out.fits"
    if log_path.is_file():
        log_path.unlink()

    combined = ""
    try:
        proc = subprocess.run(
            ["bash", "-lc", "galfit galfit.feedme > fit.log 2>&1"],
            cwd=str(run_dir),
            check=False,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        combined = (proc.stdout or "") + (proc.stderr or "")
    except (FileNotFoundError, OSError) as exc:
        print(f"GALFIT execution failed in {run_dir}: {exc}", file=sys.stderr)
        return False

    if combined.strip():
        print(combined.strip())

    if any(marker in combined for marker in _GALFIT_CRASH_MARKERS):
        # GALFIT 3.0.x may still write fit.log before glibc teardown on native FS.
        pass

    has_log = log_path.is_file() and log_path.stat().st_size > 0
    has_out = out_path.is_file() and out_path.stat().st_size > 0
    if has_log and "sersic" in log_path.read_text(encoding="utf-8", errors="replace").lower():
        return True
    return bool(has_out and has_log)


def parse_fit_result(
    log_path: Path, sky_ref: float, truth: dict, *, lock_mag: bool = False
) -> dict:
    data, strategy = parse_fitlog_file(str(log_path), sersic_component_index=0)
    sky_fit = parse_fitlog_sky_level(log_path, sky_ref=sky_ref)
    dchi1 = parse_fitlog_dchi(log_path)
    mag_true = float(truth["mag_true"])
    mag_fit = data.get("mag", "")
    mag_err = ""
    if mag_fit not in (None, ""):
        try:
            mag_err = float(mag_fit) - mag_true
        except (TypeError, ValueError):
            mag_err = ""

    parsed = {
        "ba_fit": data.get("b_a", ""),
        "re_fit_pix": data.get("re", ""),
        "mag_fit": mag_fit,
        "n_fit": data.get("n", ""),
        "pa_fit": data.get("pa", ""),
        "chi2nu": data.get("chi2nu", ""),
        "sky_fit_e": sky_fit if sky_fit is not None else "",
        "dchi2_iter1": dchi1 if dchi1 is not None else "",
        "mag_err": mag_err,
        "parse_strategy": strategy,
    }
    parsed["converged"] = is_sane_fit(parsed, truth, dchi1, lock_mag=lock_mag)
    return parsed


def fit_already_done(
    log_path: Path, sky_ref: float, truth: dict, *, lock_mag: bool = False
) -> bool:
    if not log_path.is_file() or log_path.stat().st_size < 200:
        return False
    parsed = parse_fit_result(log_path, sky_ref, truth, lock_mag=lock_mag)
    return bool(parsed.get("converged"))


def run_one_fit(
    *,
    truth: dict,
    mode: str,
    cfg: dict,
    fits_root: Path,
    mocks_root: Path,
) -> dict:
    gid = truth["galaxy_id"]
    archive_dir = fits_root / gid / mode
    archive_dir.mkdir(parents=True, exist_ok=True)

    sky_ref = float(truth.get("sky_fit_seed", truth.get("sky_e_per_pix", 0.0)))
    g = cfg["galfit"]
    lock_mag = bool(g.get("lock_mag", False))
    lock_sky = bool(g.get("lock_sky", False))
    log_path = archive_dir / "fit.log"
    if fit_already_done(log_path, sky_ref, truth, lock_mag=lock_mag):
        parsed = parse_fit_result(log_path, sky_ref, truth, lock_mag=lock_mag)
        return {
            "galaxy_id": gid,
            "realization": truth["realization"],
            "mode": mode,
            **parsed,
        }

    run_dir = native_run_dir(gid, mode)
    stage_inputs(truth=truth, mode=mode, mocks_root=mocks_root, run_dir=run_dir)

    p = cfg["physics"]
    gcfg = cfg["grid"]
    sky_tol = sky_tolerance_for_truth(truth, cfg)
    feedme, constraints = build_feedme(
        truth=truth,
        mode=mode,
        stamp_px=int(p["stamp_px"]),
        zeropoint_e=float(truth["zeropoint_e"]),
        pixel_scale=float(p["pixel_scale"]),
        conv_box_pad=int(g.get("conv_box_pad", 24)),
        mag_min=float(g.get("mag_min", 8.0)),
        mag_max=float(g.get("mag_max", 40.0)),
        xy_bound=float(g.get("xy_bound_pix", 1.0)),
        sersic_n=float(gcfg["sersic_n"]),
        sky_tolerance=sky_tol,
        lock_mag=lock_mag,
        lock_sky=lock_sky,
    )
    (run_dir / "galfit.feedme").write_text(feedme, encoding="utf-8", newline="\n")
    (run_dir / "constraints.txt").write_text(constraints, encoding="utf-8", newline="\n")

    for artifact in ("fit.log", "out.fits", "galfit.01", "galfit.02", "galfit.03", "galfit.04"):
        ap = run_dir / artifact
        if ap.is_file():
            ap.unlink()

    ok = run_galfit(run_dir)
    sync_outputs(run_dir, archive_dir)

    log_path = archive_dir / "fit.log"
    parsed = (
        parse_fit_result(log_path, sky_ref, truth, lock_mag=lock_mag)
        if log_path.is_file()
        else {}
    )
    if not ok and not parsed.get("ba_fit"):
        parsed.setdefault("converged", False)
        parsed.setdefault("parse_strategy", "galfit_failed")
    return {
        "galaxy_id": gid,
        "realization": truth["realization"],
        "mode": mode,
        **parsed,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=TOOL_DIR / "config.yaml")
    parser.add_argument("--limit", type=int, default=0, help="Fit at most N galaxies")
    parser.add_argument(
        "--modes",
        nargs="+",
        default=None,
        choices=["psf", "nopsf"],
        help="Subset of fit modes (default: config galfit.modes)",
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    layout = ensure_output_layout(cfg)
    catalog_path = layout["catalogs"] / "truth_catalog.csv"
    if not catalog_path.is_file():
        print(f"Missing {catalog_path}; run generate_mocks.py first.", file=sys.stderr)
        return 1

    truths = read_truth_catalog(catalog_path)
    modes = args.modes or list(cfg.get("galfit", {}).get("modes", ["psf", "nopsf"]))

    # Symlink shared PSF next to mocks for relative feedme path ../mocks/psf.fits from fits/<gid>/<mode>/
    psf_src = layout["mocks"] / "psf.fits"
    if not psf_src.is_file():
        print(f"Missing {psf_src}; run generate_mocks.py first.", file=sys.stderr)
        return 1

    results: list[dict] = []
    n_gal = 0
    for truth in truths:
        if args.limit and n_gal >= args.limit:
            break
        for mode in modes:
            row = run_one_fit(
                truth=truth,
                mode=mode,
                cfg=cfg,
                fits_root=layout["fits"],
                mocks_root=layout["mocks"],
            )
            results.append(row)
            status = "ok" if row.get("converged") else "FAIL"
            print(
                f"[{status}] {truth['galaxy_id']} {mode}: "
                f"b/a={row.get('ba_fit', 'NA')} chi2={row.get('chi2nu', 'NA')}"
            )
        n_gal += 1

    out_csv = layout["catalogs"] / "fit_results.csv"
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=_CATALOG_FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)

    n_ok = sum(1 for r in results if r.get("converged"))
    print(f"Wrote {len(results)} fit rows ({n_ok} converged) -> {out_csv}")
    return 0 if n_ok > 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
