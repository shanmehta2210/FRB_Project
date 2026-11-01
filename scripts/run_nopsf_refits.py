"""
No-PSF control experiment for the 56 seed hosts.

Each host is refit with the exact production GALFIT setup from
``pipeline_scripts/Output/<FRB>_all/`` (cutout, sigma, mask, feedme initial
guesses, constraints) with only ``D)`` set to ``none``. AstroPhot is rerun
through ``checks/astrophot_refit_nopsf.py`` (production check with
``psf_mode='none'``). Production outputs are only read, never written.

Stages:
  galfit     batch GALFIT no-PSF fits (parallel)
  astrophot  batch AstroPhot no-PSF fits (runs in .venv_astrophot)
  table      merge into outputs/nopsf/nopsf_refits.csv
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
VER_DIR = REPO / "pipeline_scripts" / "verification"
for p in (REPO, VER_DIR):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import vercommon as vc  # noqa: E402

SEED_CSV = REPO / "Publication" / "materials" / "host_properties" / "seed_56.csv"
OUTPUT_ROOT = REPO / "pipeline_scripts" / "Output"
NOPSF_ROOT = VER_DIR / "outputs" / "nopsf"
GALFIT_ROOT = NOPSF_ROOT / "galfit"
ASTROPHOT_ROOT = NOPSF_ROOT / "astrophot"
TABLE_CSV = NOPSF_ROOT / "nopsf_refits.csv"
ASTROPHOT_PY = REPO / ".venv_astrophot" / "Scripts" / "python.exe"
ASTROPHOT_SCRIPT = VER_DIR / "checks" / "astrophot_refit_nopsf.py"

INPUT_FILES = ("host_cutout.fits", "host_sigma.fits", "host_mask.fits", "constraints.txt")
GALFIT_TIMEOUT_S = 3600
CRASH_MARKERS = (
    "GALFIT crashed",
    "Singular Matrix",
    "Numerical Recipes run-time error",
    "Segmentation fault",
)
Q0 = 0.2

_PARAM_LINE = re.compile(r"^(\s*)(\d+|Z)\)\s*(.*?)(\s*#.*)?$")


def seed_frbs() -> list[str]:
    return pd.read_csv(SEED_CSV, dtype={"frb": str})["frb"].tolist()


def nopsf_feedme(text: str) -> str:
    out = []
    for line in text.splitlines():
        if line.lstrip().startswith("D)"):
            line = "D) none  # Input PSF file (no-PSF control)"
        out.append(line)
    return "\n".join(out) + "\n"


def seed_from_production(text: str, prod_out: Path) -> str:
    """Replace host Sersic and sky initial values with production best-fit values."""
    hdr = vc.parse_out_header(str(prod_out))
    host = hdr["components"][0]
    sky = hdr["sky"].get("level", float("nan"))
    host_vals = {
        "1": (host["xc"], host["yc"]),
        "3": (host["mag"],),
        "4": (host["re"],),
        "5": (host["n"],),
        "9": (host["q"],),
        "10": (host["pa"],),
    }
    out = []
    comp_type = ""
    sersic_seen = 0
    for line in text.splitlines():
        m = _PARAM_LINE.match(line)
        if m and m.group(2) == "0":
            comp_type = m.group(3).split()[0].lower() if m.group(3).split() else ""
            sersic_seen += comp_type == "sersic"
        elif m:
            num, toks = m.group(2), m.group(3).split()
            vals = None
            if comp_type == "sersic" and sersic_seen == 1 and num in host_vals:
                vals = host_vals[num]
            elif comp_type == "sky" and num == "1" and math.isfinite(sky):
                vals = (sky,)
            if vals is not None and all(math.isfinite(v) for v in vals):
                new = [f"{v:.6g}" for v in vals] + toks[len(vals):]
                line = f"{m.group(1)}{num}) {' '.join(new)}{m.group(4) or ''}"
        out.append(line)
    return "\n".join(out) + "\n"


def run_galfit(wkdir: Path) -> tuple[bool, str]:
    for name in ("fit.log", "out.fits"):
        (wkdir / name).unlink(missing_ok=True)
    for stale in wkdir.glob("galfit.[0-9][0-9]"):
        stale.unlink()
    try:
        proc = subprocess.run(
            ["wsl", "galfit", "galfit.feedme"],
            cwd=wkdir,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            stdin=subprocess.DEVNULL,
            timeout=GALFIT_TIMEOUT_S,
            check=False,
        )
        combined = (proc.stdout or "") + (proc.stderr or "")
    except subprocess.TimeoutExpired:
        return False, f"timeout_{GALFIT_TIMEOUT_S}s"
    (wkdir / "galfit_stdout.log").write_text(combined, encoding="utf-8")
    for marker in CRASH_MARKERS:
        if marker in combined:
            return False, f"crash: {marker}"
    out = wkdir / "out.fits"
    log = wkdir / "fit.log"
    if not (out.is_file() and out.stat().st_size > 0 and log.is_file()):
        return False, f"no_output (exit {proc.returncode})"
    return True, "ok"


def galfit_one(frb: str) -> dict:
    prod = OUTPUT_ROOT / f"{frb}_all"
    wkdir = GALFIT_ROOT / frb
    wkdir.mkdir(parents=True, exist_ok=True)
    for name in INPUT_FILES:
        src = prod / name
        if src.is_file():
            shutil.copy2(src, wkdir / name)
    feed_prod = (prod / "galfit.feedme").read_text(encoding="utf-8")
    feed = nopsf_feedme(feed_prod)
    (wkdir / "galfit.feedme").write_text(feed, encoding="utf-8", newline="\n")

    ok, reason = run_galfit(wkdir)
    status = "ok"
    if not ok:
        retry = seed_from_production(feed, prod / "out.fits")
        (wkdir / "galfit.feedme.pass1").write_text(feed, encoding="utf-8", newline="\n")
        (wkdir / "galfit.feedme").write_text(retry, encoding="utf-8", newline="\n")
        ok, reason2 = run_galfit(wkdir)
        status = "ok_retry_prod_seed" if ok else f"failed: {reason}; retry {reason2}"
    res = {"frb": frb, "status": status}
    if ok:
        hdr = vc.parse_out_header(str(wkdir / "out.fits"))
        host = hdr["components"][0]
        flags = _flagged_params(wkdir / "fit.log")
        res.update(
            {
                "q": host["q"],
                "q_err": host["q_err"],
                "n": host["n"],
                "re_px": host["re"],
                "pa": host["pa"],
                "mag": host["mag"],
                "chi2nu": hdr["chi2nu"],
                "sky": hdr["sky"].get("level"),
                "flags": flags,
            }
        )
    (wkdir / "result.json").write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
    return res


def _flagged_params(log_path: Path) -> str:
    """GALFIT marks parameters with numerical problems as ``*value*``."""
    try:
        lines = log_path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return ""
    names = ("x", "y", "mag", "re", "n", "q", "pa")
    for line in lines:
        if line.lstrip().startswith("sersic") and ":" in line:
            body = line.split(":", 1)[1]
            for ch in "(),":
                body = body.replace(ch, " ")
            return ",".join(nm for nm, t in zip(names, body.split()) if "*" in t)
    return ""


def stage_galfit(frbs: list[str], workers: int) -> None:
    GALFIT_ROOT.mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futs = {pool.submit(galfit_one, f): f for f in frbs}
        for fut in as_completed(futs):
            frb = futs[fut]
            try:
                r = fut.result()
                print(f"[galfit] {frb}: {r['status']} q={r.get('q')} n={r.get('n')}", flush=True)
            except Exception as exc:
                print(f"[galfit] {frb}: exception {type(exc).__name__}: {exc}", flush=True)


def stage_astrophot(frbs: list[str]) -> None:
    subprocess.run([str(ASTROPHOT_PY), str(ASTROPHOT_SCRIPT), *frbs], check=False)


def inclination_deg(q: float, q0: float = Q0) -> float:
    if q is None or not math.isfinite(q):
        return float("nan")
    if q >= 1.0:
        return 0.0
    if q <= q0:
        return 90.0
    return math.degrees(math.acos(math.sqrt((q * q - q0 * q0) / (1.0 - q0 * q0))))


def psf_fwhm_arcsec(frb: str) -> float:
    prod = OUTPUT_ROOT / f"{frb}_all"
    ps = vc.plate_scale_arcsec(frb)
    wcs = vc.parse_psfex_xml(str(prod / "psfex.xml")).get("FWHM_WCS_Mean", float("nan"))
    if math.isfinite(wcs) and wcs > 0:
        return wcs
    fwhm = vc.parse_psfex_xml(str(prod / "psfex.xml")).get("FWHM_FromFluxRadius_Mean", float("nan"))
    if not (math.isfinite(fwhm) and fwhm > 0):
        psf = vc._read_image(str(prod / "proto_image.fits"))
        fwhm = vc.psf_second_moments(psf).get("fwhm_px", float("nan")) if psf is not None else float("nan")
    return fwhm * ps


def _read(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def stage_table() -> pd.DataFrame:
    seed = pd.read_csv(SEED_CSV, dtype={"frb": str})
    rows = []
    for _, s in seed.iterrows():
        frb = s["frb"]
        ps = vc.plate_scale_arcsec(frb)
        g = _read(GALFIT_ROOT / frb / "result.json")
        a = _read(ASTROPHOT_ROOT / f"{frb}.json")
        nan = float("nan")
        re_px = g.get("re_px", nan)
        row = {
            "frb": frb,
            "q_galfit_psf": s["q_winning"],
            "q_galfit_nopsf": g.get("q", nan),
            "n_galfit_nopsf": g.get("n", nan),
            "re_galfit_nopsf": re_px * ps if re_px is not None else nan,
            "chi2nu_galfit_nopsf": g.get("chi2nu", nan),
            "q_astrophot_psf": s["q_astrophot"],
            "q_astrophot_nopsf": a.get("ap_q", nan),
            "n_astrophot_nopsf": a.get("ap_n", nan),
            "status_galfit": g.get("status", "not_run"),
            "status_astrophot": a.get("status", "not_run"),
            "q_galfit_production": s["q_galfit_production"],
            "q_winning_leg": s["q_winning_leg"],
            "n_galfit_psf": s["n_galfit"],
            "re_arcsec_galfit": s["re_arcsec_galfit"],
            "psf_fwhm_arcsec": psf_fwhm_arcsec(frb),
            "galfit_nopsf_flags": g.get("flags", ""),
            "re_astrophot_nopsf": a.get("ap_re_arcsec", nan),
            "chi2nu_astrophot_nopsf": a.get("ap_chi2nu", nan),
        }
        row["re_over_fwhm"] = row["re_arcsec_galfit"] / row["psf_fwhm_arcsec"]
        rows.append(row)
    df = pd.DataFrame(rows)
    for col in ("q_galfit_psf", "q_galfit_nopsf", "q_astrophot_psf", "q_astrophot_nopsf",
                "q_galfit_production"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
        df["i_" + col[2:] + "_deg"] = df[col].map(inclination_deg)
    NOPSF_ROOT.mkdir(parents=True, exist_ok=True)
    df.to_csv(TABLE_CSV, index=False, float_format="%.5g")
    print(f"Wrote {TABLE_CSV} ({len(df)} rows)")
    return df


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=("galfit", "astrophot", "table", "all"))
    parser.add_argument("--frb", nargs="*", help="Subset of FRBs (default: all of seed_56)")
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    frbs = args.frb or seed_frbs()
    if args.stage in ("galfit", "all"):
        stage_galfit(frbs, args.workers)
    if args.stage in ("astrophot", "all"):
        stage_astrophot(frbs)
    if args.stage in ("table", "all"):
        stage_table()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
