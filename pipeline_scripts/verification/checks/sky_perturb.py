"""Check 4 — sky perturbation.

Background error is the dominant systematic for faint and extended galaxies: it
moves outer-disk flux, and outer-disk flux is what sets q. This refits each host
twice with the sky held fixed one sigma above and one sigma below its best-fit
value, and reports how far q moves.

The perturbation is the **large-scale** sky uncertainty, not the per-pixel
noise. SExtractor writes a local background estimate per detected source, so the
scatter of that column across the field measures exactly the structure that
matters — gradients, scattered light, flat-field residuals. Per-pixel sky RMS is
recorded for reference but not used: it is far larger than the uncertainty on
the mean sky and would turn a realistic systematic into a worst case.

Nothing in ``pipeline_scripts/Output/`` is touched. The Phase 3a inputs are
copied into the verification tree and the feedme is edited with a minimal diff,
so everything except the sky line is byte-identical to production.
"""

from __future__ import annotations

import contextlib
import importlib.util
import logging
import math
import os
import re
import shutil
import sys

import numpy as np
from astropy.io import fits

import vercommon as vc

NAME = "sky"

_COMP_HEADER = re.compile(r"^\s*#\s*Component number:\s*(\d+)", re.IGNORECASE)
_PARAM = re.compile(r"^(\s*)(\d+)\)\s*(.*?)(\s*#.*)?$")


def _run_galfit():
    """Load production's ``run_galfit`` so the two paths share crash handling."""
    if "_ver_run_galfit_fitting" in sys.modules:
        return sys.modules["_ver_run_galfit_fitting"].run_galfit
    path = os.path.join(vc.PIPELINE_DIR, "galfit_fitting", "run_galfit_fitting.py")
    spec = importlib.util.spec_from_file_location("_ver_run_galfit_fitting", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    # 128 runs of verbatim GALFIT output would bury the batch progress; it is
    # kept per-leg in galfit_stdout.log instead.
    logging.getLogger("frb_pipeline.phase3b").setLevel(logging.ERROR)
    return module.run_galfit


def sky_sigma_from_catalog(cat_path: str) -> dict:
    """1.4826 * MAD of SExtractor's per-source BACKGROUND over clean detections."""
    if not os.path.isfile(cat_path):
        return {"sky_sigma_adu": float("nan"), "sky_sigma_source": "missing_cat"}
    try:
        with fits.open(cat_path, memmap=False) as hdul:
            tab = None
            for hdu in hdul:
                cols = getattr(getattr(hdu, "columns", None), "names", None)
                if cols and "BACKGROUND" in cols:
                    tab = hdu.data
                    break
            if tab is None:
                return {"sky_sigma_adu": float("nan"), "sky_sigma_source": "no_background_col"}
            bkg = np.asarray(tab["BACKGROUND"], dtype=float)
            flags = (np.asarray(tab["FLAGS"], dtype=float)
                     if "FLAGS" in tab.columns.names else np.zeros_like(bkg))
    except Exception as exc:
        return {"sky_sigma_adu": float("nan"),
                "sky_sigma_source": f"unreadable: {type(exc).__name__}"}

    sel = np.isfinite(bkg) & (flags == 0)
    if np.count_nonzero(sel) < 10:
        sel = np.isfinite(bkg)
    vals = bkg[sel]
    if vals.size < 5:
        return {"sky_sigma_adu": float("nan"), "sky_sigma_source": "too_few_sources"}
    med = float(np.median(vals))
    sigma = 1.4826 * float(np.median(np.abs(vals - med)))
    return {
        "sky_sigma_adu": sigma,
        "sky_sigma_source": "image.cat BACKGROUND MAD",
        "sky_cat_median_adu": med,
        "sky_cat_nsources": int(vals.size),
    }


def _edit_feedme(lines: list[str], comps: list[dict], sky_value: float,
                 sky_comp: int, fix_n: float | None = None) -> list[str]:
    """Reseed every Sersic parameter at its best fit and fix the sky.

    Reseeding at the answer isolates the sky effect: whatever moves is the sky's
    doing, not a different path through the optimizer. If the parent fit held
    host ``n`` fixed, keep it fixed here too (otherwise sky± just re-opens the
    n–sky degeneracy the re-fit closed).
    """
    out = list(lines)
    comp_id = 0
    comp_type = ""
    sersic_seen = 0
    params: dict | None = None

    for i, line in enumerate(lines):
        head = _COMP_HEADER.match(line)
        if head:
            comp_id = int(head.group(1))
            comp_type = ""
            params = None
            continue
        m = _PARAM.match(line)
        if not m:
            continue
        indent, num, body, comment = m.group(1), m.group(2), m.group(3), m.group(4) or ""
        if num == "0":
            comp_type = body.split()[0].lower() if body.split() else ""
            if comp_type == "sersic":
                params = comps[sersic_seen] if sersic_seen < len(comps) else None
                sersic_seen += 1
            continue

        if comp_type == "sky" and comp_id == sky_comp and num == "1":
            out[i] = f"{indent}1) {sky_value:.8f} 0{comment}"
        elif comp_type == "sersic" and params is not None:
            is_host = sersic_seen == 1
            if num == "5" and is_host and fix_n is not None:
                out[i] = f"{indent}5) {float(fix_n):.4f} 0{comment}"
                continue
            new = {
                "1": f"{params['xc']:.4f} {params['yc']:.4f} 1 1",
                "3": f"{params['mag']:.4f} 1",
                "4": f"{params['re']:.4f} 1",
                "5": f"{params['n']:.4f} 1",
                "9": f"{params['q']:.4f} 1",
                "10": f"{params['pa']:.4f} 1",
            }.get(num)
            if new is not None:
                out[i] = f"{indent}{num}) {new}{comment}"
    return out


def _strip_sky_constraints(src: str, dst: str, sky_comp: int) -> None:
    """Drop any constraint on the sky component so the fixed value cannot move."""
    try:
        with open(src, encoding="utf-8", errors="replace") as f:
            lines = f.read().splitlines()
    except OSError:
        lines = []
    kept = [ln for ln in lines
            if not (ln.split() and ln.split()[0] == str(sky_comp))]
    with open(dst, "w", newline="\n", encoding="utf-8") as f:
        f.write("\n".join(kept) + ("\n" if kept else ""))


def _stage(host: vc.HostData, outdir: str, label: str, sky_value: float,
           sky_comp: int, comps: list[dict]) -> str:
    wkdir = os.path.join(outdir, f"galfit_sky_{label}")
    os.makedirs(wkdir, exist_ok=True)
    for name in ("host_cutout.fits", "host_sigma.fits", "host_mask.fits",
                 "proto_image.fits"):
        src = os.path.join(host.dir, name)
        if os.path.isfile(src):
            shutil.copy2(src, os.path.join(wkdir, name))
    _strip_sky_constraints(os.path.join(host.dir, "constraints.txt"),
                           os.path.join(wkdir, "constraints.txt"), sky_comp)

    feed = vc.parse_feedme(os.path.join(host.dir, "galfit.feedme"))
    fix_n = vc.host_n_held_fixed(os.path.join(host.dir, "galfit.feedme"))
    edited = _edit_feedme(feed["lines"], comps, sky_value, sky_comp, fix_n=fix_n)
    with open(os.path.join(wkdir, "galfit.feedme"), "w", newline="\n",
              encoding="utf-8") as f:
        f.write("\n".join(edited) + "\n")
    # If n is held fixed, drop n constraint rows so they cannot fight the flag.
    if fix_n is not None:
        cpath = os.path.join(wkdir, "constraints.txt")
        try:
            with open(cpath, encoding="utf-8", errors="replace") as f:
                clines = f.read().splitlines()
        except OSError:
            clines = []
        kept = [ln for ln in clines
                if not (ln.split() and len(ln.split()) >= 2
                        and ln.split()[1].lower() == "n")]
        with open(cpath, "w", newline="\n", encoding="utf-8") as f:
            f.write("\n".join(kept) + ("\n" if kept else ""))

    for stale in ("fit.log", "out.fits"):
        path = os.path.join(wkdir, stale)
        if os.path.isfile(path):
            os.remove(path)
    return wkdir


def _read_result(wkdir: str) -> dict:
    out_path = os.path.join(wkdir, "out.fits")
    if not os.path.isfile(out_path):
        return {"status": "no_output"}
    try:
        hdr = vc.parse_out_header(out_path)
    except Exception as exc:
        return {"status": f"unreadable: {type(exc).__name__}"}
    if not hdr["components"]:
        return {"status": "no_component"}
    c = hdr["components"][0]
    return {
        "status": "ok",
        "q": c["q"], "pa": c["pa"], "re": c["re"], "n": c["n"], "mag": c["mag"],
        "q_err": c["q_err"],
        "sky": hdr["sky"].get("level", float("nan")),
    }


def run(host: vc.HostData, outdir: str) -> dict:
    out: dict = {"q_galfit": host.q, "pa_galfit_deg": host.pa,
                 "re_galfit_px": host.re, "n_galfit": host.n,
                 "mag_galfit": host.mag, "sky_best_adu": host.sky_level}
    out.update(sky_sigma_from_catalog(os.path.join(host.dir, "image.cat")))

    cal = vc.sigma_calibration_ratio(host)
    out["sky_pixel_rms_adu"] = cal.get("sky_mad_adu", float("nan"))
    sigma = out.get("sky_sigma_adu", float("nan"))
    if not (math.isfinite(sigma) and sigma > 0):
        return {**out, "status": "no_sky_sigma"}
    if not math.isfinite(host.sky_level):
        return {**out, "status": "no_best_fit_sky"}
    out["sky_sigma_over_pixel_rms"] = (
        sigma / out["sky_pixel_rms_adu"]
        if math.isfinite(out["sky_pixel_rms_adu"]) and out["sky_pixel_rms_adu"] > 0
        else float("nan")
    )

    hdr = vc.parse_out_header(os.path.join(host.dir, "out.fits"))
    comps = hdr["components"]
    sky_comp = hdr["sky"].get("comp", len(comps) + 1)

    run_galfit = _run_galfit()
    legs: dict[str, dict] = {}
    for label, sign in (("plus", +1.0), ("minus", -1.0)):
        wkdir = _stage(host, outdir, label, host.sky_level + sign * sigma,
                       sky_comp, comps)
        ok = False
        try:
            with open(os.path.join(wkdir, "galfit_stdout.log"), "w",
                      encoding="utf-8") as fh, contextlib.redirect_stdout(fh):
                ok = bool(run_galfit(wkdir))
        except Exception as exc:
            legs[label] = {"status": f"launch_failed: {type(exc).__name__}: {exc}"}
        if label not in legs:
            res = _read_result(wkdir)
            res["converged"] = ok
            legs[label] = res
        out[f"sky_{label}_adu"] = host.sky_level + sign * sigma
        for key in ("q", "pa", "re", "n", "mag", "status", "converged"):
            out[f"{key}_sky_{label}"] = legs[label].get(key, float("nan"))

    shifts: dict[str, list[float]] = {}
    for key, base in (("q", host.q), ("pa", host.pa), ("re", host.re),
                      ("n", host.n), ("mag", host.mag)):
        vals = []
        for label in ("plus", "minus"):
            v = legs[label].get(key, float("nan"))
            if isinstance(v, (int, float)) and math.isfinite(v):
                delta = (vc.wrap_pa(v - base) if key == "pa" else v - base)
                vals.append(float(delta))
        shifts[key] = vals
        out[f"d{key}_sky"] = max((abs(v) for v in vals), default=float("nan"))
        out[f"d{key}_sky_signed"] = (float(np.mean(vals)) if vals else float("nan"))

    both = all(legs[label].get("status") == "ok" for label in ("plus", "minus"))
    out["status"] = "ok" if both else "partial"
    # Report the triad explicitly: production q0 alongside the two sky legs.
    out["q_sky_0"] = float(host.q)
    out["q_sky_plus"] = float(legs["plus"].get("q", float("nan")))
    out["q_sky_minus"] = float(legs["minus"].get("q", float("nan")))
    # The comparison that matters: is the sky systematic bigger than the
    # statistical error GALFIT quotes on q?
    out["dq_sky_over_q_err"] = (out["dq_sky"] / host.q_err
                                if math.isfinite(host.q_err) and host.q_err > 0
                                else float("nan"))
    return out
