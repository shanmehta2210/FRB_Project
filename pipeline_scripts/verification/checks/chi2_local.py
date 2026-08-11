"""Check 1 — global and host-localized chi2/nu.

The global value GALFIT reports covers the whole fitting region. Since the
Re-separation ROI change the host can occupy a small part of the stamp, so sky
pixels dilute host residuals toward 1: a compact galaxy can fit badly and still
report chi2/nu ~ 1. The localized value fixes that; the sigma calibration ratio
says whether either number can be read absolutely at all.
"""

from __future__ import annotations

import math
import os
import re

import numpy as np

import vercommon as vc

NAME = "chi2"

_NDOF = re.compile(r"ndof\s*=\s*(\d+)")
_CHI2 = re.compile(r"Chi\^2\s*=\s*([-0-9.eE+]+)")


def _global_from_fitlog(path: str) -> dict:
    """Last ``Chi^2 = ..., ndof = ...`` pair written by GALFIT."""
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            text = f.read()
    except OSError:
        return {}
    ndofs = _NDOF.findall(text)
    chi2s = _CHI2.findall(text)
    out: dict = {}
    if ndofs:
        out["ndof_global"] = int(ndofs[-1])
    if chi2s:
        try:
            out["chi2_global"] = float(chi2s[-1])
        except ValueError:
            pass
    return out


def _localized(host: vc.HostData, n_re: float) -> dict:
    a, _ = vc.elliptical_coords(host.shape, host.xc, host.yc, host.q, host.pa)
    sel = vc.valid_mask(host) & (a <= n_re * host.re)
    npix = int(np.count_nonzero(sel))
    suffix = f"{n_re:g}re".replace(".", "p")
    if npix == 0:
        return {f"chi2nu_local_{suffix}": float("nan"), f"chi2_local_npix_{suffix}": 0}
    chi2 = float(np.sum((host.resid[sel] / host.sigma[sel]) ** 2))
    # 7 free Sersic parameters per component plus the single sky level.
    k = 7 * (1 + len(host.neighbours)) + 1
    nu = max(npix - k, 1)
    return {
        f"chi2nu_local_{suffix}": chi2 / nu,
        f"chi2_local_npix_{suffix}": npix,
        f"chi2_local_nu_{suffix}": nu,
    }


def run(host: vc.HostData, outdir: str) -> dict:
    out: dict = {
        "chi2nu_global": host.chi2nu_log,
        "re_px": host.re,
        "re_arcsec": host.re_arcsec,
        "psf_fwhm_px": host.psf_fwhm,
        "re_over_fwhm": (host.re / host.psf_fwhm
                         if math.isfinite(host.psf_fwhm) and host.psf_fwhm > 0
                         else float("nan")),
        "n_fit_components": 1 + len(host.neighbours),
        "residual_closure": host.residual_closure,
        "log_vs_header": host.log_vs_header,
    }
    out.update(_global_from_fitlog(os.path.join(host.dir, "fit.log")))
    out.update(_localized(host, 1.0))
    out.update(_localized(host, 2.0))
    out.update(vc.sigma_calibration_ratio(host))

    # chi2 scales as the square of any sigma miscalibration, so these are the
    # values the fit would have reported with a correctly normalized sigma map.
    ratio = out.get("sigma_calibration_ratio", float("nan"))
    if math.isfinite(ratio) and ratio > 0:
        out["chi2nu_global_corrected"] = (
            out.get("chi2nu_global", float("nan")) / ratio**2
        )
        out["chi2nu_local_2re_corrected"] = (
            out.get("chi2nu_local_2re", float("nan")) / ratio**2
        )
    else:
        out["chi2nu_global_corrected"] = float("nan")
        out["chi2nu_local_2re_corrected"] = float("nan")

    out["snr_win"] = float(host.sky_audit.get("snr_win", float("nan")) or float("nan"))
    return out
