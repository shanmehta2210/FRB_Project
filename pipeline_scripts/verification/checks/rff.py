"""Check 2 — localized residual flux fraction.

    RFF = ( sum |I - M| - 0.8 sum sigma ) / sum M_galaxy

The 0.8 is sqrt(2/pi), the expectation of |x| for zero-mean Gaussian noise.
Subtracting it means pure noise gives RFF ~ 0 whatever the depth, which is what
makes the number comparable between a bright host and a faint one. At low SNR
the subtraction is noisy and RFF can land below zero; that is expected, and is
why the analytic error is always reported with it.
"""

from __future__ import annotations

import math

import numpy as np

import vercommon as vc

NAME = "rff"

SQRT_2_OVER_PI = math.sqrt(2.0 / math.pi)
VAR_ABS = 1.0 - 2.0 / math.pi


def _rff(host: vc.HostData, sel: np.ndarray, model_gal: np.ndarray) -> dict:
    npix = int(np.count_nonzero(sel))
    if npix < 10:
        return {"rff": float("nan"), "rff_err": float("nan"), "npix": npix,
                "model_flux": float("nan")}
    absres = np.abs(host.resid[sel])
    sig = host.sigma[sel]
    denom = float(np.sum(model_gal[sel]))
    if not math.isfinite(denom) or denom <= 0:
        return {"rff": float("nan"), "rff_err": float("nan"), "npix": npix,
                "model_flux": denom}
    value = (float(np.sum(absres)) - SQRT_2_OVER_PI * float(np.sum(sig))) / denom
    err = math.sqrt(VAR_ABS * float(np.sum(sig**2))) / denom
    return {"rff": value, "rff_err": err, "npix": npix, "model_flux": denom}


def run(host: vc.HostData, outdir: str) -> dict:
    a, _ = vc.elliptical_coords(host.shape, host.xc, host.yc, host.q, host.pa)
    valid = vc.valid_mask(host)

    # The GALFIT model plane includes the sky component; the RFF denominator
    # must be galaxy flux only, or a large aperture of sky would dilute it.
    sky = host.sky_level if math.isfinite(host.sky_level) else 0.0
    model_gal = host.model - sky

    out: dict = {"sky_level_adu": sky, "re_px": host.re}
    for label, sel in (
        ("1re", valid & (a <= 1.0 * host.re)),
        ("2re", valid & (a <= 2.0 * host.re)),
        ("annulus_1_2re", valid & (a > 1.0 * host.re) & (a <= 2.0 * host.re)),
    ):
        res = _rff(host, sel, model_gal)
        out[f"rff_{label}"] = res["rff"]
        out[f"rff_{label}_err"] = res["rff_err"]
        out[f"rff_{label}_npix"] = res["npix"]
        out[f"rff_{label}_model_flux"] = res["model_flux"]

    for label in ("1re", "2re", "annulus_1_2re"):
        val, err = out[f"rff_{label}"], out[f"rff_{label}_err"]
        out[f"rff_{label}_sig"] = (val / err) if (math.isfinite(err) and err > 0) else float("nan")

    # A large annulus-minus-core difference points at background error rather
    # than a wrong model shape.
    out["rff_outer_minus_inner"] = out["rff_annulus_1_2re"] - out["rff_1re"]
    return out
