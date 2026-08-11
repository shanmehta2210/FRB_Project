"""Check 3 — Fourier decomposition of the residual.

A geometry error does not scatter residuals randomly: it writes a specific
angular signature. An axis-ratio error appears as cos(2 theta) (lobes on the
axes), a position-angle error as sin(2 theta) (lobes at 45 degrees). The two are
orthogonal, so A_2 and B_2 measure them independently — and because the
amplitude carries the model's own radial gradient, the result converts straight
back into dq and dPA in the units of the parameters themselves.

Derivation (full version in FIT_VERIFICATION_CHECKS.md):

    da/dq   = -(a / 2q) (1 - cos 2t)      =>  dq   = -2 q A_2 / (a |f'(a)|)
    da/dphi =  (a / 2) (1/q - q) sin 2t   =>  dphi = -2 B_2 / (a |f'(a)| (1/q - q))

That closed form assumes the PSF does not change the conversion, and it does:
convolution damps a quadrupole more than it damps the radial gradient, so the
analytic version under-recovers an injected dq by ~18% at this cohort's
resolution. The suite therefore calibrates the conversion numerically. The
host's Sersic is rebuilt analytically (it matches GALFIT's own model plane to
under 0.1% of peak), perturbed by +/- eps in q and in PA, convolved with the
same ``proto_image.fits``, and the response is decomposed on exactly the same
annuli, pixels and weights as the residual:

    K_q = [M_conv(q+eps) - M_conv(q-eps)] / 2 eps   ->  c_q(a) = its cos 2t term
    dq(a) = A_2(a) / c_q(a)

Injection test: dq = 0.020 -> 0.0193, dq = 0.050 -> 0.0459, against 0.0169 and
0.0402 for the analytic form. The residual few percent is genuine second-order
nonlinearity (dq = 0.05 is an 8% change in q), and shrinks with the
perturbation; dPA recovers 3.000 deg as 2.996. The analytic form is still
reported as ``*_radial`` for comparison, and is the fallback when the rebuild
is unavailable.
"""

from __future__ import annotations

import math
import os

import numpy as np

import vercommon as vc

NAME = "fourier"

M_MAX = 4
MIN_PIX = 20
EPS_Q = 0.02
EPS_PA_DEG = 2.0

# Conditions under which dq is a measurement rather than arithmetic.
MAX_RECON_ERR = 0.03      # rebuilt model within 3% of peak of GALFIT's own
MIN_ANNULI = 4            # enough rings inside 2 Re to average over
MIN_RE_OVER_FWHM = 0.8    # some resolved structure to decompose
# Crosstalk is still *reported* (fourier_kernel_crosstalk) but no longer gates
# reliability: a large PA error usually implies q is untrustworthy too, so
# excluding dq when modes mix would hide the hosts that most need a look.


def _wls(design: np.ndarray, y: np.ndarray, w: np.ndarray):
    """Inverse-variance weighted least squares.

    Least squares rather than an FFT because masking makes the azimuthal
    sampling non-uniform, which an FFT would silently alias.
    """
    xtw = design.T * w
    xtwx = xtw @ design
    try:
        cov = np.linalg.inv(xtwx)
    except np.linalg.LinAlgError:
        return None, None
    return cov @ (xtw @ y), np.sqrt(np.clip(np.diag(cov), 0.0, None))


def decompose(
    resid: np.ndarray,
    model_gal: np.ndarray,
    sigma: np.ndarray,
    valid: np.ndarray,
    *,
    xc: float,
    yc: float,
    q: float,
    pa: float,
    re: float,
    fwhm: float,
    kernel_q: np.ndarray | None = None,
    kernel_pa: np.ndarray | None = None,
    m_max: int = M_MAX,
    a_min_re: float = 0.2,
    a_max_re: float = 3.0,
) -> dict:
    """Per-annulus Fourier coefficients and the dq / dPA they imply.

    ``kernel_q`` / ``kernel_pa`` are dM/dq and dM/dPA (per degree) for the
    convolved model. When supplied they are decomposed with the identical
    design matrix and weights, so the ratio of coefficients is the parameter
    error in the exact least-squares sense. Deliberately free of ``HostData``
    so the unit test can drive it with synthetic images.
    """
    a_map, theta = vc.elliptical_coords(resid.shape, xc, yc, q, pa)
    annuli = vc.azimuthal_annuli(a_map, re, fwhm, a_min_re, a_max_re)
    n_ann = len(annuli)
    n_coef = 1 + 2 * m_max

    a_mid = np.full(n_ann, np.nan)
    coef = np.full((n_ann, n_coef), np.nan)
    coef_err = np.full((n_ann, n_coef), np.nan)
    ker_q = np.full((n_ann, n_coef), np.nan)
    ker_pa = np.full((n_ann, n_coef), np.nan)
    prof = np.full(n_ann, np.nan)
    npix = np.zeros(n_ann, dtype=int)

    for i, (lo, hi) in enumerate(annuli):
        a_mid[i] = 0.5 * (lo + hi)
        sel = valid & (a_map >= lo) & (a_map < hi)
        npix[i] = int(np.count_nonzero(sel))
        if npix[i] < MIN_PIX:
            continue
        t = theta[sel]
        design = np.empty((npix[i], n_coef))
        design[:, 0] = 1.0
        for m in range(1, m_max + 1):
            design[:, 2 * m - 1] = np.cos(m * t)
            design[:, 2 * m] = np.sin(m * t)
        w = 1.0 / sigma[sel] ** 2
        beta, err = _wls(design, resid[sel], w)
        if beta is None:
            continue
        coef[i] = beta
        coef_err[i] = err
        prof[i] = float(np.mean(model_gal[sel]))
        if kernel_q is not None:
            bq, _ = _wls(design, kernel_q[sel], w)
            if bq is not None:
                ker_q[i] = bq
        if kernel_pa is not None:
            bp, _ = _wls(design, kernel_pa[sel], w)
            if bp is not None:
                ker_pa[i] = bp

    # Numerical gradient of the azimuthally averaged convolved model, computed
    # on the finite subset so one empty annulus does not poison the rest.
    dmda = np.full(n_ann, np.nan)
    good = np.isfinite(prof) & np.isfinite(a_mid)
    if np.count_nonzero(good) >= 2:
        dmda[good] = np.gradient(prof[good], a_mid[good])

    a2, b2 = coef[:, 3], coef[:, 4]
    a2e, b2e = coef_err[:, 3], coef_err[:, 4]

    with np.errstate(divide="ignore", invalid="ignore"):
        # Analytic (PSF-blind) conversion, kept for comparison.
        scale = a_mid * np.abs(dmda)
        scale = np.where(scale > 0, scale, np.nan)
        dq_rad = -2.0 * q * a2 / scale
        dq_rad_err = np.abs(2.0 * q * a2e / scale)
        pa_lever = (1.0 / q - q) if q < 1.0 else np.nan
        dpa_rad_ = np.degrees(-2.0 * b2 / (scale * pa_lever))
        dpa_rad_err = np.abs(np.degrees(2.0 * b2e / (scale * pa_lever)))

        # Calibrated conversion.
        cq, cpa = ker_q[:, 3], ker_pa[:, 4]
        cq = np.where(np.abs(cq) > 0, cq, np.nan)
        cpa = np.where(np.abs(cpa) > 0, cpa, np.nan)
        dq_cal = a2 / cq
        dq_cal_err = np.abs(a2e / cq)
        dpa_cal = b2 / cpa
        dpa_cal_err = np.abs(b2e / cpa)
        # Fraction of K_q that leaked into sin2θ. Ideal ~0 (K_q even, K_PA odd).
        # Reported always; does not gate reliability (see FIT_VERIFICATION_CHECKS).
        crosstalk = np.abs(ker_q[:, 4] / cq)

        amp = np.full((m_max + 1, n_ann), np.nan)
        for m in range(1, m_max + 1):
            amp[m] = np.hypot(coef[:, 2 * m - 1], coef[:, 2 * m]) / np.abs(prof)

    calibrated = np.any(np.isfinite(dq_cal))
    dq = dq_cal if calibrated else dq_rad
    dq_err = dq_cal_err if calibrated else dq_rad_err
    dpa = dpa_cal if calibrated else dpa_rad_
    dpa_err = dpa_cal_err if calibrated else dpa_rad_err

    # Phase of the m=2 pattern: R_2 = C cos(2(theta - psi_2)), so psi_2 is
    # defined modulo 180 degrees and must be unwrapped on 2*psi_2.
    psi2 = 0.5 * np.arctan2(b2, a2)
    with np.errstate(divide="ignore", invalid="ignore"):
        psi2_err = 0.5 * np.sqrt(a2**2 * b2e**2 + b2**2 * a2e**2) / (a2**2 + b2**2)

    return {
        "a_mid": a_mid,
        "a_mid_re": a_mid / re,
        "npix": npix,
        "coef": coef,
        "coef_err": coef_err,
        "kernel_q_coef": ker_q,
        "kernel_pa_coef": ker_pa,
        "model_profile": prof,
        "dmda": dmda,
        "dq": dq,
        "dq_err": dq_err,
        "dpa_deg": dpa,
        "dpa_err_deg": dpa_err,
        "dq_radial": dq_rad,
        "dq_radial_err": dq_rad_err,
        "dpa_radial_deg": dpa_rad_,
        "crosstalk": crosstalk,
        "amp": amp,
        "psi2_rad": psi2,
        "psi2_err_rad": psi2_err,
        "calibrated": bool(calibrated),
    }


def response_kernels(host: vc.HostData, eps_q: float = EPS_Q,
                     eps_pa: float = EPS_PA_DEG):
    """dM/dq and dM/dPA (per degree) for the PSF-convolved host model."""
    eq = min(eps_q, 0.4 * host.q, 0.4 * (1.0 - host.q)) if host.q < 1 else eps_q
    eq = max(eq, 1e-3)
    kw = {"include_sky": False, "include_neighbours": False}
    kq = (vc.build_model(host, q=host.q + eq, **kw)
          - vc.build_model(host, q=host.q - eq, **kw)) / (2.0 * eq)
    kpa = (vc.build_model(host, pa=host.pa + eps_pa, **kw)
           - vc.build_model(host, pa=host.pa - eps_pa, **kw)) / (2.0 * eps_pa)
    return kq, kpa, eq


def _phase_winding(a_re: np.ndarray, psi2: np.ndarray, psi2_err: np.ndarray,
                   sel: np.ndarray) -> dict:
    """Linear fit of the unwrapped m=2 phase against radius.

    Flat phase with smoothly growing amplitude means a geometry error; phase
    that winds with radius is, by definition, a spiral.
    """
    nan3 = {"fourier_m2_phase_slope_deg_per_re": float("nan"),
            "fourier_m2_phase_slope_sig": float("nan"),
            "fourier_m2_phase_scatter_deg": float("nan")}
    ok = sel & np.isfinite(psi2) & np.isfinite(psi2_err) & (psi2_err > 0)
    if np.count_nonzero(ok) < 3:
        return nan3
    x = a_re[ok]
    y = np.degrees(np.unwrap(2.0 * psi2[ok]) / 2.0)
    w = 1.0 / np.degrees(psi2_err[ok]) ** 2
    sw, sx, sxx = np.sum(w), np.sum(w * x), np.sum(w * x * x)
    sy, sxy = np.sum(w * y), np.sum(w * x * y)
    det = sw * sxx - sx * sx
    if abs(det) < 1e-30:
        return nan3
    slope = (sw * sxy - sx * sy) / det
    slope_err = math.sqrt(sw / det)
    intercept = (sxx * sy - sx * sxy) / det
    return {
        "fourier_m2_phase_slope_deg_per_re": float(slope),
        "fourier_m2_phase_slope_sig": float(slope / slope_err) if slope_err > 0 else float("nan"),
        "fourier_m2_phase_scatter_deg": float(np.std(y - (slope * x + intercept))),
    }


def run(host: vc.HostData, outdir: str) -> dict:
    sky = host.sky_level if math.isfinite(host.sky_level) else 0.0
    recon = vc.model_reconstruction_error(host)
    try:
        kq, kpa, eps_q = response_kernels(host)
    except Exception as exc:
        kq = kpa = None
        eps_q = float("nan")
        recon["model_recon_status"] = f"kernels failed: {type(exc).__name__}: {exc}"

    prof = decompose(
        host.resid,
        host.model - sky,
        host.sigma,
        vc.valid_mask(host),
        xc=host.xc, yc=host.yc, q=host.q, pa=host.pa,
        re=host.re, fwhm=host.psf_fwhm,
        kernel_q=kq, kernel_pa=kpa,
    )

    a_re = prof["a_mid_re"]
    inner = np.isfinite(a_re) & (a_re <= 2.0)
    outer = np.isfinite(a_re) & (a_re > 2.0)

    dq, dq_e = vc.weighted_mean(prof["dq"][inner], prof["dq_err"][inner])
    dpa, dpa_e = vc.weighted_mean(prof["dpa_deg"][inner], prof["dpa_err_deg"][inner])
    dq_r, dq_r_e = vc.weighted_mean(prof["dq_radial"][inner], prof["dq_radial_err"][inner])

    out: dict = {
        "fourier_dq": dq,
        "fourier_dq_err": dq_e,
        "fourier_dq_sig": (dq / dq_e) if (math.isfinite(dq_e) and dq_e > 0) else float("nan"),
        "fourier_dpa_deg": dpa,
        "fourier_dpa_err_deg": dpa_e,
        "fourier_dpa_sig": (dpa / dpa_e) if (math.isfinite(dpa_e) and dpa_e > 0) else float("nan"),
        "fourier_dq_radial": dq_r,
        "fourier_dq_radial_err": dq_r_e,
        "fourier_calibrated": prof["calibrated"],
        "fourier_eps_q": eps_q,
        "fourier_n_annuli": int(np.count_nonzero(prof["npix"] >= MIN_PIX)),
        "fourier_n_annuli_inner": int(np.count_nonzero(inner & (prof["npix"] >= MIN_PIX))),
        "q_galfit": host.q,
        "q_fourier_corrected": (host.q + dq) if math.isfinite(dq) else float("nan"),
    }
    out.update(recon)
    out.update(_phase_winding(a_re, prof["psi2_rad"], prof["psi2_err_rad"], inner))

    ct = prof["crosstalk"][inner]
    ct = ct[np.isfinite(ct)]
    out["fourier_kernel_crosstalk"] = float(np.max(ct)) if ct.size else float("nan")

    # Normalized amplitudes: m=1 is lopsidedness/centroid drift, m=4 is
    # boxy/disky isophotes (a physical property, not a fit failure).
    for m in (1, 2, 3, 4):
        vals = prof["amp"][m][inner]
        vals = np.abs(vals[np.isfinite(vals)])
        out[f"fourier_m{m}_amp_max"] = float(np.max(vals)) if vals.size else float("nan")
        out[f"fourier_m{m}_amp_median"] = float(np.median(vals)) if vals.size else float("nan")

    # A_0 beyond 2 Re is an independent residual sky offset, cross-checking the
    # sky perturbation test from a completely different direction. Median over
    # annuli, not mean: one bad ring should not set the answer.
    cal = vc.sigma_calibration_ratio(host)
    sky_mad = cal.get("sky_mad_adu", float("nan"))
    a0_outer = prof["coef"][:, 0][outer]
    a0_outer = a0_outer[np.isfinite(a0_outer)]
    a0_med = float(np.median(a0_outer)) if a0_outer.size >= 2 else float("nan")
    out["fourier_a0_outer_adu"] = a0_med
    out["fourier_a0_n_annuli"] = int(a0_outer.size)
    out["fourier_a0_sky_offset_sigma"] = (
        a0_med / sky_mad
        if (math.isfinite(sky_mad) and sky_mad > 0 and math.isfinite(a0_med))
        else float("nan")
    )

    # The estimator needs resolved structure and a rebuilt model that matches
    # GALFIT. Crosstalk is recorded but does not gate: see FIT_VERIFICATION_CHECKS.
    recon_err = out.get("model_recon_max_frac", float("nan"))
    rof = host.re / host.psf_fwhm if (math.isfinite(host.psf_fwhm)
                                      and host.psf_fwhm > 0) else float("nan")
    out["re_over_fwhm"] = rof
    reasons = []
    if not prof["calibrated"]:
        reasons.append("uncalibrated")
    if not (math.isfinite(recon_err) and recon_err < MAX_RECON_ERR):
        reasons.append("model_rebuild_poor")
    if out["fourier_n_annuli_inner"] < MIN_ANNULI:
        reasons.append("too_few_annuli")
    if not (math.isfinite(rof) and rof >= MIN_RE_OVER_FWHM):
        reasons.append("unresolved")
    out["fourier_reliable"] = not reasons
    out["fourier_unreliable_reasons"] = ",".join(reasons)

    np.savez_compressed(
        os.path.join(outdir, "fourier_profiles.npz"),
        **{k: np.asarray(v) for k, v in prof.items()},
    )
    return out
