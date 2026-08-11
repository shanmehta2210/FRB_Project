"""Check 5 — free isophote fits on the data and on the PSF-convolved model.

An independent, non-parametric view of the same geometry, and the only check
that says *where in radius* a single-Sersic description starts to fail.

Both fits are run with q and PA free at every isophote, from the same starting
geometry, with the same mask, on the same radial grid. Comparing against the
**convolved** model is essential: an unconvolved model would attribute the PSF's
own rounding of the inner isophotes to a fit error. Applying the identical mask
to both sides means any quirk of the masking cancels in the difference.
"""

from __future__ import annotations

import math
import os
import warnings

import numpy as np
from photutils.isophote import Ellipse, EllipseGeometry

import vercommon as vc

NAME = "isophote"

STEP = 0.1
NCLIP = 3
MIN_RE_OVER_FWHM = 0.8


GOOD_ENOUGH = 8


def _masked(image: np.ndarray, mask: np.ndarray) -> np.ma.MaskedArray:
    return np.ma.masked_array(np.nan_to_num(image, nan=0.0), mask=mask)


def _fit_with(image, mask, st: dict, minsma: float, maxsma: float):
    """One ellipse fit at a fixed strategy. Returns ``(isolist, error)``."""
    geom = EllipseGeometry(x0=st["x0"], y0=st["y0"], sma=st["sma0"],
                           eps=st["eps"], pa=st["pa"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            iso = Ellipse(_masked(image, mask), geometry=geom).fit_image(
                minsma=minsma, maxsma=maxsma, step=st["step"], nclip=NCLIP,
                fix_center=st["fix_center"], fix_pa=False, fix_eps=False,
            )
        except Exception as exc:
            return None, f"{type(exc).__name__}: {exc}"
    return (iso, "") if len(iso) >= 3 else (None, f"only {len(iso)} isophotes")


def _search(image, mask, *, x0: float, y0: float, eps0: float, pa0: float,
            minsma: float, maxsma: float):
    """Scan starting geometries until the fit takes.

    ``fit_image`` grows outward and inward from ``sma0`` and gives up entirely
    if that first isophote does not converge, so the result is startlingly
    sensitive to it: on 20240119A, sma0 = 8 px yields 21 isophotes while 5 and
    12 yield none. A single starting guess is therefore not a fit strategy.
    The best (most isophotes) result is kept, not merely the first success.
    """
    seeds = np.geomspace(max(minsma * 1.2, minsma + 0.5), max(maxsma * 0.9,
                                                              minsma + 2.0), 12)
    best, best_st, last = None, None, ""
    for step, fix_center in ((STEP, True), (0.15, True), (STEP, False)):
        for sma0 in seeds:
            st = {"x0": x0, "y0": y0, "sma0": float(sma0), "eps": eps0,
                  "pa": pa0, "step": step, "fix_center": fix_center}
            iso, err = _fit_with(image, mask, st, minsma, maxsma)
            if iso is None:
                last = err
                continue
            if best is None or len(iso) > len(best):
                best, best_st = iso, st
            if best is not None and len(best) >= GOOD_ENOUGH:
                return best, best_st, ""
        if best is not None:
            return best, best_st, ""
    return None, None, last or "no isophotes"


def _arrays(iso) -> dict:
    eps = np.asarray(iso.eps, dtype=float)
    return {
        "sma": np.asarray(iso.sma, dtype=float),
        "q": 1.0 - eps,
        "q_err": np.asarray(iso.ellip_err, dtype=float),
        # photutils PA is measured from +x; GALFIT measures it from +y.
        "pa_deg": np.array([vc.wrap_pa(math.degrees(p) - 90.0) for p in
                            np.asarray(iso.pa, dtype=float)]),
        "pa_err_deg": np.degrees(np.asarray(iso.pa_err, dtype=float)),
        "intens": np.asarray(iso.intens, dtype=float),
        "intens_err": np.asarray(iso.int_err, dtype=float),
    }


def _interp(x_new: np.ndarray, x: np.ndarray, y: np.ndarray) -> np.ndarray:
    ok = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(ok) < 2:
        return np.full_like(x_new, np.nan, dtype=float)
    order = np.argsort(x[ok])
    return np.interp(x_new, x[ok][order], y[ok][order], left=np.nan, right=np.nan)


def run(host: vc.HostData, outdir: str) -> dict:
    sky = host.sky_level if math.isfinite(host.sky_level) else 0.0
    bad = ~vc.valid_mask(host)
    minsma = max(2.0, 0.2 * host.re)
    maxsma = min(3.0 * host.re, 0.45 * min(host.shape))
    if not math.isfinite(maxsma) or maxsma <= minsma + 2.0:
        return {"status": "stamp_too_small", "iso_minsma": minsma, "iso_maxsma": maxsma}

    eps0 = float(np.clip(1.0 - host.q, 0.02, 0.90))
    pa0 = math.radians(host.pa + 90.0)
    rof = (host.re / host.psf_fwhm if math.isfinite(host.psf_fwhm)
           and host.psf_fwhm > 0 else float("nan"))

    geom_kw = {"x0": host.xc, "y0": host.yc, "eps0": eps0, "pa0": pa0,
               "minsma": minsma, "maxsma": maxsma}
    iso_d, st_d, err_d = _search(host.data - sky, bad, **geom_kw)
    if iso_d is None:
        # A source smaller than the PSF has no isophotes to fit; that is a
        # property of the data, not a failure of the check.
        unresolved = math.isfinite(rof) and rof < MIN_RE_OVER_FWHM
        return {"status": "unresolved" if unresolved else "isophote_fit_failed",
                "iso_error_data": err_d, "iso_minsma": minsma,
                "iso_maxsma": maxsma, "re_over_fwhm": rof}

    # Fit the model with the strategy that worked on the data, so the two
    # profiles are produced identically; only fall back to its own search if
    # that fails outright.
    iso_m, err_m = _fit_with(host.model - sky, bad, st_d, minsma, maxsma)
    same_strategy = iso_m is not None
    if iso_m is None:
        iso_m, _st_m, err_m = _search(host.model - sky, bad, **geom_kw)
    if iso_m is None:
        return {"status": "isophote_fit_failed", "iso_error_model": err_m,
                "iso_minsma": minsma, "iso_maxsma": maxsma}

    d = _arrays(iso_d)
    m = _arrays(iso_m)

    # Compare on the data's radial grid; the model grid can terminate elsewhere.
    sma = d["sma"]
    m_q = _interp(sma, m["sma"], m["q"])
    m_qe = _interp(sma, m["sma"], m["q_err"])
    m_pa = _interp(sma, m["sma"], m["pa_deg"])
    m_pae = _interp(sma, m["sma"], m["pa_err_deg"])
    m_int = _interp(sma, m["sma"], m["intens"])

    dq = d["q"] - m_q
    dq_err = np.sqrt(d["q_err"] ** 2 + m_qe**2)
    dpa = np.array([vc.wrap_pa(a - b) for a, b in zip(d["pa_deg"], m_pa)])
    # Floored at 0.5 deg: photutils returns zero PA errors on isophotes it
    # considers perfectly determined, which would otherwise dominate the mean.
    dpa_err = np.sqrt(np.maximum(d["pa_err_deg"], 0.5) ** 2
                      + np.maximum(m_pae, 0.5) ** 2)

    a_re = sma / host.re
    inner = np.isfinite(a_re) & (a_re <= 2.0)
    dq_mean, dq_mean_err = vc.weighted_mean(dq[inner], dq_err[inner])
    dpa_mean, dpa_mean_err = vc.weighted_mean(dpa[inner], dpa_err[inner])

    # Radius at which the data and the model isophotes part company: where a
    # single Sersic stops describing the galaxy.
    signif = np.abs(dq) / np.where(dq_err > 0, dq_err, np.nan)
    broke = np.where(np.isfinite(signif) & (signif > 3.0))[0]
    break_re = float(a_re[broke[0]]) if broke.size else float("nan")

    mu_d = vc.counts_to_mu(d["intens"], host.magzpt, host.plate_scale)
    mu_m = vc.counts_to_mu(m_int, host.magzpt, host.plate_scale)
    with np.errstate(divide="ignore", invalid="ignore"):
        mu_err = 1.0857 * d["intens_err"] / d["intens"]
    # Negative intensity in the faint outskirts gives a negative "error"; those
    # isophotes have no surface brightness to quote, so drop them outright.
    mu_ok = (np.isfinite(mu_err) & (mu_err > 0) & (mu_err <= 0.3)
             & np.isfinite(mu_d) & (d["intens"] > 0))

    np.savez_compressed(
        os.path.join(outdir, "isophote_profiles.npz"),
        sma=sma, a_re=a_re,
        q_data=d["q"], q_data_err=d["q_err"],
        pa_data_deg=d["pa_deg"], pa_data_err_deg=d["pa_err_deg"],
        q_model=m_q, q_model_err=m_qe, pa_model_deg=m_pa,
        dq=dq, dq_err=dq_err, dpa_deg=dpa, dpa_err_deg=dpa_err,
        intens_data=d["intens"], intens_data_err=d["intens_err"], intens_model=m_int,
        mu_data=mu_d, mu_model=mu_m, mu_err=mu_err, mu_valid=mu_ok,
        sma_model=m["sma"], q_model_native=m["q"], pa_model_native_deg=m["pa_deg"],
    )

    return {
        "iso_dq_2re": dq_mean,
        "iso_dq_2re_err": dq_mean_err,
        "iso_dq_2re_sig": (dq_mean / dq_mean_err
                           if math.isfinite(dq_mean_err) and dq_mean_err > 0
                           else float("nan")),
        "iso_dpa_2re_deg": dpa_mean,
        "iso_dpa_2re_err_deg": dpa_mean_err,
        "iso_q_at_1re_data": float(_interp(np.array([host.re]), sma, d["q"])[0]),
        "iso_q_at_1re_model": float(_interp(np.array([host.re]), sma, m_q)[0]),
        "iso_q_at_2re_data": float(_interp(np.array([2 * host.re]), sma, d["q"])[0]),
        "iso_pa_at_1re_data_deg": float(_interp(np.array([host.re]), sma, d["pa_deg"])[0]),
        "iso_break_radius_re": break_re,
        "iso_dq_max_abs_inner": (float(np.nanmax(np.abs(dq[inner])))
                                 if np.any(inner) else float("nan")),
        "iso_frac_discrepant_inner": (float(np.mean(signif[inner] > 3.0))
                                      if np.any(inner) else float("nan")),
        "iso_n_data": int(len(sma)),
        "iso_n_model": int(len(m["sma"])),
        "iso_same_strategy": bool(same_strategy),
        "iso_sma0": float(st_d["sma0"]),
        "iso_step": float(st_d["step"]),
        "iso_fix_center": bool(st_d["fix_center"]),
        "iso_minsma": minsma,
        "iso_maxsma": maxsma,
        "iso_maxsma_re": maxsma / host.re,
        "iso_mu_valid_max_re": float(np.max(a_re[mu_ok])) if np.any(mu_ok) else float("nan"),
        "q_galfit": host.q,
        "pa_galfit_deg": vc.wrap_pa(host.pa),
        "fwhm_over_re": (host.psf_fwhm / host.re
                         if math.isfinite(host.psf_fwhm) and host.re > 0
                         else float("nan")),
    }
