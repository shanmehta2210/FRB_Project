"""Check 8 — independent Sersic refit with AstroPhot.

Agreement between two codes with different optimizers and different likelihood
implementations rules out a whole class of implementation-specific and
local-minimum failures that no residual metric can see.

The PSF is ``proto_image.fits`` — the exact stamp GALFIT convolved with — so
this compares two fits, not two PSF treatments. ``host_mask.fits`` is applied.
AstroPhot is left to find its own starting point for free parameters: seeding
it with GALFIT's answer would make agreement much cheaper than it looks.

Constraints match the GALFIT feedme for this host:
* sky is always subtracted at GALFIT's level (AstroPhot's flat-sky
  parametrization cannot represent a zero/negative background);
* if GALFIT held host ``n`` fixed (fit flag 0), AstroPhot locks the same ``n``.
  Leaving ``n`` free while GALFIT froze it re-opens the n–sky–size degeneracy
  and is not a fair cross-check.
"""

from __future__ import annotations

import contextlib
import io
import math
import os
import time

import numpy as np

import vercommon as vc

NAME = "astrophot"


def _value(param) -> float:
    try:
        return float(param.value.item())
    except Exception:
        try:
            return float(np.asarray(param.value).reshape(-1)[0])
        except Exception:
            return float("nan")


def _uncertainty(param) -> float:
    try:
        return float(param.uncertainty.item())
    except Exception:
        return float("nan")


def run(host: vc.HostData, outdir: str) -> dict:
    import astrophot as ap

    out: dict = {
        "q_galfit": host.q,
        "pa_galfit_deg": vc.wrap_pa(host.pa),
        "re_galfit_px": host.re,
        "n_galfit": host.n,
        "mag_galfit": host.mag,
        "chi2nu_galfit": host.chi2nu_log,
    }

    variance = np.asarray(host.sigma, dtype=np.float64) ** 2
    bad = ~np.isfinite(variance) | (variance <= 0)
    if np.any(bad):
        good = variance[~bad]
        # Down-weight rather than drop: AstroPhot needs a finite variance
        # everywhere, and the mask already removes what should be ignored.
        variance[bad] = (np.nanmax(good) * 10.0) if good.size else 1.0

    mask = host.mask | ~np.isfinite(host.data)
    psf = np.clip(np.nan_to_num(np.asarray(host.psf, dtype=np.float64)), 0.0, None)
    if psf.sum() <= 0 or psf.size < 4:
        return {**out, "status": "no_psf_stamp"}
    psf = psf / psf.sum()
    # AstroPhot wants an odd-sized PSF centred on a pixel.
    if psf.shape[0] % 2 == 0 or psf.shape[1] % 2 == 0:
        psf = psf[: psf.shape[0] - (1 - psf.shape[0] % 2),
                  : psf.shape[1] - (1 - psf.shape[1] % 2)]
        psf = psf / psf.sum()

    ps = host.plate_scale
    sky = host.sky_level if math.isfinite(host.sky_level) else 0.0
    fix_n = vc.host_n_held_fixed(os.path.join(host.dir, "galfit.feedme"))
    model_kw: dict = {
        "name": "host",
        "model_type": "sersic galaxy model",
        "psf_mode": "full",
    }
    if fix_n is not None:
        model_kw["n"] = {"value": float(fix_n), "locked": True}

    t0 = time.time()
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            target = ap.image.Target_Image(
                data=np.asarray(host.data, dtype=np.float64) - sky,
                pixelscale=ps,
                variance=variance,
                mask=mask,
                psf=ap.image.PSF_Image(data=psf, pixelscale=ps),
            )
            galaxy = ap.models.AstroPhot_Model(target=target, **model_kw)
            galaxy.initialize()
            if fix_n is not None:
                # Belt-and-suspenders: some AstroPhot versions unlock on init.
                try:
                    galaxy["n"].value = float(fix_n)
                    galaxy["n"].locked = True
                except Exception:
                    pass
            result = ap.fit.LM(galaxy, verbose=0).fit()
            # AstroPhot's LM tracks the reduced chi2, so this is directly
            # comparable to GALFIT's Chi^2/nu (do not divide by ndf again).
            chi2nu = float(result.res_loss())
    except Exception as exc:
        return {**out, "status": f"astrophot_failed: {type(exc).__name__}: {exc}",
                "ap_runtime_s": round(time.time() - t0, 2),
                "ap_n_fixed": fix_n}

    q = _value(galaxy["q"])
    # AstroPhot and GALFIT share the position-angle convention (from +y);
    # AstroPhot reports it in radians.
    pa = vc.wrap_pa(math.degrees(_value(galaxy["PA"])))
    re_arcsec = _value(galaxy["Re"])

    out.update(
        {
            "ap_q": q,
            "ap_q_err": _uncertainty(galaxy["q"]),
            "ap_pa_deg": pa,
            "ap_re_arcsec": re_arcsec,
            "ap_re_px": re_arcsec / ps if ps > 0 else float("nan"),
            "ap_n": _value(galaxy["n"]),
            "ap_n_fixed": fix_n,
            "ap_sky_fixed_adu": sky,
            "ap_chi2nu": chi2nu,
            "ap_message": str(getattr(result, "message", "")),
            "ap_runtime_s": round(time.time() - t0, 2),
            "status": "ok",
        }
    )
    out["dq_astrophot"] = out["ap_q"] - host.q
    out["dpa_astrophot_deg"] = vc.wrap_pa(pa - vc.wrap_pa(host.pa))
    out["dre_astrophot_frac"] = (out["ap_re_px"] / host.re - 1.0
                                 if host.re > 0 else float("nan"))
    out["dn_astrophot"] = out["ap_n"] - host.n
    # Scale the disagreement by the statistical error: two codes differing by
    # 0.01 in q is fine unless GALFIT claims 0.001.
    out["dq_astrophot_over_q_err"] = (
        abs(out["dq_astrophot"]) / host.q_err
        if math.isfinite(host.q_err) and host.q_err > 0 else float("nan")
    )
    return out
