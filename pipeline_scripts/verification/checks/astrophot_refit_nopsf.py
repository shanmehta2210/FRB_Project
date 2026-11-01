"""No-PSF control for Check 8 — identical to ``astrophot_refit`` except no PSF.

Same image, variance, mask, fixed sky (production GALFIT level), and the same
``n`` lock rule as the production AstroPhot refit; only ``psf_mode='none'`` and
no ``PSF_Image`` on the target. Results go to ``outputs/nopsf/astrophot/<FRB>.json``.

Run with the AstroPhot environment::

    .venv_astrophot\\Scripts\\python.exe pipeline_scripts/verification/checks/astrophot_refit_nopsf.py [FRB ...]
"""

from __future__ import annotations

import contextlib
import io
import math
import os
import sys
import time

import numpy as np

_VER_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _VER_DIR not in sys.path:
    sys.path.insert(0, _VER_DIR)

import vercommon as vc  # noqa: E402
from checks.astrophot_refit import _uncertainty, _value  # noqa: E402

NAME = "astrophot_nopsf"
OUT_DIR = os.path.join(vc.OUT_ROOT, "nopsf", "astrophot")
SEED_CSV = os.path.join(vc.REPO, "Publication", "materials", "host_properties", "seed_56.csv")


def run(host: vc.HostData, outdir: str, seed_from_galfit: bool = False) -> dict:
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
        variance[bad] = (np.nanmax(good) * 10.0) if good.size else 1.0

    mask = host.mask | ~np.isfinite(host.data)
    ps = host.plate_scale
    sky = host.sky_level if math.isfinite(host.sky_level) else 0.0
    fix_n = vc.host_n_held_fixed(os.path.join(host.dir, "galfit.feedme"))
    model_kw: dict = {
        "name": "host",
        "model_type": "sersic galaxy model",
        "psf_mode": "none",
    }
    if seed_from_galfit:
        model_kw["q"] = {"value": float(host.q)}
        model_kw["PA"] = {"value": math.radians(float(host.pa)) % math.pi}
        model_kw["Re"] = {"value": float(host.re) * ps}
        model_kw["n"] = {"value": float(host.n)}
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
            )
            galaxy = ap.models.AstroPhot_Model(target=target, **model_kw)
            galaxy.initialize()
            if fix_n is not None:
                try:
                    galaxy["n"].value = float(fix_n)
                    galaxy["n"].locked = True
                except Exception:
                    pass
            result = ap.fit.LM(galaxy, verbose=0).fit()
            chi2nu = float(result.res_loss())
    except Exception as exc:
        return {**out, "status": f"astrophot_failed: {type(exc).__name__}: {exc}",
                "ap_runtime_s": round(time.time() - t0, 2),
                "ap_n_fixed": fix_n}

    q = _value(galaxy["q"])
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
            "ap_psf_mode": "none",
            "status": "ok",
        }
    )
    out["dq_astrophot"] = out["ap_q"] - host.q
    out["dpa_astrophot_deg"] = vc.wrap_pa(pa - vc.wrap_pa(host.pa))
    out["dre_astrophot_frac"] = (out["ap_re_px"] / host.re - 1.0
                                 if host.re > 0 else float("nan"))
    out["dn_astrophot"] = out["ap_n"] - host.n
    return out


def main(frbs: list[str]) -> int:
    import pandas as pd

    if not frbs:
        frbs = pd.read_csv(SEED_CSV, dtype={"frb": str})["frb"].tolist()
    os.makedirs(OUT_DIR, exist_ok=True)
    for frb in frbs:
        path = os.path.join(OUT_DIR, f"{frb}.json")
        try:
            res = run(vc.load_host(frb), OUT_DIR)
        except Exception as exc:
            res = {"status": f"load_failed: {type(exc).__name__}: {exc}"}
        if not str(res.get("status", "")).startswith("ok"):
            first = res.get("status")
            try:
                res = run(vc.load_host(frb), OUT_DIR, seed_from_galfit=True)
            except Exception as exc:
                res = {"status": f"load_failed: {type(exc).__name__}: {exc}"}
            res["status"] = ("ok_retry_prod_seed" if res.get("status") == "ok"
                             else f"failed: {first}; retry {res.get('status')}")
        vc.write_json(path, {"_frb": frb, "_check": NAME, **res})
        print(f"[astrophot_nopsf] {frb}: {res.get('status')} q={res.get('ap_q')} n={res.get('ap_n')}",
              flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
