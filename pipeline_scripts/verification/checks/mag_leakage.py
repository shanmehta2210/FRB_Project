"""Check 7 — magnitude leakage.

A model can be wrong in a way chi2 tolerates but total flux does not. Sky
over-subtraction eats the outer disk; a runaway Sersic index inflates it. Both
land here as ``dmag = mag_GALFIT - mag_ref``.

Per-host this collects dmag and the quantities it is tested against; the trends
themselves (slope and significance against Re, n, sky offset, chi2nu) are fitted
across the cohort in ``aggregate.py``.
"""

from __future__ import annotations

import math

import numpy as np
import pandas as pd

import vercommon as vc

NAME = "mag"

# constraints.txt in production: n 0.5 to 6.0, re 1.5 to <re_max>
N_BOUNDS = (0.5, 6.0)
RE_FLOOR = 1.5
DMAG_FLAG = 0.5


def _constraint_bounds(path: str) -> dict:
    """Per-parameter bounds for component 1, as GALFIT was given them.

    Production writes ``<comp> <param> <lo> to <hi>``, e.g. ``1 re 1.5 to 100.0``.
    """
    bounds: dict[str, tuple[float, float]] = {}
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            for line in f:
                toks = line.split()
                if len(toks) >= 5 and toks[0] == "1" and toks[3].lower() == "to":
                    try:
                        bounds[toks[1].lower()] = (float(toks[2]), float(toks[4]))
                    except ValueError:
                        continue
    except OSError:
        pass
    return bounds


def _at_bound(value: float, bound: tuple[float, float] | None, frac: float = 0.05) -> bool:
    """True when a fitted value is pinned at a constraint limit.

    Measured relative to the *bound*, not to the width of the allowed range:
    with ``n`` free over 0.5-6.0, a range-relative test calls n = 0.62 pinned,
    which it plainly is not.
    """
    if bound is None or not math.isfinite(value):
        return False
    lo, hi = bound
    if math.isfinite(lo) and value <= lo * (1.0 + frac) + 1e-9:
        return True
    return math.isfinite(hi) and value >= hi * (1.0 - frac)


def run(host: vc.HostData, outdir: str) -> dict:
    row = vc.cohort("all64")
    row = row[row["frb"].astype(str) == host.frb]
    ref_mag = float("nan")
    out: dict = {}
    if len(row):
        r = row.iloc[0]
        ref_mag = pd.to_numeric(r.get("ref_mag"), errors="coerce")
        out.update(
            {
                "ref_survey": r.get("ref_survey"),
                "ref_mag": float(ref_mag) if pd.notna(ref_mag) else float("nan"),
                "ref_mag_err": float(pd.to_numeric(r.get("ref_mag_err"), errors="coerce")),
                "ref_sep_arcsec": float(pd.to_numeric(r.get("ref_sep_arcsec"), errors="coerce")),
                "zp_ok": bool(r.get("zp_ok")) if pd.notna(r.get("zp_ok")) else None,
                "mag_final_source": r.get("mag_final_source"),
                "snr_win": float(pd.to_numeric(r.get("snr_win"), errors="coerce")),
                "in_53": bool(r.get("in_53")),
            }
        )

    out["mag_galfit"] = host.mag
    out["mag_galfit_err"] = host.mag_err
    dmag = (host.mag - float(ref_mag)) if pd.notna(ref_mag) else float("nan")
    out["dmag_ref"] = float(dmag)
    out["dmag_flag"] = bool(math.isfinite(dmag) and abs(dmag) > DMAG_FLAG)

    out["re_px"] = host.re
    out["re_arcsec"] = host.re_arcsec
    out["n_sersic"] = host.n
    out["q_galfit"] = host.q
    out["chi2nu_global"] = host.chi2nu_log

    bounds = _constraint_bounds(f"{host.dir}/constraints.txt")
    out["n_bounds"] = list(bounds.get("n", N_BOUNDS))
    out["re_bounds"] = list(bounds.get("re", (RE_FLOOR, float("nan"))))
    nb = bounds.get("n", N_BOUNDS)
    out["n_at_bound"] = _at_bound(host.n, nb)
    out["n_at_ceiling"] = bool(math.isfinite(host.n) and host.n >= nb[1] * 0.95)
    out["n_at_floor"] = bool(math.isfinite(host.n) and host.n <= nb[0] * 1.05)
    out["re_at_bound"] = _at_bound(host.re, bounds.get("re"))
    out["mag_at_bound"] = _at_bound(host.mag, bounds.get("mag"))

    # Sky offset relative to the SExtractor seed, in units of the sky noise:
    # the quantity a dmag-versus-Re trend would be blaming.
    sky_ref = float(host.sky_audit.get("sky_ref_adu", float("nan")) or float("nan"))
    out["sky_fitted_adu"] = host.sky_level
    out["sky_ref_adu"] = sky_ref
    out["sky_offset_adu"] = (host.sky_level - sky_ref
                             if math.isfinite(sky_ref) and math.isfinite(host.sky_level)
                             else float("nan"))
    cal = vc.sigma_calibration_ratio(host)
    sky_mad = cal.get("sky_mad_adu", float("nan"))
    out["sky_offset_sigma"] = (out["sky_offset_adu"] / sky_mad
                               if math.isfinite(sky_mad) and sky_mad > 0
                               else float("nan"))

    # Flux the model carries beyond the stamp: a large fraction means the
    # magnitude is an extrapolation, not a measurement.
    try:
        model_in = float(np.nansum(vc.build_model(host, include_sky=False,
                                                  include_neighbours=False)))
        total = 10.0 ** ((host.magzpt - host.mag) / 2.5)
        out["flux_outside_stamp_frac"] = (1.0 - model_in / total) if total > 0 else float("nan")
    except Exception:
        out["flux_outside_stamp_frac"] = float("nan")
    return out
