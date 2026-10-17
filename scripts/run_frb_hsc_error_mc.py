"""
Protocol B Monte Carlo of the confirmed-50 vs HSC tests.

Each draw: q ~ N(q, sigma_q_B), clip to (q0, 1], Hubble, then
i += N(0, 5 deg), clip to [0, 90], convert to cos i.
HSC stays the fixed point-estimate pool. Tests are pool tests.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from plot_frb_inflated_error_cdf import _hubble_inc_deg  # noqa: E402
from plot_frb_vs_hsc_confirmed50 import OUT_DIR, load_hsc_pool  # noqa: E402
from run_frb_hsc_confirmed50_tests import _chi2_hsc_bins, _edge_on_mask  # noqa: E402

N_DRAWS = 2000
SEED = 20240924
Q0 = 0.2
FLOOR_DEG = 5.0
P_2024 = 4.0 / 23.0


def _pct(a: np.ndarray) -> tuple[float, float, float]:
    lo, mid, hi = np.percentile(a, [16, 50, 84])
    return float(lo), float(mid), float(hi)


def main() -> None:
    hosts = pd.read_csv(OUT_DIR / "confirmed50_q.csv")
    q = pd.to_numeric(hosts["q"], errors="coerce").to_numpy(dtype=float)
    sig = pd.to_numeric(hosts["sigma_q_B"], errors="coerce").to_numpy(dtype=float)
    sig = np.maximum(np.nan_to_num(sig, nan=0.0), 0.0)
    if len(q) < 1 or not np.isfinite(q).all():
        raise RuntimeError(f"expected finite q, got {np.isfinite(q).sum()} of {len(q)}")

    hsc, n_hsc = load_hsc_pool(mag_limit=22.0, q0=Q0)
    edges = np.quantile(hsc, np.linspace(0.0, 1.0, 5))
    edges[0] = -np.inf
    edges[-1] = np.inf
    p_hsc = float(_edge_on_mask(hsc, 60.0).mean())
    p_scaled = (P_2024 / 0.5) * p_hsc

    rng = np.random.default_rng(SEED)
    ad_p = np.empty(N_DRAWS)
    ks_p = np.empty(N_DRAWS)
    mw_p = np.empty(N_DRAWS)
    chi_p = np.empty(N_DRAWS)
    k60 = np.empty(N_DRAWS, dtype=int)
    med = np.empty(N_DRAWS)
    mean = np.empty(N_DRAWS)
    counts = np.empty((N_DRAWS, 4), dtype=int)

    for i in range(N_DRAWS):
        qq = np.clip(rng.normal(q, sig), Q0 + 1e-6, 1.0)
        inc = _hubble_inc_deg(qq, q0=Q0)
        inc = np.clip(inc + rng.normal(0.0, FLOOR_DEG, size=len(q)), 0.0, 90.0)
        cosi = np.cos(np.radians(inc))
        ad_p[i] = float(stats.anderson_ksamp([hsc, cosi]).pvalue)
        ks_p[i] = float(stats.ks_2samp(cosi, hsc, alternative="less", method="asymp").pvalue)
        mw_p[i] = float(stats.mannwhitneyu(cosi, hsc, alternative="greater").pvalue)
        chi = _chi2_hsc_bins(cosi, hsc, 4)
        chi_p[i] = chi["p"]
        counts[i] = np.histogram(cosi, bins=edges)[0]
        k60[i] = int(np.sum(inc > 60.0))
        med[i] = float(np.median(cosi))
        mean[i] = float(np.mean(cosi))
        if (i + 1) % 500 == 0:
            print(f"draw {i + 1}", flush=True)

    n = len(q)
    p_hsc_less = np.array([stats.binom.cdf(k, n, p_hsc) for k in k60])
    p_2024_ge = np.array([stats.binom.sf(k - 1, n, P_2024) for k in k60])
    p_scaled_ge = np.array([stats.binom.sf(k - 1, n, p_scaled) for k in k60])

    def show(name: str, a: np.ndarray) -> None:
        lo, mid, hi = _pct(a)
        print(f"{name}: 16% {lo:.4g}  median {mid:.4g}  84% {hi:.4g}")

    print(f"N_FRB {n}  N_HSC {n_hsc}  draws {N_DRAWS}  median sigma_q_B {np.median(sig):.4f}")
    print(f"p_HSC(i>60) {p_hsc:.6f}  p_2024 {P_2024:.6f}  p_scaled {p_scaled:.6f}")
    show("median cos i", med)
    show("mean cos i", mean)
    show("AD p", ad_p)
    show("KS less p", ks_p)
    show("MW greater p", mw_p)
    show("chi2 p", chi_p)
    show("k(i>60)", k60.astype(float))
    show("binom less vs HSC", p_hsc_less)
    show("P(K>=k | 4/23)", p_2024_ge)
    show("P(K>=k | scaled)", p_scaled_ge)
    print(f"fraction AD p < 0.05: {np.mean(ad_p < 0.05):.3f}")
    print(f"fraction KS p < 0.05: {np.mean(ks_p < 0.05):.3f}")
    print(f"fraction MW p < 0.05: {np.mean(mw_p < 0.05):.3f}")
    print(f"fraction chi2 p < 0.05: {np.mean(chi_p < 0.05):.3f}")
    print(f"fraction k<=15: {np.mean(k60 <= 15):.3f}")
    print(f"fraction k>=15: {np.mean(k60 >= 15):.3f}")
    print("quartile counts median", np.median(counts, axis=0).tolist())
    print("quartile counts 16%", np.percentile(counts, 16, axis=0).tolist())
    print("quartile counts 84%", np.percentile(counts, 84, axis=0).tolist())
    # point-estimate reference
    cosi0 = pd.to_numeric(hosts["cosi_point"], errors="coerce").to_numpy(dtype=float)
    print("point median cos i", float(np.median(cosi0)))
    print("point k60", int(_edge_on_mask(cosi0, 60.0).sum()))


if __name__ == "__main__":
    main()
