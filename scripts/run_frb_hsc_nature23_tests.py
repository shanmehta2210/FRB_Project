"""Point and Protocol B tests for the 23 hosts in Bhardwaj et al. 2024."""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from plot_frb_inflated_error_cdf import _hubble_inc_deg  # noqa: E402
from plot_frb_vs_hsc_confirmed50 import OUT_DIR, load_hsc_pool  # noqa: E402
from run_frb_hsc_confirmed50_tests import (  # noqa: E402
    _chi2_hsc_bins,
    _edge_on_mask,
    directional_tests,
    mc_directional_p,
    observed_tests,
)

NATURE = _SCRIPTS.parent / "Archive" / "original paper" / "frb_sample.txt"
N_MC = 10_000
N_ERR = 2_000
SEED = 20240924
Q0 = 0.2
FLOOR = 5.0
P_2024 = 4.0 / 23.0


def _pct(a: np.ndarray) -> list[float]:
    return [float(x) for x in np.percentile(a, [16, 50, 84])]


def main() -> None:
    names = []
    paper_i = []
    for line in NATURE.read_text(encoding="utf-8").splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        names.append(parts[1].strip())
        paper_i.append(int(parts[2].replace('"', "").split(",")[0]))
    paper = pd.DataFrame({"frb": names, "i_2024": paper_i})

    hosts = pd.read_csv(OUT_DIR / "confirmed50_q.csv")
    merged = paper.merge(hosts, on="frb", how="left")
    missing = merged.loc[~merged["q"].notna(), "frb"].tolist()
    sub = merged.loc[merged["q"].notna()].copy()
    frb = pd.to_numeric(sub["cosi_point"], errors="coerce").to_numpy(dtype=float)
    q = pd.to_numeric(sub["q"], errors="coerce").to_numpy(dtype=float)
    sig = np.maximum(np.nan_to_num(pd.to_numeric(sub["sigma_q_B"], errors="coerce").to_numpy(dtype=float), nan=0.0), 0.0)
    hsc, n_hsc = load_hsc_pool(mag_limit=22.0, q0=Q0)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        obs = observed_tests(frb, hsc)
        direc = directional_tests(frb, hsc)
        mc = mc_directional_p(frb, hsc, n_draws=N_MC, seed=SEED)

    # AD matched-N, same seed family as the directional MC but separate stream.
    rng = np.random.default_rng(SEED + 7)
    ad_ge = 0
    ad_obs = obs["ad_stat"]
    for _ in range(N_MC):
        draw = rng.choice(hsc, size=len(frb), replace=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ad = float(stats.anderson_ksamp([hsc, draw]).statistic)
        if ad >= ad_obs:
            ad_ge += 1
    ad_matched = (1.0 + ad_ge) / (1.0 + N_MC)

    p_hsc = float(_edge_on_mask(hsc, 60.0).mean())
    k60 = int(_edge_on_mask(frb, 60.0).sum())
    p_scaled = (P_2024 / 0.5) * p_hsc
    binom_2024 = float(stats.binom.sf(k60 - 1, len(frb), P_2024))
    binom_scaled = float(stats.binom.sf(k60 - 1, len(frb), p_scaled))
    edges = np.quantile(hsc, [0.25, 0.5, 0.75])

    # Protocol B error draws.
    rng = np.random.default_rng(SEED)
    ad_p = np.empty(N_ERR)
    ks_p = np.empty(N_ERR)
    mw_p = np.empty(N_ERR)
    chi_p = np.empty(N_ERR)
    ks_d = np.empty(N_ERR)
    k_err = np.empty(N_ERR, dtype=int)
    med = np.empty(N_ERR)
    counts = np.empty((N_ERR, 4), dtype=int)
    qedges = np.quantile(hsc, np.linspace(0.0, 1.0, 5))
    qedges[0], qedges[-1] = -np.inf, np.inf
    for i in range(N_ERR):
        qq = np.clip(rng.normal(q, sig), Q0 + 1e-6, 1.0)
        inc = np.clip(_hubble_inc_deg(qq, q0=Q0) + rng.normal(0.0, FLOOR, size=len(q)), 0.0, 90.0)
        cosi = np.cos(np.radians(inc))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ad_p[i] = float(stats.anderson_ksamp([hsc, cosi]).pvalue)
        ks = stats.ks_2samp(cosi, hsc, alternative="less", method="asymp")
        ks_p[i] = float(ks.pvalue)
        ks_d[i] = float(ks.statistic)
        mw_p[i] = float(stats.mannwhitneyu(cosi, hsc, alternative="greater").pvalue)
        chi_p[i] = _chi2_hsc_bins(cosi, hsc, 4)["p"]
        counts[i] = np.histogram(cosi, bins=qedges)[0]
        k_err[i] = int(np.sum(inc > 60.0))
        med[i] = float(np.median(cosi))

    rows = []
    for _, r in sub.iterrows():
        cosi = float(r["cosi_point"])
        rows.append({
            "frb": r["frb"],
            "q": round(float(r["q"]), 3),
            "cosi": round(cosi, 3),
            "i_deg": round(float(np.degrees(np.arccos(np.clip(cosi, 0, 1)))), 1),
            "i_2024": int(r["i_2024"]),
            "sigma_q_B": round(float(r["sigma_q_B"]), 4),
        })

    out = {
        "n_nature": len(paper),
        "n_in_confirmed50": int(len(sub)),
        "missing": missing,
        "n_hsc": int(n_hsc),
        "median_cosi": float(np.median(frb)),
        "mean_cosi": float(np.mean(frb)),
        "hsc_median_cosi": float(np.median(hsc)),
        "hsc_mean_cosi": float(np.mean(hsc)),
        "ad_stat": obs["ad_stat"],
        "ad_p_pool": obs["ad_p"],
        "ad_p_matched": ad_matched,
        "ks_d": direc["ks_less_d"],
        "ks_p_pool": direc["ks_less_p"],
        "ks_p_matched": mc["ks_less_d"],
        "mw_u": direc["mwu_greater_u"],
        "mw_p_pool": direc["mwu_greater_p"],
        "mw_p_matched": mc["mwu_greater_u"],
        "mw_rbc": obs["mwu_rbc"],
        "chi2_stat": direc["chi2_stat"],
        "chi2_p": direc["chi2_p"],
        "chi2_p_matched": mc["chi2_stat"],
        "chi2_counts": direc["chi2_counts"],
        "chi2_expected": direc["chi2_expected"],
        "k60": k60,
        "p_hsc_60": p_hsc,
        "binom_hsc_p": direc["binom"][60.0]["p"],
        "binom_hsc_matched": mc["binom60"],
        "binom_2024_p": binom_2024,
        "binom_scaled_p": binom_scaled,
        "p_scaled": p_scaled,
        "hsc_quartile_cos": [float(x) for x in edges],
        "hosts": rows,
        "err": {
            "n": N_ERR,
            "median_cosi": _pct(med),
            "ad_p": _pct(ad_p),
            "ks_p": _pct(ks_p),
            "mw_p": _pct(mw_p),
            "chi_p": _pct(chi_p),
            "k60": _pct(k_err.astype(float)),
            "frac_ad_lt_05": float(np.mean(ad_p < 0.05)),
            "frac_ks_lt_05": float(np.mean(ks_p < 0.05)),
            "frac_mw_lt_05": float(np.mean(mw_p < 0.05)),
            "counts_median": np.median(counts, axis=0).tolist(),
            "p_2024": _pct(np.array([stats.binom.sf(k - 1, len(frb), P_2024) for k in k_err])),
            "p_hsc_less": _pct(np.array([stats.binom.cdf(k, len(frb), p_hsc) for k in k_err])),
        },
    }
    dest = OUT_DIR / "nature23_tests.json"
    dest.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps({k: out[k] for k in out if k != "hosts"}, indent=2))
    print("missing", missing)
    print("wrote", dest)


if __name__ == "__main__":
    main()
