"""Appendix robustness checks for the inclination analysis (51 hosts vs HSC).

A. Control sample: cos i histogram vs uniform; Sersic-n window; magnitude cuts;
   magnitude-matched draws; colour / quiescent cut.
B. Disk classification: Gini-M20 diagram with the Lotz et al. (2008) boundaries,
   reliability flags, quiescent hosts marked.
C1. q0 = 0.10, 0.15, 0.20, 0.25 applied to both samples.
C2. Measurement-uncertainty Monte Carlo (reflected, control scattered alike).

Every test block reports the same matched-N quantities (p_m = fraction of
N_DRAWS random control subsamples of the FRB sample size whose statistic
against the full pool is at least as extreme as the FRB statistic):
two-sided AD, two-sided MW, one-sided MW (edge-on deficit: FRB cos i larger),
and the quartile chi^2 (bins = quartiles of the pool in use).

Outputs: Publication/materials/appendix/robustness/
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

_SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS))
import plot_paper_inclination_figures as P  # noqa: E402
from null_catalog_utils import REPO_ROOT, hubble_cosi_from_ba  # noqa: E402
from plot_frb_inflated_error_cdf import _hubble_inc_deg  # noqa: E402

warnings.filterwarnings("ignore")

OUT = REPO_ROOT / "Publication" / "materials" / "appendix" / "robustness"
STATMORPH_FULL = REPO_ROOT / "pipeline_scripts" / "verification" / "outputs" / "tables" / "statmorph_51_full.csv"
HSC_PARENT = P.HSC_PARENT
Q0 = P.Q0
MAIN_SEED = 42        # run_frb_hsc_confirmed50_tests.SEED: the main-text p_m run (AD stream = SEED, MW / chi^2 stream = SEED + 1)
FLOOR_DEG = 5.0
SN_MIN = 2.5          # statmorph recommended minimum <S/N> per pixel
RES_MIN = 1.0         # R_e / FWHM below which the host is marginally resolved
QUIESCENT = ["20210807D", "20220509G", "20221012A", "20190523A"]
MAG_BINS = [-np.inf, 18.0, 19.0, 20.0, 21.0, 22.0]      # for the bookkeeping table
MAG_BINS_DRAW = [-np.inf, 19.0, 20.0, 21.0, 22.0]        # bright end merged (HSC has 1 galaxy at r <= 18)


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------
def cosi_from_q(q: np.ndarray, q0: float) -> np.ndarray:
    return np.array([hubble_cosi_from_ba(float(v), q0=q0) for v in np.asarray(q, dtype=float)])


def i_deg(cosi: np.ndarray) -> np.ndarray:
    return np.degrees(np.arccos(np.clip(cosi, -1, 1)))


def reflect_deg(inc: np.ndarray) -> np.ndarray:
    """Reflect inclinations back into [0, 90]."""
    inc = np.abs(inc)
    inc = np.where(inc > 90.0, 180.0 - inc, inc)
    return np.abs(inc)


def chi2_stat(sample: np.ndarray, edges: np.ndarray) -> tuple[float, np.ndarray]:
    counts = np.histogram(sample, bins=edges)[0].astype(float)
    exp = len(sample) / (len(edges) - 1)
    return float(((counts - exp) ** 2 / exp).sum()), counts


def observed_stats(frb: np.ndarray, pool: np.ndarray, edges: np.ndarray) -> dict:
    n, m = len(frb), len(pool)
    ad = stats.anderson_ksamp([pool, frb])
    u = float(stats.mannwhitneyu(frb, pool, alternative="greater").statistic)
    chi, counts = chi2_stat(frb, edges)
    return {
        "ad_stat": float(ad.statistic), "ad_p_pool": float(ad.pvalue),
        "mw_u": u, "mw_u_dev": abs(u - n * m / 2.0),
        "mw2_p_pool": float(stats.mannwhitneyu(frb, pool, alternative="two-sided").pvalue),
        "mw1_p_pool": float(stats.mannwhitneyu(frb, pool, alternative="greater").pvalue),
        "chi2_stat": chi, "chi2_counts": counts,
        "chi2_p_pool": float(stats.chi2.sf(chi, len(edges) - 2)),
    }


def null_draws(pool: np.ndarray, n: int, edges: np.ndarray, rng, n_draws: int, draw_fn=None) -> dict:
    """Statistics of n_draws subsamples of size n (draw_fn(rng) -> indices) vs the full pool.

    `rng` is either one Generator (shared stream, used by the error MC) or a pair (rng_ad, rng_dir):
    the AD draws and the MW / chi^2 draws then come from the same two streams as the main-text run
    (`run_frb_hsc_confirmed50_tests.py`: `mc_tail_p` seed SEED for AD, `mc_directional_p` SEED + 1 for MW / chi^2), so the
    adopted row of every table reproduces the main-text p_m values exactly.
    """
    m = len(pool)
    rng_ad, rng_dir = rng if isinstance(rng, tuple) else (rng, rng)
    ad = np.empty(n_draws); u = np.empty(n_draws); chi = np.empty(n_draws)
    for k in range(n_draws):
        idx = draw_fn(rng_ad) if draw_fn is not None else rng_ad.choice(m, n, replace=False)
        ad[k] = stats.anderson_ksamp([pool, pool[idx]]).statistic
        if rng_dir is not rng_ad:
            idx = draw_fn(rng_dir) if draw_fn is not None else rng_dir.choice(m, n, replace=False)
        d = pool[idx]
        u[k] = stats.mannwhitneyu(d, pool, alternative="greater").statistic
        chi[k] = chi2_stat(d, edges)[0]
    return {"ad": ad, "u": u, "u_dev": np.abs(u - n * m / 2.0), "chi2": chi}


def pm_from_null(obs: dict, null: dict) -> dict:
    # (1 + k) / (1 + N) as in the main-text run.
    n = len(null["ad"])
    p = lambda a, s: float((1.0 + np.sum(a >= s)) / (1.0 + n))
    return {
        "ad_pm": p(null["ad"], obs["ad_stat"]),
        "mw2_pm": p(null["u_dev"], obs["mw_u_dev"]),
        "mw1_pm": p(null["u"], obs["mw_u"]),
        "chi2_pm": p(null["chi2"], obs["chi2_stat"]),
    }


def matched_tests(frb: np.ndarray, pool: np.ndarray, rng, n_draws: int, draw_fn=None) -> dict:
    edges = np.quantile(pool, [0, 0.25, 0.5, 0.75, 1.0]); edges[0] = -np.inf; edges[-1] = np.inf
    obs = observed_stats(frb, pool, edges)
    # Fresh, fixed streams per variant (same seeds as the main-text run); `rng` is ignored here.
    streams = (np.random.default_rng(MAIN_SEED), np.random.default_rng(MAIN_SEED + 1))
    null = null_draws(pool, len(frb), edges, streams, n_draws, draw_fn)
    pm = pm_from_null(obs, null)
    out = {
        "N_frb": int(len(frb)), "N_pool": int(len(pool)),
        "median_cosi_frb": float(np.median(frb)), "median_cosi_pool": float(np.median(pool)),
        "median_i_frb": float(i_deg(np.median(frb))), "median_i_pool": float(i_deg(np.median(pool))),
        "frac_i_gt60_frb": float(np.mean(i_deg(frb) > 60)), "frac_i_gt60_pool": float(np.mean(i_deg(pool) > 60)),
        "quartile_counts": ",".join(str(int(c)) for c in obs["chi2_counts"][::-1]),  # face-on -> edge-on
        "ad_p_pool": obs["ad_p_pool"], "ad_pm": pm["ad_pm"],
        "mw2_p_pool": obs["mw2_p_pool"], "mw2_pm": pm["mw2_pm"],
        "mw1_p_pool": obs["mw1_p_pool"], "mw1_pm": pm["mw1_pm"],
        "chi2_stat": obs["chi2_stat"], "chi2_p_pool": obs["chi2_p_pool"], "chi2_pm": pm["chi2_pm"],
    }
    return out


def load_parent() -> pd.DataFrame:
    cols = ["rmag", "gmag", "imag", "ba", "n_sersic", "goodfits_flag", "quiescent_flag"]
    df = pd.read_csv(HSC_PARENT, usecols=cols)
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    ok = (df["goodfits_flag"] == 1) & np.isfinite(df["rmag"]) & np.isfinite(df["ba"]) & (df["ba"] > 0) & (df["ba"] <= 1) & np.isfinite(df["n_sersic"])
    return df.loc[ok].reset_index(drop=True)


def select(parent: pd.DataFrame, *, n_lo=0.4, n_hi=1.5, mag=22.0, q0=Q0, extra=None) -> pd.DataFrame:
    ok = (parent["n_sersic"] > n_lo) & (parent["n_sersic"] < n_hi) & (parent["rmag"] <= mag) & (parent["ba"] > q0)
    if extra is not None:
        ok &= extra
    return parent.loc[ok].reset_index(drop=True)


def f2(x): return f"{x:.2f}"
def f3(x): return f"{x:.3f}"


ROW_HDR = "| sample | $N_{\\rm FRB}$ | $N_{\\rm pool}$ | median $\\cos i$ FRB / pool | quartile counts | AD $p_{\\rm m}$ | MW 2-sided $p_{\\rm m}$ | MW 1-sided $p_{\\rm m}$ | $\\chi^2$ $p_{\\rm m}$ |"
ROW_SEP = "|---|---|---|---|---|---|---|---|---|"


def md_row(label: str, r: dict) -> str:
    return (f"| {label} | {r['N_frb']} | {r['N_pool']:,} | {f3(r['median_cosi_frb'])} / {f3(r['median_cosi_pool'])} | "
            f"{r['quartile_counts']} | {f2(r['ad_pm'])} | {f2(r['mw2_pm'])} | {f2(r['mw1_pm'])} | {f2(r['chi2_pm'])} |")


# ----------------------------------------------------------------------------
# A1. control histogram
# ----------------------------------------------------------------------------
def fig_control_hist(pool_cosi: np.ndarray) -> dict:
    inc = i_deg(pool_cosi)
    med_i = float(np.median(inc)); f60 = float(np.mean(inc > 60))
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6))
    ax = axes[0]
    ax.hist(pool_cosi, bins=40, range=(0, 1), density=True, color="0.6", edgecolor="white", label=f"HSC disks, $N={len(pool_cosi):,}$")
    ax.axhline(1.0, color="k", ls="--", lw=1.2, label="uniform in $\\cos i$ (random orientation)")
    ax.axvline(np.median(pool_cosi), color="C3", lw=1.2, label=f"median $\\cos i={np.median(pool_cosi):.3f}$ ($i={med_i:.1f}^\\circ$)")
    ax.set_xlabel("$\\cos i$"); ax.set_ylabel("density"); ax.set_xlim(0, 1); ax.legend(fontsize=7.5, loc="upper left")
    ax = axes[1]
    ax.hist(inc, bins=36, range=(0, 90), density=True, color="0.6", edgecolor="white")
    th = np.linspace(0, 90, 400)
    ax.plot(th, np.sin(np.radians(th)) * np.pi / 180, "k--", lw=1.2, label="random orientation ($\\propto \\sin i$)")
    ax.axvline(60, color="C0", lw=1.0, ls=":", label=f"$i>60^\\circ$: {100*f60:.1f}% (uniform: 50%)")
    ax.axvline(med_i, color="C3", lw=1.2, label=f"median $i={med_i:.1f}^\\circ$ (uniform: $60^\\circ$)")
    ax.set_xlabel("$i$ [deg]"); ax.set_xlim(0, 90); ax.legend(fontsize=7.5, loc="upper left")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"draft_figA_control_hist.{ext}", dpi=200)
    plt.close(fig)
    ks = stats.kstest(pool_cosi, "uniform")
    return {"N": int(len(pool_cosi)), "median_cosi": float(np.median(pool_cosi)), "median_i": med_i,
            "frac_i_gt60": f60, "mean_cosi": float(pool_cosi.mean()),
            "ks_vs_uniform_D": float(ks.statistic), "ks_vs_uniform_p": float(ks.pvalue)}


# ----------------------------------------------------------------------------
# B. Gini-M20
# ----------------------------------------------------------------------------
def fig_gini_m20(hosts: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    sm = pd.read_csv(STATMORPH_FULL, dtype={"frb": str})
    sm = sm.merge(hosts[["frb", "mag"]], on="frb", how="inner")
    assert len(sm) == 51, len(sm)
    sm["reliable"] = (sm["sm_flag"] == 0) & (sm["sm_sn_per_pixel"] >= SN_MIN) & (sm["re_over_fwhm"] >= RES_MIN)
    sm["why_unreliable"] = [
        "; ".join([s for s, c in (("flag", r.sm_flag != 0), (f"S/N<{SN_MIN}", r.sm_sn_per_pixel < SN_MIN),
                                  (f"Re/FWHM<{RES_MIN}", r.re_over_fwhm < RES_MIN)) if c])
        for r in sm.itertuples()]
    sm["quiescent_lit"] = sm["frb"].isin(QUIESCENT)

    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    m = np.linspace(-3.0, 0.0, 200)
    ax.plot(m, -0.14 * m + 0.33, "k-", lw=1.0)
    mm = m[m <= (0.33 - 0.80) / (0.14 + 0.14)]  # early/late line only below the merger line
    ax.plot(mm, 0.14 * mm + 0.80, "k-", lw=1.0)
    ax.text(-2.75, 0.73, "mergers", fontsize=9); ax.text(-2.75, 0.56, "E/S0/Sa", fontsize=9); ax.text(-1.05, 0.40, "Sb–Irr", fontsize=9)
    rel = sm[sm.reliable]; unr = sm[~sm.reliable]
    ax.scatter(rel.sm_m20, rel.sm_gini, s=34, c="C0", edgecolor="k", lw=0.5, label=f"reliable ({len(rel)})", zorder=3)
    ax.scatter(unr.sm_m20, unr.sm_gini, s=34, facecolor="none", edgecolor="C0", lw=1.0, label=f"flagged / marginal ({len(unr)})", zorder=3)
    qq = sm[sm.quiescent_lit]
    ax.scatter(qq.sm_m20, qq.sm_gini, s=120, marker="s", facecolor="none", edgecolor="C3", lw=1.2, label="quiescent in literature", zorder=2)
    lab = sm[(sm.lotz2008 != "late") | sm.quiescent_lit].sort_values("sm_gini")
    for k, r in enumerate(lab.itertuples()):
        dy = 4 if k % 2 == 0 else -9
        ax.annotate(r.frb, (r.sm_m20, r.sm_gini), xytext=(5, dy), textcoords="offset points", fontsize=6.5)
    ax.set_xlim(-3.0, 0.0); ax.set_ylim(0.3, 0.8)
    ax.set_xlabel("$M_{20}$"); ax.set_ylabel("Gini")
    ax.legend(fontsize=8, loc="lower left")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"draft_figB_gini_m20.{ext}", dpi=200)
    plt.close(fig)

    cols = ["frb", "mag", "re_over_fwhm", "sm_sn_per_pixel", "sm_flag", "sm_flag_sersic", "sm_gini", "sm_m20",
            "sm_F_bulge", "sm_C", "sm_A", "lotz2008", "reliable", "why_unreliable", "quiescent_lit", "lit_note"]
    cols = [c for c in cols if c in sm.columns]
    tab = sm[cols].sort_values("frb")
    tab.to_csv(OUT / "gini_m20_51.csv", index=False)
    summ = {
        "lotz_counts_all": tab.lotz2008.value_counts().to_dict(),
        "lotz_counts_reliable": tab[tab.reliable].lotz2008.value_counts().to_dict(),
        "n_reliable": int(tab.reliable.sum()), "n_unreliable": int((~tab.reliable).sum()),
        "n_flag": int((tab.sm_flag != 0).sum()), "n_lowsn": int((tab.sm_sn_per_pixel < SN_MIN).sum()),
        "n_marginal": int((tab.re_over_fwhm < RES_MIN).sum()),
        "not_late": tab[tab.lotz2008 != "late"][["frb", "lotz2008", "reliable", "sm_gini", "sm_m20"]].to_dict("records"),
        "quiescent": tab[tab.quiescent_lit][["frb", "lotz2008", "reliable", "sm_gini", "sm_m20", "sm_F_bulge", "why_unreliable"]].to_dict("records"),
    }
    return tab, summ


# ----------------------------------------------------------------------------
# C2. error Monte Carlo
# ----------------------------------------------------------------------------
def error_mc(hosts: pd.DataFrame, pool_q: np.ndarray, rng: np.random.Generator, n_mc: int, n_draws: int) -> dict:
    q = hosts["q"].to_numpy(dtype=float)
    sig = hosts["sigma_q_B"].to_numpy(dtype=float)
    pool_i0 = _hubble_inc_deg(pool_q, q0=Q0)
    n = len(q)

    def scatter_pool() -> np.ndarray:
        inc = reflect_deg(pool_i0 + rng.normal(0.0, FLOOR_DEG, size=len(pool_i0)))
        return np.cos(np.radians(inc))

    # Null distributions: draws from scattered pools (a fresh scattered pool every 500 draws).
    t0 = time.time()
    null = {"ad": [], "u": [], "u_dev": [], "chi2": []}
    block = 500
    for b in range(0, n_draws, block):
        pool_s = scatter_pool()
        edges = np.quantile(pool_s, [0, 0.25, 0.5, 0.75, 1.0]); edges[0] = -np.inf; edges[-1] = np.inf
        nd = null_draws(pool_s, n, edges, rng, min(block, n_draws - b))
        for k in null:
            null[k].append(nd[k])
    null = {k: np.concatenate(v) for k, v in null.items()}
    print(f"  [C2] null built in {time.time()-t0:.0f}s")

    keys = ["ad_pm", "mw2_pm", "mw1_pm", "chi2_pm", "ad_p_pool", "mw2_p_pool", "mw1_p_pool", "chi2_p_pool", "median_cosi", "n_gt60"]
    rec = {k: np.empty(n_mc) for k in keys}
    counts = np.empty((n_mc, 4))
    for k in range(n_mc):
        qq = rng.normal(q, sig)
        bad = (qq <= Q0) | (qq > 1.0)
        while bad.any():  # truncate, do not clip
            qq[bad] = rng.normal(q[bad], sig[bad]); bad = (qq <= Q0) | (qq > 1.0)
        inc = reflect_deg(_hubble_inc_deg(qq, q0=Q0) + rng.normal(0.0, FLOOR_DEG, size=n))
        frb = np.cos(np.radians(inc))
        pool_s = scatter_pool()
        edges = np.quantile(pool_s, [0, 0.25, 0.5, 0.75, 1.0]); edges[0] = -np.inf; edges[-1] = np.inf
        obs = observed_stats(frb, pool_s, edges)
        pm = pm_from_null(obs, null)
        for kk in ("ad_pm", "mw2_pm", "mw1_pm", "chi2_pm"):
            rec[kk][k] = pm[kk]
        for kk in ("ad_p_pool", "mw2_p_pool", "mw1_p_pool", "chi2_p_pool"):
            rec[kk][k] = obs[kk]
        rec["median_cosi"][k] = np.median(frb); rec["n_gt60"][k] = (inc > 60).sum()
        counts[k] = obs["chi2_counts"][::-1]
    pct = lambda a: [float(np.percentile(a, p)) for p in (16, 50, 84)]
    out = {"n_mc": n_mc, "n_draws_null": n_draws, "floor_deg": FLOOR_DEG, "N_frb": n, "N_pool": int(len(pool_q)),
           "median_sigma_q": float(np.median(sig)),
           "percentiles_16_50_84": {k: pct(v) for k, v in rec.items()},
           "frac_pm_lt_0.05": {k: float(np.mean(rec[k] < 0.05)) for k in ("ad_pm", "mw2_pm", "mw1_pm", "chi2_pm")},
           "frac_ppool_lt_0.05": {k: float(np.mean(rec[k] < 0.05)) for k in ("ad_p_pool", "mw2_p_pool", "mw1_p_pool", "chi2_p_pool")},
           "quartile_counts_median": [float(x) for x in np.median(counts, axis=0)]}
    out["max_frac_pm_lt_0.05"] = max(out["frac_pm_lt_0.05"].values())
    pd.DataFrame(rec).to_csv(OUT / "error_mc_realizations.csv", index=False)
    return out


# ----------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=10_000)
    ap.add_argument("--mc", type=int, default=2_000)
    ap.add_argument("--skip", default="", help="comma list of blocks to skip: A1,A2,A3,A4,A5,B,C1,C2")
    args = ap.parse_args()
    skip = set(s.strip() for s in args.skip.split(",") if s.strip())
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(P.SEED + 101)

    hosts = P.load_hosts(); assert len(hosts) == 51
    frb = hosts["cosi_point"].to_numpy(dtype=float)
    parent = load_parent()
    base = select(parent)
    base_cosi = cosi_from_q(base["ba"].to_numpy(), Q0)
    assert len(base_cosi) == 24_450, len(base_cosi)
    res: dict = {"n_draws": args.draws, "n_mc": args.mc}
    md: list[str] = []

    # ---- A1 ------------------------------------------------------------------
    if "A1" not in skip:
        res["A1_control"] = fig_control_hist(base_cosi)
        a = res["A1_control"]
        md += ["## A1. The control sample is not uniform in $\\cos i$", "",
               f"`figA_control_hist.pdf` (final version: `scripts/appendix_figures.py`): HSC disk pool ($N = {a['N']:,}$, $r \\le 22$, $0.4 < n < 1.5$, $q > 0.2$). "
               f"Median $\\cos i = {a['median_cosi']:.3f}$ ($i = {a['median_i']:.1f}^\\circ$); fraction with $i > 60^\\circ$: {100*a['frac_i_gt60']:.1f}% "
               f"(random orientation: $60^\\circ$ and 50%). KS against uniform in $\\cos i$: $D = {a['ks_vs_uniform_D']:.3f}$, $p = {a['ks_vs_uniform_p']:.1e}$. "
               "The pool is depleted at both ends: at $\\cos i \\to 0$ because the $q > q_0$ cut, the magnitude limit and dust remove the most inclined disks, and because the "
               "Hubble formula maps $q$ just above $q_0$ onto a narrow range of $\\cos i$; at $\\cos i \\to 1$ because measured axis ratios rarely reach $q = 1$ "
               "(intrinsic disk ellipticity, non-axisymmetric structure, PSF and noise). Tests must therefore be made against this pool, not against a uniform distribution, "
               "and the quartile bins of Figure 2 are defined on it.", ""]
        print("[A1]", a)

    # ---- A2 Sersic window ----------------------------------------------------
    if "A2" not in skip:
        rows = []
        for lo, hi, lab in [(0.4, 1.5, "$0.4 < n < 1.5$ (adopted)"), (0.4, 2.0, "$0.4 < n < 2.0$"),
                            (0.4, 2.5, "$0.4 < n < 2.5$"), (0.0, 2.5, "$n < 2.5$"), (0.0, 10.0, "all $n$ (no morphology cut)")]:
            sel = select(parent, n_lo=lo, n_hi=hi)
            r = {"label": lab, "n_lo": lo, "n_hi": hi, "median_n_pool": float(np.median(sel["n_sersic"]))}
            r.update(matched_tests(frb, cosi_from_q(sel["ba"].to_numpy(), Q0), rng, args.draws)); rows.append(r)
            print("[A2]", lab, {k: r[k] for k in ("N_pool", "ad_pm", "mw2_pm", "mw1_pm", "chi2_pm")})
        df = pd.DataFrame(rows); df.to_csv(OUT / "A2_sersic_window.csv", index=False); res["A2"] = rows
        md += ["## A2. Sérsic-index window of the control", "", "FRB sample fixed (51). Pool re-selected; quartile bins recomputed from each pool.", "",
               ROW_HDR, ROW_SEP] + [md_row(r["label"] + f" (median $n$ = {r['median_n_pool']:.2f})", r) for r in rows] + [""]

    # ---- A3 magnitude cuts ---------------------------------------------------
    if "A3" not in skip:
        rows = []
        for mag in (22.0, 21.0, 20.0):
            sel = select(parent, mag=mag)
            pool = cosi_from_q(sel["ba"].to_numpy(), Q0)
            r = {"label": f"control $r \\le {mag:.0f}$, all 51 hosts", "mag": mag, "frb_cut": False}
            r.update(matched_tests(frb, pool, rng, args.draws)); rows.append(r)
            print("[A3]", r["label"], {k: r[k] for k in ("N_pool", "ad_pm", "mw2_pm", "mw1_pm", "chi2_pm")})
            if mag < 22:
                f2m = hosts.loc[hosts["mag"] <= mag, "cosi_point"].to_numpy(dtype=float)
                r = {"label": f"both samples $r \\le {mag:.0f}$", "mag": mag, "frb_cut": True}
                r.update(matched_tests(f2m, pool, rng, args.draws)); rows.append(r)
                print("[A3]", r["label"], {k: r[k] for k in ("N_frb", "N_pool", "ad_pm", "mw2_pm", "mw1_pm", "chi2_pm")})
        df = pd.DataFrame(rows); df.to_csv(OUT / "A3_magnitude_cuts.csv", index=False); res["A3"] = rows
        md += ["## A3. Magnitude limit of the control", "", ROW_HDR, ROW_SEP] + [md_row(r["label"], r) for r in rows] + [""]

    # ---- A4 magnitude-matched draws -----------------------------------------
    if "A4" not in skip:
        fmag = hosts["mag"].to_numpy(dtype=float)
        pmag = base["rmag"].to_numpy(dtype=float)
        fbin = np.digitize(fmag, MAG_BINS[1:-1]); pbin = np.digitize(pmag, MAG_BINS[1:-1])
        nb = len(MAG_BINS) - 1
        bins_tab = []
        for b in range(nb):
            lo, hi = MAG_BINS[b], MAG_BINS[b + 1]
            bins_tab.append({"bin": f"{'' if np.isinf(lo) else f'{lo:.0f} <'} r <= {hi:.0f}", "n_frb": int((fbin == b).sum()),
                             "n_control": int((pbin == b).sum()),
                             "frb_hosts": ", ".join(hosts.loc[fbin == b, "frb"])})
        res["A4_bins"] = bins_tab
        for bt in bins_tab:
            print("[A4] bin", bt)
        # draws use the merged binning
        fbd = np.digitize(fmag, MAG_BINS_DRAW[1:-1]); pbd = np.digitize(pmag, MAG_BINS_DRAW[1:-1])
        nbd = len(MAG_BINS_DRAW) - 1
        idx_by_bin = [np.flatnonzero(pbd == b) for b in range(nbd)]
        draw_bins = [{"bin": f"{'' if np.isinf(MAG_BINS_DRAW[b]) else f'{MAG_BINS_DRAW[b]:.0f} <'} r <= {MAG_BINS_DRAW[b+1]:.0f}",
                      "n_frb": int((fbd == b).sum()), "n_control": int(len(idx_by_bin[b]))} for b in range(nbd)]
        res["A4_draw_bins"] = draw_bins

        def run_matched(mask_frb: np.ndarray, label: str) -> dict:
            nbb = [int((fbd[mask_frb] == b).sum()) for b in range(nbd)]
            def draw(rg):
                return np.concatenate([rg.choice(idx_by_bin[b], nbb[b], replace=False) for b in range(nbd) if nbb[b] > 0])
            r = {"label": label, "n_frb_per_bin": ",".join(map(str, nbb)),
                 "n_control_per_bin": ",".join(str(len(idx_by_bin[b])) for b in range(nbd))}
            r.update(matched_tests(frb[mask_frb], base_cosi, rng, args.draws, draw_fn=draw))
            return r
        rows = [run_matched(np.ones(51, bool), "all 51, bright end merged into one $r \\le 19$ bin"),
                run_matched(fmag > 19.0, "hosts with $r > 19$ only")]
        for r in rows:
            print("[A4]", r["label"], {k: r[k] for k in ("N_frb", "ad_pm", "mw2_pm", "mw1_pm", "chi2_pm")})
        pd.DataFrame(rows).to_csv(OUT / "A4_magnitude_matched.csv", index=False); res["A4"] = rows
        md += ["## A4. Magnitude-matched control draws", "",
               "Each null draw takes, from every $r$ bin, the same number of control galaxies as there are FRB hosts in that bin (without replacement); "
               "$p_{\\rm m}$ is the fraction of such magnitude-matched draws that differ from the full pool at least as much as the FRB hosts. "
               "The FRB statistic itself is unchanged (51 hosts vs the full pool); only the null changes.", "",
               "| $r$ bin | FRB hosts | control galaxies | hosts |", "|---|---|---|---|"]
        md += [f"| {bt['bin']} | {bt['n_frb']} | {bt['n_control']:,} | {bt['frb_hosts']} |" for bt in bins_tab]
        md += ["", "The HSC pool is essentially empty brighter than $r = 18$ (bright galaxies saturate / are excluded from the Kawinwanichakij et al. fits), so a bin-by-bin match "
               "of the ten $r \\le 18$ hosts is impossible. Two runs are made: the bright end merged into a single $r \\le 19$ bin "
               f"({draw_bins[0]['n_frb']} hosts drawn from {draw_bins[0]['n_control']} control galaxies, so those draws are strongly repeated), and hosts with $r > 19$ only.", "",
               ROW_HDR, ROW_SEP] + [md_row(r["label"], r) for r in rows] + [""]

    # ---- A5 colour / quiescent -----------------------------------------------
    if "A5" not in skip:
        rows = []
        gi = parent["gmag"] - parent["imag"]
        for lab, extra in [("star-forming only (`quiescent_flag = 0`)", parent["quiescent_flag"] == 0),
                           ("blue, $g - i < 1.2$", gi < 1.2), ("red, $g - i \\ge 1.2$", gi >= 1.2)]:
            sel = select(parent, extra=extra)
            r = {"label": lab}
            r.update(matched_tests(frb, cosi_from_q(sel["ba"].to_numpy(), Q0), rng, args.draws)); rows.append(r)
            print("[A5]", lab, {k: r[k] for k in ("N_pool", "ad_pm", "mw2_pm", "mw1_pm", "chi2_pm")})
        pd.DataFrame(rows).to_csv(OUT / "A5_colour_cut.csv", index=False); res["A5"] = rows
        md += ["## A5. Colour / star-formation cut on the control", "",
               "HSC $g - i$ is observed-frame (no $k$-correction); `quiescent_flag` is the Kawinwanichakij et al. (2021) UVJ-based flag. The FRB sample is unchanged.", "",
               ROW_HDR, ROW_SEP] + [md_row(r["label"], r) for r in rows] + [""]

    # ---- B Gini-M20 ----------------------------------------------------------
    if "B" not in skip:
        tab, summ = fig_gini_m20(hosts); res["B"] = summ
        print("[B]", json.dumps(summ, indent=1, default=str)[:1500])
        lc, lr = summ["lotz_counts_all"], summ["lotz_counts_reliable"]
        md += ["## B. Gini–$M_{20}$ (statmorph on the GALFIT cutouts, PSF and masks)", "",
               f"`figB_gini_m20.pdf` (final version: `scripts/appendix_figures.py`). Lotz et al. (2008) boundaries: mergers $G > -0.14\\,M_{{20}} + 0.33$; E/S0/Sa $G > 0.14\\,M_{{20}} + 0.80$ (below the merger line); Sb–Irr otherwise.", "",
               f"All 51: late {lc.get('late',0)}, early {lc.get('early',0)}, merger {lc.get('merger',0)}. "
               f"Reliable only (`flag = 0`, $\\langle S/N\\rangle_{{\\rm pix}} \\ge {SN_MIN}$, $R_e/{{\\rm FWHM}} \\ge {RES_MIN}$; $N = {summ['n_reliable']}$): "
               f"late {lr.get('late',0)}, early {lr.get('early',0)}, merger {lr.get('merger',0)}. "
               f"Marked unreliable: {summ['n_unreliable']} (flag $\\ne 0$: {summ['n_flag']}; low S/N: {summ['n_lowsn']}; marginally resolved: {summ['n_marginal']}).", "",
               "Hosts outside the late-type region:", "", "| FRB | class | reliable | $G$ | $M_{20}$ |", "|---|---|---|---|---|"]
        md += [f"| {r['frb']} | {r['lotz2008']} | {r['reliable']} | {r['sm_gini']:.3f} | {r['sm_m20']:.2f} |" for r in summ["not_late"]]
        md += ["", "Hosts described as quiescent in the literature:", "", "| FRB | class | reliable | $G$ | $M_{20}$ | $F(G,M_{20})$ | note |", "|---|---|---|---|---|---|---|"]
        md += [f"| {r['frb']} | {r['lotz2008']} | {r['reliable']} | {r['sm_gini']:.3f} | {r['sm_m20']:.2f} | {r['sm_F_bulge']:+.2f} | {r['why_unreliable']} |" for r in summ["quiescent"]]
        md += ["", "Per-host values: `gini_m20_51.csv`.", ""]

    # ---- C1 q0 ---------------------------------------------------------------
    if "C1" not in skip:
        qser = P.magnitude_q(hosts)
        par_q = parent.loc[(parent["n_sersic"] > 0.4) & (parent["n_sersic"] < 1.5) & (parent["rmag"] <= 22.0), "ba"].to_numpy(dtype=float)
        rows = []
        for q0 in (0.10, 0.15, 0.20, 0.25):
            fq = qser[qser > q0]
            added = qser[(qser > q0) & (qser <= Q0)]; removed = qser[(qser <= q0) & (qser > Q0)]
            r = {"q0": q0, "added_vs_0.2": ", ".join(f"{i} ($q={v:.2f}$)" for i, v in added.items()),
                 "removed_vs_0.2": ", ".join(f"{i} ($q={v:.2f}$)" for i, v in removed.items())}
            r.update(matched_tests(cosi_from_q(fq.to_numpy(), q0), cosi_from_q(par_q[par_q > q0], q0), rng, args.draws)); rows.append(r)
            print("[C1]", q0, {k: r[k] for k in ("N_frb", "N_pool", "ad_pm", "mw2_pm", "mw1_pm", "chi2_pm", "added_vs_0.2", "removed_vs_0.2")})
        pd.DataFrame(rows).to_csv(OUT / "C1_q0.csv", index=False); res["C1"] = rows
        md += ["## C1. Intrinsic thickness $q_0$", "",
               "$\\cos^2 i = (q^2 - q_0^2)/(1 - q_0^2)$ and the cut $q > q_0$ applied to both samples. Hosts with $q \\le 0.2$ that pass the magnitude and verification cuts re-enter at lower $q_0$.", "",
               "| $q_0$ | hosts added | hosts removed |", "|---|---|---|"]
        md += [f"| {r['q0']:.2f} | {r['added_vs_0.2'] or '—'} | {r['removed_vs_0.2'] or '—'} |" for r in rows]
        md += ["", ROW_HDR.replace("sample", "$q_0$"), ROW_SEP] + [md_row(f"{r['q0']:.2f}", r) for r in rows] + [""]

    # ---- C2 error MC ---------------------------------------------------------
    if "C2" not in skip:
        res["C2"] = error_mc(hosts, base["ba"].to_numpy(dtype=float), rng, args.mc, args.draws)
        c = res["C2"]; pc = c["percentiles_16_50_84"]
        print("[C2]", json.dumps(c, indent=1)[:2500])
        md += ["## C2. Measurement uncertainties (Monte Carlo)", "",
               f"{c['n_mc']} realizations. Each host: $q \\sim \\mathcal{{N}}(q, \\sigma_q)$ truncated to $q_0 < q \\le 1$ (redrawn, not clipped; median $\\sigma_q = {c['median_sigma_q']:.3f}$), "
               f"converted to $i$, plus $\\mathcal{{N}}(0, {FLOOR_DEG:.0f}^\\circ)$ reflected at $0^\\circ$ and $90^\\circ$. "
               f"The same $5^\\circ$ reflected scatter is applied to every control galaxy in every realization, so both samples are treated alike. "
               f"$p_{{\\rm m}}$ in each realization is read against a null of {c['n_draws_null']:,} draws of 51 from scattered control pools.", "",
               "| quantity | 16% | median | 84% | fraction of realizations with $p_{\\rm m} < 0.05$ |", "|---|---|---|---|---|"]
        for k, lab in (("ad_pm", "AD $p_{\\rm m}$"), ("mw2_pm", "MW two-sided $p_{\\rm m}$"), ("mw1_pm", "MW one-sided $p_{\\rm m}$"), ("chi2_pm", "quartile $\\chi^2$ $p_{\\rm m}$")):
            md.append(f"| {lab} | {pc[k][0]:.2f} | {pc[k][1]:.2f} | {pc[k][2]:.2f} | {100*c['frac_pm_lt_0.05'][k]:.1f}% |")
        md += [f"| median $\\cos i$ | {pc['median_cosi'][0]:.3f} | {pc['median_cosi'][1]:.3f} | {pc['median_cosi'][2]:.3f} | |",
               f"| $N(i > 60^\\circ)$ | {pc['n_gt60'][0]:.0f} | {pc['n_gt60'][1]:.0f} | {pc['n_gt60'][2]:.0f} | |",
               "", f"Largest fraction of realizations in which any test reaches $p_{{\\rm m}} < 0.05$: **{100*c['max_frac_pm_lt_0.05']:.1f}%**. "
               f"Median quartile counts (face-on to edge-on): {', '.join(f'{x:.0f}' for x in c['quartile_counts_median'])}.", ""]

    if skip and (OUT / "robustness_results.json").exists():  # keep previously computed blocks
        old = json.loads((OUT / "robustness_results.json").read_text(encoding="utf-8"))
        res = {**old, **res}
    (OUT / "robustness_results.json").write_text(json.dumps(res, indent=1, default=str), encoding="utf-8")
    head = ["# Appendix robustness checks (51 hosts vs HSC)", "",
            f"Source: `scripts/appendix_robustness.py` (draws per null: {args.draws:,}; error-MC realizations: {args.mc:,}; matched-N draws use the main-text streams, seed {MAIN_SEED} for AD and {MAIN_SEED + 1} for MW / $\\chi^2$, "
            f"so the adopted row equals the main-text run; error-MC seed {P.SEED + 101}). "
            "Hosts: `plots/plots_null/v2/frb_vs_hsc_confirmed50/confirmed50_q.csv` (51, winning-leg $q$). Control: Kawinwanichakij et al. (2021) HSC-SSP parent, `goodfits_flag = 1`, "
            "$r \\le 22$, $0.4 < n < 1.5$, $q > 0.2$ unless stated. Inclinations from $\\cos^2 i = (q^2 - q_0^2)/(1 - q_0^2)$, $q_0 = 0.2$.", "",
            "$p_{\\rm m}$ = fraction of random control subsamples of the FRB sample size whose statistic against the full pool is at least as extreme as the FRB value "
            "(two-sided AD statistic; two-sided MW $|U - nm/2|$; one-sided MW $U$ in the edge-on-deficit direction, i.e. FRB $\\cos i$ larger; quartile $\\chi^2$ with bins at the pool quartiles). "
            "Quartile counts are listed face-on to edge-on.", ""]
    (OUT / "ROBUSTNESS.md").write_text("\n".join(head + md), encoding="utf-8")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
