"""Correlations among redshift, stellar mass, SFR, DM_host and host inclination, within tertiles of one quantity.

Grid A (12 panels): for each binning quantity C in (z, M*, SFR, DM_host), the three pairs of the other three.
Grid B (16 panels): for each C, each of the four quantities against inclination.

Each panel uses every host with the values it needs. Per tertile of C: Spearman r_S, bootstrap 95% interval,
delete-1 and delete-3 jackknife ranges, and Benjamini-Hochberg q over the tertile correlations of the grid.
Per panel: the pooled within-tertile correlation (Fisher z, var = 1.06 / (n - 3)) with Cochran's Q test of
whether the tertiles agree, and Kendall's partial tau of the pair given C (p from shuffling within strata of C).
Alongside every Spearman test, a two-sample Anderson-Darling test of whether y has the same distribution in the
lower and upper halves of x (permutation p for the measured values, asymptotic p inside the Monte Carlo).

Measurement errors give the error bars and 2,000 Monte Carlo draws in which every value is redrawn and the
tertiles are re-formed before the correlations are recomputed:
  z        spectroscopic redshifts fixed; the photometric redshift of 20240229A with sigma_z = 0.03
  M*, SFR  published asymmetric errors as split normals in dex, at least 0.1 dex; missing SFR errors get
           the sample median
  DM_host  DM_cosmic from the Macquart et al. (2020) distribution and DM_MW with a 20% error
  i        q ~ N(q, sigma_q) redrawn outside (q0, 1], Hubble formula, plus N(0, 5 deg) reflected into [0, 90]

Output: Publication/materials/appendix/binned_correlations/
  RESULTS.md, binned_correlations.csv, grid.{pdf,png}, grid_inclination.{pdf,png},
  heatmap.{pdf,png}, heatmap_inclination.{pdf,png}, binned_by_<C>__<y>_vs_<x>.png
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dm_host_physics as D  # noqa: E402
import robust_correlations as RC  # noqa: E402
import robust_extras as E  # noqa: E402

OUT = D.ROOT / "Publication" / "materials" / "appendix" / "binned_correlations"
PROPS = D.ROOT / "Publication" / "materials" / "host_properties" / "host_properties_57.csv"
DM = D.DM
INC = "i_deg"
VARS = ["z", "logmstar", "log_sfr", DM]
NAMES = {"z": "z", "logmstar": "M*", "log_sfr": "SFR", DM: "DM_host", INC: "i"}
AXIS = {"z": "Redshift", "logmstar": r"log $M_\star$ (M$_\odot$)", "log_sfr": r"log SFR (M$_\odot$ yr$^{-1}$)",
        DM: r"DM$_{\rm host}$ (rest frame, pc cm$^{-3}$)", INC: r"Inclination $i$ (deg)"}
SHORT = {"z": "z", "logmstar": r"log $M_\star$", "log_sfr": "log SFR", DM: r"DM$_{\rm host}$", INC: r"$i$"}
COLORS = ["#3b6fb6", "#e0a526", "#c0392b"]
N_BOOT = 2_000
N_MC = D.N_MC
Q0 = 0.2
SYS_I_DEG = 5.0
LOG_ERR_FLOOR = 0.1
PHOTO_Z = {"20240229A": 0.03}


# ---------------------------------------------------------------- measurement errors
def load_hosts() -> pd.DataFrame:
    h = D.load().reset_index(drop=True)
    sfr_err = pd.read_csv(PROPS, dtype={"frb": str})[["frb", "sfr_ep", "sfr_em"]]
    return h.merge(sfr_err, on="frb", how="left")


def split_normal(mid, lo, hi, rng):
    u = rng.standard_normal((N_MC, len(mid)))
    return mid + u * np.where(u > 0, hi, lo)


def log_sfr_errors(h):
    hi = np.log10(1 + h.sfr_ep / h.sfr)
    lo = np.where(h.sfr_em < h.sfr, -np.log10(1 - h.sfr_em / h.sfr), hi)
    lo, hi = (np.maximum(np.asarray(e, float), LOG_ERR_FLOOR) for e in (lo, hi))
    med = np.nanmedian(np.r_[lo, hi])
    return np.where(np.isfinite(lo), lo, med), np.where(np.isfinite(hi), hi, med)


def inclination_draws(q, sig, rng):
    q, sig = np.broadcast_to(q, (N_MC, len(q))), np.broadcast_to(sig, (N_MC, len(q)))
    qq = rng.normal(q, sig)
    bad = (qq <= Q0) | (qq > 1)
    while bad.any():
        qq[bad] = rng.normal(q[bad], sig[bad])
        bad = (qq <= Q0) | (qq > 1)
    i = np.degrees(np.arccos(np.sqrt((qq ** 2 - Q0 ** 2) / (1 - Q0 ** 2)))) + rng.normal(0, SYS_I_DEG, qq.shape)
    i = np.abs(i)
    return np.where(i > 90, 180 - i, i)


def measurement_draws(h, rng) -> dict[str, np.ndarray]:
    """(N_MC, N_hosts) draws of every quantity."""
    z_err = h.frb.map(PHOTO_Z).fillna(0.0).to_numpy()
    m_lo = np.maximum(h.logmstar_em.abs().fillna(LOG_ERR_FLOOR), LOG_ERR_FLOOR).to_numpy()
    m_hi = np.maximum(h.logmstar_ep.abs().fillna(LOG_ERR_FLOOR), LOG_ERR_FLOOR).to_numpy()
    return {"z": np.clip(h.z.to_numpy() + rng.normal(0, 1, (N_MC, len(h))) * z_err, 1e-3, None),
            "logmstar": split_normal(h.logmstar.to_numpy(), m_lo, m_hi, rng),
            "log_sfr": split_normal(h.log_sfr.to_numpy(), *log_sfr_errors(h), rng),
            DM: D.dm_host_draws(h, rng),
            INC: inclination_draws(h.q.to_numpy(), h.sigma_q_B.fillna(0).to_numpy(), rng)}


def error_bars(h, draws) -> dict[str, np.ndarray]:
    """(2, N_hosts) distances from each value to the 16th and 84th percentiles of its draws."""
    out = {}
    for v, d in draws.items():
        p16, p84 = np.nanpercentile(d, [16, 84], axis=0)
        mid = h[v].to_numpy(float)
        out[v] = np.nan_to_num(np.clip([mid - p16, p84 - mid], 0, None))
    return out


# ---------------------------------------------------------------- statistics
def tertiles(v):
    return np.digitize(v, np.percentile(v, [100 / 3, 200 / 3]))


def fisher_pool(r, n):
    """Pooled correlation over bins, its p-value, and Cochran's Q p-value for equal correlations."""
    r, n = np.asarray(r, float), np.asarray(n, float)
    z, w = np.arctanh(np.clip(r, -0.999, 0.999)), (n - 3) / 1.06
    zbar = np.sum(w * z) / np.sum(w)
    p = 2 * stats.norm.sf(abs(zbar) * np.sqrt(np.sum(w)))
    q = np.sum(w * (z - zbar) ** 2)
    return float(np.tanh(zbar)), float(p), float(stats.chi2.sf(q, len(r) - 1))


def bin_stats(x, y, rng):
    r, p = stats.spearmanr(x, y)
    p_ad = RC.ad_split(x, y)
    boot = [stats.spearmanr(x[i], y[i])[0] for i in (rng.integers(0, len(x), len(x)) for _ in range(N_BOOT))]
    jack = [stats.spearmanr(np.delete(x, k), np.delete(y, k))[0] for k in range(len(x))]
    lo, hi = np.nanpercentile(boot, [2.5, 97.5])
    jk3_min, jk3_max, _, _, jk3_fsig = E.jackknife(x, y, 3)
    return dict(N=len(x), r=r, p=p, p_ad=p_ad, lo=lo, hi=hi, jk_min=min(jack), jk_max=max(jack),
                jk3_min=jk3_min, jk3_max=jk3_max, jk3_fsig=jk3_fsig)


def monte_carlo(X, Y, C, partial):
    r_bin, p_bin = np.full((N_MC, 3), np.nan), np.full((N_MC, 3), np.nan)
    ad_bin = np.full((N_MC, 3), np.nan)
    r_all, p_all, ad_all, ptau = np.empty(N_MC), np.empty(N_MC), np.empty(N_MC), np.full(N_MC, np.nan)
    for d in range(N_MC):
        t = tertiles(C[d])
        for k in range(3):
            m = t == k
            r_bin[d, k], p_bin[d, k] = stats.spearmanr(X[d, m], Y[d, m])
            ad_bin[d, k] = RC.ad_split(X[d, m], Y[d, m], exact=False)
        r_all[d], p_all[d] = stats.spearmanr(X[d], Y[d])
        ad_all[d] = RC.ad_split(X[d], Y[d], exact=False)
        if partial:
            ptau[d] = D.partial_tau(X[d], Y[d], C[d])
    pct = lambda a: np.nanpercentile(a, [50, 16, 84], axis=0)
    return dict(r_bin=pct(r_bin), fsig_bin=np.mean(p_bin < 0.05, axis=0), r_all=pct(r_all),
                fsig_all=float(np.mean(p_all < 0.05)), ad_fsig_bin=np.mean(ad_bin < 0.05, axis=0),
                ad_fsig_all=float(np.mean(ad_all < 0.05)), ptau=pct(ptau) if partial else np.full(3, np.nan))


def analyse(h, draws, c, xc, yc, rng) -> dict:
    cols = list(dict.fromkeys([c, xc, yc]))
    idx = h.index[h[cols].notna().all(axis=1)].to_numpy()
    x, y, cv = (h.loc[idx, v].to_numpy(float) for v in (xc, yc, c))
    tert = tertiles(cv)
    r_all, p_all = stats.spearmanr(x, y)
    row = dict(C=c, x=xc, y=yc, N=len(idx), r_all=r_all, p_all=p_all, p_ad_all=RC.ad_split(x, y),
               idx=idx, tert=tert)
    row["bins"] = []
    for k in range(3):
        m = tert == k
        row["bins"].append(dict(bin_stats(x[m], y[m], rng), c_lo=cv[m].min(), c_hi=cv[m].max()))
    row["r_pool"], row["p_pool"], row["p_Q"] = fisher_pool([b["r"] for b in row["bins"]],
                                                           [b["N"] for b in row["bins"]])
    partial = c not in (xc, yc)
    row["ptau"], row["ptau_lo"], row["ptau_hi"], row["p_ptau"] = (
        D.partial_tau_test(x, y, cv, rng) if partial else (np.nan,) * 4)
    row["mc"] = monte_carlo(draws[xc][:, idx], draws[yc][:, idx], draws[c][:, idx], partial)
    return row


# ---------------------------------------------------------------- figures
def fmt_p(p):
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}" if p < 0.1 else f"p = {p:.2f}"


def fmt_c(c, v):
    return f"{v:.0f}" if c == DM else f"{v:.2f}"


def fmt_mc(a):
    return f"{a[0]:+.2f} [{a[1]:+.2f}, {a[2]:+.2f}]"


def draw(ax, h, err, row, small=False):
    c, xc, yc, idx, tert = row["C"], row["x"], row["y"], row["idx"], row["tert"]
    x, y = h.loc[idx, xc].to_numpy(float), h.loc[idx, yc].to_numpy(float)
    ex, ey = err[xc][:, idx], err[yc][:, idx]
    for k, b in enumerate(row["bins"]):
        m = tert == k
        lab = (f"{SHORT[c]} {fmt_c(c, b['c_lo'])}–{fmt_c(c, b['c_hi'])}: $r_S$ = {b['r']:+.2f} ({fmt_p(b['p'])}); "
               f"AD {fmt_p(b['p_ad'])}\n"
               f"MC $r_S$ {fmt_mc(row['mc']['r_bin'][:, k])}; p < 0.05 in {row['mc']['fsig_bin'][k]:.0%} "
               f"(AD {row['mc']['ad_fsig_bin'][k]:.0%})")
        ax.errorbar(x[m], y[m], xerr=ex[:, m], yerr=ey[:, m], fmt="none", ecolor=COLORS[k], elinewidth=0.6,
                    alpha=0.5, zorder=2)
        ax.scatter(x[m], y[m], s=13 if small else 22, color=COLORS[k], edgecolor="k", lw=0.3, label=lab, zorder=3)
        slope, icpt, _, _ = stats.theilslopes(y[m], x[m])
        xx = np.array([x[m].min(), x[m].max()])
        ax.plot(xx, icpt + slope * xx, color=COLORS[k], lw=1.4, zorder=4)
    fs = 6 if small else 7.5
    ax.legend(fontsize=fs - 0.7, loc="upper left", frameon=False, handletextpad=0.2, borderaxespad=0.2)
    ax.set_xlabel(AXIS[xc], fontsize=fs + 1)
    ax.set_ylabel(AXIS[yc], fontsize=fs + 1)
    partial = (f"partial $\\tau$ | {SHORT[c]} = {row['ptau']:+.2f} ({fmt_p(row['p_ptau'])})"
               if np.isfinite(row["ptau"]) else "partial $\\tau$: binned on the y quantity")
    mc_partial = f"\nMC partial $\\tau$ {fmt_mc(row['mc']['ptau'])}" if np.isfinite(row["ptau"]) else ""
    ax.text(0.5, -0.17 if small else -0.15,
            f"all {row['N']}: $r_S$ = {row['r_all']:+.2f} ({fmt_p(row['p_all'])}); AD {fmt_p(row['p_ad_all'])};  "
            f"{partial}\n"
            f"pooled within tertiles: $r$ = {row['r_pool']:+.2f} ({fmt_p(row['p_pool'])});  "
            f"tertiles agree: {fmt_p(row['p_Q'])}\n"
            f"through the errors (MC median [16, 84]%): all $r_S$ {fmt_mc(row['mc']['r_all'])}, "
            f"p < 0.05 in {row['mc']['fsig_all']:.0%} (AD {row['mc']['ad_fsig_all']:.0%}){mc_partial}",
            transform=ax.transAxes, ha="center", va="top", fontsize=fs)
    ax.tick_params(labelsize=fs)
    lo, hi = np.nanmin(y - ey[0]), np.nanmax(y + ey[1])
    ax.set_ylim(lo - 0.04 * (hi - lo), hi + 0.75 * (hi - lo))
    if xc == INC:
        ax.set_xlim(0, 90)


def title(row):
    return f"{NAMES[row['y']]} vs {NAMES[row['x']]}, binned by {NAMES[row['C']]}"


def grid(h, err, rows, shape, name, suptitle):
    fig, axes = plt.subplots(*shape, figsize=(4.6 * shape[1], 4.2 * shape[0]))
    for ax, row in zip(axes.flat, rows):
        draw(ax, h, err, row, small=True)
        ax.set_title(title(row), fontsize=8.5)
    fig.suptitle(suptitle, fontsize=11, y=1.0)
    fig.tight_layout(h_pad=3.4)
    fig.text(0.5, -0.01,
             "How to read each legend entry: tertile range of the binning quantity; Spearman rank correlation $r_S$ of "
             "the plotted pair in that tertile, with its two-sided p-value (null hypothesis: $r_S$ = 0).\n"
             "AD: two-sample Anderson–Darling p-value that the y quantity has the same distribution for hosts below and "
             "above the median of the x quantity (5,000 permutations); it also detects non-monotonic dependence.\n"
             "MC $r_S$: the same Spearman coefficient recomputed in 2,000 draws in which every value is resampled within "
             "its error bars and the tertiles are re-formed; median [16th, 84th percentile] of $r_S$ (not a p-value), "
             "and the fraction of draws with Spearman p < 0.05 (in brackets: with AD p < 0.05).",
             ha="center", va="top", fontsize=9)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=170, bbox_inches="tight")
    plt.close(fig)


def heatmap(rows, name, suptitle):
    cols = ["low tertile", "mid tertile", "high tertile", "all", "pooled", "partial τ"]
    row_h = 0.8
    fig, ax = plt.subplots(figsize=(10.5, row_h * len(rows) + 1.6))
    short = lambda p: fmt_p(p).replace("p = ", "").replace("p ", "")
    for j, row in enumerate(rows):
        cells = [(b["r"], b["p"], b["p_ad"], row["mc"]["r_bin"][0, k], b["q"] < 0.05)
                 for k, b in enumerate(row["bins"])]
        cells += [(row["r_all"], row["p_all"], row["p_ad_all"], row["mc"]["r_all"][0], False),
                  (row["r_pool"], row["p_pool"], np.nan, np.nan, False),
                  (row["ptau"], row["p_ptau"], np.nan, row["mc"]["ptau"][0], False)]
        for i, (r, p, p_ad, mc, strong) in enumerate(cells):
            color = plt.cm.RdBu_r(0.5 + 0.5 * np.clip(r, -1, 1)) if np.isfinite(r) else "#eeeeee"
            ax.add_patch(Rectangle((i, j), 1, 1, color=color, ec="white"))
            if strong:
                ax.add_patch(Rectangle((i + 0.03, j + 0.05), 0.94, 0.9, fill=False, ec="k", lw=2))
            if np.isfinite(r):
                txt = f"{r:+.2f} ({short(p)})"
                txt += f"\nAD {short(p_ad)}" if np.isfinite(p_ad) else ""
                txt += f"\nMC {mc:+.2f}" if np.isfinite(mc) else ""
                ax.text(i + 0.5, j + 0.5, txt, ha="center", va="center", fontsize=7.5,
                        weight="bold" if p < 0.05 else "normal")
            else:
                ax.text(i + 0.5, j + 0.5, "n/a", ha="center", va="center", fontsize=7.5, color="0.4")
    ax.set_xlim(0, len(cols))
    ax.set_ylim(len(rows), 0)
    ax.set_xticks(np.arange(len(cols)) + 0.5, cols, fontsize=8.5)
    ax.xaxis.tick_top()
    ax.set_yticks(np.arange(len(rows)) + 0.5, [title(r) + f"  (N = {r['N']})" for r in rows], fontsize=8)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(suptitle, fontsize=10, pad=24)
    fig.text(0.5, 0.01, "Cells: Spearman $r_S$ (p) for the tertiles, all hosts and the Fisher-pooled within-tertile "
             "value; Kendall partial $\\tau$ given the binning quantity.\nAD: two-sample Anderson–Darling p that the y "
             "quantity has the same distribution below and above the median of the x quantity.\nMC: median $r_S$ over "
             "2,000 draws through the measurement errors with tertiles re-formed. Bold: Spearman p < 0.05. Black frame: "
             "Benjamini–Hochberg q < 0.05 over the tertiles of the grid.", ha="center", va="bottom", fontsize=7.5)
    height = row_h * len(rows) + 1.6
    fig.subplots_adjust(left=0.3, right=0.99, top=1 - 1.1 / height, bottom=1.0 / height)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=170)
    plt.close(fig)


def figures(h, err, rows_a, rows_b):
    for row in rows_a + rows_b:
        fig, ax = plt.subplots(figsize=(5.4, 4.6))
        draw(ax, h, err, row)
        ax.set_title(title(row), fontsize=10)
        fig.savefig(OUT / f"binned_by_{row['C']}__{row['y']}_vs_{row['x']}.png", dpi=200, bbox_inches="tight")
        plt.close(fig)
    grid(h, err, rows_a, (4, 3), "grid", "Pairs of z, M*, SFR and DM_host within tertiles of the fourth quantity")
    grid(h, err, rows_b, (4, 4), "grid_inclination",
         "z, M*, SFR and DM_host against inclination within tertiles of each quantity")
    heatmap(rows_a, "heatmap", "Correlations within tertiles: z, M*, SFR, DM_host")
    heatmap(rows_b, "heatmap_inclination", "Correlations with inclination within tertiles of z, M*, SFR, DM_host")


# ---------------------------------------------------------------- report
def table(rows, label) -> list[str]:
    md = []
    for c in VARS:
        md += ["", f"### {label}: binned by {NAMES[c]}", "",
               "| pair | tertile | range of " + NAMES[c] + " | N | r_S (p) | BH q | AD p (BH q) | bootstrap 95% | "
               "delete-1 range | delete-3 range (frac. p < 0.05) | MC r_S [16, 84] (frac. p < 0.05; AD) |",
               "|---|---|---|---|---|---|---|---|---|---|---|"]
        for row in (r for r in rows if r["C"] == c):
            pair = f"{NAMES[row['y']]} vs {NAMES[row['x']]}"
            for k, b in enumerate(row["bins"]):
                md.append(f"| {pair if k == 0 else ''} | {['low', 'mid', 'high'][k]} | "
                          f"{fmt_c(c, b['c_lo'])} to {fmt_c(c, b['c_hi'])} | {b['N']} | {b['r']:+.2f} ({fmt_p(b['p'])}) | "
                          f"{b['q']:.2f} | {fmt_p(b['p_ad'])} ({b['q_ad']:.2f}) | [{b['lo']:+.2f}, {b['hi']:+.2f}] | "
                          f"{b['jk_min']:+.2f} to {b['jk_max']:+.2f} | "
                          f"{b['jk3_min']:+.2f} to {b['jk3_max']:+.2f} ({b['jk3_fsig']:.2f}) | "
                          f"{fmt_mc(row['mc']['r_bin'][:, k])} ({row['mc']['fsig_bin'][k]:.2f}; "
                          f"{row['mc']['ad_fsig_bin'][k]:.2f}) |")
            partial = (f"partial tau given {NAMES[c]}: {row['ptau']:+.2f} [{row['ptau_lo']:+.2f}, {row['ptau_hi']:+.2f}] "
                       f"({fmt_p(row['p_ptau'])}); MC {fmt_mc(row['mc']['ptau'])}" if np.isfinite(row["ptau"]) else "")
            md.append(f"| | **summary** | | {row['N']} | all: {row['r_all']:+.2f} ({fmt_p(row['p_all'])}) | | "
                      f"all: {fmt_p(row['p_ad_all'])} | "
                      f"pooled: {row['r_pool']:+.2f} ({fmt_p(row['p_pool'])}); tertiles agree: {fmt_p(row['p_Q'])} | "
                      f"{partial} | | MC all: {fmt_mc(row['mc']['r_all'])} ({row['mc']['fsig_all']:.2f}; "
                      f"{row['mc']['ad_fsig_all']:.2f}) |")
    return md


def report(h, rows_a, rows_b, pairs):
    md = ["# Correlations of z, M*, SFR, DM_host and inclination within tertiles of one quantity", "",
          "Generated by `scripts/binned_correlations.py`. Each panel uses every host with the values it needs. "
          "DM_host is rest frame. BH q: Benjamini-Hochberg over the tertile correlations of each grid (36 in grid A, "
          "48 in grid B), separately for Spearman and AD. AD: two-sample Anderson-Darling p that the y quantity has the "
          "same distribution below and above the median of the x quantity (5,000 permutations). MC: 2,000 draws "
          "through the measurement errors, tertiles re-formed in each draw; median and 16-84% range of r_S, and the "
          "fractions of draws with Spearman p < 0.05 and with AD p < 0.05.", "",
          "## Whole-sample correlations", "",
          "| y vs x | N | r_S (p) | AD p | MC r_S [16, 84] (frac. p < 0.05; AD) |", "|---|---|---|---|---|"]
    md += [f"| {NAMES[a]} vs {NAMES[b]} | {n} | {r:+.2f} ({fmt_p(p)}) | {fmt_p(p_ad)} | {fmt_mc(mc)} ({f:.2f}; {f_ad:.2f}) |"
           for a, b, n, r, p, p_ad, mc, f, f_ad in pairs]
    md += ["", "## Grid A: pairs of z, M*, SFR and DM_host"] + table(rows_a, "A")
    md += ["", "## Grid B: each quantity against inclination"] + table(rows_b, "B")
    (OUT / "RESULTS.md").write_text("\n".join(md) + "\n", encoding="utf-8")


def whole_sample(h, draws):
    out = []
    for a, b in combinations(VARS + [INC], 2):
        idx = h.index[h[[a, b]].notna().all(axis=1)].to_numpy()
        xc, yc = (b, a) if b == INC else (a, b)
        x, y = h.loc[idx, xc].to_numpy(float), h.loc[idx, yc].to_numpy(float)
        r, p = stats.spearmanr(x, y)
        mc = [stats.spearmanr(draws[a][d, idx], draws[b][d, idx]) for d in range(N_MC)]
        ad = [RC.ad_split(draws[xc][d, idx], draws[yc][d, idx], exact=False) for d in range(N_MC)]
        rr = np.array([m[0] for m in mc])
        out.append((yc, xc, len(idx), r, p, RC.ad_split(x, y), np.percentile(rr, [50, 16, 84]),
                    np.mean([m[1] < 0.05 for m in mc]), np.mean(np.array(ad) < 0.05)))
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(5)
    h = load_hosts()
    draws = measurement_draws(h, rng)
    err = error_bars(h, draws)
    rows_a = [analyse(h, draws, c, xc, yc, rng)
              for c in VARS for xc, yc in combinations([v for v in VARS if v != c], 2)]
    print("grid A done", flush=True)
    rows_b = [analyse(h, draws, c, INC, yc, rng) for c in VARS for yc in VARS]
    print("grid B done", flush=True)
    for rows in (rows_a, rows_b):
        bins = [b for r in rows for b in r["bins"]]
        for b, q, q_ad in zip(bins, stats.false_discovery_control([b["p"] for b in bins]),
                              stats.false_discovery_control([b["p_ad"] for b in bins])):
            b["q"], b["q_ad"] = q, q_ad
    flat = []
    for grid_name, rows in (("A", rows_a), ("B", rows_b)):
        for r in rows:
            base = {k: v for k, v in r.items() if k not in ("bins", "idx", "tert", "mc")}
            base.update(grid=grid_name, mc_r_all=r["mc"]["r_all"][0], mc_ptau=r["mc"]["ptau"][0],
                        mc_ad_fsig_all=r["mc"]["ad_fsig_all"])
            for k, b in enumerate(r["bins"]):
                flat.append(dict(base, tertile=k, mc_r=r["mc"]["r_bin"][0, k], mc_r16=r["mc"]["r_bin"][1, k],
                                 mc_r84=r["mc"]["r_bin"][2, k], mc_fsig=r["mc"]["fsig_bin"][k],
                                 mc_ad_fsig=r["mc"]["ad_fsig_bin"][k],
                                 **{f"bin_{kk}": vv for kk, vv in b.items()}))
    pd.DataFrame(flat).to_csv(OUT / "binned_correlations.csv", index=False)
    figures(h, err, rows_a, rows_b)
    report(h, rows_a, rows_b, whole_sample(h, draws))
    print("done", flush=True)


if __name__ == "__main__":
    main()
