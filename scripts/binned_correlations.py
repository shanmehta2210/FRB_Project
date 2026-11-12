"""Pairwise correlations among DM_host, redshift, stellar mass and SFR, each pair tested within bins of one of the
other two quantities left over (4 binning quantities x 3 pairs = 12 combinations).

For each binning quantity C and each pair (x, y) of the remaining three, on the hosts with all four values:
  * Spearman r_S in tertiles of C, with bootstrap 95% intervals, delete-1 and delete-3 jackknife ranges, and
    Benjamini-Hochberg q-values over all 36 tertile correlations;
  * the pooled within-tertile correlation (inverse-variance mean of Fisher z, var = 1.06 / (n - 3)) and
    Cochran's Q test of whether the tertiles agree;
  * Kendall's partial tau of x and y given C (Akritas & Siebert 1996), p from shuffling y within strata of C;
  * for pairs with DM_host and C != DM_host, the median r_S over 2,000 draws of DM_cosmic (Macquart et al. 2020)
    and DM_MW (NE2001, 20%).

Output: Publication/materials/appendix/binned_correlations/{binned_correlations.csv, RESULTS.md,
        binned_by_<C>__<y>_vs_<x>.png, grid.pdf, grid.png}
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dm_host_physics as D  # noqa: E402
import robust_extras as E  # noqa: E402

OUT = D.ROOT / "Publication" / "materials" / "appendix" / "binned_correlations"
DM = D.DM
VARS = ["z", "logmstar", "log_sfr", DM]
NAMES = {"z": "z", "logmstar": "M*", "log_sfr": "SFR", DM: "DM_host"}
AXIS = {"z": "Redshift", "logmstar": r"log $M_\star$ (M$_\odot$)", "log_sfr": r"log SFR (M$_\odot$ yr$^{-1}$)",
        DM: r"DM$_{\rm host}$ (rest frame, pc cm$^{-3}$)"}
SHORT = {"z": "z", "logmstar": r"log $M_\star$", "log_sfr": "log SFR", DM: r"DM$_{\rm host}$"}
COLORS = ["#3b6fb6", "#e0a526", "#c0392b"]
N_BOOT = 2_000


def fisher_pool(r, n):
    """Pooled correlation over bins, its p-value, and Cochran's Q p-value for equal correlations."""
    r, n = np.asarray(r), np.asarray(n)
    z, w = np.arctanh(np.clip(r, -0.999, 0.999)), (n - 3) / 1.06
    zbar = np.sum(w * z) / np.sum(w)
    p = 2 * stats.norm.sf(abs(zbar) * np.sqrt(np.sum(w)))
    q = np.sum(w * (z - zbar) ** 2)
    return float(np.tanh(zbar)), float(p), float(stats.chi2.sf(q, len(r) - 1))


def bin_stats(x, y, rng):
    r, p = stats.spearmanr(x, y)
    boot = [stats.spearmanr(x[i], y[i])[0] for i in (rng.integers(0, len(x), len(x)) for _ in range(N_BOOT))]
    jack = [stats.spearmanr(np.delete(x, k), np.delete(y, k))[0] for k in range(len(x))]
    lo, hi = np.nanpercentile(boot, [2.5, 97.5])
    jk3_min, jk3_max, _, _, jk3_fsig = E.jackknife(x, y, 3)
    return dict(N=len(x), r=r, p=p, lo=lo, hi=hi, jk_min=min(jack), jk_max=max(jack),
                jk3_min=jk3_min, jk3_max=jk3_max, jk3_fsig=jk3_fsig)


def analyse(h, draws, c, xc, yc, rng) -> dict:
    edges = np.percentile(h[c], [100 / 3, 200 / 3])
    tert = np.digitize(h[c], edges)
    x, y, cv = h[xc].to_numpy(float), h[yc].to_numpy(float), h[c].to_numpy(float)
    r_all, p_all = stats.spearmanr(x, y)
    row = dict(C=c, x=xc, y=yc, N=len(h), r_all=r_all, p_all=p_all, edges=edges)
    bins = []
    for k in range(3):
        m = tert == k
        b = bin_stats(x[m], y[m], rng)
        b.update(c_lo=cv[m].min(), c_hi=cv[m].max())
        if DM in (xc, yc) and c != DM:
            other = x[m] if yc == DM else y[m]
            r_mc, _ = D.spearman_rows(other, draws[:, h.index[m]])
            b["r_mc"] = float(np.median(r_mc))
        bins.append(b)
    row["bins"] = bins
    row["r_pool"], row["p_pool"], row["p_Q"] = fisher_pool([b["r"] for b in bins], [b["N"] for b in bins])
    row["ptau"], row["ptau_lo"], row["ptau_hi"], row["p_ptau"] = D.partial_tau_test(x, y, cv, rng)
    if DM in (xc, yc) and c != DM:
        other = x if yc == DM else y
        r_mc, _ = D.spearman_rows(other, draws[:, h.index])
        row["r_all_mc"] = float(np.median(r_mc))
    return row


def fmt_p(p):
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}" if p < 0.1 else f"p = {p:.2f}"


def fmt_c(c, v):
    return f"{v:.2f}" if c != DM else f"{v:.0f}"


def draw(ax, h, row, small=False):
    c, xc, yc = row["C"], row["x"], row["y"]
    tert = np.digitize(h[c], row["edges"])
    for k, b in enumerate(row["bins"]):
        s = h[tert == k]
        lab = (f"{SHORT[c]} {fmt_c(c, b['c_lo'])}–{fmt_c(c, b['c_hi'])}: "
               f"$r_S$ = {b['r']:+.2f} ({fmt_p(b['p'])})")
        ax.scatter(s[xc], s[yc], s=14 if small else 22, color=COLORS[k], edgecolor="k", lw=0.3, label=lab, zorder=3)
        slope, icpt, _, _ = stats.theilslopes(s[yc], s[xc])
        xx = np.array([s[xc].min(), s[xc].max()])
        ax.plot(xx, icpt + slope * xx, color=COLORS[k], lw=1.4, zorder=2)
    fs = 6 if small else 7.5
    ax.legend(fontsize=fs - 0.5, loc="upper left", frameon=False, handletextpad=0.2, borderaxespad=0.2)
    ax.set_xlabel(AXIS[xc], fontsize=fs + 1)
    ax.set_ylabel(AXIS[yc], fontsize=fs + 1)
    ax.text(0.5, -0.17 if small else -0.15,
            f"all {row['N']}: $r_S$ = {row['r_all']:+.2f} ({fmt_p(row['p_all'])});  "
            f"partial $\\tau$ | {SHORT[c]} = {row['ptau']:+.2f} ({fmt_p(row['p_ptau'])})\n"
            f"pooled within tertiles: $r$ = {row['r_pool']:+.2f} ({fmt_p(row['p_pool'])});  "
            f"tertiles agree: {fmt_p(row['p_Q'])}",
            transform=ax.transAxes, ha="center", va="top", fontsize=fs)
    ax.tick_params(labelsize=fs)
    lo, hi = h[yc].min(), h[yc].max()
    ax.set_ylim(lo - 0.05 * (hi - lo), hi + 0.3 * (hi - lo))


def figures(h, rows):
    for row in rows:
        fig, ax = plt.subplots(figsize=(5.0, 4.2))
        draw(ax, h, row)
        ax.set_title(f"{NAMES[row['y']]} vs {NAMES[row['x']]}, binned by {NAMES[row['C']]}", fontsize=10)
        fig.savefig(OUT / f"binned_by_{row['C']}__{row['y']}_vs_{row['x']}.png", dpi=200, bbox_inches="tight")
        plt.close(fig)
    fig, axes = plt.subplots(4, 3, figsize=(13.5, 15.5))
    for ax, row in zip(axes.flat, rows):
        draw(ax, h, row, small=True)
        ax.set_title(f"binned by {NAMES[row['C']]}", fontsize=8.5)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"grid.{ext}", dpi=170, bbox_inches="tight")
    plt.close(fig)


def report(h, rows, pairs_full):
    md = ["# Pairwise correlations of DM_host, z, M* and SFR within bins of the other quantities", "",
          f"Generated by `scripts/binned_correlations.py`. {len(h)} hosts have all four quantities; tertiles of the "
          "binning quantity hold 14-16 hosts each. DM_host is rest frame, as in the paper. BH q: Benjamini-Hochberg "
          "over all 36 tertile correlations.", "",
          "## Whole-sample correlations", "",
          "| pair | N (all with both) | r_S (p) | N (all four) | r_S (p) |", "|---|---|---|---|---|"]
    md += [f"| {NAMES[a]} - {NAMES[b]} | {n1} | {r1:+.2f} ({fmt_p(p1)}) | {n2} | {r2:+.2f} ({fmt_p(p2)}) |"
           for a, b, n1, r1, p1, n2, r2, p2 in pairs_full]
    for c in VARS:
        md += ["", f"## Binned by {NAMES[c]}", "",
               "| pair | tertile | range of " + NAMES[c] + " | N | r_S (p) | BH q | bootstrap 95% | delete-1 range | "
               "delete-3 range (frac. p < 0.05) | DM_cosmic MC r_S |", "|---|---|---|---|---|---|---|---|---|---|"]
        for row in (r for r in rows if r["C"] == c):
            pair = f"{NAMES[row['y']]} vs {NAMES[row['x']]}"
            for k, b in enumerate(row["bins"]):
                mc = f"{b['r_mc']:+.2f}" if "r_mc" in b else ""
                md.append(f"| {pair if k == 0 else ''} | {['low', 'mid', 'high'][k]} | "
                          f"{fmt_c(c, b['c_lo'])} to {fmt_c(c, b['c_hi'])} | {b['N']} | {b['r']:+.2f} ({fmt_p(b['p'])}) | "
                          f"{b['q']:.2f} | [{b['lo']:+.2f}, {b['hi']:+.2f}] | {b['jk_min']:+.2f} to {b['jk_max']:+.2f} | "
                          f"{b['jk3_min']:+.2f} to {b['jk3_max']:+.2f} ({b['jk3_fsig']:.2f}) | {mc} |")
            mc = f"; DM_cosmic MC r_S {row['r_all_mc']:+.2f}" if "r_all_mc" in row else ""
            md.append(f"| | **summary** | | {row['N']} | all: {row['r_all']:+.2f} ({fmt_p(row['p_all'])}) | | "
                      f"pooled within tertiles: {row['r_pool']:+.2f} ({fmt_p(row['p_pool'])}); tertiles agree: "
                      f"{fmt_p(row['p_Q'])} | partial tau given {NAMES[c]}: {row['ptau']:+.2f} "
                      f"[{row['ptau_lo']:+.2f}, {row['ptau_hi']:+.2f}] ({fmt_p(row['p_ptau'])}){mc} | | |")
    (OUT / "RESULTS.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    full = D.load().reset_index(drop=True)
    rng = np.random.default_rng(5)
    draws = D.dm_host_draws(full, rng)
    h = full[np.all([np.isfinite(full[v]) for v in VARS], axis=0)]
    pairs_full = []
    for a, b in combinations(VARS, 2):
        s = full[np.isfinite(full[a]) & np.isfinite(full[b])]
        pairs_full.append((a, b, len(s), *stats.spearmanr(s[a], s[b]), len(h), *stats.spearmanr(h[a], h[b])))
    rows = []
    for c in VARS:
        for xc, yc in combinations([v for v in VARS if v != c], 2):
            rows.append(analyse(h, draws, c, xc, yc, rng))
            print(f"by {c}: {yc} vs {xc} done", flush=True)
    q = stats.false_discovery_control([b["p"] for r in rows for b in r["bins"]])
    for b, qq in zip((b for r in rows for b in r["bins"]), q):
        b["q"] = qq
    flat = []
    for r in rows:
        base = {k: v for k, v in r.items() if k not in ("bins", "edges")}
        for k, b in enumerate(r["bins"]):
            flat.append(dict(base, tertile=k, **{f"bin_{kk}": vv for kk, vv in b.items()}))
    pd.DataFrame(flat).to_csv(OUT / "binned_correlations.csv", index=False)
    figures(h, rows)
    report(h, rows, pairs_full)


if __name__ == "__main__":
    main()
