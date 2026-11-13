"""Are the few significant tertile correlations of the binned grids more than chance?

The tertile tests of each grid are strongly dependent (the same pair appears under several binnings and
the panels share hosts), so the null distributions are built by permutation, which keeps that dependence:
  grid B (48 tertiles)               inclination shuffled across hosts: i independent of z, M*, SFR,
                                     DM_host, with every relation among those four kept
  grid A, DM_host vs z/M*/SFR (18)   DM_host shuffled across hosts: DM_host independent of z, M*, SFR,
                                     with the selection trends among those three kept
  grid A, z/M*/SFR pairs (18)        y shuffled within the tertiles of the binning quantity: no
                                     association at fixed C
  grid A, all 36                     both at once
For each family and for Spearman and Anderson-Darling (median split, asymptotic p): the number of tertiles
with p < 0.05, the smallest p, and Fisher's combined statistic, each with its permutation p-value, and
Westfall-Young single-step minP family-wise adjusted p-values for every tertile.

Output: Publication/materials/appendix/binned_correlations/
  GLOBAL_TESTS.md, global_tests.csv, global_tests_tertiles.csv, global_tests.{pdf,png}
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
import binned_correlations as B  # noqa: E402
import robust_correlations as RC  # noqa: E402

N_PERM = 5_000
ALPHA = 0.05
TERT = ["low", "mid", "high"]

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "stix", "font.serif": ["STIXGeneral"]})


def panels():
    a = [("A", c, xc, yc) for c in B.VARS for xc, yc in combinations([v for v in B.VARS if v != c], 2)]
    b = [("B", c, B.INC, yc) for c in B.VARS for yc in B.VARS]
    return a + b


def prepare(h, specs):
    out = []
    for g, c, xc, yc in specs:
        idx = h.index[h[list(dict.fromkeys([c, xc, yc]))].notna().all(axis=1)].to_numpy()
        null = "inc" if g == "B" else "dm" if B.DM in (xc, yc) else "within"
        out.append(dict(grid=g, C=c, x=xc, y=yc, idx=idx, null=null))
    return out


def tertile_p(cols, pan, shuffled=None, rng=None):
    """Spearman and AD p of every tertile. cols maps variable -> full-length array; shuffled maps each null
    type ('dm', 'inc') to the permuted columns, and 'within' panels get y permuted inside each tertile."""
    ps, pa = [], []
    for p in pan:
        src = cols if shuffled is None or p["null"] == "within" else shuffled[p["null"]]
        x, y, cv = (src[v][p["idx"]] for v in (p["x"], p["y"], p["C"]))
        t = B.tertiles(cv)
        for k in range(3):
            m = np.flatnonzero(t == k)
            yy = y[m]
            if shuffled is not None and p["null"] == "within":
                yy = rng.permutation(yy)
            ps.append(stats.spearmanr(x[m], yy)[1])
            pa.append(RC.ad_split(x[m], yy, exact=False))
    return np.array(ps), np.array(pa)


def family_stats(p):
    return dict(n_sig=int(np.sum(p < ALPHA)), min_p=float(np.min(p)), fisher=float(-2 * np.sum(np.log(p))))


def main() -> None:
    rng = np.random.default_rng(17)
    h = B.load_hosts()
    pan = prepare(h, panels())
    labels = [(p["grid"], p["C"], p["x"], p["y"], k, p["null"]) for p in pan for k in range(3)]
    lab = pd.DataFrame(labels, columns=["grid", "C", "x", "y", "tertile", "null"])
    cols = {v: h[v].to_numpy(float) for v in B.VARS + [B.INC]}
    obs_s, obs_a = tertile_p(cols, pan)

    fams = {"A, all 36": (lab.grid == "A").to_numpy(),
            "A, DM_host vs z/M*/SFR (18)": (lab.null == "dm").to_numpy(),
            "A, z/M*/SFR pairs (18)": (lab.null == "within").to_numpy(),
            "B, all 48": (lab.grid == "B").to_numpy()}

    def shuffle(v):
        out, ok = dict(cols), np.isfinite(cols[v])
        out[v] = cols[v].copy()
        out[v][ok] = rng.permutation(cols[v][ok])
        return out

    null_s, null_a = np.empty((N_PERM, len(lab))), np.empty((N_PERM, len(lab)))
    for n in range(N_PERM):
        null_s[n], null_a[n] = tertile_p(cols, pan, {"dm": shuffle(B.DM), "inc": shuffle(B.INC)}, rng)
        if n % 500 == 0:
            print(n, flush=True)

    rows, adj = [], {}
    for name, m in fams.items():
        for test, obs, null in (("Spearman", obs_s, null_s), ("AD", obs_a, null_a)):
            o = family_stats(obs[m])
            nl = [family_stats(r[m]) for r in null]
            ns, mn, fi = (np.array([d[k] for d in nl]) for k in ("n_sig", "min_p", "fisher"))
            rows.append(dict(family=name, test=test, n_tests=int(m.sum()), n_sig=o["n_sig"],
                             expected_indep=ALPHA * m.sum(), null_mean_n_sig=ns.mean(),
                             null_p95_n_sig=float(np.percentile(ns, 95)),
                             p_count=(1 + np.sum(ns >= o["n_sig"])) / (1 + N_PERM),
                             p_count_binomial=float(stats.binom.sf(o["n_sig"] - 1, m.sum(), ALPHA)),
                             min_p=o["min_p"], p_minp=(1 + np.sum(mn <= o["min_p"])) / (1 + N_PERM),
                             fisher=o["fisher"], p_fisher=(1 + np.sum(fi >= o["fisher"])) / (1 + N_PERM)))
            adj[(name, test)] = (np.array([(1 + np.sum(mn <= p)) / (1 + N_PERM) for p in obs[m]]), ns)
    summary = pd.DataFrame(rows)

    t = lab.copy()
    t["p_S"], t["p_AD"] = obs_s, obs_a
    for test, col in (("Spearman", "p_S"), ("AD", "p_AD")):
        for fam, out in (("A, all 36", f"wy_{col}_gridA"), ("A, DM_host vs z/M*/SFR (18)", f"wy_{col}_A_DM"),
                         ("B, all 48", f"wy_{col}_gridB")):
            t[out] = np.nan
            t.loc[fams[fam], out] = adj[(fam, test)][0]
    t["tertile"] = t.tertile.map(dict(enumerate(TERT)))
    out = B.OUT
    summary.to_csv(out / "global_tests.csv", index=False)
    t.to_csv(out / "global_tests_tertiles.csv", index=False)
    figure(summary, adj, out)
    report(summary, t, out)
    print(summary.round(4).to_string())


def figure(summary, adj, out):
    fams = ["A, all 36", "A, DM_host vs z/M*/SFR (18)", "A, z/M*/SFR pairs (18)", "B, all 48"]
    fig, axes = plt.subplots(len(fams), 2, figsize=(10, 3.0 * len(fams)))
    for i, fam in enumerate(fams):
        for j, test in enumerate(("Spearman", "AD")):
            ax = axes[i, j]
            ns = adj[(fam, test)][1]
            r = summary[(summary.family == fam) & (summary.test == test)].iloc[0]
            bins = np.arange(-0.5, max(ns.max(), r.n_sig) + 1.5)
            ax.hist(ns, bins=bins, color="0.75", ec="w", density=True, label="permutation null")
            ax.axvline(r.n_sig, color="#D62728", lw=2, label=f"observed: {r.n_sig}")
            ax.axvline(r.expected_indep, color="k", ls="--", lw=1, label=f"5% of {r.n_tests} = {r.expected_indep:.1f}")
            ax.text(0.98, 0.62, f"P(count $\\geq$ observed) = {r.p_count:.3f}\n"
                    f"smallest p = {r.min_p:.4f}; family-wise p = {r.p_minp:.3f}\n"
                    f"Fisher combined: p = {r.p_fisher:.3f}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=7.5)
            ax.set_title(f"Grid {fam}: {test}", fontsize=9.5)
            ax.set_xlabel("number of tertiles with p < 0.05", fontsize=8.5)
            ax.set_ylabel("fraction of permutations", fontsize=8.5)
            ax.tick_params(labelsize=7.5)
            ax.legend(fontsize=7, frameon=False, loc="upper right")
    fig.suptitle("Global tests of the binned correlations against permutation nulls "
                 f"({N_PERM:,} permutations)", fontsize=11)
    fig.text(0.5, -0.01, "Grid A: DM$_{\\rm host}$ shuffled across hosts (keeps the selection trends among z, M*, SFR); "
             "the panels pairing z, M*, SFR: y shuffled within tertiles.\nGrid B: inclination shuffled "
             "across hosts. Family-wise p: Westfall–Young minP, the chance that any tertile in the family has a p this "
             "small.", ha="center", va="top", fontsize=8)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(out / f"global_tests.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def report(summary, t, out):
    n = B.NAMES
    md = ["# Global multiple-testing checks of the binned grids", "",
          f"Generated by `scripts/binned_global_tests.py`; {N_PERM:,} permutations. Spearman and Anderson-Darling "
          "(median split, asymptotic p clipped to [0.001, 0.25]) tertile p-values as in `RESULTS.md`.", "",
          "Null hypotheses, built by permutation so that the dependence between the tests is kept:", "",
          "- Grid B: inclination shuffled across hosts (i independent of z, M*, SFR and DM_host).",
          "- Grid A, DM_host vs z, M* or SFR (6 panels): DM_host shuffled across hosts (independent of z, M*, SFR; "
          "their mutual selection trends kept).",
          "- Grid A, pairs among z, M*, SFR (6 panels, three of them binned by DM_host): y shuffled within the "
          "tertiles of the binning quantity.", "",
          "| family | test | tests | observed p < 0.05 | 5% of tests | permutation mean (95th pct.) | P(count >= obs.), "
          "permutation | P(count >= obs.), binomial | smallest p | family-wise P(min p <= obs.) | Fisher p |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in summary.itertuples():
        md.append(f"| {r.family} | {r.test} | {r.n_tests} | {r.n_sig} | {r.expected_indep:.1f} | "
                  f"{r.null_mean_n_sig:.1f} ({r.null_p95_n_sig:.0f}) | {r.p_count:.3f} | {r.p_count_binomial:.3f} | "
                  f"{r.min_p:.4f} | {r.p_minp:.3f} | {r.p_fisher:.3f} |")
    md += ["", "## Tertiles with p < 0.05 in either test, with Westfall-Young family-wise adjusted p", "",
           "| grid | pair | binned by | tertile | Spearman p | AD p | adjusted Spearman (grid) | adjusted AD (grid) | "
           "adjusted Spearman (A, DM_host vs z/M*/SFR) | adjusted AD (A, DM_host vs z/M*/SFR) |",
           "|---|---|---|---|---|---|---|---|---|---|"]
    hits = t[(t.p_S < ALPHA) | (t.p_AD < ALPHA)]
    for r in hits.itertuples():
        g = "A" if r.grid == "A" else "B"
        md.append(f"| {g} | {n[r.y]} vs {n[r.x]} | {n[r.C]} | {r.tertile} | {r.p_S:.4f} | {r.p_AD:.4f} | "
                  f"{getattr(r, f'wy_p_S_grid{g}'):.3f} | {getattr(r, f'wy_p_AD_grid{g}'):.3f} | "
                  + (f"{r.wy_p_S_A_DM:.3f} | {r.wy_p_AD_A_DM:.3f} |" if np.isfinite(r.wy_p_S_A_DM) else "– | – |"))
    (out / "GLOBAL_TESTS.md").write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
