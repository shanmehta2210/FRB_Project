"""Outlier-robust correlations of host properties with axis ratio q (equivalently cos i).

For each property: Spearman r_S (Fisher 95% interval), Kendall tau_b, bootstrap 95% intervals of both,
leave-one-out jackknife range of r_S (and the most influential host), the Theil-Sen slope with its
95% interval, and a two-sample Anderson-Darling test of whether the property has the same distribution
below and above the median of q (permutation p). DM_host is also tested after removing the foreground-halo DM of intervening galaxies
(scripts/foreground_halo_dm.py: corrected version of the Bhardwaj, Snelders et al. 2025 method, which is kept
as the "published method" comparison).

Inputs:  plots/plots_null/v2/host_correlations/host_props.csv
         Publication/materials/foreground_halos/foreground_dm_per_host.csv (+ foreground_dm_controls.csv)
Outputs: Publication/materials/appendix/robust_correlations/{robust_correlations.csv, ROBUST.md}
"""
from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
HOSTS = ROOT / "plots" / "plots_null" / "v2" / "host_correlations" / "host_props.csv"
FG = ROOT / "Publication" / "materials" / "foreground_halos"
OUT = ROOT / "Publication" / "materials" / "appendix" / "robust_correlations"
N_BOOT = 10_000
N_PERM_AD = 5_000
SEED = 42
AD_RNG = np.random.default_rng(11)
BULGE_FITS = ["20190608B", "20220725A", "20231226A"]

PROPS = [
    ("dm_host_rest", r"DM$_{\rm host}$ (rest frame)"),
    ("dm_host_fg", r"DM$_{\rm host}$, foreground halos removed"),
    ("dm_host_fg_excess", r"DM$_{\rm host}$, foreground excess over random sightlines removed"),
    ("dm_host_fg_cond", r"DM$_{\rm host}$, foreground removed (all photo-$z$ galaxies in front)"),
    ("dm_host_fg_spec", r"DM$_{\rm host}$, spectroscopic foreground only"),
    ("dm_host_fg_excess_published", r"DM$_{\rm host}$, foreground excess removed (published method)"),
    ("log_abs_rm", r"$\log|{\rm RM}|$"),
    ("log_tau", r"$\log\tau_{\rm scatt}$"),
    ("lin_pol", "Linear polarization"),
    ("av_old", r"$A_V$, old stars"),
]


def fisher_ci(r, n):
    z, s = np.arctanh(r), 1.96 / np.sqrt(n - 3)
    return float(np.tanh(z - s)), float(np.tanh(z + s))


def load() -> pd.DataFrame:
    """Host table with the foreground-corrected DM_host columns. dm_host_fg_excess removes only the excess
    of the sightline's halo DM over the random-sightline mean, which <DM_cosmic> already contains; the
    full subtractions (dm_host_fg, _cond, _spec) count that mean twice and are bounds."""
    h = pd.read_csv(HOSTS)
    f = pd.read_csv(FG / "foreground_dm_per_host.csv")[
        ["frb", "in_ls", "dm_fg", "dm_fg_conditional", "dm_fg_specz_only", "dm_fg_expected",
         "dm_fg_published", "dm_fg_expected_published"]]
    h = h.merge(f, on="frb", how="left")
    ls = h.in_ls.astype(str).isin(["dr9", "dr10"])
    for col, src in (("dm_host_fg", "dm_fg"), ("dm_host_fg_cond", "dm_fg_conditional"),
                     ("dm_host_fg_spec", "dm_fg_specz_only"), ("dm_host_fg_published", "dm_fg_published")):
        h[col] = np.where(ls, h.dm_host_rest - h[src].fillna(0) * (1 + h.z), np.nan)
    for col, fg, ex in (("dm_host_fg_excess", "dm_fg", "dm_fg_expected"),
                        ("dm_host_fg_excess_published", "dm_fg_published", "dm_fg_expected_published")):
        h[col] = np.where(ls, h.dm_host_rest - (h[fg] - h[ex]) * (1 + h.z), np.nan)
    return h


def ad_split(x, y, exact=True) -> float:
    """Two-sample Anderson-Darling p that y has the same distribution below and above the median of x.
    exact: permutation p; otherwise asymptotic p, which scipy clips to [0.001, 0.25]."""
    lower = np.zeros(len(x), bool)
    lower[np.argsort(x, kind="stable")[:len(x) // 2]] = True
    method = stats.PermutationMethod(n_resamples=N_PERM_AD, rng=AD_RNG) if exact else None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return float(stats.anderson_ksamp([y[lower], y[~lower]], method=method).pvalue)


def robust_row(frb, x, y, rng) -> dict:
    ok = np.isfinite(x) & np.isfinite(y)
    frb, x, y = frb[ok], x[ok], y[ok]
    n = len(x)
    rs, ps = stats.spearmanr(x, y)
    p_ad = ad_split(x, y)
    tau = stats.kendalltau(x, y)
    lo, hi = fisher_ci(rs, n)
    idx = rng.integers(0, n, size=(N_BOOT, n))
    bs = np.array([stats.spearmanr(x[i], y[i])[0] for i in idx])
    bt = np.array([stats.kendalltau(x[i], y[i])[0] for i in idx])
    jk = np.array([stats.spearmanr(np.delete(x, k), np.delete(y, k))[0] for k in range(n)])
    jk_se = np.sqrt((n - 1) / n * np.sum((jk - jk.mean()) ** 2))
    k_inf = int(np.argmax(np.abs(jk - rs)))
    ts = stats.theilslopes(y, x, alpha=0.95)
    return dict(N=n, r_S=rs, p_S=ps, p_AD=p_ad, rS_fisher_lo=lo, rS_fisher_hi=hi,
                rS_boot_lo=np.nanpercentile(bs, 2.5), rS_boot_hi=np.nanpercentile(bs, 97.5),
                tau=tau.statistic, p_tau=tau.pvalue,
                tau_boot_lo=np.nanpercentile(bt, 2.5), tau_boot_hi=np.nanpercentile(bt, 97.5),
                rS_jack_min=jk.min(), rS_jack_max=jk.max(), rS_jack_se=jk_se,
                most_influential=frb[k_inf], rS_without_it=jk[k_inf],
                theil_sen_slope=ts.slope, ts_lo=ts.low_slope, ts_hi=ts.high_slope,
                frac_boot_le_m044=float(np.mean(bs <= -0.44)))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    h = load()
    rng = np.random.default_rng(SEED)
    rows = []
    for col, lab in PROPS:
        if col not in h or h[col].notna().sum() < 6:
            continue
        r = robust_row(h.frb.to_numpy(), h.q.to_numpy(float), h[col].to_numpy(float), rng)
        rows.append(dict(property=col, label=lab, **r))
        print(f"{col:20s} N={r['N']:2d} rS={r['r_S']:+.2f} [{r['rS_boot_lo']:+.2f},{r['rS_boot_hi']:+.2f}] "
              f"tau={r['tau']:+.2f} (p={r['p_tau']:.2f}) jack [{r['rS_jack_min']:+.2f},{r['rS_jack_max']:+.2f}] "
              f"TS={r['theil_sen_slope']:.1f} [{r['ts_lo']:.1f},{r['ts_hi']:.1f}] P(rS<=-0.44)={r['frac_boot_le_m044']:.3f}")
    subs = [("m_r < 21", h.m_r < 21), ("n < 2", h.sersic_n < 2), ("bulge fits removed", ~h.frb.isin(BULGE_FITS))]
    for col in ("dm_host_rest", "dm_host_fg_excess"):
        for lab, sel in subs:
            s = h[sel]
            r = robust_row(s.frb.to_numpy(), s.q.to_numpy(float), s[col].to_numpy(float), rng)
            rows.append(dict(property=f"{col} [{lab}]", label=f"{dict(PROPS)[col]}, {lab}", **r))
            print(f"{col} [{lab}] N={r['N']} rS={r['r_S']:+.2f} [{r['rS_boot_lo']:+.2f},{r['rS_boot_hi']:+.2f}] tau p={r['p_tau']:.2f}")
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "robust_correlations.csv", index=False)
    sanity = {}
    for col in ("log_tau", "log_abs_rm"):
        ok = h.dm_host_rest.notna() & h[col].notna()
        rs = stats.spearmanr(h.loc[ok, "dm_host_rest"], h.loc[ok, col])
        sanity[col] = (int(ok.sum()), float(rs[0]), float(rs[1]))
    print("DM_host sanity correlations:", sanity)
    pd.DataFrame([dict(property=k, N=v[0], r_S=v[1], p=v[2]) for k, v in sanity.items()]).to_csv(
        OUT / "dm_host_sanity.csv", index=False)
    ls = h.in_ls.astype(str).isin(["dr9", "dr10"])
    fg = h.loc[ls, "dm_fg"] * (1 + h.loc[ls, "z"])
    imp = dict(n_ls=int(ls.sum()), n_out=int((~ls & h.z.notna()).sum()), median=float(fg.median()),
               p84=float(fg.quantile(0.84)), max=float(fg.max()), n_gt50=int((fg > 50).sum()),
               n_gt100=int((fg > 100).sum()),
               frac_median=float((fg / h.loc[ls, "dm_host_rest"].where(h.loc[ls, "dm_host_rest"] > 0)).median()),
               n_neg_after=int((h.loc[ls, "dm_host_fg"] < 0).sum()),
               n_neg_after_excess=int((h.loc[ls, "dm_host_fg_excess"] < 0).sum()),
               median_published=float((h.loc[ls, "dm_fg_published"] * (1 + h.loc[ls, "z"])).median()),
               r_fg_q=float(stats.spearmanr(h.loc[ls, "q"], h.loc[ls, "dm_fg"])[0]),
               p_fg_q=float(stats.spearmanr(h.loc[ls, "q"], h.loc[ls, "dm_fg"])[1]),
               p_fg_q_ad=ad_split(h.loc[ls, "q"].to_numpy(float), h.loc[ls, "dm_fg"].to_numpy(float)),
               expected_median=float((h.loc[ls, "dm_fg_expected"] * (1 + h.loc[ls, "z"])).median()))
    pd.Series(imp).to_csv(OUT / "foreground_impact.csv", header=["value"])
    print(imp)
    md = ["# Outlier-robust correlations with axis ratio $q$", "",
          f"Bootstrap: {N_BOOT:,} resamples (seed {SEED}). Jackknife: leave one host out. "
          "Theil-Sen slope in property units per unit $q$, 95% interval. "
          "$P(r_S \\le -0.44)$: bootstrap fraction at or below the Marnoch et al. value. "
          "AD $p$: two-sample Anderson-Darling test of the property below vs above the median $q$ "
          f"({N_PERM_AD:,} permutations).", "",
          "| property | $N$ | $r_S$ (Fisher 95%) | $r_S$ $p$ | AD $p$ | $r_S$ bootstrap 95% | Kendall $\\tau$ ($p$) | $\\tau$ bootstrap 95% | jackknife $r_S$ range | most influential host ($r_S$ without) | Theil-Sen slope (95%) | $P(r_S\\le-0.44)$ |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r['label']} | {r['N']} | {r['r_S']:+.2f} ({r['rS_fisher_lo']:+.2f}, {r['rS_fisher_hi']:+.2f}) | "
                  f"{r['p_S']:.3f} | {r['p_AD']:.3f} | "
                  f"({r['rS_boot_lo']:+.2f}, {r['rS_boot_hi']:+.2f}) | {r['tau']:+.2f} ({r['p_tau']:.2f}) | "
                  f"({r['tau_boot_lo']:+.2f}, {r['tau_boot_hi']:+.2f}) | {r['rS_jack_min']:+.2f} to {r['rS_jack_max']:+.2f} | "
                  f"{r['most_influential']} ({r['rS_without_it']:+.2f}) | {r['theil_sen_slope']:.3g} ({r['ts_lo']:.3g}, {r['ts_hi']:.3g}) | "
                  f"{r['frac_boot_le_m044']:.3f} |")
    md += ["", "## Foreground-halo DM over the sample (rest frame of the host)", "",
           f"LS footprint: {imp['n_ls']} hosts; outside: {imp['n_out']}. Median {imp['median']:.0f} "
           f"(published method: {imp['median_published']:.0f}), 84th percentile {imp['p84']:.0f}, "
           f"maximum {imp['max']:.0f} pc cm$^{{-3}}$; {imp['n_gt50']} hosts above 50 and {imp['n_gt100']} above 100. "
           f"Median fraction of DM$_{{\\rm host}}$: {imp['frac_median']:.2f}. Hosts negative after removing the excess: "
           f"{imp['n_neg_after_excess']}; after the full subtraction: {imp['n_neg_after']}. "
           f"Spearman of DM$_{{\\rm fg}}$ with $q$: $r_S = {imp['r_fg_q']:+.2f}$ ($p = {imp['p_fg_q']:.2f}$; "
           f"AD $p = {imp['p_fg_q_ad']:.2f}$). "
           f"Median expected (random-sightline) term: {imp['expected_median']:.0f} pc cm$^{{-3}}$."]
    (OUT / "ROBUST.md").write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
