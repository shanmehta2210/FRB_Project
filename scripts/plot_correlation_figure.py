"""Figure 3: correlations of DM_host and burst properties with host axis ratio.

(a) Rest-frame DM_host versus axis ratio q, before (filled) and after (open) removing the foreground-halo
    DM of intervening galaxies, with the Theil-Sen line and its 95% slope interval.
(b) Spearman r_S with bootstrap 95% intervals for DM_host (all hosts, foreground-corrected, subsamples)
    and for the propagation-related burst properties plus old-star extinction.

Inputs: scripts/robust_correlations.py outputs (and its load() for the per-host values)
Output: Publication/materials/figures/fig3_correlations.{pdf,png}
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_paper_inclination_figures as P  # noqa: E402  (style)
import robust_correlations as R  # noqa: E402

OUT = P.REPO_ROOT / "Publication" / "materials" / "figures"

DM_ROWS = [("dm_host_rest", "All hosts"), ("dm_host_fg", "Foreground removed"),
           ("dm_host_fg_excess", "Foreground excess removed"), ("dm_host_rest [m_r < 21]", r"$m_r < 21$"),
           ("dm_host_rest [n < 2]", r"$n < 2$"), ("dm_host_rest [bulge fits removed]", "Bulge fits removed")]
PROP_ROWS = [("log_abs_rm", r"$|{\rm RM}|$"), ("log_tau", "Scattering time"),
             ("lin_pol", "Linear polarization"), ("av_old", r"$A_V$, old stars")]


def panel_a(ax, h: pd.DataFrame) -> None:
    d = h[h.dm_host_rest.notna()]
    x, y = d.q.to_numpy(float), d.dm_host_rest.to_numpy(float)
    yc = d.dm_host_fg.to_numpy(float)
    has = np.isfinite(yc) & (np.abs(yc - y) > 1)
    ax.vlines(x[has], yc[has], y[has], color=P.FRB_C, lw=0.5, alpha=0.5, zorder=2)
    ax.scatter(x[has], yc[has], s=12, facecolor="white", edgecolor=P.FRB_C, lw=0.7, zorder=3,
               label="Foreground halos removed")
    ax.scatter(x, y, s=14, color=P.FRB_C, lw=0, zorder=4, label=r"DM$_{\rm host}$")
    ts = stats.theilslopes(y, x, alpha=0.95)
    xx = np.linspace(0.2, 1.0, 50)
    xm, ym = np.median(x), np.median(y)
    lo, hi = ym + ts.low_slope * (xx - xm), ym + ts.high_slope * (xx - xm)
    ax.fill_between(xx, np.minimum(lo, hi), np.maximum(lo, hi), color=P.FRB_C, alpha=0.15, lw=0, zorder=1)
    ax.plot(xx, ts.intercept + ts.slope * xx, color=P.FRB_C, lw=1.0, zorder=2, label="Theil--Sen fit")
    rs = stats.spearmanr(x, y)
    tau = stats.kendalltau(x, y)
    ax.text(0.97, 0.97, f"$N = {len(x)}$\n$r_S = {rs[0]:.2f}$, $p = {rs[1]:.2f}$\n"
            rf"$\tau = {tau.statistic:.2f}$, $p = {tau.pvalue:.2f}$",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.5)
    ax.set_xlim(0.2, 1.0)
    ax.set_ylim(-450, 1150)
    ax.set_xlabel(r"Axis ratio $q$")
    ax.set_ylabel(r"DM$_{\rm host}$ (pc cm$^{-3}$, rest frame)")
    ax.legend(loc="upper left", frameon=False, fontsize=6)
    ax.text(0.03, 0.03, "(a)", transform=ax.transAxes, fontsize=9)


def panel_b(ax, rob: pd.DataFrame) -> None:
    ext = pd.read_csv(R.OUT / "robust_extras.csv").set_index("label")
    rob = rob.copy()
    rob["p_AD"] = rob.label.map(ext.p_AD)
    for k, lab in (("dm_host_rest [m_r < 21]", "m_r < 21"), ("dm_host_rest [n < 2]", "n < 2"),
                   ("dm_host_rest [bulge fits removed]", "bulge fits removed")):
        rob.loc[rob.property == k, "p_AD"] = ext.loc[rf"DM$_{{\rm host}}$, {lab}", "p_AD"]
    r = rob.set_index("property")
    rows = [(lab, r.loc[k], P.FRB_C) for k, lab in DM_ROWS if k in r.index]
    n_dm = len(rows)
    rows += [(lab, r.loc[k], P.HSC_C) for k, lab in PROP_ROWS]
    y = np.arange(len(rows))[::-1].astype(float)
    y[n_dm:] -= 0.6
    for yy, (lab, s, col) in zip(y, rows):
        filled = s.p_S < 0.05
        ax.errorbar(s.r_S, yy, xerr=[[s.r_S - s.rS_boot_lo], [s.rS_boot_hi - s.r_S]], fmt="o", ms=3.5, color=col,
                    mfc=col if filled else "white", lw=0.9, capsize=1.5, zorder=3)
        ax.text(1.03, yy, f"{int(s.N)}", transform=ax.get_yaxis_transform(), va="center", fontsize=6)
        for xt, p in ((1.14, s.p_S), (1.27, s.p_AD)):
            ax.text(xt, yy, f"{p:.2f}" if p >= 0.01 else f"{p:.3f}", transform=ax.get_yaxis_transform(),
                    va="center", fontsize=6)
    ax.axvline(0, color="0.7", lw=0.5, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([row[0] for row in rows], fontsize=6.5)
    ax.text(1.03, y.max() + 0.9, "$N$", transform=ax.get_yaxis_transform(), va="bottom", fontsize=6)
    ax.text(1.14, y.max() + 0.9, "$p_S$", transform=ax.get_yaxis_transform(), va="bottom", fontsize=6)
    ax.text(1.27, y.max() + 0.9, r"$p_{\rm AD}$", transform=ax.get_yaxis_transform(), va="bottom", fontsize=6)
    ax.set_ylim(y.min() - 0.8, y.max() + 1.6)
    ax.set_xlim(-1, 1)
    ax.set_xlabel(r"Spearman $r_S$ with $q$ (or $\cos i$)")
    ax.text(-0.98, y[0] + 0.75, r"DM$_{\rm host}$", fontsize=6, color="0.3", va="center")
    ax.text(-0.98, y[n_dm] + 0.75, "Other properties", fontsize=6, color="0.3", va="center")
    ax.tick_params(axis="y", length=0)
    ax.text(0.03, 0.03, "(b)", transform=ax.transAxes, fontsize=9)


def main() -> None:
    h = R.load()
    rob = pd.read_csv(R.OUT / "robust_correlations.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.0), gridspec_kw={"width_ratios": [1.0, 1.0], "wspace": 0.95})
    panel_a(axes[0], h)
    panel_b(axes[1], rob)
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig3_correlations.{ext}", dpi=250, bbox_inches="tight")
    plt.close(fig)
    print("wrote fig3_correlations")


if __name__ == "__main__":
    main()
