"""Host stellar mass versus inclination for the confirmed sample.

Input:  plots/plots_null/v2/host_correlations/host_props.csv
Output: plots/plots_null/v2/host_correlations/mass_vs_inclination.{pdf,png}
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

HC = P.REPO_ROOT / "plots" / "plots_null" / "v2" / "host_correlations"
Q0 = 0.2
SYS_I = 5.0  # deg, systematic term added in quadrature as in the paper


def main() -> None:
    allh = pd.read_csv(HC / "host_props.csv", dtype={"frb": str}).dropna(subset=["i_deg"])
    errs = pd.read_csv(P.REPO_ROOT / "Publication" / "materials" / "host_properties" / "host_properties_57.csv",
                       dtype={"frb": str})[["frb", "logmstar_ep", "logmstar_em"]]
    allh = allh.merge(errs, on="frb", how="left")
    q = allh.q.to_numpy(float)
    c = np.sqrt((q ** 2 - Q0 ** 2) / (1 - Q0 ** 2))
    didq = np.degrees(q / ((1 - Q0 ** 2) * c * np.sqrt(1 - c ** 2)))
    allh["sig_i"] = np.sqrt((allh.sigma_q_B.to_numpy(float) * didq) ** 2 + SYS_I ** 2)
    d = allh.dropna(subset=["logmstar"])
    nomass = allh[allh.logmstar.isna()]
    rs, ps = stats.spearmanr(d.i_deg, d.logmstar)
    fig, ax = plt.subplots(figsize=(3.5, 3.0))
    yerr = np.vstack([d.logmstar_em.abs().fillna(0), d.logmstar_ep.abs().fillna(0)])
    ax.errorbar(d.i_deg, d.logmstar, xerr=d.sig_i, yerr=yerr, fmt="o", ms=3.5, color=P.FRB_C,
                elinewidth=0.5, capsize=0, alpha=0.9, label=f"FRB hosts (N = {len(d)})")
    ymin = d.logmstar.min() - (yerr[0].max() if len(d) else 0) - 0.3
    ax.errorbar(nomass.i_deg, [ymin + 0.1] * len(nomass), xerr=nomass.sig_i, fmt="v", ms=4, color="k",
                elinewidth=0.5, capsize=0, label=rf"No published $M_*$ (N = {len(nomass)})")
    ax.set_ylim(bottom=ymin)
    ax.axvline(60, color="k", ls="--", lw=0.6)
    ax.set_xlim(0, 90)
    ax.set_xlabel(r"Inclination $i$ (deg)")
    ax.set_ylabel(r"$\log(M_*/M_\odot)$")
    ax.legend(loc="upper left", frameon=False, fontsize=6.5)
    ax.set_ylim(top=d.logmstar.max() + 0.3)
    ax.text(0.97, 0.97, f"$r_S = {rs:+.2f}$ ($p = {ps:.2f}$)", transform=ax.transAxes,
            ha="right", va="top", fontsize=6.5)
    for ext in ("pdf", "png"):
        fig.savefig(HC / f"mass_vs_inclination.{ext}", dpi=250, bbox_inches="tight")
    print(len(d), rs, ps)


if __name__ == "__main__":
    main()
