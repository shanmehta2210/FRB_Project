"""Rest-frame DM_host versus host SFR and versus host stellar mass.

Inputs: plots/plots_null/v2/host_correlations/host_props.csv,
        Publication/materials/host_properties/host_properties_57.csv (published errors, SFR limits)
Output: plots/plots_null/v2/host_correlations/{dm_host_vs_sfr, dm_host_vs_mass}.{pdf,png}
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
HP = P.REPO_ROOT / "Publication" / "materials" / "host_properties" / "host_properties_57.csv"


def load() -> pd.DataFrame:
    h = pd.read_csv(HC / "host_props.csv", dtype={"frb": str})
    e = pd.read_csv(HP, dtype={"frb": str})[["frb", "logmstar_ep", "logmstar_em", "sfr_ep", "sfr_em", "sfr_limit"]]
    h = h.merge(e, on="frb", how="left")
    lim = pd.to_numeric(h.sfr_limit.astype(str).str.replace("<", "", regex=False), errors="coerce")
    h["is_limit"] = lim.notna() & h.sfr.isna()
    h.loc[h.is_limit, "sfr"] = lim[h.is_limit]
    return h.dropna(subset=["dm_host_rest"])


def finish(ax, x, y, stem: str, xlabel: str) -> None:
    rs, ps = stats.spearmanr(x, y)
    ax.axhline(0, color="k", ls="--", lw=0.6)
    ax.set_ylim(top=1050)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(r"DM$_{\rm host}$ (pc cm$^{-3}$, rest frame)")
    ax.legend(loc="upper left", frameon=False, fontsize=6.5)
    ax.text(0.97, 0.97, f"$r_S = {rs:+.2f}$ ($p = {ps:.2f}$)", transform=ax.transAxes,
            ha="right", va="top", fontsize=6.5)
    for ext in ("pdf", "png"):
        ax.figure.savefig(HC / f"{stem}.{ext}", dpi=250, bbox_inches="tight")
    plt.close(ax.figure)
    print(stem, len(x), round(rs, 3), round(ps, 3))


def main() -> None:
    h = load()

    d = h.dropna(subset=["sfr"])
    d = d[d.sfr > 0]
    det, lim = d[~d.is_limit], d[d.is_limit]
    fig, ax = plt.subplots(figsize=(3.5, 3.0))
    xerr = np.vstack([det.sfr_em.abs().fillna(0).clip(upper=det.sfr * 0.999), det.sfr_ep.abs().fillna(0)])
    ax.errorbar(det.sfr, det.dm_host_rest, xerr=xerr, fmt="o", ms=3.5, color=P.FRB_C, elinewidth=0.5,
                capsize=0, label=f"FRB hosts (N = {len(det)})")
    if len(lim):
        ax.errorbar(lim.sfr, lim.dm_host_rest, xerr=0.3 * lim.sfr, xuplims=True, fmt="o", ms=3.5,
                    mfc="white", color=P.FRB_C, elinewidth=0.5, label=f"SFR upper limit (N = {len(lim)})")
    ax.set_xscale("log")
    finish(ax, det.sfr.to_numpy(float), det.dm_host_rest.to_numpy(float), "dm_host_vs_sfr",
           r"SFR ($M_\odot$ yr$^{-1}$)")

    d = h.dropna(subset=["logmstar"])
    fig, ax = plt.subplots(figsize=(3.5, 3.0))
    xerr = np.vstack([d.logmstar_em.abs().fillna(0), d.logmstar_ep.abs().fillna(0)])
    ax.errorbar(d.logmstar, d.dm_host_rest, xerr=xerr, fmt="o", ms=3.5, color=P.FRB_C, elinewidth=0.5,
                capsize=0, label=f"FRB hosts (N = {len(d)})")
    finish(ax, d.logmstar.to_numpy(float), d.dm_host_rest.to_numpy(float), "dm_host_vs_mass",
           r"$\log(M_*/M_\odot)$")


if __name__ == "__main__":
    main()
