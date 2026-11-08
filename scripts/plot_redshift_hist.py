"""Redshift distribution of the FRB hosts.

Input:  plots/plots_null/v2/host_correlations/host_props.csv,
        Publication/materials/host_properties/host_properties_57.csv (z_type)
Output: plots/plots_null/v2/host_correlations/redshift_hist.{pdf,png}
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_paper_inclination_figures as P  # noqa: E402  (style)

HC = P.REPO_ROOT / "plots" / "plots_null" / "v2" / "host_correlations"
HP = P.REPO_ROOT / "Publication" / "materials" / "host_properties" / "host_properties_57.csv"


def main() -> None:
    h = pd.read_csv(HC / "host_props.csv", dtype={"frb": str}).dropna(subset=["z"])
    zt = pd.read_csv(HP, dtype={"frb": str})[["frb", "z_type"]]
    h = h.merge(zt, on="frb", how="left")
    phot = h.z_type.astype(str).str.lower().str.startswith("phot")
    bins = np.arange(0, h.z.max() + 0.1, 0.05)
    fig, ax = plt.subplots(figsize=(3.5, 3.0))
    ax.hist([h.z[~phot], h.z[phot]], bins=bins, stacked=True, color=[P.FRB_C, "k"],
            label=[f"Spectroscopic (N = {(~phot).sum()})", f"Photometric (N = {phot.sum()})"])
    ax.axvline(h.z.median(), color="k", ls="--", lw=0.6)
    ax.text(h.z.median() + 0.01, 0.97, f"median $z = {h.z.median():.2f}$", transform=ax.get_xaxis_transform(),
            va="top", fontsize=6.5)
    ax.set_xlim(0, bins[-1])
    ax.set_xlabel("Redshift $z$")
    ax.set_ylabel("Number of FRB hosts")
    ax.legend(loc="upper right", frameon=False, fontsize=6.5, bbox_to_anchor=(1.0, 0.9))
    ax.tick_params(top=False)
    for ext in ("pdf", "png"):
        fig.savefig(HC / f"redshift_hist.{ext}", dpi=250, bbox_inches="tight")
    print(len(h), (~phot).sum(), phot.sum(), h.z.min(), h.z.median(), h.z.max())


if __name__ == "__main__":
    main()
