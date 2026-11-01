"""Figure 3: host inclination versus DM_host and other host and burst properties.

(a) Rest-frame DM_host versus GALFIT axis ratio q, hosts shared with Marnoch et al.
    (2026) outlined, their best-fitting relation overlaid.
(b) Spearman r_S against cos i with 95% Fisher intervals for burst and host
    properties; the vertical line marks the r_S = -0.44 of Marnoch et al. (2026).

Inputs: plots/plots_null/v2/host_correlations/{host_props.csv, correlations.csv},
        Publication/materials/marnoch_recheck/marnoch2026_table1.csv
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

CORR = P.REPO_ROOT / "plots" / "plots_null" / "v2" / "host_correlations"
MARNOCH = P.REPO_ROOT / "Publication" / "materials" / "marnoch_recheck" / "marnoch2026_table1.csv"
OUT = P.REPO_ROOT / "Publication" / "materials" / "figures"
R_MARNOCH = -0.44

BURST = [("dm_host_rest", r"DM$_{\rm host}$"), ("log_abs_rm", r"$|{\rm RM}|$"),
         ("abs_rm_host_rest", r"$|{\rm RM}_{\rm host}|$"), ("b_par_host", r"$B_{\parallel,{\rm host}}$"),
         ("log_tau", r"Scattering time"), ("lin_pol", "Linear polarization"),
         ("circ_pol", "Circular polarization"), ("width_ms", "Width"), ("fluence_jy_ms", "Fluence"),
         ("snr", "S/N"), ("scint_bw", "Scintillation bandwidth")]
HOST = [("av_old", r"$A_V$, old stars"), ("av", r"$A_V$, young stars"), ("sersic_n", "Sérsic $n$"),
        ("re_kpc", r"$R_e$"), ("logmstar", r"$M_\star$"), ("sfr", "SFR"), ("log_ssfr", "sSFR"),
        ("offset_normalized", r"Offset$/R_e$"), ("u_r_color", r"$u - r$")]


def fisher_ci(r: float, n: int) -> tuple[float, float]:
    z, s = np.arctanh(r), 1.0 / np.sqrt(n - 3)
    return float(np.tanh(z - 1.96 * s)), float(np.tanh(z + 1.96 * s))


def panel_a(ax, hp: pd.DataFrame, shared: set[str]) -> None:
    d = hp[["frb", "q", "dm_host_rest"]].dropna()
    r, p = stats.spearmanr(d.q, d.dm_host_rest)
    sh = d.frb.isin(shared)
    ax.scatter(d.q[~sh], d.dm_host_rest[~sh], s=16, color=P.FRB_C, lw=0, label="This work", zorder=3)
    ax.scatter(d.q[sh], d.dm_host_rest[sh], s=16, color=P.FRB_C, edgecolor="k", lw=0.7,
               label=f"Shared with Marnoch et al. ({sh.sum()})", zorder=4)
    x = np.linspace(0.2, 1.0, 50)
    ax.plot(x, -309.3 * (x - 0.5) + 155.9, "k--", lw=0.9, label="Marnoch et al. (2026) fit", zorder=2)
    ax.axhline(0, color=P.GREY, lw=0.4, zorder=1)
    ax.set_xlim(0.2, 1.0)
    ax.set_xlabel(r"Axis ratio $q$")
    ax.set_ylabel(r"DM$_{\rm host}$ (pc cm$^{-3}$, rest frame)")
    ax.text(0.04, 0.78, f"$N = {len(d)}$\n$r_S = {r:.2f}$, $p = {p:.2f}$", transform=ax.transAxes,
            ha="left", va="top", fontsize=7)
    ax.legend(loc="upper left", frameon=False, fontsize=6.5)
    ax.text(0.03, 0.03, "(a)", transform=ax.transAxes, fontsize=9)


def panel_b(ax, corr: pd.DataFrame) -> None:
    c = corr.set_index("property")
    rows = [("burst", k, lab) for k, lab in BURST] + [("host", k, lab) for k, lab in HOST]
    y = np.arange(len(rows))[::-1].astype(float)
    y[len(BURST):] -= 0.8
    for yy, (grp, k, lab) in zip(y, rows):
        r, n, p = float(c.loc[k, "spearman"]), int(c.loc[k, "N"]), float(c.loc[k, "spearman_p"])
        lo, hi = fisher_ci(r, n)
        col = P.FRB_C if grp == "burst" else P.HSC_C
        ax.errorbar(r, yy, xerr=[[r - lo], [hi - r]], fmt="o", ms=3.5, color=col,
                    mfc=col if p < 0.05 else "white", lw=0.9, capsize=1.5, zorder=3)
        ax.text(1.03, yy, f"{n}", transform=ax.get_yaxis_transform(), va="center", fontsize=6)
    ax.axvline(0, color=P.GREY, lw=0.5, zorder=1)
    ax.axvline(R_MARNOCH, color="k", ls="--", lw=0.9, zorder=2)
    ax.text(R_MARNOCH, y.max() + 1.0, "Marnoch et al.", ha="center", va="bottom", fontsize=6)
    ax.set_yticks(y)
    ax.set_yticklabels([lab for _, _, lab in rows], fontsize=6.5)
    ax.set_ylim(y.min() - 0.8, y.max() + 1.8)
    ax.set_xlim(-1, 1)
    ax.set_xlabel(r"Spearman $r_S$ with $\cos i$")
    ax.text(1.03, y.max() + 1.0, "$N$", transform=ax.get_yaxis_transform(), va="bottom", fontsize=6)
    ax.tick_params(axis="y", length=0)
    ax.text(0.03, 0.03, "(b)", transform=ax.transAxes, fontsize=9)


def main() -> None:
    hp = pd.read_csv(CORR / "host_props.csv", dtype={"frb": str})
    corr = pd.read_csv(CORR / "correlations.csv")
    shared = set(pd.read_csv(MARNOCH, dtype={"frb": str}).frb) & set(hp.frb)
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.6), gridspec_kw={"width_ratios": [1.1, 1.0], "wspace": 0.55})
    panel_a(axes[0], hp, shared)
    panel_b(axes[1], corr)
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig3_correlations.{ext}", dpi=250, bbox_inches="tight")
    plt.close(fig)
    print("wrote fig3_correlations; shared hosts:", len(shared))


if __name__ == "__main__":
    main()
