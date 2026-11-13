"""Per-FRB diagnostic of the DM from intervening (foreground) galaxy halos, after the audit.

FRBs are ordered by host redshift. Foreground values are the corrected estimate of
scripts/foreground_audit.py (halo mass from P(Mh | M*), one halo per group, host environment excluded),
which matches the halo-model expectation along random sightlines; the published pipeline value is shown
for comparison. The Macquart mean <DM_cosmic> already contains the average halo gas, so only the
sightline's excess over the random-sightline expectation <DM_fg>(z) is removed from DM_host.
Panels:
  (a) observer-frame DM budget: DM_MW, diffuse cosmic gas <DM_cosmic> - <DM_fg>(z), this sightline's
      foreground halos, and what is left for the host
  (b) foreground-halo DM: corrected and published estimate, and the random-sightline expectation
  (c) rest-frame DM_host: standard (median DM_cosmic) and with the foreground accounted for consistently,
      with the predicted DM of the host's own halo
  (d) what the correction removed from the published estimate
Hosts outside the Legacy Survey footprint have no foreground estimate and are shaded.

Inputs:  Publication/materials/foreground_halos/audit/{budget, sightlines}.csv
         Publication/materials/appendix/dm_host_physics/halo_budget.csv
Output:  Publication/materials/foreground_halos/foreground_diagnostic.{pdf,png}
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import robust_correlations as RC  # noqa: E402

FG = RC.FG
AUDIT = FG / "audit"
HALO = RC.ROOT / "Publication" / "materials" / "appendix" / "dm_host_physics" / "halo_budget.csv"
COL = {"mw": "0.6", "cosmic": "#9ecae1", "fg": "#e6550d", "host": "#31a354", "pipe": "#D62728", "exp": "#3182bd"}

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "stix", "font.serif": ["STIXGeneral"]})


def load() -> pd.DataFrame:
    h = RC.load()
    b = pd.read_csv(AUDIT / "budget.csv", dtype={"frb": str})[
        ["frb", "mean_c", "fg_pipe", "fg_corr", "fg_expected_corr", "host_obs_corr", "dm_host_corr", "sigma_corr"]]
    s = pd.read_csv(AUDIT / "sightlines.csv", dtype={"sightline": str})
    s = s[s.kind == "frb"].rename(columns={"sightline": "frb"})[["frb", "dm", "dm_cond", "cond_central", "hostenv_cond"]]
    halo = pd.read_csv(HALO, dtype={"frb": str})[["frb", "dm_halo_pred"]]
    h = h.merge(b, on="frb", how="left").merge(s, on="frb", how="left").merge(halo, on="frb", how="left")
    h = h[h.dm_host_rest.notna()].sort_values("z").reset_index(drop=True)
    h["in_fp"] = h.fg_corr.notna()
    return h


def shade_missing(ax, h):
    for k in np.flatnonzero(~h.in_fp):
        ax.axvspan(k - 0.5, k + 0.5, color="0.92", zorder=0, lw=0)


def budget(ax, h, x):
    m = h.in_fp
    fg = h.fg_corr.where(m, 0)
    diffuse = np.where(m, h.mean_c - h.fg_expected_corr, h.dm_exgal - h.dm_host_obs)
    host = np.where(m, h.host_obs_corr, h.dm_host_obs)
    base = np.zeros(len(h))
    for v, c, lab in ((h.dm_mw, COL["mw"], r"DM$_{\rm MW}$"),
                      (diffuse, COL["cosmic"], r"diffuse cosmic gas, $\langle$DM$_{\rm cosmic}\rangle - \langle$DM$_{\rm fg}\rangle(z)$"),
                      (fg, COL["fg"], "foreground halos on this sightline")):
        ax.bar(x, v, bottom=base, color=c, width=0.8, label=lab, zorder=2)
        base = base + np.asarray(v)
    ax.bar(x, np.clip(host, 0, None), bottom=base, color=COL["host"], width=0.8, label="left for the host", zorder=2)
    ax.bar(x, np.clip(host, None, 0), bottom=0, color=COL["host"], width=0.8, alpha=0.5, hatch="///", ec="w", lw=0,
           label="host term < 0", zorder=2)
    ax.scatter(x, h.dm_obs, marker="_", s=60, color="k", lw=1.4, label="observed DM", zorder=3)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_ylabel("DM (pc cm$^{-3}$, observer frame)")
    ax.legend(ncol=3, fontsize=7.5, loc="upper left", frameon=False)
    ax.set_title("(a) Where the observed DM comes from (shaded columns: standard split with the median DM$_{\\rm cosmic}$)",
                 fontsize=9.5, loc="left")


def foreground(ax, h, x):
    m = h.in_fp
    ax.errorbar(x[m], h.fg_corr[m], yerr=0.5 * h.fg_corr[m], fmt="o", ms=4, color=COL["fg"], elinewidth=0.8,
                capsize=1.5, label="corrected estimate (50% systematic)", zorder=3)
    ax.scatter(x[m], h.fg_pipe[m], marker="o", s=16, facecolor="none", edgecolor=COL["pipe"], lw=0.8,
               label="published pipeline", zorder=4)
    ax.plot(x[m], h.fg_expected_corr[m], color=COL["exp"], lw=1.2, drawstyle="steps-mid",
            label=r"expected for a random sightline to the same $z$, $\langle$DM$_{\rm fg}\rangle(z)$", zorder=2)
    ax.text(0.5, 0.95, f"median over {m.sum()} FRBs: corrected {h.fg_corr[m].median():.0f}, published "
            f"{h.fg_pipe[m].median():.0f} pc cm$^{{-3}}$;  > 100 pc cm$^{{-3}}$: {(h.fg_corr[m] > 100).sum()} corrected, "
            f"{(h.fg_pipe[m] > 100).sum()} published;  largest share of the observed DM: "
            f"{(h.fg_corr / h.dm_obs)[m].max():.0%} corrected, {(h.fg_pipe / h.dm_obs)[m].max():.0%} published",
            transform=ax.transAxes, ha="center", va="top", fontsize=7.5)
    ax.set_ylabel("DM$_{\\rm fg}$ (pc cm$^{-3}$, observer frame)")
    ax.set_ylim(bottom=-8)
    ax.legend(fontsize=7.5, loc="upper left", frameon=False, bbox_to_anchor=(0, 0.88))
    ax.set_title("(b) DM from intervening galaxy halos", fontsize=9.5, loc="left")


def host(ax, h, x):
    m = h.in_fp
    corr = h.dm_host_corr
    ax.vlines(x[m], np.minimum(h.dm_host_rest, corr)[m], np.maximum(h.dm_host_rest, corr)[m], color=COL["fg"],
              lw=1.2, zorder=2)
    ax.scatter(x, h.dm_host_rest, s=20, color="k", label=r"standard: DM $-$ DM$_{\rm MW}$ $-$ median DM$_{\rm cosmic}$", zorder=3)
    ax.scatter(x[m], corr[m], s=20, facecolor="w", edgecolor=COL["host"], lw=1.2,
               label=r"DM $-$ DM$_{\rm MW}$ $-$ diffuse cosmic gas $-$ foreground halos", zorder=4)
    ax.scatter(x, h.dm_halo_pred, marker="_", s=50, color=COL["exp"], lw=1.4, label="predicted DM of the host's own halo",
               zorder=4)
    ax.axhline(0, color="k", lw=0.6)
    worst = (h.host_obs_corr / h.sigma_corr)[m].min()
    ax.text(0.5, 0.95, f"negative: standard {(h.dm_host_rest[m] < 0).sum()}, with the foreground "
            f"{(corr[m] < 0).sum()} of {m.sum()}; the most negative is {worst:.1f}$\\sigma$ from zero (audit/AUDIT.md)",
            transform=ax.transAxes, ha="center", va="top", fontsize=7.5)
    ax.set_ylabel(r"DM$_{\rm host}$ (pc cm$^{-3}$, rest frame)")
    ax.legend(fontsize=7.5, loc="upper left", frameon=False, bbox_to_anchor=(0, 0.88))
    ax.set_title("(c) Host DM, standard and with the foreground halos accounted for", fontsize=9.5, loc="left")


def removed(ax, h, x):
    m = h.in_fp
    parts = ((h.cond_central, COL["fg"], "kept (corrected estimate)"),
             (h.hostenv_cond, "#756bb1", "host environment (at the host redshift, inside its halo)"),
             (h.dm_cond - h.cond_central - h.hostenv_cond, "#bcbddc", "satellites: second halo for the same group"),
             (h.dm - h.dm_cond, "#fdd0a2", r"halo mass from inverting $M_*(M_h)$ instead of $P(M_h|M_*)$"))
    base = np.zeros(m.sum())
    for v, c, lab in parts:
        ax.bar(x[m], v[m].clip(lower=0), bottom=base, color=c, width=0.8, label=lab, zorder=2)
        base = base + v[m].clip(lower=0).to_numpy()
    ax.set_ylabel("DM$_{\\rm fg}$ (pc cm$^{-3}$)")
    ax.legend(fontsize=7.5, loc="upper left", frameon=False, ncol=2)
    ax.set_title("(d) Breakdown of the published estimate", fontsize=9.5, loc="left")


def main() -> None:
    h = load()
    x = np.arange(len(h))
    fig, axes = plt.subplots(4, 1, figsize=(15, 14), sharex=True, gridspec_kw=dict(height_ratios=[1.2, 1, 1.1, 0.8]))
    for ax in axes:
        shade_missing(ax, h)
        ax.tick_params(labelsize=7.5)
        ax.grid(axis="y", color="0.9", lw=0.6, zorder=0)
    budget(axes[0], h, x)
    foreground(axes[1], h, x)
    host(axes[2], h, x)
    removed(axes[3], h, x)
    axes[-1].set_xticks(x, [f"{f}  ({z:.2f})" for f, z in zip(h.frb, h.z)], rotation=90, fontsize=7)
    axes[-1].set_xlim(-0.7, len(h) - 0.3)
    axes[-1].set_xlabel("FRB (host redshift), ordered by redshift; grey columns: outside the Legacy Survey footprint, "
                        "no foreground estimate", fontsize=8.5)
    fig.suptitle(f"Intervening-halo DM for the {len(h)} FRBs with DM$_{{\\rm host}}$", fontsize=11, y=0.995)
    fig.tight_layout(h_pad=0.6)
    for ext in ("pdf", "png"):
        fig.savefig(FG / f"foreground_diagnostic.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"{len(h)} FRBs, {h.in_fp.sum()} in the footprint")


if __name__ == "__main__":
    main()
