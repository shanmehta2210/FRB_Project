"""Per-FRB diagnostic of the DM from intervening (foreground) galaxy halos.

FRBs are ordered by host redshift. Panels:
  (a) observer-frame DM budget: DM_MW, <DM_cosmic>, foreground halos, and the remaining host term
  (b) foreground-halo DM (observer frame) with its error, the spectroscopic-only lower bound, and the
      mean of the random sightlines 0.75 deg from the host
  (c) rest-frame DM_host before and after removing the foreground halos, with the predicted DM of the
      host's own halo
  (d) probability that the sightline crosses at least one foreground halo, and the number of likely
      crossings (p_fg >= 0.84)
Hosts outside the Legacy Survey footprint have no foreground estimate and are shaded.

Inputs:  Publication/materials/foreground_halos/{foreground_dm_per_host, foreground_dm_controls, sightline_counts}.csv
         Publication/materials/appendix/dm_host_physics/halo_budget.csv
Output:  Publication/materials/foreground_halos/foreground_diagnostic.{pdf,png}
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
import robust_correlations as RC  # noqa: E402

FG = RC.FG
HALO = RC.ROOT / "Publication" / "materials" / "appendix" / "dm_host_physics" / "halo_budget.csv"
COL = {"mw": "0.6", "cosmic": "#9ecae1", "fg": "#e6550d", "host": "#31a354", "rand": "#3182bd", "spec": "k"}

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "stix", "font.serif": ["STIXGeneral"]})


def load() -> pd.DataFrame:
    h = RC.load()
    fg = pd.read_csv(FG / "foreground_dm_per_host.csv", dtype={"frb": str})[
        ["frb", "dm_fg_err", "n_cand"]]
    ctrl = pd.read_csv(FG / "foreground_dm_controls.csv", dtype={"host": str})
    ctrl = ctrl[ctrl.in_ls.astype(str).isin(["dr9", "dr10"])].groupby("host").dm_fg.mean().rename("dm_fg_random")
    sl = pd.read_csv(FG / "sightline_counts.csv", dtype={"frb": str})[["frb", "p_any", "front", "specz"]]
    halo = pd.read_csv(HALO, dtype={"frb": str})[["frb", "dm_halo_pred"]]
    h = (h.merge(fg, on="frb", how="left").merge(sl, on="frb", how="left").merge(halo, on="frb", how="left")
         .merge(ctrl, left_on="frb", right_index=True, how="left"))
    h = h[h.dm_host_rest.notna()].sort_values("z").reset_index(drop=True)
    h["in_fp"] = h.in_ls.astype(str).isin(["dr9", "dr10"])
    h["dm_cosmic"] = h.dm_exgal - h.dm_host_obs
    h["fg_rest"] = h.dm_fg * (1 + h.z)
    return h


def shade_missing(ax, h):
    for k in np.flatnonzero(~h.in_fp):
        ax.axvspan(k - 0.5, k + 0.5, color="0.92", zorder=0, lw=0)


def budget(ax, h, x):
    fg = h.dm_fg.where(h.in_fp, 0).fillna(0)
    host = h.dm_host_obs - fg
    base = np.zeros(len(h))
    for v, c, lab in ((h.dm_mw, COL["mw"], r"DM$_{\rm MW}$"), (h.dm_cosmic, COL["cosmic"], r"$\langle$DM$_{\rm cosmic}\rangle$"),
                      (fg, COL["fg"], "foreground halos")):
        ax.bar(x, v, bottom=base, color=c, width=0.8, label=lab, zorder=2)
        base = base + v.to_numpy()
    pos = host.clip(lower=0)
    ax.bar(x, pos, bottom=base, color=COL["host"], width=0.8, label="host (after foreground)", zorder=2)
    ax.bar(x, host.clip(upper=0), bottom=0, color=COL["host"], width=0.8, alpha=0.5, hatch="///", ec="w", lw=0,
           label="host term < 0", zorder=2)
    ax.scatter(x, h.dm_obs, marker="_", s=60, color="k", lw=1.4, label=r"observed DM", zorder=3)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_ylabel("DM (pc cm$^{-3}$, observer frame)")
    ax.legend(ncol=6, fontsize=7.5, loc="upper left", frameon=False)
    ax.set_title("(a) Where the observed DM comes from", fontsize=9.5, loc="left")


def foreground(ax, h, x):
    m = h.in_fp
    ax.errorbar(x[m], h.dm_fg[m], yerr=h.dm_fg_err[m], fmt="o", ms=4, color=COL["fg"], ecolor=COL["fg"],
                elinewidth=0.8, capsize=1.5, label="foreground halos (photo-$z$ + spec-$z$)", zorder=3)
    ax.scatter(x[m], h.dm_fg_specz_only[m], marker="x", s=18, color=COL["spec"], lw=0.8,
               label="spectroscopic redshifts only (lower bound)", zorder=4)
    ax.scatter(x[m], h.dm_fg_random[m], marker="D", s=14, facecolor="none", edgecolor=COL["rand"], lw=0.9,
               label="random sightlines 0.75° from the host (mean of about 2)", zorder=4)
    med_f, med_r = h.dm_fg[m].median(), h.dm_fg_random[m].median()
    ax.text(0.5, 0.95, f"median over FRBs: own sightline {med_f:.0f}, mean of its random sightlines {med_r:.0f} "
            f"pc cm$^{{-3}}$;  > 100 pc cm$^{{-3}}$: {(h.dm_fg[m] > 100).sum()} of {m.sum()} FRBs\n"
            f"in the host rest frame, DM$_{{\\rm fg}}(1+z)$: median {h.fg_rest[m].median():.0f}, "
            f"> 100 for {(h.fg_rest[m] > 100).sum()} of {m.sum()}",
            transform=ax.transAxes, ha="center", va="top", fontsize=7.5)
    ax.set_ylabel("DM$_{\\rm fg}$ (pc cm$^{-3}$, observer frame)")
    ax.set_ylim(bottom=-8)
    ax.legend(fontsize=7.5, loc="upper left", frameon=False)
    ax.set_title("(b) DM from intervening galaxy halos", fontsize=9.5, loc="left")


def host(ax, h, x):
    m = h.in_fp
    corr = h.dm_host_rest - h.fg_rest
    ax.vlines(x[m], np.minimum(h.dm_host_rest, corr)[m], np.maximum(h.dm_host_rest, corr)[m], color=COL["fg"],
              lw=1.2, zorder=2)
    ax.scatter(x, h.dm_host_rest, s=20, color="k", label=r"DM$_{\rm host}$ (standard)", zorder=3)
    ax.scatter(x[m], corr[m], s=20, facecolor="w", edgecolor=COL["host"], lw=1.2,
               label=r"DM$_{\rm host}$, foreground halos removed", zorder=4)
    ax.scatter(x, h.dm_halo_pred, marker="_", s=50, color=COL["rand"], lw=1.4, label="predicted DM of the host's own halo",
               zorder=4)
    ax.axhline(0, color="k", lw=0.6)
    frac = (h.fg_rest / h.dm_host_rest).where(m & (h.dm_host_rest > 0))
    ax.text(0.5, 0.95, f"median foreground share of DM$_{{\\rm host}}$: {frac.median():.0%};  "
            f"negative after removal: {(corr[m] < 0).sum()} of {m.sum()}",
            transform=ax.transAxes, ha="center", va="top", fontsize=7.5)
    ax.set_ylabel(r"DM$_{\rm host}$ (pc cm$^{-3}$, rest frame)")
    ax.legend(fontsize=7.5, loc="upper left", frameon=False)
    ax.set_title("(c) Host DM before and after removing the foreground halos", fontsize=9.5, loc="left")


def crossings(ax, h, x):
    m = h.in_fp
    ax.bar(x[m], h.p_any[m], width=0.8, color="#fdae6b", zorder=2, label="P(at least one foreground halo crossed)")
    for k in np.flatnonzero(m):
        ax.text(x[k], h.p_any[k] + 0.03, f"{int(h.front[k])}", ha="center", va="bottom", fontsize=6.5)
    ax.set_ylim(0, 1.25)
    ax.set_yticks([0, 0.5, 1])
    ax.set_ylabel("probability")
    ax.legend(handles=[Patch(color="#fdae6b", label="P(at least one foreground halo crossed)"),
                       Patch(color="w", label="numbers: likely crossings with $p_{\\rm fg} \\geq 0.84$")],
              fontsize=7.5, loc="upper left", frameon=False, ncol=2)
    ax.set_title(f"(d) Halo crossings; expected number of sightlines with at least one: {h.p_any[m].sum():.1f} of {m.sum()}",
                 fontsize=9.5, loc="left")


def main() -> None:
    h = load()
    x = np.arange(len(h))
    fig, axes = plt.subplots(4, 1, figsize=(15, 13), sharex=True, gridspec_kw=dict(height_ratios=[1.2, 1, 1.1, 0.6]))
    for ax in axes:
        shade_missing(ax, h)
        ax.tick_params(labelsize=7.5)
        ax.grid(axis="y", color="0.9", lw=0.6, zorder=0)
    budget(axes[0], h, x)
    foreground(axes[1], h, x)
    host(axes[2], h, x)
    crossings(axes[3], h, x)
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
