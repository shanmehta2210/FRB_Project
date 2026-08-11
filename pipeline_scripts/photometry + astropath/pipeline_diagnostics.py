"""Phase 2 diagnostic plots: separation vs galaxy angular size (PATH geometry)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

# Viability lines: sep = x_max * R_eff (shape_r)
DEFAULT_X_MAX_LINES: tuple[tuple[float, str, str], ...] = (
    (6, "C1", "-"),
    (12, "C2", "--"),
    (24, "C3", ":"),
    (36, "C4", "-."),
)


def add_sep_arcsec(df: pd.DataFrame, frb_coord: SkyCoord) -> pd.DataFrame:
    """Add separation from FRB (arcsec) using ra/dec columns."""
    out = df.copy()
    coords = SkyCoord(ra=out["ra"].to_numpy(), dec=out["dec"].to_numpy(), unit="deg")
    out["sep_arcsec"] = frb_coord.separation(coords).to(u.arcsec).value
    return out


def plot_candidate_geometry(
    posteriors: pd.DataFrame,
    frb_coord: SkyCoord,
    theta_max: float,
    out_dir: Path | str = ".",
    *,
    title_prefix: str = "",
    show_title: bool = True,
    concise_legend: bool = False,
    base_fontsize: float = 11,
    best_objid: int | None = None,
    x_max_lines: tuple[tuple[float, str, str], ...] | None = None,
) -> tuple[Path, Path]:
    """
    Write sep_vs_shape_r.png and sep_vs_x_max_reff.png.

    posteriors must have columns: ra, dec, ang_size, mag; posterior_O optional.
    theta_max is the x_max used for the viability panel (PATH theta_max).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = add_sep_arcsec(posteriors, frb_coord)
    df = df.copy()
    df["x_max_reff_arcsec"] = theta_max * df["ang_size"]
    if x_max_lines is None:
        x_max_lines = DEFAULT_X_MAX_LINES

    if best_objid is None:
        if "posterior_O" in df.columns and df["posterior_O"].notna().any():
            best = df.loc[df["posterior_O"].idxmax()]
        else:
            best = df.loc[df["sep_arcsec"].idxmin()]
    else:
        best = df[df["objid"] == best_objid].iloc[0]

    prefix = f"{title_prefix}: " if title_prefix else ""
    legend_fs = base_fontsize if concise_legend else base_fontsize - 2
    tick_fs = base_fontsize if concise_legend else base_fontsize - 1
    plt.rcParams.update(
        {
            "font.size": base_fontsize,
            "axes.labelsize": base_fontsize,
            "axes.titlesize": base_fontsize + 1,
            "legend.fontsize": legend_fs,
            "xtick.labelsize": tick_fs,
            "ytick.labelsize": tick_fs,
        }
    )

    # --- sep vs ang_size (shape_r / FLUX_RADIUS) ---
    fig, ax = plt.subplots(figsize=(9, 7))
    sc = ax.scatter(
        df["sep_arcsec"],
        df["ang_size"],
        c=df["mag"],
        cmap="viridis_r",
        s=28,
        alpha=0.75,
        edgecolors="k",
        linewidths=0.2,
    )
    cbar = plt.colorbar(sc, ax=ax, label=r"$m_r$ (AB)")
    cbar.ax.tick_params(labelsize=tick_fs)
    best_label = (
        f"best host (objid {int(best['objid'])})"
        if not concise_legend
        else f"M49 (objid {int(best['objid'])})"
    )
    ax.scatter(
        best["sep_arcsec"],
        best["ang_size"],
        s=120,
        facecolors="none",
        edgecolors="red",
        linewidths=2,
        label=best_label,
        zorder=5,
    )

    sep_max = max(float(df["sep_arcsec"].max()) * 1.02, 1.0)
    sep_grid = np.linspace(0, sep_max, 200)
    for xmax, color, ls in x_max_lines:
        line_label = (
            rf"$\mathrm{{sep}}={xmax:g}\,R_{{\mathrm{{eff}}}}$"
            if concise_legend
            else rf"$x_{{\max}}={xmax:g}$ ($\mathrm{{sep}}={xmax:g}\,R_{{\mathrm{{eff}}}}$)"
        )
        ax.plot(
            sep_grid,
            sep_grid / xmax,
            ls=ls,
            color=color,
            lw=1.2,
            alpha=0.7,
            label=line_label,
        )

    ax.set_xlabel("Separation from FRB position (arcsec)")
    ax.set_ylabel(r"$R_{\mathrm{eff}}$ (arcsec)")
    if show_title:
        ax.set_title(f"{prefix}offset vs. galaxy size")
    ax.legend(loc="upper right", fontsize=legend_fs)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    fig.tight_layout()
    out1 = out_dir / "sep_vs_shape_r.png"
    fig.savefig(out1, dpi=150)
    plt.close(fig)

    # --- sep vs x_max * R_eff viability ---
    fig2, ax2 = plt.subplots(figsize=(9, 6))
    viable = df["sep_arcsec"] <= theta_max * df["ang_size"]
    ax2.scatter(
        df.loc[~viable, "sep_arcsec"],
        df.loc[~viable, "x_max_reff_arcsec"],
        c="0.6",
        s=24,
        alpha=0.5,
        label=rf"excluded ($x_{{\max}}$={theta_max:g})",
    )
    ax2.scatter(
        df.loc[viable, "sep_arcsec"],
        df.loc[viable, "x_max_reff_arcsec"],
        c=df.loc[viable, "mag"],
        cmap="viridis_r",
        s=40,
        edgecolors="k",
        linewidths=0.3,
        label=rf"viable ($x_{{\max}}$={theta_max:g})",
    )
    lim = max(float(df["sep_arcsec"].max()), float(df["x_max_reff_arcsec"].max())) * 1.05
    ax2.plot([0, lim], [0, lim], "k--", lw=1, alpha=0.5, label=rf"sep = {theta_max:g}$\times R_{{\mathrm{{eff}}}}$")
    ax2.scatter(
        best["sep_arcsec"],
        best["x_max_reff_arcsec"],
        s=120,
        facecolors="none",
        edgecolors="red",
        linewidths=2,
        zorder=5,
    )
    ax2.set_xlabel("Separation from FRB (arcsec)")
    ax2.set_ylabel(r"$x_{\max} \times R_{\mathrm{eff}}$ (arcsec)")
    if show_title:
        ax2.set_title(f"{prefix}PATH viability at $x_{{\\max}}$ = {theta_max:g}")
    ax2.legend(loc="upper left", fontsize=9)
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, lim)
    ax2.set_ylim(0, lim)
    fig2.tight_layout()
    out2 = out_dir / "sep_vs_x_max_reff.png"
    fig2.savefig(out2, dpi=150)
    plt.close(fig2)

    return out1, out2
