#!/usr/bin/env python3
"""COSMOS HST vs SDSS b/a–magnitude audit plots (pass 1).

Run from repo root::

    python scripts/plot_cosmos_ba_mag_audit.py
    python scripts/plot_cosmos_ba_mag_audit.py --build
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS))

from cosmos_ba_audit_utils import (  # noqa: E402
    COSMOS_BA_MIN,
    COSMOS_PLOTS,
    HST_BA_COL,
    HST_CSV,
    HST_MAG_COL,
    HST_MAG_MAX_RELIABLE,
    MIN_N_BIN,
    SDSS_BA_COL,
    SDSS_CSV,
    SDSS_MAG_COL,
    mag_bin_table,
    spearman_mag_ba,
)

PYTHON = sys.executable


def plot_joint_panel(
    bins: pd.DataFrame,
    *,
    out_png: Path,
    ylabel: str,
    title: str,
    mag_label: str,
    pool_n: int,
) -> None:
    fig, (ax_ba, ax_n) = plt.subplots(
        2,
        1,
        figsize=(9, 7),
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1], "hspace": 0.08},
    )
    x = 0.5 * (bins["mag_lo"] + bins["mag_hi"])
    ax_ba.plot(x, bins["median_b_a"], "o-", color="#377eb8", lw=2, ms=5)
    ax_ba.set_ylabel(ylabel)
    ax_ba.set_ylim(0, 1.02)
    ax_ba.axhline(COSMOS_BA_MIN, color="0.6", ls=":", lw=1, label=f"b/a floor {COSMOS_BA_MIN:g}")
    ax_ba.legend(loc="upper right", fontsize=8)
    ax_ba.grid(True, alpha=0.3)
    ax_ba.set_title(f"{title}\nN={pool_n:,}")

    ax_n.plot(x, bins["n"], "o-", color="#377eb8", lw=1.5, ms=4)
    ax_n.set_xlabel(f"{mag_label} (0.5 mag bins)")
    ax_n.set_ylabel(r"$N$ per mag bin")
    ax_n.set_yscale("log")
    ax_n.grid(True, alpha=0.3, which="both")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(hspace=0.12)
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _bins_from_df(
    df: pd.DataFrame,
    mag_col: str,
    ba_col: str,
    *,
    mag_max: float | None,
    min_n: int,
) -> pd.DataFrame:
    mag = pd.to_numeric(df[mag_col], errors="coerce").to_numpy(dtype=float)
    ba = pd.to_numeric(df[ba_col], errors="coerce").to_numpy(dtype=float)
    if mag_max is not None:
        keep = np.isfinite(mag) & (mag <= mag_max)
        mag, ba = mag[keep], ba[keep]
        pool_n = int(keep.sum())
    else:
        pool_n = len(df)
    return mag_bin_table(mag, ba, pool_n=pool_n, mag_col_label=mag_col, min_n=min_n)


def plot_overlay(
    hst_bins: pd.DataFrame,
    sdss_bins: pd.DataFrame,
    *,
    out_png: Path,
    n_hst: int,
    n_sdss: int,
    mag_max: float | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(9, 5.5))
    xh = 0.5 * (hst_bins["mag_lo"] + hst_bins["mag_hi"])
    xs = 0.5 * (sdss_bins["mag_lo"] + sdss_bins["mag_hi"])
    ax.plot(
        xh,
        hst_bins["median_b_a"],
        "o-",
        color="#e41a1c",
        lw=2,
        ms=5,
        label=f"HST Zurich GIM2D (N={n_hst:,})",
    )
    ax.plot(
        xs,
        sdss_bins["median_b_a"],
        "s-",
        color="#377eb8",
        lw=2,
        ms=5,
        label=f"SDSS DR17 modelMag$_r$ (N={n_sdss:,})",
    )
    ax.axhline(COSMOS_BA_MIN, color="0.6", ls=":", lw=1)
    ax.set_xlabel("Magnitude (survey-native; 0.5 mag bins)")
    ax.set_ylabel("median projected b/a")
    ax.set_ylim(0, 1.02)
    ax.legend(loc="upper right")
    ax.grid(True, alpha=0.3)
    subtitle = f"strict b/a > {COSMOS_BA_MIN:g}; no colour cut; bands differ (ACS I vs r)"
    if mag_max is not None:
        subtitle = f"mag <= {mag_max:g}; " + subtitle
    ax.set_title("COSMOS ACS footprint: median b/a vs mag\n" + subtitle)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def run_build() -> None:
    for cmd in (
        [PYTHON, str(_SCRIPTS / "build_cosmos_hst_zurich_catalog.py")],
        [PYTHON, str(_SCRIPTS / "build_cosmos_sdss_catalog.py"), "--no-color-cut"],
    ):
        print(f"[*] Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hst-csv", type=Path, default=HST_CSV)
    parser.add_argument("--sdss-csv", type=Path, default=SDSS_CSV)
    parser.add_argument("--out-dir", type=Path, default=COSMOS_PLOTS)
    parser.add_argument("--build", action="store_true", help="Run catalog builders first.")
    parser.add_argument("--min-n", type=int, default=MIN_N_BIN)
    args = parser.parse_args()

    if args.build or not args.sdss_csv.is_file() or not args.hst_csv.is_file():
        run_build()

    sdss = pd.read_csv(args.sdss_csv)
    hst = pd.read_csv(args.hst_csv)
    print(f"[*] HST pool N={len(hst):,}; SDSS pool N={len(sdss):,}")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    s_mag = pd.to_numeric(sdss[SDSS_MAG_COL], errors="coerce").to_numpy(dtype=float)
    s_ba = pd.to_numeric(sdss[SDSS_BA_COL], errors="coerce").to_numpy(dtype=float)

    sdss_bins = mag_bin_table(
        s_mag, s_ba, pool_n=len(sdss), mag_col_label=SDSS_MAG_COL, min_n=args.min_n
    )
    plot_joint_panel(
        sdss_bins,
        out_png=args.out_dir / "sdss_ba_mag_joint_panel.png",
        ylabel="median expAB$_r$",
        title="SDSS DR17 COSMOS — median expAB_r vs modelMag_r (strict)",
        mag_label="modelMag$_r$",
        pool_n=len(sdss),
    )

    h_mag = pd.to_numeric(hst[HST_MAG_COL], errors="coerce").to_numpy(dtype=float)
    h_ba = pd.to_numeric(hst[HST_BA_COL], errors="coerce").to_numpy(dtype=float)
    hst_bins = mag_bin_table(
        h_mag, h_ba, pool_n=len(hst), mag_col_label=HST_MAG_COL, min_n=args.min_n
    )
    plot_joint_panel(
        hst_bins,
        out_png=args.out_dir / "hst_ba_mag_joint_panel.png",
        ylabel="median b/a (1 - ELL_GIM2D)",
        title="Zurich GIM2D — median b/a vs ACS_MAG_AUTO (strict)",
        mag_label="ACS_MAG_AUTO",
        pool_n=len(hst),
    )
    plot_overlay(
        hst_bins,
        sdss_bins,
        out_png=args.out_dir / "hst_sdss_ba_mag_overlay.png",
        n_hst=len(hst),
        n_sdss=len(sdss),
    )

    hst_bins_225 = _bins_from_df(
        hst, HST_MAG_COL, HST_BA_COL, mag_max=HST_MAG_MAX_RELIABLE, min_n=args.min_n
    )
    sdss_bins_225 = _bins_from_df(
        sdss, SDSS_MAG_COL, SDSS_BA_COL, mag_max=HST_MAG_MAX_RELIABLE, min_n=args.min_n
    )
    n_h_225 = int((pd.to_numeric(hst[HST_MAG_COL], errors="coerce") <= HST_MAG_MAX_RELIABLE).sum())
    n_s_225 = int((pd.to_numeric(sdss[SDSS_MAG_COL], errors="coerce") <= HST_MAG_MAX_RELIABLE).sum())
    if len(hst_bins_225) and len(sdss_bins_225):
        plot_overlay(
            hst_bins_225,
            sdss_bins_225,
            out_png=args.out_dir / "hst_sdss_ba_mag_overlay_mag22p5.png",
            n_hst=n_h_225,
            n_sdss=n_s_225,
            mag_max=HST_MAG_MAX_RELIABLE,
        )

    rho_h, p_h = spearman_mag_ba(h_mag, h_ba)
    rho_s, p_s = spearman_mag_ba(s_mag, s_ba)
    print(
        f"[*] Spearman: HST rho={rho_h:.3f}, SDSS rho={rho_s:.3f}\n"
        f"[*] Plots written to {args.out_dir}"
    )


if __name__ == "__main__":
    main()
