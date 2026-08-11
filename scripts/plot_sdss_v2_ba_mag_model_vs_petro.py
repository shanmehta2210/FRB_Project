"""
Median expAB_r vs mag for SDSS v2: modelMag_r vs petroMag_r.

Rebuilds the full-catalog ``ba_mag_joint_panel`` style plot using Petrosian r,
then overlays both binning choices so the difference is visible.

Outputs under plots/plots_null/v2/sdss_audit/formal/:
  - ba_mag_joint_panel_petro.png
  - ba_mag_joint_panel_model_vs_petro_overlay.png
  - ba_per_mag_bin_petro.csv

Run from repo root::

    python scripts/plot_sdss_v2_ba_mag_model_vs_petro.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from audit_sdss_v2_cosi_formal import mag_bin_table  # noqa: E402
from pipeline_null_plot_utils import DEFAULT_SDSS_V2, REPO_ROOT  # noqa: E402

OUT_DIR = REPO_ROOT / "plots" / "plots_null" / "v2" / "sdss_audit" / "formal"


def bin_ba_by_mag(mag: np.ndarray, ba: np.ndarray) -> pd.DataFrame:
    ok = np.isfinite(mag) & np.isfinite(ba) & (ba >= 0.0) & (ba <= 1.0)
    return mag_bin_table(
        mag[ok],
        ba[ok],
        median_col="median_expAB_r",
        mean_col="mean_expAB_r",
    )


def plot_overlay(model_bins: pd.DataFrame, petro_bins: pd.DataFrame, out_png: Path) -> None:
    xm = 0.5 * (model_bins["mag_lo"] + model_bins["mag_hi"])
    xp = 0.5 * (petro_bins["mag_lo"] + petro_bins["mag_hi"])

    fig, ax1 = plt.subplots(figsize=(9.5, 5.2))
    ax1.plot(
        xm,
        model_bins["median_expAB_r"],
        "o-",
        color="#377eb8",
        lw=2.0,
        ms=5,
        label="median expAB_r vs modelMag_r",
    )
    ax1.plot(
        xp,
        petro_bins["median_expAB_r"],
        "s--",
        color="#e41a1c",
        lw=2.0,
        ms=5,
        label="median expAB_r vs petroMag_r",
    )
    ax1.set_xlabel(r"$m_r$ (bin center; model or Petrosian)")
    ax1.set_ylabel("median expAB_r")
    ax1.set_ylim(0.0, 1.05)
    ax1.grid(True, alpha=0.3)

    ax2 = ax1.twinx()
    ax2.plot(
        xm,
        model_bins["frac_pool"],
        ":",
        color="#377eb8",
        alpha=0.55,
        lw=1.4,
        label="N(m) frac (modelMag_r)",
    )
    ax2.plot(
        xp,
        petro_bins["frac_pool"],
        ":",
        color="#e41a1c",
        alpha=0.55,
        lw=1.4,
        label="N(m) frac (petroMag_r)",
    )
    ax2.set_ylabel("pool fraction per bin")

    # Combined legend
    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="upper right", fontsize=8)

    ax1.set_title(
        "SDSS v2 full catalog: median expAB_r — modelMag_r bins vs petroMag_r bins"
    )
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sdss-csv", type=Path, default=DEFAULT_SDSS_V2)
    p.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[*] Loading {args.sdss_csv} ...", flush=True)
    df = pd.read_csv(
        args.sdss_csv,
        usecols=["modelMag_r", "petroMag_r", "expAB_r"],
    )
    model = pd.to_numeric(df["modelMag_r"], errors="coerce").to_numpy()
    petro = pd.to_numeric(df["petroMag_r"], errors="coerce").to_numpy()
    ba = pd.to_numeric(df["expAB_r"], errors="coerce").to_numpy()

    model_bins = bin_ba_by_mag(model, ba)
    petro_bins = bin_ba_by_mag(petro, ba)
    petro_bins.insert(0, "pool", "full_catalog_petroMag_r")
    petro_bins.insert(1, "pool_n", int(np.isfinite(petro).sum()))
    petro_bins.to_csv(args.out_dir / "ba_per_mag_bin_petro.csv", index=False)

    # Standalone petro panel (same dual-axis style as ba_mag_joint_panel.png)
    x = 0.5 * (petro_bins["mag_lo"] + petro_bins["mag_hi"])
    fig, ax1 = plt.subplots(figsize=(9, 5))
    ax1.plot(x, petro_bins["median_expAB_r"], "o-", color="C0", label="median expAB_r")
    ax1.set_xlabel("petroMag_r")
    ax1.set_ylabel("median expAB_r", color="C0")
    ax1.tick_params(axis="y", labelcolor="C0")
    ax2 = ax1.twinx()
    ax2.plot(x, petro_bins["frac_pool"], "s--", color="C1", alpha=0.8, label="N(m) fraction")
    ax2.set_ylabel("pool fraction per bin", color="C1")
    ax2.tick_params(axis="y", labelcolor="C1")
    ax1.set_title("median expAB_r vs petroMag_r (full v2 catalog, no cuts)")
    fig.tight_layout()
    fig.savefig(args.out_dir / "ba_mag_joint_panel_petro.png", dpi=200, bbox_inches="tight")
    plt.close(fig)

    plot_overlay(
        model_bins,
        petro_bins,
        args.out_dir / "ba_mag_joint_panel_model_vs_petro_overlay.png",
    )
    print(f"[*] Wrote petro panel + overlay -> {args.out_dir}", flush=True)


if __name__ == "__main__":
    main()
