"""
Side-by-side b/a vs mag hexbin scatters: LS EXP | DES Y1 | HSC.

Writes ``plots/plots_null/v2/comparisons/ba_vs_mag_scatter_trio.png``.

Run from repo root::

    python scripts/plot_ba_mag_scatter_trio.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (  # noqa: E402
    DES_Y1_MORPH_SAMPLE_DEFAULT,
    HSC_KAWIN_SAMPLE_DEFAULT,
    LS_CATALOG_V2_EXP_DEFAULT,
)
from pipeline_null_plot_utils import REPO_ROOT  # noqa: E402

OUT_DIR = REPO_ROOT / "plots" / "plots_null" / "v2" / "comparisons"
PLOT_MAG_MAX = 26.0
GRIDSIZE = 70


def q_from_e1e2(e1: np.ndarray, e2: np.ndarray) -> np.ndarray:
    eabs = np.hypot(np.asarray(e1, dtype=float), np.asarray(e2, dtype=float))
    q = np.full(eabs.shape, np.nan, dtype=float)
    good = np.isfinite(eabs) & (eabs < 1.0)
    q[good] = (1.0 - eabs[good]) / (1.0 + eabs[good])
    return q


def load_ls(path: Path) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(path, usecols=["modelMag_r", "shape_e1", "shape_e2"])
    mag = pd.to_numeric(df["modelMag_r"], errors="coerce").to_numpy(dtype=float)
    ba = q_from_e1e2(
        pd.to_numeric(df["shape_e1"], errors="coerce").to_numpy(dtype=float),
        pd.to_numeric(df["shape_e2"], errors="coerce").to_numpy(dtype=float),
    )
    ok = np.isfinite(mag) & np.isfinite(ba) & (ba >= 0.0) & (ba <= 1.0) & (mag <= PLOT_MAG_MAX)
    return mag[ok], ba[ok]


def load_ba_mag(path: Path, mag_col: str, ba_col: str) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(path, usecols=[mag_col, ba_col])
    mag = pd.to_numeric(df[mag_col], errors="coerce").to_numpy(dtype=float)
    ba = pd.to_numeric(df[ba_col], errors="coerce").to_numpy(dtype=float)
    ok = np.isfinite(mag) & np.isfinite(ba) & (ba >= 0.0) & (ba <= 1.0) & (mag <= PLOT_MAG_MAX)
    return mag[ok], ba[ok]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=OUT_DIR / "ba_vs_mag_scatter_trio.png")
    args = p.parse_args()

    panels = [
        ("LS EXP", REPO_ROOT / LS_CATALOG_V2_EXP_DEFAULT, "ls"),
        ("DES Y1", REPO_ROOT / DES_Y1_MORPH_SAMPLE_DEFAULT, "des"),
        ("HSC", REPO_ROOT / HSC_KAWIN_SAMPLE_DEFAULT, "hsc"),
    ]

    data: list[tuple[str, np.ndarray, np.ndarray]] = []
    for label, path, kind in panels:
        print(f"[*] Loading {path.name} ...", flush=True)
        if kind == "ls":
            mag, ba = load_ls(path)
        elif kind == "des":
            mag, ba = load_ba_mag(path, "mag_r", "ba_r")
        else:
            mag, ba = load_ba_mag(path, "mag", "ba")
        print(f"    {label}: N={len(mag):,}", flush=True)
        data.append((label, mag, ba))

    # Shared log color scale across panels
    # Estimate max hex count roughly from densest panel via 2d hist
    vmax = 1.0
    for _, mag, ba in data:
        H, _, _ = np.histogram2d(mag, ba, bins=[GRIDSIZE, GRIDSIZE])
        vmax = max(vmax, float(H.max()))
    norm = LogNorm(vmin=1, vmax=max(vmax, 10))

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    last_hb = None
    for ax, (label, mag, ba) in zip(axes, data):
        xmin = float(np.nanpercentile(mag, 0.5))
        last_hb = ax.hexbin(
            mag,
            ba,
            gridsize=GRIDSIZE,
            cmap="viridis",
            norm=norm,
            mincnt=1,
            extent=(xmin, PLOT_MAG_MAX, 0.0, 1.0),
        )
        ax.set_xlim(xmin, PLOT_MAG_MAX)
        ax.set_ylim(0.0, 1.0)
        ax.set_xlabel(r"$m_r$")
        ax.set_title(f"{label}\nN={len(mag):,}")
        ax.grid(True, alpha=0.25)

    axes[0].set_ylabel(r"$b/a$")
    fig.subplots_adjust(right=0.90, wspace=0.12)
    cax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
    cb = fig.colorbar(last_hb, cax=cax)
    cb.set_label(r"$\log_{10}(N)$ per hex")
    fig.suptitle(r"$b/a$ vs magnitude", fontsize=13, y=1.02)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[*] Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
