"""
Side-by-side median b/a vs mag: LS EXP | LS EXP (scaled) | HSC | DES Y1.

Writes ``plots/plots_null/v2/comparisons/median_ba_vs_mag_quad.png``.

LS scaled: keep ``b/a ≤ 0.8``, plot median of ``(b/a)/0.8`` (same as audit).

Run from repo root::

    python scripts/plot_median_ba_vs_mag_quad.py
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

from null_catalog_utils import (  # noqa: E402
    DES_Y1_MORPH_SAMPLE_DEFAULT,
    HSC_KAWIN_SAMPLE_DEFAULT,
    LS_CATALOG_V2_EXP_DEFAULT,
)
from pipeline_null_plot_utils import REPO_ROOT  # noqa: E402

OUT_DIR = REPO_ROOT / "plots" / "plots_null" / "v2" / "comparisons"
MAG_CLIP = 15.0
BIN_WIDTH = 0.5
PLOT_MAG_MAX = 26.0
MIN_BIN_N = 1
BA_FACE_CAP = 0.8


def q_from_e1e2(e1: np.ndarray, e2: np.ndarray) -> np.ndarray:
    eabs = np.hypot(np.asarray(e1, dtype=float), np.asarray(e2, dtype=float))
    q = np.full(eabs.shape, np.nan, dtype=float)
    good = np.isfinite(eabs) & (eabs < 1.0)
    q[good] = (1.0 - eabs[good]) / (1.0 + eabs[good])
    return q


def mag_bin_edges(mag: np.ndarray, *, clip: float, width: float) -> np.ndarray:
    finite = mag[np.isfinite(mag)]
    hi = float(np.ceil(finite.max() / width) * width) if len(finite) else clip + width
    edges: list[float] = [-np.inf, clip]
    m = clip
    while m < hi:
        m += width
        edges.append(m)
    return np.asarray(edges, dtype=float)


def bin_mask(mag: np.ndarray, lo: float, hi: float) -> np.ndarray:
    if not np.isfinite(lo):
        return np.isfinite(mag) & (mag <= hi)
    return np.isfinite(mag) & (mag > lo) & (mag <= hi)


def make_bin_summary(mag: np.ndarray, ba: np.ndarray, edges: np.ndarray) -> pd.DataFrame:
    rows: list[dict] = []
    for b in range(len(edges) - 1):
        lo, hi = edges[b], edges[b + 1]
        mask = bin_mask(mag, lo, hi)
        n = int(np.sum(mask))
        rows.append(
            {
                "mag_lo": float(lo) if np.isfinite(lo) else np.nan,
                "mag_hi": float(hi),
                "n_galaxies": n,
                "median_ba": float(np.median(ba[mask])) if n else np.nan,
            }
        )
    return pd.DataFrame(rows)


def display_summary(summary: pd.DataFrame) -> pd.DataFrame:
    hi_ok = summary["mag_hi"] <= PLOT_MAG_MAX + 1e-9
    n_ok = summary["n_galaxies"] >= MIN_BIN_N
    plot = summary.loc[hi_ok & n_ok & summary["median_ba"].notna()].copy()
    centers = []
    for _, r in plot.iterrows():
        if not np.isfinite(r["mag_lo"]):
            centers.append(r["mag_hi"] - 0.25)
        else:
            centers.append(0.5 * (r["mag_lo"] + r["mag_hi"]))
    plot["mag_center"] = np.asarray(centers, dtype=float)
    return plot


def load_ls(path: Path) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(path, usecols=["modelMag_r", "shape_e1", "shape_e2"])
    mag = pd.to_numeric(df["modelMag_r"], errors="coerce").to_numpy(dtype=float)
    ba = q_from_e1e2(
        pd.to_numeric(df["shape_e1"], errors="coerce").to_numpy(dtype=float),
        pd.to_numeric(df["shape_e2"], errors="coerce").to_numpy(dtype=float),
    )
    ok = np.isfinite(mag) & np.isfinite(ba) & (ba >= 0.0) & (ba <= 1.0)
    return mag[ok], ba[ok]


def load_ba_mag(path: Path, mag_col: str, ba_col: str) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(path, usecols=[mag_col, ba_col])
    mag = pd.to_numeric(df[mag_col], errors="coerce").to_numpy(dtype=float)
    ba = pd.to_numeric(df[ba_col], errors="coerce").to_numpy(dtype=float)
    ok = np.isfinite(mag) & np.isfinite(ba) & (ba >= 0.0) & (ba <= 1.0)
    return mag[ok], ba[ok]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out", type=Path, default=OUT_DIR / "median_ba_vs_mag_quad.png")
    args = p.parse_args()

    print(f"[*] Loading LS ...", flush=True)
    mag_ls, ba_ls = load_ls(REPO_ROOT / LS_CATALOG_V2_EXP_DEFAULT)
    keep = ba_ls <= BA_FACE_CAP
    mag_ls_s, ba_ls_s = mag_ls[keep], ba_ls[keep] / BA_FACE_CAP
    print(f"    LS: N={len(mag_ls):,}  |  LS scaled: N={len(mag_ls_s):,}", flush=True)

    print(f"[*] Loading HSC ...", flush=True)
    mag_hsc, ba_hsc = load_ba_mag(
        REPO_ROOT / HSC_KAWIN_SAMPLE_DEFAULT, "mag", "ba"
    )
    print(f"    HSC: N={len(mag_hsc):,}", flush=True)

    print(f"[*] Loading DES ...", flush=True)
    mag_des, ba_des = load_ba_mag(
        REPO_ROOT / DES_Y1_MORPH_SAMPLE_DEFAULT, "mag_r", "ba_r"
    )
    print(f"    DES: N={len(mag_des):,}", flush=True)

    panels: list[tuple[str, np.ndarray, np.ndarray]] = [
        ("LS EXP", mag_ls, ba_ls),
        ("LS EXP (scaled)", mag_ls_s, ba_ls_s),
        ("HSC", mag_hsc, ba_hsc),
        ("DES Y1", mag_des, ba_des),
    ]

    colors = ["#377eb8", "#4daf4a", "#984ea3", "#e41a1c"]
    fig, axes = plt.subplots(1, 4, figsize=(16, 4.5), sharey=True)

    for ax, (title, mag, ba), color in zip(axes, panels, colors):
        edges = mag_bin_edges(mag, clip=MAG_CLIP, width=BIN_WIDTH)
        plot = display_summary(make_bin_summary(mag, ba, edges))
        ax.plot(
            plot["mag_center"],
            plot["median_ba"],
            "o-",
            color=color,
            markersize=4,
            linewidth=1.4,
        )
        ax.set_xlim(MAG_CLIP - 0.5, PLOT_MAG_MAX)
        ax.set_ylim(0.0, 1.0)
        ax.set_xlabel(r"$m_r$ (bin center)")
        ax.set_title(f"{title}\nN={len(mag):,}", fontsize=11)
        ax.grid(True, alpha=0.25)

    axes[0].set_ylabel(r"median $b/a$  (scaled: $(b/a)/0.8$)")

    fig.suptitle(r"Median $b/a$ vs magnitude", fontsize=13, y=1.03)
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[*] Wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
