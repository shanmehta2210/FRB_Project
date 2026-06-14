#!/usr/bin/env python3
"""
Histogram of SDSS v2 null-pool galaxy counts vs modelMag_r (0.5 mag bins).

Production pool (no mag limit): u-r < 2.3, lnL exp-wins, **expAB_r > 0.2** only
(exponential profile axis ratio — not best_model_ba_r).

Binning: all modelMag_r < 15 in one bin; then 0.5 mag steps.

Outputs under plots/plots_null/v2_sdss_audit/:
  mag_histogram_expAB_r.png
  mag_histogram_expAB_r.csv

Run from repo root::

    python scripts/plot_sdss_v2_mag_histogram.py
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
    Q0,
    SDSS_Q_COLUMN_CDF,
    SDSS_UR_MAX_CDF,
    prepare_null_strict_color_base,
    read_sdss_null_catalog,
)
from pipeline_null_plot_utils import DEFAULT_SDSS_V2, REPO_ROOT  # noqa: E402

OUT_ROOT = REPO_ROOT / "plots" / "plots_null" / "v2_sdss_audit"
MAG_CLIP = 15.0
BIN_WIDTH = 0.5


def mag_bin_edges(mag_values: np.ndarray, *, clip: float, width: float) -> np.ndarray:
    """Edges: (-inf, clip], then (clip, clip+width], ... through max mag."""
    finite = mag_values[np.isfinite(mag_values)]
    hi = float(np.ceil(finite.max() / width) * width) if len(finite) else clip + width
    edges: list[float] = [-np.inf, clip]
    m = clip
    while m < hi:
        m += width
        edges.append(m)
    return np.asarray(edges, dtype=float)


def mag_bin_labels(edges: np.ndarray) -> list[str]:
    labels: list[str] = []
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        if not np.isfinite(lo):
            labels.append(f"< {hi:g}")
        elif i == len(edges) - 2:
            labels.append(f"{lo:g}–{hi:g}")
        else:
            labels.append(f"{lo:g}–{hi:g}")
    return labels


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sdss-csv", type=Path, default=DEFAULT_SDSS_V2)
    parser.add_argument("--out-root", type=Path, default=OUT_ROOT)
    parser.add_argument("--mag-clip", type=float, default=MAG_CLIP)
    parser.add_argument("--bin-width", type=float, default=BIN_WIDTH)
    parser.add_argument("--q0", type=float, default=Q0)
    parser.add_argument("--ur-max", type=float, default=SDSS_UR_MAX_CDF)
    args = parser.parse_args()

    if SDSS_Q_COLUMN_CDF != "expAB_r":
        raise RuntimeError(f"Expected expAB_r for CDF pools, got {SDSS_Q_COLUMN_CDF!r}")

    print(f"[*] Loading {args.sdss_csv} ...")
    raw = read_sdss_null_catalog(args.sdss_csv)
    pool = prepare_null_strict_color_base(
        raw,
        mag_column="modelMag_r",
        q0=args.q0,
        q_column=SDSS_Q_COLUMN_CDF,
        is_legacy=False,
        sdss_ur_max=args.ur_max,
        sdss_exp_winner_only=True,
    )
    if "expAB_r" not in pool.columns:
        raise KeyError("Pool missing expAB_r after strict base preparation")
    if (pool.columns == "best_model_ba_r").any():
        raise RuntimeError("Pool must not use best_model_ba_r for shape cuts")

    mag = pd.to_numeric(pool["modelMag_r"], errors="coerce").to_numpy(dtype=float)
    expab = pd.to_numeric(pool["expAB_r"], errors="coerce").to_numpy(dtype=float)
    n_pool = len(pool)
    print(
        f"[*] Production pool N={n_pool:,} "
        f"(u-r<{args.ur_max:g}, lnL exp-wins, expAB_r>{args.q0:g})"
    )

    edges = mag_bin_edges(mag, clip=args.mag_clip, width=args.bin_width)
    labels = mag_bin_labels(edges)
    counts = np.zeros(len(edges) - 1, dtype=int)
    median_ba = np.full(len(edges) - 1, np.nan, dtype=float)

    for b in range(len(edges) - 1):
        lo, hi = edges[b], edges[b + 1]
        if not np.isfinite(lo):
            mask = np.isfinite(mag) & (mag <= hi)
        elif b < len(edges) - 2:
            mask = np.isfinite(mag) & (mag > lo) & (mag <= hi)
        else:
            mask = np.isfinite(mag) & (mag > lo) & (mag <= hi)
        counts[b] = int(np.sum(mask))
        if counts[b]:
            median_ba[b] = float(np.median(expab[mask]))

    summary = pd.DataFrame(
        {
            "bin_label": labels,
            "mag_lo": [edges[i] if np.isfinite(edges[i]) else np.nan for i in range(len(labels))],
            "mag_hi": edges[1:],
            "n_galaxies": counts,
            "frac_of_pool_pct": 100.0 * counts / max(1, n_pool),
            "median_expAB_r": median_ba,
        }
    )

    args.out_root.mkdir(parents=True, exist_ok=True)
    csv_path = args.out_root / "mag_histogram_expAB_r.csv"
    summary.to_csv(csv_path, index=False)
    print(f"[*] Wrote {csv_path}")

    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(len(labels))
    bars = ax.bar(x, counts, color="#4daf4a", edgecolor="white", linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Number of galaxies")
    ax.set_xlabel(r"modelMag$_r$ bin")
    ax.set_title(
        f"SDSS v2 strict null pool vs magnitude\n"
        f"N={n_pool:,}  |  shape: expAB$_r$ > {args.q0:g} (exp profile, lnL exp-wins)  |  "
        f"u-r < {args.ur_max:g}"
    )
    for bar, n in zip(bars, counts):
        if n > 0:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height(),
                f"{n:,}",
                ha="center",
                va="bottom",
                fontsize=6,
                rotation=0,
            )
    ax.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()
    png_path = args.out_root / "mag_histogram_expAB_r.png"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[*] Saved {png_path}")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
