"""
Legacy Survey DR10 v2 EXP catalog — magnitude vs axis-ratio audit.

Pool: full ``catalog/LS_catalog_v2_fullsky_exp.csv`` (Tractor ``type=EXP`` only).
Axis ratio is recomputed from Tractor ellipticity via

    |e| = sqrt(e1^2 + e2^2),   b/a = (1 - |e|) / (1 + |e|)

(no SDSS ``expAB_r`` column; no Hubble transform).

Outputs under ``plots/plots_null/v2/ls_audit/``:
  - mag_histogram_ba.png / .csv
  - median_ba_vs_mag.png / .csv
  - ba_vs_mag_scatter.png
  - README.md

Scaled pool (same face-on cap as CDF ``strict_scaled``):
  keep ``b/a ≤ 0.8`` (including ``b/a < q0``); plot ``b/a' = (b/a)/0.8``
  so the face-on edge maps to 1. Strict ``b/a > q0`` remains CDF-only.
  → ``plots/plots_null/v2/ls_audit/scaled/``

Run from repo root::

    python scripts/audit_ls_v2_mag_ba.py
    python scripts/audit_ls_v2_mag_ba.py --scaled
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

from null_catalog_utils import LS_CATALOG_V2_EXP_DEFAULT, Q0  # noqa: E402
from pipeline_null_plot_utils import REPO_ROOT  # noqa: E402

OUT_ROOT = REPO_ROOT / "plots" / "plots_null" / "v2" / "ls_audit"
DEFAULT_CSV = REPO_ROOT / LS_CATALOG_V2_EXP_DEFAULT

USECOLS = (
    "ls_id",
    "modelMag_r",
    "shape_e1",
    "shape_e2",
    "RA_ICRS",
    "DE_ICRS",
)

MAG_CLIP = 15.0
BIN_WIDTH = 0.5
# Display cut: bins fainter than this stay in CSV but not plots.
PLOT_MAG_MAX = 26.0
MIN_BIN_N = 1
# Scatter: hexbin of all finite rows; optional thin overlay subsample for markers.
SCATTER_OVERLAY_N = 40_000
BA_FACE_CAP = 0.8  # scaled mode: drop face-on vicinity (matches CDF strict_scaled)


def q_from_e1e2(e1: np.ndarray, e2: np.ndarray) -> np.ndarray:
    """Tractor ellipticity → axis ratio b/a = (1-|e|)/(1+|e|)."""
    e1a = np.asarray(e1, dtype=float)
    e2a = np.asarray(e2, dtype=float)
    eabs = np.hypot(e1a, e2a)
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


def mag_bin_labels(edges: np.ndarray) -> list[str]:
    labels: list[str] = []
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        if not np.isfinite(lo):
            labels.append(f"< {hi:g}")
        else:
            labels.append(f"{lo:g}–{hi:g}")
    return labels


def bin_mask(mag: np.ndarray, lo: float, hi: float, *, last: bool) -> np.ndarray:
    if not np.isfinite(lo):
        return np.isfinite(mag) & (mag <= hi)
    if last:
        return np.isfinite(mag) & (mag > lo) & (mag <= hi)
    return np.isfinite(mag) & (mag > lo) & (mag <= hi)


def load_pool(path: Path) -> pd.DataFrame:
    print(f"[*] Loading {path} ...", flush=True)
    df = pd.read_csv(path, usecols=list(USECOLS))
    e1 = pd.to_numeric(df["shape_e1"], errors="coerce").to_numpy(dtype=float)
    e2 = pd.to_numeric(df["shape_e2"], errors="coerce").to_numpy(dtype=float)
    mag = pd.to_numeric(df["modelMag_r"], errors="coerce").to_numpy(dtype=float)
    ba = q_from_e1e2(e1, e2)
    ok = (
        np.isfinite(mag)
        & np.isfinite(ba)
        & (ba >= 0.0)
        & (ba <= 1.0)
        & np.isfinite(e1)
        & np.isfinite(e2)
    )
    out = pd.DataFrame(
        {
            "ls_id": df.loc[ok, "ls_id"].to_numpy(),
            "modelMag_r": mag[ok],
            "ba": ba[ok],
            "shape_e1": e1[ok],
            "shape_e2": e2[ok],
            "RA_ICRS": pd.to_numeric(df.loc[ok, "RA_ICRS"], errors="coerce").to_numpy(),
            "DE_ICRS": pd.to_numeric(df.loc[ok, "DE_ICRS"], errors="coerce").to_numpy(),
        }
    )
    print(
        f"[*] Finite mag + b/a(e1,e2): N={len(out):,}  "
        f"(dropped {(~ok).sum():,})",
        flush=True,
    )
    return out


def make_bin_summary(mag: np.ndarray, ba: np.ndarray, edges: np.ndarray) -> pd.DataFrame:
    labels = mag_bin_labels(edges)
    rows: list[dict] = []
    n_tot = len(mag)
    for b in range(len(edges) - 1):
        lo, hi = edges[b], edges[b + 1]
        mask = bin_mask(mag, lo, hi, last=(b == len(edges) - 2))
        n = int(np.sum(mask))
        rows.append(
            {
                "bin_label": labels[b],
                "mag_lo": float(lo) if np.isfinite(lo) else np.nan,
                "mag_hi": float(hi),
                "n_galaxies": n,
                "frac_of_pool_pct": round(100.0 * n / max(1, n_tot), 4),
                "median_ba": round(float(np.median(ba[mask])), 4) if n else np.nan,
                "mean_ba": round(float(np.mean(ba[mask])), 4) if n else np.nan,
            }
        )
    return pd.DataFrame(rows)


def plot_display_summary(summary: pd.DataFrame, *, mag_max: float, min_n: int) -> pd.DataFrame:
    """Bins kept for histogram / median plots (CSV retains the full table)."""
    hi_ok = summary["mag_hi"] <= mag_max + 1e-9
    n_ok = summary["n_galaxies"] >= min_n
    return summary.loc[hi_ok & n_ok].copy()


def plot_mag_histogram(
    summary: pd.DataFrame,
    n_pool: int,
    out_png: Path,
    *,
    subtitle: str,
) -> None:
    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(len(summary))
    counts = summary["n_galaxies"].to_numpy()
    bars = ax.bar(x, counts, color="#377eb8", edgecolor="white", linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(summary["bin_label"], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Number of galaxies")
    ax.set_xlabel(r"modelMag$_r$ bin")
    ax.set_title(
        f"Legacy v2 EXP catalog vs magnitude\nN={n_pool:,}  |  {subtitle}"
        if subtitle
        else "Legacy v2 EXP (scaled)"
    )
    for bar, n in zip(bars, counts):
        if n > 0:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height(),
                f"{n:,}",
                ha="center",
                va="bottom",
                fontsize=5.5,
            )
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_median_ba_vs_mag(
    summary: pd.DataFrame,
    n_pool: int,
    out_png: Path,
    *,
    subtitle: str,
    ba_ylim: tuple[float, float] = (0.0, 1.0),
    ylabel: str = r"median $b/a$ from $e_1,e_2$",
) -> None:
    plot = summary.loc[summary["n_galaxies"] > 0].copy()
    centers = []
    for _, r in plot.iterrows():
        if not np.isfinite(r["mag_lo"]):
            centers.append(r["mag_hi"] - 0.25)
        else:
            centers.append(0.5 * (r["mag_lo"] + r["mag_hi"]))
    centers = np.asarray(centers, dtype=float)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(
        centers,
        plot["median_ba"],
        "o-",
        color="#377eb8",
        markersize=5,
        linewidth=1.5,
        label=f"median b/a (N={n_pool:,})",
    )
    ax.set_xlabel(r"modelMag$_r$ (bin center)")
    ax.set_ylabel(ylabel)
    ax.set_title(
        "Legacy v2 EXP: median axis ratio vs magnitude\n"
        rf"$b/a=(1-|e|)/(1+|e|)$  |  {subtitle}"
        if subtitle
        else "Legacy v2 EXP (scaled) — median $(b/a)'$"
    )
    ax.set_ylim(*ba_ylim)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_ba_vs_mag_scatter(
    mag: np.ndarray,
    ba: np.ndarray,
    n_pool: int,
    out_png: Path,
    *,
    overlay_n: int,
    mag_max: float,
    subtitle: str,
    ba_ylim: tuple[float, float] = (0.0, 1.0),
    ylabel: str = r"$b/a$ from $e_1,e_2$",
    seed: int = 42,
) -> None:
    # Restrict hexbin x-range to the well-populated faint limit (CSV uses full sample).
    in_range = np.isfinite(mag) & (mag <= mag_max)
    mag_p = mag[in_range]
    ba_p = ba[in_range]

    fig, ax = plt.subplots(figsize=(10, 6))
    hb = ax.hexbin(
        mag_p,
        ba_p,
        gridsize=80,
        bins="log",
        cmap="viridis",
        mincnt=1,
        extent=(
            float(np.nanpercentile(mag_p, 0.1)),
            float(mag_max),
            float(ba_ylim[0]),
            float(ba_ylim[1]),
        ),
    )
    cb = fig.colorbar(hb, ax=ax, pad=0.02)
    cb.set_label(r"log$_{10}$(N) per hex")

    rng = np.random.default_rng(seed)
    if len(mag_p) > overlay_n:
        idx = rng.choice(len(mag_p), size=overlay_n, replace=False)
        ax.scatter(
            mag_p[idx],
            ba_p[idx],
            s=2,
            c="white",
            alpha=0.08,
            linewidths=0,
            rasterized=True,
            zorder=2,
        )

    ax.set_xlabel(r"modelMag$_r$")
    ax.set_ylabel(ylabel)
    ax.set_ylim(*ba_ylim)
    ax.set_xlim(right=mag_max)
    ax.set_title(
        f"Legacy v2 EXP: $b/a$ vs magnitude (hexbin + {overlay_n:,} overlay pts)\n"
        f"N={n_pool:,} shown for modelMag$_r$ ≤ {mag_max:g}  |  {subtitle}"
        if subtitle
        else "Legacy v2 EXP (scaled)"
    )
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_readme(
    out_root: Path,
    *,
    n_pool: int,
    ra_span: tuple[float, float],
    dec_span: tuple[float, float],
    median_ba: float,
    mean_ba: float,
    mag_med: float,
    scaled: bool,
) -> None:
    if scaled:
        text = f"""# Legacy Survey v2 EXP — scaled b/a audit

Catalog: `catalog/LS_catalog_v2_fullsky_exp.csv` (Tractor `type=EXP`).

Face-on cap matches CDF `strict_scaled`: keep **b/a ≤ {BA_FACE_CAP:g}**
(including b/a < q0={Q0:g}; strict q0 cut is CDF-only). Plot the stretched
axis ratio

$$(b/a)' = (b/a) / {BA_FACE_CAP:g}$$

so the former face-on edge at b/a={BA_FACE_CAP:g} maps to 1.

| Metric | Value |
|--------|-------|
| Rows (finite mag + b/a ≤ {BA_FACE_CAP:g}) | {n_pool:,} |
| RA span | {ra_span[0]:.3f} … {ra_span[1]:.3f} |
| Dec span | {dec_span[0]:.3f} … {dec_span[1]:.3f} |
| median modelMag_r | {mag_med:.3f} |
| median (b/a)' | {median_ba:.4f} |
| mean (b/a)' | {mean_ba:.4f} |

## Plots

| File | Content |
|------|---------|
| `mag_histogram_ba.png` | Counts per 0.5 mag bin |
| `median_ba_vs_mag.png` | Median (b/a)' per mag bin |
| `ba_vs_mag_scatter.png` | (b/a)' vs mag hexbin (y→1 at former ba={BA_FACE_CAP:g}) |

## Regenerate

```bash
python scripts/audit_ls_v2_mag_ba.py --scaled
```
"""
    else:
        text = f"""# Legacy Survey v2 EXP catalog audit

Catalog: `catalog/LS_catalog_v2_fullsky_exp.csv` (Tractor `type=EXP` only, full LS footprint).

| Metric | Value |
|--------|-------|
| Rows (finite mag + b/a) | {n_pool:,} |
| RA span | {ra_span[0]:.3f} … {ra_span[1]:.3f} |
| Dec span | {dec_span[0]:.3f} … {dec_span[1]:.3f} |
| median modelMag_r | {mag_med:.3f} |
| median b/a | {median_ba:.4f} |
| mean b/a | {mean_ba:.4f} |

## Axis ratio

All galaxies are already EXP, so shape is taken from Tractor ellipticity:

$$|e| = \\sqrt{{e_1^2 + e_2^2}},\\qquad b/a = (1 - |e|) / (1 + |e|)$$

No Hubble transform; no color or `b/a > q_0` cuts in this audit.

## Plots

| File | Content |
|------|---------|
| `mag_histogram_ba.png` | Galaxy counts per 0.5 mag bin (`modelMag_r`; bright open bin `<15`; display through mag ≤ 26) |
| `median_ba_vs_mag.png` | Median b/a per mag bin (same display cut) |
| `ba_vs_mag_scatter.png` | b/a vs mag (log hexbin + thin scatter overlay; mag ≤ 26) |

CSV companions retain the full magnitude range: `mag_histogram_ba.csv`, `median_ba_vs_mag.csv`.

Scaled pool (same ba cut as CDF `strict_scaled`): see `scaled/`.

```bash
python scripts/audit_ls_v2_mag_ba.py --scaled
```

## Mag-cut cos(i) CDFs (`cdfs/`)

Hubble cos(i) from b/a(e1,e2) with q0=0.2 and b/a > q0. No color cut.

| File | Cut |
|------|-----|
| `mag20.png` | modelMag_r ≤ 20 |
| `mag21.png` | modelMag_r ≤ 21 |
| `mag22.png` | modelMag_r ≤ 22 |
| `overlay.png` | all three overlaid |
| `summary.csv` | N + median cos(i) |

```bash
python scripts/plot_ls_v2_mag_cut_cdfs.py
python scripts/plot_ls_v2_mag_cut_cdfs.py --mode strict_scaled
```

## Regenerate

```bash
python scripts/audit_ls_v2_mag_ba.py
```
"""
    (out_root / "README.md").write_text(text, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--out-root", type=Path, default=None)
    parser.add_argument("--mag-clip", type=float, default=MAG_CLIP)
    parser.add_argument("--bin-width", type=float, default=BIN_WIDTH)
    parser.add_argument("--overlay-n", type=int, default=SCATTER_OVERLAY_N)
    parser.add_argument(
        "--scaled",
        action="store_true",
        help=f"Keep b/a ≤ {BA_FACE_CAP:g} (incl. b/a < q0), plot (b/a)/{BA_FACE_CAP:g} "
        f"so face-on edge maps to 1; write under v2/ls_audit/scaled/",
    )
    args = parser.parse_args()

    out_root = args.out_root
    if out_root is None:
        out_root = OUT_ROOT / "scaled" if args.scaled else OUT_ROOT

    pool = load_pool(args.csv)
    if args.scaled:
        before = len(pool)
        # Face-on cap only: keep ba < q0 in scatter (strict q0 is CDF-only).
        keep = pool["ba"] <= BA_FACE_CAP
        pool = pool.loc[keep].reset_index(drop=True)
        # Stretch so former ba=BA_FACE_CAP maps to 1.
        pool = pool.copy()
        pool["ba"] = pool["ba"].to_numpy(dtype=float) / BA_FACE_CAP
        n_below_q0 = int(np.sum(pool["ba"] * BA_FACE_CAP <= Q0))
        print(
            f"[*] scaled: kept b/a <= {BA_FACE_CAP:g}: "
            f"{len(pool):,} / {before:,} (dropped {before - len(pool):,} with ba>{BA_FACE_CAP:g}); "
            f"of which native ba<=q0={Q0:g}: {n_below_q0:,}",
            flush=True,
        )
        print(
            f"[*] plotted ba' = (b/a)/{BA_FACE_CAP:g}  "
            f"(median ba'={float(pool['ba'].median()):.4f})",
            flush=True,
        )
        subtitle = ""  # scaled plots: short title only; details live in axis labels / README
        ba_ylim = (0.0, 1.0)
        ylabel = r"$(b/a)/0.8$"
    else:
        subtitle = r"type=EXP  |  b/a from e1,e2  |  no color / q cuts"
        ba_ylim = (0.0, 1.0)
        ylabel = r"$b/a$ from $e_1,e_2$"

    mag = pool["modelMag_r"].to_numpy(dtype=float)
    ba = pool["ba"].to_numpy(dtype=float)
    n_pool = len(pool)

    edges = mag_bin_edges(mag, clip=args.mag_clip, width=args.bin_width)
    summary = make_bin_summary(mag, ba, edges)

    out_root.mkdir(parents=True, exist_ok=True)

    hist_csv = out_root / "mag_histogram_ba.csv"
    med_csv = out_root / "median_ba_vs_mag.csv"
    summary.to_csv(hist_csv, index=False)
    summary.to_csv(med_csv, index=False)

    display = plot_display_summary(
        summary, mag_max=PLOT_MAG_MAX, min_n=MIN_BIN_N
    )
    plot_mag_histogram(
        display, n_pool, out_root / "mag_histogram_ba.png", subtitle=subtitle
    )
    plot_median_ba_vs_mag(
        display,
        n_pool,
        out_root / "median_ba_vs_mag.png",
        subtitle=subtitle,
        ba_ylim=ba_ylim,
        ylabel=f"median {ylabel}",
    )
    plot_ba_vs_mag_scatter(
        mag,
        ba,
        n_pool,
        out_root / "ba_vs_mag_scatter.png",
        overlay_n=args.overlay_n,
        mag_max=PLOT_MAG_MAX,
        subtitle=subtitle,
        ba_ylim=ba_ylim,
        ylabel=ylabel,
    )

    write_readme(
        out_root,
        n_pool=n_pool,
        ra_span=(float(pool["RA_ICRS"].min()), float(pool["RA_ICRS"].max())),
        dec_span=(float(pool["DE_ICRS"].min()), float(pool["DE_ICRS"].max())),
        median_ba=round(float(np.median(ba)), 4),
        mean_ba=round(float(np.mean(ba)), 4),
        mag_med=round(float(np.median(mag)), 3),
        scaled=args.scaled,
    )

    print(f"[*] Wrote {out_root}", flush=True)
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
