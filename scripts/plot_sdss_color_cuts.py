#!/usr/bin/env python3
"""
SDSS-only u-r color-cut mag vs b/a plots (streaming; bounded RAM).

Reads the v1 SDSS CSV in small chunks, keeps galaxies passing each u-r cut,
then plots with matplotlib hexbin (same as audit diagnostics — not square bins).

Peak RAM: ~few MB per cut for stored (mag, b/a) float32 arrays of passing rows only,
not the full 500k-row catalog.

Run from repo root:
    python scripts/plot_sdss_color_cuts.py
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
    Q0,
    SDSS_BA_FLOOR_MIN,
    UR_CUTS_DEFAULT,
    ur_cut_folder_tag,
)
from pipeline_null_plot_utils import (  # noqa: E402
    DEFAULT_SDSS,
    PLOTS_NULL,
    save_figure,
)

DIAG_ROOT = PLOTS_NULL / "v1_null_cdf_inclination" / "diagnostics"

CHUNKSIZE = 8_000
SCATTER_CAP = 15_000
SCATTER_SEED = 42
HEXBIN_GRIDSIZE = 50
PLOT_USECOLS = ("modelMag_r", "best_model_ba_r", "u_r")


def _reservoir_add(
    mag_buf: np.ndarray,
    ba_buf: np.ndarray,
    n_stored: int,
    n_seen: int,
    mag_new: np.ndarray,
    ba_new: np.ndarray,
    k: int,
    rng: np.random.Generator,
) -> tuple[int, int]:
    for m, b in zip(mag_new, ba_new):
        n_seen += 1
        if n_stored < k:
            mag_buf[n_stored] = m
            ba_buf[n_stored] = b
            n_stored += 1
        else:
            j = int(rng.integers(0, n_seen))
            if j < k:
                mag_buf[j] = m
                ba_buf[j] = b
    return n_stored, n_seen


def _stream_pass(
    csv_path: Path,
    ur_cuts: tuple[float, ...],
    *,
    chunksize: int = CHUNKSIZE,
    scatter_cap: int = SCATTER_CAP,
    limit_chunks: int | None = None,
) -> dict[float, dict]:
    """
    One streaming read: per u-r cut, collect all passing (mag, ba) for hexbin
    plus a capped reservoir for the scatter panel.
    """
    n_cuts = len(ur_cuts)
    mag_parts: list[list[np.ndarray]] = [[] for _ in ur_cuts]
    ba_parts: list[list[np.ndarray]] = [[] for _ in ur_cuts]
    mag_bufs = [np.empty(scatter_cap, dtype=np.float64) for _ in ur_cuts]
    ba_bufs = [np.empty(scatter_cap, dtype=np.float64) for _ in ur_cuts]
    n_stored = [0] * n_cuts
    n_seen = [0] * n_cuts
    n_pass = dict.fromkeys(ur_cuts, 0)
    rngs = [np.random.default_rng(SCATTER_SEED + i) for i in range(n_cuts)]

    chunk_i = 0
    for chunk in pd.read_csv(csv_path, usecols=list(PLOT_USECOLS), chunksize=chunksize):
        chunk_i += 1
        if limit_chunks is not None and chunk_i > limit_chunks:
            break

        mr = pd.to_numeric(chunk["modelMag_r"], errors="coerce").to_numpy(dtype=np.float64)
        ba = pd.to_numeric(chunk["best_model_ba_r"], errors="coerce").to_numpy(dtype=np.float64)
        ur = pd.to_numeric(chunk["u_r"], errors="coerce").to_numpy(dtype=np.float64)

        ok = np.isfinite(mr) & np.isfinite(ba) & np.isfinite(ur) & (ba >= 0) & (ba <= 1)
        if not ok.any():
            del chunk
            continue

        mr_ok = mr[ok]
        ba_ok = ba[ok]
        ur_ok = ur[ok]

        for i, ur_max in enumerate(ur_cuts):
            sel = ur_ok < ur_max
            if not sel.any():
                continue
            m = mr_ok[sel].astype(np.float32, copy=False)
            b = ba_ok[sel].astype(np.float32, copy=False)
            n_pass[ur_max] += int(len(m))
            mag_parts[i].append(m)
            ba_parts[i].append(b)
            n_stored[i], n_seen[i] = _reservoir_add(
                mag_bufs[i],
                ba_bufs[i],
                n_stored[i],
                n_seen[i],
                m.astype(np.float64),
                b.astype(np.float64),
                scatter_cap,
                rngs[i],
            )

        del chunk, mr, ba, ur, mr_ok, ba_ok, ur_ok

    out: dict[float, dict] = {}
    for i, ur_max in enumerate(ur_cuts):
        mp, bp = mag_parts[i], ba_parts[i]
        if mp:
            mag_all = np.concatenate(mp).astype(np.float64, copy=False)
            ba_all = np.concatenate(bp).astype(np.float64, copy=False)
        else:
            mag_all = np.array([], dtype=np.float64)
            ba_all = np.array([], dtype=np.float64)

        out[ur_max] = {
            "mag_hex": mag_all,
            "ba_hex": ba_all,
            "mag_scatter": mag_bufs[i][: n_stored[i]],
            "ba_scatter": ba_bufs[i][: n_stored[i]],
            "n_pass": n_pass[ur_max],
        }
    return out


def plot_mag_vs_ba_sdss_cut(
    *,
    mag_hex: np.ndarray,
    ba_hex: np.ndarray,
    mag_scatter: np.ndarray,
    ba_scatter: np.ndarray,
    title: str,
    out_stem: Path,
    n_total: int,
) -> None:
    """Two-panel layout matching audit plot_mag_vs_ba (hexbin + scatter)."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax in axes:
        ax.axhspan(
            0,
            SDSS_BA_FLOOR_MIN,
            color="#d62728",
            alpha=0.12,
            zorder=0,
            label=f"PhotoObj $b/a$ floor ($\\leq {SDSS_BA_FLOOR_MIN}$)",
        )

    # mincnt=1: draw only bins with >=1 object; LogNorm(vmin=1): color scale starts at 1 count
    hb = axes[0].hexbin(
        mag_hex,
        ba_hex,
        gridsize=HEXBIN_GRIDSIZE,
        cmap="viridis",
        mincnt=1,
        norm=LogNorm(vmin=1),
        linewidths=0.0,
    )
    axes[0].axhline(Q0, color="#d62728", linewidth=1.5, linestyle="--", label=f"strict cut $q_0={Q0}$")
    axes[0].set_xlabel("modelMag_r")
    axes[0].set_ylabel("best_model_ba_r")
    axes[0].set_ylim(0, 1)
    axes[0].set_title(f"{title} — hex density (N={n_total:,})")
    axes[0].legend(loc="upper right", fontsize=7)
    fig.colorbar(hb, ax=axes[0], label="count")

    n = len(mag_scatter)
    axes[1].scatter(
        mag_scatter,
        ba_scatter,
        s=4,
        alpha=0.2,
        c="#4daf4a",
        linewidths=0,
        rasterized=True,
    )
    axes[1].axhline(Q0, color="#d62728", linewidth=1.5, linestyle="--", label=f"$q_0={Q0}$")
    axes[1].set_xlabel("modelMag_r")
    axes[1].set_ylabel("best_model_ba_r")
    axes[1].set_ylim(0, 1)
    axes[1].legend(loc="upper right", fontsize=7)
    axes[1].set_title(f"{title} — subsample (n={n:,})")

    if n_total > 0 and n > 0:
        n_floor = int((ba_scatter <= SDSS_BA_FLOOR_MIN).sum())
        fig.text(
            0.5,
            0.01,
            f"Red band: SDSS PhotoObj floor at b/a={SDSS_BA_FLOOR_MIN} "
            f"({n_floor:,} in subsample, {100 * n_floor / n:.1f}% of subsample). "
            "Included in SDSS null CDF pools (not excluded).",
            ha="center",
            fontsize=8,
            color="0.35",
        )

    plt.tight_layout(rect=[0, 0.04, 1, 1])
    save_figure(fig, out_stem)
    plt.close(fig)
    del mag_hex, ba_hex, mag_scatter, ba_scatter


def plot_mag_vs_ba_sdss_cuts_combined(
    by_cut: dict[float, dict],
    *,
    ur_cuts: tuple[float, ...],
    out_stem: Path,
    roles: dict[float, str] | None = None,
) -> None:
    """Single figure: hexbin mag vs b/a per u-r cut (same hexbin logic as single-cut plots)."""
    del roles
    n = len(ur_cuts)
    fig, axes = plt.subplots(1, n, figsize=(5.5 * n, 5.5), sharey=True)
    if n == 1:
        axes = [axes]

    last_hb = None
    for ax, ur_max in zip(axes, ur_cuts):
        d = by_cut[ur_max]
        mag_hex = d["mag_hex"]
        ba_hex = d["ba_hex"]
        n_total = d["n_pass"]
        ax.axhspan(
            0,
            SDSS_BA_FLOOR_MIN,
            color="#d62728",
            alpha=0.12,
            zorder=0,
        )
        if len(mag_hex) > 0:
            # Same as plot_mag_vs_ba_sdss_cut hex panel: autoscale per cut, LogNorm(vmin=1).
            last_hb = ax.hexbin(
                mag_hex,
                ba_hex,
                gridsize=HEXBIN_GRIDSIZE,
                cmap="viridis",
                mincnt=1,
                norm=LogNorm(vmin=1),
                linewidths=0.0,
            )
        ax.axhline(Q0, color="#d62728", linewidth=1.5, linestyle="--")
        ax.set_xlabel("modelMag_r")
        ax.set_ylim(0, 1)
        ax.set_title(f"$u-r < {ur_max:g}$\n$N={n_total:,}$", fontsize=10)
    axes[0].set_ylabel("best_model_ba_r")

    fig.subplots_adjust(left=0.07, right=0.84, wspace=0.22, top=0.88, bottom=0.12)
    if last_hb is not None:
        cax = fig.add_axes([0.86, 0.14, 0.025, 0.72])
        fig.colorbar(last_hb, cax=cax, label="count")
    save_figure(fig, out_stem)
    plt.close(fig)


def run_sdss_color_cuts(
    sdss_csv: Path,
    out_root: Path,
    *,
    ur_cuts: tuple[float, ...] = UR_CUTS_DEFAULT,
    limit_chunks: int | None = None,
) -> None:
    out_root.mkdir(parents=True, exist_ok=True)
    print(
        f"[*] Streaming {sdss_csv.name} (chunksize={CHUNKSIZE}, "
        f"hexbin gridsize={HEXBIN_GRIDSIZE}, scatter cap={SCATTER_CAP:,})"
    )
    if limit_chunks:
        print(f"[!] --limit-chunks {limit_chunks} (test mode)")

    by_cut = _stream_pass(sdss_csv, ur_cuts, limit_chunks=limit_chunks)
    counts = {ur: by_cut[ur]["n_pass"] for ur in ur_cuts}
    _write_color_cuts_md(out_root / "COLOR_CUTS.md", ur_cuts, counts)

    plot_mag_vs_ba_sdss_cuts_combined(
        by_cut,
        ur_cuts=ur_cuts,
        out_stem=out_root / "mag_vs_ba_sdss_color_cuts_combined",
    )
    print(
        f"[*] Combined color-cut panel: "
        f"{out_root / 'mag_vs_ba_sdss_color_cuts_combined.png'}"
    )

    rows = []
    for ur_max in ur_cuts:
        tag = ur_cut_folder_tag(ur_max)
        sub = out_root / tag
        sub.mkdir(parents=True, exist_ok=True)
        d = by_cut[ur_max]
        plot_mag_vs_ba_sdss_cut(
            mag_hex=d["mag_hex"],
            ba_hex=d["ba_hex"],
            mag_scatter=d["mag_scatter"],
            ba_scatter=d["ba_scatter"],
            title=f"SDSS null ($u-r < {ur_max:g}$)",
            out_stem=sub / "mag_vs_ba_sdss",
            n_total=d["n_pass"],
        )
        rows.append(
            {
                "ur_max": ur_max,
                "n_sdss": d["n_pass"],
                "n_sdss_scatter": len(d["mag_scatter"]),
                "folder": tag,
            }
        )
        print(
            f"[*] {tag}: N={d['n_pass']:,} hexbin "
            f"(scatter n={len(d['mag_scatter']):,})"
        )
        del by_cut[ur_max]

    pd.DataFrame(rows).to_csv(out_root / "color_cut_counts.csv", index=False)


_ROLES = {
    3.5: "Broad blue / star-forming population",
    2.2: "Near Strateva et al. (2001) u*-r* ~ 2.22 separator",
    1.5: "Strongly blue tail",
}


def _write_color_cuts_md(path: Path, ur_cuts: tuple[float, ...], counts: dict[float, int]) -> None:
    lines = [
        "# Color cuts for null-catalog diagnostics (SDSS only)",
        "",
        "Uses `u_r` from the augmented v1 SDSS catalog.",
        "Density panel: **matplotlib hexbin** (`gridsize=50`, `mincnt=1`, counts ≥1 only).",
        "See **`HEXBIN_AND_SDSS_BANDS_AUDIT.md`** for bin geometry and horizontal bands.",
        "",
        "| u-r cut | Role | N in cut |",
        "|---------|------|----------|",
    ]
    for ur in ur_cuts:
        role = _ROLES.get(ur, "")
        lines.append(f"| < {ur:g} | {role} | {counts.get(ur, 0):,} |")
    lines.extend(["", "Legacy g-r proxy plots: **not generated** (deferred).", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sdss-csv", type=Path, default=DEFAULT_SDSS)
    parser.add_argument("--out-dir", type=Path, default=DIAG_ROOT / "color_cuts")
    parser.add_argument(
        "--limit-chunks",
        type=int,
        default=None,
        help="Process only first N chunks (smoke test; not for production).",
    )
    args = parser.parse_args()

    if not args.sdss_csv.is_file():
        raise SystemExit(f"Missing SDSS catalog: {args.sdss_csv}")

    run_sdss_color_cuts(
        args.sdss_csv,
        args.out_dir,
        limit_chunks=args.limit_chunks,
    )
    print(f"[*] Done: {args.out_dir}")


if __name__ == "__main__":
    main()
