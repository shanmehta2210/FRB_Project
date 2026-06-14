#!/usr/bin/env python3
"""
Legacy Survey g-r color-cut mag vs Tractor b/a plots (streaming; bounded RAM).

Reads the v1 Legacy CSV in chunks, keeps galaxies passing g-r < gr_max
(plus optional REX exclusion and strict b/a > q0 to match null CDF pools).

Peak RAM: stores only passing (mag, b/a) rows as float32 arrays, not the full catalog.

Run from repo root:
    python scripts/plot_legacy_color_cuts.py
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
    LEGACY_GR_MAX_CDF,
    Q0,
    gr_cut_folder_tag,
)
from pipeline_null_plot_utils import (  # noqa: E402
    DEFAULT_LEGACY,
    PLOTS_NULL,
    save_figure,
)

DIAG_ROOT = PLOTS_NULL / "v1_null_cdf_inclination" / "diagnostics" / "color_cuts"

CHUNKSIZE = 8_000
SCATTER_CAP = 15_000
SCATTER_SEED = 42
HEXBIN_GRIDSIZE = 50
PLOT_USECOLS = ("tractor_mag_r", "rmag", "gmag", "expAB_r", "tractor_type")


def _reservoir_add(
    mag_buf: np.ndarray,
    ba_buf: np.ndarray,
    rex_buf: np.ndarray,
    n_stored: int,
    n_seen: int,
    mag_new: np.ndarray,
    ba_new: np.ndarray,
    rex_new: np.ndarray,
    k: int,
    rng: np.random.Generator,
) -> tuple[int, int]:
    for m, b, r in zip(mag_new, ba_new, rex_new):
        n_seen += 1
        if n_stored < k:
            mag_buf[n_stored] = m
            ba_buf[n_stored] = b
            rex_buf[n_stored] = r
            n_stored += 1
        else:
            j = int(rng.integers(0, n_seen))
            if j < k:
                mag_buf[j] = m
                ba_buf[j] = b
                rex_buf[j] = r
    return n_stored, n_seen


def _stream_legacy(
    csv_path: Path,
    *,
    gr_max: float,
    match_cdf_cuts: bool,
    exclude_rex: bool,
    q0: float,
    chunksize: int = CHUNKSIZE,
    scatter_cap: int = SCATTER_CAP,
    limit_chunks: int | None = None,
) -> dict:
    """Collect (mag, Tractor b/a) for Legacy rows passing g-r and optional CDF cuts."""
    mag_parts: list[np.ndarray] = []
    ba_parts: list[np.ndarray] = []
    mag_buf = np.empty(scatter_cap, dtype=np.float64)
    ba_buf = np.empty(scatter_cap, dtype=np.float64)
    rex_buf = np.empty(scatter_cap, dtype=bool)
    n_stored = 0
    n_seen = 0
    n_pass = 0
    rng = np.random.default_rng(SCATTER_SEED)

    chunk_i = 0
    for chunk in pd.read_csv(csv_path, usecols=list(PLOT_USECOLS), chunksize=chunksize):
        chunk_i += 1
        if limit_chunks is not None and chunk_i > limit_chunks:
            break

        mag = pd.to_numeric(chunk["tractor_mag_r"], errors="coerce").to_numpy(dtype=np.float64)
        if not np.isfinite(mag).any():
            mag = pd.to_numeric(chunk["rmag"], errors="coerce").to_numpy(dtype=np.float64)
        ba = pd.to_numeric(chunk["expAB_r"], errors="coerce").to_numpy(dtype=np.float64)
        gmag = pd.to_numeric(chunk["gmag"], errors="coerce").to_numpy(dtype=np.float64)
        rmag = pd.to_numeric(chunk["rmag"], errors="coerce").to_numpy(dtype=np.float64)
        gr = gmag - rmag
        types = chunk["tractor_type"].astype(str).str.upper().to_numpy()

        ok = (
            np.isfinite(mag)
            & np.isfinite(ba)
            & np.isfinite(gr)
            & (ba >= 0)
            & (ba <= 1)
            & (gr < gr_max)
        )
        if match_cdf_cuts:
            ok &= ba > q0
            if exclude_rex:
                ok &= types != "REX"

        if not ok.any():
            del chunk
            continue

        m = mag[ok].astype(np.float32, copy=False)
        b = ba[ok].astype(np.float32, copy=False)
        rex = types[ok] == "REX"
        n_pass += int(len(m))
        mag_parts.append(m)
        ba_parts.append(b)
        n_stored, n_seen = _reservoir_add(
            mag_buf,
            ba_buf,
            rex_buf,
            n_stored,
            n_seen,
            m.astype(np.float64),
            b.astype(np.float64),
            rex,
            scatter_cap,
            rng,
        )
        del chunk, mag, ba, gmag, rmag, gr, types, m, b, rex

    if mag_parts:
        mag_all = np.concatenate(mag_parts).astype(np.float64, copy=False)
        ba_all = np.concatenate(ba_parts).astype(np.float64, copy=False)
    else:
        mag_all = np.array([], dtype=np.float64)
        ba_all = np.array([], dtype=np.float64)

    return {
        "mag_hex": mag_all,
        "ba_hex": ba_all,
        "mag_scatter": mag_buf[:n_stored],
        "ba_scatter": ba_buf[:n_stored],
        "rex_scatter": rex_buf[:n_stored],
        "n_pass": n_pass,
    }


def plot_mag_vs_ba_legacy(
    *,
    mag_hex: np.ndarray,
    ba_hex: np.ndarray,
    mag_scatter: np.ndarray,
    ba_scatter: np.ndarray,
    rex_scatter: np.ndarray,
    title: str,
    out_stem: Path,
    n_total: int,
    gr_max: float,
    match_cdf_cuts: bool,
    q0: float,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    hb = axes[0].hexbin(
        mag_hex,
        ba_hex,
        gridsize=HEXBIN_GRIDSIZE,
        cmap="viridis",
        mincnt=1,
        norm=LogNorm(vmin=1),
        linewidths=0.0,
    )
    axes[0].axhline(q0, color="#d62728", linewidth=1.5, linestyle="--", label=f"strict cut $q_0={q0}$")
    axes[0].set_xlabel("tractor_mag_r")
    axes[0].set_ylabel("Tractor $b/a$ (from $e_1$, $e_2$)")
    axes[0].set_ylim(0, 1)
    axes[0].set_title(f"{title} — hex density (N={n_total:,})")
    axes[0].legend(loc="upper right", fontsize=7)
    fig.colorbar(hb, ax=axes[0], label="count")

    n = len(mag_scatter)
    if n > 0 and rex_scatter.size == n:
        not_rex = ~rex_scatter
        if not_rex.any():
            axes[1].scatter(
                mag_scatter[not_rex],
                ba_scatter[not_rex],
                s=4,
                alpha=0.25,
                c="#377eb8",
                linewidths=0,
                label="non-REX",
                rasterized=True,
            )
        if (~not_rex).any():
            axes[1].scatter(
                mag_scatter[~not_rex],
                ba_scatter[~not_rex],
                s=6,
                alpha=0.45,
                c="#ff7f00",
                linewidths=0,
                label="REX",
                rasterized=True,
            )
        axes[1].legend(loc="upper right", fontsize=8)
    else:
        axes[1].scatter(
            mag_scatter,
            ba_scatter,
            s=4,
            alpha=0.2,
            c="#377eb8",
            linewidths=0,
            rasterized=True,
        )

    axes[1].axhline(q0, color="#d62728", linewidth=1.5, linestyle="--", label=f"$q_0={q0}$")
    axes[1].set_xlabel("tractor_mag_r")
    axes[1].set_ylabel("Tractor $b/a$")
    axes[1].set_ylim(0, 1)
    axes[1].set_title(f"{title} — subsample (n={n:,})")

    extra = ""
    if match_cdf_cuts:
        extra = f"; also $b/a > {q0}$ and no REX (matches Legacy null CDF cuts except mag limit)"
    fig.text(
        0.5,
        0.01,
        f"Color cut: $g-r < {gr_max:g}$ on Tractor magnitudes.{extra}",
        ha="center",
        fontsize=8,
        color="0.35",
    )

    plt.tight_layout(rect=[0, 0.04, 1, 1])
    save_figure(fig, out_stem)
    plt.close(fig)


def run_legacy_color_cuts(
    legacy_csv: Path,
    out_root: Path,
    *,
    gr_max: float = LEGACY_GR_MAX_CDF,
    match_cdf_cuts: bool = True,
    exclude_rex: bool = True,
    q0: float = Q0,
    limit_chunks: int | None = None,
) -> None:
    out_root.mkdir(parents=True, exist_ok=True)
    print(
        f"[*] Streaming {legacy_csv.name} (chunksize={CHUNKSIZE}, "
        f"g-r < {gr_max:g}, match_cdf_cuts={match_cdf_cuts})"
    )
    if limit_chunks:
        print(f"[!] --limit-chunks {limit_chunks} (test mode)")

    data = _stream_legacy(
        legacy_csv,
        gr_max=gr_max,
        match_cdf_cuts=match_cdf_cuts,
        exclude_rex=exclude_rex,
        q0=q0,
        limit_chunks=limit_chunks,
    )
    n_total = data["n_pass"]
    print(f"[*] Passing rows: N={n_total:,}")

    tag = gr_cut_folder_tag(gr_max)
    sub = out_root / tag
    sub.mkdir(parents=True, exist_ok=True)

    title = f"Legacy ($g-r < {gr_max:g}$)"
    if match_cdf_cuts:
        title += f", $b/a > {q0}$, no REX"

    plot_mag_vs_ba_legacy(
        mag_hex=data["mag_hex"],
        ba_hex=data["ba_hex"],
        mag_scatter=data["mag_scatter"],
        ba_scatter=data["ba_scatter"],
        rex_scatter=data["rex_scatter"],
        title=title,
        out_stem=sub / "mag_vs_ba_legacy",
        n_total=n_total,
        gr_max=gr_max,
        match_cdf_cuts=match_cdf_cuts,
        q0=q0,
    )
    print(f"[*] {sub / 'mag_vs_ba_legacy.png'}")

    pd.DataFrame(
        [
            {
                "gr_max": gr_max,
                "n_legacy": n_total,
                "n_scatter": len(data["mag_scatter"]),
                "match_cdf_cuts": match_cdf_cuts,
                "exclude_rex": exclude_rex,
                "folder": tag,
            }
        ]
    ).to_csv(out_root / "legacy_color_cut_counts.csv", index=False)

    _append_legacy_color_cuts_md(out_root / "COLOR_CUTS.md", gr_max, n_total, match_cdf_cuts)


def _append_legacy_color_cuts_md(
    path: Path,
    gr_max: float,
    n_pass: int,
    match_cdf_cuts: bool,
) -> None:
    if path.is_file():
        text = path.read_text(encoding="utf-8")
        text = text.replace(
            "Legacy g-r hexbin diagnostics: **not generated** (deferred).",
            "",
        )
        marker = "## Legacy Survey"
        if marker in text:
            head, _, _tail = text.partition(marker)
            base = head.rstrip() + "\n"
        else:
            base = text.rstrip() + "\n"
    else:
        base = (
            "# Color cuts for null-catalog diagnostics\n\n"
            "See **`HEXBIN_AND_SDSS_BANDS_AUDIT.md`** for hexbin settings.\n\n"
        )

    cdf_note = (
        f", $b/a > {Q0}$, exclude REX (matches null CDF)"
        if match_cdf_cuts
        else " (color cut only)"
    )
    legacy_block = (
        f"## Legacy Survey (`plot_legacy_color_cuts.py`)\n\n"
        f"Tractor $g-r = \\mathrm{{gmag}} - \\mathrm{{rmag}}$; axis ratio = Tractor $b/a$ "
        f"from $e_1$, $e_2$ (column `expAB_r`).\n\n"
        f"| g-r cut | Role | N in cut |\n"
        f"|---------|------|----------|\n"
        f"| < {gr_max:g} | Diagnostic pool{cdf_note} (no mag limit) | {n_pass:,} |\n\n"
        f"Output: `color_cuts/{gr_cut_folder_tag(gr_max)}/mag_vs_ba_legacy.png`\n"
    )
    path.write_text(base + "\n" + legacy_block, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-csv", type=Path, default=DEFAULT_LEGACY)
    parser.add_argument("--out-dir", type=Path, default=DIAG_ROOT)
    parser.add_argument("--gr-max", type=float, default=LEGACY_GR_MAX_CDF)
    parser.add_argument(
        "--color-cut-only",
        action="store_true",
        help="Only apply g-r cut (do not also require b/a > q0 or exclude REX).",
    )
    parser.add_argument(
        "--include-rex",
        action="store_true",
        help="Keep REX rows when --color-cut-only (ignored if match_cdf_cuts).",
    )
    parser.add_argument("--q0", type=float, default=Q0)
    parser.add_argument("--limit-chunks", type=int, default=None)
    args = parser.parse_args()

    if not args.legacy_csv.is_file():
        raise SystemExit(f"Missing Legacy catalog: {args.legacy_csv}")

    match_cdf = not args.color_cut_only
    run_legacy_color_cuts(
        args.legacy_csv,
        args.out_dir,
        gr_max=args.gr_max,
        match_cdf_cuts=match_cdf,
        exclude_rex=not args.include_rex,
        q0=args.q0,
        limit_chunks=args.limit_chunks,
    )
    print(f"[*] Done: {args.out_dir}")


if __name__ == "__main__":
    main()
