#!/usr/bin/env python3
"""
Audit v1 null catalogs: cut funnel table + diagnostic figures.

Outputs under plots/plots_null/v1_null_cdf_inclination/diagnostics/

RAM: default does NOT load v1 null CSVs. Full funnel + mag vs b/a needs --full
(~500k SDSS + ~500k Legacy rows in memory).

SDSS u-r color-cut plots (low RAM):
    python scripts/plot_sdss_color_cuts.py

Full audit (high RAM):
    python scripts/audit_and_plot_null_v1_diagnostics.py --full
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
    LEGACY_CDF_TYPE_EXCLUDE,
    LEGACY_GR_MAX_CDF,
    Q0,
    SDSS_BA_FLOOR_MIN,
    SDSS_Q_COLUMN_CDF,
    SDSS_UR_MAX_CDF,
    cut_funnel_rows,
    ensure_sdss_colors,
    resolve_mag_column,
    resolve_q_column,
)
from pipeline_null_plot_utils import (  # noqa: E402
    DEFAULT_LEGACY,
    DEFAULT_SDSS,
    PLOTS_NULL,
    save_figure,
)

DIAG_ROOT = PLOTS_NULL / "v1_null_cdf_inclination" / "diagnostics"
MAG_LIMITS_AUDIT = [21.0, 19, 17, 16, 15]
MAG_BIN_EDGES = np.arange(14.0, 22.5, 1.0)
SUBSAMPLE_N = 15_000
SUBSAMPLE_SEED = 42


def build_cut_funnel(
    legacy: pd.DataFrame,
    sdss: pd.DataFrame,
    mag_limits: list[float],
    q0: float,
    exclude_types: str,
) -> pd.DataFrame:
    rows: list[dict] = []
    for mag_limit in mag_limits:
        rows.extend(
            cut_funnel_rows(
                legacy,
                survey="legacy",
                mag_limit=mag_limit,
                sample_mode="strict",
                mag_column="rmag",
                q_column="expAB_r",
                q0=q0,
                exclude_legacy_types=exclude_types,
                is_legacy=True,
                legacy_gr_max=LEGACY_GR_MAX_CDF,
                legacy_spiral_morph_only=True,
                legacy_cdf_type_exclude=LEGACY_CDF_TYPE_EXCLUDE,
            )
        )
        rows.extend(
            cut_funnel_rows(
                sdss,
                survey="sdss",
                mag_limit=mag_limit,
                sample_mode="strict",
                mag_column="modelMag_r",
                q_column=SDSS_Q_COLUMN_CDF,
                q0=q0,
                is_legacy=False,
                exclude_sdss_ba_floor=False,
                sdss_ur_max=SDSS_UR_MAX_CDF,
                sdss_exp_winner_only=True,
                legacy_spiral_morph_only=True,
                legacy_cdf_type_exclude=LEGACY_CDF_TYPE_EXCLUDE,
            )
        )
    return pd.DataFrame(rows)


def plot_mag_vs_ba(
    df: pd.DataFrame,
    *,
    mag_col: str,
    ba_col: str,
    title: str,
    out_stem: Path,
    type_col: str | None = None,
    sdss_floor_min: float | None = None,
) -> None:
    mag_c = resolve_mag_column(df, mag_col)
    ba_c = resolve_q_column(df, ba_col)
    mag = pd.to_numeric(df[mag_c], errors="coerce").to_numpy(dtype=float)
    ba = pd.to_numeric(df[ba_c], errors="coerce").to_numpy(dtype=float)
    ok = np.isfinite(mag) & np.isfinite(ba) & (ba >= 0) & (ba <= 1)
    mag, ba = mag[ok], ba[ok]

    if type_col and type_col in df.columns:
        types = df.loc[ok, type_col].astype(str).str.upper().to_numpy()
        is_rex = types == "REX"
    else:
        is_rex = None

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    if sdss_floor_min is not None:
        for ax in axes:
            ax.axhspan(
                0,
                sdss_floor_min,
                color="#d62728",
                alpha=0.12,
                zorder=0,
                label=f"PhotoObj $b/a$ floor ($\\leq {sdss_floor_min}$)",
            )

    hb = axes[0].hexbin(
        mag,
        ba,
        gridsize=50,
        cmap="viridis",
        mincnt=1,
        norm=LogNorm(vmin=1),
        linewidths=0.0,
    )
    axes[0].axhline(Q0, color="#d62728", linewidth=1.5, linestyle="--", label=f"strict cut $q_0={Q0}$")
    axes[0].set_xlabel(mag_c)
    axes[0].set_ylabel(ba_c)
    axes[0].set_ylim(0, 1)
    axes[0].set_title(f"{title} — hex density (N={len(mag):,})")
    axes[0].legend(loc="upper right", fontsize=7)
    fig.colorbar(hb, ax=axes[0], label="count")

    rng = np.random.default_rng(SUBSAMPLE_SEED)
    n = min(SUBSAMPLE_N, len(mag))
    idx = rng.choice(len(mag), size=n, replace=False) if len(mag) > n else np.arange(len(mag))
    if is_rex is not None:
        axes[1].scatter(
            mag[idx][~is_rex[idx]],
            ba[idx][~is_rex[idx]],
            s=4,
            alpha=0.25,
            c="#377eb8",
            linewidths=0,
            label="non-REX",
            rasterized=True,
        )
        axes[1].scatter(
            mag[idx][is_rex[idx]],
            ba[idx][is_rex[idx]],
            s=6,
            alpha=0.45,
            c="#ff7f00",
            linewidths=0,
            label="REX",
            rasterized=True,
        )
        axes[1].legend(loc="upper right", fontsize=8)
    else:
        axes[1].scatter(mag[idx], ba[idx], s=4, alpha=0.2, c="#4daf4a", linewidths=0, rasterized=True)
    axes[1].axhline(Q0, color="#d62728", linewidth=1.5, linestyle="--", label=f"$q_0={Q0}$")
    axes[1].set_xlabel(mag_c)
    axes[1].set_ylabel(ba_c)
    axes[1].set_ylim(0, 1)
    axes[1].legend(loc="upper right", fontsize=7)
    axes[1].set_title(f"{title} — subsample (n={n:,})")

    if sdss_floor_min is not None:
        n_floor = int((ba <= sdss_floor_min).sum())
        fig.text(
            0.5,
            0.01,
            f"Red band: SDSS PhotoObj floor at b/a={sdss_floor_min} ({n_floor:,} objects, "
            f"{100 * n_floor / len(ba):.1f}% of catalog). "
            f"Included in SDSS null CDF pools (not excluded).",
            ha="center",
            fontsize=8,
            color="0.35",
        )

    plt.tight_layout(rect=[0, 0.04, 1, 1])
    save_figure(fig, out_stem)
    plt.close(fig)


def plot_rex_fraction(legacy: pd.DataFrame, out_stem: Path) -> pd.DataFrame:
    mag_c = resolve_mag_column(legacy, "rmag")
    mag = pd.to_numeric(legacy[mag_c], errors="coerce")
    types = legacy["tractor_type"].astype(str).str.upper()
    ok = mag.notna() & np.isfinite(mag)
    mag = mag.loc[ok]
    types = types.loc[ok]

    bins = MAG_BIN_EDGES
    records = []
    for i in range(len(bins) - 1):
        lo, hi = bins[i], bins[i + 1]
        mask = (mag >= lo) & (mag < hi)
        n_tot = int(mask.sum())
        n_rex = int((types[mask] == "REX").sum())
        frac = n_rex / n_tot if n_tot > 0 else float("nan")
        records.append(
            {
                "mag_bin_lo": lo,
                "mag_bin_hi": hi,
                "n_total": n_tot,
                "n_rex": n_rex,
                "n_non_rex": n_tot - n_rex,
                "frac_rex": frac,
            }
        )
    rex_df = pd.DataFrame(records)

    fig, ax = plt.subplots(figsize=(9, 5))
    centers = 0.5 * (rex_df["mag_bin_lo"] + rex_df["mag_bin_hi"])
    ax.bar(
        centers,
        rex_df["frac_rex"] * 100,
        width=0.85,
        color="#ff7f00",
        edgecolor="white",
        linewidth=0.6,
    )
    ax.set_xlabel(f"{mag_c} (1 mag bins)")
    ax.set_ylabel("REX fraction (%)")
    ax.set_title("Legacy Tractor: fraction of type=REX vs magnitude (full catalog)")
    ax.set_xlim(bins[0] - 0.5, bins[-1] + 0.5)
    ymax = float(np.nanmax(rex_df["frac_rex"]) * 100 * 1.15) if rex_df["frac_rex"].notna().any() else 100
    ax.set_ylim(0, max(10, ymax))
    for x, row in zip(centers, rex_df.itertuples()):
        if row.n_total > 0:
            ax.text(
                x,
                row.frac_rex * 100 + 0.5,
                f"{row.n_rex}/{row.n_total}",
                ha="center",
                va="bottom",
                fontsize=7,
            )
    plt.tight_layout()
    save_figure(fig, out_stem)
    plt.close(fig)
    return rex_df


def plot_frac_q_le_q0(
    sdss: pd.DataFrame,
    legacy: pd.DataFrame,
    out_stem: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    for df, mag_col, ba_col, color, label in (
        (sdss, "modelMag_r", "expAB_r", "#4daf4a", "SDSS expAB_r"),
        (legacy, "rmag", "expAB_r", "#377eb8", "Legacy Tractor"),
    ):
        mag_c = resolve_mag_column(df, mag_col)
        ba_c = resolve_q_column(df, ba_col)
        mag = pd.to_numeric(df[mag_c], errors="coerce")
        ba = pd.to_numeric(df[ba_c], errors="coerce")
        ok = mag.notna() & ba.notna()
        mag = mag.loc[ok]
        ba = ba.loc[ok]
        fracs, centers = [], []
        for i in range(len(MAG_BIN_EDGES) - 1):
            lo, hi = MAG_BIN_EDGES[i], MAG_BIN_EDGES[i + 1]
            mask = (mag >= lo) & (mag < hi)
            n = int(mask.sum())
            if n == 0:
                continue
            fracs.append(float((ba[mask] <= Q0).sum()) / n)
            centers.append(0.5 * (lo + hi))
        ax.plot(centers, np.array(fracs) * 100, "o-", color=color, label=label, linewidth=2)

    ax.axhline(Q0 * 100, color="0.4", linestyle=":", linewidth=1, alpha=0.5)
    ax.set_xlabel("r-band magnitude (1 mag bins)")
    ax.set_ylabel(f"Fraction with $b/a \\leq q_0$ ({Q0}) (%)")
    ax.set_title(f"Edge-on floor population vs magnitude ($q_0={Q0}$)")
    ax.legend(fontsize=9)
    ax.set_xlim(MAG_BIN_EDGES[0], MAG_BIN_EDGES[-1])
    plt.tight_layout()
    save_figure(fig, out_stem)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-csv", type=Path, default=DEFAULT_LEGACY)
    parser.add_argument("--sdss-csv", type=Path, default=DEFAULT_SDSS)
    parser.add_argument("--out-dir", type=Path, default=DIAG_ROOT)
    parser.add_argument("--q0", type=float, default=Q0)
    parser.add_argument(
        "--exclude-types",
        default="REX",
        help="Type exclusion for mag-vs-b/a diagnostics (REX visible); funnel uses CDF cuts.",
    )
    parser.add_argument(
        "--mag-limits",
        type=float,
        nargs="+",
        default=MAG_LIMITS_AUDIT,
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="Load full v1 Legacy+SDSS CSVs and regenerate funnel, mag vs b/a, REX plots (high RAM).",
    )
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if not args.full:
        print(
            "Nothing to do: this script no longer runs by default (avoids loading ~1M null rows).\n"
            "  SDSS u-r color cuts (low RAM):  python scripts/plot_sdss_color_cuts.py\n"
            "  Full funnel + mag vs b/a:     python scripts/audit_and_plot_null_v1_diagnostics.py --full\n"
            "  Mag-cut CDFs:                 python scripts/plot_null_mag_cut_cdfs.py --no-clear\n"
            "See plots/plots_null/v1_null_cdf_inclination/diagnostics/MEMORY_SAFE_NULL_WORK.md"
        )
        return

    legacy = pd.read_csv(args.legacy_csv)
    sdss = pd.read_csv(args.sdss_csv)
    print(f"Legacy rows: {len(legacy):,}  SDSS rows: {len(sdss):,}")

    funnel = build_cut_funnel(
        legacy,
        sdss,
        mag_limits=args.mag_limits,
        q0=args.q0,
        exclude_types=args.exclude_types,
    )
    funnel_path = args.out_dir / "cut_funnel.csv"
    funnel.to_csv(funnel_path, index=False)
    print(f"Wrote {funnel_path} ({len(funnel)} rows)")

    sdss = ensure_sdss_colors(sdss)

    plot_mag_vs_ba(
        sdss,
        mag_col="modelMag_r",
        ba_col="best_model_ba_r",
        title="SDSS DR16 null catalog (full, incl. floor rung)",
        out_stem=args.out_dir / "mag_vs_ba_sdss",
        sdss_floor_min=SDSS_BA_FLOOR_MIN,
    )
    plot_mag_vs_ba(
        legacy,
        mag_col="rmag",
        ba_col="expAB_r",
        title="Legacy DR10 Tractor null catalog",
        out_stem=args.out_dir / "mag_vs_ba_legacy",
        type_col="tractor_type",
    )

    rex_df = plot_rex_fraction(legacy, args.out_dir / "rex_fraction_vs_mag")
    rex_csv = args.out_dir / "rex_fraction_by_mag.csv"
    rex_df.to_csv(rex_csv, index=False)
    print(f"Wrote {rex_csv}")

    plot_frac_q_le_q0(sdss, legacy, args.out_dir / "frac_q_le_q0_vs_mag")

    print(f"[*] Diagnostics complete: {args.out_dir}")
    print("[*] SDSS u-r color cuts (separate, low RAM): python scripts/plot_sdss_color_cuts.py")


if __name__ == "__main__":
    main()
