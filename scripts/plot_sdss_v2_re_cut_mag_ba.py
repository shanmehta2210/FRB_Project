#!/usr/bin/env python3
"""
Large-galaxy Re cut on SDSS v2 catalog — mag composition vs median b/a.

Keeps rows with finite modelMag_r, expAB_r, and expRad_r (SDSS exponential-
profile effective radius in arcsec).

Compares:
  - Full sample (optionally strict production cuts on both overlays)
  - Subsample with expRad_r > re_min (same cuts + Re threshold)

Use ``--strict`` for production pool: u-r < 2.3, lnL exp-wins, expAB_r > 0.2.

Outputs under plots/plots_null/v2/sdss_audit/re_cut/:
  re_cut_mag_ba_summary[_strict]_re{N}.csv
  re_cut_mag_ba_panel[_strict]_re{N}.png

Run from repo root::

    python scripts/plot_sdss_v2_re_cut_mag_ba.py --re-min 4
    python scripts/plot_sdss_v2_re_cut_mag_ba.py --re-min 3 --strict
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
    apply_strict_q_cut,
    ensure_sdss_colors,
    filter_sdss_drop_dev_winners,
    filter_sdss_ur,
)
from pipeline_null_plot_utils import DEFAULT_SDSS_V2, REPO_ROOT  # noqa: E402

OUT_DIR = REPO_ROOT / "plots" / "plots_null" / "v2" / "sdss_audit" / "re_cut"
USECOLS = ("modelMag_r", "expAB_r", "expRad_r", "best_model_re_r")
USECOLS_STRICT = USECOLS + (
    "modelMag_u",
    "modelMag_g",
    "lnLExp_r",
    "lnLDeV_r",
)
MAG_STEP = 0.5
MAG_LO_START = 15.0
MIN_N_FULL = 30
MIN_N_RE = 5


def mag_bin_table(
    mag: np.ndarray,
    ba: np.ndarray,
    *,
    pool_n: int,
    step: float = MAG_STEP,
    min_n: int = 30,
) -> pd.DataFrame:
    rows: list[dict] = []
    lo = MAG_LO_START
    while lo < 28.0:
        hi = lo + step
        if lo == MAG_LO_START:
            mask = np.isfinite(mag) & (mag <= hi)
        else:
            mask = np.isfinite(mag) & (mag > lo) & (mag <= hi)
        n = int(mask.sum())
        if n < min_n:
            lo = hi
            continue
        rows.append(
            {
                "mag_lo": lo,
                "mag_hi": hi,
                "n": n,
                "frac_pool_pct": 100.0 * n / max(1, pool_n),
                "median_expAB_r": float(np.median(ba[mask])),
                "mean_expAB_r": float(np.mean(ba[mask])),
            }
        )
        lo = hi
    return pd.DataFrame(rows)


def apply_strict_cuts(df: pd.DataFrame) -> pd.DataFrame:
    """Production strict pool: u-r, lnL exp-wins, expAB_r > q0."""
    if SDSS_Q_COLUMN_CDF != "expAB_r":
        raise RuntimeError(f"Expected expAB_r for strict pools, got {SDSS_Q_COLUMN_CDF!r}")
    out = filter_sdss_ur(df, SDSS_UR_MAX_CDF)
    out = filter_sdss_drop_dev_winners(out)
    out = apply_strict_q_cut(out, q_col="expAB_r", q0=Q0)
    return out.reset_index(drop=True)


def load_catalog(path: Path, *, strict: bool) -> pd.DataFrame:
    usecols = USECOLS_STRICT if strict else USECOLS
    print(f"[*] Loading {path} ({', '.join(usecols)}) ...")
    df = pd.read_csv(path, usecols=list(usecols))
    if strict:
        df = ensure_sdss_colors(df)
    mag = pd.to_numeric(df["modelMag_r"], errors="coerce")
    ba = pd.to_numeric(df["expAB_r"], errors="coerce")
    re = pd.to_numeric(df["expRad_r"], errors="coerce")
    ok = mag.notna() & ba.notna() & re.notna() & (ba >= 0) & (ba <= 1) & (re > 0)
    out = df.loc[ok].copy()
    out["modelMag_r"] = mag[ok].to_numpy()
    out["expAB_r"] = ba[ok].to_numpy()
    out["expRad_r"] = re[ok].to_numpy()
    print(f"[*] Finite mag + expAB_r + expRad_r: N={len(out):,}")
    if strict:
        out = apply_strict_cuts(out)
        print(
            f"[*] After strict cuts (u-r<{SDSS_UR_MAX_CDF:g}, lnL exp-wins, "
            f"expAB_r>{Q0:g}): N={len(out):,}"
        )
    return out


def plot_panel(
    full_bins: pd.DataFrame,
    re_bins: pd.DataFrame,
    *,
    re_min: float,
    n_full: int,
    n_re: int,
    out_png: Path,
    strict: bool,
) -> None:
    x_full = 0.5 * (full_bins["mag_lo"] + full_bins["mag_hi"])
    x_re = 0.5 * (re_bins["mag_lo"] + re_bins["mag_hi"])

    fig, ax_ba = plt.subplots(figsize=(10, 5.5))
    ax_pct = ax_ba.twinx()

    ax_ba.plot(
        x_full,
        full_bins["median_expAB_r"],
        "o--",
        color="#377eb8",
        alpha=0.85,
        lw=1.5,
        ms=5,
        label=f"{'Strict pool' if strict else 'Full catalog'} median expAB$_r$ (N={n_full:,})",
    )
    ax_ba.plot(
        x_re,
        re_bins["median_expAB_r"],
        "s-",
        color="#e41a1c",
        lw=2,
        ms=6,
        label=(
            f"strict + expRad$_r$ > {re_min:g}\" median expAB$_r$ (N={n_re:,})"
            if strict
            else f"expRad$_r$ > {re_min:g}\" median expAB$_r$ (N={n_re:,})"
        ),
    )

    ax_pct.plot(
        x_re,
        re_bins["frac_pool_pct"],
        "^-",
        color="#4daf4a",
        lw=1.8,
        ms=5,
        alpha=0.9,
        label=f"expRad$_r$ > {re_min:g}\" % per mag bin",
    )

    ax_ba.set_xlabel(r"modelMag$_r$ (0.5 mag bins)")
    ax_ba.set_ylabel("median expAB$_r$ (raw SDSS)", color="#333333")
    ax_pct.set_ylabel(f"% of Re > {re_min:g}\" sample per bin", color="#4daf4a")
    ax_pct.tick_params(axis="y", labelcolor="#4daf4a")
    ax_ba.set_ylim(0, 1.02)
    ax_ba.grid(True, alpha=0.3)

    lines_ba, labels_ba = ax_ba.get_legend_handles_labels()
    lines_pct, labels_pct = ax_pct.get_legend_handles_labels()
    ax_ba.legend(lines_ba + lines_pct, labels_ba + labels_pct, loc="upper right", fontsize=8)

    cut_note = (
        f"Strict: u-r < {SDSS_UR_MAX_CDF:g}, lnL exp-wins, expAB$_r$ > {Q0:g}"
        if strict
        else "No u-r / lnL / b/a cuts"
    )
    ax_ba.set_title(
        f"SDSS v2 {'strict pool' if strict else 'full catalog'}: Re cut vs median b/a\n"
        f"{cut_note}  |  Re = PhotoObj expRad$_r$ (arcsec)  |  shape = expAB$_r$"
    )
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sdss-csv", type=Path, default=DEFAULT_SDSS_V2)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--re-min", type=float, default=4.0, help="expRad_r cut (arcsec).")
    parser.add_argument("--re-column", default="expRad_r", help="Re column (arcsec).")
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Apply production strict cuts to full overlay and Re subsample.",
    )
    args = parser.parse_args()

    tag = "strict" if args.strict else "raw"
    df = load_catalog(args.sdss_csv, strict=args.strict)
    if args.re_column not in df.columns:
        raise KeyError(f"Missing {args.re_column!r} in catalog")

    re = pd.to_numeric(df[args.re_column], errors="coerce").to_numpy()
    mag = df["modelMag_r"].to_numpy(dtype=float)
    ba = df["expAB_r"].to_numpy(dtype=float)

    n_full = len(df)
    re_mask = re > args.re_min
    n_re = int(re_mask.sum())
    pass_frac = 100.0 * n_re / max(1, n_full)
    print(
        f"[*] {args.re_column} > {args.re_min:g}\" arcsec: N={n_re:,} "
        f"({pass_frac:.2f}% of finite sample)"
    )
    print(
        f"[*] {args.re_column} percentiles (full): "
        f"p50={np.percentile(re, 50):.3f}\", "
        f"p90={np.percentile(re, 90):.3f}\", "
        f"p99={np.percentile(re, 99):.3f}\""
    )

    full_bins = mag_bin_table(mag, ba, pool_n=n_full, min_n=MIN_N_FULL)
    re_bins = mag_bin_table(mag[re_mask], ba[re_mask], pool_n=n_re, min_n=MIN_N_RE)

    merged = full_bins.merge(
        re_bins,
        on=["mag_lo", "mag_hi"],
        how="outer",
        suffixes=("_full", "_re_cut"),
    ).sort_values("mag_lo")
    merged["re_min_arcsec"] = args.re_min
    merged["re_column"] = args.re_column
    merged["sample_mode"] = tag
    merged["n_full_catalog"] = n_full
    merged["n_re_cut"] = n_re

    args.out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"re_cut_mag_ba_{tag}_re{args.re_min:g}"
    csv_path = args.out_dir / f"{stem}_summary.csv"
    merged.to_csv(csv_path, index=False)
    print(f"[*] Wrote {csv_path}")

    png_path = args.out_dir / f"{stem}_panel.png"
    plot_panel(
        full_bins,
        re_bins,
        re_min=args.re_min,
        n_full=n_full,
        n_re=n_re,
        out_png=png_path,
        strict=args.strict,
    )
    print(f"[*] Saved {png_path}")

    print("\n--- Re-cut sample mag composition ---")
    print(re_bins[["mag_lo", "mag_hi", "n", "frac_pool_pct", "median_expAB_r"]].to_string(index=False))
    print("\n--- Median b/a delta (Re-cut minus full) ---")
    both = merged.dropna(subset=["median_expAB_r_full", "median_expAB_r_re_cut"])
    both = both.assign(
        delta_median_ba=both["median_expAB_r_re_cut"] - both["median_expAB_r_full"]
    )
    print(both[["mag_lo", "mag_hi", "median_expAB_r_full", "median_expAB_r_re_cut", "delta_median_ba"]].to_string(index=False))


if __name__ == "__main__":
    main()
