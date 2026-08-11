#!/usr/bin/env python3
"""Merge truth + GALFIT results, aggregate realizations, plot b/a vs magnitude."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sim_utils import TOOL_DIR, ensure_output_layout, load_config


def _to_numeric(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for c in cols:
        if c in out.columns:
            out[c] = pd.to_numeric(out[c], errors="coerce")
    return out


def merge_catalogs(truth_path: Path, fit_path: Path) -> pd.DataFrame:
    truth = pd.read_csv(truth_path)
    fits = pd.read_csv(fit_path)
    merged = fits.merge(
        truth,
        on=["galaxy_id", "realization"],
        how="left",
        suffixes=("", "_truth"),
    )
    num_cols = [
        "ba_fit",
        "re_fit_pix",
        "mag_fit",
        "n_fit",
        "pa_fit",
        "chi2nu",
        "sky_fit_e",
        "ba_true",
        "re_arcsec_true",
        "re_pix_true",
        "mag_true",
        "pa_true",
        "snr_at_re",
    ]
    return _to_numeric(merged, num_cols)


def aggregate_realizations(df: pd.DataFrame) -> pd.DataFrame:
    group_cols = ["mode", "re_arcsec_true", "ba_true", "mag_true"]
    rows = []
    for keys, grp in df.groupby(group_cols, sort=True):
        mode, re_arc, ba, mag = keys
        ba_vals = grp["ba_fit"].dropna()
        if ba_vals.empty:
            continue
        q16, q50, q84 = np.percentile(ba_vals, [16, 50, 84])
        rows.append(
            {
                "mode": mode,
                "re_arcsec_true": re_arc,
                "ba_true": ba,
                "mag_true": mag,
                "n_realizations": len(ba_vals),
                "ba_fit_median": round(float(q50), 4),
                "ba_fit_std": round(float(ba_vals.std(ddof=0)), 4) if len(ba_vals) > 1 else 0.0,
                "ba_fit_p16": round(float(q16), 4),
                "ba_fit_p84": round(float(q84), 4),
            }
        )
    return pd.DataFrame(rows)


def build_summary(agg: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (mode, re_arc, ba), grp in agg.groupby(["mode", "re_arcsec_true", "ba_true"], sort=True):
        grp = grp.sort_values("mag_true")
        below = grp.loc[grp["ba_fit_median"] < 0.5, "mag_true"]
        rows.append(
            {
                "mode": mode,
                "re_arcsec_true": re_arc,
                "ba_true": ba,
                "n_mag_bins": len(grp),
                "n_realizations_max": int(grp["n_realizations"].max()),
                "ba_min_median": round(float(grp["ba_fit_median"].min()), 4),
                "mag_first_below_0p5": round(float(below.iloc[0]), 2) if len(below) else "",
            }
        )
    return pd.DataFrame(rows)


def pivot_psf_nopsf(agg: pd.DataFrame) -> pd.DataFrame:
    """One row per (Re, ba_true, mag) with PSF and no-PSF columns side by side."""
    psf = agg[agg["mode"] == "psf"].rename(
        columns={
            "ba_fit_median": "ba_psf",
            "ba_fit_p16": "ba_psf_p16",
            "ba_fit_p84": "ba_psf_p84",
        }
    )
    nopsf = agg[agg["mode"] == "nopsf"].rename(
        columns={
            "ba_fit_median": "ba_nopsf",
            "ba_fit_p16": "ba_nopsf_p16",
            "ba_fit_p84": "ba_nopsf_p84",
        }
    )
    keys = ["re_arcsec_true", "ba_true", "mag_true"]
    wide = psf.merge(nopsf, on=keys, how="outer", suffixes=("_psf", "_nopsf"))
    wide["ba_true"] = wide["ba_true"].astype(float)
    wide["re_arcsec_true"] = wide["re_arcsec_true"].astype(float)
    wide["mag_true"] = wide["mag_true"].astype(float)
    return wide.sort_values(keys)


def plot_truth_psf_nopsf(wide: pd.DataFrame, out_path: Path, cfg: dict | None = None) -> None:
    """
    3x3 grid: rows = intrinsic b/a, cols = Re.
    Each panel shows three curves vs magnitude:
      - true b/a (horizontal)
      - GALFIT with PSF
      - GALFIT without PSF
    """
    ba_values = sorted(wide["ba_true"].unique())
    re_values = sorted(wide["re_arcsec_true"].unique())
    n_ba, n_re = len(ba_values), len(re_values)

    fig, axes = plt.subplots(
        n_ba,
        n_re,
        figsize=(4.2 * n_re, 3.6 * n_ba),
        sharex=True,
        sharey=True,
        squeeze=False,
    )

    for i, ba in enumerate(ba_values):
        for j, re_arc in enumerate(re_values):
            ax = axes[i, j]
            sub = wide[
                (wide["ba_true"] == ba) & (wide["re_arcsec_true"] == re_arc)
            ].sort_values("mag_true")
            if sub.empty:
                ax.set_visible(False)
                continue

            x = sub["mag_true"].values

            # True b/a: flat reference (same for all mags in this panel).
            ax.axhline(
                ba,
                color="black",
                linestyle="--",
                linewidth=2.0,
                label=rf"True $b/a={ba:.1f}$",
                zorder=1,
            )

            psf_y = sub["ba_psf"].values
            ax.plot(
                x,
                psf_y,
                color="C0",
                marker="o",
                markersize=4,
                linewidth=1.8,
                label="GALFIT + PSF",
                zorder=3,
            )

            nopsf_y = sub["ba_nopsf"].values
            ax.plot(
                x,
                nopsf_y,
                color="C3",
                marker="s",
                markersize=3.5,
                linewidth=1.8,
                label="GALFIT, no PSF",
                zorder=2,
            )

            ax.set_title(rf"True $b/a={ba:.1f}$, $R_e={re_arc:.1f}''$", fontsize=10)
            if j == 0:
                ax.set_ylabel(r"Reported $b/a$")
            if i == n_ba - 1:
                ax.set_xlabel(r"Input $m_r$ (bright $\leftarrow$ faint)")

            ax.set_xlim(16.5, 24.5)
            ax.set_ylim(0.0, 1.0)
            ax.grid(True, alpha=0.3)
            ax.invert_xaxis()  # bright (17) on the left

            if i == 0 and j == n_re - 1:
                ax.legend(fontsize=8, loc="upper right")

    fig.suptitle(
        "True b/a vs GALFIT-reported b/a (with PSF deconvolution vs without)"
        + (" [mag+sky locked to truth]" if (cfg or {}).get("galfit", {}).get("lock_mag") else ""),
        fontsize=13,
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=TOOL_DIR / "config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    layout = ensure_output_layout(cfg)
    truth_path = layout["catalogs"] / "truth_catalog.csv"
    fit_path = layout["catalogs"] / "fit_results.csv"
    if not truth_path.is_file() or not fit_path.is_file():
        print("Missing truth_catalog.csv or fit_results.csv; run prior steps first.")
        return 1

    merged = merge_catalogs(truth_path, fit_path)
    merged["converged"] = merged["converged"].astype(str).str.lower().eq("true")
    merged_path = layout["catalogs"] / "merged.csv"
    merged.to_csv(merged_path, index=False)

    plot_df = merged[merged["converged"]].copy()
    if plot_df.empty:
        print("No converged fits to plot; run audit_fits.py and re-fit with fixed config.")
        return 1

    agg = aggregate_realizations(plot_df)
    agg_path = layout["catalogs"] / "aggregated.csv"
    agg.to_csv(agg_path, index=False)

    summary = build_summary(agg)
    summary_path = layout["catalogs"] / "summary.csv"
    summary.to_csv(summary_path, index=False)

    wide = pivot_psf_nopsf(agg)
    wide_path = layout["catalogs"] / "wide_psf_nopsf.csv"
    wide.to_csv(wide_path, index=False)

    plot_truth_psf_nopsf(wide, layout["plots"] / "ba_vs_mag_combined.png", cfg)

    print(f"Wrote merged -> {merged_path}")
    print(f"Wrote aggregated -> {agg_path}")
    print(f"Wrote wide_psf_nopsf -> {wide_path}")
    print(f"Wrote summary -> {summary_path}")
    print(f"Plots -> {layout['plots']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
