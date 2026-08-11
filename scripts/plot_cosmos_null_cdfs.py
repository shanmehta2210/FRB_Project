#!/usr/bin/env python3
"""COSMOS field-level cos(i) CDF comparison: HST Zurich vs SDSS.

Output layout::

    plots/plots_null/v2/sdss_audit/COSMOS/plots/cdfs/mag20/strict.png
    plots/plots_null/v2/sdss_audit/COSMOS/plots/cdfs/mag20/entire.png
    plots/plots_null/v2/sdss_audit/COSMOS/plots/cdfs/mag21/...
    plots/plots_null/v2/sdss_audit/COSMOS/plots/cdfs/mag22/...

Run from repo root::

    python scripts/plot_cosmos_null_cdfs.py
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

_SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS))

from cosmos_ba_audit_utils import (  # noqa: E402
    COSMOS_BA_MIN,
    COSMOS_CDFS,
    HST_BA_COL,
    HST_DISK_CSV,
    HST_DISK_ENTIRE_CSV,
    HST_MAG_COL,
    SDSS_BA_COL,
    SDSS_DISK_CSV,
    SDSS_DISK_ENTIRE_CSV,
    SDSS_MAG_COL,
)
from null_catalog_utils import Q0, cosi_array_from_df  # noqa: E402
from pipeline_null_plot_utils import (  # noqa: E402
    add_inclination_top_axis,
    cdf_envelope,
    default_font,
    save_figure,
)

MAG_CUTS: tuple[tuple[str, float], ...] = (
    ("mag20", 20.0),
    ("mag21", 21.0),
    ("mag22", 22.0),
)
HST_COLOR = "#e41a1c"
SDSS_COLOR = "#377eb8"


def slice_by_mag(df: pd.DataFrame, mag_col: str, mag_max: float) -> pd.DataFrame:
    mag = pd.to_numeric(df[mag_col], errors="coerce")
    return df.loc[mag.notna() & (mag <= mag_max)].copy()


def plot_dual_cosi_cdf(
    *,
    hst_cosi: np.ndarray,
    sdss_cosi: np.ndarray,
    mag_cut: float,
    mode_label: str,
    n_sample: int,
    n_draws: int,
    out_stem: Path,
) -> tuple[float, float, int, int]:
    font_prop = default_font()
    fig, ax = plt.subplots(figsize=(8, 8))

    for label, color, cosi in (
        ("HST GIM2D disk (TYPE=2)", HST_COLOR, hst_cosi),
        ("SDSS DR17 disk (exp wins)", SDSS_COLOR, sdss_cosi),
    ):
        x_n, mn, lo, hi = cdf_envelope(cosi, n_sample=n_sample, n_draws=n_draws)
        ax.plot(x_n, mn, color=color, linewidth=2.0, label=f"{label} (N={len(cosi):,})")
        ax.fill_between(x_n, lo, hi, color=color, alpha=0.18)

    ks = stats.ks_2samp(hst_cosi, sdss_cosi)
    ba_note = f"b/a > {COSMOS_BA_MIN:g}" if mode_label == "strict" else "all b/a"
    ax.plot((0, 1), (0, 1), "k--", linewidth=1.0, label="Uniform")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("cos(i)", fontproperties=font_prop, fontsize=11)
    ax.set_ylabel("Cumulative distribution", fontproperties=font_prop, fontsize=11)
    ax.set_title(
        f"COSMOS ACS — cos(i) CDFs ({mode_label}, mag <= {mag_cut:g})\n"
        f"disk only; {ba_note}; no colour cut; KS D={ks.statistic:.3f}, p={ks.pvalue:.2e}",
        fontsize=11,
    )
    ax.legend(fontsize=8, loc="upper left")
    add_inclination_top_axis(ax, font_prop)
    plt.tight_layout()
    save_figure(fig, out_stem)
    plt.close(fig)
    return float(ks.statistic), float(ks.pvalue), len(hst_cosi), len(sdss_cosi)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hst-entire-csv", type=Path, default=HST_DISK_ENTIRE_CSV)
    parser.add_argument("--hst-strict-csv", type=Path, default=HST_DISK_CSV)
    parser.add_argument("--sdss-entire-csv", type=Path, default=SDSS_DISK_ENTIRE_CSV)
    parser.add_argument("--sdss-strict-csv", type=Path, default=SDSS_DISK_CSV)
    parser.add_argument("--out-dir", type=Path, default=COSMOS_CDFS)
    parser.add_argument("--mc-draws", type=int, default=5000)
    parser.add_argument("--q0", type=float, default=Q0)
    parser.add_argument("--clean", action="store_true", help="Remove old files under out-dir first.")
    args = parser.parse_args()

    hst_entire = pd.read_csv(args.hst_entire_csv)
    hst_strict = pd.read_csv(args.hst_strict_csv)
    sdss_entire = pd.read_csv(args.sdss_entire_csv)
    sdss_strict = pd.read_csv(args.sdss_strict_csv)

    if args.clean and args.out_dir.exists():
        shutil.rmtree(args.out_dir)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    for folder, mag_cut in MAG_CUTS:
        mag_dir = args.out_dir / folder
        mag_dir.mkdir(parents=True, exist_ok=True)

        for mode, hst_df, sdss_df in (
            ("entire", hst_entire, sdss_entire),
            ("strict", hst_strict, sdss_strict),
        ):
            h_sub = slice_by_mag(hst_df, HST_MAG_COL, mag_cut)
            s_sub = slice_by_mag(sdss_df, SDSS_MAG_COL, mag_cut)

            h_cosi = cosi_array_from_df(h_sub, q_col=HST_BA_COL, q0=args.q0)
            s_cosi = cosi_array_from_df(s_sub, q_col=SDSS_BA_COL, q0=args.q0)

            n_sample = min(len(h_cosi), len(s_cosi))
            if n_sample < 50:
                print(
                    f"[!] {folder}/{mode}: pool too small "
                    f"(HST={len(h_cosi)}, SDSS={len(s_cosi)})"
                )
                continue

            out_stem = mag_dir / mode
            ks_d, ks_p, n_h, n_s = plot_dual_cosi_cdf(
                hst_cosi=h_cosi,
                sdss_cosi=s_cosi,
                mag_cut=mag_cut,
                mode_label=mode,
                n_sample=n_sample,
                n_draws=args.mc_draws,
                out_stem=out_stem,
            )
            rows.append(
                {
                    "mag_folder": folder,
                    "mag_cut": mag_cut,
                    "mode": mode,
                    "n_hst": n_h,
                    "n_sdss": n_s,
                    "n_sample_mc": n_sample,
                    "ks_D": ks_d,
                    "ks_pvalue": ks_p,
                }
            )
            print(f"[*] {folder}/{mode}: HST N={n_h:,}, SDSS N={n_s:,}, KS D={ks_d:.3f}")

    if rows:
        pd.DataFrame(rows).to_csv(args.out_dir / "ks_summary.csv", index=False)
    print(f"[*] CDF plots written to {args.out_dir}")


if __name__ == "__main__":
    main()
