#!/usr/bin/env python3
"""
Strict SDSS v2 pool: median expAB_r vs mag **within fixed Re slices**.

Avoids the Re-threshold confound (conditioning on Re > X at fixed mag selects a
different population). Instead, slice Re first, then track b/a vs mag inside each
comparable size bin.

Default Re slices (arcsec): [1, 2), [2, 3), [3, ∞).

Outputs under plots/plots_null/v2/sdss_audit/re_cut/:
  re_slice_mag_ba_summary.csv
  re_slice_mag_ba_panel.png

Run from repo root::

    python scripts/plot_sdss_v2_re_slice_mag_ba.py
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

from plot_sdss_v2_re_cut_mag_ba import (  # noqa: E402
    MAG_LO_START,
    MAG_STEP,
    USECOLS_STRICT,
    apply_strict_cuts,
    mag_bin_table,
)
from null_catalog_utils import Q0, SDSS_UR_MAX_CDF, ensure_sdss_colors  # noqa: E402
from pipeline_null_plot_utils import DEFAULT_SDSS_V2, REPO_ROOT  # noqa: E402

OUT_DIR = REPO_ROOT / "plots" / "plots_null" / "v2" / "sdss_audit" / "re_cut"

# (lo inclusive, hi exclusive or inf, label)
DEFAULT_RE_SLICES: tuple[tuple[float, float, str], ...] = (
    (1.0, 2.0, "1–2"),
    (2.0, 3.0, "2–3"),
    (3.0, float("inf"), "3+"),
)


def load_strict_pool(path: Path) -> pd.DataFrame:
    print(f"[*] Loading {path} ...")
    df = pd.read_csv(path, usecols=list(USECOLS_STRICT))
    df = ensure_sdss_colors(df)
    mag = pd.to_numeric(df["modelMag_r"], errors="coerce")
    ba = pd.to_numeric(df["expAB_r"], errors="coerce")
    re = pd.to_numeric(df["expRad_r"], errors="coerce")
    ok = mag.notna() & ba.notna() & re.notna() & (ba >= 0) & (ba <= 1) & (re > 0)
    out = df.loc[ok].copy()
    out["modelMag_r"] = mag[ok].to_numpy()
    out["expAB_r"] = ba[ok].to_numpy()
    out["expRad_r"] = re[ok].to_numpy()
    out = apply_strict_cuts(out)
    print(
        f"[*] Strict pool N={len(out):,} "
        f"(u-r<{SDSS_UR_MAX_CDF:g}, lnL exp-wins, expAB_r>{Q0:g})"
    )
    return out


def slice_re(df: pd.DataFrame, lo: float, hi: float) -> pd.DataFrame:
    re = df["expRad_r"].to_numpy(dtype=float)
    if np.isfinite(hi):
        mask = (re >= lo) & (re < hi)
    else:
        mask = re >= lo
    return df.loc[mask].copy()


def parse_re_slices(spec: str) -> tuple[tuple[float, float, str], ...]:
    """Parse '1:2,2:3,3:inf' into slice tuples."""
    slices: list[tuple[float, float, str]] = []
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        lo_s, hi_s = part.split(":")
        lo = float(lo_s)
        hi = float("inf") if hi_s.lower() in ("inf", "infinity") else float(hi_s)
        label = f"{lo:g}–{hi_s}" if np.isfinite(hi) else f"{lo:g}+"
        slices.append((lo, hi, label))
    return tuple(slices)


def plot_slices(
    tables: dict[str, pd.DataFrame],
    pool_counts: dict[str, int],
    *,
    full_bins: pd.DataFrame | None,
    n_strict: int,
    out_png: Path,
) -> None:
    colors = ("#377eb8", "#ff7f00", "#e41a1c", "#984ea3", "#4daf4a")
    fig, (ax_ba, ax_n) = plt.subplots(
        2,
        1,
        figsize=(10, 8),
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1], "hspace": 0.08},
    )

    if full_bins is not None and len(full_bins):
        x = 0.5 * (full_bins["mag_lo"] + full_bins["mag_hi"])
        ax_ba.plot(
            x,
            full_bins["median_expAB_r"],
            "--",
            color="0.45",
            lw=1.8,
            label=f"All strict (N={n_strict:,})",
            zorder=1,
        )

    for i, (label, bins) in enumerate(tables.items()):
        if bins.empty:
            continue
        x = 0.5 * (bins["mag_lo"] + bins["mag_hi"])
        c = colors[i % len(colors)]
        ax_ba.plot(
            x,
            bins["median_expAB_r"],
            "o-",
            color=c,
            lw=2,
            ms=5,
            label=f"Re {label}\" (N={pool_counts[label]:,})",
            zorder=2 + i,
        )
        ax_n.plot(x, bins["n"], "o-", color=c, lw=1.5, ms=4, label=f"Re {label}\"")

    ax_ba.set_ylabel("median expAB$_r$")
    ax_ba.set_ylim(0, 1.02)
    ax_ba.legend(loc="upper right", fontsize=8)
    ax_ba.grid(True, alpha=0.3)
    ax_ba.set_title(
        "Strict pool: median b/a vs mag within fixed Re slices\n"
        f"u-r < {SDSS_UR_MAX_CDF:g}, lnL exp-wins, expAB$_r$ > {Q0:g}  |  "
        "Re = PhotoObj expRad$_r$ (arcsec)"
    )

    ax_n.set_xlabel(r"modelMag$_r$ (0.5 mag bins)")
    ax_n.set_ylabel(r"$N$ per mag bin")
    ax_n.set_yscale("log")
    ax_n.grid(True, alpha=0.3, which="both")
    ax_n.legend(loc="upper right", fontsize=7)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(hspace=0.12)
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sdss-csv", type=Path, default=DEFAULT_SDSS_V2)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument(
        "--re-slices",
        default="1:2,2:3,3:inf",
        help="Re bin edges as lo:hi,... with hi=inf for open top (default: 1:2,2:3,3:inf).",
    )
    parser.add_argument("--min-n", type=int, default=25, help="Min galaxies per mag bin.")
    parser.add_argument(
        "--no-full-reference",
        action="store_true",
        help="Omit all-strict reference curve.",
    )
    args = parser.parse_args()

    re_slices = parse_re_slices(args.re_slices)
    pool = load_strict_pool(args.sdss_csv)

    re_vals = pool["expRad_r"].to_numpy(dtype=float)
    print(
        "[*] Strict pool expRad_r percentiles: "
        f"p25={np.percentile(re_vals, 25):.3f}\", "
        f"p50={np.percentile(re_vals, 50):.3f}\", "
        f"p75={np.percentile(re_vals, 75):.3f}\", "
        f"p90={np.percentile(re_vals, 90):.3f}\""
    )

    tables: dict[str, pd.DataFrame] = {}
    pool_counts: dict[str, int] = {}
    long_rows: list[pd.DataFrame] = []

    for lo, hi, label in re_slices:
        sub = slice_re(pool, lo, hi)
        n_sub = len(sub)
        pool_counts[label] = n_sub
        pct = 100.0 * n_sub / max(1, len(pool))
        hi_str = f"{hi:g}" if np.isfinite(hi) else "inf"
        print(f"[*] Re [{lo:g}, {hi_str}) arcsec: N={n_sub:,} ({pct:.1f}% of strict pool)")

        mag = sub["modelMag_r"].to_numpy(dtype=float)
        ba = sub["expAB_r"].to_numpy(dtype=float)
        bins = mag_bin_table(mag, ba, pool_n=n_sub, min_n=args.min_n)
        bins.insert(0, "re_slice", label)
        bins.insert(1, "re_lo", lo)
        bins.insert(2, "re_hi", hi if np.isfinite(hi) else np.nan)
        bins.insert(3, "slice_n", n_sub)
        tables[label] = bins
        long_rows.append(bins)

    full_bins = None
    if not args.no_full_reference:
        mag_all = pool["modelMag_r"].to_numpy(dtype=float)
        ba_all = pool["expAB_r"].to_numpy(dtype=float)
        full_bins = mag_bin_table(mag_all, ba_all, pool_n=len(pool), min_n=args.min_n)

    summary = pd.concat(long_rows, ignore_index=True)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.out_dir / "re_slice_mag_ba_summary.csv"
    summary.to_csv(csv_path, index=False)
    print(f"[*] Wrote {csv_path}")

    png_path = args.out_dir / "re_slice_mag_ba_panel.png"
    plot_slices(
        tables,
        pool_counts,
        full_bins=full_bins,
        n_strict=len(pool),
        out_png=png_path,
    )
    print(f"[*] Saved {png_path}")

    print("\n--- Per-slice median b/a at mag 20–21 (where comparable) ---")
    for label, bins in tables.items():
        sub = bins[(bins["mag_lo"] >= 20.0) & (bins["mag_hi"] <= 21.0)]
        if sub.empty:
            print(f"  Re {label}\": (no bins with n>={args.min_n})")
        else:
            print(
                f"  Re {label}\": "
                + ", ".join(
                    f"{r.mag_lo:g}–{r.mag_hi:g} med b/a={r.median_expAB_r:.3f} (n={int(r.n)})"
                    for r in sub.itertuples()
                )
            )


if __name__ == "__main__":
    main()
