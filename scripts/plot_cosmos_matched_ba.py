#!/usr/bin/env python3
"""Pass-2 matched-sample Δ(b/a) audit: HST Zurich vs SDSS in COSMOS ACS footprint.

Cross-match radius 0.5 arcsec. Plot SDSS expAB_r − HST b_a vs HST ACS_MAG_AUTO.

Run from repo root::

    python scripts/plot_cosmos_matched_ba.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord, match_coordinates_sky
import astropy.units as u

_SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS))

from cosmos_ba_audit_utils import (  # noqa: E402
    COSMOS_BA_MIN,
    COSMOS_PLOTS,
    HST_BA_COL,
    HST_CSV,
    HST_MAG_COL,
    MIN_N_BIN,
    SDSS_BA_COL,
    SDSS_CSV,
    SDSS_MAG_COL,
    mag_bin_table,
)

MATCH_RADIUS_ARCSEC = 0.5


def cross_match(hst: pd.DataFrame, sdss: pd.DataFrame, radius_arcsec: float) -> pd.DataFrame:
    h_ra = pd.to_numeric(hst["ra"], errors="coerce").to_numpy(dtype=float)
    h_dec = pd.to_numeric(hst["dec"], errors="coerce").to_numpy(dtype=float)
    s_ra = pd.to_numeric(sdss["ra"], errors="coerce").to_numpy(dtype=float)
    s_dec = pd.to_numeric(sdss["dec"], errors="coerce").to_numpy(dtype=float)

    h_ok = np.isfinite(h_ra) & np.isfinite(h_dec)
    s_ok = np.isfinite(s_ra) & np.isfinite(s_dec)

    h_coords = SkyCoord(ra=h_ra[h_ok] * u.deg, dec=h_dec[h_ok] * u.deg)
    s_coords = SkyCoord(ra=s_ra[s_ok] * u.deg, dec=s_dec[s_ok] * u.deg)

    idx, sep, _ = match_coordinates_sky(h_coords, s_coords)
    match_mask = sep <= radius_arcsec * u.arcsec

    h_idx = np.where(h_ok)[0][match_mask]
    s_idx = np.where(s_ok)[0][idx[match_mask].astype(int)]
    sep_matched = sep[match_mask].arcsec

    rows = []
    for hi, si, sep_arcsec in zip(h_idx, s_idx, sep_matched):
        rows.append(
            {
                "hst_ra": h_ra[hi],
                "hst_dec": h_dec[hi],
                "sdss_ra": s_ra[si],
                "sdss_dec": s_dec[si],
                "sep_arcsec": float(sep_arcsec),
                HST_MAG_COL: pd.to_numeric(hst.iloc[hi][HST_MAG_COL], errors="coerce"),
                HST_BA_COL: pd.to_numeric(hst.iloc[hi][HST_BA_COL], errors="coerce"),
                SDSS_MAG_COL: pd.to_numeric(sdss.iloc[si][SDSS_MAG_COL], errors="coerce"),
                SDSS_BA_COL: pd.to_numeric(sdss.iloc[si][SDSS_BA_COL], errors="coerce"),
            }
        )

    matched = pd.DataFrame(rows)
    if matched.empty:
        return matched

    matched["delta_b_a"] = matched[SDSS_BA_COL] - matched[HST_BA_COL]
    return matched


def plot_delta_ba(matched: pd.DataFrame, *, out_png: Path, out_csv: Path, min_n: int) -> None:
    mag = matched[HST_MAG_COL].to_numpy(dtype=float)
    delta = matched["delta_b_a"].to_numpy(dtype=float)

    bins = mag_bin_table(mag, delta, pool_n=len(matched), mag_col_label=HST_MAG_COL, min_n=min_n)
    summary = bins.copy()
    summary["median_delta_b_a"] = summary["median_b_a"]
    summary["mean_delta_b_a"] = summary["mean_b_a"]
    summary = summary.drop(columns=["median_b_a", "mean_b_a"])
    summary.to_csv(out_csv, index=False)

    fig, (ax_main, ax_n) = plt.subplots(
        2,
        1,
        figsize=(9, 7),
        sharex=True,
        gridspec_kw={"height_ratios": [2.2, 1], "hspace": 0.08},
    )

    ok = np.isfinite(mag) & np.isfinite(delta)
    ax_main.scatter(mag[ok], delta[ok], s=8, alpha=0.15, color="0.5", label="matched pairs")
    if len(bins):
        x = 0.5 * (bins["mag_lo"] + bins["mag_hi"])
        ax_main.plot(x, bins["median_b_a"], "o-", color="#377eb8", lw=2, ms=5, label="median Δ(b/a)")
    ax_main.axhline(0, color="k", ls="--", lw=1)
    ax_main.set_ylabel(r"$\Delta(b/a)$ = SDSS $-$ HST")
    ax_main.legend(loc="upper right", fontsize=8)
    ax_main.grid(True, alpha=0.3)
    ax_main.set_title(
        f"COSMOS matched sample (r <= {MATCH_RADIUS_ARCSEC} arcsec, N={len(matched):,})\n"
        f"b/a > {COSMOS_BA_MIN:g}, no colour cut"
    )

    if len(bins):
        ax_n.plot(x, bins["n"], "o-", color="#377eb8", lw=1.5, ms=4)
    ax_n.set_xlabel(f"HST {HST_MAG_COL} (0.5 mag bins)")
    ax_n.set_ylabel(r"$N$ per mag bin")
    ax_n.set_yscale("log")
    ax_n.grid(True, alpha=0.3, which="both")

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(hspace=0.12)
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hst-csv", type=Path, default=HST_CSV)
    parser.add_argument("--sdss-csv", type=Path, default=SDSS_CSV)
    parser.add_argument("--out-dir", type=Path, default=COSMOS_PLOTS)
    parser.add_argument("--radius-arcsec", type=float, default=MATCH_RADIUS_ARCSEC)
    parser.add_argument("--min-n", type=int, default=MIN_N_BIN)
    args = parser.parse_args()

    hst = pd.read_csv(args.hst_csv)
    sdss = pd.read_csv(args.sdss_csv)
    matched = cross_match(hst, sdss, args.radius_arcsec)
    print(f"[*] Matched N={len(matched):,} within {args.radius_arcsec} arcsec")

    if matched.empty:
        print("[!] No matches found")
        return

    args.out_dir.mkdir(parents=True, exist_ok=True)
    matched.to_csv(args.out_dir / "matched_pairs.csv", index=False)
    plot_delta_ba(
        matched,
        out_png=args.out_dir / "matched_delta_ba_vs_mag.png",
        out_csv=args.out_dir / "matched_ba_delta_summary.csv",
        min_n=args.min_n,
    )
    print(f"[*] Wrote matched plots to {args.out_dir}")


if __name__ == "__main__":
    main()
