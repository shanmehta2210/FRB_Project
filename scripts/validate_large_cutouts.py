#!/usr/bin/env python3
"""Sanity-check large_cutouts flux+invvar before running the pipeline."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

REPO = Path(__file__).resolve().parents[1]
CUTOUT_DIR = REPO / "large_cutouts"
LOC_CSV = REPO / "master_frb_localization.csv"

# Legacy-like nanomaggie cutouts (53 original + DR10)
FLUX_MED_MIN = 1e-7
FLUX_MED_MAX = 1.0
INVVAR_MED_MIN = 1e3
INVVAR_GOOD_FRAC = 0.5


def validate_pair(frb: str) -> dict:
    flux_p = CUTOUT_DIR / f"{frb}_flux.fits"
    inv_p = CUTOUT_DIR / f"{frb}_invvar.fits"
    row = {"frb": frb, "ok": False, "issues": []}
    if not flux_p.is_file() or not inv_p.is_file():
        row["issues"].append("missing_flux_or_invvar")
        return row
    with fits.open(flux_p) as hf, fits.open(inv_p) as hi:
        flux = np.squeeze(hf[0].data).astype(float)
        inv = np.squeeze(hi[0].data).astype(float)
        row["shape"] = flux.shape
        row["flux_med"] = float(np.nanmedian(flux))
        row["inv_med"] = float(np.nanmedian(inv[inv > 0])) if np.any(inv > 0) else 0.0
        row["inv_good_frac"] = float(np.mean(inv > 0))
    if flux.shape != (2290, 2290):
        row["issues"].append(f"bad_shape_{flux.shape}")
    if not np.isfinite(row["flux_med"]) or row["flux_med"] < FLUX_MED_MIN or row["flux_med"] > FLUX_MED_MAX:
        row["issues"].append(f"flux_median_out_of_range({row['flux_med']:.3g})")
    if row["inv_good_frac"] < INVVAR_GOOD_FRAC:
        row["issues"].append(f"invvar_mostly_zero(frac={row['inv_good_frac']:.3f})")
    if row["inv_med"] < INVVAR_MED_MIN:
        row["issues"].append(f"invvar_median_low({row['inv_med']:.3g})")
    row["ok"] = len(row["issues"]) == 0
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frb", nargs="*", help="Only these FRBs")
    parser.add_argument("--list", type=Path)
    parser.add_argument("--csv", type=Path, default=LOC_CSV)
    args = parser.parse_args()

    names = list(args.frb or [])
    if args.list and args.list.is_file():
        names.extend(
            ln.strip()
            for ln in args.list.read_text().splitlines()
            if ln.strip() and not ln.startswith("#")
        )
    if not names:
        loc = pd.read_csv(args.csv)
        names = loc["frb"].tolist()

    rows = [validate_pair(f) for f in names]
    df = pd.DataFrame(rows)
    out = CUTOUT_DIR / "cutout_validation.csv"
    df.to_csv(out, index=False)
    n_ok = int(df["ok"].sum())
    print(f"[validate] {n_ok}/{len(df)} passed -> {out}")
    bad = df[~df["ok"]]
    if len(bad):
        print(bad[["frb", "issues", "flux_med", "inv_med", "inv_good_frac"]].to_string(index=False))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
