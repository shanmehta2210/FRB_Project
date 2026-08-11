#!/usr/bin/env python3
"""Merge lnL patch cache into SDSS v1 CSV with robust sky matching."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import assign_sdss_profile_winner_columns, sdss_exp_wins_lnl_mask

COORD_DECIMALS = 6
MAG_DECIMALS = 3
LNL_COLS = ("lnLDeV_r", "lnLExp_r", "deVMag_r", "expMag_r", "modelMag_r")


def _add_keys(
    df: pd.DataFrame,
    *,
    ra_col: str,
    dec_col: str,
    mag_col: str | None,
) -> pd.DataFrame:
    out = df.copy()
    out["_ra_k"] = pd.to_numeric(out[ra_col], errors="coerce").round(COORD_DECIMALS)
    out["_dec_k"] = pd.to_numeric(out[dec_col], errors="coerce").round(COORD_DECIMALS)
    if mag_col:
        out["_mag_k"] = pd.to_numeric(out[mag_col], errors="coerce").round(MAG_DECIMALS)
    return out


def merge_patch(base: pd.DataFrame, patch: pd.DataFrame) -> pd.DataFrame:
    """Two-pass merge: (ra,dec,cmodelMag) then (ra,dec) for remaining."""
    for c in LNL_COLS:
        if c in patch.columns:
            patch[c] = pd.to_numeric(patch[c], errors="coerce")

    b = _add_keys(base, ra_col="RA_ICRS", dec_col="DE_ICRS", mag_col="rmag")
    p = _add_keys(patch, ra_col="ra", dec_col="dec", mag_col="cmodelMag_r")

    keep_cols = ["_ra_k", "_dec_k", "_mag_k"] + [c for c in LNL_COLS if c in p.columns]
    p3 = p[keep_cols].drop_duplicates(subset=["_ra_k", "_dec_k", "_mag_k"], keep="first")

    out = b.merge(p3, on=["_ra_k", "_dec_k", "_mag_k"], how="left")

    miss = out["lnLExp_r"].isna() if "lnLExp_r" in out.columns else pd.Series(True, index=out.index)
    if miss.any():
        p2 = p[["_ra_k", "_dec_k"] + [c for c in LNL_COLS if c in p.columns]].drop_duplicates(
            subset=["_ra_k", "_dec_k"], keep="first"
        )
        fill = b.loc[miss, ["_ra_k", "_dec_k"]].merge(p2, on=["_ra_k", "_dec_k"], how="left")
        for c in LNL_COLS:
            if c in fill.columns:
                out.loc[miss, c] = fill[c].to_numpy()

    return out.drop(columns=["_ra_k", "_dec_k", "_mag_k"], errors="ignore")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--in-csv", default="catalog/SDSS_catalog_v1_allsky_modelmr.csv")
    parser.add_argument("--patch-csv", default="catalog/SDSS_lnl_patch_cache.csv")
    parser.add_argument("--out-csv", default="catalog/SDSS_catalog_v1_allsky_modelmr.csv")
    args = parser.parse_args()

    base = pd.read_csv(args.in_csv)
    patch = pd.read_csv(args.patch_csv)
    merged = merge_patch(base, patch)
    n_lnl = int(pd.to_numeric(merged["lnLExp_r"], errors="coerce").notna().sum())
    print(f"Matched lnL: {n_lnl} / {len(base)} ({n_lnl / len(base):.1%})")

    merged = assign_sdss_profile_winner_columns(merged)
    finite = pd.to_numeric(merged["lnLExp_r"], errors="coerce").notna() & pd.to_numeric(
        merged["lnLDeV_r"], errors="coerce"
    ).notna()
    wins = sdss_exp_wins_lnl_mask(merged)
    print(
        f"exp-winner fraction (finite lnL only): "
        f"{wins.sum() / max(int(finite.sum()), 1):.3f}"
    )

    merged.to_csv(args.out_csv, index=False)
    print(f"Wrote {args.out_csv}")


if __name__ == "__main__":
    main()
