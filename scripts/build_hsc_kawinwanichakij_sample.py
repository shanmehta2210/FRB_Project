"""
Build a random sample (up to 500k) from Kawinwanichakij+2021 HSC structural catalogs.

Inputs (downloaded from https://member.ipmu.jp/john.silverman/HSC/):
  - Kawinwanichakij2021_dud_pdr2_final.fits
  - Kawinwanichakij2021_wide_efeds_pdr2_final.fits

Run from repo root::

    python scripts/build_hsc_kawinwanichakij_sample.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.table import Table

DEFAULT_DIR = Path("catalog_downloads/kawinwanichakij")
DEFAULT_OUT = "catalog/HSC_kawinwanichakij_sample_500k.csv"
DEFAULT_TARGET = 500_000
DEFAULT_SEED = 42

FILES = (
    "Kawinwanichakij2021_dud_pdr2_final.fits",
    "Kawinwanichakij2021_wide_efeds_pdr2_final.fits",
)

# Prefer these columns when present (names vary slightly across releases).
PREFERRED = [
    "object_id",
    "Object_id",
    "ra",
    "dec",
    "ira",
    "idec",
    "fitted_q",
    "fitted_sersic",
    "fitted_reff",
    "fitted_mag",
    "fitted_flux",
    "fitted_redchi2",
    "corrected_q",
    "corrected_sersic",
    "corrected_reff",
    "corrected_mag",
    "rmag",
    "imag",
    "gmag",
    "zmag",
    "ymag",
    "photoz_best",
    "stellar_mass",
    "goodfits_flag",
    "use_flag",
    "calib_flag",
    "quiescent_flag",
]


def load_one(path: Path) -> pd.DataFrame:
    print(f"[*] Reading {path} ...", flush=True)
    tbl = Table.read(path, format="fits")
    df = tbl.to_pandas()
    # Normalize column names to lower where helpful for object_id
    rename = {}
    if "Object_id" in df.columns and "object_id" not in df.columns:
        rename["Object_id"] = "object_id"
    if rename:
        df = df.rename(columns=rename)
    df["source_file"] = path.name
    print(f"    rows={len(df):,} cols={len(df.columns)}", flush=True)
    return df


def select_columns(df: pd.DataFrame) -> pd.DataFrame:
    keep = [c for c in PREFERRED if c in df.columns]
    # Always keep identifiers / coords / key morph if present under any case
    extra = []
    for c in df.columns:
        cl = c.lower()
        if cl in {
            "object_id",
            "ra",
            "dec",
            "fitted_q",
            "fitted_sersic",
            "fitted_reff",
            "fitted_mag",
            "corrected_q",
            "corrected_sersic",
            "corrected_reff",
            "goodfits_flag",
            "use_flag",
        }:
            if c not in keep:
                extra.append(c)
    cols = list(dict.fromkeys(keep + extra + ["source_file"]))
    # If still too few morph cols, keep all
    morphish = [c for c in cols if "fitted" in c.lower() or "corrected" in c.lower()]
    if len(morphish) < 3:
        print("[!] Few morph columns matched; writing full column set.", flush=True)
        return df
    return df[cols].copy()


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--data-dir", type=Path, default=DEFAULT_DIR)
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    p.add_argument("--target-rows", type=int, default=DEFAULT_TARGET)
    p.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = p.parse_args()

    frames: list[pd.DataFrame] = []
    for name in FILES:
        path = args.data_dir / name
        if not path.is_file():
            print(f"[!] Missing {path}", flush=True)
            continue
        frames.append(load_one(path))

    if not frames:
        raise FileNotFoundError(f"No Kawinwanichakij FITS found in {args.data_dir}")

    cat = pd.concat(frames, ignore_index=True)
    # Dedup on object_id if present
    id_col = "object_id" if "object_id" in cat.columns else None
    if id_col:
        before = len(cat)
        cat = cat.drop_duplicates(subset=[id_col], keep="first")
        print(f"[*] Dedup {id_col}: {before:,} -> {len(cat):,}", flush=True)

    cat = select_columns(cat)

    # Standard aliases for our pipeline
    if "fitted_q" in cat.columns:
        cat["ba"] = pd.to_numeric(cat["fitted_q"], errors="coerce")
    if "corrected_q" in cat.columns:
        cat["ba_corr"] = pd.to_numeric(cat["corrected_q"], errors="coerce")
    if "fitted_sersic" in cat.columns:
        cat["n_sersic"] = pd.to_numeric(cat["fitted_sersic"], errors="coerce")
    if "fitted_reff" in cat.columns:
        cat["re_arcsec"] = pd.to_numeric(cat["fitted_reff"], errors="coerce")
    if "fitted_mag" in cat.columns:
        cat["mag"] = pd.to_numeric(cat["fitted_mag"], errors="coerce")
    for src, dst in (("ra", "RA_ICRS"), ("ira", "RA_ICRS"), ("dec", "DE_ICRS"), ("idec", "DE_ICRS")):
        if src in cat.columns and dst not in cat.columns:
            cat[dst] = pd.to_numeric(cat[src], errors="coerce")

    n = len(cat)
    if n > args.target_rows:
        cat = cat.sample(n=args.target_rows, random_state=args.seed).reset_index(drop=True)
        print(f"[*] Random sample {args.target_rows:,} of {n:,}", flush=True)
    else:
        print(f"[*] Catalog has only {n:,} rows (< {args.target_rows:,}); writing all.", flush=True)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    cat.to_csv(args.out, index=False)
    print(f"Wrote {args.out} ({len(cat):,} rows, {len(cat.columns)} cols)", flush=True)
    if "ba" in cat.columns:
        ba = pd.to_numeric(cat["ba"], errors="coerce")
        print(f"  median fitted_q/ba={ba.median():.4f}", flush=True)
    if "RA_ICRS" in cat.columns:
        print(
            f"  RA=[{cat['RA_ICRS'].min():.3f},{cat['RA_ICRS'].max():.3f}]  "
            f"Dec=[{cat['DE_ICRS'].min():.3f},{cat['DE_ICRS'].max():.3f}]",
            flush=True,
        )


if __name__ == "__main__":
    main()
