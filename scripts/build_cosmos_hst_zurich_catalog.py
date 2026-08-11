#!/usr/bin/env python3
"""Build HST pools from COSMOS Zurich morphology catalog (GIM2D).

Writes both entire (no b/a floor) and strict (b/a > 0.2) CSVs.

Run from repo root::

    python scripts/build_cosmos_hst_zurich_catalog.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS))

from cosmos_ba_audit_utils import (  # noqa: E402
    COSMOS_BA_MIN,
    GIM2D_DISK_TYPE,
    HST_CSV,
    HST_DISK_CSV,
    HST_DISK_ENTIRE_CSV,
    HST_ENTIRE_CSV,
    ZURICH_TBL,
    in_footprint_mask,
)

SENTINEL = -999999.0

HST_OUT_COLS = (
    "SequentialID",
    "ra",
    "dec",
    "ACS_MAG_AUTO",
    "b_a",
    "Re_arcsec",
    "ELL_GIM2D",
    "TYPE",
    "STELLARITY",
)


def _finite_series(s: pd.Series) -> pd.Series:
    v = pd.to_numeric(s, errors="coerce")
    return v.replace(SENTINEL, np.nan)


def read_zurich(path: Path) -> pd.DataFrame:
    print(f"[*] Reading {path} ...")
    with path.open(encoding="utf-8", errors="replace") as f:
        header = [c.strip() for c in f.readline().strip("| \n").split("|") if c.strip()]
        for _ in range(3):
            f.readline()
        df = pd.read_csv(f, sep=r"\s+", names=header, engine="python")
    df = df.rename(columns={"RA": "ra", "DEC": "dec"})
    print(f"[*] Loaded N={len(df):,} Zurich rows")
    return df


def apply_zurich_cuts(
    df: pd.DataFrame, *, strict_ba: bool, disk_only: bool
) -> tuple[pd.DataFrame, list[dict]]:
    funnel: list[dict] = []

    ra = _finite_series(df["ra"]).to_numpy(dtype=float)
    dec = _finite_series(df["dec"]).to_numpy(dtype=float)
    fp = in_footprint_mask(ra, dec)
    out = df.loc[fp].copy()
    funnel.append({"stage": "footprint_acs_box", "n_remaining": len(out)})

    mu = _finite_series(out["ACS_MU_CLASS"])
    out = out.loc[mu == 1].copy()
    funnel.append({"stage": "ACS_MU_CLASS_galaxy", "n_remaining": len(out)})

    st = _finite_series(out["STELLARITY"])
    out = out.loc[st == 0].copy()
    funnel.append({"stage": "STELLARITY_galaxy", "n_remaining": len(out)})

    jf = _finite_series(out["JUNKFLAG"])
    out = out.loc[jf == 0].copy()
    funnel.append({"stage": "JUNKFLAG_clean", "n_remaining": len(out)})

    mag = _finite_series(out["ACS_MAG_AUTO"])
    flux = _finite_series(out["FLUX_GIM2D"])
    ell = _finite_series(out["ELL_GIM2D"])

    gim2d_ok = mag.notna() & (flux > 0) & ell.notna() & (ell >= 0) & (ell < 1)
    out = out.loc[gim2d_ok].copy()
    funnel.append({"stage": "finite_gim2d_mag_shape", "n_remaining": len(out)})

    out["ACS_MAG_AUTO"] = _finite_series(out["ACS_MAG_AUTO"])
    out["ELL_GIM2D"] = _finite_series(out["ELL_GIM2D"])
    out["b_a"] = 1.0 - out["ELL_GIM2D"]
    out["Re_arcsec"] = _finite_series(out["R_0P5_GIM2D"])
    out["ra"] = _finite_series(out["ra"])
    out["dec"] = _finite_series(out["dec"])

    out = out.loc[out["b_a"].notna() & (out["b_a"] > 0) & (out["b_a"] <= 1.0)].copy()
    funnel.append({"stage": "finite_b_a", "n_remaining": len(out)})

    if disk_only:
        typ = _finite_series(out["TYPE"])
        out = out.loc[typ == GIM2D_DISK_TYPE].copy()
        funnel.append({"stage": f"TYPE_eq_{GIM2D_DISK_TYPE}_disk", "n_remaining": len(out)})

    if strict_ba:
        out = out.loc[out["b_a"] > COSMOS_BA_MIN].copy()
        funnel.append({"stage": f"b_a_gt_{COSMOS_BA_MIN:g}", "n_remaining": len(out)})

    return out.reset_index(drop=True), funnel


def write_pool(pool: pd.DataFrame, path: Path) -> None:
    keep = [c for c in HST_OUT_COLS if c in pool.columns]
    path.parent.mkdir(parents=True, exist_ok=True)
    pool[keep].to_csv(path, index=False)
    print(f"[*] Wrote {path} (N={len(pool):,})")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--zurich-tbl", type=Path, default=ZURICH_TBL)
    parser.add_argument("--out-entire", type=Path, default=HST_ENTIRE_CSV)
    parser.add_argument("--out-strict", type=Path, default=HST_CSV)
    args = parser.parse_args()

    raw = read_zurich(args.zurich_tbl)

    entire, funnel_entire = apply_zurich_cuts(raw, strict_ba=False, disk_only=False)
    strict, funnel_strict = apply_zurich_cuts(raw, strict_ba=True, disk_only=False)
    disk_entire, funnel_de = apply_zurich_cuts(raw, strict_ba=False, disk_only=True)
    disk_strict, funnel_ds = apply_zurich_cuts(raw, strict_ba=True, disk_only=True)

    write_pool(entire, args.out_entire)
    write_pool(strict, args.out_strict)
    write_pool(disk_entire, HST_DISK_ENTIRE_CSV)
    write_pool(disk_strict, HST_DISK_CSV)

    pd.DataFrame(funnel_entire).to_csv(
        args.out_entire.with_name("cosmos_hst_entire_cut_funnel.csv"), index=False
    )
    pd.DataFrame(funnel_strict).to_csv(
        args.out_strict.with_name("cosmos_hst_cut_funnel.csv"), index=False
    )
    pd.DataFrame(funnel_de).to_csv(
        HST_DISK_ENTIRE_CSV.with_name("cosmos_hst_disk_entire_cut_funnel.csv"), index=False
    )
    pd.DataFrame(funnel_ds).to_csv(
        HST_DISK_CSV.with_name("cosmos_hst_disk_cut_funnel.csv"), index=False
    )

    print(
        f"[*] disk entire N={len(disk_entire):,}, disk strict N={len(disk_strict):,}; "
        f"median ACS_MAG_AUTO disk entire={disk_entire['ACS_MAG_AUTO'].median():.3f}"
    )


if __name__ == "__main__":
    main()
