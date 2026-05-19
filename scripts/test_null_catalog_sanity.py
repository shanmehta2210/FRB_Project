"""
Sanity checks for v1 null catalogs (Legacy + SDSS).

Exit 0 if both strict pools have >= min_pool rows and basic column/footprint checks pass.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (
    JOINT_DEC_MAX,
    JOINT_DEC_MIN,
    MAG_LIMIT,
    Q0,
    footprint_summary,
    prepare_null_sample,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_LEGACY = REPO_ROOT / "LS_catalog_v1_allsky_modelmr.csv"
DEFAULT_SDSS = REPO_ROOT / "SDSS_catalog_v1_allsky_modelmr.csv"


def _check_catalog(path: Path, label: str, mag_col: str, is_legacy: bool, min_pool: int) -> bool:
    if not path.is_file():
        print(f"[FAIL] {label}: missing {path}")
        return False

    df = pd.read_csv(path)
    footprint_summary(df, label)

    required = {"RA_ICRS", "DE_ICRS", "expAB_r", mag_col}
    missing = required - set(df.columns)
    if missing:
        print(f"[FAIL] {label}: missing columns {missing}")
        return False

    ra = pd.to_numeric(df["RA_ICRS"], errors="coerce")
    dec = pd.to_numeric(df["DE_ICRS"], errors="coerce")
    if not (ra.notna().all() and dec.notna().all()):
        print(f"[FAIL] {label}: non-finite RA/Dec")
        return False

    dec_ok = (dec >= JOINT_DEC_MIN - 0.5) & (dec <= JOINT_DEC_MAX + 0.5)
    if dec_ok.mean() < 0.99:
        print(f"[FAIL] {label}: only {dec_ok.mean():.1%} rows in joint Dec footprint")
        return False

    q = pd.to_numeric(df["expAB_r"], errors="coerce")
    q_ok = np.isfinite(q) & (q >= 0) & (q <= 1)
    if q_ok.mean() < 0.99:
        print(f"[FAIL] {label}: only {q_ok.mean():.1%} rows with expAB_r in [0,1]")
        return False

    if is_legacy and "tractor_mag_r" in df.columns and "rmag" in df.columns:
        diff = (pd.to_numeric(df["tractor_mag_r"]) - pd.to_numeric(df["rmag"])).abs()
        med = float(diff.median())
        print(f"[OK] {label}: median |tractor_mag_r - rmag| = {med:.6f}")
        if med > 1e-6:
            print(f"[WARN] {label}: tractor_mag_r and rmag differ (expected identical)")

    strict = prepare_null_sample(
        df,
        sample_mode="strict",
        mag_column=mag_col,
        mag_limit=MAG_LIMIT,
        q0=Q0,
        is_legacy=is_legacy,
    )
    print(f"[OK] {label}: strict pool n={len(strict)} (min {min_pool})")
    if len(strict) < min_pool:
        print(f"[FAIL] {label}: strict pool too small")
        return False

    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--legacy-csv", default=str(DEFAULT_LEGACY))
    parser.add_argument("--sdss-csv", default=str(DEFAULT_SDSS))
    parser.add_argument("--min-pool", type=int, default=10_000)
    args = parser.parse_args()

    ok_legacy = _check_catalog(
        Path(args.legacy_csv),
        "Legacy",
        mag_col="tractor_mag_r",
        is_legacy=True,
        min_pool=args.min_pool,
    )
    ok_sdss = _check_catalog(
        Path(args.sdss_csv),
        "SDSS",
        mag_col="rmag",
        is_legacy=False,
        min_pool=args.min_pool,
    )

    if ok_legacy and ok_sdss:
        print("[PASS] All null catalog sanity checks passed.")
        return 0
    print("[FAIL] One or more checks failed.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
