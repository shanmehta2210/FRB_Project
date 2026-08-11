"""

Sanity checks for v1/v2 null catalogs (Legacy + SDSS).



Exit 0 if pools meet size thresholds and basic column/footprint checks pass.

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

    LEGACY_CDF_TYPE_EXCLUDE,

    LEGACY_GR_MAX_CDF,

    MAG_LIMIT,

    Q0,

    SDSS_MAG20_LIMIT,

    SDSS_MIN_STRICT_MAG20_POOL,

    SDSS_UR_MAX_CDF,

    footprint_summary,

    prepare_null_strict_color_base,

    slice_null_base_by_mag,

)



REPO_ROOT = Path(__file__).resolve().parents[1]

DEFAULT_LEGACY = REPO_ROOT / "catalog/LS_catalog_v1_allsky_modelmr.csv"

DEFAULT_SDSS = REPO_ROOT / "catalog/SDSS_catalog_v1_allsky_modelmr.csv"

DEFAULT_SDSS_V2 = REPO_ROOT / "catalog/SDSS_catalog_v2_fullsky_modelmr.csv"

# Post-morphology Legacy CDF pool at m<21 is ~4.4k (see legacy_morphology diagnostic).

LEGACY_MIN_STRICT_POOL = 4_000





def _check_sdss_v2(path: Path, label: str) -> bool:

    if not path.is_file():

        print(f"[FAIL] {label}: missing {path}")

        return False



    df = pd.read_csv(path)

    footprint_summary(df, label)



    required = {"objID", "RA_ICRS", "DE_ICRS", "expAB_r", "modelMag_r", "lnLExp_r"}

    missing = required - set(df.columns)

    if missing:

        print(f"[FAIL] {label}: missing columns {missing}")

        return False



    if df["objID"].duplicated().any():

        n_dup = int(df["objID"].duplicated().sum())

        print(f"[FAIL] {label}: {n_dup} duplicate objID values")

        return False

    print(f"[OK] {label}: objID unique ({len(df):,} rows)")



    ra = pd.to_numeric(df["RA_ICRS"], errors="coerce")

    dec = pd.to_numeric(df["DE_ICRS"], errors="coerce")

    if not (ra.notna().all() and dec.notna().all()):

        print(f"[FAIL] {label}: non-finite RA/Dec")

        return False



    # v2 should extend beyond joint footprint (southern stripes).

    frac_joint = ((dec >= JOINT_DEC_MIN) & (dec <= JOINT_DEC_MAX)).mean()

    print(f"[OK] {label}: {frac_joint:.1%} rows inside joint Dec [{JOINT_DEC_MIN}, {JOINT_DEC_MAX}]")

    if dec.min() > JOINT_DEC_MIN + 1.0:

        print(

            f"[WARN] {label}: dec_min={dec.min():.2f} — no southern-stripe coverage "

            f"below joint clip ({JOINT_DEC_MIN}°)?"

        )



    lnl_frac = float(pd.to_numeric(df["lnLExp_r"], errors="coerce").notna().mean())

    print(f"[OK] {label}: lnL coverage {lnl_frac:.1%}")

    if lnl_frac < 0.95:

        print(f"[WARN] {label}: lnL coverage < 95%; run patch_sdss_profile_winner.py --footprint full")



    base = prepare_null_strict_color_base(

        df,

        mag_column="modelMag_r",

        q0=Q0,

        q_column="expAB_r",

        is_legacy=False,

        sdss_ur_max=SDSS_UR_MAX_CDF,

        sdss_exp_winner_only=True,

    )

    strict_mag20 = slice_null_base_by_mag(

        base, mag_column="modelMag_r", mag_limit=SDSS_MAG20_LIMIT

    )

    strict_mag21 = slice_null_base_by_mag(base, mag_column="modelMag_r", mag_limit=MAG_LIMIT)



    print(

        f"[OK] {label}: strict+color+morph @ m<{SDSS_MAG20_LIMIT:g} "

        f"n={len(strict_mag20)} (min {SDSS_MIN_STRICT_MAG20_POOL})"

    )

    print(f"[OK] {label}: strict+color+morph @ m<{MAG_LIMIT:g} n={len(strict_mag21)}")



    if len(strict_mag20) < SDSS_MIN_STRICT_MAG20_POOL:

        print(f"[FAIL] {label}: strict mag20 pool too small")

        return False



    return True





def _check_catalog(

    path: Path,

    label: str,

    mag_col: str,

    is_legacy: bool,

    min_pool: int,

    *,

    check_joint_footprint: bool = True,

) -> bool:

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



    if check_joint_footprint:

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



    if is_legacy:

        base = prepare_null_strict_color_base(

            df,

            mag_column="tractor_mag_r",

            q0=Q0,

            q_column="expAB_r",

            exclude_legacy_types=LEGACY_CDF_TYPE_EXCLUDE,

            is_legacy=True,

            legacy_gr_max=LEGACY_GR_MAX_CDF,

            legacy_spiral_morph_only=True,

        )

        strict = slice_null_base_by_mag(base, mag_column="tractor_mag_r", mag_limit=MAG_LIMIT)

        min_required = LEGACY_MIN_STRICT_POOL

    else:

        if "lnLExp_r" not in df.columns:

            print(f"[FAIL] {label}: missing lnLExp_r (run patch_sdss_profile_winner.py)")

            return False

        base = prepare_null_strict_color_base(

            df,

            mag_column="modelMag_r",

            q0=Q0,

            q_column="expAB_r",

            is_legacy=False,

            sdss_ur_max=SDSS_UR_MAX_CDF,

            sdss_exp_winner_only=True,

        )

        strict = slice_null_base_by_mag(base, mag_column="modelMag_r", mag_limit=MAG_LIMIT)

        min_required = min_pool



    print(f"[OK] {label}: strict+color+morph @ m<{MAG_LIMIT} n={len(strict)} (min {min_required})")

    if len(strict) < min_required:

        print(f"[FAIL] {label}: strict pool too small")

        return False



    return True





def main() -> int:

    parser = argparse.ArgumentParser()

    parser.add_argument("--legacy-csv", default=str(DEFAULT_LEGACY))

    parser.add_argument("--sdss-csv", default=str(DEFAULT_SDSS))

    parser.add_argument("--sdss-v2-csv", default=str(DEFAULT_SDSS_V2))

    parser.add_argument(

        "--survey-version",

        choices=("v1", "v2", "both"),

        default="v1",

        help="v1: Legacy+v1 SDSS; v2: SDSS v2 only; both: all three.",

    )

    parser.add_argument("--min-pool", type=int, default=10_000)

    args = parser.parse_args()



    ok = True

    if args.survey_version in ("v1", "both"):

        ok_legacy = _check_catalog(

            Path(args.legacy_csv),

            "Legacy",

            mag_col="tractor_mag_r",

            is_legacy=True,

            min_pool=args.min_pool,

        )

        ok_sdss = _check_catalog(

            Path(args.sdss_csv),

            "SDSS v1",

            mag_col="rmag",

            is_legacy=False,

            min_pool=args.min_pool,

        )

        ok = ok and ok_legacy and ok_sdss



    if args.survey_version in ("v2", "both"):

        ok_v2 = _check_sdss_v2(Path(args.sdss_v2_csv), "SDSS v2")

        ok = ok and ok_v2



    if ok:

        print("[PASS] All null catalog sanity checks passed.")

        return 0

    print("[FAIL] One or more checks failed.")

    return 1





if __name__ == "__main__":

    sys.exit(main())


