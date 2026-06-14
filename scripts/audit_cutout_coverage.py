#!/usr/bin/env python3

"""Probe Legacy / PS1 / DES r-band coverage for FRBs missing large cutouts."""



from __future__ import annotations



import argparse

from pathlib import Path



import pandas as pd



from cutout_fetch_common import (

    DEC_SPLIT,

    LS_LAYERS_GLOBAL,

    LS_LAYERS_NORTH,

    PS1_DEC_MIN,

    legacy_layers_for_dec,

    probe_des,

    probe_legacy,

    get_ps1_r_filename,

)



REPO_ROOT = Path(__file__).resolve().parents[1]

LOC_CSV = REPO_ROOT / "master_frb_localization.csv"

CUTOUT_DIR = REPO_ROOT / "large_cutouts"

DEFAULT_OUT = CUTOUT_DIR / "coverage_audit.csv"





def missing_frbs(loc: pd.DataFrame, only_missing: bool) -> pd.DataFrame:

    have = {p.stem.replace("_flux", "") for p in CUTOUT_DIR.glob("*_flux.fits")}

    if only_missing:

        return loc[~loc["frb"].isin(have)].copy()

    return loc.copy()





def recommend_tier(row: dict) -> str:

    if row["legacy_global_ok"]:

        return "legacy_global"

    if row["legacy_hemisphere_ok"]:

        return "legacy_hemisphere"

    if row["ps1_ok"]:

        return "ps1"

    if row["des_ok"]:

        return "des"

    return "none"





def audit_row(frb: str, ra: float, dec: float, survey: str) -> dict:

    leg_g = probe_legacy(ra, dec, LS_LAYERS_GLOBAL[0])

    hem_layers = legacy_layers_for_dec(dec)

    leg_h = probe_legacy(ra, dec, hem_layers[0]) if not leg_g else True

    ps1_ok = dec > PS1_DEC_MIN and get_ps1_r_filename(ra, dec) is not None

    des_ok = probe_des(ra, dec)



    row = {

        "frb": frb,

        "ra_deg": ra,

        "dec_deg": dec,

        "survey_csv": survey,

        "legacy_global_ok": leg_g,

        "legacy_hemisphere_ok": leg_h,

        "hemisphere_layers": hem_layers[0] if dec < DEC_SPLIT else LS_LAYERS_NORTH[0],

        "ps1_ok": ps1_ok,

        "des_ok": des_ok,

        "notes": "",

    }

    row["recommended_tier"] = recommend_tier(row)

    if dec <= PS1_DEC_MIN:

        row["notes"] = "dec<=PS1 limit; use Legacy or DES"

    return row





def main():

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--csv", type=Path, default=LOC_CSV)

    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)

    parser.add_argument("--frb", nargs="*", help="Restrict to these FRB names")

    parser.add_argument("--all", action="store_true", help="Audit every CSV row, not only missing cutouts")

    args = parser.parse_args()



    loc = pd.read_csv(args.csv)

    targets = missing_frbs(loc, only_missing=not args.all)

    if args.frb:

        targets = targets[targets["frb"].isin(args.frb)]



    CUTOUT_DIR.mkdir(parents=True, exist_ok=True)

    rows = []

    n = len(targets)

    print(f"[audit] {n} FRB(s) — probes: Legacy, PS1, DES (no SkyMapper)")

    for i, (_, r) in enumerate(targets.iterrows(), 1):

        frb = r["frb"]

        ra, dec = float(r["ra_deg"]), float(r["dec_deg"])

        survey = str(r.get("survey", ""))

        print(f"  [{i}/{n}] {frb} dec={dec:.2f} ...", flush=True)

        rows.append(audit_row(frb, ra, dec, survey))



    out_df = pd.DataFrame(rows)

    out_df.to_csv(args.out, index=False)

    print(f"[audit] wrote {args.out}")

    if len(out_df):

        print(out_df["recommended_tier"].value_counts().to_string())





if __name__ == "__main__":

    main()

