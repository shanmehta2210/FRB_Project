"""Seed table for the 56-host property sweep.

Collects what the repository already holds per host (localization, DM, z,
repeater flag, literature cite, pipeline photometry, GALFIT / AstroPhot /
isophote axis ratios, imaging survey) into one CSV that the literature
sweep agents extend. Output:
    Publication/materials/host_properties/seed_56.csv
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Publication" / "materials" / "host_properties"
Q0 = 0.2


def hubble_i(q: float) -> float:
    if q is None or not (Q0 < q <= 1):
        return float("nan")
    return math.degrees(math.acos(math.sqrt((q * q - Q0 * Q0) / (1 - Q0 * Q0))))


def load_json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8")) if p.exists() else {}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    win = pd.read_csv(ROOT / "plots/plots_null/v2/frb_vs_hsc_confirmed50/confirmed50_q.csv", dtype={"frb": str})
    win["i_deg"] = [hubble_i(float(q)) for q in win.q]
    loc = pd.read_csv(ROOT / "master_frb_localization.csv", dtype={"frb": str}).set_index("frb")
    prod = pd.read_csv(ROOT / "pipeline_galfit_results.csv", dtype={"frb": str}).set_index("frb")
    audit = pd.read_csv(ROOT / "production67_lit_astropath_audit.csv", dtype={"frb": str}).set_index("frb")
    n23 = pd.read_csv(ROOT / "Archive/original paper/frb_sample.txt", sep="\t", dtype={"FRB": str}).set_index("FRB")

    rows = []
    for _, w in win.iterrows():
        f = w.frb
        L = loc.loc[f]
        P = prod.loc[f]
        hdr = fits.getheader(ROOT / f"large_cutouts/{f}_flux.fits")
        survey = str(hdr.get("SURVEY", "")).strip()
        version = str(hdr.get("VERSION", "")).strip()
        if "HSC" in survey:
            imaging = "HSC-SSP PDR2"
        elif "Pan-STARRS" in survey or "PSLIB_V" in hdr:
            imaging = "PS1 DR1"
        else:
            imaging = f"LS {version}" if version else "LS"
        ver = ROOT / "pipeline_scripts/verification/outputs/per_host" / f
        ap = load_json(ver / "astrophot.json")
        iso = load_json(ver / "isophote.json")
        ap_q = ap.get("ap_q")
        iso_q = iso.get("iso_q_at_1re_data")
        r = {
            "frb": f,
            "ra_deg": L.ra_deg,
            "dec_deg": L.dec_deg,
            "ra_hms": L.ra_hms,
            "dec_dms": L.dec_dms,
            "loc_survey": L.survey,
            "z": L.z,
            "DM_obs": L.DM,
            "DM_MW_NE2001": L.DM_MW,
            "DM_exgal": L.DM_exgal,
            "repeater": L.repeater,
            "repeater_source": L.repeater_source,
            "host_citation": audit.loc[f, "citation"] if f in audit.index else "",
            "in_bhardwaj2024": f in n23.index,
            "imaging": imaging,
            "m_r_galfit": P.mag,
            "m_r_galfit_err": P.mag_err,
            "ref_survey": P.ref_survey,
            "m_r_ref": P.ref_mag,
            "re_arcsec_galfit": P.re,
            "n_galfit": P.n,
            "q_galfit_production": P.b_a,
            "q_winning": w.q,
            "q_winning_err": w.q_err,
            "q_winning_leg": w.leg,
            "i_galfit_deg": w.i_deg,
            "q_astrophot": ap_q,
            "i_astrophot_deg": hubble_i(ap_q) if ap_q else float("nan"),
            "q_isophote_1re": iso_q,
            "i_isophote_deg": hubble_i(iso_q) if iso_q else float("nan"),
            "i_bhardwaj2024": n23.loc[f, "Inclination angle (deg)"] if f in n23.index else "",
            "SFR_bhardwaj2024": n23.loc[f, "SFR"] if f in n23.index else "",
            "logsSFR_bhardwaj2024": n23.loc[f, "log(sSFR)"] if f in n23.index else "",
            "Mstar_1e10_bhardwaj2024": n23.loc[f, "Mass"] if f in n23.index else "",
            "DM_host_bhardwaj2024": n23.loc[f, "DM_Host"] if f in n23.index else "",
            "m_r_bhardwaj2024": n23.loc[f, "Host m_r"] if f in n23.index else "",
        }
        rows.append(r)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "seed_56.csv", index=False)
    print(f"wrote {OUT / 'seed_56.csv'} ({len(df)} rows)")


if __name__ == "__main__":
    main()
