#!/usr/bin/env python3
"""Build root ``repeater_localizations.csv`` — all localized repeaters in this work.

Merges:
  * ``CHIME/repeater_localizations.csv`` (16 CHIME-discovered hosts)
  * ``master_frb_localization.csv`` ``repeater=yes`` rows not already in CHIME
    (20121102A, 20180301A, 20230814B)
  * ``20190520B`` (FAST; documented in SOURCES_AUDIT / Gordon et al. 2023)

Enrichment from ``CHIME/chime_host_magnitudes.csv``, ``CHIME/chime_hosts_inclination.csv``,
and Gordon ``Lz_host_data.csv`` / literature for non-CHIME hosts.

Run from repo root:
    python scripts/build_root_repeater_localizations.py
"""
from __future__ import annotations

import csv
from pathlib import Path

import pandas as pd
from astropy.coordinates import SkyCoord

REPO = Path(__file__).resolve().parents[1]
CHIME = REPO / "CHIME"
OUT = REPO / "repeater_localizations.csv"

# Base localization columns (shared with CHIME / master schemas)
BASE_COLS = [
    "frb",
    "ra_deg",
    "dec_deg",
    "ra_hms",
    "dec_dms",
    "coord_semantics",
    "xmin",
    "xmax",
    "ymin",
    "ymax",
    "survey",
    "ra_err_as",
    "dec_err_as",
    "major_sigma_as",
    "minor_sigma_as",
    "pa_deg",
    "z",
    "DM",
    "DM_MW",
    "DM_exgal",
    "status",
    "repeater",
    "repeater_source",
]

EXTRA_COLS = [
    "discovery_facility",
    "localization_facility",
    "catalog_origin",
    "mag_r",
    "mag_r_err",
    "mag_band",
    "mag_source",
    "inc_deg",
    "inc_err_deg",
    "b_a",
    "notes",
]

ALL_COLS = BASE_COLS + EXTRA_COLS

# CHIME-discovered: discovery always CHIME; localization facility by source
CHIME_LOCALIZATION = {
    "20180814A": "CHIME baseband",
    "20180916B": "EVN (Marcote et al. 2020)",
    "20181030A": "CHIME / optical host",
    "20190110C": "CHIME (Ibik et al. 2024)",
    "20190208A": "EVN PRECISE (Hewitt et al. 2024)",
    "20190303A": "CHIME Outriggers (via 20231204A)",
    "20190417A": "EVN PRECISE (Kirsten et al. 2025)",
    "20190711A": "ASKAP/CRAFT (Heintz et al. 2020)",
    "20191106C": "CHIME Outriggers (via 20231128A)",
    "20200120E": "EVN (Kirsten et al. 2022)",
    "20200223B": "CHIME (Ibik et al. 2024)",
    "20201124A": "ASKAP + EVN VLBI (Fong/Nimmo)",
    "20220912A": "DSA-110 (Ravi et al. 2023)",
    "20240114A": "MeerKAT + EVN (Snelders et al. 2025)",
    "20240209A": "CHIME Outriggers (Shah et al. 2024)",
    "20251229A": "CHIME Outriggers (ATel #17709)",
}

# Non-CHIME hosts already in master_frb_localization.csv
MASTER_EXTRA_META = {
    "20121102A": {
        "discovery_facility": "Arecibo",
        "localization_facility": "VLA/EVN VLBI",
        "catalog_origin": "master_frb_localization.csv",
        "mag_r": 23.73,
        "mag_r_err": "",
        "mag_band": "GMOS-N r",
        "mag_source": "Gordon et al. 2023 / Tendulkar et al. 2017",
        "inc_deg": "",
        "inc_err_deg": "",
        "b_a": "",
        "notes": "First localized repeater; PRS; excluded from CHIME-only catalog",
        "status": "secure",
    },
    "20180301A": {
        "discovery_facility": "Parkes",
        "localization_facility": "RealFAST/VLA",
        "catalog_origin": "master_frb_localization.csv",
        "mag_r": 21.21,
        "mag_r_err": "",
        "mag_band": "NOT r",
        "mag_source": "Gordon et al. 2023 / Bhandari et al. 2022",
        "inc_deg": "",
        "inc_err_deg": "",
        "b_a": "",
        "notes": "Parkes discovery; excluded from CHIME-only catalog",
        "status": "secure",
    },
    "20230814B": {
        "discovery_facility": "DSA-110",
        "localization_facility": "DSA-110",
        "catalog_origin": "master_frb_localization.csv",
        "mag_r": "",
        "mag_r_err": "",
        "mag_band": "",
        "mag_source": "",
        "inc_deg": "",
        "inc_err_deg": "",
        "b_a": "",
        "notes": "DSA designation FRB 20230814A (ATel #16191); pipeline fit_log missing",
        "status": "secure",
    },
}

# FAST repeater recognized in SOURCES_AUDIT / Gordon but not in master CSV
FAST_20190520B = {
    "frb": "20190520B",
    "ra_deg": 240.51779166666665,
    "dec_deg": -11.28813888888889,
    "ra_hms": "16:02:04.27",
    "dec_dms": "-11:17:17.3",
    "coord_semantics": "host",
    "xmin": "",
    "xmax": "",
    "ymin": "",
    "ymax": "",
    "survey": "FAST/VLA",
    "ra_err_as": "",
    "dec_err_as": "",
    "major_sigma_as": "",
    "minor_sigma_as": "",
    "pa_deg": "",
    "z": 0.2418,
    "DM": 1204.7,
    "DM_MW": 113.0,
    "DM_exgal": 1091.7,
    "status": "secure",
    "repeater": "yes",
    "repeater_source": (
        "Niu et al. 2022, Nature 606, 873 (doi:10.1038/s41586-022-04755-5); "
        "Gordon et al. 2023 host photometry"
    ),
    "discovery_facility": "FAST",
    "localization_facility": "VLA (Niu et al. 2022)",
    "catalog_origin": "literature (Gordon/SOURCES_AUDIT; not in master CSV)",
    "mag_r": 22.16,
    "mag_r_err": "",
    "mag_band": "SOAR r",
    "mag_source": "Gordon et al. 2023 / Niu et al. 2022",
    "inc_deg": "",
    "inc_err_deg": "",
    "b_a": "",
    "notes": (
        "PRS host; DM_MW≈113 (Niu+22 / Simha+23); excluded from CHIME-only catalog"
    ),
}


def _blank_base() -> dict:
    return {c: "" for c in ALL_COLS}


def _fmt_hms_dms(ra: float, dec: float) -> tuple[str, str]:
    c = SkyCoord(ra=ra, dec=dec, unit="deg")
    # Prefer colon HMS/DMS matching existing CSVs
    ra_hms = c.ra.to_string(unit="hour", sep=":", precision=2, pad=True)
    dec_dms = c.dec.to_string(unit="deg", sep=":", precision=1, pad=True, alwayssign=True)
    return ra_hms, dec_dms


def main() -> None:
    chime = pd.read_csv(CHIME / "repeater_localizations.csv", dtype=str).fillna("")
    master = pd.read_csv(REPO / "master_frb_localization.csv", dtype=str).fillna("")
    mags = pd.read_csv(CHIME / "chime_host_magnitudes.csv", dtype=str).fillna("")
    incs = pd.read_csv(CHIME / "chime_hosts_inclination.csv", dtype=str).fillna("")

    mag_by = {r["frb"]: r for _, r in mags.iterrows()}
    inc_by = {r["frb"]: r for _, r in incs.iterrows()}

    rows: list[dict] = []
    seen: set[str] = set()

    # --- 16 CHIME hosts ---
    for _, r in chime.iterrows():
        frb = str(r["frb"])
        out = _blank_base()
        for c in BASE_COLS:
            out[c] = r.get(c, "")
        out["discovery_facility"] = "CHIME"
        out["localization_facility"] = CHIME_LOCALIZATION.get(frb, "CHIME")
        out["catalog_origin"] = "CHIME/repeater_localizations.csv"
        if frb in mag_by:
            out["mag_r"] = mag_by[frb].get("mag_r", "")
            out["mag_r_err"] = mag_by[frb].get("mag_r_err", "")
            out["mag_band"] = mag_by[frb].get("mag_band", "")
            out["mag_source"] = mag_by[frb].get("source", "")
            note = mag_by[frb].get("notes", "")
            if note:
                out["notes"] = note
        if frb in inc_by:
            out["inc_deg"] = inc_by[frb].get("inc", "")
            out["inc_err_deg"] = inc_by[frb].get("inc_err", "")
            out["b_a"] = inc_by[frb].get("b_a", "")
        if frb == "20190208A":
            out["notes"] = (
                "True host r≈27.3 absent from Legacy cutout; no usable GALFIT inclination"
            )
        if frb == "20251229A":
            out["notes"] = "PRELIMINARY — ATel-only (no arXiv/refereed paper yet)"
        rows.append(out)
        seen.add(frb)

    # --- master repeater=yes not already in CHIME ---
    master_rep = master.loc[master["repeater"].str.lower() == "yes"]
    for _, r in master_rep.iterrows():
        frb = str(r["frb"])
        if frb in seen:
            continue
        meta = MASTER_EXTRA_META.get(frb, {})
        out = _blank_base()
        for c in BASE_COLS:
            out[c] = r.get(c, "")
        if not out["status"] or out["status"] == "pending":
            out["status"] = meta.get("status", "secure")
        for k, v in meta.items():
            if k == "status":
                continue
            out[k] = v
        # Prefer Gordon HMS if blank
        if not out["ra_hms"] or not out["dec_dms"]:
            try:
                hms, dms = _fmt_hms_dms(float(out["ra_deg"]), float(out["dec_deg"]))
                out["ra_hms"], out["dec_dms"] = hms, dms
            except (TypeError, ValueError):
                pass
        rows.append(out)
        seen.add(frb)

    # --- FAST 20190520B ---
    if "20190520B" not in seen:
        out = _blank_base()
        out.update(FAST_20190520B)
        rows.append(out)
        seen.add("20190520B")

    rows.sort(key=lambda d: d["frb"])

    with OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=ALL_COLS, extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow({c: row.get(c, "") for c in ALL_COLS})

    n_chime = sum(1 for r in rows if r["discovery_facility"] == "CHIME")
    n_other = len(rows) - n_chime
    print(f"Wrote {len(rows)} rows to {OUT}")
    print(f"  CHIME-discovered: {n_chime}")
    print(f"  Other facilities: {n_other} ({', '.join(sorted(seen - set(chime['frb'])))})")


if __name__ == "__main__":
    main()
