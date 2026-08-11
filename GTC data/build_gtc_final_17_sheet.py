#!/usr/bin/env python3
"""Build GTC proposal sheet for the final 17 FRBs."""
import json
from pathlib import Path

import pandas as pd
from astropy import units as u
from astropy.coordinates import Angle

REPO = Path(__file__).resolve().parents[1]
OUT_DIR = REPO / "GTC data"
GALFIT = REPO / "pipeline_galfit_results.csv"

# User-selected 17 FRBs (Package 1: 8 production/science; Package 2: 9 trial hosts)
FRB_ORDER = [
    # Package 1 — visible from user's pipeline list + science-flagged (excl. 20220319D)
    "20220105A",
    "20221219A",
    "20240119A",
    "20190523A",
    "20250518",
    "20240304B",
    "20220307B",
    "20220912A",
    # Package 2 — 6 host-missing + 3 photometry-problematic (excl. 20230930A)
    "20210117A",
    "20220204A",
    "20220506D",
    "20230501A",
    "20230521B",
    "20230814B",
    "20230913",
    "20240203",
    "20221116A",
]

# Redshift overrides when master_frb_localization.csv is blank but literature table has z
Z_OVERRIDE = {
    "20220912A": 0.0771,  # new_confident_hosts.txt / Ravi et al. 2023
}

TONIGHT = "2026-06-24"


def fmt_coords(ra_hms, dec_dms) -> str:
    dec = str(dec_dms).strip()
    if dec and not dec.startswith(("+", "-")):
        dec = "+" + dec
    return f"{ra_hms} {dec}"


def fmt_coords_deg(ra_deg: float, dec_deg: float) -> str:
    ra = Angle(ra_deg * u.deg).to_string(unit=u.hour, sep=" ", precision=2, pad=True)
    dec = Angle(dec_deg * u.deg).to_string(sep=" ", precision=1, alwayssign=True, pad=True)
    return f"{ra} {dec}"


def fmt_z(v) -> str:
    if pd.isna(v) or str(v).strip() == "":
        return ""
    return f"{float(v):.4f}".rstrip("0").rstrip(".")


def galfit_output_dir(frb: str) -> Path:
    return REPO / "pipeline_scripts" / "Output" / f"{frb}_all"


def resolve_host_number(frb: str, galfit: pd.DataFrame) -> int | None:
    """Host component ID when GALFIT ran; None if fitting never completed."""
    if frb in galfit.index and pd.notna(galfit.loc[frb, "host_number"]):
        return int(galfit.loc[frb, "host_number"])

    out_dir = galfit_output_dir(frb)
    if not (out_dir / "fit.log").is_file():
        return None

    meta_path = out_dir / "cutout_meta.json"
    if meta_path.is_file():
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        if meta.get("host_number") is not None:
            return int(meta["host_number"])

    audit_path = out_dir / "sky_fit_audit.json"
    if audit_path.is_file():
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        if audit.get("host_number") is not None:
            return int(audit["host_number"])

    return None


def galfit_host_coords(frb: str, galfit: pd.DataFrame) -> tuple[float | None, float | None]:
    """SExtractor ALPHAWIN/DELTAWIN for the host component GALFIT fits."""
    host_num = resolve_host_number(frb, galfit)
    if host_num is None:
        return None, None

    hc_path = galfit_output_dir(frb) / "host_components.csv"
    if not hc_path.is_file():
        return None, None

    hc = pd.read_csv(hc_path)
    row = hc.loc[hc["NUMBER"] == host_num]
    if row.empty:
        return None, None
    r = row.iloc[0]
    return float(r["ALPHAWIN_J2000"]), float(r["DELTAWIN_J2000"])


def tonight_label(vis_row) -> str:
    if vis_row is None or vis_row.empty:
        return "Unknown"
    r = vis_row.iloc[0]
    if bool(r.get("rigorous_science_pass")):
        start = r.get("best_start_utc", "")
        end = r.get("best_end_utc", "")
        if pd.notna(start) and pd.notna(end):
            return f"Yes ({start} – {end} UTC)"
        return "Yes"
    reasons = []
    if not bool(r.get("gate1_mechanical_pass", True)):
        reasons.append("below GTC elevation limit")
    if not bool(r.get("gate2_airmass_pass", True)):
        reasons.append("airmass")
    if not bool(r.get("duration_pass", True)):
        reasons.append("no >=30 min window")
    if not reasons:
        reasons.append("rigorous gates failed")
    return "No (" + "; ".join(reasons) + ")"


def main() -> None:
    master = pd.read_csv(REPO / "master_frb_localization.csv").set_index("frb")
    galfit = pd.read_csv(GALFIT).set_index("frb")
    vis = pd.read_csv(REPO / f"GTC data/visibility/nightly/gtc_visibility_{TONIGHT}.csv")

    rows = []
    for frb in FRB_ORDER:
        if frb not in master.index:
            raise KeyError(f"{frb} missing from master_frb_localization.csv")
        m = master.loc[frb]
        v = vis[vis["frb"] == frb]

        z = Z_OVERRIDE.get(frb, m.get("z"))
        loc_coords = fmt_coords(m["ra_hms"], m["dec_dms"])
        loc_ra, loc_dec = float(m["ra_deg"]), float(m["dec_deg"])

        host_ra, host_dec = galfit_host_coords(frb, galfit)
        host_coords = fmt_coords_deg(host_ra, host_dec) if host_ra is not None else ""

        rows.append(
            {
                "FRB": frb,
                "Coordinates": loc_coords,
                "Host_galaxy_coordinates": host_coords,
                "Redshift": fmt_z(z),
                f"GTC_observable_{TONIGHT}": tonight_label(v),
                "RA_deg": loc_ra,
                "Dec_deg": loc_dec,
                "Host_RA_deg": host_ra if host_ra is not None else "",
                "Host_Dec_deg": host_dec if host_dec is not None else "",
            }
        )

    df = pd.DataFrame(rows)

    out_cols = [
        "FRB",
        "Coordinates",
        "Host_galaxy_coordinates",
        "Redshift",
        f"GTC_observable_{TONIGHT}",
    ]
    csv_path = OUT_DIR / "gtc_final_17_proposal.csv"
    xlsx_path = OUT_DIR / "gtc_final_17_proposal.xlsx"

    df[out_cols].to_csv(csv_path, index=False)

    portal_cols = [
        "FRB",
        "RA_deg",
        "Dec_deg",
        "Host_RA_deg",
        "Host_Dec_deg",
    ] + out_cols[1:]
    portal = df[portal_cols].copy()

    with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
        df[out_cols].to_excel(writer, sheet_name="GTC_17_proposal", index=False)
        portal.to_excel(writer, sheet_name="with_RA_Dec_deg", index=False)

    print(f"Wrote {csv_path}")
    print(f"Wrote {xlsx_path}")
    print(df[out_cols].to_string(index=False))


if __name__ == "__main__":
    main()
