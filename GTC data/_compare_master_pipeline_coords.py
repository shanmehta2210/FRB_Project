"""Compare master FRB coords vs pipeline host coords for the GTC-17 list."""
from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path(__file__).resolve().parents[1]
FRBS = [
    "20220105A", "20221219A", "20240119A", "20190523A", "20250518", "20240304B",
    "20220307B", "20220912A", "20210117A", "20220204A", "20220506D", "20230501A",
    "20230521B", "20230814B", "20230913", "20240203", "20221116A",
]

master = pd.read_csv(ROOT / "master_frb_localization.csv").set_index("frb")
galfit = pd.read_csv(ROOT / "pipeline_galfit_results.csv").set_index("frb")


def sep_as(ra1, dec1, ra2, dec2) -> float:
    c1 = SkyCoord(ra=ra1 * u.deg, dec=dec1 * u.deg)
    c2 = SkyCoord(ra=ra2 * u.deg, dec=dec2 * u.deg)
    return c1.separation(c2).arcsec


def pipeline_input_coords(frb: str) -> tuple[float | None, float | None]:
    log = ROOT / "pipeline_scripts" / "Output" / f"{frb}_all" / "master_run.log"
    if not log.is_file():
        return None, None
    text = log.read_text(encoding="utf-8", errors="replace")
    m = re.search(r"--ra'\s+'([0-9.+-]+)'\s+'--dec'\s+'([0-9.+-]+)'", text)
    if m:
        return float(m.group(1)), float(m.group(2))
    m = re.search(r"localization host \(--ra/--dec\) -> RA=([0-9.+-]+), Dec=([0-9.+-]+)", text)
    if m:
        return float(m.group(1)), float(m.group(2))
    return None, None


def galfit_host_coords(frb: str, host_num: int) -> tuple[float | None, float | None]:
    hc = ROOT / "pipeline_scripts" / "Output" / f"{frb}_all" / "host_components.csv"
    if not hc.is_file():
        return None, None
    h = pd.read_csv(hc)
    row = h.loc[h["NUMBER"] == host_num]
    if row.empty:
        return None, None
    r = row.iloc[0]
    return float(r["ALPHAWIN_J2000"]), float(r["DELTAWIN_J2000"])


def astropath_best(frb: str) -> tuple[float | None, float | None, float | None]:
    ap = ROOT / "pipeline_scripts" / "Output" / f"{frb}_all" / "astropath_posteriors.csv"
    if not ap.is_file():
        return None, None, None
    d = pd.read_csv(ap)
    if d.empty:
        return None, None, None
    best = d.loc[d["posterior_O"].idxmax()]
    return float(best["ra"]), float(best["dec"]), float(best["posterior_O"])


rows = []
for frb in FRBS:
    m = master.loc[frb]
    mra, mdec = float(m.ra_deg), float(m.dec_deg)
    out_dir = ROOT / "pipeline_scripts" / "Output" / f"{frb}_all"
    has_out = out_dir.is_dir()
    pin_ra, pin_dec = pipeline_input_coords(frb) if has_out else (None, None)

    gf_sep = ap_sep = pin_sep = None
    gf_ra = gf_dec = None
    ap_ra = ap_dec = ap_p = None

    if frb in galfit.index:
        hn = int(galfit.loc[frb, "host_number"])
        gf_ra, gf_dec = galfit_host_coords(frb, hn)
        if gf_ra is not None:
            gf_sep = sep_as(mra, mdec, gf_ra, gf_dec)

    if has_out:
        ap_ra, ap_dec, ap_p = astropath_best(frb)
        if ap_ra is not None:
            ap_sep = sep_as(mra, mdec, ap_ra, ap_dec)
        if pin_ra is not None:
            pin_sep = sep_as(mra, mdec, pin_ra, pin_dec)

    rows.append(
        {
            "frb": frb,
            "coord_semantics": m.coord_semantics,
            "survey": m.survey,
            "master_ra": mra,
            "master_dec": mdec,
            "pipeline_output": has_out,
            "galfit_in_results": frb in galfit.index,
            "pipeline_input_sep_as": pin_sep,
            "galfit_host_sep_as": gf_sep,
            "astropath_best_sep_as": ap_sep,
            "astropath_P_O": ap_p,
            "galfit_host_ra": gf_ra,
            "galfit_host_dec": gf_dec,
            "astropath_best_ra": ap_ra,
            "astropath_best_dec": ap_dec,
        }
    )

df = pd.DataFrame(rows)
out = ROOT / "GTC data" / "coord_compare_17.csv"
df.to_csv(out, index=False, float_format="%.8f")

print("frb,pipeline,galfit,d_master_vs_galfit_host,d_master_vs_astropath_best,d_master_vs_pipeline_input")
for _, r in df.iterrows():
    def fmt(v):
        return "NA" if v is None or (isinstance(v, float) and pd.isna(v)) else f"{v:.3f}"

    print(
        f"{r.frb},{r.pipeline_output},{r.galfit_in_results},"
        f"{fmt(r.galfit_host_sep_as)},{fmt(r.astropath_best_sep_as)},{fmt(r.pipeline_input_sep_as)}"
    )
