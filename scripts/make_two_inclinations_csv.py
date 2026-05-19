from pathlib import Path

import pandas as pd
from astropy.io import fits

src = pd.read_csv("legacy_vs_galfit_inclination_comparison.csv")

# Include size diagnostics if already merged into the main comparison table.
has_reff = {"re_psf_arcsec", "shape_r_ls_arcsec"}.issubset(src.columns)
if has_reff and "delta_reff_arcsec" not in src.columns:
    src["delta_reff_arcsec"] = src["shape_r_ls_arcsec"] - src["re_psf_arcsec"]

base_cols = ["FRB", "galfit_inc_psf_deg", "ls_inc_deg", "ls_inc_err_deg", "type_ls", "sersic_ls_fit"]
if has_reff:
    base_cols += ["re_psf_arcsec", "shape_r_ls_arcsec", "delta_reff_arcsec"]

simp = src.loc[
    src["ls_inc_deg"].notna(),
    base_cols,
].copy()
simp["sersic_n_fit"] = pd.NA
is_ser = simp["type_ls"].astype(str).str.upper().eq("SER")
simp.loc[is_ser, "sersic_n_fit"] = simp.loc[is_ser, "sersic_ls_fit"]
simp = simp.drop(columns=["sersic_ls_fit"])
simp = simp.sort_values("FRB")
out_file = "legacy_vs_galfit_two_inclinations.csv"
simp.to_csv(out_file, index=False)

master = pd.read_csv("master_frb_summary.csv")
frbs = set(master["FRB"].astype(str))
survey_map = {}
for p in Path("large_cutouts").glob("*_flux.fits"):
    frb = p.name.split("_")[0]
    if frb not in frbs:
        continue
    try:
        with fits.open(p) as h:
            surv = str(h[0].header.get("SURVEY", "")).strip()
    except Exception:
        surv = ""
    survey_map[frb] = surv

panstarrs = sorted([f for f, s in survey_map.items() if "Pan-STARRS" in s])
legacy = sorted([f for f, s in survey_map.items() if "Legacy Survey" in s])

print("Saved:", out_file)
print("Rows in simplified file:", len(simp))
print("Pan-STARRS hosts in cutout headers:", len(panstarrs))
print("Pan-STARRS list:", ", ".join(panstarrs) if panstarrs else "(none)")
print("Legacy hosts in cutout headers:", len(legacy))
print("Source columns used for query: ['FRB','RA_deg','DEC_deg','inc_psf']")
