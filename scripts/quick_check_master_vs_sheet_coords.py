import numpy as np
import pandas as pd

master = pd.read_csv("master_frb_summary.csv")
sheet = pd.read_excel("FRB_Inclination angle estimate sheet.xlsx")

# Flexible column detection for the sheet
frb_col = next((c for c in sheet.columns if str(c).strip().lower() in {"frb", "frb name", "frb_name", "name"}), None)
ra_col = next((c for c in sheet.columns if "ra" in str(c).strip().lower() and "deg" in str(c).strip().lower()), None)
dec_col = next((c for c in sheet.columns if ("dec" in str(c).strip().lower() or str(c).strip().lower().startswith("de")) and "deg" in str(c).strip().lower()), None)
if ra_col is None:
    ra_col = next((c for c in sheet.columns if "ra" in str(c).strip().lower()), None)
if dec_col is None:
    dec_col = next((c for c in sheet.columns if "dec" in str(c).strip().lower() or str(c).strip().lower().startswith("de")), None)

if frb_col is None or ra_col is None or dec_col is None:
    raise RuntimeError(f"Could not detect sheet columns. Found: {list(sheet.columns)}")

m = master[["FRB", "RA_deg", "DEC_deg"]].copy()
s = sheet[[frb_col, ra_col, dec_col]].copy()
s.columns = ["FRB", "RA_sheet", "DEC_sheet"]
s["FRB"] = s["FRB"].astype(str).str.strip()

j = m.merge(s, on="FRB", how="inner")

# Angular separation approximation in arcsec
j["dRA_arcsec"] = (j["RA_deg"] - j["RA_sheet"]) * 3600.0 * np.cos(np.radians(j["DEC_deg"]))
j["dDEC_arcsec"] = (j["DEC_deg"] - j["DEC_sheet"]) * 3600.0
j["sep_arcsec"] = np.hypot(j["dRA_arcsec"], j["dDEC_arcsec"])

print("Detected sheet columns:", frb_col, ra_col, dec_col)
print("Matched FRBs:", len(j), "/", len(m))
print("Median sep arcsec:", round(float(j["sep_arcsec"].median()), 4))
print("Mean sep arcsec:", round(float(j["sep_arcsec"].mean()), 4))
print("Max sep arcsec:", round(float(j["sep_arcsec"].max()), 4))
print("FRBs with sep <= 1 arcsec:", int((j["sep_arcsec"] <= 1).sum()))
print("FRBs with sep <= 5 arcsec:", int((j["sep_arcsec"] <= 5).sum()))

print("\nTop separations:")
print(
    j[["FRB", "RA_deg", "RA_sheet", "DEC_deg", "DEC_sheet", "sep_arcsec"]]
    .sort_values("sep_arcsec", ascending=False)
    .head(10)
    .to_string(index=False)
)
