"""
Integrates 12 new FRB GALFIT results into:
1. galfit_vs_legacy_master.csv  (with Legacy Survey TAP query via pyvo)
2. galfit_vs_legacy_quick_read.csv
3. master_frb_summary.csv  (galfit fields + new n_err column)

Uses the same pyvo TAP approach as compare_galfit_vs_tractor_inclination.py.
Inclination from b/a: Hubble formula with q0=0.2
  cos^2(i) = (q^2 - q0^2) / (1 - q0^2)
  if q <= q0: i = 90 deg
"""

import math, time
import numpy as np
import pandas as pd
import pyvo

# ─── Config ─────────────────────────────────────────────────────────────────────
import os
BASE = r"C:\Users\lenovo\Desktop\Bhardwajetal_2024_nature_inclination_angle-main"

EXCLUDE = {"20190611B", "20190711A", "20230526A", "20240310A"}
ALL_16 = [
    "20190611B","20190711A","20200430A","20220105A","20220725A","20221106A",
    "20230526A","20230708A","20230902A","20231226A","20240201A","20240208A",
    "20240210A","20240304A","20240310A","20240318A"
]
FRB12 = [f for f in ALL_16 if f not in EXCLUDE]
Q0 = 0.2

# ─── Hubble inclination helpers (identical to compare_galfit_vs_tractor) ────────
def incl_from_q(q, q0=Q0):
    if not np.isfinite(q): return np.nan
    if q <= q0: return 90.0
    val = (q*q - q0*q0) / (1.0 - q0*q0)
    val = min(1.0, max(0.0, val))
    return math.degrees(math.acos(math.sqrt(val)))

def q_from_e1e2(e1, e2):
    eabs = float(np.hypot(e1, e2))
    if eabs >= 1.0: return np.nan, eabs
    return (1.0 - eabs) / (1.0 + eabs), eabs

def sigma_from_ivar(ivar):
    if not np.isfinite(ivar) or ivar <= 0: return np.nan
    return float(1.0 / np.sqrt(ivar))

def incl_from_q_array(q, q0=Q0):
    q = np.asarray(q, dtype=float)
    val = (q*q - q0*q0) / (1.0 - q0*q0)
    val = np.clip(val, 0.0, 1.0)
    return np.degrees(np.arccos(np.sqrt(val)))

def incl_err_mc(e1, e2, e1_ivar, e2_ivar, n_draws=4000):
    s1, s2 = sigma_from_ivar(e1_ivar), sigma_from_ivar(e2_ivar)
    if not (np.isfinite(s1) and np.isfinite(s2)): return np.nan
    rng = np.random.default_rng(12345)
    e1_d = rng.normal(float(e1), s1, n_draws)
    e2_d = rng.normal(float(e2), s2, n_draws)
    eabs = np.hypot(e1_d, e2_d)
    valid = eabs < 1.0
    if np.count_nonzero(valid) < 50: return np.nan
    q = (1.0 - eabs[valid]) / (1.0 + eabs[valid])
    return float(np.nanstd(incl_from_q_array(q), ddof=1))

def ba_to_inc(ba, q0=Q0):
    return incl_from_q(float(ba), q0) if ba is not None else None

def ba_err_to_inc_err(ba, ba_err, q0=Q0):
    """Symmetric finite difference. Does NOT zero out error at edge-on boundary."""
    try:
        ba, ba_err = float(ba), float(ba_err)
        q_up   = min(1.0, ba + ba_err)
        q_down = max(0.0, ba - ba_err)
        return abs(incl_from_q(q_up, q0) - incl_from_q(q_down, q0)) / 2.0
    except: return 0.0

def _ra_clause(ra, dra):
    lo, hi = ra - dra, ra + dra
    if lo < 0:   return f"(ra > {lo+360:.8f} OR ra < {hi:.8f})"
    if hi > 360: return f"(ra > {lo:.8f} OR ra < {hi-360:.8f})"
    return f"ra > {lo:.8f} AND ra < {hi:.8f}"

def query_nearest(svc, ra, dec, radius_arcsec=10.0):
    dec_clip = max(-85.0, min(85.0, dec))
    dra  = (radius_arcsec / 3600.0) / math.cos(math.radians(dec_clip))
    ddec = radius_arcsec / 3600.0
    q = f"""
    SELECT TOP 300 objid, ra, dec, type, flux_r, sersic, shape_e1, shape_e2,
                   shape_e1_ivar, shape_e2_ivar, shape_r
    FROM ls_dr10.tractor
    WHERE {_ra_clause(ra, dra)}
      AND dec > {dec-ddec:.8f} AND dec < {dec+ddec:.8f}
      AND flux_r > 0
    """
    tab = svc.search(q).to_table()
    if len(tab) == 0: return None

    ra_arr  = np.array(tab["ra"],  dtype=float)
    dec_arr = np.array(tab["dec"], dtype=float)
    dra_as  = (ra_arr  - ra)  * np.cos(np.radians(dec)) * 3600.0
    ddec_as = (dec_arr - dec) * 3600.0
    sep = np.hypot(dra_as, ddec_as)

    types = np.array(tab["type"]).astype(str)
    non_psf = np.where(types != "PSF")[0]
    idx = int(non_psf[np.argmin(sep[non_psf])]) if len(non_psf) > 0 else int(np.argmin(sep))
    r = tab[idx]
    return {
        "objid":         int(r["objid"]),
        "ra_ls":         float(r["ra"]),
        "dec_ls":        float(r["dec"]),
        "type_ls":       str(r["type"]),
        "flux_r_ls":     float(r["flux_r"]),
        "sersic_ls_fit": float(r["sersic"]),
        "shape_e1":      float(r["shape_e1"]),
        "shape_e2":      float(r["shape_e2"]),
        "shape_e1_ivar": float(r["shape_e1_ivar"]),
        "shape_e2_ivar": float(r["shape_e2_ivar"]),
        "sep_arcsec":    float(sep[idx]),
        "n_candidates":  int(len(tab)),
        "shape_r_ls_arcsec": float(r["shape_r"]) * 0.262 if r["shape_r"] else None,
    }

# ─── Load data ──────────────────────────────────────────────────────────────────
results = pd.read_csv(os.path.join(BASE, "new_16_frbs_galfit_results.csv"))
master  = pd.read_csv(os.path.join(BASE, "galfit_vs_legacy_master.csv"))
quick   = pd.read_csv(os.path.join(BASE, "galfit_vs_legacy_quick_read.csv"))
summary = pd.read_csv(os.path.join(BASE, "master_frb_summary.csv"))

# Remove any previous entries for these 12 FRBs before appending fresh ones
master = master[~master["FRB"].isin(FRB12)].copy()
quick  = quick[~quick["FRB"].isin(FRB12)].copy()

results12 = results[results["FRB"].isin(FRB12)].copy()
print(f"Processing {len(results12)} FRBs: {results12['FRB'].tolist()}")

# ─── TAP queries ────────────────────────────────────────────────────────────────
svc = pyvo.dal.TAPService("https://datalab.noirlab.edu/tap")

ls_data = {}
for frb in FRB12:
    smask = summary["FRB"] == frb
    if not smask.any():
        print(f"  [SKIP] {frb}: not in summary")
        continue
    srow = summary[smask].iloc[0]
    ra, dec = float(srow["RA_deg"]), float(srow["DEC_deg"])
    print(f"  [{frb}] querying at host RA={ra:.5f}, DEC={dec:.5f} ...", end=" ", flush=True)

    match = None
    last_err = None
    for attempt in range(3):
        try:
            match = query_nearest(svc, ra, dec)
            break
        except Exception as e:
            last_err = str(e)
            time.sleep(2.0)

    if match is None:
        print(f"NO MATCH (err: {last_err})")
        ls_data[frb] = {"match_found": False, "query_error": last_err}
    else:
        q_ls, eabs = q_from_e1e2(match["shape_e1"], match["shape_e2"])
        match["ellipticity_abs"] = eabs
        match["q_ls_from_e"]     = q_ls
        match["shape_e1_sigma"]  = sigma_from_ivar(match["shape_e1_ivar"])
        match["shape_e2_sigma"]  = sigma_from_ivar(match["shape_e2_ivar"])
        match["ls_inc_deg"]      = incl_from_q(q_ls)
        match["ls_inc_err_deg"]  = incl_err_mc(match["shape_e1"], match["shape_e2"],
                                                match["shape_e1_ivar"], match["shape_e2_ivar"])
        match["sersic_n_ls_fit"] = match["sersic_ls_fit"]
        match["match_found"]     = True
        ls_data[frb] = match
        print(f"OK (sep={match['sep_arcsec']:.2f}\", q={q_ls:.3f}, inc={match['ls_inc_deg']:.1f}°)")

def safe(v, default=None):
    try: return float(v)
    except: return default

# ─── Build new rows ─────────────────────────────────────────────────────────────
new_master_rows, new_quick_rows = [], []

for _, row in results12.iterrows():
    frb = row["FRB"]
    ba_psf  = safe(row.get("b_a_psf"))
    bae_psf = safe(row.get("b_a_err_psf"), 0.0)
    n_psf   = safe(row.get("n_psf"))
    re_psf  = safe(row.get("re_psf"))
    chi2_psf = safe(row.get("chi2nu_psf"))

    smask = summary["FRB"] == frb
    ra  = float(summary[smask]["RA_deg"].iloc[0])  if smask.any() else None
    dec = float(summary[smask]["DEC_deg"].iloc[0]) if smask.any() else None

    galfit_inc = ba_to_inc(ba_psf)
    galfit_err = ba_err_to_inc_err(ba_psf, bae_psf) if ba_psf else None
    re_arcsec  = round(re_psf * 0.262, 4) if re_psf else None

    ls = ls_data.get(frb, {"match_found": False})

    new_master_rows.append({
        "FRB": frb, "RA_deg": ra, "DEC_deg": dec,
        "galfit_inc_psf_deg": round(galfit_inc, 4) if galfit_inc is not None else None,
        "match_found": ls.get("match_found", False),
        "objid": ls.get("objid"), "ra_ls": ls.get("ra_ls"), "dec_ls": ls.get("dec_ls"),
        "type_ls": ls.get("type_ls"), "flux_r_ls": ls.get("flux_r_ls"),
        "sersic_ls_fit": ls.get("sersic_ls_fit"),
        "shape_e1": ls.get("shape_e1"), "shape_e2": ls.get("shape_e2"),
        "shape_e1_ivar": ls.get("shape_e1_ivar"), "shape_e2_ivar": ls.get("shape_e2_ivar"),
        "sep_arcsec": ls.get("sep_arcsec"), "n_candidates": ls.get("n_candidates"),
        "ellipticity_abs": ls.get("ellipticity_abs"), "q_ls_from_e": ls.get("q_ls_from_e"),
        "shape_e1_sigma": ls.get("shape_e1_sigma"), "shape_e2_sigma": ls.get("shape_e2_sigma"),
        "ls_inc_deg": ls.get("ls_inc_deg"), "ls_inc_err_deg": ls.get("ls_inc_err_deg"),
        "delta_deg_ls_minus_galfit": round(float(ls["ls_inc_deg"]) - galfit_inc, 4)
            if ls.get("ls_inc_deg") and galfit_inc else None,
        "re_psf_arcsec": re_arcsec, "shape_r_ls_arcsec": ls.get("shape_r_ls_arcsec"),
        "sep_arcsec_reff": ls.get("sep_arcsec"), "tractor_objid_reff": ls.get("objid"),
        "delta_reff_arcsec": round(abs(re_arcsec - float(ls["shape_r_ls_arcsec"])), 4)
            if re_arcsec and ls.get("shape_r_ls_arcsec") else None,
        "tractor_objid": ls.get("objid"), "sersic_n_ls_fit": ls.get("sersic_n_ls_fit"),
        "query_error": ls.get("query_error"),
        "galfit_chi2nu_local": chi2_psf,
        "ls_chi2nu_local": None, "chi2_superior_model": None,
        "fit_region_found": None, "fit_region_pixels": None,
        "galfit_chi2nu_fitregion": None, "legacy_chi2nu_fitregion": None,
        "chi2nu_ratio_legacy_over_galfit_fitregion": None,
        "chi2nu_ratio_assessment_fitregion": None,
        "chi2nu_functionally_equivalent_fitregion": None,
        "sersic_n_fit": round(n_psf, 4) if n_psf else None,
    })

    new_quick_rows.append({
        "FRB": frb,
        "galfit_n":    round(n_psf, 2) if n_psf else None,
        "galfit_inc":  round(galfit_inc, 2) if galfit_inc else None,
        "galfit_err":  round(galfit_err, 2) if galfit_err else None,
        "legacy_type": ls.get("type_ls"),
        "legacy_n":    round(float(ls["sersic_n_ls_fit"]), 2) if ls.get("sersic_n_ls_fit") else None,
        "legacy_inc":  round(float(ls["ls_inc_deg"]), 2)      if ls.get("ls_inc_deg")      else None,
        "legacy_err":  round(float(ls["ls_inc_err_deg"]), 2)  if ls.get("ls_inc_err_deg")  else None,
    })

# ─── Save CSVs ──────────────────────────────────────────────────────────────────
master_out = pd.concat([master, pd.DataFrame(new_master_rows, columns=master.columns)], ignore_index=True)
master_out.to_csv(os.path.join(BASE, "galfit_vs_legacy_master.csv"), index=False)
print(f"\nMaster CSV: {len(master_out)} rows")

quick_out = pd.concat([quick, pd.DataFrame(new_quick_rows, columns=quick.columns)], ignore_index=True)
quick_out.to_csv(os.path.join(BASE, "galfit_vs_legacy_quick_read.csv"), index=False)
print(f"Quick read CSV: {len(quick_out)} rows")

# ─── Update master_frb_summary ──────────────────────────────────────────────────
for col in ["n_err_psf", "n_err_nopsf"]:
    if col not in summary.columns: summary[col] = ""

for _, row in results12.iterrows():
    frb = row["FRB"]
    idx = summary.index[summary["FRB"] == frb]
    if len(idx) == 0: continue
    i = idx[0]
    ba_psf   = safe(row.get("b_a_psf"))
    ba_nopsf = safe(row.get("b_a_nopsf"))
    summary.at[i, "chi2nu_psf"]   = safe(row.get("chi2nu_psf"))
    summary.at[i, "n_psf"]        = safe(row.get("n_psf"))
    summary.at[i, "n_err_psf"]    = safe(row.get("n_err_psf"))
    summary.at[i, "b_a_psf"]      = ba_psf
    summary.at[i, "inc_psf"]      = round(ba_to_inc(ba_psf), 4) if ba_psf else None
    summary.at[i, "chi2nu_nopsf"] = safe(row.get("chi2nu_nopsf"))
    summary.at[i, "n_nopsf"]      = safe(row.get("n_nopsf"))
    summary.at[i, "n_err_nopsf"]  = safe(row.get("n_err_nopsf"))
    summary.at[i, "b_a_nopsf"]    = ba_nopsf
    summary.at[i, "inc_nopsf"]    = round(ba_to_inc(ba_nopsf), 4) if ba_nopsf else None
    summary.at[i, "Notes"]        = str(summary.at[i, "Notes"]) + "; Free-n GALFIT (n_err available)"

summary.to_csv(os.path.join(BASE, "master_frb_summary.csv"), index=False)
print(f"master_frb_summary updated for {len(FRB12)} FRBs")
print("\n=== DONE ===")
