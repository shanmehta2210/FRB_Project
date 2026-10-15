"""Merge the literature sweep into one wide host table and emit LaTeX.

Inputs
    Publication/materials/host_properties/seed_51.csv
    Publication/materials/host_properties/group_*.csv   (long format, see SWEEP_PROTOCOL.md)
    Publication/materials/host_properties/group_*.bib

Outputs (same folder)
    host_properties_long.csv     every literature row, deduplicated
    host_properties_wide.csv     one row per host, preferred value per quantity
    host_table.tex               AASTeX deluxetable mimicking Bhardwaj+2024 Table 1
    host_table_extended.tex      landscape table with every column that has data
    host_properties.bib          merged BibTeX
    COVERAGE.md                  which quantity is filled for which host, and from where
"""
from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import quad

ROOT = Path(__file__).resolve().parents[1]
HP = ROOT / "Publication" / "materials" / "host_properties"
Q0 = 0.2

# Preferred source order per quantity family. Earlier = preferred.
# Generic fallback (host identity, offsets, photometry, burst properties).
SOURCE_PRIORITY = [
    "gordon2025", "sharma2024", "shannon2025", "connor2025", "law2024", "gordon2023",
    "hussaini2025", "pastormarazuela2026", "li2025", "leung2025", "bhandari2022",
    "heintz2020", "connor2023", "ravi2025", "ravi2023", "caleb2023", "glowacki2025",
    "scott2025", "verdi2026", "bhardwaj2024",
]
# Stellar mass and SFR: one SED code per survey where possible (Prospector for
# DSA-110 via Sharma 2024, Prospector for CRAFT via Gordon 2023), then other
# SED fits, then H-alpha. K-band M/L masses (Gordon 2025) and secondary
# compilations are excluded from the preferred pick but stay in the long file.
SED_PRIORITY = [
    "sharma2024", "gordon2023", "law2024", "muller2026", "leung2025", "pastormarazuela2026",
    "caleb2023", "ravi2025", "connor2023", "sharma2023", "faber2024", "ravi2023",
    "ravi2019", "heintz2020", "bhandari2022", "bhandari2020a", "bhandari2020b",
    "prochaska2019", "bannister2019",
    # Li 2025 Bagpipes re-fits of Legacy photometry: primary dedicated SED papers win; used where none exists
    # (the five Connor 2025 hosts and the 20220912A SFR).
    "li2025", "leewaddell2023", "mahony2018", "muller2026",
]
SED_EXCLUDE_SOURCES = {"bhardwaj2024", "gordon2025", "glowacki2025", "glowacki2023",
                       "athukoralalage2026", "bernalescortes2025", "mannings2021", "chittidi2021"}
SED_EXCLUDE_METHOD_RE = re.compile(r"kband|k-band|compiled|from_|literature|radio|wise|galex|"
                                   r"colour|color|photometric|NED-LVS|published", re.I)
SFR_100MYR_RE = re.compile(r"100\s*myr|0-100", re.I)
# z / DM / DM_MW: localisation catalogues first.
BURST_PRIORITY = [
    "shannon2025", "sharma2024", "law2024", "caleb2023", "connor2025", "pastormarazuela2026",
    "ravi2023", "ravi2025", "ravi2019", "gordon2025", "gordon2023", "hussaini2025",
    "heintz2020", "bhandari2022", "bhandari2020a", "bhandari2020b", "bannister2019",
    "prochaska2019", "macquart2020", "leewaddell2023", "mahony2018", "scott2025", "day2020",
]
# ---- Round-2 audit decisions (audit_group_*.csv / *_notes.md) -------------
# Rows removed from the long file before picking: (frb, quantity, source).
AUDIT_EXCLUDE = {
    ("20190102C", "m_r", "gordon2023"),     # 20.77 is VLT/FORS2 I band (Gordon 2023 Tab. 1 filter column)
    ("20190102C", "m_r", "bhardwaj2024"),   # copied from the above
    ("20220912A", "sfr", "bhardwaj2024"),   # 2.5 is the Hewitt 2024 radio-implied SFR, not Ravi 2023
    ("20200906A", "logmstar", "bhardwaj2024"),  # 10.12 is Bhandari 2022 CIGALE, mis-cited to Gordon 2023
    ("20211212A", "z", "leung2025"),        # 0.0715 from James 2022 'in prep'; Shannon 2025 / Gordon 2023 give 0.0707
    ("20211212A", "z", "roxburgh2025"),
    ("20190523A", "sfr", "li2025"),         # 30.8 +/- 48.4: unconstrained Bagpipes fit; keep the Ravi 2019 limit (audit D)
}
# Rows whose source label was wrong: (frb, quantity, source, value) -> correct source.
AUDIT_RELABEL = {
    ("20220914A", "m_r", "sharma2024", "20.05"): "law2024",   # Law 2024 PS1 values mislabelled; Sharma BASS r below
    ("20220920A", "m_r", "sharma2024", "19.46"): "law2024",
    ("20221012A", "m_r", "sharma2024", "19.67"): "law2024",
}
# Extra rows established by the audit: (frb, quantity, value, ep, em, unit, method, source, location, notes).
AUDIT_ADD = [
    ("20220914A", "m_r", "20.20", "0.19", "0.19", "mag", "BASS r", "sharma2024", "ED Table 1", "audit F"),
    ("20220920A", "m_r", "19.43", "0.20", "0.20", "mag", "BASS r", "sharma2024", "ED Table 1", "audit F"),
    ("20221012A", "m_r", "19.64", "0.19", "0.19", "mag", "BASS r", "sharma2024", "ED Table 1", "audit F"),
]
# Final pick forced for a cell: (frb, quantity) -> (value, source). Used for audited conflicts.
PICK_OVERRIDE = {
    ("20210410D", "dm_obs"): ("578.78", "caleb2023"),   # S/N-maximising DM, Caleb 2023 Tab. 1 (Connor 575 is the structure DM)
    ("20221012A", "dm_obs"): ("441.2", "sherman2024"),  # Sherman 2024 / Law 441.08; Connor 442.2 is older
    ("20171020A", "m_r"): ("15.141", "leewaddell2023"),  # PS1 Petrosian r; 14.97 is photographic SuperCOSMOS R
    # Shannon 2025 Table 1 NE2001 misprints (audit B): recomputed NE2001 at the tabulated coordinates,
    # confirmed by Hussaini 2025 DM_exgal where available.
    ("20230526A", "dm_mw"): ("31.9", "this_audit"),
    ("20231226A", "dm_mw"): ("38.07", "hussaini2025"),
    ("20240208A", "dm_mw"): ("39.0", "this_audit"),
    ("20240304A", "dm_mw"): ("73.76", "hussaini2025"),
    # Audit D: Leung 2025 "i ~ 81" is an undocumented text estimate; Li 2025 tabulates a Sersic b/a value.
    ("20231120A", "inclination_deg"): ("71.74", "li2025"),
    # Poorly localised ASKAP burst: Heintz 2020 Tab. 3 gives 1.7 +/- 19.2 kpc, Macquart 2020 / Bhandari 2020a
    # 3.1 (+15.7 / -3.1) kpc. Both are genuine; the localisation paper's asymmetric value is adopted.
    ("20181112A", "offset_kpc"): ("3.1", "macquart2020"),
}
# Published inclinations / axis ratios: deep VLT GALFIT (Glowacki 2025) over PetroFit on Legacy imaging (Li 2025).
INCL_PRIORITY = ["glowacki2025", "chittidi2021", "yamanaka2024", "bhardwaj2024", "leung2025", "li2025", "glowacki2023"]

# Host-specific overrides after reading the group notes.
Z_OVERRIDE = {
    # Sharma 2024 LRIS spectrum and Verdi 2026 give 0.0368; Connor 2025 / Hussaini 2025
    # list 0.0700 (NED photometric 0.0664). Spectroscopic value adopted, flagged.
    "20231120A": ("0.0368", "sharma2024"),
    # Connor 2025 Table 2 / Verdi 2026 spectroscopic value over the Table 1 photo-z.
    "20240119A": ("0.3760", "connor2025"),
    # Sharma 2024 (measurement paper) 0.2621 +/- 0.0001; Connor / Hussaini / Verdi list 0.2625 (audit D).
    "20231123B": ("0.2621", "sharma2024"),
    # Connor 2025 Table 1 dagger is stale (ESM lists a Keck/LRIS spectrum, 2024 May 3); Hussaini is a compilation.
    "20240229A": ("0.2870", "connor2025"),
}

# Planck 2018 for <DM_cosmic>(z)
H0 = 67.4  # km/s/Mpc
OM = 0.315
OB = 0.0493
FD = 0.844  # diffuse baryon fraction, Macquart+2020


def mean_dm_cosmic(z: float) -> float:
    """Macquart et al. (2020) mean cosmic DM in pc cm^-3."""
    if not np.isfinite(z) or z <= 0:
        return float("nan")
    # 3 c H0 Omega_b / (8 pi G m_p) in pc cm^-3 units: 935.8 * (H0/67.4) * (Ob/0.0493) approx
    pref = 935.8 * (H0 / 67.4) * (OB / 0.0493) * FD / 0.844

    def integrand(zz: float) -> float:
        return (1 + zz) / math.sqrt(OM * (1 + zz) ** 3 + (1 - OM))

    val, _ = quad(integrand, 0, z)
    return pref * val


def median_dm_cosmic(z: float, F: float = 0.2) -> float:
    """Median of the Macquart et al. (2020) DM_cosmic distribution.

    p(Delta) = A Delta^-beta exp(-(Delta^-alpha - C0)^2 / (2 alpha^2 sigma^2)),
    alpha = beta = 3, sigma = F / sqrt(z), C0 fixed so that <Delta> = 1.
    """
    mean = mean_dm_cosmic(z)
    if not np.isfinite(mean):
        return float("nan")
    sigma = F / math.sqrt(z)
    alpha = beta = 3.0
    d = np.linspace(1e-3, 6.0, 60000)

    def pdf(c0: float) -> np.ndarray:
        p = d ** (-beta) * np.exp(-((d ** (-alpha) - c0) ** 2) / (2 * alpha ** 2 * sigma ** 2))
        return p / np.trapezoid(p, d)

    lo, hi = -5.0, 5.0
    for _ in range(60):  # bisection on C0 so that the mean is unity
        mid = 0.5 * (lo + hi)
        m = np.trapezoid(d * pdf(mid), d)
        if m > 1:  # mean falls as C0 rises
            lo = mid
        else:
            hi = mid
    p = pdf(0.5 * (lo + hi))
    cdf = np.cumsum(p) * (d[1] - d[0])
    med = float(d[np.searchsorted(cdf, 0.5)])
    return mean * med


def source_rank(s: str, order: list[str] | None = None) -> int:
    order = order or SOURCE_PRIORITY
    s = str(s).lower().replace(" ", "")
    for k, name in enumerate(order):
        if s.startswith(name):
            return k
    return len(order)


def to_num(x):
    try:
        return float(str(x).replace("<", "").replace(">", "").strip())
    except Exception:
        return float("nan")


def load_long() -> pd.DataFrame:
    frames = []
    for p in sorted(HP.glob("group_*.csv")):
        try:
            d = pd.read_csv(p, dtype=str, encoding="utf-8").fillna("")
        except UnicodeDecodeError:
            d = pd.read_csv(p, dtype=str, encoding="latin-1").fillna("")
        d["group_file"] = p.name
        frames.append(d)
    if not frames:
        return pd.DataFrame(columns=["frb", "quantity", "value", "err_plus", "err_minus", "unit",
                                     "method", "source_short", "arxiv_id", "doi", "location", "notes"])
    L = pd.concat(frames, ignore_index=True)
    L["frb"] = L["frb"].str.strip()
    L["quantity"] = L["quantity"].str.strip().str.lower()
    L = L.drop_duplicates(subset=["frb", "quantity", "value", "source_short"])
    return apply_audits(L)


def apply_audits(L: pd.DataFrame) -> pd.DataFrame:
    """Fold the round-2 audit files into the long table.

    * ``corrected`` rows: the (frb, quantity, current_source) rows they replace are
      dropped, the verified value is appended with its own source.
    * ``new`` / ``conflict`` rows with a value are appended (they are alternatives
      the pick logic can choose from; conflicts are resolved by priority or
      PICK_OVERRIDE).
    * Hard exclusions / relabels / additions from the module-level tables.
    """
    cols = ["frb", "quantity", "value", "err_plus", "err_minus", "unit", "method", "source_short",
            "arxiv_id", "doi", "location", "notes"]
    for c in cols:
        if c not in L.columns:
            L[c] = ""
    extra = []
    for p in sorted(HP.glob("audit_group_*.csv")):
        try:
            a = pd.read_csv(p, dtype=str, encoding="utf-8").fillna("")
        except UnicodeDecodeError:
            a = pd.read_csv(p, dtype=str, encoding="latin-1").fillna("")
        a["frb"] = a.frb.str.strip()
        a["quantity"] = a.quantity.str.strip().str.lower()
        a["verdict"] = a.verdict.str.strip().str.lower()
        for _, r in a.iterrows():
            if r.verdict == "corrected" and r.current_source:
                drop = (L.frb == r.frb) & (L.quantity == r.quantity) & (L.source_short == r.current_source)
                if r.current_value:
                    cv = str(r.current_value).strip()
                    vals = L.value.astype(str).str.strip()
                    same = vals == cv
                    if np.isfinite(to_num(cv)):
                        same |= np.isclose(vals.map(to_num).astype(float), to_num(cv), rtol=0, atol=1e-9)
                    drop &= same
                L = L[~drop]
            if r.verdict in ("corrected", "new", "conflict") and str(r.verified_value).strip() and r.source_short:
                extra.append({"frb": r.frb, "quantity": r.quantity, "value": str(r.verified_value).strip(),
                              "err_plus": r.err_plus, "err_minus": r.err_minus, "unit": r.unit, "method": r.method,
                              "source_short": r.source_short, "arxiv_id": r.arxiv_id, "doi": r.doi,
                              "location": r.location, "notes": f"[{p.stem}:{r.verdict}] {r.notes}",
                              "group_file": p.name})
    for frb, q, v, ep, em, unit, meth, src, loc, note in AUDIT_ADD:
        extra.append({"frb": frb, "quantity": q, "value": v, "err_plus": ep, "err_minus": em, "unit": unit,
                      "method": meth, "source_short": src, "arxiv_id": "", "doi": "", "location": loc,
                      "notes": note, "group_file": "AUDIT_ADD"})
    if extra:
        L = pd.concat([L, pd.DataFrame(extra)], ignore_index=True)
    for (frb, q, src, val), new_src in AUDIT_RELABEL.items():
        m = (L.frb == frb) & (L.quantity == q) & (L.source_short == src) & (L.value.astype(str).str.strip() == val)
        L.loc[m, "source_short"] = new_src
    for frb, q, src in AUDIT_EXCLUDE:
        L = L[~((L.frb == frb) & (L.quantity == q) & (L.source_short == src))]
    L = L.drop_duplicates(subset=["frb", "quantity", "value", "source_short"])
    return L.reset_index(drop=True)


def pick(L: pd.DataFrame, frb: str, q: str):
    """Return (value, err_plus, err_minus, source, method) for the preferred row."""
    sub = L[(L.frb == frb) & (L.quantity == q) & (L.value != "")].copy()
    if sub.empty:
        return None
    src_l = sub.source_short.str.lower().str.replace(" ", "")
    meth = sub.method.fillna("")
    if q in ("logmstar", "sfr"):
        keep = ~src_l.isin(SED_EXCLUDE_SOURCES) & ~meth.str.contains(SED_EXCLUDE_METHOD_RE)
        keep &= ~sub.value.str.strip().str.startswith(("<", ">"))
        keep &= ~meth.str.contains(r"limit", case=False)
        if q == "sfr":
            # Tier 0: explicit 0-100 Myr SED SFR; 1: other SED SFR; 2: H-alpha; 3: short-timescale SED.
            is100 = meth.str.contains(SFR_100MYR_RE)
            is_sed = meth.str.contains(r"sed|prospector|cigale|bagpipes|prospect", case=False)
            is20 = meth.str.contains(r"20\s*myr|30\s*myr", case=False)
            is_ha = meth.str.contains(r"halpha|h_alpha|hbeta", case=False)
            tier = np.where(is100 & ~is20, 0, np.where(is_sed & ~is20, 1, np.where(is_ha, 2, np.where(is20, 3, 4))))
            sub = sub.assign(tier=tier)
        else:
            sub = sub.assign(tier=0)
        sub = sub[keep]
        if sub.empty:
            return None
        sub["rank"] = sub.source_short.map(lambda s: source_rank(s, SED_PRIORITY))
        r = sub.sort_values(["tier", "rank"]).iloc[0]
    elif q in ("z", "dm_obs", "dm_mw", "z_type", "dm_exgal"):
        if q == "dm_mw":
            sub = sub[~meth.str.contains(r"ymw|ne2025|halo|cgm", case=False)]
            if sub.empty:
                return None
        if q == "z":
            sub = sub[~meth.str.contains(r"photo", case=False) | (meth == "")]
            if sub.empty:
                return None
        sub["rank"] = sub.source_short.map(lambda s: source_rank(s, BURST_PRIORITY))
        r = sub.sort_values("rank").iloc[0]
    elif q in ("inclination_deg", "axis_ratio", "position_angle_deg"):
        sub["rank"] = sub.source_short.map(lambda s: source_rank(s, INCL_PRIORITY))
        r = sub.sort_values("rank").iloc[0]
    elif q == "tau_scatt_ms":
        # Homogenise the reference frequency: drop compilations rescaled to 600 MHz (Li 2025:
        # 20210410D 226.9 ms is Caleb's 29.4 ms at 1 GHz x (1/0.6)^4) and prefer values quoted at 1 GHz.
        sub = sub[~meth.str.contains(r"compiled|scaled_to_600", case=False)]
        if sub.empty:
            return None
        txt = (sub.method.fillna("") + " " + sub.notes.fillna("")).str.lower()
        sub = sub.assign(tier=np.where(txt.str.contains(r"1\s*ghz"), 0, 1),
                         rank=sub.source_short.map(source_rank))
        r = sub.sort_values(["tier", "rank"]).iloc[0]
    else:
        sub["rank"] = sub.source_short.map(source_rank)
        r = sub.sort_values("rank").iloc[0]
    return r.value, r.err_plus, r.err_minus, r.source_short, r.method


def fmt(v, nd=2):
    x = to_num(v)
    if not np.isfinite(x):
        return "\\nodata"
    return f"{x:.{nd}f}"


def fmt_err(v, ep, em, nd=2, nd_max=6):
    x = to_num(v)
    if not np.isfinite(x):
        return "\\nodata"
    p, m = to_num(ep), to_num(em)
    if np.isfinite(p) and p > 0:
        # Decimals so that the (smaller) error keeps at least one significant figure, never fewer than nd.
        e_min = min(p, m) if (np.isfinite(m) and m > 0) else p
        nd = max(nd, int(math.ceil(-math.log10(e_min))))
        nd = min(nd, nd_max)
        if p < 0.5 * 10 ** (-nd):
            p = m = float("nan")  # error below the displayed precision: drop it rather than print 0
    if np.isfinite(p) and np.isfinite(m) and p > 0 and m > 0 and abs(p - m) >= 0.5 * 10 ** (-nd):
        return f"${x:.{nd}f}^{{+{p:.{nd}f}}}_{{-{m:.{nd}f}}}$"
    if np.isfinite(p) and p > 0:
        return f"${x:.{nd}f} \\pm {p:.{nd}f}$"
    return f"{x:.{nd}f}"


def bin_of(i: float) -> str:
    if not np.isfinite(i):
        return "\\nodata"
    edges = [41.4, 60.0, 75.5]  # equal cos i quartiles: cos i = 0.75, 0.5, 0.25
    return f"Bin{1 + sum(i > e for e in edges)}"


EXTRA_HOSTS = ["20200430A", "20220105A", "20220717A", "20230708A", "20240310A", "20240201A"]  # b/a <= 0.2


def load_seed() -> pd.DataFrame:
    """51 working-sample hosts plus the six pre-cut hosts (q <= 0.2), flagged."""
    seed = pd.read_csv(HP / "seed_51.csv", dtype={"frb": str})
    seed["in_working_sample"] = True
    # Adopted sigma_q (GALFIT formal error in quadrature with the sky term) from the statistics run.
    cq = pd.read_csv(ROOT / "plots" / "plots_null" / "v2" / "frb_vs_hsc_confirmed50" / "confirmed50_q.csv",
                     dtype={"frb": str}).set_index("frb")
    seed["q_winning_err"] = seed["frb"].map(cq["sigma_q_B"]).fillna(seed["q_winning_err"])
    loc = pd.read_csv(ROOT / "master_frb_localization.csv", dtype={"frb": str}).set_index("frb")
    prod = pd.read_csv(ROOT / "pipeline_galfit_results.csv", dtype={"frb": str}).set_index("frb")
    extra = []
    for f in EXTRA_HOSTS:
        r = {c: np.nan for c in seed.columns}
        r.update({"frb": f, "in_working_sample": False, "in_bhardwaj2024": False})
        if f in loc.index:
            for src, dst in (("z", "z"), ("DM_obs", "DM_obs"), ("DM_MW_NE2001", "DM_MW_NE2001"),
                             ("loc_survey", "loc_survey"), ("repeater", "repeater")):
                if src in loc.columns:
                    r[dst] = loc.loc[f, src]
        if f in prod.index:
            r["m_r_galfit"] = prod.loc[f, "mag"] if "mag" in prod.columns else np.nan
            r["q_galfit_production"] = prod.loc[f, "b_a"] if "b_a" in prod.columns else np.nan
            r["imaging"] = prod.loc[f, "imaging"] if "imaging" in prod.columns else np.nan
        extra.append(r)
    return pd.concat([seed, pd.DataFrame(extra)], ignore_index=True)


def main() -> None:
    seed = load_seed()
    L = load_long()
    L.to_csv(HP / "host_properties_long.csv", index=False, encoding="utf-8")

    quantities = sorted(L.quantity.unique()) if not L.empty else []
    wide_rows, cov_rows = [], []
    for _, s in seed.iterrows():
        f = s.frb
        row = {"frb": f}
        for q in quantities:
            got = pick(L, f, q)
            if got:
                v, ep, em, src, meth = got
                row[q] = v
                row[q + "_ep"] = ep
                row[q + "_em"] = em
                row[q + "_src"] = src
                row[q + "_method"] = meth
                cov_rows.append((f, q, src))
        if f in Z_OVERRIDE:
            row["z"], row["z_src"] = Z_OVERRIDE[f]
            row["z_method"] = "override_spectroscopic (see notes)"
        for (of, oq), (ov, osrc) in PICK_OVERRIDE.items():
            if of != f:
                continue
            hit = L[(L.frb == f) & (L.quantity == oq) & (L.source_short == osrc) &
                    (L.value.astype(str).str.strip() == ov)]
            row[oq], row[oq + "_src"] = ov, osrc
            row[oq + "_ep"] = hit.err_plus.iloc[0] if len(hit) else ""
            row[oq + "_em"] = hit.err_minus.iloc[0] if len(hit) else ""
            row[oq + "_method"] = (hit.method.iloc[0] if len(hit) else "") + " [audit override]"
        # SFR upper limit when no measured SFR exists (e.g. 20190523A: Ravi 2019 < 1.3).
        if not np.isfinite(to_num(row.get("sfr", ""))):
            lim = L[(L.frb == f) & (L.quantity == "sfr") & L.value.astype(str).str.strip().str.startswith("<")
                    & ~L.source_short.str.lower().isin(SED_EXCLUDE_SOURCES)].copy()
            if len(lim):
                lim["rank"] = lim.source_short.map(lambda s_: source_rank(s_, SED_PRIORITY))
                lr = lim.sort_values("rank").iloc[0]
                row["sfr_limit"] = str(lr.value).strip()
                row["sfr_limit_src"] = lr.source_short
                row.pop("sfr", None)
        # fall back to repo values
        z = to_num(row.get("z")) if "z" in row else float("nan")
        if not np.isfinite(z):
            z = to_num(s.z)
            if np.isfinite(z):
                row["z"], row["z_src"] = z, "repo:master_frb_localization"
        dm = to_num(row.get("dm_obs")) if "dm_obs" in row else float("nan")
        if not np.isfinite(dm):
            dm = to_num(s.DM_obs)
            if np.isfinite(dm):
                row["dm_obs"], row["dm_obs_src"] = dm, "repo:master_frb_localization"
        dmmw = to_num(row.get("dm_mw")) if "dm_mw" in row else float("nan")
        if not np.isfinite(dmmw):
            dmmw = to_num(s.DM_MW_NE2001)
            if np.isfinite(dmmw):
                row["dm_mw"], row["dm_mw_src"] = dmmw, "repo:NE2001"
        row["dm_cosmic_mean"] = mean_dm_cosmic(z)
        row["dm_cosmic_median"] = median_dm_cosmic(z)
        ok = np.isfinite(dm) and np.isfinite(dmmw) and np.isfinite(z)
        row["dm_host_obsframe_calc"] = dm - dmmw - row["dm_cosmic_median"] if ok else float("nan")
        row["dm_host_restframe_calc"] = row["dm_host_obsframe_calc"] * (1 + z) if ok else float("nan")
        # our measurements
        row.update({
            "i_galfit_deg": s.i_galfit_deg,
            "i_astrophot_deg": s.i_astrophot_deg,
            "i_isophote_deg": s.i_isophote_deg,
            "q_winning": s.q_winning,
            "q_winning_err": s.q_winning_err,
            "q_galfit_production": s.q_galfit_production,
            "m_r_galfit": s.m_r_galfit,
            "imaging": s.imaging,
            "loc_survey": s.loc_survey,
            "host_citation": s.host_citation,
            "in_bhardwaj2024": s.in_bhardwaj2024,
            "in_working_sample": bool(s.in_working_sample),
            "repeater": s.repeater,
        })
        # sSFR always from the adopted SFR and M* so the three columns are consistent.
        m, sf = to_num(row.get("logmstar", "")), to_num(row.get("sfr", ""))
        if np.isfinite(m) and np.isfinite(sf) and sf > 0:
            row["log_ssfr"] = f"{math.log10(sf) - m:.2f}"
            row["log_ssfr_src"] = f"calc:{row.get('sfr_src', '')}/{row.get('logmstar_src', '')}"
            row["log_ssfr_ep"] = row["log_ssfr_em"] = ""
        else:
            # Never mix a published sSFR with a missing SFR or mass column.
            for k in ("log_ssfr", "log_ssfr_ep", "log_ssfr_em", "log_ssfr_src", "log_ssfr_method"):
                row.pop(k, None)
        wide_rows.append(row)
    W_all = pd.DataFrame(wide_rows)
    W_all.to_csv(HP / "host_properties_wide.csv", index=False, encoding="utf-8")
    W = W_all[W_all.in_working_sample].copy()
    W_extra = W_all[~W_all.in_working_sample].copy()

    # ---- coverage -------------------------------------------------------
    key = ["z", "dm_obs", "dm_mw", "logmstar", "sfr", "log_ssfr", "m_r", "offset_kpc", "re_kpc",
           "dm_host_restframe_calc"]
    lines = ["# Coverage of the 51-host property sweep (+6 pre-cut hosts at the end)", "",
             "| FRB | " + " | ".join(key) + " |", "|---|" + "---|" * len(key)]
    for _, r in pd.concat([W, W_extra]).iterrows():
        cells = []
        for q in key:
            v = r.get(q, "")
            cells.append("x" if (isinstance(v, str) and v != "") or (isinstance(v, float) and np.isfinite(v)) else "")
        tag = "" if r.in_working_sample else " (pre-cut)"
        lines.append(f"| {r.frb}{tag} | " + " | ".join(cells) + " |")
    lines += ["", "Counts (51 working sample):"]
    for q in key:
        n = sum(1 for _, r in W.iterrows() if (q in W.columns) and str(r.get(q, "")) not in ("", "nan"))
        lines.append(f"- {q}: {n}/51")
    (HP / "COVERAGE.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    # ---- bib ---------------------------------------------------------------
    bibs = []
    seen = set()
    for p in sorted(HP.glob("group_*.bib")):
        txt = p.read_text(encoding="utf-8", errors="replace")
        for entry in re.split(r"(?=@\w+\{)", txt):
            m = re.match(r"@\w+\{([^,]+),", entry)
            if m and m.group(1) not in seen and entry.strip():
                seen.add(m.group(1))
                bibs.append(entry.strip())
    (HP / "host_properties.bib").write_text("\n\n".join(bibs) + "\n", encoding="utf-8")

    # ---- Table 1 clone -----------------------------------------------------
    def citekey(src: str) -> str:
        s = str(src).lower().replace(" ", "")
        return re.sub(r"[^a-z0-9]", "", s) or "repo"

    W_sorted = W.copy()
    W_sorted["cosi"] = np.cos(np.radians(W_sorted.i_galfit_deg))
    W_sorted = W_sorted.sort_values("frb")
    t1 = []
    t1.append(r"\begin{longrotatetable}")
    t1.append(r"\begin{deluxetable*}{clccccccccl}")
    t1.append(r"\tabletypesize{\scriptsize}")
    t1.append(r"\tablecaption{Sample of 51 FRB host galaxies employed in this work\label{tab:sample}}")
    t1.append(r"\tablehead{")
    t1.append(r"\colhead{Number\tablenotemark{a}} & \colhead{FRB} & \colhead{Inclination angle\tablenotemark{b}} & \colhead{Bin \#} & "
              r"\colhead{SFR$_{0-100\,\rm Myr}$\tablenotemark{c}} & \colhead{log(sSFR)} & \colhead{Mass} & \colhead{$z$} & "
              r"\colhead{DM$_{\rm host}$\tablenotemark{d}} & \colhead{Host $m_r$\tablenotemark{e}} & \colhead{Reference} \\")
    t1.append(r"\colhead{} & \colhead{} & \colhead{Degrees} & \colhead{} & \colhead{M$_\odot$\,yr$^{-1}$} & \colhead{log(yr$^{-1}$)} & "
              r"\colhead{M$_\odot \times 10^{10}$} & \colhead{} & \colhead{pc\,cm$^{-3}$} & \colhead{AB mag} & \colhead{}")
    t1.append(r"}")
    t1.append(r"\startdata")
    for k, (_, r) in enumerate(W_sorted.iterrows(), 1):
        num = f"\\textbf{{{k}}}" if "PS1" in str(r.imaging) or "HSC" in str(r.imaging) else str(k)
        inc = f"{to_num(r.i_galfit_deg):.0f}, {to_num(r.i_astrophot_deg):.0f}" if np.isfinite(to_num(r.i_astrophot_deg)) else f"{to_num(r.i_galfit_deg):.0f}, \\nodata"
        sfr = fmt(r.get("sfr"), 2)
        sm = str(r.get("sfr_method", "")).lower()
        if sfr != "\\nodata":
            if re.search(r"halpha|hbeta", sm):
                sfr += "$^{\\dagger}$"
            elif re.search(r"20\s*myr|30\s*myr", sm):
                sfr += "$^{\\ddagger}$"
        ssfr = fmt(r.get("log_ssfr"), 2)
        lm = to_num(r.get("logmstar"))
        lim = str(r.get("sfr_limit", "")).strip()
        if sfr == "\\nodata" and lim.startswith("<") and np.isfinite(to_num(lim[1:])):
            sfr = f"$<{to_num(lim[1:]):.2f}$"
            if np.isfinite(lm):
                ssfr = f"$<{math.log10(to_num(lim[1:])) - lm:.2f}$"
        mass = f"{10 ** lm / 1e10:.2f}" if np.isfinite(lm) else "\\nodata"
        z = fmt(r.get("z"), 4)
        dmh = to_num(r.get("dm_host_restframe_calc"))
        dmh_s = f"{dmh:.0f}" if np.isfinite(dmh) else "\\nodata"
        mr = fmt(r.m_r_galfit, 2)
        refs = []
        for c in ("logmstar_src", "sfr_src", "sfr_limit_src", "z_src"):
            v = str(r.get(c, "")).strip()
            if v and not v.startswith("repo") and v != "nan":
                k = citekey(v)
                if k not in refs:
                    refs.append(k)
        ref = "; ".join(f"\\citealt{{{x}}}" for x in refs) if refs else "\\nodata"
        t1.append(f"{num} & {r.frb} & {inc} & {bin_of(to_num(r.i_galfit_deg))} & {sfr} & {ssfr} & {mass} & {z} & {dmh_s} & {mr} & {ref} \\\\")
    t1.append(r"\enddata")
    t1.append(r"\tablenotetext{a}{Boldface indicates FRB hosts for which Pan-STARRS1 DR1 (or, for FRB 20240208A, HSC-SSP PDR2) $r$-band images were used. For all other hosts we use Legacy Surveys DR10 $r$-band images (DR9 for FRBs 20190523A and 20230907D).}")
    t1.append(r"\tablenotetext{b}{The first inclination angle in each row is from the GALFIT single-S\'ersic fit; the second is from the independent AstroPhot fit of the same image. Both use $\cos^2 i = (q^2 - q_0^2)/(1 - q_0^2)$ with $q_0 = 0.2$. Bins are equal-width in $\cos i$ (edges at $i = 41.4^\circ$, $60.0^\circ$, $75.5^\circ$).}")
    t1.append(r"\tablenotetext{c}{SFR$_{0-100\,\rm Myr}$ and stellar mass are SED-fit values taken from the references in the last column (Prospector for the DSA-110 and CRAFT hosts; CIGALE or Bagpipes where no Prospector fit exists). log(sSFR) is computed from the two tabulated values. $\dagger$: H$\alpha$-based SFR (no SED SFR published); $\ddagger$: SED SFR over the last 20 Myr only. Every published alternative is listed in Tables~\ref{tab:extended_a} and \ref{tab:extended_b}.}")
    t1.append(r"\tablenotetext{d}{Rest-frame host contribution, $[{\rm DM} - {\rm DM}_{\rm MW,NE2001} - \widetilde{\rm DM}_{\rm cosmic}(z)]\,(1+z)$, where $\widetilde{\rm DM}_{\rm cosmic}$ is the median of the \citet{macquart2020} distribution ($F = 0.2$, Planck 2018). No halo term is subtracted, following \citet{bhardwaj2024}.}")
    t1.append(r"\tablenotetext{e}{GALFIT $r$-band magnitude from this work (AB).}")
    t1.append(r"\end{deluxetable*}")
    t1.append(r"\end{longrotatetable}")
    (HP / "host_table.tex").write_text("\n".join(t1) + "\n", encoding="utf-8")

    # ---- 57-host table in exactly the Bhardwaj et al. (2024) Table 1 layout (no errors) ----
    # Bins: quartiles of the control sample's cos i (same pool as the statistics and Figure 2),
    # numbered face-on (1) to edge-on (4), assigned from the unrounded cos i of the winning q.
    import sys
    sys.path.insert(0, str(ROOT / "scripts"))
    from plot_frb_vs_hsc_confirmed50 import load_hsc_pool  # noqa: E402
    hsc_cosi, n_hsc = load_hsc_pool(mag_limit=22.0, q0=Q0)
    q_edges = np.quantile(hsc_cosi, [0.25, 0.5, 0.75])          # cos i, increasing
    i_edges = np.degrees(np.arccos(q_edges))[::-1]              # deg, increasing: 41.0, 54.9, 67.1

    def cosi_of(q) -> float:
        q = to_num(q)
        if not np.isfinite(q) or q <= Q0:
            return float("nan")
        return math.sqrt((q * q - Q0 ** 2) / (1 - Q0 ** 2))

    def quartile_bin(cosi: float) -> str:
        if not np.isfinite(cosi):
            return "\\nodata"
        return str(1 + int(np.sum(cosi < q_edges)))              # cos i above all edges -> 1 (face-on)

    def t57_row(k: int, r, working: bool) -> str:
        num = f"\\textbf{{{k}}}" if "PS1" in str(r.imaging) or "HSC" in str(r.imaging) else str(k)
        cq = cosi_of(r.q_winning if working else r.get("q_galfit_production"))
        if working and np.isfinite(cq):
            inc = f"{math.degrees(math.acos(cq)):.2f}"
            bn = quartile_bin(cq)
        else:
            inc = f"$q = {to_num(r.get('q_galfit_production')):.2f}$"
            bn = "\\nodata"
        sfr = fmt(r.get("sfr"), 2)
        sm = str(r.get("sfr_method", ""))
        # Mark SED SFRs that are not a 0-100 Myr average (instantaneous, or timescale not stated).
        is100 = re.search(r"100\s*myr|0-100", sm, re.I) and not re.search(
            r"instantaneous|not stated|not a 100", sm, re.I)
        if sfr != "\\nodata" and not is100:
            sfr += "$^{\\dagger}$"
        ssfr = fmt(r.get("log_ssfr"), 2)
        lm = to_num(r.get("logmstar"))
        lim = str(r.get("sfr_limit", "")).strip()
        if sfr == "\\nodata" and lim.startswith("<") and np.isfinite(to_num(lim[1:])):
            sfr = f"$<{to_num(lim[1:]):.2f}$"
            if np.isfinite(lm):
                ssfr = f"$<{math.log10(to_num(lim[1:])) - lm:.2f}$"
        mass = f"{10 ** lm / 1e10:.2f}" if np.isfinite(lm) else "\\nodata"
        dmh = to_num(r.get("dm_host_restframe_calc"))
        refs = []
        for c in ("logmstar_src", "sfr_src", "sfr_limit_src", "z_src"):
            v = str(r.get(c, "")).strip()
            if v and not v.startswith("repo") and v != "nan" and citekey(v) not in refs:
                refs.append(citekey(v))
        ref = "; ".join(f"\\citealt{{{x}}}" for x in refs) if refs else "\\nodata"
        return (f"{num} & {r.frb} & {inc} & {bn} & {sfr} & {ssfr} & {mass} & {fmt(r.get('z'), 4)} & "
                f"{f'{dmh:.0f}' if np.isfinite(dmh) else chr(92) + 'nodata'} & {fmt(r.m_r_galfit, 2)} & {ref} \\\\")

    # Inclination uncertainties for note b: sigma_q (formal + sky) propagated through the Hubble relation,
    # then the 5 deg floor added in quadrature (Section on errors).
    qw = pd.to_numeric(W_sorted.q_winning, errors="coerce").to_numpy()
    sq = pd.to_numeric(W_sorted.q_winning_err, errors="coerce").to_numpy()
    cw = np.sqrt((qw ** 2 - Q0 ** 2) / (1 - Q0 ** 2))
    didq = np.degrees(qw / ((1 - Q0 ** 2) * cw * np.sqrt(1 - cw ** 2)))
    sig_formal = sq * didq
    sig_tot = np.sqrt(sig_formal ** 2 + 5.0 ** 2)
    sig_note = (f"Propagated $q$ uncertainties (formal plus sky term) give a median $\\sigma_i = "
                f"{np.nanmedian(sig_formal):.1f}^\\circ$; with the $5^\\circ$ systematic term added in quadrature, "
                f"the adopted uncertainties range from ${np.nanmin(sig_tot):.1f}^\\circ$ to "
                f"${np.nanmax(sig_tot):.1f}^\\circ$ (median ${np.nanmedian(sig_tot):.1f}^\\circ$).")

    t57 = [r"\begin{longrotatetable}", r"\begin{deluxetable*}{clccccccccl}", r"\tabletypesize{\scriptsize}",
           r"\tablecaption{FRB host galaxies considered in this work\label{tab:sample57}}",
           r"\tablehead{",
           r"\colhead{Number\tablenotemark{a}} & \colhead{FRB} & \colhead{Inclination angle\tablenotemark{b}} & "
           r"\colhead{Bin \#} & \colhead{SFR$_{0-100\,\rm Myr}$\tablenotemark{c}} & \colhead{log(sSFR)} & "
           r"\colhead{Mass} & \colhead{$z$} & \colhead{DM$_{\rm host}$} & \colhead{Host m$_r$\tablenotemark{d}} & "
           r"\colhead{Reference} \\",
           r"\colhead{} & \colhead{} & \colhead{Degrees} & \colhead{} & \colhead{M$_\odot$\,yr$^{-1}$} & "
           r"\colhead{log(yr$^{-1}$)} & \colhead{M$_\odot \times 10^{10}$} & \colhead{} & \colhead{pc\,cm$^{-3}$} & "
           r"\colhead{AB mag} & \colhead{}",
           r"}", r"\startdata"]
    k = 0
    for _, r in W_sorted.iterrows():
        k += 1
        t57.append(t57_row(k, r, True))
    t57.append(r"\cutinhead{Hosts with $q \leq q_0$, excluded from the analysis (Section~\ref{sec:inclination})}")
    for _, r in W_extra.sort_values("frb").iterrows():
        k += 1
        t57.append(t57_row(k, r, False))
    e1, e2, e3 = i_edges
    t57 += [r"\enddata",
            r"\tablenotetext{a}{Boldface indicates FRB hosts for which Pan-STARRS1 DR1 (or, for FRB 20240208A, "
            r"HSC-SSP PDR2) $r$-band images were used. For all other hosts we use Legacy Surveys DR10 $r$-band images "
            r"(DR9 for FRBs 20190523A and 20230907D).}",
            r"\tablenotetext{b}{Inclination from our \texttt{GALFIT} single-S\'ersic axis ratio $q$ "
            r"(Equation~\ref{eq:HubbleFormula}, $q_0 = 0.2$). Bins are the quartiles of the control sample's "
            r"$\cos i$ distribution, numbered from face-on to edge-on (edges at "
            f"$i = {e1:.2f}^\\circ$, ${e2:.2f}^\\circ$, and ${e3:.2f}^\\circ$; Figure~\\ref{{fig:bins}}). "
            r"For the six hosts with $q \leq q_0$, the inclination is undefined and we list $q$: five are marginally "
            r"resolved, and their fits collapse to unphysically thin profiles, while FRB 20240201A is a resolved "
            r"edge-on disk. " + sig_note + "}",
            r"\tablenotetext{c}{SFR and stellar mass are SED-fit values from the references in the last column; "
            r"$\log({\rm sSFR})$ is computed from the two tabulated values. SFRs are averaged over the last "
            r"0--100\,Myr, except "
            r"those marked $\dagger$, for which the published SED SFR is instantaneous or its timescale is not stated. "
            r"DM$_{\rm host}$ is the rest-frame host contribution "
            r"$[{\rm DM} - {\rm DM}_{\rm MW} - \widetilde{\rm DM}_{\rm cosmic}(z)]\,(1+z)$, with DM$_{\rm MW}$ from "
            r"NE2001 \citep{cordes2002} and $\widetilde{\rm DM}_{\rm cosmic}$ the median of the \citet{macquart2020} "
            r"distribution ($F = 0.2$, Planck 2018 cosmology; \citealt{planck2020}); no halo term is subtracted, "
            r"similar to \citet{Bhardwaj2024}. Negative values reflect the scatter of DM$_{\rm cosmic}$ about its "
            r"median.}",
            r"\tablenotetext{d}{\texttt{GALFIT} $r$-band magnitude from this work.}",
            r"\end{deluxetable*}", r"\end{longrotatetable}"]
    (HP / "host_table_57.tex").write_text("\n".join(t57) + "\n", encoding="utf-8")

    # Portrait version: same rows and notes, compact headers, \startlongtable instead of longrotatetable,
    # numbered references (order of first appearance) collected in \tablerefs, negatives in math mode.
    body_start = t57.index(r"\startdata")
    body_end = t57.index(r"\enddata")
    ref_order: list[str] = []

    def portrait_row(line: str) -> str:
        if not re.search(r"^(\\textbf\{)?\d+\}? & 20", line):
            return line
        cells = [c.strip() for c in line[:-3].split(" & ")]  # strip trailing ' \\'
        keys = re.findall(r"\\citealt\{([^}]+)\}", cells[-1])
        nums = []
        for kk in keys:
            if kk not in ref_order:
                ref_order.append(kk)
            nums.append(str(ref_order.index(kk) + 1))
        cells[-1] = ", ".join(nums) if nums else "\\nodata"
        for j in (5, 8):  # log(sSFR), DM_host
            if re.fullmatch(r"-[\d.]+", cells[j]):
                cells[j] = f"${cells[j]}$"
        return " & ".join(cells) + " \\\\"

    body = [portrait_row(ln) for ln in t57[body_start:body_end + 1]]
    t57p = [r"\startlongtable", r"\begin{deluxetable*}{rlccccccccl}", r"\tabletypesize{\scriptsize}",
            r"\tablecaption{FRB host galaxies considered in this work\label{tab:sample57}}",
            r"\tablehead{",
            r"\colhead{No.\tablenotemark{a}} & \colhead{FRB} & \colhead{$i$\tablenotemark{b}} & \colhead{Bin} & "
            r"\colhead{SFR\tablenotemark{c}} & \colhead{$\log({\rm sSFR})$} & \colhead{Mass} & \colhead{$z$} & "
            r"\colhead{DM$_{\rm host}$} & \colhead{Host $m_r$\tablenotemark{d}} & \colhead{Ref.} \\",
            r"\colhead{} & \colhead{} & \colhead{(deg)} & \colhead{} & \colhead{($M_\odot$\,yr$^{-1}$)} & "
            r"\colhead{(yr$^{-1}$)} & \colhead{($10^{10}\,M_\odot$)} & \colhead{} & \colhead{(pc\,cm$^{-3}$)} & "
            r"\colhead{(AB mag)} & \colhead{}",
            r"}"]
    t57p += body
    notes = [ln for ln in t57[body_end + 1:] if ln not in (r"\end{longrotatetable}", r"\end{deluxetable*}")]
    t57p += notes
    t57p.append(r"\tablerefs{" + "; ".join(f"({n}) \\citet{{{kk}}}" for n, kk in enumerate(ref_order, 1)) + ".}")
    t57p.append(r"\end{deluxetable*}")
    (HP / "host_table_57_portrait.tex").write_text("\n".join(t57p) + "\n", encoding="utf-8")
    bins57 = [quartile_bin(cosi_of(q)) for q in W_sorted.q_winning]
    print("57-table bins (face-on -> edge-on):", [bins57.count(str(b)) for b in (1, 2, 3, 4)],
          "edges", np.round(i_edges, 1), f"HSC N={n_hsc}")

    # ---- Extended table ----------------------------------------------------
    def has_data(c: str) -> bool:
        return c in W.columns and W[c].astype(str).str.strip().replace("nan", "").ne("").any()

    # Quantities we measure ourselves (m_r, q, i, DM_host) are always our values; the literature is used
    # only for what we do not measure.
    ext_parts = {
        "A": ("Burst and dispersion-measure properties of the 51 hosts", "tab:extended_a",
              [c for c in ["host_name", "z", "dm_obs", "dm_mw", "dm_host_restframe_calc", "tau_scatt_ms",
                           "rm_rad_m2", "path_posterior", "offset_kpc", "offset_normalized"] if has_data(c)],
              r"\tablecomments{$z$, DM, DM$_{\rm MW}$ (NE2001), $\tau$, RM, the PATH posterior and the projected "
              r"offset are literature values, one preferred source per cell (source, method and every alternative "
              r"value are listed in \texttt{host\_properties\_wide.csv} and \texttt{host\_properties\_long.csv}). "
              r"DM$_{\rm host}$ is our homogeneous rest-frame estimate, "
              r"$[{\rm DM} - {\rm DM}_{\rm MW} - \widetilde{\rm DM}_{\rm cosmic}(z)](1+z)$, with "
              r"$\widetilde{\rm DM}_{\rm cosmic}$ the median of the \citet{macquart2020} distribution ($F = 0.2$, "
              r"Planck 2018), as in Table~\ref{tab:sample}. $\tau$ is quoted at 1\,GHz where the source gives it and "
              r"otherwise at the band centre of the detection; compilations rescaled to other frequencies are not used.}"),
        "B": ("Stellar-population and structural properties of the 51 hosts", "tab:extended_b",
              [c for c in ["m_r_galfit", "q_winning", "i_galfit_deg", "logmstar", "sfr", "tm_gyr",
                           "metallicity_logz", "av", "u_r_color", "re_kpc", "galaxy_type"] if has_data(c)],
              r"\tablecomments{$m_r$, $q = b/a$ and $i$ are our GALFIT measurements (Section~\ref{sec:methods}); "
              r"$i$ follows from $q$ with $q_0 = 0.2$. The stellar-population quantities, $R_e$ and the host type are "
              r"literature values, one preferred source per cell (source, method and every alternative value in "
              r"\texttt{host\_properties\_wide.csv} and \texttt{host\_properties\_long.csv}). Types are normalised to "
              r"star-forming / transition / quiescent; qualifiers from the source are kept in parentheses.}"),
    }
    hdr = {"host_name": "Host", "z": "$z$", "dm_obs": "DM", "dm_mw": "DM$_{\\rm MW}$",
           "dm_host_restframe_calc": "DM$_{\\rm host}$ (this work)",
           "tau_scatt_ms": "$\\tau$ (ms)", "rm_rad_m2": "RM", "m_r_galfit": "$m_r$ (this work)",
           "q_winning": "$q$ (this work)", "i_galfit_deg": "$i$ (this work)", "logmstar": "log $M_\\star$",
           "sfr": "SFR", "log_ssfr": "log sSFR", "age_mass_weighted_gyr": "$t_m$ (Gyr)", "tm_gyr": "$t_m$ (Gyr)",
           "metallicity_logz": "log $Z$", "av": "$A_V$", "u_r_color": "$u-r$", "re_kpc": "$R_e$ (kpc)",
           "offset_kpc": "Offset (kpc)", "offset_normalized": "$\\delta R/R_e$", "path_posterior": "$P(O|x)$",
           "galaxy_type": "Type", "environment": "Env."}

    def tex_text(v) -> str:
        s = str(v)
        if s in ("", "nan"):
            return "\\nodata"
        return s.replace("_", "\\_").replace("&", "\\&").replace("%", "\\%")

    def norm_type(v) -> str:
        s = str(v).strip()
        if s in ("", "nan"):
            return "\\nodata"
        low = s.lower()
        if low in ("sf", "star-forming", "star forming"):
            return "star-forming"
        if low in ("sf (y)",):
            return "star-forming"
        if low in ("transition", "transitioning"):
            return "transition"
        if low in ("quiescent",):
            return "quiescent"
        if "quiescent" in low:
            rest = re.sub(r"\s*\(quiescent-like\)|quiescent", "", s, flags=re.I).strip(" ,;")
            return "quiescent" + (f" ({tex_text(rest)})" if rest else "")
        if low.startswith(("sf", "star-forming", "star forming")):
            rest = re.sub(r"^(sf|star[- ]forming)[,\s]*", "", s, flags=re.I).strip()
            return "star-forming" + (f" ({tex_text(rest)})" if rest else "")
        return tex_text(s)

    ex: list[str] = []
    for part, (caption, label, cols, comment) in ext_parts.items():
        ex += [r"\begin{longrotatetable}", r"\begin{deluxetable*}{l" + "c" * len(cols) + "}",
               r"\tabletypesize{\scriptsize}",
               f"\\tablecaption{{{caption}\\label{{{label}}}}}",
               r"\tablehead{\colhead{FRB} & " + " & ".join(f"\\colhead{{{hdr.get(c, c)}}}" for c in cols) + "}",
               r"\startdata"]
        for _, r in W_sorted.iterrows():
            cells = []
            for c in cols:
                v = r.get(c, "")
                if c == "galaxy_type":
                    cells.append(norm_type(v))
                elif c in ("host_name", "environment"):
                    cells.append(tex_text(v))
                elif c == "dm_host_restframe_calc":
                    cells.append(fmt(v, 0))
                elif c == "q_winning":
                    cells.append(fmt_err(v, r.get("q_winning_err", ""), r.get("q_winning_err", ""), 2, nd_max=3))
                elif c == "i_galfit_deg":
                    cells.append(fmt(v, 0))
                elif c == "m_r_galfit":
                    cells.append(fmt(v, 2))
                else:
                    nd = 4 if c == "z" else 2
                    cells.append(fmt_err(v, r.get(c + "_ep", ""), r.get(c + "_em", ""), nd))
            ex.append(f"{r.frb} & " + " & ".join(cells) + " \\\\")
        ex += [r"\enddata", comment, r"\end{deluxetable*}", r"\end{longrotatetable}", ""]
    (HP / "host_table_extended.tex").write_text("\n".join(ex) + "\n", encoding="utf-8")

    # ---- the six hosts removed by the b/a > 0.2 cut ------------------------
    sx = [r"\begin{deluxetable*}{lcccccccl}", r"\tabletypesize{\scriptsize}",
          r"\tablecaption{Hosts that pass the magnitude and verification cuts but have $b/a \le 0.2$"
          r"\label{tab:precut}}",
          r"\tablehead{\colhead{FRB} & \colhead{$q_{\rm GALFIT}$} & \colhead{$z$} & \colhead{DM} & "
          r"\colhead{DM$_{\rm MW}$} & \colhead{log $M_\star$} & \colhead{SFR} & \colhead{$(b/a)_{\rm lit}$ / $i_{\rm lit}$} & "
          r"\colhead{Reference} \\ \colhead{} & \colhead{} & \colhead{} & \colhead{pc\,cm$^{-3}$} & "
          r"\colhead{pc\,cm$^{-3}$} & \colhead{M$_\odot$} & \colhead{M$_\odot$\,yr$^{-1}$} & \colhead{} & \colhead{}}",
          r"\startdata"]
    for _, r in W_extra.sort_values("frb").iterrows():
        q = fmt(r.get("q_galfit_production"), 2)
        lit = []
        if str(r.get("axis_ratio", "")) not in ("", "nan"):
            lit.append(f"{fmt(r.axis_ratio, 2)}")
        if str(r.get("inclination_deg", "")) not in ("", "nan"):
            lit.append(f"{fmt(r.inclination_deg, 0)}$^\\circ$")
        refs = []
        for c in ("logmstar_src", "sfr_src", "sfr_limit_src", "z_src"):
            v = str(r.get(c, "")).strip()
            if v and not v.startswith("repo") and v != "nan" and citekey(v) not in refs:
                refs.append(citekey(v))
        ref = "; ".join(f"\\citealt{{{x}}}" for x in refs) if refs else "\\nodata"
        sx.append(f"{r.frb} & {q} & {fmt(r.get('z'), 4)} & {fmt(r.get('dm_obs'), 1)} & {fmt(r.get('dm_mw'), 0)} & "
                  f"{fmt(r.get('logmstar'), 2)} & {fmt(r.get('sfr'), 2)} & {' / '.join(lit) if lit else chr(92) + 'nodata'} & {ref} \\\\")
    sx += [r"\enddata",
           r"\tablecomments{Five of these are unresolved fits in which GALFIT collapses to a sliver (20200430A, "
           r"20220105A, 20220717A, 20230708A, 20240310A); FRB 20240201A is a resolved edge-on disk ($q = 0.18$) that "
           r"re-enters the sample for $q_0 \le 0.15$ (Appendix). Literature axis ratios are the VLT/FORS2 GALFIT values "
           r"of \citet{glowacki2025}; literature inclinations from \citet{li2025} or \citet{leung2025}.}",
           r"\end{deluxetable*}"]
    (HP / "host_table_precut.tex").write_text("\n".join(sx) + "\n", encoding="utf-8")

    # ---- Accounting table: the 57 hosts with the Table-1 quantities and their sources ----
    acc_cols = ["frb", "in_working_sample", "in_bhardwaj2024", "repeater", "imaging", "q_winning", "q_winning_err",
                "i_galfit_deg", "i_astrophot_deg", "z", "z_src", "z_type", "dm_obs", "dm_obs_src", "dm_mw", "dm_mw_src",
                "dm_cosmic_median", "dm_host_restframe_calc", "logmstar", "logmstar_ep", "logmstar_em", "logmstar_src",
                "sfr", "sfr_ep", "sfr_em", "sfr_src", "sfr_method", "sfr_limit", "sfr_limit_src", "log_ssfr",
                "m_r_galfit", "m_r", "m_r_src", "offset_kpc", "offset_kpc_src", "inclination_deg", "inclination_deg_src"]
    acc = W_all.reindex(columns=[c for c in acc_cols if c in W_all.columns])
    acc = acc.sort_values(["in_working_sample", "frb"], ascending=[False, True])
    acc.to_csv(HP / "host_properties_57.csv", index=False, encoding="utf-8")

    # ---- ApJL wrapper: all 57 hosts (51 working + 6 pre-cut), no compile needed ----
    apjl = r"""% ApJL-format host tables: 51 working-sample hosts (Table 1 layout of Bhardwaj et al. 2024),
% the 6 hosts removed by the b/a <= 0.2 cut, and the extended published-property tables.
% Generated by scripts/build_host_property_table.py from host_properties_wide.csv.
% \input{} the individual table files into the manuscript, or compile this file as-is:
%   pdflatex host_tables_apjl ; bibtex host_tables_apjl ; pdflatex host_tables_apjl x2
\documentclass[apjl, linenumbers]{aastex701}
\usepackage{amsmath}

\begin{document}

\title{Host-galaxy properties of the FRB sample}

\author{}
\affiliation{}
\email{}

\section{Working sample}
% Table 1: the 51 hosts in the working sample.
\input{host_table}

\section{Hosts removed by the axis-ratio cut}
% The 6 hosts that pass the magnitude and verification cuts but have b/a <= 0.2.
\input{host_table_precut}

\section{Extended published properties}
% Tables A (burst and DM properties) and B (SED and morphology) for the 51 hosts.
\input{host_table_extended}

\bibliographystyle{aasjournalv7}
\bibliography{host_properties}
\end{document}
"""
    (HP / "host_tables_apjl.tex").write_text(apjl, encoding="utf-8")
    apjl57 = r"""% ApJL-format table of all 57 hosts (51 working sample + 6 with b/a <= 0.2) in exactly the
% Table 1 layout of Bhardwaj et al. (2024): no uncertainties; inclination, bin and m_r are our
% GALFIT measurements; DM_host is our homogeneous rest-frame estimate.
% Generated by scripts/build_host_property_table.py.
\documentclass[apjl, linenumbers]{aastex701}
\usepackage{amsmath}

\begin{document}

\title{Sample of 57 FRB host galaxies}

\author{}
\affiliation{}
\email{}

\input{host_table_57}

\bibliographystyle{aasjournalv7}
\bibliography{host_properties}
\end{document}
"""
    (HP / "host_table_57_apjl.tex").write_text(apjl57, encoding="utf-8")
    print(f"long rows {len(L)}, quantities {len(quantities)}, wide cols {W.shape[1]}")


if __name__ == "__main__":
    main()
