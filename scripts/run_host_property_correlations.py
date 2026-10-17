"""
Host and burst properties vs inclination for the confirmed FRB-host sample.

The property list follows the relations tested in the literature:
  Marnoch+2026      DM_host (Macquart, obs / rest frame) vs q and i; m<21, n<2 subsamples;
                    DM_host vs offset
  Bhardwaj+2024     DM_host vs i
  Glowacki+2025     tau vs i, stellar surface density, mass-weighted age, metallicity
  Leung+2025        DM_host vs M*
  Li+2025           DM_exc vs SFR, sSFR, M*, inclination
  Verdi+2026        |RM|, DM_host vs z
  Day+2020          DM vs tau
plus selection / systematics checks (z, m_r, R_e, Sersic n, S/N, fluence, width, PATH P,
RA, Dec, |b|, DM_MW, survey).

Inclination: winning-leg q -> Hubble cos i (q0 = 0.2). Uncertainty: Protocol B Monte Carlo
(q ~ N(q, sigma_q_B), i += N(0, 5 deg)) on the Spearman coefficient.

Writes plots/plots_null/v2/host_correlations/{host_props.csv, correlations.csv,
CORRELATIONS.md, dmhost_vs_inclination.png}.

Run from repo root::

    python scripts/run_host_property_correlations.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy import stats

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

import build_host_property_table as B  # noqa: E402
from plot_frb_inflated_error_cdf import _hubble_inc_deg  # noqa: E402

ROOT = _SCRIPTS.parent
SAMPLE = ROOT / "plots" / "plots_null" / "v2" / "frb_vs_hsc_confirmed50" / "confirmed50_q.csv"
OUT = ROOT / "plots" / "plots_null" / "v2" / "host_correlations"
Q0 = 0.2
FLOOR_DEG = 5.0
N_DRAWS = 2000
SEED = 20261018

SURVEY_FILL = {
    "20231210F": "MeerKAT", "20240203D": "DSA-110", "20240224A": "DSA-110", "20240312D": "ASKAP",
    "20240525A": "ASKAP", "20240615B": "ASKAP", "20250316A": "CHIME",
}


def num(x) -> float:
    s = str(x).strip()
    if s.startswith(("<", ">")):
        return float("nan")
    return B.to_num(s)


def picked(L: pd.DataFrame, frb: str, q: str):
    if (frb, q) in B.PICK_OVERRIDE:
        v, src = B.PICK_OVERRIDE[(frb, q)]
        return v, src
    if q == "z" and frb in B.Z_OVERRIDE:
        return B.Z_OVERRIDE[frb]
    got = B.pick(L, frb, q)
    return (got[0], got[3]) if got else (None, None)


def any_value(L: pd.DataFrame, frb: str, q: str):
    """Fallback for quantities the SED pick rejects (K-band / compiled masses)."""
    sub = L[(L.frb == frb) & (L.quantity == q)]
    sub = sub[np.isfinite(sub.value.map(num).astype(float))]
    if sub.empty:
        return None, None
    sub = sub.assign(rank=sub.source_short.map(B.source_rank)).sort_values("rank")
    return sub.value.iloc[0], sub.source_short.iloc[0] + " (fallback)"


# Same quantity under different names in the group files; first source by priority wins.
MERGED = {
    "av": ["av"],
    "av_old": ["av_old"],
    "g_r_color": ["g_r_color"],
    "rm_host_rest": ["rm_host_rad_m2", "rm_host"],
    "b_par_host": ["b_parallel_host_ug"],
    "lin_pol": ["linear_pol_frac", "pol_frac_linear", "pol_linear_frac"],
    "circ_pol": ["circular_pol_frac", "pol_frac_circular"],
    "metal_gas": ["metallicity_12logoh"],
    "halpha_ew": ["halpha_ew_a"],
    "q_pah": ["q_pah"],
    "log_tau_rest": ["log_tau_1ghz_rest_ms", "log_tau_rest_1ghz_ms"],
    "scint_bw": ["scint_bandwidth_mhz"],
}


def merged_value(L: pd.DataFrame, frb: str, names: list[str]):
    sub = L[(L.frb == frb) & L.quantity.isin(names)]
    sub = sub[np.isfinite(sub.value.map(num).astype(float))]
    if sub.empty:
        return float("nan"), None
    sub = sub.assign(rank=sub.source_short.map(B.source_rank)).sort_values("rank")
    return num(sub.value.iloc[0]), sub.source_short.iloc[0]


def survey_group(s: str) -> str:
    s = str(s)
    for k in ("ASKAP", "DSA", "MeerKAT", "CHIME"):
        if k.lower() in s.lower():
            return "DSA-110" if k == "DSA" else k
    return "other"


def build_table() -> pd.DataFrame:
    hosts = pd.read_csv(SAMPLE, dtype={"frb": str})
    L = B.load_long()
    master = pd.read_csv(ROOT / "master_frb_localization.csv", dtype={"frb": str}).drop_duplicates("frb").set_index("frb")
    prod = pd.read_csv(ROOT / "pipeline_galfit_results.csv", dtype={"frb": str}).drop_duplicates("frb").set_index("frb")
    rows = []
    for _, h in hosts.iterrows():
        f = h.frb
        r = {"frb": f, "q": float(h.q), "sigma_q_B": float(h.sigma_q_B), "m_r": float(h.mag)}
        r["cosi"] = math.sqrt(max((r["q"] ** 2 - Q0 ** 2) / (1 - Q0 ** 2), 0.0))
        r["i_deg"] = math.degrees(math.acos(r["cosi"]))
        r["sersic_n"] = num(prod.loc[f, "n"]) if f in prod.index else float("nan")
        if h.leg in ("n1", "n1_sky"):
            r["sersic_n"] = 1.0
        for q in ("z", "dm_obs", "dm_host", "tau_scatt_ms", "rm_rad_m2", "offset_kpc", "offset_normalized",
                  "re_kpc", "path_posterior", "tm_gyr", "metallicity_logz", "u_r_color", "snr",
                  "fluence_jy_ms", "width_ms", "sfr", "logmstar", "telescope_localization"):
            v, src = picked(L, f, q)
            if q == "logmstar" and v is None:
                v, src = any_value(L, f, q)
            r[q] = v if q == "telescope_localization" else (num(v) if v is not None else float("nan"))
            r[q + "_src"] = src
        for k, names in MERGED.items():
            r[k], r[k + "_src"] = merged_value(L, f, names)
        for k in ("lin_pol", "circ_pol"):
            if r[k] > 1.5:  # Sherman+2024 quote percent; others fractions
                r[k] /= 100.0
        m = master.loc[f]
        r["ra_deg"], r["dec_deg"] = float(m.ra_deg), float(m.dec_deg)
        r["gal_abs_b"] = abs(float(SkyCoord(r["ra_deg"] * u.deg, r["dec_deg"] * u.deg).galactic.b.deg))
        r["dm_mw"] = num(m.DM_MW)
        if not np.isfinite(r["dm_obs"]):
            r["dm_obs"] = num(m.DM)
        r["survey"] = survey_group(r["telescope_localization"] or SURVEY_FILL.get(f, m.survey))
        z = r["z"]
        ok = np.isfinite(z) and np.isfinite(r["dm_obs"]) and np.isfinite(r["dm_mw"])
        r["dm_exgal"] = r["dm_obs"] - r["dm_mw"] if np.isfinite(r["dm_obs"]) and np.isfinite(r["dm_mw"]) else float("nan")
        dmc = B.median_dm_cosmic(z) if ok else float("nan")
        r["dm_host_obs"] = r["dm_obs"] - r["dm_mw"] - dmc if ok else float("nan")
        r["dm_host_rest"] = r["dm_host_obs"] * (1 + z) if ok else float("nan")
        if not np.isfinite(r["offset_normalized"]) and np.isfinite(r["offset_kpc"]) and np.isfinite(r["re_kpc"]):
            r["offset_normalized"] = r["offset_kpc"] / r["re_kpc"]
        r["abs_rm"] = abs(r["rm_rad_m2"]) if np.isfinite(r["rm_rad_m2"]) else float("nan")
        r["log_tau"] = math.log10(r["tau_scatt_ms"]) if r["tau_scatt_ms"] > 0 else float("nan")
        r["log_abs_rm"] = math.log10(r["abs_rm"]) if r["abs_rm"] > 0 else float("nan")
        r["log_ssfr"] = (math.log10(r["sfr"]) - r["logmstar"]) if r["sfr"] > 0 and np.isfinite(r["logmstar"]) else float("nan")
        r["log_sigma_star"] = (r["logmstar"] - math.log10(2 * math.pi * r["re_kpc"] ** 2)
                               if np.isfinite(r["logmstar"]) and r["re_kpc"] > 0 else float("nan"))
        r["abs_rm_host_rest"] = abs(r["rm_host_rest"]) if np.isfinite(r["rm_host_rest"]) else float("nan")
        r["b_par_proxy"] = (1.232 * r["abs_rm"] * (1 + z) ** 2 / r["dm_host_rest"]
                            if np.isfinite(r["abs_rm"]) and r["dm_host_rest"] > 0 else float("nan"))
        rows.append(r)
    return pd.DataFrame(rows)


def mc_cosi(T: pd.DataFrame, rng: np.random.Generator) -> np.ndarray:
    q = T.q.to_numpy(float)
    s = np.nan_to_num(T.sigma_q_B.to_numpy(float), nan=0.0)
    out = np.empty((N_DRAWS, len(q)))
    for k in range(N_DRAWS):
        qq = np.clip(rng.normal(q, s), Q0 + 1e-6, 1.0)
        inc = np.clip(_hubble_inc_deg(qq, q0=Q0) + rng.normal(0, FLOOR_DEG, len(q)), 0.0, 90.0)
        out[k] = np.cos(np.radians(inc))
    return out


def corr(x: np.ndarray, y: np.ndarray) -> dict:
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 6:
        return {"N": int(len(x))}
    rs = stats.spearmanr(x, y)
    kt = stats.kendalltau(x, y)
    pr = stats.pearsonr(x, y)
    return {"N": int(len(x)), "spearman": rs.statistic, "spearman_p": rs.pvalue,
            "kendall": kt.statistic, "kendall_p": kt.pvalue, "pearson": pr.statistic, "pearson_p": pr.pvalue}


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[float, float, int]:
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    if ok.sum() < 8:
        return float("nan"), float("nan"), int(ok.sum())
    rx, ry, rz = (stats.rankdata(a[ok]) for a in (x, y, z))
    ex = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
    ey = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    r, _ = stats.pearsonr(ex, ey)
    n = int(ok.sum())
    t = r * math.sqrt((n - 3) / max(1 - r * r, 1e-12))
    return float(r), float(2 * stats.t.sf(abs(t), n - 3)), n


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, float)
    q = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return q
    pv = p[ok]
    order = np.argsort(pv)
    ranked = pv[order] * len(pv) / (np.arange(len(pv)) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(pv)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


INCL_PROPS = [
    ("dm_host_rest", "DM_host rest (Macquart median, NE2001)", "Marnoch+26; Bhardwaj+24; Li+25"),
    ("dm_host_obs", "DM_host observer frame", "Marnoch+26"),
    ("dm_host", "DM_host as published", "Marnoch+26; Verdi+26"),
    ("dm_exgal", "DM - DM_MW", "Li+25"),
    ("dm_obs", "DM_obs", "Bhardwaj+24"),
    ("log_tau", "log tau (scattering)", "Glowacki+25; Ocker+22"),
    ("log_abs_rm", "log |RM|", "Verdi+26; Mannings+23"),
    ("offset_kpc", "projected offset [kpc]", "Marnoch+26; Mannings+21"),
    ("offset_normalized", "offset / R_e", "Marnoch+26; Gordon+25"),
    ("logmstar", "log M*", "Leung+25; Li+25"),
    ("sfr", "SFR", "Li+25"),
    ("log_ssfr", "log sSFR", "Li+25"),
    ("log_sigma_star", "log Sigma* (M*/2 pi R_e^2)", "Glowacki+25"),
    ("tm_gyr", "mass-weighted age", "Glowacki+25"),
    ("metallicity_logz", "stellar metallicity", "Glowacki+25"),
    ("u_r_color", "u - r colour", "Marnoch+26 (colour)"),
    ("re_kpc", "R_e [kpc]", "Li+25 (projected area)"),
    ("av", "A_V young stars (birth clouds + diffuse; Prospector)", "dust"),
    ("av_old", "A_V old stars (diffuse ISM; Prospector)", "positive control: diffuse disc dust rises edge-on"),
    ("g_r_color", "g - r colour", "positive control (reddening)"),
    ("abs_rm_host_rest", "|RM_host| rest (published)", "Verdi+26; Mannings+23"),
    ("b_par_host", "B_par host (published)", "Verdi+26"),
    ("b_par_proxy", "1.23 |RM| (1+z)^2 / DM_host,rest", "Verdi+26 style, RM_MW not removed"),
    ("lin_pol", "linear polarisation fraction", "depolarisation through the disc"),
    ("circ_pol", "circular polarisation fraction", "propagation"),
    ("log_tau_rest", "log tau rest 1 GHz (published)", "Glowacki+25"),
    ("metal_gas", "gas metallicity 12+log O/H", "Glowacki+25"),
    ("halpha_ew", "H-alpha EW", "Glowacki+25"),
    ("q_pah", "PAH fraction q_PAH", "dust"),
    ("scint_bw", "scintillation bandwidth", "null check (Milky Way screen)"),
    ("z", "redshift", "selection"),
    ("m_r", "m_r (ours)", "selection; Marnoch m<21 cut"),
    ("sersic_n", "Sersic n (ours)", "selection; Marnoch n<2 cut"),
    ("path_posterior", "PATH P(O|x)", "association"),
    ("snr", "S/N", "selection"),
    ("fluence_jy_ms", "fluence", "selection"),
    ("width_ms", "width", "selection; Glowacki+25"),
    ("dm_mw", "DM_MW (NE2001)", "null check"),
    ("gal_abs_b", "|b| Galactic", "null check"),
    ("dec_deg", "Dec", "null check (survey footprint)"),
    ("ra_deg", "RA", "null check"),
]

LIT_PAIRS = [
    ("logmstar", "dm_host_rest", "Leung+25: DM_host anticorrelates with M*"),
    ("log_ssfr", "dm_exgal", "Li+25: DM_exc tightly correlates with sSFR"),
    ("sfr", "dm_exgal", "Li+25: DM_exc vs SFR"),
    ("log_sigma_star", "log_tau", "Glowacki+25: tau correlates with Sigma*"),
    ("tm_gyr", "log_tau", "Glowacki+25: tau correlates with mass-weighted age"),
    ("metallicity_logz", "log_tau", "Glowacki+25: weak tau-metallicity"),
    ("dm_host_rest", "log_tau", "Cordes+22 / Day+20: tau vs DM"),
    ("dm_host_rest", "log_abs_rm", "Verdi+26: DM_host vs RM_host"),
    ("offset_normalized", "dm_host_rest", "Marnoch+26 / Mannings+21: DM_host vs offset"),
    ("z", "dm_host_rest", "Verdi+26: no z evolution of DM_host"),
]


def fmt(x, nd=3):
    if x is None or not np.isfinite(x):
        return "—"
    return f"{x:.2e}" if abs(x) < 1e-3 and x != 0 else f"{x:.{nd}f}"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    T = build_table()
    T.to_csv(OUT / "host_props.csv", index=False)
    rng = np.random.default_rng(SEED)
    C = mc_cosi(T, rng)
    cosi = T.cosi.to_numpy(float)

    rows = []
    for col, label, ref in INCL_PROPS:
        x = T[col].to_numpy(float)
        r = {"property": col, "label": label, "literature": ref}
        r.update(corr(cosi, x))
        ok = np.isfinite(x)
        if ok.sum() >= 6:
            rho = np.array([stats.spearmanr(c[ok], x[ok]).statistic for c in C])
            pp = np.array([stats.spearmanr(c[ok], x[ok]).pvalue for c in C[:500]])
            r["mc_rho_16"], r["mc_rho_50"], r["mc_rho_84"] = np.percentile(rho, [16, 50, 84])
            r["mc_frac_p05"] = float(np.mean(pp < 0.05))
            edge = T.i_deg.to_numpy(float) > 60
            a, b = x[ok & edge], x[ok & ~edge]
            if len(a) >= 3 and len(b) >= 3:
                r["N_edge"], r["N_face"] = len(a), len(b)
                r["med_edge"], r["med_face"] = float(np.median(a)), float(np.median(b))
                r["mwu_edge_face_p"] = float(stats.mannwhitneyu(a, b, alternative="two-sided").pvalue)
        rows.append(r)
    R = pd.DataFrame(rows)
    R["spearman_q_BH"] = bh(R.get("spearman_p", pd.Series(dtype=float)).to_numpy(float))
    R.to_csv(OUT / "correlations.csv", index=False)

    # Marnoch-style DM_host block
    mrows = []
    masks = {
        "all": np.ones(len(T), bool),
        "m_r < 21": T.m_r.to_numpy(float) < 21,
        "n < 2": T.sersic_n.to_numpy(float) < 2,
        "m_r < 21 & n < 2": (T.m_r.to_numpy(float) < 21) & (T.sersic_n.to_numpy(float) < 2),
        "spectroscopic z <= 0.3": T.z.to_numpy(float) <= 0.3,
    }
    for dmcol in ("dm_host_rest", "dm_host_obs", "dm_host"):
        for xname, xv in (("q", T.q.to_numpy(float)), ("i", T.i_deg.to_numpy(float))):
            for mname, mk in masks.items():
                d = corr(np.where(mk, xv, np.nan), T[dmcol].to_numpy(float))
                pr, pp, pn = partial_spearman(np.where(mk, xv, np.nan), T[dmcol].to_numpy(float), T.z.to_numpy(float))
                d.update({"dm": dmcol, "x": xname, "subsample": mname, "partial_rho_z": pr, "partial_p_z": pp})
                mrows.append(d)
    M = pd.DataFrame(mrows)
    M.to_csv(OUT / "marnoch_dmhost.csv", index=False)

    lrows = []
    for a, b, claim in LIT_PAIRS:
        d = corr(T[a].to_numpy(float), T[b].to_numpy(float))
        d.update({"x": a, "y": b, "claim": claim})
        lrows.append(d)
    Lr = pd.DataFrame(lrows)
    Lr.to_csv(OUT / "literature_relations.csv", index=False)

    groups = {g: T.loc[T.survey == g, "cosi"].to_numpy(float) for g in sorted(T.survey.unique())}
    big = {g: v for g, v in groups.items() if len(v) >= 5}
    kw = stats.kruskal(*big.values()) if len(big) >= 2 else None

    # figure
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4))
    for ax, (col, lab) in zip(axes, (("dm_host_rest", r"DM$_{\rm host}$ rest [pc cm$^{-3}$]"),
                                     ("log_tau", r"$\log\,\tau$ [ms]"), ("log_abs_rm", r"$\log|{\rm RM}|$"))):
        ok = np.isfinite(T[col])
        ax.scatter(T.cosi[ok], T[col][ok], s=16, c="k")
        rr = R.set_index("property").loc[col]
        ax.set_title(fr"$r_S$={rr.spearman:.2f}, p={rr.spearman_p:.2f}, N={int(rr.N)}", fontsize=9)
        ax.set_xlabel(r"$\cos i$")
        ax.set_ylabel(lab)
    fig.tight_layout()
    fig.savefig(OUT / "dmhost_tau_rm_vs_cosi.png", dpi=150)

    n = len(T)
    md = [f"# Host and burst properties vs inclination (confirmed sample, N = {n})", "",
          "Source: `scripts/run_host_property_correlations.py`. Inclination is the winning-leg $q$ through "
          "$\\cos^2 i = (q^2 - q_0^2)/(1 - q_0^2)$, $q_0 = 0.2$. Correlations are against $\\cos i$ (Spearman "
          "and Kendall are identical in magnitude, opposite in sign, against $i$). "
          f"MC: {N_DRAWS} Protocol B draws ($q \\sim N(q, \\sigma_{{q,B}})$, $i$ += $N(0, 5^\\circ)$); "
          "`frac p<0.05` uses the first 500 draws. BH $q$ is the Benjamini–Hochberg false-discovery value across "
          "the rows of Table 1. Edge-on = $i > 60^\\circ$.", "",
          "DM$_{\\rm host}$ = (DM − DM$_{\\rm MW,NE2001}$ − median Macquart DM$_{\\rm cosmic}(z)$) [× (1+z) for rest frame], "
          "no halo term, spectroscopic $z$ only. Literature values are picked with the priority rules of "
          "`build_host_property_table.py`; masses fall back to K-band / compiled values where no SED mass exists.", "",
          "## 1. Each property vs $\\cos i$", "",
          "| property | literature | N | Spearman $r_S$ (p) | Kendall $\\tau$ (p) | Pearson $r$ (p) | MC $r_S$ 16/50/84 | "
          "MC frac p<0.05 | BH q | median edge / face (MWU p) |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in R.iterrows():
        if not np.isfinite(r.get("spearman", np.nan)):
            md.append(f"| {r.label} | {r.literature} | {int(r.N)} | — | — | — | — | — | — | — |")
            continue
        ef = (f"{fmt(r.med_edge, 3)} / {fmt(r.med_face, 3)} ({fmt(r.mwu_edge_face_p)}; {int(r.N_edge)}/{int(r.N_face)})"
              if np.isfinite(r.get("mwu_edge_face_p", np.nan)) else "—")
        md.append(f"| {r.label} | {r.literature} | {int(r.N)} | {fmt(r.spearman)} ({fmt(r.spearman_p)}) | "
                  f"{fmt(r.kendall)} ({fmt(r.kendall_p)}) | {fmt(r.pearson)} ({fmt(r.pearson_p)}) | "
                  f"{fmt(r.mc_rho_16, 2)} / {fmt(r.mc_rho_50, 2)} / {fmt(r.mc_rho_84, 2)} | {fmt(r.mc_frac_p05, 2)} | "
                  f"{fmt(r.spearman_q_BH)} | {ef} |")
    md += ["", "## 2. Marnoch et al. (2026) DM$_{\\rm host}$ test", "",
           "Marnoch+2026 (CRAFT, N = 34): $r_S \\approx -0.44$, $p \\approx 0.01$ between rest-frame DM$_{\\rm host}$ "
           "and $q$ (more DM in edge-on hosts). Partial $r_S$ controls for $z$.", "",
           "| DM_host | vs | subsample | N | Spearman (p) | Pearson (p) | partial Spearman given z (p) |",
           "|---|---|---|---|---|---|---|"]
    for _, r in M.iterrows():
        md.append(f"| {r.dm} | {r.x} | {r.subsample} | {int(r.N)} | {fmt(r.get('spearman'))} ({fmt(r.get('spearman_p'))}) | "
                  f"{fmt(r.get('pearson'))} ({fmt(r.get('pearson_p'))}) | {fmt(r.partial_rho_z)} ({fmt(r.partial_p_z)}) |")
    md += ["", "## 3. Literature host-property relations on this sample (no inclination)", "",
           "| x | y | claim | N | Spearman (p) | Pearson (p) |", "|---|---|---|---|---|---|"]
    for _, r in Lr.iterrows():
        md.append(f"| {r.x} | {r.y} | {r.claim} | {int(r.N)} | {fmt(r.get('spearman'))} ({fmt(r.get('spearman_p'))}) | "
                  f"{fmt(r.get('pearson'))} ({fmt(r.get('pearson_p'))}) |")
    md += ["", "## 4. $\\cos i$ by survey", "",
           "| survey | N | median $\\cos i$ |", "|---|---|---|"]
    for g, v in groups.items():
        md.append(f"| {g} | {len(v)} | {np.median(v):.3f} |")
    if kw is not None:
        md += ["", f"Kruskal–Wallis over surveys with N ≥ 5 ({', '.join(big)}): H = {kw.statistic:.2f}, p = {kw.pvalue:.3f}."]
    miss = T.loc[~np.isfinite(T.dm_host_rest), "frb"].tolist()
    md += ["", f"Hosts without DM$_{{\\rm host}}$: {', '.join(miss)}.", ""]
    (OUT / "CORRELATIONS.md").write_text("\n".join(md), encoding="utf-8")
    sys.stdout.reconfigure(encoding="utf-8")
    print("\n".join(md))


if __name__ == "__main__":
    main()
