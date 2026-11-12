"""Checks of foreground_halo_dm.py: agreement with the original notebook code (Estimate_DM_halo.ipynb, cell 26)
and science consistency.

Output: Publication/materials/foreground_halos/CODE_CHECKS.md
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import constants as const
from astropy import units as u
from astropy.coordinates import SkyCoord
from colossus.cosmology import cosmology
from colossus.halo import concentration
from scipy import integrate, optimize, stats
from scipy.special import hyp2f1

sys.path.insert(0, str(Path(__file__).resolve().parent))
import foreground_halo_dm as F  # noqa: E402
import robust_correlations as RC  # noqa: E402

cosmology.setCosmology("planck18")
OUT = F.OUT / "CODE_CHECKS.md"


# ---------------------------------------------------------------- notebook reference (cell 26, deterministic)
def nb_logmh(logms, z):
    return optimize.fsolve(lambda x: F.moster13_logms(x, z) - logms, logms + 2.0)[0]


def nb_concentration(m_msun, z, h_units=False):
    """Colossus expects M in Msun/h; the notebook passes Msun (h_units=False)."""
    return concentration.concentration(m_msun * (F.cosmo.h if h_units else 1.0), "200c", z, model="dutton14")


def nb_dm(logmh, z, theta_arcsec, num=1000, h_units=False):
    """Notebook ModifiedNFW.compute_DM / (1 + z), with c from Colossus as in the notebook."""
    m = 10 ** logmh * u.Msun
    c = nb_concentration(m.value, z, h_units)
    rhoc = F.cosmo.critical_density(z)
    r200 = ((3 * m / (4 * np.pi * 200 * rhoc)) ** (1 / 3)).to(u.kpc)
    rho0 = (200 * rhoc / 3) * c ** 3 / (np.log(1 + c) - c / (1 + c))
    d = theta_arcsec * u.arcsec * F.cosmo.kpc_proper_per_arcmin(z).to(u.kpc / u.arcsec)
    if d >= r200:
        return 0.0
    chord = 2 * np.sqrt((r200 ** 2 - d ** 2).value)
    x = np.linspace(-chord / 2, chord / 2, num) * u.kpc
    y = c * np.sqrt(x.value ** 2 + d.value ** 2) / r200.value
    rho = 0.75 * F.FB * rho0 / (y ** (1 - 2) * (2 + y) ** 4)
    ne = (rho / (1.3 * const.m_p)).to(u.cm ** -3).value * 1.1667
    return float(np.trapezoid(ne, x.to(u.pc).value) / (1 + z))


def notebook_agreement():
    rows = []
    for logms in (8.0, 9.0, 10.0, 11.0, 11.5):
        for z in (0.02, 0.1, 0.3, 0.6):
            lmh_nb, lmh_ours = nb_logmh(logms, z), F.logmh_from_logms(logms, z)[0]
            m = 10 ** lmh_ours
            r200_arcsec = F.r200_kpc(m, z) / F.kpc_per_arcsec(z)
            for frac in (0.02, 0.2, 0.5, 0.9):
                theta = frac * r200_arcsec
                rows.append(dict(logms=logms, z=z, frac_r200=frac, dlogmh=lmh_ours - lmh_nb,
                                 c_ratio=F.dutton14_c200c(m, z) / nb_concentration(m, z),
                                 c_ratio_h=F.dutton14_c200c(m, z) / nb_concentration(m, z, True),
                                 dm_ours=F.dm_mnfw(lmh_ours, z, theta)[0], dm_nb=nb_dm(lmh_nb, z, theta),
                                 dm_nb_h=nb_dm(lmh_nb, z, theta, h_units=True)))
    d = pd.DataFrame(rows)
    d["rel"] = d.dm_ours / d.dm_nb - 1
    d["rel_h"] = d.dm_ours / d.dm_nb_h - 1
    return d


def cosmic_budget():
    """Foreground-halo DM along the random control sightlines as a fraction of the Macquart mean DM_cosmic."""
    import build_host_property_table as BT
    c = pd.read_csv(F.OUT / "foreground_dm_controls.csv")
    c = c[c.in_ls.astype(str).isin(["dr9", "dr10"])]
    frac = c.dm_fg / c.z.map(BT.mean_dm_cosmic)
    return len(c), frac.median(), frac.mean()


# ---------------------------------------------------------------- science checks
def fy_dm(c):
    return np.log(1 + c) - c / (1 + c)


def fy_b(c, alpha=2.0, y0=2.0):
    """Enclosed baryon mass function of the modified NFW (frb package, Prochaska & Zheng 2019)."""
    return (c / (y0 + c)) ** (1 + alpha) * (y0 ** -alpha * (y0 + c) ** (1 + alpha)
                                            * hyp2f1(1 + alpha, 1 + alpha, 2 + alpha, -c / y0) - y0) / (1 + alpha) / y0


def mass_norm_factor(m_msun, z):
    """DM ratio (gas-mass-normalized / dark-matter-normalized) for the same profile shape."""
    c = F.dutton14_c200c(m_msun, z)
    return fy_dm(c) / fy_b(c)


def mw_benchmark():
    """Radial DM from the Sun (8 kpc) to r200 for a 1.5e12 Msun halo, c = 7.67 (Prochaska & Zheng 2019)."""
    m, c, z = 1.5e12, 7.67, 0.0
    r200 = F.r200_kpc(m, z)
    rho0 = (200 * F.rho_crit(z) / 3) * c ** 3 / fy_dm(c)
    r = np.linspace(8.0, r200, 20000)
    y = c * r / r200
    ne = 0.75 * F.FB * rho0 / (y ** -1 * (2 + y) ** 4) / (1.3 * F.M_P) * 1.1667
    dm = np.trapezoid(ne, r * 1e3)
    return dm, dm * fy_dm(c) / fy_b(c), r200


def integration_accuracy():
    out = []
    for logmh, z, frac in ((11.0, 0.05, 0.05), (12.0, 0.2, 0.3), (13.0, 0.4, 0.8)):
        m = 10 ** logmh
        c, r200 = F.dutton14_c200c(m, z), F.r200_kpc(m, z)
        rho0 = (200 * F.rho_crit(z) / 3) * c ** 3 / fy_dm(c)
        b = frac * r200

        def ne(s):
            y = c * np.sqrt(s ** 2 + b ** 2) / r200
            return 0.75 * F.FB * rho0 / (y ** -1 * (2 + y) ** 4) / (1.3 * F.M_P) * 1.1667

        half = np.sqrt(r200 ** 2 - b ** 2)
        exact = 2 * integrate.quad(ne, 0, half, limit=200)[0] * 1e3 / (1 + z)
        ours = F.dm_mnfw(logmh, z, b / F.kpc_per_arcsec(z))[0]
        out.append(ours / exact - 1)
    return np.array(out)


def bell_error_coefficient():
    """Median sigma(log M*) with the coded r-band coefficient (0.659) and the analytic one (1.497)."""
    sig = []
    for f in F.CACHE.glob("*_ls.csv"):
        c = pd.read_csv(f)[["flux_g", "flux_ivar_g", "flux_r", "flux_ivar_r"]].apply(pd.to_numeric, errors="coerce")
        if not len(c):
            continue
        ge, re_ = F.mag_err(c.flux_g.to_numpy(), c.flux_ivar_g.to_numpy()), F.mag_err(c.flux_r.to_numpy(), c.flux_ivar_r.to_numpy())
        sig.append(np.c_[np.hypot(np.hypot(1.097 * ge, 0.659 * re_), 0.1), np.hypot(np.hypot(1.097 * ge, 1.497 * re_), 0.1)])
    s = np.vstack(sig)
    return np.median(s, 0), np.percentile(s, 90, axis=0)


def z_ref(gal, hosts):
    """Redshift the code uses for each galaxy: spec-z (already capped at z_host), else min(z_phot, z_host)."""
    z_host = gal.frb.map(hosts.set_index("frb").z).to_numpy()
    return np.where(np.isfinite(gal.z_spec), gal.z_spec, np.minimum(gal.z_phot, z_host))


def kcorrection_effect(gal, hosts):
    d = gal[gal.type != "SGA"].copy()
    z = z_ref(d, hosts)
    no_k = 1.097 * (d.g - d.r) - 0.406 - 0.4 * (d.r - F.cosmo.distmod(z).value - 4.65)
    return np.nanmedian(d.logms_at_zref - no_k), np.nanpercentile(d.logms_at_zref - no_k, [16, 84])


def local_group(hosts):
    """FRBs inside the projected r200 of M31 and M33 (excluded by the SGA z_leda > 0.0005 cut)."""
    frb = SkyCoord(hosts.ra_deg.to_numpy() * u.deg, hosts.dec_deg.to_numpy() * u.deg)
    out = []
    for name, ra, dec, d_kpc, logmh in (("M31", 10.6847, 41.2690, 765, 12.1), ("M33", 23.4621, 30.6599, 840, 11.4)):
        r200 = F.r200_kpc(10 ** logmh, 0.0)
        radius = np.degrees(r200 / d_kpc)
        sep = frb.separation(SkyCoord(ra * u.deg, dec * u.deg)).deg
        out.append((name, radius, int((sep < radius).sum()), float(sep.min())))
    return out


def sga_radius_check():
    """Largest angular r200 among SGA galaxies (does the 30 deg search radius bind?)."""
    best = 0.0
    for f in F.CACHE.glob("*_sga.csv"):
        s = pd.read_csv(f)
        if not len(s):
            continue
        s = s.apply(pd.to_numeric, errors="coerce")
        g, r = F.ab_mag(s.flux_g, s.mw_transmission_g), F.ab_mag(s.flux_r, s.mw_transmission_r)
        ok = (np.isfinite(g) & np.isfinite(r)).to_numpy()
        if not ok.any():
            continue
        g, r, z = g[ok].to_numpy(), r[ok].to_numpy(), s.z_leda[ok].to_numpy()
        logms = 1.097 * (g - r) - 0.406 - 0.4 * (r - F.cosmo.distmod(z).value - 4.65)
        z_bin = np.round(z, 2)
        lmh = np.empty_like(logms)
        for zb in np.unique(z_bin):
            k = z_bin == zb
            lmh[k] = np.interp(logms[k], F.moster13_logms(F.LOGMH_GRID, max(zb, 1e-3)), F.LOGMH_GRID)
        best = max(best, np.max(F.r200_kpc(10 ** lmh, z) / F.kpc_per_arcsec(z) / 3600))
    return best


def renormalized_results(gal, hosts_fg):
    """Re-derive the foreground-corrected DM_host with gas-mass-normalized halos and re-test the key correlations."""
    g = gal.copy()
    z = z_ref(g, hosts_fg)
    logms = np.array([m if np.isfinite(m) else float(F.logms_bell03_gr(gg, rr, zz))
                      for m, gg, rr, zz in zip(g.logms_at_zref, g.g, g.r, z)])
    g["factor"] = [mass_norm_factor(10 ** F.logmh_from_logms(m, zz)[0], zz) for m, zz in zip(logms, z)]
    g["dm_w"] = g.dm_mean * g.p_fg
    per = g.groupby("frb").apply(lambda t: pd.Series(dict(dm_fg=t.dm_w.sum(), dm_fg_mass=(t.dm_w * t.factor).sum())),
                                 include_groups=False).reset_index()
    h = RC.load().merge(per, on="frb", how="left", suffixes=("", "_sum"))
    ls = h.in_ls.astype(str).isin(["dr9", "dr10"])
    k = np.nanmedian(g.factor)
    h["dm_host_fg_mass"] = np.where(ls, h.dm_host_rest - h.dm_fg_mass.fillna(0) * (1 + h.z), np.nan)
    h["dm_host_fg_excess_mass"] = np.where(ls & h.dm_fg_expected.notna(), h.dm_host_rest
                                           - (h.dm_fg_mass.fillna(0) - k * h.dm_fg_expected) * (1 + h.z), np.nan)
    tests = []
    for col in ("dm_host_rest", "dm_host_fg", "dm_host_fg_mass", "dm_host_fg_excess", "dm_host_fg_excess_mass"):
        for other, lab in (("q", "q"), ("log_abs_rm", r"log \|RM\|"), ("log_tau", "log tau")):
            d = h[np.isfinite(h[col]) & np.isfinite(h[other])]
            r, p = stats.spearmanr(d[other], d[col])
            tests.append(dict(dm=col, x=lab, N=len(d), r=r, p=p))
    return per, k, pd.DataFrame(tests), g.factor.describe()


def main():
    gal = pd.read_csv(F.OUT / "foreground_galaxies.csv", dtype={"frb": str})
    hosts = pd.read_csv(F.OUT / "foreground_dm_per_host.csv", dtype={"frb": str})
    host_props = pd.read_csv(F.HOSTS, dtype={"frb": str})

    nb = notebook_agreement()
    val = F.validate()
    nonzero = nb[nb.dm_nb > 0]
    mw_dm, mw_dm_mass, mw_r200 = mw_benchmark()
    integ = integration_accuracy()
    print("notebook and integration checks done", flush=True)
    (s_med, s_med_fix), (s_90, s_90_fix) = bell_error_coefficient()
    dk, dk_range = kcorrection_effect(gal, hosts)
    n_ctrl, budget_med, budget_mean = cosmic_budget()
    lg = local_group(host_props)
    per, k, tests, fac = renormalized_results(gal, hosts)
    print("renormalization done", flush=True)
    sga_max = sga_radius_check()
    fyb_numeric = integrate.quad(lambda y: y ** 3 / (2 + y) ** 4, 0, 7.67)[0]
    ls_hosts = hosts[hosts.in_ls.astype(str).isin(["dr9", "dr10"])]
    cond_ratio = (ls_hosts.dm_fg_conditional / ls_hosts.dm_fg).replace([np.inf], np.nan)
    cs = np.array([5.0, 7.0, 10.0, 15.0])

    lines = [
        "# Checks of the foreground-halo DM code", "",
        "Generated by `scripts/check_foreground_halo_dm.py`.", "",
        "## 1. Agreement with the original code", "",
        f"- **Notebook validation (FRB 20240114A, Estimate_DM_halo.ipynb cell 26):** notebook medians "
        f"{F.NOTEBOOK_VALIDATION}; this code {[float(v) for v in val]} (about 5% higher; explained below).",
        f"- **Deterministic grid** ({len(nb)} cases: log M* 8-11.5, z 0.02-0.6, impact 0.02-0.9 r200) against a "
        f"re-implementation of the notebook's `ModifiedNFW.compute_DM` with Colossus concentrations and `fsolve` "
        f"halo masses:",
        f"  - exactly as in the notebook: this code is higher by median {np.median(nonzero.rel) * 100:.1f}% "
        f"(max {nonzero.rel.max() * 100:.1f}%);",
        f"  - with the mass passed to Colossus in Msun/h, as Colossus expects: median "
        f"{np.median(np.abs(nonzero.rel_h)) * 100:.3f}% (max {np.abs(nonzero.rel_h).max() * 100:.3f}%).",
        f"  - So the whole difference is the notebook passing M in Msun to Colossus, which lowers c by "
        f"{(1 - 1 / nb.c_ratio.min()) * 100:.1f}-{(1 - 1 / nb.c_ratio.max()) * 100:.1f}%. This code's analytic "
        f"Dutton & Maccio (2014) c matches Colossus with Msun/h to a ratio of "
        f"{nb.c_ratio_h.min():.4f}-{nb.c_ratio_h.max():.4f}.",
        f"  - Halo-mass inversion (grid vs `fsolve`): max |difference| {np.abs(nb.dlogmh).max():.4f} dex.",
        f"- **Integration accuracy** (801-point trapezoid vs adaptive quadrature): "
        f"max |difference| {np.abs(integ).max() * 100:.3f}%.",
        "- **Refactor:** the tidied script reproduced `foreground_dm_per_host.csv`, `foreground_galaxies.csv` and "
        "`validation_20240114A.csv` exactly (commit d788368).", "",
        "Deliberate differences from the notebooks (all as in the Bhardwaj+2025 text unless noted):", "",
        "| Step | Notebooks / paper | This code | Why |", "|---|---|---|---|",
        "| Stellar mass | Taylor+11 g-i (SDSS) | Bell+03 g-r, Chabrier (-0.406), as in atel_analysis.ipynb | LS DR9 has no i band |",
        f"| Colours | observed | K-corrected to rest frame | median shift {dk:+.2f} dex in log M* (16-84%: {dk_range[0]:+.2f} to {dk_range[1]:+.2f}) |",
        "| Photo-z galaxies | truncated Gaussian, mean DM (in front assumed) | same (`dm_fg_conditional`); main column also weights by P(in front) | conditional sum overcounts galaxies that are probably behind |",
        f"| | | | conditional / weighted total: median {np.nanmedian(cond_ratio):.2f} |",
        "| Local-volume check | HECATE within 30 deg | SGA-2020 within 30 deg | SGA is in the same Data Lab service |",
        "| Halo mass | `fsolve` per draw | interpolation on a 0.005 dex grid | speed |", "",
        "## 2. Science consistency", "",
        "Correct as implemented:", "",
        "- Moster+13 parameters (N, M1, beta, gamma and their z/(1+z) evolution), Dutton & Maccio (2014) c200c with "
        "M in h^-1 Msun, r200 at 200 rho_c(z), proper impact parameter, n_e/n_H = 1.1667 (X = 0.75, fully ionized "
        "H and He), path in pc so DM is in pc cm^-3, and the observer-frame factor 1/(1 + z_gal).",
        "- K-correction sign: for AB magnitudes M_rest = m(lambda_rest (1+z)) - DM_dist + 2.5 log(1+z).",
        "- The 0.018 size ratio is applied as R_eff / theta > 0.018 (the paper text writes it inverted; the "
        "notebook and this code use R_eff / theta).",
        "- Subtraction: foreground DM is observer-frame, removed from DM_obs before the (1 + z_host) factor.",
        "- Foreground halos are also part of the Macquart mean DM_cosmic; the 'excess over random sightlines' "
        "column removes only the excess.", "",
        "Issues found:", "",
        f"1. **Gas normalization.** Eq. B3 of Bhardwaj+2025 (and the notebook, and this code) normalizes the "
        f"modified profile with the dark-matter NFW rho_0. The gas mass inside r200 is then f_gas f_b M200 x "
        f"f_b(c)/f_dm(c), only {', '.join(f'{fy_b(c) / fy_dm(c):.2f}' for c in cs)} of the intended value for "
        f"c = {', '.join(f'{c:g}' for c in cs)}. The `frb` package (Prochaska & Zheng 2019's own code) normalizes "
        f"by gas mass. Milky Way benchmark (1.5e12 Msun, c = 7.67, Sun to r200 = {mw_r200:.0f} kpc): "
        f"{mw_dm:.0f} pc cm^-3 as coded vs {mw_dm_mass:.0f} gas-mass-normalized; Prochaska & Zheng (2019) quote "
        f"50-80. The coded halo DMs are a factor {fac['25%']:.1f}-{fac['75%']:.1f} below the gas-mass-normalized "
        f"profile (interquartile over the foreground galaxies, median {fac['50%']:.1f}); equivalently, the code "
        f"behaves like f_hot = 0.75 x f_b(c)/f_dm(c), i.e. "
        f"{0.75 * fy_b(cs[0]) / fy_dm(cs[0]):.2f}-{0.75 * fy_b(cs[-1]) / fy_dm(cs[-1]):.2f}. (Check of the analytic f_b: "
        f"{fy_b(7.67):.4f} vs numerical {fyb_numeric:.4f}.)",
        f"   - **Cosmic budget test.** Along the {n_ctrl} random control sightlines, the coded foreground halos "
        f"supply a mean {budget_mean * 100:.0f}% (median {budget_med * 100:.0f}%) of the Macquart mean DM_cosmic. "
        f"Gas-mass normalization with f_hot = 0.75 would raise this to about {budget_mean * k * 100:.0f}%, from "
        f"LS-detected galaxies alone, leaving too little for the intergalactic medium. f_hot = 0.75 for every halo "
        f"is therefore an upper bound, and the coded convention is a reasonable, if low-side, choice.",
        f"2. **r-band error coefficient.** sigma(log M*) uses 0.659 for the r term; the derivative of the Bell "
        f"relation is 1.097 + 0.4 = 1.497. Median sigma(log M*) {s_med:.3f} vs {s_med_fix:.3f} dex "
        f"(90th percentile {s_90:.3f} vs {s_90_fix:.3f}); the 0.10 dex systematic dominates.",
        f"3. **Local Group.** SGA `z_leda > 0.0005` excludes M31 and M33. FRBs inside their projected r200: "
        + "; ".join(f"{n} (radius {rad:.0f} deg): {k_} (closest {sep:.0f} deg)" for n, rad, k_, sep in lg) + ".",
        f"4. **30 deg search radius.** Largest angular r200 among cached SGA galaxies: {sga_max:.1f} deg, so 30 deg "
        "does not bind.", "",
        "## 3. Effect of the gas normalization on our results", "",
        f"Each galaxy's DM rescaled by f_dm(c)/f_b(c) at its median mass; random-sightline expectation scaled by the "
        f"median factor ({k:.2f}). Per-host foreground DM (probability-weighted, observer frame): median "
        f"{per.dm_fg.median():.1f} -> {per.dm_fg_mass.median():.1f} pc cm^-3. Rows ending in `_mass` are the "
        f"upper bound from the cosmic budget test; the coded rows are the baseline.", "",
        "| DM_host version | vs | N | r_S | p |", "|---|---|---|---|---|"]
    for t in tests.itertuples():
        lines.append(f"| {t.dm} | {t.x} | {t.N} | {t.r:+.2f} | {t.p:.3f} |")
    OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
