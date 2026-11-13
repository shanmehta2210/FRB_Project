"""Test of the halo-mass inversion behind the foreground-halo DM of Bhardwaj et al. (2025, arXiv:2506.11915).

Everything is computed with the functions of Estimate_DM_halo.ipynb (cell 26) executed unchanged:
stellarmass_from_halomass (Moster+13), halomass_from_stellarmass, concentration_index, ModifiedNFW,
and the Taylor stellar masses and Monte Carlo loop of its main block. The only change is the line that
assigns a halo mass to each stellar-mass draw:

  inverse      logMh = halomass_from_stellarmass(logM*, z)            (the notebook and the paper)
  conditional  logMh drawn from P(Mh | M*) ~ dn/dlnMh * N(logM*; Moster13(Mh), 0.15)

Moster+13 fit the mean stellar mass at fixed halo mass with 0.15 dex scatter. The conditional
distribution is what that model, combined with the halo mass function, predicts for a galaxy of given
stellar mass; the inverse is the halo whose average galaxy has that stellar mass. They agree below
log M* ~ 10.5 and separate above it, where M* ~ Mh^0.3 and massive halos are rare.

Outputs: Publication/materials/foreground_halos/mass_inversion/
  conditional_mass.csv, dm_profiles.csv, frb20240114a_galaxies.csv, frb20240114a_budget.csv,
  scatter_sensitivity.csv, mass_inversion.{pdf,png}
Run: python scripts/mass_inversion_test.py
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from colossus.cosmology import cosmology
from colossus.halo import mass_adv
from colossus.lss import mass_function

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK = ROOT / "Estimate_DM_halo.ipynb"
OUT = ROOT / "Publication/materials/foreground_halos/mass_inversion"
SCATTER = 0.15
N_MC = 5000
rng = np.random.default_rng(11)
cosmology.setCosmology("planck18")
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "stix", "font.serif": ["STIXGeneral"]})

# notebook main-block inputs: host, central galaxy, Table 5 galaxies 4, 6, 7, 8
GAL = pd.DataFrame({
    "name": ["host", "central", "4", "6", "7", "8"],
    "g": [22.361, 20.687, 19.238, 16.487, 17.767, 16.431],
    "g_err": [0.137, 0.052, 0.019, 0.006, 0.009, 0.007],
    "i": [21.7388, 19.4518, 18.1228, 14.6448, 16.4618, 15.7618],
    "i_err": [0.134, 0.040, 0.013, 0.003, 0.006, 0.006],
    "z": [0.13, 0.1296, 0.09, 0.090, 0.09, 0.029],
    "theta": [0.4042, 35.4207, 87.17, 281.2489, 123.5773, 499.2689],
})

# paper, Sec. III.2 and Appendix B
DM_OBS, DM_DISK, DM_MWHALO, DM_IGM, DM_ICM = 527.72, 45.0, 45.0, 116.0, 60.0
SIG_DISK, SIG_MWHALO, SIG_IGM = 9.0, 22.5, 61.0
DM_PHOTZ = {"2": 9.0, "3": 3.0, "5": 3.0}
DM_CENTRAL_PAPER = 11.0
HALPHA = (126.0, 26.0, 81.0)
CLUSTER = dict(name="J212719.9+042225", ra="21h27m19.9s", dec="+04d22m25s", z=0.0903, m500=6.8e13)
GAL6_POS = ("21h27m25.44s", "+04d22m46.5s")


def notebook_code():
    src = "".join(json.load(open(NOTEBOOK, encoding="utf-8"))["cells"][26]["source"])
    if not hasattr(np, "trapz"):
        np.trapz = np.trapezoid  # removed in numpy 2.4; same function
    ns = {"__name__": "notebook"}
    exec(src, ns)
    return ns


NB = notebook_code()


def taylor_logms(row):
    d_l = NB["cosmo"].luminosity_distance(row.z).value
    m_i = row.i - 5 * (np.log10(d_l * 1e6) - 1)
    logm = 1.15 + 0.70 * (row.g - row.i) - 0.4 * m_i
    stat = np.sqrt((0.70 * row.g_err) ** 2 + (1.1 * row.i_err) ** 2)
    return logm, np.sqrt(stat ** 2 + 0.1 ** 2)


def halo_dm(logmh, z, theta):
    """The notebook's per-draw DM for a given halo mass (observer frame)."""
    c = NB["concentration_index"](10 ** logmh, z)
    halo = NB["ModifiedNFW"](logmh, c, z)
    d_proj = theta * u.arcsec * NB["kpc_per_arcsec"](z)
    if d_proj >= halo.r200:
        return 0.0
    return (halo.compute_DM(theta) / (1 + z)).value


LMH = np.linspace(10.0, 16.0, 1201)


def conditional_weights(logms, z, scatter=SCATTER):
    """P(log Mh | log M*) on LMH for each log M* in logms."""
    h = NB["cosmo"].h
    dn = mass_function.massFunction(10 ** LMH * h, z, mdef="200c", model="tinker08", q_out="dndlnM")
    mu = NB["stellarmass_from_halomass"](LMH, z=z)
    w = dn[None, :] * np.exp(-0.5 * ((np.atleast_1d(logms)[:, None] - mu[None, :]) / scatter) ** 2)
    return w / w.sum(1, keepdims=True)


def draw_conditional(logms, z, scatter=SCATTER):
    cdf = np.cumsum(conditional_weights(logms, z, scatter), 1)
    u_ = rng.random(len(cdf))[:, None]
    return LMH[np.minimum((cdf < u_).sum(1), len(LMH) - 1)]


def inverse(logms, z):
    return np.array([NB["halomass_from_stellarmass"](m, z=z) for m in np.atleast_1d(logms)])


def cluster_m200c():
    h = NB["cosmo"].h
    m200, r200, c200 = mass_adv.changeMassDefinitionCModel(CLUSTER["m500"] * h, CLUSTER["z"], "500c", "200c",
                                                           c_model="dutton14")
    return np.log10(m200 / h)


# ---------------------------------------------------------------- analyses
def conditional_table():
    rows = []
    for z in (0.09, 0.3):
        ms = np.arange(9.0, 12.01, 0.05)
        w = conditional_weights(ms, z)
        cdf = np.cumsum(w, 1)
        q = {p: LMH[(cdf < p).sum(1)] for p in (0.16, 0.5, 0.84)}
        inv = inverse(ms, z)
        for k, m in enumerate(ms):
            rows.append(dict(z=z, logms=m, logmh_inverse=inv[k], logmh_mean=w[k] @ LMH, logmh_p16=q[0.16][k],
                             logmh_median=q[0.5][k], logmh_p84=q[0.84][k]))
    return pd.DataFrame(rows)


def dm_profiles(logms, z=0.09, n=300):
    thetas = np.geomspace(20, 1500, 30)
    inv = inverse(logms, z)[0]
    cond = draw_conditional(np.full(n, logms), z)
    lmc = cluster_m200c()
    rows = []
    for t in thetas:
        d = np.array([halo_dm(m, z, t) for m in cond])
        rows.append(dict(logms=logms, theta=t, dm_inverse=halo_dm(inv, z, t), dm_cond_mean=d.mean(),
                         dm_cond_p16=np.percentile(d, 16), dm_cond_p84=np.percentile(d, 84),
                         dm_cluster_halo=halo_dm(lmc, z, t)))
    return pd.DataFrame(rows)


def galaxy_mc():
    rows, samples = [], {}
    for _, g in GAL.iterrows():
        logm, sig = taylor_logms(g)
        ls = rng.normal(logm, sig, N_MC)
        mh = {"inverse": inverse(ls, g.z), "conditional": draw_conditional(ls, g.z)}
        for k, m in mh.items():
            d = np.array([halo_dm(x, g.z, g.theta) for x in m])
            samples[(g["name"], k)] = d
            rows.append(dict(galaxy=g["name"], z=g.z, theta=g.theta, logms=logm, sigma_logms=sig, method=k,
                             logmh_median=np.median(m), dm_median=np.median(d), dm_mean=d.mean(),
                             dm_p16=np.percentile(d, 16), dm_p84=np.percentile(d, 84), f_inside_r200=(d > 0).mean()))
        print(g["name"], f"log M* = {logm:.2f}", *(f"{r['method']}: Mh {r['logmh_median']:.2f}, DM {r['dm_median']:.1f}"
                                                   for r in rows[-2:]))
    return pd.DataFrame(rows), samples


def scatter_sensitivity(n=2000):
    """Galaxy 6 under different intrinsic scatters of M* at fixed Mh (0 is the inverse)."""
    g = GAL[GAL.name == "6"].iloc[0]
    logm, sig = taylor_logms(g)
    ls = rng.normal(logm, sig, n)
    rows = []
    for s in (0.0, 0.05, 0.10, 0.15, 0.20, 0.25):
        mh = inverse(ls, g.z) if s == 0 else draw_conditional(ls, g.z, s)
        d = np.array([halo_dm(x, g.z, g.theta) for x in mh])
        rows.append(dict(scatter_dex=s, logmh_median=np.median(mh), dm_median=np.median(d)))
    return pd.DataFrame(rows)


def budget(gal):
    med = gal.pivot(index="galaxy", columns="method", values="dm_median")
    photz = sum(DM_PHOTZ.values())
    cases = {
        "paper": dict(halos=135.0, note="Table 5 (124) + central galaxy (11)"),
        "notebook, inverse": dict(halos=photz + DM_CENTRAL_PAPER + med.loc[["4", "6", "7", "8"], "inverse"].sum(),
                                  note="this run of the notebook code, Table 5 galaxies 4, 6, 7, 8"),
        "conditional": dict(halos=photz + DM_CENTRAL_PAPER + med.loc[["4", "6", "7", "8"], "conditional"].sum(),
                            note="P(Mh|M*) for galaxies 4, 6, 7, 8"),
        "conditional, galaxy 6 in ICM": dict(halos=photz + DM_CENTRAL_PAPER + med.loc[["4", "7", "8"], "conditional"].sum(),
                                             note="galaxy 6 is a member of the cluster whose gas DM_ICM already counts"),
    }
    rows = []
    for k, c in cases.items():
        fhalo = DM_ICM + c["halos"]
        sig = np.sqrt(SIG_DISK ** 2 + SIG_MWHALO ** 2 + SIG_IGM ** 2 + (0.5 * c["halos"]) ** 2)
        host = DM_OBS - (DM_DISK + DM_MWHALO + DM_IGM + fhalo)
        rows.append(dict(case=k, dm_c_halo=c["halos"], dm_fhalo=fhalo, dm_host_obs=host, sigma=sig,
                         halpha_minus_budget_sigma=(HALPHA[0] - host) / np.hypot(sig, HALPHA[2] if host > HALPHA[0] else HALPHA[1]),
                         note=c["note"]))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- figure
def figure(cond, prof, samples, bud, gal):
    fig, ax = plt.subplots(2, 2, figsize=(11, 8.6))
    a = ax[0, 0]
    c = cond[cond.z == 0.09]
    a.fill_between(c.logms, c.logmh_p16, c.logmh_p84, color="C0", alpha=0.25, lw=0)
    a.plot(c.logms, c.logmh_median, "C0", label=r"$P(M_h\,|\,M_*)$ median, 16–84%")
    a.plot(c.logms, c.logmh_inverse, "C3", lw=2, label="inverse (notebook)")
    g6 = gal[(gal.galaxy == "6")].set_index("method")
    a.errorbar(g6.logms.iloc[0], g6.loc["inverse", "logmh_median"], fmt="*", ms=13, color="C3", mec="k")
    a.errorbar(g6.logms.iloc[0], g6.loc["conditional", "logmh_median"], fmt="*", ms=13, color="C0", mec="k")
    a.annotate("galaxy 6", (g6.logms.iloc[0], g6.loc["inverse", "logmh_median"]), xytext=(-62, 4),
               textcoords="offset points", fontsize=9)
    lmc = cluster_m200c()
    a.axhline(lmc, color="k", ls="--", lw=1)
    a.text(9.05, lmc + 0.07, f"{CLUSTER['name']}: whole cluster, $M_{{200c}}$ from $M_{{500}}$", fontsize=8)
    a.set(xlim=(9, 12), ylim=(10.5, 16), xlabel=r"$\log M_*$ [$M_\odot$]", ylabel=r"$\log M_{200c}$ [$M_\odot$]",
          title=r"(a) Halo mass assigned to a galaxy of stellar mass $M_*$, $z=0.09$")
    a.legend(loc="upper left", fontsize=8.5, frameon=False)

    a = ax[0, 1]
    p = prof
    a.plot(p.theta, p.dm_inverse, "C3", lw=2, label="inverse (notebook)")
    a.fill_between(p.theta, p.dm_cond_p16, p.dm_cond_p84, color="C0", alpha=0.25, lw=0)
    a.plot(p.theta, p.dm_cond_mean, "C0", label=r"$P(M_h\,|\,M_*)$ mean, 16–84%")
    a.plot(p.theta, p.dm_cluster_halo, "k--", lw=1, label="halo of the whole cluster's mass")
    a.axvline(281.2, color="0.5", lw=0.8)
    a.text(290, 0.75, "galaxy 6\nsightline", fontsize=8, color="0.3")
    at = {k: np.interp(281.2, p.theta, p[k]) for k in ("dm_inverse", "dm_cond_mean", "dm_cluster_halo")}
    a.text(25, 1500, f"at 281\": inverse {at['dm_inverse']:.0f}, conditional {at['dm_cond_mean']:.0f}, "
                     f"cluster-mass halo {at['dm_cluster_halo']:.0f}", fontsize=8)
    a.set(xscale="log", yscale="log", ylim=(0.5, 3000), xlabel=r"impact angle $\theta$ [arcsec]",
          ylabel=r"DM$_{\rm halo}$ (observer) [pc cm$^{-3}$]",
          title=rf"(b) Notebook mNFW DM for $\log M_*={p.logms.iloc[0]:.2f}$ at $z=0.09$")
    a.legend(fontsize=8.5, frameon=False, loc="lower left")

    a = ax[1, 0]
    bins = np.linspace(0, 260, 53)
    for k, col in (("inverse", "C3"), ("conditional", "C0")):
        d = samples[("6", k)]
        a.hist(d, bins, color=col, alpha=0.5, label=f"{k}: median {np.median(d):.0f}, "
                                                    f"{100 * (d > 0).mean():.0f}% inside $r_{{200}}$")
    a.axvline(104, color="k", lw=1)
    a.text(108, a.get_ylim()[1] * 0.6, "paper: 104", fontsize=8.5)
    a.set(xlabel=r"DM$_{\rm halo}$ of galaxy 6 (observer) [pc cm$^{-3}$]", ylabel=f"draws (of {N_MC})",
          title=r"(c) Galaxy 6, notebook Monte Carlo, one line changed")
    a.legend(fontsize=8.5, frameon=False, loc="upper right")

    a = ax[1, 1]
    y = np.arange(len(bud))[::-1]
    a.barh(y, bud.dm_fhalo, color="C1", alpha=0.7, label=r"DM$_{\rm fhalo}$ = ICM (60) + halos")
    a.errorbar(bud.dm_host_obs, y, xerr=bud.sigma, fmt="o", color="k", capsize=3, label=r"budget DM$_{\rm host,obs}$")
    a.axvspan(HALPHA[0] - HALPHA[1], HALPHA[0] + HALPHA[2], color="C2", alpha=0.2, lw=0)
    a.axvline(HALPHA[0], color="C2", lw=1.5, label=r"H$\alpha$ DM$_{\rm host,obs}$ (paper)")
    for yy, r in zip(y, bud.itertuples()):
        a.text(4, yy + 0.27, f"fhalo {r.dm_fhalo:.0f},  host {r.dm_host_obs:.0f}±{r.sigma:.0f}", fontsize=8)
    a.set(yticks=y, yticklabels=bud.case, xlim=(0, 400), xlabel=r"pc cm$^{-3}$ (observer frame)",
          title="(d) FRB 20240114A DM budget")
    a.legend(fontsize=8, frameon=False, loc="lower right")
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"mass_inversion.{ext}", dpi=160)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cond = conditional_table()
    cond.to_csv(OUT / "conditional_mass.csv", index=False)
    gal, samples = galaxy_mc()
    gal.to_csv(OUT / "frb20240114a_galaxies.csv", index=False)
    g6 = gal[gal.galaxy == "6"].logms.iloc[0]
    prof = dm_profiles(g6)
    prof.to_csv(OUT / "dm_profiles.csv", index=False)
    bud = budget(gal)
    bud.to_csv(OUT / "frb20240114a_budget.csv", index=False)
    sens = scatter_sensitivity()
    sens.to_csv(OUT / "scatter_sensitivity.csv", index=False)
    print(sens.round(2).to_string(index=False))
    sep = SkyCoord(*GAL6_POS).separation(SkyCoord(CLUSTER["ra"], CLUSTER["dec"])).arcsec
    kpc = sep * NB["kpc_per_arcsec"](CLUSTER["z"]).value
    print(f"cluster M200c = 10^{cluster_m200c():.2f}; galaxy 6 is {sep:.0f}\" = {kpc:.0f} kpc from the cluster centre")
    print(bud.to_string())
    figure(cond, prof, samples, bud, gal)


if __name__ == "__main__":
    main()
