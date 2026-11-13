"""Audit of the foreground-halo DM estimate and of the DM budget it implies.

1. Calibration. Along random sightlines the mean foreground DM must equal the mean electron density in
   halo gas, which a halo mass function fixes for the profile used by foreground_halo_dm.py (same mNFW
   profile, f_gas, concentration, r200 truncation). That is an upper bound: it counts every halo above
   1e10 Msun once, while the survey only sees some of them.
2. Per-galaxy recomputation of every FRB and control sightline from the cached catalogues, with
   variants that isolate the suspected sources of excess:
     noscat   no stellar-mass scatter
     cond     halo mass drawn from P(Mh | M*), the halo-mass distribution of galaxies of that stellar
              mass, instead of inverting the mean M*(Mh) relation. Above log M* ~ 11 the inverse runs
              0.4-1 dex high (M* ~ Mh^0.3 there, and massive halos are rare), up to 1e15.5 Msun
     central  one halo per group: a galaxy inside r200 of a more massive candidate at a consistent
              redshift is a satellite and its halo is not added again
     cond_central  both, with host-environment galaxies excluded; the corrected estimate, which
              foreground_halo_dm.py now uses (the "pipeline" column is its published-method total)
     hostenv  galaxies at the host redshift inside the host's own r200 are host environment, not
              foreground (FRB sightlines only)
3. DM budget. With the halo term of the sightline known, the consistent subtraction is
     DM_host,obs = DM - DM_MW - (<DM_cosmic> - <DM_fg>_random(z)) - DM_fg,
   because the Macquart mean already contains the average halo gas. Each FRB's residual is compared with
   its uncertainty to test whether the foreground estimate fits in the observed DM.

Outputs: Publication/materials/foreground_halos/audit/
  sightlines.csv, galaxies.csv, calibration.csv, budget.csv, audit.{pdf,png}, AUDIT.md
Run: python scripts/foreground_audit.py   (add --recompute to redo the per-galaxy step)
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_host_property_table as BT  # noqa: E402
import foreground_halo_dm as F  # noqa: E402
import robust_correlations as RC  # noqa: E402
from foreground_halo_dm import (ZGRID, conditional_logmh, flag_groups, halo_electron_fraction,  # noqa: E402
                                halo_model_mean)

OUT = F.OUT / "audit"
N_MC = 500
SIG_MW = 0.2              # fractional uncertainty of NE2001
rng = np.random.default_rng(5)
plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "stix", "font.serif": ["STIXGeneral"]})


# ---------------------------------------------------------------- 2. per-galaxy recomputation
def draws(g, z_host, z_spec):
    if np.isfinite(z_spec):
        return np.full(N_MC, z_spec), 1.0
    sz = max(g.z_phot_std, 0.01)
    lo, hi = (1e-3 - g.z_phot_mean) / sz, (z_host - g.z_phot_mean) / sz
    return (stats.truncnorm.rvs(lo, hi, loc=g.z_phot_mean, scale=sz, size=N_MC, random_state=rng),
            float(stats.norm.cdf(hi) - stats.norm.cdf(lo)))


def recompute_sightline(name, ra, dec, z_host, seed_row=None, host_logms=np.nan):
    dr = F.ls_footprint(name, ra, dec)
    if dr is None:
        return None, []
    cand = F.load_candidates(name, ra, dec, z_host, dr, seed_row, {})
    rows = []
    for _, g in cand.iterrows():
        z_spec = min(float(g.z_best), z_host) if np.isfinite(g.z_best) else np.nan
        z_ref = z_spec if np.isfinite(z_spec) else min(g.z_phot_mean, z_host)
        w1_ok = g.flux_w1 * np.sqrt(max(g.flux_ivar_w1, 0)) > 3
        mags = dict(imag=g.dered_mag_i if g.dered_mag_i > 0 else np.nan,
                    zmag=g.dered_mag_z if g.dered_mag_z > 0 else np.nan, w1=g.dered_mag_w1 if w1_ok else np.nan)
        if not np.isfinite(z_spec):
            m = dict(g=g.dered_mag_g, r=g.dered_mag_r, i=mags["imag"], z=mags["zmag"], w1=mags["w1"])
            if F.rest_mag(m, z_ref, F.BAND_WAVELENGTH["r"]) - F.cosmo.distmod(z_ref).value < F.MR_BRIGHTEST:
                continue
        zs, p_fg = draws(g, z_host, z_spec)
        ms0 = F.logms_bell03_gr(g.dered_mag_g, g.dered_mag_r, zs, **mags)
        sig = np.hypot(np.hypot(1.097 * F.mag_err(g.flux_g, g.flux_ivar_g), 1.497 * F.mag_err(g.flux_r, g.flux_ivar_r)),
                       F.MASS_SYS_DEX)
        ms = ms0 + rng.normal(0, sig, N_MC)
        dm_full = F.dm_mnfw(F.logmh_from_logms(ms, zs), zs, g.theta)
        dm_nos = F.dm_mnfw(F.logmh_from_logms(ms0, zs), zs, g.theta)
        dm_cond = F.dm_mnfw(conditional_logmh(ms, zs), zs, g.theta)
        logms_ref = float(F.logms_bell03_gr(g.dered_mag_g, g.dered_mag_r, z_ref, **mags))
        rows.append(galaxy_row(name, g.ls_id, g.ra, g.dec, g.theta, g.dered_mag_r, z_ref, z_spec,
                               0.0 if np.isfinite(z_spec) else max(g.z_phot_std, 0.01), p_fg, logms_ref,
                               p_fg * dm_full.mean(), p_fg * dm_nos.mean(), p_fg * dm_cond.mean()))
    sga = F.sga_wide(name, ra, dec, z_host)
    for _, s in sga[sga.theta_deg > F.SEARCH_DEG].iterrows():
        gm, rm = F.ab_mag(s.flux_g, s.mw_transmission_g), F.ab_mag(s.flux_r, s.mw_transmission_r)
        if not (np.isfinite(gm) and np.isfinite(rm)):
            continue
        logms = float(F.logms_bell03_gr(gm, rm, s.z_leda))
        logmh = float(F.logmh_from_logms(logms, s.z_leda)[0])
        theta = s.theta_deg * 3600
        if theta * F.kpc_per_arcsec(s.z_leda) >= F.r200_kpc(10 ** logmh, s.z_leda):
            continue
        ms = logms + rng.normal(0, np.hypot(0.05 * 1.5, F.MASS_SYS_DEX), N_MC)
        zs = np.full(N_MC, s.z_leda)
        rows.append(galaxy_row(name, f"SGA {s.galaxy}", s.ra, s.dec, theta, rm, s.z_leda, s.z_leda, 0.0, 1.0, logms,
                               F.dm_mnfw(F.logmh_from_logms(ms, zs), zs, theta).mean(),
                               F.dm_mnfw(logmh, s.z_leda, theta)[0],
                               F.dm_mnfw(conditional_logmh(ms, zs), zs, theta).mean()))
    gal = pd.DataFrame(rows)
    if len(gal):
        flag_groups(gal, z_host, host_logms)
    return dr, gal


def galaxy_row(name, ls_id, ra, dec, theta, r, z_ref, z_spec, z_err, p_fg, logms, dm, dm_noscat, dm_cond):
    logmh = float(conditional_logmh(logms, z_ref, draw=False))
    return dict(sightline=name, ls_id=str(ls_id), ra=ra, dec=dec, theta_arcsec=theta, r=r, z_ref=z_ref, z_spec=z_spec,
                z_err=z_err, p_fg=p_fg, logms=logms, logmh_inverse=float(F.logmh_from_logms(logms, z_ref)[0]),
                logmh=logmh, r200_kpc=float(F.r200_kpc(10 ** logmh, z_ref)),
                dm=float(dm), dm_noscat=float(dm_noscat), dm_cond=float(dm_cond))


def control_positions(hosts, done):
    """The positions drawn by foreground_halo_dm.controls (same seed and order)."""
    r = np.random.default_rng(7)
    pos = {}
    for _, row in hosts[hosts.frb.isin(done[done.in_ls.astype(str).isin(["dr9", "dr10"])].index)].iterrows():
        for k in range(2):
            pa = r.uniform(0, 2 * np.pi)
            dec = row.dec_deg + 0.75 * np.cos(pa)
            pos[f"{row.frb}_ctrl{k}"] = ((row.ra_deg + 0.75 * np.sin(pa) / np.cos(np.radians(dec))) % 360, dec)
    return pos


def recompute_all():
    hosts = pd.read_csv(F.HOSTS)
    seed = pd.read_csv(F.SEED).set_index("frb")
    ctrl = pd.read_csv(F.OUT / "foreground_dm_controls.csv")
    lines, gals = [], []

    def add(kind, name, host, z, dr, gal, pipeline):
        tot = dict.fromkeys(["dm", "dm_noscat", "dm_cond", "central", "cond_central", "hostenv", "hostenv_cond"], 0.0)
        if len(gal):
            keep = ~gal.satellite & ~gal.hostenv
            tot.update(dm=gal.dm.sum(), dm_noscat=gal.dm_noscat.sum(), dm_cond=gal.dm_cond.sum(),
                       central=gal.dm[keep].sum(), cond_central=gal.dm_cond[keep].sum(),
                       hostenv=gal.dm[gal.hostenv].sum(), hostenv_cond=gal.dm_cond[gal.hostenv].sum())
            gals.append(gal.assign(kind=kind, host=host))
        lines.append(dict(kind=kind, sightline=name, host=host, z=z, dr=dr, n_gal=len(gal), pipeline=pipeline,
                          **{k: float(v) for k, v in tot.items()}))
        print(f"{name}: pipeline {pipeline:6.1f}  recomputed {tot['dm']:6.1f}  noscat {tot['dm_noscat']:6.1f}  "
              f"cond {tot['dm_cond']:6.1f}  central {tot['central']:6.1f}  cond+central {tot['cond_central']:6.1f}",
              flush=True)

    done = pd.read_csv(F.OUT / "foreground_dm_per_host.csv", dtype={"frb": str}).set_index("frb")
    for _, row in hosts[hosts.z.notna()].iterrows():
        sr = seed.loc[row.frb] if row.frb in seed.index else None
        dr, gal = recompute_sightline(row.frb, row.ra_deg, row.dec_deg, row.z, sr, row.logmstar)
        if dr is not None:
            add("frb", row.frb, row.frb, row.z, dr, gal, done.loc[row.frb, "dm_fg_published"])
    pos = control_positions(hosts, done)
    for _, c in ctrl.iterrows():
        dr, gal = recompute_sightline(c.frb, *pos[c.frb], c.z)
        if dr is not None:
            add("control", c.frb, c.host, c.z, dr, gal, c.dm_fg_published)
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(lines).to_csv(OUT / "sightlines.csv", index=False)
    pd.concat(gals, ignore_index=True).to_csv(OUT / "galaxies.csv", index=False)


# ---------------------------------------------------------------- 3. calibration and budget
VARIANTS = {"pipeline": "pipeline (as published)", "dm_noscat": "no mass scatter",
            "central": "one halo per group", "dm_cond": r"$M_h$ from $P(M_h|M_*)$",
            "cond_central": r"$P(M_h|M_*)$ + one halo per group"}
MW_HALO, MW_HALO_ERR = 40.0, 20.0


def macquart_sd(z, F_=0.2):
    """Half the 16-84% range of DM_cosmic under the Macquart et al. (2020) distribution (its variance
    diverges logarithmically with the Delta^-3 tail)."""
    mean = BT.mean_dm_cosmic(z)
    sigma, d = F_ / np.sqrt(z), np.linspace(1e-3, 6.0, 60000)

    def pdf(c0):
        p = d ** -3 * np.exp(-((d ** -3 - c0) ** 2) / (18 * sigma ** 2))
        return p / np.trapezoid(p, d)

    lo, hi = -5.0, 5.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if np.trapezoid(d * pdf(mid), d) > 1 else (lo, mid)
    cdf = np.cumsum(pdf(0.5 * (lo + hi))) * (d[1] - d[0])
    return mean * 0.5 * (d[np.searchsorted(cdf, 0.84)] - d[np.searchsorted(cdf, 0.16)])


def calibration(lines, frac):
    c = lines[lines.kind == "control"].copy()
    c["model"] = halo_model_mean(c.z.to_numpy(), frac)
    rows = []
    for v in VARIANTS:
        boot = [c[v].iloc[i].sum() / c.model.iloc[i].sum() for i in rng.integers(0, len(c), (2000, len(c)))]
        rows.append(dict(variant=v, label=VARIANTS[v], ratio_to_halo_model=c[v].sum() / c.model.sum(),
                         ratio_lo=np.percentile(boot, 16), ratio_hi=np.percentile(boot, 84)))
    return pd.DataFrame(rows), c


def frb_vs_random(lines, frac):
    """Foreground of FRB sightlines against random ones, each normalized by the halo-model mean at its z."""
    s = lines[lines.n_gal.notna() & (lines.z > 0.05)].copy()
    s["model"] = halo_model_mean(s.z.to_numpy(), frac)
    rows = []
    for v in ("pipeline", "cond_central"):
        x = (s[v] / s.model).to_numpy()
        is_frb = (s.kind == "frb").to_numpy()
        obs = x[is_frb].mean() - x[~is_frb].mean()
        perm = np.array([(lambda m: x[m].mean() - x[~m].mean())(rng.permutation(is_frb)) for _ in range(10_000)])
        rows.append(dict(variant=v, frb_over_model=x[is_frb].mean(), random_over_model=x[~is_frb].mean(),
                         p_perm=(1 + np.sum(perm >= obs)) / 10_001,
                         p_ad=stats.anderson_ksamp([x[is_frb], x[~is_frb]]).pvalue))
    return pd.DataFrame(rows)


def budget(lines, frac, cal, ctrl):
    h = RC.load()
    f = lines[lines.kind == "frb"].set_index("sightline")
    h = h[h.frb.isin(f.index) & h.dm_host_rest.notna()].copy()
    h["mean_c"] = h.z.map(BT.mean_dm_cosmic)
    h["med_c"] = h.z.map(BT.median_dm_cosmic)
    h["sd_c"] = h.z.map(macquart_sd)
    model = halo_model_mean(h.z.to_numpy(), frac)
    for v, tag in (("pipeline", "pipe"), ("cond_central", "corr")):
        ratio = float(cal.set_index("variant").loc[v, "ratio_to_halo_model"])
        fg = h.frb.map(f[v]).to_numpy(float)
        expect = ratio * model
        resid_sd = np.interp(h.z, *zip(*[(ctrl[(ctrl.z > a) & (ctrl.z <= b)].z.mean(),
                                          ctrl[(ctrl.z > a) & (ctrl.z <= b)][v].std())
                                         for a, b in ((0, .15), (.15, .3), (.3, .45), (.45, .7))]))
        sd_diffuse = np.sqrt(np.clip(h.sd_c ** 2 - resid_sd ** 2, (0.1 * h.mean_c) ** 2, None))
        other = (SIG_MW * h.dm_mw) ** 2 + (0.5 * fg) ** 2 + MW_HALO_ERR ** 2
        h[f"fg_{tag}"] = fg
        h[f"fg_expected_{tag}"] = expect
        h[f"host_obs_{tag}"] = h.dm_obs - h.dm_mw - (h.mean_c - expect) - fg
        h[f"host_obs_{tag}_mwhalo"] = h[f"host_obs_{tag}"] - MW_HALO
        h[f"sigma_{tag}"] = np.sqrt(other + sd_diffuse ** 2)
        h[f"zscore_{tag}"] = h[f"host_obs_{tag}_mwhalo"] / h[f"sigma_{tag}"]
        h[f"zscore_{tag}_floor"] = h[f"host_obs_{tag}_mwhalo"] / np.sqrt(other + (0.1 * h.mean_c) ** 2)
        h[f"dm_host_{tag}"] = h[f"host_obs_{tag}"] * (1 + h.z)
        h[f"dm_host_{tag}_excess"] = h.dm_host_rest - (fg - expect) * (1 + h.z)
    h["host_obs_standard"] = h.dm_obs - h.dm_mw - h.med_c
    h["host_obs_full_subtraction"] = h.host_obs_standard - h.fg_pipe
    return h.sort_values("z")


def correlation_impact(b):
    rows = []
    tests = [("q", "all"), ("log_sfr", "all"), ("log_sfr", "z > 0.29"), ("log_abs_rm", "all"), ("logmstar", "all"),
             ("z", "all")]
    b = b.assign(log_sfr=np.log10(b.sfr.where(b.sfr > 0)))
    for x, sample in tests:
        s = b if sample == "all" else b[b.z > 0.29]
        for col, lab in (("dm_host_rest", "standard (paper)"), ("dm_host_fg_published", "published foreground fully removed"),
                         ("dm_host_pipe", "pipeline, consistent subtraction"),
                         ("dm_host_corr_excess", "standard minus corrected excess"),
                         ("dm_host_corr", "corrected foreground, consistent subtraction")):
            ok = s[x].notna() & s[col].notna()
            r, p = stats.spearmanr(s.loc[ok, x], s.loc[ok, col])
            rows.append(dict(x=x, sample=sample, dm_host=lab, N=int(ok.sum()), r_S=r, p_S=p,
                             p_AD=RC.ad_split(s.loc[ok, x].to_numpy(float), s.loc[ok, col].to_numpy(float))))
    return pd.DataFrame(rows)


def figure(cal_ctrl, frac, b, gal):
    fig, ax = plt.subplots(2, 2, figsize=(14, 10.5))
    a = ax[0, 0]
    zz = np.linspace(0.01, 0.75, 60)
    a.plot(zz, halo_model_mean(zz, frac), "k-", lw=2, label=r"halo model, every halo $>10^{10}\,M_\odot$ (upper bound)")
    for v, col, mk in (("pipeline", "#D62728", "o"), ("cond_central", "#1F77B4", "s")):
        a.scatter(cal_ctrl.z, cal_ctrl[v], s=10, color=col, alpha=0.35, marker=mk)
        edges = [0, .15, .3, .45, .7]
        mids, means, errs = [], [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            s = cal_ctrl[(cal_ctrl.z > lo) & (cal_ctrl.z <= hi)][v]
            mids.append(cal_ctrl[(cal_ctrl.z > lo) & (cal_ctrl.z <= hi)].z.mean())
            means.append(s.mean())
            errs.append(s.std() / np.sqrt(len(s)))
        a.errorbar(mids, means, errs, color=col, marker=mk, ms=7, lw=1.8, capsize=3,
                   label=f"random sightlines, {VARIANTS[v]}: mean")
    a.set_yscale("symlog", linthresh=10)
    a.set_xlabel("redshift of the sightline end")
    a.set_ylabel(r"foreground halo DM (pc cm$^{-3}$, observer frame)")
    a.set_title("(a) Random sightlines against the most the halos can supply", fontsize=10)
    a.legend(fontsize=7.5, frameon=False, loc="upper left")

    a = ax[0, 1]
    ms = np.linspace(9.5, 12.2, 100)
    for z, ls in ((0.1, "-"), (0.5, "--")):
        a.plot(ms, F.logmh_from_logms(ms, z), color="#D62728", ls=ls, label=f"inverse of mean $M_*(M_h)$, z = {z}")
        a.plot(ms, conditional_logmh(ms, np.full_like(ms, z), draw=False), color="#1F77B4", ls=ls,
               label=f"mean $M_h$ at fixed $M_*$, z = {z}")
    g = gal[(gal.kind == "frb") & (gal.dm > 5)]
    a.scatter(g.logms, g.logmh_inverse, s=12, color="#D62728", alpha=0.4, label="FRB foreground galaxies with DM > 5")
    a.set_xlabel(r"$\log M_*/M_\odot$")
    a.set_ylabel(r"$\log M_h/M_\odot$")
    a.set_title("(b) The halo mass assigned to massive galaxies", fontsize=10)
    a.legend(fontsize=7.5, frameon=False, loc="upper left")

    a = ax[1, 0]
    x = np.arange(len(b))
    a.axhline(0, color="k", lw=0.8)
    a.scatter(x - 0.2, b.host_obs_full_subtraction - MW_HALO, s=16, color="0.6", marker="x",
              label="standard DM$_{\\rm host}$ minus full pipeline foreground (double counted)")
    a.errorbar(x, b.host_obs_pipe_mwhalo, b.sigma_pipe, fmt="o", ms=4, color="#D62728", lw=0.8, capsize=0,
               label="pipeline foreground, consistent subtraction")
    a.errorbar(x + 0.2, b.host_obs_corr_mwhalo, b.sigma_corr, fmt="s", ms=4, color="#1F77B4", lw=0.8, capsize=0,
               label="corrected foreground, consistent subtraction")
    a.set_xticks(x)
    a.set_xticklabels([f"{r.frb} ({r.z:.2f})" for r in b.itertuples()], rotation=90, fontsize=6)
    a.set_ylabel(r"DM$_{\rm host}$/(1+z) left for the host (pc cm$^{-3}$, observer frame)")
    a.set_title(r"(c) DM budget closure per FRB (DM$_{\rm MW,halo}$ = 40 subtracted; 1$\sigma$ errors)", fontsize=10)
    a.legend(fontsize=7.5, frameon=False, loc="upper left")

    a = ax[1, 1]
    bins = np.linspace(-4, 6.5, 22)
    a.hist(b.zscore_pipe.clip(upper=6.2), bins, histtype="step", lw=2, color="#D62728",
           label=f"pipeline: {(b.zscore_pipe < -2).sum()} below $-2\\sigma$, {(b.host_obs_pipe_mwhalo < 0).sum()} negative")
    a.hist(b.zscore_corr.clip(upper=6.2), bins, histtype="step", lw=2, color="#1F77B4",
           label=f"corrected: {(b.zscore_corr < -2).sum()} below $-2\\sigma$, {(b.host_obs_corr_mwhalo < 0).sum()} negative")
    a.axvline(-2, color="k", ls="--", lw=1)
    a.set_xlabel(r"(DM left for the host) / $\sigma$  (last bin: $\geq 6$)")
    a.set_ylabel("FRBs")
    a.set_title("(d) Does the foreground fit inside the observed DM?", fontsize=10)
    a.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"audit.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)


def analyse():
    lines = pd.read_csv(OUT / "sightlines.csv", dtype={"sightline": str, "host": str})
    gal = pd.read_csv(OUT / "galaxies.csv", dtype={"sightline": str, "host": str})
    frac = np.array([halo_electron_fraction(z) for z in ZGRID])
    cal, ctrl = calibration(lines, frac)
    b = budget(lines, frac, cal, ctrl)
    imp = correlation_impact(b)
    fvr = frb_vs_random(lines, frac)
    fvr.to_csv(OUT / "frb_vs_random.csv", index=False)
    print(fvr.round(3).to_string(index=False))
    cal.to_csv(OUT / "calibration.csv", index=False)
    b.to_csv(OUT / "budget.csv", index=False)
    imp.to_csv(OUT / "correlation_impact.csv", index=False)
    figure(ctrl, frac, b, gal)
    pd.set_option("display.width", 250)
    print(cal.round(3).to_string(index=False))
    print(b[["frb", "z", "dm_obs", "mean_c", "fg_pipe", "fg_corr", "fg_expected_pipe", "fg_expected_corr",
             "host_obs_standard", "host_obs_full_subtraction", "host_obs_pipe_mwhalo", "host_obs_corr_mwhalo",
             "sigma_pipe", "sigma_corr", "zscore_pipe", "zscore_corr", "zscore_corr_floor"]].round(1).to_string(index=False))
    print(imp.round(3).to_string(index=False))


if __name__ == "__main__":
    if "--recompute" in sys.argv or not (OUT / "sightlines.csv").exists():
        recompute_all()
    analyse()
