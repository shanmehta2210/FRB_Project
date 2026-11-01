"""Foreground (intervening) galaxy-halo DM for the 56-host sample, following
Bhardwaj, Snelders et al. (2025, arXiv:2506.11915) Appendix B and the root notebooks
(Estimate_DM_halo.ipynb cell 26, Galaxy_offset_code.ipynb, atel_analysis.ipynb).

Selection (Legacy Surveys DR9 only):
  * Tractor galaxies (type not PSF/DUP) within 0.5 deg of the FRB,
  * LS DR9 photo-z (Zhou+2021) with z_phot - 2 sigma_z <= z_host,
  * R_eff / theta > 0.018 (stellar-to-halo radius ratio),
  * host galaxy (and its Tractor shreds) excluded.
Spec-z from the LS photo-z compilation, DESI DR1, SDSS DR17 and SGA-2020 (HyperLEDA);
galaxies with z_spec > z_host are dropped. Photo-z-only galaxies use a Gaussian
N(z_phot, sigma_z) truncated to z <= z_host (paper convention).
Mass: Bell+2003 g-r M/L (Chabrier zero-point -0.406, as in atel_analysis.ipynb), since
LS DR9 has no i band; Moster+2013 SHMR; Dutton & Maccio 2014 c(M); modified NFW gas
profile (Prochaska & Zheng 2019; alpha=y0=2, f_gas=0.75); DM_obs = DM/(1+z_gal).
A 50% one-sigma systematic is applied to each host total.
SGA-2020 galaxies with z < z_host out to 30 deg are checked for halos that reach the
sightline from beyond 0.5 deg (the paper's HECATE check).
"""
from __future__ import annotations

import io
import time
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import constants as const
from astropy import units as u
from astropy.cosmology import Planck18 as cosmo
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Publication" / "materials" / "foreground_halos"
CACHE = OUT / "ls_cache"
HOSTS = ROOT / "plots" / "plots_null" / "v2" / "host_correlations" / "host_props.csv"
SEED = ROOT / "Publication" / "materials" / "host_properties" / "seed_56.csv"

SEARCH_DEG = 0.5
SRHR = 0.018
SGA_DEG = 30.0
N_MC = 2000
SYS_DEX = 0.10
SYS_FRAC = 0.50
MR_MIN = -24.0
RNG = np.random.default_rng(42)

H = cosmo.h
FB = cosmo.Ob0 / cosmo.Om0
MP = const.m_p.cgs.value
KPC_CM = u.kpc.to(u.cm)
KPC_PC = 1e3
LOGMH_GRID = np.linspace(9.0, 15.5, 1301)


# ---------------------------------------------------------------- halo physics
def moster13_logms(logmh, z):
    zf = z / (1 + z)
    n = 0.0351 - 0.0247 * zf
    beta = 1.376 - 0.826 * zf
    gamma = 0.608 + 0.329 * zf
    logm1 = 11.59 + 1.195 * zf
    x = 10 ** (logmh - logm1)
    return logmh + np.log10(2 * n) - np.log10(x ** (-beta) + x ** gamma)


def logmh_from_logms(logms, z):
    logms, z = np.broadcast_arrays(np.atleast_1d(logms), np.atleast_1d(z))
    grid = moster13_logms(LOGMH_GRID[None, :], z.ravel()[:, None])
    out = np.array([np.interp(t, g, LOGMH_GRID) for t, g in zip(logms.ravel(), grid)])
    return out.reshape(logms.shape)


def dutton14_c200c(m_msun, z):
    b = -0.101 + 0.026 * z
    a = 0.520 + (0.905 - 0.520) * np.exp(-0.617 * z ** 1.21)
    return 10 ** (a + b * np.log10(m_msun * H / 1e12))


def rho_crit(z):
    return cosmo.critical_density(z).to(u.g / u.cm ** 3).value


def r200_kpc(m_msun, z):
    m_g = m_msun * const.M_sun.cgs.value
    return (3 * m_g / (4 * np.pi * 200 * rho_crit(z))) ** (1 / 3) / KPC_CM


def kpc_per_arcsec(z):
    return cosmo.kpc_proper_per_arcmin(z).to(u.kpc / u.arcsec).value


def dm_mnfw(logmh, z, theta_arcsec, alpha=2.0, y0=2.0, f_gas=0.75, mu_h=1.3, wt=1.1667, npts=801):
    """Observer-frame DM of a modified-NFW halo crossed at angular offset theta (vectorised)."""
    logmh, z, theta = np.broadcast_arrays(*(np.atleast_1d(np.asarray(v, float)) for v in (logmh, z, theta_arcsec)))
    m = 10 ** logmh
    c = dutton14_c200c(m, z)
    r200 = r200_kpc(m, z)
    d = theta * kpc_per_arcsec(z)
    fc = np.log(1 + c) - c / (1 + c)
    rho0 = (200 * rho_crit(z) / 3) * c ** 3 / fc
    inside = d < r200
    half = np.sqrt(np.clip(r200 ** 2 - d ** 2, 0, None))
    t = np.linspace(-1, 1, npts)
    x = half[:, None] * t[None, :]
    r = np.sqrt(x ** 2 + d[:, None] ** 2)
    y = c[:, None] * r / r200[:, None]
    rho = f_gas * FB * rho0[:, None] / (y ** (1 - alpha) * (y0 + y) ** (2 + alpha))
    ne = rho / (mu_h * MP) * wt
    dm = np.trapezoid(ne, x * KPC_PC, axis=1)
    return np.where(inside, dm / (1 + z), 0.0)


# ---------------------------------------------------------------- stellar mass
LAM = {"g": 4810.0, "r": 6440.0, "i": 7840.0, "z": 9260.0, "w1": 33680.0}


def rest_mag(mags, z, lam_rest):
    """Observed AB SED (log-lambda interpolation) evaluated at lam_rest*(1+z); K = -2.5 log(1+z)."""
    bands = [b for b in ("g", "r", "i", "z", "w1") if np.isfinite(mags.get(b, np.nan))]
    xp = np.log10([LAM[b] for b in bands])
    fp = np.array([mags[b] for b in bands])
    return np.interp(np.log10(lam_rest * (1 + np.asarray(z))), xp, fp) + 2.5 * np.log10(1 + np.asarray(z))


def logms_bell03_gr(g, r, z, imag=np.nan, zmag=np.nan, w1=np.nan):
    mags = dict(g=g, r=r, i=imag, z=zmag, w1=w1)
    g0, r0 = rest_mag(mags, z, LAM["g"]), rest_mag(mags, z, LAM["r"])
    mr = r0 - cosmo.distmod(z).value
    return 1.097 * (g0 - r0) - 0.406 - 0.4 * (mr - 4.65)


def logms_taylor11_gi(g, i, z):
    mi = i - cosmo.distmod(z).value
    return 1.15 + 0.70 * (g - i) - 0.4 * mi


def galaxy_dm_samples(g, r, g_err, r_err, theta, z_host, z_spec=np.nan, z_phot=np.nan, z_phot_err=np.nan,
                      n=N_MC, mass_fn=logms_bell03_gr, **mags):
    """MC samples of a galaxy's observer-frame DM; returns (samples, P_foreground)."""
    if np.isfinite(z_spec):
        zs = np.full(n, z_spec)
        p_fg = 1.0
    else:
        lo, hi = (1e-3 - z_phot) / z_phot_err, (z_host - z_phot) / z_phot_err
        zs = stats.truncnorm.rvs(lo, hi, loc=z_phot, scale=z_phot_err, size=n, random_state=RNG)
        p_fg = float(stats.norm.cdf(hi) - stats.norm.cdf(lo))
    stat = np.hypot(1.097 * g_err, 0.659 * r_err) if mass_fn is logms_bell03_gr else np.hypot(0.70 * g_err, 1.1 * r_err)
    sig = np.hypot(stat, SYS_DEX)
    logms = mass_fn(g, r, zs, **mags) + RNG.normal(0, sig, n)
    logmh = logmh_from_logms(logms, zs)
    return dm_mnfw(logmh, zs, theta), p_fg


# ---------------------------------------------------------------- catalogue access
def dl_query(sql, retries=3):
    from dl import queryClient as qc
    for k in range(retries):
        try:
            return pd.read_csv(io.StringIO(qc.query(sql=sql, fmt="csv", timeout=600)))
        except Exception as e:  # noqa: BLE001
            if k == retries - 1:
                raise
            print(f"   retry ({e.__class__.__name__}: {str(e)[:80]})")
            time.sleep(5)


def ls_candidates(frb, ra, dec, z_host, dr):
    f = CACHE / f"{frb}_ls.csv"
    if f.exists():
        return pd.read_csv(f)
    imag = "t.dered_mag_i," if dr == "dr10" else ""
    sql = f"""
    SELECT t.ls_id, t.ra, t.dec, t.type, t.shape_r, t.dered_mag_g, t.dered_mag_r, {imag} t.dered_mag_z, t.dered_mag_w1,
           t.flux_g, t.flux_r, t.flux_ivar_g, t.flux_ivar_r, t.flux_w1, t.flux_ivar_w1,
           p.z_phot_mean, p.z_phot_std, p.z_spec,
           q3c_dist(t.ra, t.dec, {ra}, {dec})*3600 AS theta
    FROM ls_{dr}.tractor AS t JOIN ls_{dr}.photo_z AS p ON t.ls_id = p.ls_id
    WHERE q3c_radial_query(t.ra, t.dec, {ra}, {dec}, {SEARCH_DEG})
      AND t.type <> 'PSF' AND t.type <> 'DUP' AND t.brick_primary = 1
      AND t.shape_r > {SRHR} * q3c_dist(t.ra, t.dec, {ra}, {dec})*3600
      AND p.z_phot_mean - 2*p.z_phot_std <= {z_host}
    """
    df = dl_query(sql)
    df.to_csv(f, index=False)
    return df


def ls_footprint(frb, ra, dec):
    """'dr9' if the field is in LS DR9 (as in the notebooks), else 'dr10' if in DR10, else None."""
    f = CACHE / f"{frb}_footprint.csv"
    if f.exists():
        d = pd.read_csv(f).dr.iloc[0]
        return d if isinstance(d, str) else None
    dr = None
    for rel in ("dr9", "dr10"):
        n = dl_query(f"""SELECT COUNT(*) AS n FROM ls_{rel}.tractor AS t JOIN ls_{rel}.photo_z AS p ON t.ls_id = p.ls_id
                         WHERE q3c_radial_query(t.ra, t.dec, {ra}, {dec}, 0.02)""").n.iloc[0]
        if n > 0:
            dr = rel
            break
    pd.DataFrame({"dr": [dr]}).to_csv(f, index=False)
    return dr


def specz(frb, ra, dec):
    f = CACHE / f"{frb}_specz.csv"
    if f.exists():
        return pd.read_csv(f)
    parts = []
    for src, sql in {
        "DESI DR1": f"""SELECT mean_fiber_ra AS ra, mean_fiber_dec AS dec, z FROM desi_dr1.zpix
            WHERE q3c_radial_query(mean_fiber_ra, mean_fiber_dec, {ra}, {dec}, {SEARCH_DEG})
            AND zwarn = 0 AND spectype = 'GALAXY' AND zcat_primary""",
        "SDSS DR17": f"""SELECT ra, dec, z FROM sdss_dr17.specobj
            WHERE q3c_radial_query(ra, dec, {ra}, {dec}, {SEARCH_DEG}) AND zwarning = 0 AND class = 'GALAXY'""",
        "SGA-2020": f"""SELECT ra, dec, z_leda AS z FROM sga2020.ellipse
            WHERE q3c_radial_query(ra, dec, {ra}, {dec}, {SEARCH_DEG}) AND z_leda > 0""",
    }.items():
        try:
            d = dl_query(sql)
            d["src"] = src
            parts.append(d)
        except Exception as e:  # noqa: BLE001
            print(f"   {src} spec-z query failed: {str(e)[:100]}")
    df = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["ra", "dec", "z", "src"])
    df.to_csv(f, index=False)
    return df


def sga_wide(frb, ra, dec, z_host):
    f = CACHE / f"{frb}_sga.csv"
    if f.exists():
        return pd.read_csv(f)
    sql = f"""
    SELECT e.sga_id, e.galaxy, e.ra, e.dec, e.z_leda, e.d26,
           t.flux_g, t.flux_r, t.mw_transmission_g, t.mw_transmission_r,
           q3c_dist(e.ra, e.dec, {ra}, {dec}) AS theta_deg
    FROM sga2020.ellipse AS e JOIN sga2020.tractor AS t ON e.sga_id = t.sga_id
    WHERE q3c_radial_query(e.ra, e.dec, {ra}, {dec}, {SGA_DEG})
      AND e.z_leda > 0.0005 AND e.z_leda <= {z_host}
    """
    df = dl_query(sql)
    df.to_csv(f, index=False)
    return df


# ---------------------------------------------------------------- per host
def match_specz(cand, sz, tol_arcsec=1.5):
    cand = cand.copy()
    cand["z_best"], cand["z_best_src"] = np.nan, ""
    ok = cand.z_spec.to_numpy(float) > 0
    cand.loc[ok, "z_best"] = cand.loc[ok, "z_spec"]
    cand.loc[ok, "z_best_src"] = "LS photo-z compilation"
    if len(sz) and len(cand):
        dra = (cand.ra.to_numpy()[:, None] - sz.ra.to_numpy()[None, :]) * np.cos(np.radians(cand.dec.to_numpy()))[:, None]
        sep = np.hypot(dra, cand.dec.to_numpy()[:, None] - sz.dec.to_numpy()[None, :]) * 3600
        j = sep.argmin(1)
        tol = np.maximum(tol_arcsec, 0.2 * np.nan_to_num(cand.shape_r.to_numpy(float)))
        hit = (sep[np.arange(len(cand)), j] < tol) & ~np.isfinite(cand.z_best.to_numpy(float))
        cand.loc[hit, "z_best"] = sz.z.to_numpy()[j[hit]]
        cand.loc[hit, "z_best_src"] = sz.src.to_numpy()[j[hit]]
    return cand


def identify_host(cand, m_r_host, re_host):
    """Host = LS source within max(5", 3 R_e) of the FRB closest in r magnitude to the GALFIT host."""
    rad = max(5.0, 3 * re_host) if np.isfinite(re_host) else 5.0
    near = cand[cand.theta < rad]
    if not len(near):
        return None
    if np.isfinite(m_r_host):
        dm = (near.dered_mag_r - m_r_host).abs()
        if dm.min() < 1.5:
            return near.loc[dm.idxmin()]
    return near.loc[near.theta.idxmin()]


def mag_err(flux, ivar):
    with np.errstate(divide="ignore", invalid="ignore"):
        e = 1.0857 / (flux * np.sqrt(ivar))
    return np.where(np.isfinite(e) & (e > 0), np.minimum(e, 1.0), 0.3)


def run_host(row, seed_row, n_mc=N_MC):
    frb, ra, dec, zh = row.frb, row.ra_deg, row.dec_deg, row.z
    res = dict(frb=frb, z=zh, in_ls=False, n_cand=0, n_specz=0, host_ls_id="", dm_fg=0.0, dm_fg_conditional=0.0,
               dm_fg_specz_only=0.0, dm_fg_sga_wide=0.0, dm_fg_mc_sd=0.0, top_gal="", top_dm=0.0)
    dr = ls_footprint(frb, ra, dec)
    if dr is None:
        return res, []
    res["in_ls"] = dr
    cand = ls_candidates(frb, ra, dec, zh, dr)
    if "dered_mag_i" not in cand:
        cand["dered_mag_i"] = np.nan
    for c in cand.columns.drop(["ls_id", "type"]):
        cand[c] = pd.to_numeric(cand[c], errors="coerce")
    cand = match_specz(cand, specz(frb, ra, dec))
    re_host = float(seed_row.re_arcsec_galfit) if seed_row is not None else np.nan
    m_r_host = float(seed_row.m_r_galfit) if seed_row is not None else np.nan
    host = identify_host(cand, m_r_host, re_host) if seed_row is not None else None
    if host is not None:
        res["host_ls_id"] = str(host.ls_id)
        hra, hdec = host.ra, host.dec
        rex = max(2 * (re_host if np.isfinite(re_host) else 1.0), 2 * float(host.shape_r), 2.0)
        sep_h = np.hypot((cand.ra - hra) * np.cos(np.radians(hdec)), cand.dec - hdec) * 3600
        zb = cand.z_best.to_numpy(float)
        shred = (sep_h < rex) & ~(np.isfinite(zb) & (np.abs(zb - zh) > 0.01))
        cand = cand[~shred | (cand.ls_id == host.ls_id)]
        cand = cand[cand.ls_id != host.ls_id]
    zb = cand.z_best.to_numpy(float)
    cand = cand[~(np.isfinite(zb) & (zb > zh + 0.001))]
    cand = cand[(cand.dered_mag_g > 0) & (cand.dered_mag_r > 0) & np.isfinite(cand.dered_mag_r)]
    res["n_cand"] = len(cand)
    res["n_specz"] = int(np.isfinite(cand.z_best.to_numpy(float)).sum())

    rows, tot, tot_c, tot_s = [], np.zeros(n_mc), 0.0, 0.0
    n_bad = 0
    for _, g in cand.iterrows():
        zsp = min(float(g.z_best), zh) if np.isfinite(g.z_best) else np.nan
        w1_ok = g.flux_w1 * np.sqrt(max(g.flux_ivar_w1, 0)) > 3
        extra = dict(imag=g.dered_mag_i if g.dered_mag_i > 0 else np.nan,
                     zmag=g.dered_mag_z if g.dered_mag_z > 0 else np.nan, w1=g.dered_mag_w1 if w1_ok else np.nan)
        if not np.isfinite(zsp):
            zr = min(g.z_phot_mean, zh)
            mags = dict(g=g.dered_mag_g, r=g.dered_mag_r, i=extra["imag"], z=extra["zmag"], w1=extra["w1"])
            if rest_mag(mags, zr, LAM["r"]) - cosmo.distmod(zr).value < MR_MIN:
                n_bad += 1
                continue
        sam, pfg = galaxy_dm_samples(
            g.dered_mag_g, g.dered_mag_r, mag_err(g.flux_g, g.flux_ivar_g), mag_err(g.flux_r, g.flux_ivar_r),
            g.theta, zh, z_spec=zsp, z_phot=g.z_phot_mean, z_phot_err=max(g.z_phot_std, 0.01), n=n_mc, **extra)
        m = float(sam.mean())
        tot += sam * (RNG.random(n_mc) < pfg)
        tot_c += m
        if np.isfinite(zsp):
            tot_s += m
        zref = zsp if np.isfinite(zsp) else min(g.z_phot_mean, zh)
        logms = logms_bell03_gr(g.dered_mag_g, g.dered_mag_r, zref, **extra)
        rows.append(dict(frb=frb, ls_id=g.ls_id, ra=g.ra, dec=g.dec, type=g.type, theta_arcsec=g.theta,
                         shape_r=g.shape_r, g=g.dered_mag_g, r=g.dered_mag_r, logms_at_zref=float(logms),
                         z_phot=g.z_phot_mean,
                         z_phot_err=g.z_phot_std, z_spec=zsp, z_src=g.z_best_src, p_fg=pfg,
                         dm_mean=m, dm_p16=float(np.percentile(sam, 16)), dm_p84=float(np.percentile(sam, 84))))
    sga = sga_wide(frb, ra, dec, zh)
    for _, s in sga[sga.theta_deg > SEARCH_DEG].iterrows():
        with np.errstate(divide="ignore", invalid="ignore"):
            gm = 22.5 - 2.5 * np.log10(s.flux_g / s.mw_transmission_g)
            rm = 22.5 - 2.5 * np.log10(s.flux_r / s.mw_transmission_r)
        if not (np.isfinite(gm) and np.isfinite(rm)):
            continue
        logmh = logmh_from_logms(logms_bell03_gr(gm, rm, s.z_leda), s.z_leda)[0]
        if s.theta_deg * 3600 * kpc_per_arcsec(s.z_leda) >= r200_kpc(10 ** logmh, s.z_leda):
            continue
        sam, _ = galaxy_dm_samples(gm, rm, 0.05, 0.05, s.theta_deg * 3600, zh, z_spec=s.z_leda, n=n_mc)
        m = float(sam.mean())
        tot += sam
        tot_c += m
        tot_s += m
        res["dm_fg_sga_wide"] += m
        rows.append(dict(frb=frb, ls_id=f"SGA {s.galaxy}", ra=s.ra, dec=s.dec, type="SGA",
                         theta_arcsec=s.theta_deg * 3600, g=gm, r=rm, z_spec=s.z_leda, z_src="SGA-2020",
                         p_fg=1.0, dm_mean=m, dm_p16=float(np.percentile(sam, 16)),
                         dm_p84=float(np.percentile(sam, 84))))
    res["n_photoz_implausible"] = n_bad
    res["dm_fg"] = float(tot.mean())
    res["dm_fg_mc_sd"] = float(tot.std())
    res["dm_fg_conditional"] = float(tot_c)
    res["dm_fg_specz_only"] = float(tot_s)
    if rows:
        top = max(rows, key=lambda d: d["dm_mean"])
        res["top_gal"], res["top_dm"] = str(top["ls_id"]), top["dm_mean"]
    return res, rows


# ---------------------------------------------------------------- validation
def validate():
    """Notebook cell 26 inputs (FRB 20240114A, SDSS g,i; Taylor+2011): 15.9, 4.6, 0, 104.4, 5.1, 0."""
    g = np.array([22.361, 20.687, 19.238, 16.487, 17.767, 16.431])
    i = np.array([21.7388, 19.4518, 18.1228, 14.6448, 16.4618, 15.7618])
    ge = np.array([0.137, 0.052, 0.019, 0.006, 0.009, 0.007])
    ie = np.array([0.134, 0.040, 0.013, 0.003, 0.006, 0.006])
    zs = np.array([0.13, 0.1296, 0.09, 0.090, 0.09, 0.029])
    th = np.array([0.4042, 35.4207, 87.17, 281.2489, 123.5773, 499.2689])
    out = [np.median(galaxy_dm_samples(g[k], i[k], ge[k], ie[k], th[k], 0.13, z_spec=zs[k],
                                       mass_fn=logms_taylor11_gi)[0]) for k in range(6)]
    return np.round(out, 1)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    val = validate()
    print("validation vs notebook [15.9, 4.6, 0, 104.4, 5.1, 0]:", val)
    hosts = pd.read_csv(HOSTS)
    seed = pd.read_csv(SEED).set_index("frb")
    res, gal = [], []
    for _, row in hosts[hosts.z.notna()].iterrows():
        t0 = time.time()
        r, g = run_host(row, seed.loc[row.frb] if row.frb in seed.index else None)
        res.append(r)
        gal += g
        print(f"{row.frb}: z={row.z:.4f} LS={r['in_ls']} ncand={r['n_cand']} DM_fg={r['dm_fg']:.1f} "
              f"(cond {r['dm_fg_conditional']:.1f}, spec {r['dm_fg_specz_only']:.1f}) [{time.time() - t0:.0f}s]")
    res = pd.DataFrame(res)
    res["dm_fg_err"] = np.hypot(SYS_FRAC * res.dm_fg, res.dm_fg_mc_sd)
    res.to_csv(OUT / "foreground_dm_per_host.csv", index=False)
    pd.DataFrame(gal).to_csv(OUT / "foreground_galaxies.csv", index=False)
    pd.DataFrame({"galaxy": np.arange(1, 7), "notebook": [15.9, 4.6, 0, 104.4, 5.1, 0], "this_code": val}).to_csv(
        OUT / "validation_20240114A.csv", index=False)


def controls(n_per_host=2, offset_deg=0.75, n_mc=500):
    """Random sightlines offset from each LS host at the same z_host: the typical halo term
    that the Macquart DM_cosmic already contains."""
    hosts = pd.read_csv(HOSTS)
    done = pd.read_csv(OUT / "foreground_dm_per_host.csv")
    ls = done[done.in_ls.astype(str).isin(["dr9", "dr10"])].frb
    hosts = hosts[hosts.frb.isin(ls)]
    rng = np.random.default_rng(7)
    out = []
    for _, row in hosts.iterrows():
        for k in range(n_per_host):
            pa = rng.uniform(0, 2 * np.pi)
            dec = row.dec_deg + offset_deg * np.cos(pa)
            ra = row.ra_deg + offset_deg * np.sin(pa) / np.cos(np.radians(dec))
            r = pd.Series(dict(frb=f"{row.frb}_ctrl{k}", ra_deg=ra % 360, dec_deg=dec, z=row.z))
            t0 = time.time()
            res, _ = run_host(r, None, n_mc=n_mc)
            res["host"] = row.frb
            out.append(res)
            print(f"{r.frb}: LS={res['in_ls']} DM_fg={res['dm_fg']:.1f} [{time.time() - t0:.0f}s]", flush=True)
            pd.DataFrame(out).to_csv(OUT / "foreground_dm_controls.csv", index=False)


if __name__ == "__main__":
    import sys
    controls() if "--controls" in sys.argv else main()
