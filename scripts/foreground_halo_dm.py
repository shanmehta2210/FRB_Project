"""Dispersion measure contributed by foreground (intervening) galaxy halos along each FRB sightline.

Method: Bhardwaj, Snelders et al. (2025, arXiv:2506.11915), Appendix B.

Candidates (Legacy Surveys DR9, or DR10 where DR9 has no coverage):
  * extended Tractor sources within 0.5 deg of the FRB,
  * photometric redshift consistent with the foreground (z_phot - 2 sigma_z <= z_host),
  * angular size R_eff / theta > 0.018,
  * the host galaxy and its shreds removed.
Spectroscopic redshifts come from the LS compilation, DESI DR1, SDSS DR17 and SGA-2020; galaxies behind the
host are dropped. Photo-z-only galaxies are drawn from N(z_phot, sigma_z) truncated at z_host.

Each galaxy: stellar mass from g-r colour (Bell+2003) -> halo mass (Moster+2013) -> concentration
(Dutton & Maccio 2014) -> modified-NFW gas profile (Prochaska & Zheng 2019) -> DM_obs = DM / (1 + z_gal).
SGA-2020 galaxies out to 30 deg are also checked for large halos reaching the sightline from outside 0.5 deg.
A 50% systematic is added to each host total.

Usage: python foreground_halo_dm.py              (all hosts)
       python foreground_halo_dm.py --controls   (random offset sightlines)
"""
from __future__ import annotations

import io
import sys
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

SEARCH_DEG = 0.5          # candidate search radius
SIZE_RATIO_MIN = 0.018    # minimum R_eff / theta
SGA_DEG = 30.0            # search radius for very large nearby galaxies
N_MC = 2000               # Monte Carlo draws per galaxy
MASS_SYS_DEX = 0.10       # systematic scatter on stellar mass
DM_SYS_FRAC = 0.50        # systematic on each host's total
MR_BRIGHTEST = -24.0      # brighter than this at z_phot means the photo-z is implausible
RNG = np.random.default_rng(42)

FB = cosmo.Ob0 / cosmo.Om0
M_P = const.m_p.cgs.value
KPC_CM = u.kpc.to(u.cm)
LOGMH_GRID = np.linspace(9.0, 15.5, 1301)
BAND_WAVELENGTH = {"g": 4810.0, "r": 6440.0, "i": 7840.0, "z": 9260.0, "w1": 33680.0}


# ---------------------------------------------------------------- halo physics
def moster13_logms(logmh, z):
    """Stellar mass from halo mass (Moster+2013)."""
    zf = z / (1 + z)
    n = 0.0351 - 0.0247 * zf
    beta = 1.376 - 0.826 * zf
    gamma = 0.608 + 0.329 * zf
    logm1 = 11.59 + 1.195 * zf
    x = 10 ** (logmh - logm1)
    return logmh + np.log10(2 * n) - np.log10(x ** (-beta) + x ** gamma)


def logmh_from_logms(logms, z):
    """Invert Moster+2013 numerically."""
    logms, z = np.broadcast_arrays(np.atleast_1d(logms), np.atleast_1d(z))
    grid = moster13_logms(LOGMH_GRID[None, :], z.ravel()[:, None])
    out = np.array([np.interp(m, g, LOGMH_GRID) for m, g in zip(logms.ravel(), grid)])
    return out.reshape(logms.shape)


def dutton14_c200c(m_msun, z):
    b = -0.101 + 0.026 * z
    a = 0.520 + (0.905 - 0.520) * np.exp(-0.617 * z ** 1.21)
    return 10 ** (a + b * np.log10(m_msun * cosmo.h / 1e12))


def rho_crit(z):
    return cosmo.critical_density(z).to(u.g / u.cm ** 3).value


def r200_kpc(m_msun, z):
    m_g = m_msun * const.M_sun.cgs.value
    return (3 * m_g / (4 * np.pi * 200 * rho_crit(z))) ** (1 / 3) / KPC_CM


def kpc_per_arcsec(z):
    return cosmo.kpc_proper_per_arcmin(z).to(u.kpc / u.arcsec).value


def dm_mnfw(logmh, z, theta_arcsec, alpha=2.0, y0=2.0, f_gas=0.75, mu_h=1.3, ne_per_nh=1.1667, npts=801):
    """Observer-frame DM of a modified-NFW halo crossed at angular offset theta; zero outside r200."""
    logmh, z, theta = np.broadcast_arrays(*(np.atleast_1d(np.asarray(v, float)) for v in (logmh, z, theta_arcsec)))
    m = 10 ** logmh
    c = dutton14_c200c(m, z)
    r200 = r200_kpc(m, z)
    impact = theta * kpc_per_arcsec(z)
    rho0 = (200 * rho_crit(z) / 3) * c ** 3 / (np.log(1 + c) - c / (1 + c))

    half_chord = np.sqrt(np.clip(r200 ** 2 - impact ** 2, 0, None))
    s = half_chord[:, None] * np.linspace(-1, 1, npts)[None, :]
    y = c[:, None] * np.sqrt(s ** 2 + impact[:, None] ** 2) / r200[:, None]
    rho_gas = f_gas * FB * rho0[:, None] / (y ** (1 - alpha) * (y0 + y) ** (2 + alpha))
    n_e = rho_gas / (mu_h * M_P) * ne_per_nh
    dm = np.trapezoid(n_e, s * 1e3, axis=1)
    return np.where(impact < r200, dm / (1 + z), 0.0)


# ---------------------------------------------------------------- stellar mass
def rest_mag(mags, z, lam_rest):
    """Rest-frame AB magnitude at lam_rest, interpolating the observed SED in log wavelength."""
    bands = [b for b in BAND_WAVELENGTH if np.isfinite(mags.get(b, np.nan))]
    xp = np.log10([BAND_WAVELENGTH[b] for b in bands])
    fp = np.array([mags[b] for b in bands])
    z = np.asarray(z)
    return np.interp(np.log10(lam_rest * (1 + z)), xp, fp) + 2.5 * np.log10(1 + z)


def logms_bell03_gr(g, r, z, imag=np.nan, zmag=np.nan, w1=np.nan):
    """Bell+2003 g-r mass-to-light ratio with the Chabrier zero-point."""
    mags = dict(g=g, r=r, i=imag, z=zmag, w1=w1)
    g0 = rest_mag(mags, z, BAND_WAVELENGTH["g"])
    r0 = rest_mag(mags, z, BAND_WAVELENGTH["r"])
    abs_r = r0 - cosmo.distmod(z).value
    return 1.097 * (g0 - r0) - 0.406 - 0.4 * (abs_r - 4.65)


def logms_taylor11_gi(g, i, z):
    """Taylor+2011 g-i relation; used only to validate against the original notebook."""
    abs_i = i - cosmo.distmod(z).value
    return 1.15 + 0.70 * (g - i) - 0.4 * abs_i


def galaxy_dm_samples(g, r, g_err, r_err, theta, z_host, z_spec=np.nan, z_phot=np.nan, z_phot_err=np.nan,
                      n=N_MC, mass_fn=logms_bell03_gr, **mags):
    """Monte Carlo DM samples for one galaxy, and the probability that it lies in front of the host."""
    if np.isfinite(z_spec):
        zs = np.full(n, z_spec)
        p_fg = 1.0
    else:
        lo, hi = (1e-3 - z_phot) / z_phot_err, (z_host - z_phot) / z_phot_err
        zs = stats.truncnorm.rvs(lo, hi, loc=z_phot, scale=z_phot_err, size=n, random_state=RNG)
        p_fg = float(stats.norm.cdf(hi) - stats.norm.cdf(lo))
    if mass_fn is logms_bell03_gr:
        mass_err = np.hypot(1.097 * g_err, 0.659 * r_err)
    else:
        mass_err = np.hypot(0.70 * g_err, 1.1 * r_err)
    sigma = np.hypot(mass_err, MASS_SYS_DEX)
    logms = mass_fn(g, r, zs, **mags) + RNG.normal(0, sigma, n)
    return dm_mnfw(logmh_from_logms(logms, zs), zs, theta), p_fg


# ---------------------------------------------------------------- catalogue queries (cached)
def dl_query(sql, retries=3):
    from dl import queryClient as qc
    for attempt in range(retries):
        try:
            return pd.read_csv(io.StringIO(qc.query(sql=sql, fmt="csv", timeout=600)))
        except Exception as e:  # noqa: BLE001
            if attempt == retries - 1:
                raise
            print(f"   retry ({e.__class__.__name__}: {str(e)[:80]})")
            time.sleep(5)


def cached(name, frb, fetch):
    f = CACHE / f"{frb}_{name}.csv"
    if f.exists():
        return pd.read_csv(f)
    df = fetch()
    df.to_csv(f, index=False)
    return df


def ls_footprint(frb, ra, dec):
    """'dr9' if the field is in LS DR9, else 'dr10' if in DR10, else None."""
    def fetch():
        for dr in ("dr9", "dr10"):
            n = dl_query(f"""SELECT COUNT(*) AS n FROM ls_{dr}.tractor AS t JOIN ls_{dr}.photo_z AS p
                             ON t.ls_id = p.ls_id WHERE q3c_radial_query(t.ra, t.dec, {ra}, {dec}, 0.02)""").n.iloc[0]
            if n > 0:
                return pd.DataFrame({"dr": [dr]})
        return pd.DataFrame({"dr": [None]})

    dr = cached("footprint", frb, fetch).dr.iloc[0]
    return dr if isinstance(dr, str) else None


def ls_candidates(frb, ra, dec, z_host, dr):
    i_band = "t.dered_mag_i," if dr == "dr10" else ""
    sql = f"""
    SELECT t.ls_id, t.ra, t.dec, t.type, t.shape_r, t.dered_mag_g, t.dered_mag_r, {i_band} t.dered_mag_z, t.dered_mag_w1,
           t.flux_g, t.flux_r, t.flux_ivar_g, t.flux_ivar_r, t.flux_w1, t.flux_ivar_w1,
           p.z_phot_mean, p.z_phot_std, p.z_spec,
           q3c_dist(t.ra, t.dec, {ra}, {dec})*3600 AS theta
    FROM ls_{dr}.tractor AS t JOIN ls_{dr}.photo_z AS p ON t.ls_id = p.ls_id
    WHERE q3c_radial_query(t.ra, t.dec, {ra}, {dec}, {SEARCH_DEG})
      AND t.type <> 'PSF' AND t.type <> 'DUP' AND t.brick_primary = 1
      AND t.shape_r > {SIZE_RATIO_MIN} * q3c_dist(t.ra, t.dec, {ra}, {dec})*3600
      AND p.z_phot_mean - 2*p.z_phot_std <= {z_host}
    """
    return cached("ls", frb, lambda: dl_query(sql))


def spectroscopic_redshifts(frb, ra, dec):
    queries = {
        "DESI DR1": f"""SELECT mean_fiber_ra AS ra, mean_fiber_dec AS dec, z FROM desi_dr1.zpix
            WHERE q3c_radial_query(mean_fiber_ra, mean_fiber_dec, {ra}, {dec}, {SEARCH_DEG})
            AND zwarn = 0 AND spectype = 'GALAXY' AND zcat_primary""",
        "SDSS DR17": f"""SELECT ra, dec, z FROM sdss_dr17.specobj
            WHERE q3c_radial_query(ra, dec, {ra}, {dec}, {SEARCH_DEG}) AND zwarning = 0 AND class = 'GALAXY'""",
        "SGA-2020": f"""SELECT ra, dec, z_leda AS z FROM sga2020.ellipse
            WHERE q3c_radial_query(ra, dec, {ra}, {dec}, {SEARCH_DEG}) AND z_leda > 0""",
    }

    def fetch():
        parts = []
        for src, sql in queries.items():
            try:
                parts.append(dl_query(sql).assign(src=src))
            except Exception as e:  # noqa: BLE001
                print(f"   {src} spec-z query failed: {str(e)[:100]}")
        return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["ra", "dec", "z", "src"])

    return cached("specz", frb, fetch)


def sga_wide(frb, ra, dec, z_host):
    sql = f"""
    SELECT e.sga_id, e.galaxy, e.ra, e.dec, e.z_leda, e.d26,
           t.flux_g, t.flux_r, t.mw_transmission_g, t.mw_transmission_r,
           q3c_dist(e.ra, e.dec, {ra}, {dec}) AS theta_deg
    FROM sga2020.ellipse AS e JOIN sga2020.tractor AS t ON e.sga_id = t.sga_id
    WHERE q3c_radial_query(e.ra, e.dec, {ra}, {dec}, {SGA_DEG})
      AND e.z_leda > 0.0005 AND e.z_leda <= {z_host}
    """
    return cached("sga", frb, lambda: dl_query(sql))


# ---------------------------------------------------------------- candidate cleaning
def separation_arcsec(ra1, dec1, ra2, dec2):
    return np.hypot((ra1 - ra2) * np.cos(np.radians(dec1)), dec1 - dec2) * 3600


def attach_specz(cand, specz, tol_arcsec=1.5):
    """Best redshift per candidate: LS compilation first, then the nearest external spec-z match."""
    cand = cand.copy()
    cand["z_best"], cand["z_best_src"] = np.nan, ""
    has_ls = cand.z_spec.to_numpy(float) > 0
    cand.loc[has_ls, "z_best"] = cand.loc[has_ls, "z_spec"]
    cand.loc[has_ls, "z_best_src"] = "LS photo-z compilation"
    if len(specz) and len(cand):
        sep = separation_arcsec(cand.ra.to_numpy()[:, None], cand.dec.to_numpy()[:, None],
                                specz.ra.to_numpy()[None, :], specz.dec.to_numpy()[None, :])
        nearest = sep.argmin(1)
        tol = np.maximum(tol_arcsec, 0.2 * np.nan_to_num(cand.shape_r.to_numpy(float)))
        hit = (sep[np.arange(len(cand)), nearest] < tol) & ~np.isfinite(cand.z_best.to_numpy(float))
        cand.loc[hit, "z_best"] = specz.z.to_numpy()[nearest[hit]]
        cand.loc[hit, "z_best_src"] = specz.src.to_numpy()[nearest[hit]]
    return cand


def identify_host(cand, m_r_host, re_host):
    """LS source within max(5", 3 R_e) of the FRB closest in r magnitude to our GALFIT host."""
    radius = max(5.0, 3 * re_host) if np.isfinite(re_host) else 5.0
    near = cand[cand.theta < radius]
    if not len(near):
        return None
    if np.isfinite(m_r_host):
        dmag = (near.dered_mag_r - m_r_host).abs()
        if dmag.min() < 1.5:
            return near.loc[dmag.idxmin()]
    return near.loc[near.theta.idxmin()]


def remove_host(cand, host, re_host, z_host):
    """Drop the host and its shreds (sources within 2 R_e of it not known to be at another redshift)."""
    radius = max(2 * (re_host if np.isfinite(re_host) else 1.0), 2 * float(host.shape_r), 2.0)
    sep = separation_arcsec(host.ra, host.dec, cand.ra, cand.dec)
    z = cand.z_best.to_numpy(float)
    shred = (sep < radius) & ~(np.isfinite(z) & (np.abs(z - z_host) > 0.01))
    return cand[~shred & (cand.ls_id != host.ls_id)]


def mag_err(flux, ivar):
    with np.errstate(divide="ignore", invalid="ignore"):
        e = 1.0857 / (flux * np.sqrt(ivar))
    return np.where(np.isfinite(e) & (e > 0), np.minimum(e, 1.0), 0.3)


def ab_mag(flux, mw_transmission):
    with np.errstate(divide="ignore", invalid="ignore"):
        return 22.5 - 2.5 * np.log10(flux / mw_transmission)


def sample_summary(sam):
    return dict(dm_mean=float(sam.mean()), dm_p16=float(np.percentile(sam, 16)), dm_p84=float(np.percentile(sam, 84)))


# ---------------------------------------------------------------- per host
def load_candidates(frb, ra, dec, z_host, dr, seed_row, res):
    cand = ls_candidates(frb, ra, dec, z_host, dr)
    if "dered_mag_i" not in cand:
        cand["dered_mag_i"] = np.nan
    for col in cand.columns.drop(["ls_id", "type"]):
        cand[col] = pd.to_numeric(cand[col], errors="coerce")
    cand = attach_specz(cand, spectroscopic_redshifts(frb, ra, dec))

    if seed_row is not None:
        re_host = float(seed_row.re_arcsec_galfit)
        host = identify_host(cand, float(seed_row.m_r_galfit), re_host)
        if host is not None:
            res["host_ls_id"] = str(host.ls_id)
            cand = remove_host(cand, host, re_host, z_host)

    z = cand.z_best.to_numpy(float)
    cand = cand[~(np.isfinite(z) & (z > z_host + 0.001))]
    return cand[(cand.dered_mag_g > 0) & (cand.dered_mag_r > 0) & np.isfinite(cand.dered_mag_r)]


def run_host(row, seed_row, n_mc=N_MC):
    """Total foreground DM for one sightline. Returns (per-host summary, per-galaxy rows)."""
    frb, ra, dec, z_host = row.frb, row.ra_deg, row.dec_deg, row.z
    res = dict(frb=frb, z=z_host, in_ls=False, n_cand=0, n_specz=0, host_ls_id="", dm_fg=0.0,
               dm_fg_conditional=0.0, dm_fg_specz_only=0.0, dm_fg_sga_wide=0.0, dm_fg_mc_sd=0.0,
               top_gal="", top_dm=0.0)
    dr = ls_footprint(frb, ra, dec)
    if dr is None:
        return res, []
    res["in_ls"] = dr
    cand = load_candidates(frb, ra, dec, z_host, dr, seed_row, res)
    res["n_cand"] = len(cand)
    res["n_specz"] = int(np.isfinite(cand.z_best.to_numpy(float)).sum())

    total = np.zeros(n_mc)       # MC total, each photo-z galaxy included with probability p_fg
    total_conditional = 0.0      # sum of means, assuming every candidate is in front
    total_specz = 0.0            # sum of means, spectroscopic galaxies only
    rows, n_implausible = [], 0

    for _, g in cand.iterrows():
        z_spec = min(float(g.z_best), z_host) if np.isfinite(g.z_best) else np.nan
        z_ref = z_spec if np.isfinite(z_spec) else min(g.z_phot_mean, z_host)
        w1_detected = g.flux_w1 * np.sqrt(max(g.flux_ivar_w1, 0)) > 3
        extra = dict(imag=g.dered_mag_i if g.dered_mag_i > 0 else np.nan,
                     zmag=g.dered_mag_z if g.dered_mag_z > 0 else np.nan,
                     w1=g.dered_mag_w1 if w1_detected else np.nan)
        if not np.isfinite(z_spec):
            mags = dict(g=g.dered_mag_g, r=g.dered_mag_r, i=extra["imag"], z=extra["zmag"], w1=extra["w1"])
            if rest_mag(mags, z_ref, BAND_WAVELENGTH["r"]) - cosmo.distmod(z_ref).value < MR_BRIGHTEST:
                n_implausible += 1
                continue

        sam, p_fg = galaxy_dm_samples(
            g.dered_mag_g, g.dered_mag_r, mag_err(g.flux_g, g.flux_ivar_g), mag_err(g.flux_r, g.flux_ivar_r),
            g.theta, z_host, z_spec=z_spec, z_phot=g.z_phot_mean, z_phot_err=max(g.z_phot_std, 0.01), n=n_mc, **extra)
        mean = float(sam.mean())
        total += sam * (RNG.random(n_mc) < p_fg)
        total_conditional += mean
        if np.isfinite(z_spec):
            total_specz += mean
        logms = logms_bell03_gr(g.dered_mag_g, g.dered_mag_r, z_ref, **extra)
        rows.append(dict(frb=frb, ls_id=g.ls_id, ra=g.ra, dec=g.dec, type=g.type, theta_arcsec=g.theta,
                         shape_r=g.shape_r, g=g.dered_mag_g, r=g.dered_mag_r, logms_at_zref=float(logms),
                         z_phot=g.z_phot_mean, z_phot_err=g.z_phot_std, z_spec=z_spec, z_src=g.z_best_src,
                         p_fg=p_fg, **sample_summary(sam)))

    sga = sga_wide(frb, ra, dec, z_host)
    for _, s in sga[sga.theta_deg > SEARCH_DEG].iterrows():
        g_mag, r_mag = ab_mag(s.flux_g, s.mw_transmission_g), ab_mag(s.flux_r, s.mw_transmission_r)
        if not (np.isfinite(g_mag) and np.isfinite(r_mag)):
            continue
        theta = s.theta_deg * 3600
        logmh = logmh_from_logms(logms_bell03_gr(g_mag, r_mag, s.z_leda), s.z_leda)[0]
        if theta * kpc_per_arcsec(s.z_leda) >= r200_kpc(10 ** logmh, s.z_leda):
            continue
        sam, _ = galaxy_dm_samples(g_mag, r_mag, 0.05, 0.05, theta, z_host, z_spec=s.z_leda, n=n_mc)
        mean = float(sam.mean())
        total += sam
        total_conditional += mean
        total_specz += mean
        res["dm_fg_sga_wide"] += mean
        rows.append(dict(frb=frb, ls_id=f"SGA {s.galaxy}", ra=s.ra, dec=s.dec, type="SGA", theta_arcsec=theta,
                         g=g_mag, r=r_mag, z_spec=s.z_leda, z_src="SGA-2020", p_fg=1.0, **sample_summary(sam)))

    res["n_photoz_implausible"] = n_implausible
    res["dm_fg"] = float(total.mean())
    res["dm_fg_mc_sd"] = float(total.std())
    res["dm_fg_conditional"] = float(total_conditional)
    res["dm_fg_specz_only"] = float(total_specz)
    if rows:
        top = max(rows, key=lambda d: d["dm_mean"])
        res["top_gal"], res["top_dm"] = str(top["ls_id"]), top["dm_mean"]
    return res, rows


# ---------------------------------------------------------------- driver
NOTEBOOK_VALIDATION = [15.9, 4.6, 0, 104.4, 5.1, 0]


def validate():
    """Reproduce the notebook's six galaxies for FRB 20240114A (SDSS g, i; Taylor+2011 masses)."""
    g = np.array([22.361, 20.687, 19.238, 16.487, 17.767, 16.431])
    i = np.array([21.7388, 19.4518, 18.1228, 14.6448, 16.4618, 15.7618])
    g_err = np.array([0.137, 0.052, 0.019, 0.006, 0.009, 0.007])
    i_err = np.array([0.134, 0.040, 0.013, 0.003, 0.006, 0.006])
    z = np.array([0.13, 0.1296, 0.09, 0.090, 0.09, 0.029])
    theta = np.array([0.4042, 35.4207, 87.17, 281.2489, 123.5773, 499.2689])
    out = [np.median(galaxy_dm_samples(g[k], i[k], g_err[k], i_err[k], theta[k], 0.13, z_spec=z[k],
                                       mass_fn=logms_taylor11_gi)[0]) for k in range(6)]
    return np.round(out, 1)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    val = validate()
    print(f"validation vs notebook {NOTEBOOK_VALIDATION}:", val)
    hosts = pd.read_csv(HOSTS)
    seed = pd.read_csv(SEED).set_index("frb")
    summaries, galaxies = [], []
    for _, row in hosts[hosts.z.notna()].iterrows():
        t0 = time.time()
        res, rows = run_host(row, seed.loc[row.frb] if row.frb in seed.index else None)
        summaries.append(res)
        galaxies += rows
        print(f"{row.frb}: z={row.z:.4f} LS={res['in_ls']} ncand={res['n_cand']} DM_fg={res['dm_fg']:.1f} "
              f"(cond {res['dm_fg_conditional']:.1f}, spec {res['dm_fg_specz_only']:.1f}) [{time.time() - t0:.0f}s]")
    summaries = pd.DataFrame(summaries)
    summaries["dm_fg_err"] = np.hypot(DM_SYS_FRAC * summaries.dm_fg, summaries.dm_fg_mc_sd)
    summaries.to_csv(OUT / "foreground_dm_per_host.csv", index=False)
    pd.DataFrame(galaxies).to_csv(OUT / "foreground_galaxies.csv", index=False)
    pd.DataFrame({"galaxy": np.arange(1, 7), "notebook": NOTEBOOK_VALIDATION, "this_code": val}).to_csv(
        OUT / "validation_20240114A.csv", index=False)


def controls(n_per_host=2, offset_deg=0.75, n_mc=500):
    """Random sightlines 0.75 deg from each host at the same z: the typical halo DM already in DM_cosmic."""
    hosts = pd.read_csv(HOSTS)
    done = pd.read_csv(OUT / "foreground_dm_per_host.csv")
    hosts = hosts[hosts.frb.isin(done[done.in_ls.astype(str).isin(["dr9", "dr10"])].frb)]
    rng = np.random.default_rng(7)
    out = []
    for _, row in hosts.iterrows():
        for k in range(n_per_host):
            pa = rng.uniform(0, 2 * np.pi)
            dec = row.dec_deg + offset_deg * np.cos(pa)
            ra = row.ra_deg + offset_deg * np.sin(pa) / np.cos(np.radians(dec))
            sightline = pd.Series(dict(frb=f"{row.frb}_ctrl{k}", ra_deg=ra % 360, dec_deg=dec, z=row.z))
            t0 = time.time()
            res, _ = run_host(sightline, None, n_mc=n_mc)
            res["host"] = row.frb
            out.append(res)
            print(f"{sightline.frb}: LS={res['in_ls']} DM_fg={res['dm_fg']:.1f} [{time.time() - t0:.0f}s]", flush=True)
            pd.DataFrame(out).to_csv(OUT / "foreground_dm_controls.csv", index=False)


if __name__ == "__main__":
    controls() if "--controls" in sys.argv else main()
