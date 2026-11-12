# Walkthrough of `scripts/foreground_halo_dm.py`

What the code does, block by block, with line numbers. It estimates the dispersion measure (DM) that galaxy halos in front of each FRB host add to the burst, following Bhardwaj, Snelders et al. (2025, arXiv:2506.11915), Appendix B, and the original notebooks (`Estimate_DM_halo.ipynb`, `Galaxy_offset_code.ipynb`, `atel_analysis.ipynb`).

Checks against the notebooks and of the physics are in [CODE_CHECKS.md](CODE_CHECKS.md); a summary is at the end of this file.

## The idea in one paragraph

Every galaxy sits in a dark-matter halo filled with hot ionized gas. If the FRB sightline passes within a foreground galaxy's halo radius $r_{200}$, the burst picks up extra DM. For each FRB, the code finds candidate galaxies near the sightline, estimates each one's stellar mass from its colour, converts that to a halo mass, models the halo gas, and integrates the electron density along the sightline. Redshift and mass uncertainties are propagated with 2,000 Monte Carlo draws per galaxy.

## Setup (lines 1–54)

- **1–20:** docstring: method and usage (`python scripts/foreground_halo_dm.py`, or `--controls` for random sightlines).
- **36–39:** paths. Outputs go to `Publication/materials/foreground_halos/`; catalogue queries are cached in `ls_cache/` so reruns need no internet.
- **41–48:** analysis constants:
  - search radius 0.5 deg;
  - size-ratio cut 0.018;
  - 30 deg radius for very large nearby galaxies;
  - 2,000 draws;
  - 0.10 dex mass systematic;
  - 50% systematic on each host's total;
  - $M_r < -24$ flags an implausible photo-z;
  - fixed random seed 42.
- **50–54:** derived constants:
  - baryon fraction $f_b = \Omega_b/\Omega_m$ (Planck18);
  - proton mass;
  - the halo-mass grid ($\log M_h$ = 9 to 15.5 in 0.005 dex steps);
  - effective wavelengths of the g, r, i, z, W1 bands.

## Halo physics (lines 57–111)

- **58–66, `moster13_logms`:** the Moster et al. (2013) stellar-to-halo mass relation,
  $$M_\star = 2N M_h \left[(M_h/M_1)^{-\beta} + (M_h/M_1)^{\gamma}\right]^{-1},$$
  with $N$, $M_1$, $\beta$, $\gamma$ each evolving linearly in $z/(1+z)$.
- **69–74, `logmh_from_logms`:** inverts that relation. $M_\star$ is evaluated on the halo-mass grid at the galaxy's redshift and interpolated back, which works because the relation is monotonic.
- **77–80, `dutton14_c200c`:** halo concentration from Dutton & Macciò (2014), $\log_{10} c = a(z) + b(z)\log_{10}[M_h/(10^{12}\,h^{-1}M_\odot)]$.
- **83–89, `rho_crit`, `r200_kpc`:** critical density at $z$, and $r_{200} = [3M_h/(4\pi \cdot 200\,\rho_c(z))]^{1/3}$.
- **92–93, `kpc_per_arcsec`:** proper kpc per arcsec, so impact parameters are physical distances.
- **96–111, `dm_mnfw`:** the DM through one halo. It is vectorized, so all 2,000 draws are computed at once.
  1. Impact parameter $b = \theta \times$ (kpc per arcsec).
  2. Central density, as in Eq. B3 of the paper:
     $$\rho_0 = \frac{200\rho_c}{3}\,\frac{c^3}{\ln(1+c) - c/(1+c)}.$$
  3. 801 points along the chord through the halo, from $-\sqrt{r_{200}^2 - b^2}$ to $+\sqrt{r_{200}^2 - b^2}$.
  4. Modified-NFW gas density, Eq. B2:
     $$\rho_{\rm gas} = \frac{f_{\rm gas} f_b \rho_0}{y^{1-\alpha}(y_0+y)^{2+\alpha}}, \qquad y = c\,r/r_{200},$$
     with $\alpha = y_0 = 2$ and $f_{\rm gas} = 0.75$.
  5. Electron density $n_e = 1.1667\,\rho_{\rm gas}/(1.3\,m_p)$, for fully ionized hydrogen and helium.
  6. $\mathrm{DM} = \int n_e\,ds$, with $ds$ in pc. The result is zero if $b \ge r_{200}$, and divided by $(1+z_{\rm gal})$ to give the observer frame.

## Stellar mass (lines 114–155)

- **115–121, `rest_mag`:** rest-frame magnitude at a chosen wavelength. It interpolates the observed g, r, (i), z, (W1) magnitudes in log wavelength at $\lambda_{\rm rest}(1+z)$ and adds $2.5\log_{10}(1+z)$ (the AB K-correction term).
- **124–130, `logms_bell03_gr`:** Bell et al. (2003) colour mass-to-light ratio with a Chabrier IMF, as in `atel_analysis.ipynb`:
  $$\log M_\star = 1.097\,(g-r)_0 - 0.406 - 0.4\,(M_r - 4.65).$$
  The paper used Taylor+11 g−i from SDSS; LS DR9 has no i band.
- **133–136, `logms_taylor11_gi`:** the paper's Taylor et al. (2011) relation. It is used only to reproduce the notebook's validation numbers.
- **139–155, `galaxy_dm_samples`:** the Monte Carlo for one galaxy.
  - **Spectroscopic redshift:** all draws use it, and $P_{\rm fg} = 1$.
  - **Photometric redshift only:** draws come from $\mathcal{N}(z_{\rm phot}, \sigma_z)$ truncated to $0.001 \le z \le z_{\rm host}$, as in the paper. $P_{\rm fg} = \Phi(\text{upper}) - \Phi(\text{lower})$ is the probability that the galaxy is actually in front.
  - **Mass draws:** stellar-mass scatter = photometric error (propagated from the g and r errors) combined in quadrature with 0.10 dex. Each draw is mapped to a halo mass and then to a DM.
  - **Returns:** the 2,000 DM samples and $P_{\rm fg}$.

## Catalogue queries (lines 158–242)

All queries go to NOIRLab's Astro Data Lab and are cached as `ls_cache/<FRB>_<name>.csv`.

- **159–168, `dl_query`:** runs SQL with up to 3 retries.
- **171–177, `cached`:** returns the cached CSV if present, otherwise queries and saves.
- **180–191, `ls_footprint`:** is the field in LS DR9? If not, is it in DR10? If neither, the FRB is skipped.
- **194–207, `ls_candidates`:** the paper's selection, done in SQL. It takes Tractor sources joined to LS photo-z, keeping those that:
  - lie within 0.5 deg;
  - are extended (not PSF or duplicate) and primary;
  - satisfy $R_{\rm eff} > 0.018\,\theta$;
  - satisfy $z_{\rm phot} - 2\sigma_z \le z_{\rm host}$.

  The size cut works because galaxy size is about 0.018 of its halo radius, so a sightline can only be inside $r_{200}$ if $R_{\rm eff}/\theta > 0.018$.
- **210–230, `spectroscopic_redshifts`:** spectroscopic redshifts within 0.5 deg from DESI DR1 (good redshifts, galaxies, primary), SDSS DR17 (good redshifts, galaxies) and SGA-2020 (HyperLEDA).
- **233–242, `sga_wide`:** SGA-2020 galaxies within 30 deg with $0.0005 < z \le z_{\rm host}$. These are nearby galaxies whose halos could reach the sightline from beyond 0.5 deg (the paper did this check with HECATE).

## Cleaning the candidate list (lines 245–302)

- **250–265, `attach_specz`:** best redshift for each candidate. It uses the LS compiled spec-z if there is one; otherwise the nearest external spec-z within the larger of 1.5″ and $0.2\,R_{\rm eff}$.
- **268–278, `identify_host`:** finds the host in the LS catalogue. Among sources within the larger of 5″ and $3R_e$ of the FRB, it takes the one closest in r magnitude to our GALFIT host (if within 1.5 mag), otherwise the nearest.
- **281–287, `remove_host`:** drops the host and its shreds. A shred is any source within the largest of $2R_e$, $2\times$ the LS size and 2″ of the host, unless its spectroscopic redshift differs from the host's by more than 0.01.
- **290–293, `mag_err`:** magnitude error $= 1.0857/(\text{S/N})$, capped at 1 mag, or 0.3 mag if undefined.
- **296–298, `ab_mag`:** SGA flux (nanomaggies) to Galactic-extinction-corrected AB magnitude.
- **301–302, `sample_summary`:** mean, 16th and 84th percentile of a galaxy's DM samples.

## Per FRB (lines 305–398)

- **306–323, `load_candidates`:** query, convert columns to numbers, attach spectroscopic redshifts, remove the host. It then drops galaxies with a spectroscopic redshift behind the host ($z > z_{\rm host} + 0.001$) and any without valid g and r magnitudes.
- **326–398, `run_host`:** the main loop.
  - **330–338:** set up the result row; skip if outside LS; count candidates and spectroscopic redshifts.
  - **340–343:** three running totals:
    - `total`: the Monte Carlo sum, where each draw includes a photo-z galaxy with probability $P_{\rm fg}$;
    - `total_conditional`: the paper's convention, every candidate assumed in front;
    - `total_specz`: spectroscopic galaxies only.
  - **345–370:** for each LS candidate:
    1. clip a spectroscopic $z$ at $z_{\rm host}$;
    2. skip photo-z galaxies that would be brighter than $M_r = -24$ (a wrong photo-z);
    3. draw DM samples;
    4. add them to the three totals;
    5. record one row in `foreground_galaxies.csv`.
  - **372–388:** for each SGA galaxy beyond 0.5 deg: compute its magnitudes, and check whether the sightline is inside $r_{200}$ at its nominal mass. If it is, draw samples (0.05 mag errors, HyperLEDA redshift) and add them to all totals.
  - **390–398:** the per-host summary: mean and spread of the Monte Carlo total, the conditional and spec-z-only sums, and the largest single contributor.

## Driver (lines 401–460)

- **405–415, `validate`:** reruns the notebook's six galaxies for FRB 20240114A (SDSS g, i; Taylor+11 masses). The expected medians are 15.9, 4.6, 0, 104.4, 5.1, 0 pc cm⁻³.
- **418–438, `main`:** validation, then every host with a redshift. It writes:
  - `foreground_dm_per_host.csv`, with $\sigma = \sqrt{(0.5\,\mathrm{DM_{fg}})^2 + \sigma_{\rm MC}^2}$;
  - `foreground_galaxies.csv`;
  - `validation_20240114A.csv`.
- **441–460, `controls`:** two random sightlines 0.75 deg from each host, at the host's redshift and with no host removal. This measures the typical foreground DM that is already part of the average cosmic DM, written to `foreground_dm_controls.csv`.

## How the outputs are used

`robust_correlations.load` (lines 46–61 of `scripts/robust_correlations.py`) builds the corrected host DMs, in the rest frame:

- `dm_host_fg` $= \mathrm{DM_{host}} - \mathrm{DM_{fg}}(1+z_{\rm host})$;
- `dm_host_fg_excess` $= \mathrm{DM_{host}} - (\mathrm{DM_{fg}} - \langle \mathrm{DM_{fg}}\rangle_{\rm random})(1+z_{\rm host})$, which removes only the excess over a random sightline, because the Macquart mean already includes typical halos.

`scripts/foreground_sightline_counts.py` counts intersected halos per sightline (`SIGHTLINES.md`).

## Summary of the checks

Full numbers are in [CODE_CHECKS.md](CODE_CHECKS.md), produced by `scripts/check_foreground_halo_dm.py`.

**Agreement with the original code.** The physics is the same as `Estimate_DM_halo.ipynb`: Moster+13 parameters, the Eq. B3 profile normalization, $f_{\rm gas}=0.75$, $\alpha=y_0=2$, $\mu_H=1.3$ and $n_e/n_H=1.1667$. On an 80-case grid the halo masses agree to 0.0000 dex and the integration to 0.000%. The DMs are 3–6% higher than the notebook's because the notebook passes $M$ in $M_\odot$ to Colossus, which expects $M_\odot/h$, so its concentrations come out about 3.5% low. With that unit fixed the two agree to 0.001%. The same unit slip explains the validation offset for FRB 20240114A (15.9 and 104.4 in the notebook, 16.9 and 109.7 here). The deliberate differences (Bell g−r masses because DR9 has no i band, K-corrections, $P_{\rm fg}$ weighting, SGA instead of HECATE, grid inversion) are listed in the table in CODE_CHECKS.md.

**Science consistency.** Units, cosmology, impact parameters, the observer-frame $1/(1+z)$ factor and the order of the subtraction are all correct. Issues found:

1. **Gas normalization.** Eq. B3 (paper, notebook and this code) normalizes the gas profile with the dark-matter $\rho_0$. The gas mass inside $r_{200}$ is then only $f_b(c)/f_{\rm dm}(c)\approx0.2$–$0.35$ of $f_{\rm gas}f_bM_{200}$, so halo DMs are about $3.8\times$ lower than the gas-mass-normalized profile of the `frb` package (Milky Way: 14 vs 61 pc cm$^{-3}$). This is not clearly wrong. Along random sightlines the coded halos already supply about 17% of the mean $\rm DM_{cosmic}$, and the gas-mass-normalized version would supply about 64% from LS galaxies alone, which is too much. The coded convention therefore behaves like an effective hot-gas fraction of 0.13–0.26, and $f_{\rm hot}=0.75$ is an upper bound. With the upper bound, DM$_{\rm host}$ vs $q$ is still null, and DM–$|{\rm RM}|$ weakens to $r_S=+0.34$ ($p=0.045$).
2. **Mass error coefficient.** The r-band term in $\sigma(\log M_*)$ uses 0.659 instead of $1.097+0.4=1.497$. The median $\sigma$ changes from 0.126 to 0.147 dex, smaller than the 0.10 dex systematic that is already added.
3. **Paper text typos (code is correct).** The intercept of Eq. B1, and the 0.018 size ratio, which is written inverted.
4. **Coverage.** M31 and M33 are not in SGA, but no FRB lies within their $r_{200}$. The 30° SGA search radius never binds (the largest angular $r_{200}$ is 5.5°).
