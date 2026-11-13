# Walkthrough of `scripts/foreground_halo_dm.py`

What the code does, block by block with line numbers, and why each step is done that way. It estimates the dispersion measure (DM) that galaxy halos in front of each FRB host add to the burst. It follows Bhardwaj, Snelders et al. (2025, arXiv:2506.11915), Appendix B, and the original notebooks `Estimate_DM_halo.ipynb`, `Galaxy_offset_code.ipynb` and `atel_analysis.ipynb`, with two corrections found by the audit: the halo mass of each galaxy and the treatment of groups. Both versions are computed from the same draws; the published one is kept as `dm_fg_published`.

Checks against the notebooks and of the physics are in [CODE_CHECKS.md](CODE_CHECKS.md); the audit is in [audit/AUDIT.md](audit/AUDIT.md); the test of the halo-mass step on the paper's own galaxies is in [mass_inversion/MASS_INVERSION.md](mass_inversion/MASS_INVERSION.md). A summary is at the end of this file. The outputs of the published method before the corrections are archived in `Publication/materials/archive/foreground_published_method/`.

## The idea

The observed DM of an FRB is the sum of several pieces: the Milky Way, the intergalactic medium (DM$_{\rm cosmic}$), the host galaxy, and any galaxy halos the sightline happens to pass through on the way. We want DM$_{\rm host}$, which is normally obtained by subtracting the Milky Way term and the mean cosmic term (the Macquart relation) from the observed DM. Anything left over, including the gas of intervening halos, is wrongly assigned to the host. If such halos are common, they add noise to DM$_{\rm host}$ that could hide or fake a correlation with host inclination. This code estimates that halo contribution for each sightline so it can be removed.

Every galaxy sits in a dark-matter halo that contains hot ionized gas (the circumgalactic medium, CGM). If the FRB sightline passes within a foreground galaxy's halo radius $r_{200}$, the burst picks up extra DM. For each FRB the code:

1. finds galaxies near the sightline that could be in front of the host and close enough for their halo to reach it;
2. estimates each galaxy's stellar mass from its colour and brightness;
3. assigns a halo mass, then a halo size and concentration;
4. models the halo gas density and integrates the electron density along the sightline;
5. repeats this 2,000 times per galaxy with different masses and redshifts drawn from their uncertainties, so the result carries an error bar;
6. counts each group of galaxies once, and leaves out galaxies that belong to the host's own halo.

Steps 2–4 are a chain of empirical relations, each with scatter, and the DM is very non-linear in mass (it drops to zero when $r_{200}$ no longer reaches the sightline). That is why the uncertainties are propagated by Monte Carlo rather than analytically.

## Where each step is in the original notebooks

| Step | This code | Original |
|---|---|---|
| Stellar mass → halo mass | `conditional_logmh` (corrected); `moster13_logms`, `logmh_from_logms` (published) | `stellarmass_from_halomass`, `halomass_from_stellarmass` (`Estimate_DM_halo.ipynb`, cell 26) |
| Concentration | `dutton14_c200c` | `concentration_index` (Colossus `dutton14`, cell 26) |
| $r_{200}$ | `r200_kpc` | `virial_radius` (cell 26) |
| Gas profile and DM through one halo | `dm_mnfw` | `ModifiedNFW` (`_setup`, `rho_gas`, `n_e`, `compute_DM`), cell 26 |
| Stellar mass from colour | `logms_bell03_gr`; `logms_taylor11_gi` for validation | Bell+03 g−r in `atel_analysis.ipynb`; Taylor+11 g−i in cell 26 |
| Monte Carlo over mass (and photo-z) | `galaxy_dm_samples` | main loop of cell 26; truncated-Gaussian photo-z from the Appendix B text |
| Candidate selection | `ls_candidates` | `process_desi_bricks`, `check_intersecting_halo` and the photo-z cut (`Galaxy_offset_code.ipynb`) |
| Large nearby galaxies | `sga_wide` | `search_hecate_galaxies` (`Galaxy_offset_code.ipynb`) |
| Groups and host environment | `flag_groups` | none (every candidate summed) |
| Random-sightline expectation | `controls`, `calibrate` | none |
| Validation | `validate` | cell 26 inputs and output for FRB 20240114A |

## Setup (lines 1–76)

- **1–31, docstring:** a summary of the method, the two estimates and usage. Run `python scripts/foreground_halo_dm.py` for all hosts, `--controls` for the random comparison sightlines (followed by the calibration), or `--calibrate` for the calibration alone.
- **53–57, paths:** inputs are our host table (`host_props.csv`: FRB positions, redshifts, DM$_{\rm host}$, host stellar mass) and our GALFIT fits (`seed_56.csv`: host magnitude and size, used to find and remove the host). Outputs go to `Publication/materials/foreground_halos/`. Every catalogue query is cached in `ls_cache/`, so reruns are reproducible and need no internet.
- **59–70, analysis constants:**
  - `SEARCH_DEG = 0.5`: search radius, as in the paper. Beyond 0.5° only very nearby galaxies can still reach the sightline, and those are handled separately (`SGA_DEG`).
  - `SIZE_RATIO_MIN = 0.018`: the size cut explained under `ls_candidates`.
  - `SGA_DEG = 30`: radius for the large-nearby-galaxy check, as in the paper (which used HECATE).
  - `N_MC = 2000`: draws per galaxy. The Monte Carlo error on a galaxy's mean DM is then about $1/\sqrt{2000} \approx 2\%$ of its spread, negligible next to the 50% systematic.
  - `MASS_SYS_DEX = 0.10`: systematic scatter added to every stellar mass, as in cell 26. It covers the intrinsic scatter of colour–mass relations, which photometric errors alone miss.
  - `DM_SYS_FRAC = 0.50`: the paper's 50% systematic on each total. It covers the choices of colour–mass relation, stellar-to-halo mass relation, concentration and gas profile, none of which the Monte Carlo varies.
  - `MR_BRIGHTEST = -24`: photo-z sanity limit (see `run_host`).
  - `SHMR_SCATTER = 0.15`: scatter in $\log M_\star$ at fixed halo mass, the value Moster et al. (2013) adopt.
  - `DZ_SPEC = 0.005`: redshift tolerance for group membership between spectroscopic redshifts, $(1+z) \times 1500$ km s$^{-1}$, a cluster-scale velocity spread.
  - `LOGM_MIN`, `ZGRID`: lowest halo mass and redshift grid of the halo-model expectation.
  - `RNG` seed 42: a fixed seed makes every run give identical numbers.
- **72–76, derived constants:**
  - $f_b = \Omega_b/\Omega_m \approx 0.158$ (astropy Planck18), the cosmic baryon fraction: the fraction of a halo's mass expected to be gas plus stars.
  - The proton mass, and the kpc-to-cm conversion.
  - `LOGMH_GRID`: halo masses from $10^9$ to $10^{15.5}\,M_\odot$ in 0.005 dex steps, used to invert the stellar-to-halo relation (published method).
  - `BAND_WAVELENGTH`: effective wavelengths of the g, r, i, z and W1 bands, used for the K-correction.

## Halo physics (lines 79–165)

- **80–88, `moster13_logms`:** the Moster et al. (2013) stellar-to-halo mass relation,
  $$M_\star = 2N M_h \left[(M_h/M_1)^{-\beta} + (M_h/M_1)^{\gamma}\right]^{-1},$$
  with $N$, $M_1$, $\beta$ and $\gamma$ each evolving linearly in $z/(1+z)$. The parameters are identical to cell 26.
  - **Why:** we observe stars but need the halo, because the halo mass sets how far the gas extends ($r_{200}$) and how much gas there is. The relation comes from abundance matching and gives the *mean stellar mass of the galaxies in halos of mass* $M_h$, with 0.15 dex scatter around it.
- **91–96, `logmh_from_logms` (published method):** inverts the relation on the halo-mass grid. This gives the halo whose average galaxy has stellar mass $M_\star$, not the typical halo of a galaxy of that stellar mass.
- **99–128, `conditional_logmh_table`, `conditional_logmh` (corrected):** the halo-mass distribution of galaxies of stellar mass $M_\star$,
  $$P(M_h \mid M_\star) \propto \frac{dn}{d\ln M_h}\,\mathcal{N}\!\left(\log M_\star;\ \log M_\star^{\rm Moster}(M_h),\ 0.15\right),$$
  with the Tinker et al. (2008) halo mass function (Colossus, $M_{200c}$). Its mean and spread are tabulated once on a grid in $z$ (0–1) and $\log M_\star$ (7–12.8); each draw takes a Gaussian with that mean and spread.
  - **Why:** above the peak of the relation, $M_\star \propto M_h^{0.37}$, so 0.15 dex of scatter in $M_\star$ is about 0.4 dex in $M_h$, and the halo mass function falls steeply there. Most galaxies of a given high stellar mass therefore live in lower-mass halos that scattered up. Ignoring this, the inverse runs 0.3 dex high at $\log M_\star = 11$ and 0.8 dex at 11.8, up to $10^{15.5}\,M_\odot$, which would be the halo of a rich cluster. Below $\log M_\star \approx 10.5$ the two agree. This is the known difference between $\langle \log M_\star\rangle(M_h)$ and $\langle M_h\rangle(M_\star)$ (Behroozi et al. 2010; Dutton et al. 2010). `mass_inversion/` shows it with the notebook's own functions: for Galaxy 6 of FRB 20240114A the inverse gives $10^{15.3}\,M_\odot$, 23 times the measured mass of the cluster the galaxy belongs to.
- **131–134, `dutton14_c200c`:** halo concentration from Dutton & Macciò (2014),
  $$\log_{10} c = a + b\log_{10}\left[\frac{M_h}{10^{12}\,h^{-1}M_\odot}\right],\quad a = 0.520 + 0.385\,e^{-0.617 z^{1.21}},\quad b = -0.101 + 0.026z.$$
  - **Why:** the concentration $c = r_{200}/r_s$ sets how centrally peaked the halo is. A more concentrated halo puts more gas near the centre, which raises the DM for sightlines passing close to the galaxy. Lower-mass halos are more concentrated, because they formed earlier, when the Universe was denser.
  - This is the same relation that Colossus `dutton14` implements and the notebook calls. The relation is defined with mass in $h^{-1}M_\odot$, hence the factor `cosmo.h`. The notebook passes $M_\odot$ instead (see the summary).
- **137–143, `rho_crit`, `r200_kpc`:** the critical density $\rho_c(z)$ and
  $$r_{200} = \left[\frac{3M_h}{4\pi \cdot 200\,\rho_c(z)}\right]^{1/3},$$
  the radius inside which the mean density is 200 times critical.
  - **Why:** $r_{200}$ is the conventional halo boundary, and the paper (like the notebook) assumes the hot gas extends to $r_{200}$ and no further. A sightline passing outside $r_{200}$ contributes nothing. This is a convention; real gas extends somewhat beyond.
- **146–147, `kpc_per_arcsec`:** proper (physical) kpc per arcsec at redshift $z$. A halo is a bound object that does not expand with the Universe, and $r_{200}$ is a physical radius, so the angular offset has to become a physical impact parameter.
- **150–165, `dm_mnfw`:** the DM through one halo, vectorized over all 2,000 draws.
  1. **Halo properties per draw:** $c$, $r_{200}$, and the impact parameter $b = \theta \times$ (kpc per arcsec) at the draw's redshift.
  2. **Central density,** as in Eq. B3 of the paper and in cell 26:
     $$\rho_0 = \frac{200\rho_c}{3}\,\frac{c^3}{\ln(1+c) - c/(1+c)}.$$
     This is the normalization of an ordinary NFW dark-matter halo with mass $M_{200}$. The checks discuss its consequence: the gas inside $r_{200}$ adds up to less than $f_{\rm gas} f_b M_{200}$.
  3. **Sightline sampling:** 801 points along the chord through the sphere of radius $r_{200}$. The odd number puts one point exactly at closest approach, where the density peaks. The trapezoid rule agrees with adaptive integration to 0.000%.
  4. **Gas density,** the modified NFW profile of Prochaska & Zheng (2019), paper Eq. B2:
     $$\rho_{\rm gas} = \frac{f_{\rm gas} f_b \rho_0}{y^{1-\alpha}(y_0+y)^{2+\alpha}}, \qquad y = c\,r/r_{200},$$
     with $\alpha = y_0 = 2$ and $f_{\rm gas} = 0.75$. A pure NFW cusp would give unrealistically large DMs near the galaxy; feedback and pressure flatten the gas core.
  5. **Electron density:** $n_e = 1.1667\,\rho_{\rm gas}/(1.3\,m_p)$: hydrogen density from the gas mass per hydrogen atom, times free electrons per hydrogen atom for fully ionized hydrogen plus helium.
  6. **DM:** $\mathrm{DM} = \int n_e\,ds$, with $ds$ in pc. Draws with $b \ge r_{200}$ give exactly zero, as in the notebook.
  7. **Observer frame:** divided by $(1+z_{\rm gal})$, so contributions at different redshifts can be added.

## Stellar mass (lines 168–211)

- **169–175, `rest_mag`:** rest-frame AB magnitude at a chosen wavelength, interpolating the observed magnitudes (g, r, plus i, z and W1 where available) linearly in $\log\lambda$ at $\lambda_{\rm rest}(1+z)$, plus $2.5\log_{10}(1+z)$.
  - **Why (the K-correction):** colour–mass relations are calibrated on rest-frame colours. Observed colours of a $z = 0.3$ galaxy would make it look redder and brighter, and its mass too high. The notebook used observed colours, which is fine for its $z \approx 0.1$ galaxies; our hosts reach $z \approx 1$, and the correction lowers $\log M_\star$ by a median 0.29 dex.
- **178–184, `logms_bell03_gr`:** Bell et al. (2003) colour–mass-to-light relation in g−r,
  $$\log M_\star = 1.097\,(g-r)_0 - 0.406 - 0.4\,(M_r - 4.65),$$
  as used in `atel_analysis.ipynb`. −0.406 is Bell's zero-point shifted by −0.1 dex to a Chabrier IMF. The paper used Taylor et al. (2011) in g−i with SDSS photometry; Legacy Surveys DR9 has no i band.
- **187–190, `logms_taylor11_gi`:** the paper's Taylor et al. (2011) relation, $\log M_\star = 1.15 + 0.70(g-i) - 0.4M_i$, used only by `validate`.
- **193–211, `galaxy_dm_samples`:** the Monte Carlo for one galaxy.
  - **Spectroscopic redshift:** all draws use it, and the galaxy is certainly in front ($P_{\rm fg} = 1$).
  - **Photometric redshift only:** draws come from $\mathcal{N}(z_{\rm phot}, \sigma_z)$ truncated to $0.001 \le z \le z_{\rm host}$, as in the paper. $P_{\rm fg} = \Phi(\text{upper}) - \Phi(\text{lower})$ is the probability that the galaxy really is in front. The paper implicitly sets it to 1; over many high-redshift hosts that would let probable background galaxies count as if they sat just in front of the host.
  - **Mass draws:** Gaussian in $\log M_\star$ with $\sigma = \sqrt{(1.097\,\sigma_g)^2 + (1.497\,\sigma_r)^2 + 0.10^2}$, the derivatives of the Bell relation plus the systematic.
  - **Two halo masses per draw:** a draw from $P(M_h \mid M_\star)$ (corrected) and the inverse (published), each passed through `dm_mnfw`. Both use the same stellar-mass and redshift draws, so their difference is the halo-mass step alone.
  - **Returns:** `{"corrected": ..., "published": ...}`, 2,000 DM samples each (observer frame), and $P_{\rm fg}$.

## Catalogue queries (lines 214–298)

All queries go to NOIRLab's Astro Data Lab and are cached as `ls_cache/<FRB>_<name>.csv`.

- **215–224, `dl_query`:** runs one SQL query, retrying up to 3 times.
- **227–233, `cached`:** returns the cached CSV if present, otherwise queries and saves, so results depend only on files in the repository.
- **236–247, `ls_footprint`:** checks whether the field is in LS DR9, else DR10. If neither covers the field, the FRB is skipped. This is why 49 of our 55 hosts have an estimate.
- **250–263, `ls_candidates`:** the paper's candidate selection, in SQL. It keeps sources within 0.5° that are extended (`type` not `PSF` or `DUP`), `brick_primary` (one copy of sources in overlapping sky tiles), with $R_{\rm eff} > 0.018\,\theta$ and $z_{\rm phot} - 2\sigma_z \le z_{\rm host}$.

  **Why the size cut works:** galaxy size scales with halo size, $R_{\rm eff} \approx 0.018\,r_{200}$ (Somerville et al. 2018). The sightline is inside the halo when $\theta < r_{200}$ in angle, i.e. $R_{\rm eff}/\theta > 0.018$. Both are angles, so the cut needs no redshift. The exact test ($b < r_{200}$ for each draw) happens later in `dm_mnfw`. (The paper writes the ratio inverted; the notebook and this code use $R_{\rm eff}/\theta$.)
- **266–286, `spectroscopic_redshifts`:** spectroscopic redshifts within 0.5° from DESI DR1, SDSS DR17 and SGA-2020 (HyperLEDA), with quality flags. A spectroscopic redshift fixes the distance, the mass and whether the galaxy is in front.
- **289–298, `sga_wide`:** Siena Galaxy Atlas galaxies within 30° with $0.0005 < z \le z_{\rm host}$. A very nearby galaxy's halo can span degrees and reach the sightline from far outside 0.5°. The paper did this check with HECATE.

## Cleaning the candidate list (lines 301–391)

- **302–303, `separation_arcsec`:** flat-sky angular separation.
- **306–321, `attach_specz`:** best redshift per candidate: the LS compiled spectroscopic redshift, else the nearest external one within the larger of 1.5″ and $0.2\,R_{\rm eff}$.
- **324–334, `identify_host`:** finds the FRB host in the LS catalogue: among sources within the larger of 5″ and $3R_e$ of the FRB, the one closest in r magnitude to our GALFIT host (within 1.5 mag), otherwise the nearest.
- **337–343, `remove_host`:** drops the host and its "shreds" (sources within the largest of $2R_e$, twice the LS size and 2″ of the host, unless their spectroscopic redshift differs by more than 0.01). Tractor often splits a large galaxy into several sources, each of which would otherwise get its own halo right on the sightline.
- **346–354, `mag_err`, `ab_mag`:** magnitude errors from flux and inverse variance (capped at 1 mag, 0.3 mag where undefined), and SGA fluxes to dust-corrected AB magnitudes.
- **357–359, `sample_summary`:** mean and 16–84% of the corrected DM samples, and the mean of the published ones.
- **362–363, `same_z`:** two galaxies are at a consistent redshift if they differ by less than their combined redshift error, or by less than `DZ_SPEC`$\,(1+z)$ when both are spectroscopic.
- **366–384, `flag_groups`:** marks galaxies whose halo is not added to the total.
  - **Satellites:** galaxies are visited from the most to the least massive halo. A galaxy that lies within $r_{200}$ of a more massive candidate already kept, at a consistent redshift, is a satellite.
  - **Why:** the members of a group or cluster share one halo, whose mass the halo mass function counts once. Giving every member its own Moster halo adds the same gas several times. The group's halo is represented by its most massive member, the central.
  - **Host environment:** a galaxy within the host's own $r_{200}$ (from the host stellar mass in `host_props.csv`), at a redshift consistent with the host's, belongs to the host system. Its gas is part of what DM$_{\rm host}$ measures, not foreground.
- **387–391, `galaxy_row`:** one row per galaxy for `foreground_galaxies.csv`: position, redshift and its error, stellar mass, conditional and inverse halo masses, $r_{200}$, $P_{\rm fg}$, and both DM summaries.

## Per FRB (lines 394–492)

- **395–412, `load_candidates`:** queries the candidates, converts columns to numbers, attaches spectroscopic redshifts and removes the host. It then drops galaxies whose spectroscopic redshift puts them behind the host ($z > z_{\rm host} + 0.001$), and any without valid g and r magnitudes. The 0.001 margin keeps group members at essentially the host redshift; they are placed at $z_{\rm host}$.
- **415–492, `run_host`:** the main loop for one sightline.
  - **418–428:** set up the result row; skip the FRB if it is outside the LS footprint; count candidates and spectroscopic redshifts.
  - **431–452:** for each LS candidate:
    1. Cap a spectroscopic $z$ at $z_{\rm host}$.
    2. Skip photo-z galaxies that would be brighter than $M_r = -24$ if placed at $\min(z_{\rm phot}, z_{\rm host})$. Even the brightest cluster galaxies are around $-23.5$, so the photo-z is wrong (usually a star or an AGN); kept, it would get an enormous halo on the sightline. `n_photoz_implausible` counts them.
    3. Draw DM samples ($\sigma_z$ floored at 0.01) and record the galaxy's row.
  - **454–466:** for each SGA galaxy beyond 0.5°, test whether its inverse-mass halo reaches the sightline (an upper bound on its size, so no galaxy is missed), and if so draw samples (0.05 mag errors, HyperLEDA redshift).
  - **471–479:** flag satellites and host environment, then build two Monte Carlo totals. In each draw a photo-z galaxy is included with probability $P_{\rm fg}$, so each total's mean is the expected foreground DM and its spread includes the uncertainty about which galaxies are in front.
    - `dm_fg`: corrected halo masses, counted galaxies only;
    - `dm_fg_published`: inverse halo masses, every galaxy (the published method).
  - **480–492:** the per-host summary: both totals and their spreads; `dm_fg_conditional` (every counted candidate assumed in front, the paper's photo-z convention) and `dm_fg_specz_only`; the DM and number of satellites and host-environment galaxies left out; and the largest single contributor (`top_gal`, `top_dm`).

## Driver (lines 495–556)

- **499–509, `validate`:** reruns the notebook's six galaxies for FRB 20240114A with the notebook's SDSS g and i photometry, spectroscopic redshifts, offsets and Taylor+11 masses, using the published (inverse) halo masses. The notebook printed medians of 15.9, 4.6, 0, 104.4, 5.1 and 0 pc cm⁻³; this code gives 16.9, 4.7, 0, 110.3, 5.4 and 0. The difference is the units passed to Colossus in the notebook (see the summary).
- **512–534, `main`:** runs the validation, then every host with a redshift. It writes:
  - `foreground_dm_per_host.csv`, with total uncertainty $\sigma = \sqrt{(0.5\,\mathrm{DM_{fg}})^2 + \sigma_{\rm MC}^2}$ for both methods;
  - `foreground_galaxies.csv`, one row per galaxy, with `satellite`, `hostenv` and `counted` flags;
  - `validation_20240114A.csv`.
- **537–556, `controls`:** the same calculation along two random sightlines 0.75° from each host, at the host's redshift, with no host removal and 500 draws per galaxy. Output: `foreground_dm_controls.csv`.

## Random-sightline expectation (lines 559–619)

The Macquart mean DM$_{\rm cosmic}$ is an average over all sightlines and already includes the gas in galaxy halos. Only the excess of a sightline's foreground over that average may be removed from DM$_{\rm host}$.

- **560–563, `gas_fraction`:** gas mass inside $r_{200}$ of the coded profile, in units of $f_{\rm gas} f_b M_{200}$: $I(c)/m(c)$ with $I = \int_0^c y^3/(2+y)^4\,dy$.
- **566–572, `halo_electron_fraction`:** the mean electron density in halo gas (every halo from $10^{10}$ to $10^{15.5}\,M_\odot$, Tinker+08) over the mean diffuse electron density of the Macquart relation. It is 7–9%.
- **575–583, `halo_model_mean`:** the mean observer-frame DM that halo gas adds out to redshift $z$, integrating that fraction with the Macquart weighting. This is the most halos can supply with this gas profile.
- **586–609, `calibrate`:** the ratio of the control sightlines to `halo_model_mean`, per method, with a bootstrap interval (`foreground_expectation.csv`). `dm_fg_expected` $=$ ratio $\times$ `halo_model_mean`$(z_{\rm host})$ is the smooth random-sightline expectation for each host; `dm_fg_expected_published` is the same for the published method. The corrected method gives 1.01 (0.91–1.12) times the halo model, i.e. it saturates the physical bound; the published method gives 2.24 (2.00–2.47) times, more than the halos can hold.

## How the outputs are used

`robust_correlations.load` builds the corrected host DMs in the host's rest frame:

- `dm_host_fg_excess` $= \mathrm{DM_{host}} - (\mathrm{DM_{fg}} - \mathrm{DM_{fg,expected}})(1+z_{\rm host})$: only the excess over random sightlines removed. This is the consistent correction and the main one.
- `dm_host_fg` $= \mathrm{DM_{host}} - \mathrm{DM_{fg}}(1+z_{\rm host})$: the full foreground removed. It counts the mean halo gas twice and is kept as a bound.
- `dm_host_fg_published`, `dm_host_fg_excess_published`: the same with the published method, for comparison.

The factor $(1+z_{\rm host})$ converts the observer-frame foreground DM to the host rest frame. `scripts/foreground_sightline_counts.py` counts the halos intersected per sightline (`SIGHTLINES.md`); `scripts/foreground_audit.py` checks the estimate and the DM budget (`audit/`).

## Summary of the checks

**Halo masses (`mass_inversion/MASS_INVERSION.md`, `scripts/mass_inversion_test.py`).** Run with the notebook's own functions on the paper's own galaxies, changing only the halo-mass line. The inverse gives Galaxy 6 of FRB 20240114A $10^{15.3}\,M_\odot$, 23 times the measured mass of the cluster it belongs to. With $P(M_h \mid M_\star)$ its DM falls from 104 to 61 pc cm$^{-3}$. Only zero scatter reproduces 104.

**Audit (`audit/AUDIT.md`, `scripts/foreground_audit.py`).** Along random sightlines the published method returns 2.2 times the mean halo DM that the halo mass function allows for the same gas profile. Three causes:
- **Halo masses of massive galaxies,** as above.
- **Group galaxies:** every member of a group gets its own halo.
- **Host environment:** host-system galaxies are counted as foreground.

With all three corrected, as this code now does, the random sightlines sit on the halo-model value (ratio 1.01). The audit also shows that the full subtraction (`dm_host_fg`) double counts the halo gas already in $\langle{\rm DM_{cosmic}}\rangle$, and that no FRB's budget is short by more than $1\sigma$.

Full numbers are in [CODE_CHECKS.md](CODE_CHECKS.md), produced by `scripts/check_foreground_halo_dm.py`.

**Agreement with the original code.** The physics is the same as `Estimate_DM_halo.ipynb` cell 26: Moster+13 parameters, the Eq. B3 profile normalization, $f_{\rm gas}=0.75$, $\alpha=y_0=2$, $\mu_H=1.3$, $n_e/n_H=1.1667$, zero DM outside $r_{200}$ and the $1/(1+z)$ factor. The candidate selection is the same as `Galaxy_offset_code.ipynb` and the Appendix B text.

- **End-to-end test on the notebook's own inputs,** with the published halo masses: the notebook re-run reproduces its medians exactly (15.9, 4.6, 0, 104.4, 5.1, 0 pc cm⁻³). This code gives about 6% more, entirely because the notebook passes the halo mass to Colossus in $M_\odot$ where Colossus expects $h^{-1}M_\odot$ (concentrations about 3.5% low). With the unit corrected, the notebook gives this code's numbers.
- **Component checks.** On an 80-case grid, inverse halo masses agree to 0.0000 dex and the integration to 0.000%.
- **Deliberate differences.** Halo masses from $P(M_h \mid M_\star)$, one halo per group and host environment excluded (the published method is kept alongside); Bell g−r masses (DR9 has no i band); K-corrections; the mean rather than the median per galaxy (means add); $P_{\rm fg}$ weighting; removal of brick duplicates; SGA instead of HECATE; no optional faint-galaxy cut. All are listed with their effects in CODE_CHECKS.md.

**Science consistency.** Units, cosmology, impact parameters, the observer-frame $1/(1+z)$ factor and the order of the subtraction are all correct. Issues found:

1. **Gas normalization.** Eq. B3 (paper, notebook and this code) normalizes the gas profile with the dark-matter $\rho_0$. The gas mass inside $r_{200}$ is then only $f_b(c)/f_{\rm dm}(c)\approx0.2$–$0.35$ of $f_{\rm gas}f_bM_{200}$. With halo masses from $P(M_h \mid M_\star)$ and one halo per group, the coded profile already saturates the halo-model bound along random sightlines, so the gas-mass-normalized profile would put more gas in halos than there is. The coded convention is kept.
2. **Mass error coefficient (fixed).** The r-band term in $\sigma(\log M_*)$ used 0.659 instead of $1.097+0.4=1.497$.
3. **Paper text vs code (the code is right).** The intercept of Eq. B1 is wrong as printed. The 0.018 size ratio is written inverted. The text cites Moster et al. (2010), but cell 26 uses the Moster et al. (2013) parameters.
4. **Coverage.** M31 and M33 are not in SGA, but no FRB lies within their $r_{200}$. The 30° SGA search radius never binds.

## In plain English

A fast radio burst travels billions of light years to reach us. On the way, its signal passes through free electrons, and each electron delays the low radio frequencies slightly more than the high ones. The total delay, the dispersion measure (DM), counts how many electrons the burst passed through. We want to know how many of those electrons belong to the galaxy the burst came from. In our paper the question is whether bursts from edge-on galaxies, which have to cross more of their galaxy's gas disc, pick up more electrons than bursts from face-on galaxies.

To get the host galaxy's share, we subtract what we know about everything else: our own Galaxy, and the average amount of gas between galaxies for a burst at that distance. A burst can also pass close to other galaxies that lie in front of its host. Each galaxy is surrounded by a large, thin cloud of hot gas, much bigger than the visible galaxy, and if the burst goes through one, its extra electrons get counted as the host's.

This code estimates that extra contribution for every burst. It finds the catalogued galaxies near the burst's position that are in front of the host and close enough to the line of sight. It estimates how heavy each one is from its brightness and colour, works out how big and dense its gas cloud is, and adds up the electrons the burst meets. Galaxies in the same group share one cloud, so a group is counted once. Because distances and masses are uncertain, the calculation is repeated thousands of times, which gives an error bar. The same calculation along random directions tells us how much of this every burst meets anyway, which the average we already subtracted includes.

What we found: most bursts pass through at least one foreground gas cloud, but the extra contribution is small for most of them, a median of 8% of the host's apparent share. Only the part above the random-direction average is removed. Removing it does not change the paper's conclusion: the host galaxy's share of the electrons shows no dependence on how tilted the galaxy is. The original method gave the heaviest galaxies halos up to cluster size and counted every member of a group separately, which made the foreground about twice what the halos can physically hold. Both are corrected here; the original numbers are kept alongside.
