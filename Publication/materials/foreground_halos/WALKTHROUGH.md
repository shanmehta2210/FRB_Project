# Walkthrough of `scripts/foreground_halo_dm.py`

What the code does, block by block with line numbers, and why each step is done that way. It estimates the dispersion measure (DM) that galaxy halos in front of each FRB host add to the burst. It follows Bhardwaj, Snelders et al. (2025, arXiv:2506.11915), Appendix B, and the original notebooks `Estimate_DM_halo.ipynb`, `Galaxy_offset_code.ipynb` and `atel_analysis.ipynb`.

Checks against the notebooks and of the physics are in [CODE_CHECKS.md](CODE_CHECKS.md); a summary is at the end of this file.

## The idea

The observed DM of an FRB is the sum of several pieces: the Milky Way, the intergalactic medium (DM$_{\rm cosmic}$), the host galaxy, and any galaxy halos the sightline happens to pass through on the way. We want DM$_{\rm host}$, which is normally obtained by subtracting the Milky Way term and the mean cosmic term (the Macquart relation) from the observed DM. Anything left over, including the gas of intervening halos, is wrongly assigned to the host. If such halos are common, they add noise to DM$_{\rm host}$ that could hide or fake a correlation with host inclination. This code estimates that halo contribution for each sightline so it can be removed.

Every galaxy sits in a dark-matter halo that contains hot ionized gas (the circumgalactic medium, CGM). If the FRB sightline passes within a foreground galaxy's halo radius $r_{200}$, the burst picks up extra DM. For each FRB the code:

1. finds galaxies near the sightline that could be in front of the host and close enough for their halo to reach it;
2. estimates each galaxy's stellar mass from its colour and brightness;
3. converts stellar mass to halo mass, then to halo size and concentration;
4. models the halo gas density and integrates the electron density along the sightline;
5. repeats this 2,000 times per galaxy with different masses and redshifts drawn from their uncertainties, so the result carries an error bar.

Steps 2–4 are a chain of empirical relations, each with scatter, and the DM is very non-linear in mass (it drops to zero when $r_{200}$ no longer reaches the sightline). That is why the uncertainties are propagated by Monte Carlo rather than analytically.

## Where each step is in the original notebooks

| Step | This code | Original |
|---|---|---|
| Stellar mass → halo mass (Moster+13) | `moster13_logms`, `logmh_from_logms` | `stellarmass_from_halomass`, `halomass_from_stellarmass` (`Estimate_DM_halo.ipynb`, cell 26) |
| Concentration | `dutton14_c200c` | `concentration_index` (Colossus `dutton14`, cell 26) |
| $r_{200}$ | `r200_kpc` | `virial_radius` (cell 26) |
| Gas profile and DM through one halo | `dm_mnfw` | `ModifiedNFW` (`_setup`, `rho_gas`, `n_e`, `compute_DM`), cell 26 |
| Stellar mass from colour | `logms_bell03_gr`; `logms_taylor11_gi` for validation | Bell+03 g−r in `atel_analysis.ipynb`; Taylor+11 g−i in cell 26 |
| Monte Carlo over mass (and photo-z) | `galaxy_dm_samples` | main loop of cell 26; truncated-Gaussian photo-z from the Appendix B text |
| Candidate selection | `ls_candidates` | `process_desi_bricks`, `check_intersecting_halo` and the photo-z cut (`Galaxy_offset_code.ipynb`) |
| Large nearby galaxies | `sga_wide` | `search_hecate_galaxies` (`Galaxy_offset_code.ipynb`) |
| Validation | `validate` | cell 26 inputs and output for FRB 20240114A |

## Setup (lines 1–54)

- **1–20, docstring:** a summary of the method and usage. Run `python scripts/foreground_halo_dm.py` for all hosts, or add `--controls` for the random comparison sightlines.
- **35–39, paths:** inputs are our host table (`host_props.csv`: FRB positions, redshifts, DM$_{\rm host}$) and our GALFIT fits (`seed_56.csv`: host magnitude and size, used to find and remove the host). Outputs go to `Publication/materials/foreground_halos/`. Every catalogue query is cached in `ls_cache/`, so reruns are reproducible and need no internet.
- **41–48, analysis constants:**
  - `SEARCH_DEG = 0.5`: search radius, as in the paper. Beyond 0.5° only very nearby galaxies can still reach the sightline, and those are handled separately (`SGA_DEG`).
  - `SIZE_RATIO_MIN = 0.018`: the size cut explained under `ls_candidates`.
  - `SGA_DEG = 30`: radius for the large-nearby-galaxy check, as in the paper (which used HECATE).
  - `N_MC = 2000`: draws per galaxy. The Monte Carlo error on a galaxy's mean DM is then about $1/\sqrt{2000} \approx 2\%$ of its spread, negligible next to the 50% systematic.
  - `MASS_SYS_DEX = 0.10`: systematic scatter added to every stellar mass, as in cell 26. It covers the intrinsic scatter of colour–mass relations, which photometric errors alone miss.
  - `DM_SYS_FRAC = 0.50`: the paper's 50% systematic on each total. It covers the choices of colour–mass relation, stellar-to-halo mass relation, concentration and gas profile, none of which the Monte Carlo varies.
  - `MR_BRIGHTEST = -24`: photo-z sanity limit (see `run_host`).
  - `RNG` seed 42: a fixed seed makes every run give identical numbers.
- **50–54, derived constants:**
  - $f_b = \Omega_b/\Omega_m \approx 0.158$ (astropy Planck18), the cosmic baryon fraction: the fraction of a halo's mass expected to be gas plus stars.
  - The proton mass, and the kpc-to-cm conversion.
  - `LOGMH_GRID`: halo masses from $10^9$ to $10^{15.5}\,M_\odot$ in 0.005 dex steps, used to invert the stellar-to-halo relation.
  - `BAND_WAVELENGTH`: effective wavelengths of the g, r, i, z and W1 bands, used for the K-correction.

## Halo physics (lines 57–111)

- **58–66, `moster13_logms`:** the Moster et al. (2013) stellar-to-halo mass relation,
  $$M_\star = 2N M_h \left[(M_h/M_1)^{-\beta} + (M_h/M_1)^{\gamma}\right]^{-1},$$
  with $N$, $M_1$, $\beta$ and $\gamma$ each evolving linearly in $z/(1+z)$. The parameters are identical to cell 26.
  - **Why:** we observe stars but need the halo, because the halo mass sets how far the gas extends ($r_{200}$) and how much gas there is. This relation comes from abundance matching (ranking galaxies by stellar mass and halos by mass, then pairing them), so it describes the average halo of a galaxy of given stellar mass at given redshift.
  - The relation has a peak efficiency near $M_h \approx 10^{12}\,M_\odot$. Above it, stellar mass grows slowly with halo mass, so a small error in the stellar mass of a massive galaxy becomes a large error in its halo mass. This is one reason massive galaxies dominate the uncertainty.
- **69–74, `logmh_from_logms`:** inverts the relation. It evaluates $M_\star(M_h)$ on the halo-mass grid at the galaxy's redshift and interpolates back.
  - **Why a grid:** the relation has no closed-form inverse. The notebook solves it with `fsolve` once per draw, which is slow for tens of thousands of draws. The relation is monotonic, so interpolation is exact to within the grid step; the checks find 0.0000 dex difference from `fsolve`.
- **77–80, `dutton14_c200c`:** halo concentration from Dutton & Macciò (2014),
  $$\log_{10} c = a + b\log_{10}\left[\frac{M_h}{10^{12}\,h^{-1}M_\odot}\right],\quad a = 0.520 + 0.385\,e^{-0.617 z^{1.21}},\quad b = -0.101 + 0.026z.$$
  - **Why:** the concentration $c = r_{200}/r_s$ sets how centrally peaked the halo is. A more concentrated halo puts more gas near the centre, which raises the DM for sightlines passing close to the galaxy. Lower-mass halos are more concentrated, because they formed earlier, when the Universe was denser.
  - This is the same relation that Colossus `dutton14` implements and the notebook calls. The relation is defined with mass in $h^{-1}M_\odot$, hence the factor `cosmo.h`. The notebook passes $M_\odot$ instead, the only difference between the two codes (see the summary).
- **83–89, `rho_crit`, `r200_kpc`:** the critical density $\rho_c(z)$ and
  $$r_{200} = \left[\frac{3M_h}{4\pi \cdot 200\,\rho_c(z)}\right]^{1/3},$$
  the radius inside which the mean density is 200 times critical.
  - **Why:** $r_{200}$ is the conventional halo boundary, and the paper (like the notebook) assumes the hot gas extends to $r_{200}$ and no further. A sightline passing outside $r_{200}$ contributes nothing. This is a convention; real gas extends somewhat beyond.
  - $\rho_c$ is evaluated at the galaxy's redshift because halos that formed earlier are denser and therefore smaller for the same mass.
- **92–93, `kpc_per_arcsec`:** proper (physical) kpc per arcsec at redshift $z$.
  - **Why proper rather than comoving:** a halo is a gravitationally bound object that does not expand with the Universe, and $r_{200}$ is a physical radius. The angular offset therefore has to become a physical impact parameter to be compared with it.
- **96–111, `dm_mnfw`:** the DM through one halo. It is vectorized: all 2,000 draws (different masses and redshifts) are computed in one array operation, which makes the Monte Carlo fast.
  1. **Halo properties per draw:** $c$, $r_{200}$, and the impact parameter $b = \theta \times$ (kpc per arcsec) at the draw's redshift.
  2. **Central density,** as in Eq. B3 of the paper and in cell 26:
     $$\rho_0 = \frac{200\rho_c}{3}\,\frac{c^3}{\ln(1+c) - c/(1+c)}.$$
     This is the normalization of an ordinary NFW dark-matter halo with mass $M_{200}$. The checks discuss its consequence: the gas inside $r_{200}$ adds up to less than $f_{\rm gas} f_b M_{200}$.
  3. **Sightline sampling:** 801 points along the chord through the sphere of radius $r_{200}$, from $-\sqrt{r_{200}^2 - b^2}$ to $+\sqrt{r_{200}^2 - b^2}$. The odd number puts one point exactly at closest approach, where the density peaks. The trapezoid rule on 801 points agrees with adaptive integration to 0.000%.
  4. **Gas density,** the modified NFW profile of Prochaska & Zheng (2019), paper Eq. B2:
     $$\rho_{\rm gas} = \frac{f_{\rm gas} f_b \rho_0}{y^{1-\alpha}(y_0+y)^{2+\alpha}}, \qquad y = c\,r/r_{200},$$
     with $\alpha = y_0 = 2$ and $f_{\rm gas} = 0.75$.
     - **Why modified:** a pure NFW profile rises as $1/r$ toward the centre. Hot gas does not follow that, because feedback and pressure flatten its core, and a $1/r$ cusp would give unrealistically large DMs for sightlines near the galaxy. With $\alpha = y_0 = 2$ the profile is flat in the core and follows NFW at large radius.
     - $f_{\rm gas} = 0.75$ is the assumed fraction of the halo's baryons that are in the hot CGM; the rest are in stars, cold gas or expelled.
  5. **Electron density:** $n_e = 1.1667\,\rho_{\rm gas}/(1.3\,m_p)$. Dividing by $1.3\,m_p$ (the gas mass per hydrogen atom, including helium) gives the hydrogen density. Multiplying by 1.1667 gives free electrons per hydrogen atom for fully ionized hydrogen plus helium ($1 + 2 \times 0.083$). The hot CGM ($\sim 10^6$ K) is fully ionized, so every electron counts.
  6. **DM:** $\mathrm{DM} = \int n_e\,ds$, with $ds$ in pc so DM is in pc cm⁻³. Draws with $b \ge r_{200}$ give exactly zero, as in the notebook.
  7. **Observer frame:** the result is divided by $(1+z_{\rm gal})$. Dispersion delay scales as $\nu^{-2}$, and the burst arrives at the halo at a frequency $(1+z)$ times higher than we observe, while the delay is stretched by $(1+z)$. The net effect is that a column $\int n_e\,ds$ at redshift $z$ adds $\int n_e\,ds/(1+z)$ to the observed DM. Expressing everything in the observer frame is what lets contributions at different redshifts be added.

## Stellar mass (lines 114–155)

- **115–121, `rest_mag`:** rest-frame AB magnitude at a chosen wavelength $\lambda_{\rm rest}$. It interpolates the observed magnitudes (g, r, plus i, z and W1 where available) linearly in $\log\lambda$, evaluates them at $\lambda_{\rm rest}(1+z)$, and adds $2.5\log_{10}(1+z)$.
  - **Why (the K-correction):** colour–mass relations are calibrated on rest-frame colours. At $z = 0.3$, the observed r band samples rest-frame light near 4,950 Å, close to the rest-frame g band. Using observed colours would make galaxies look redder and brighter than they are, and their masses too high. Interpolating the observed spectral energy distribution (SED) at the redshifted wavelength recovers the rest-frame magnitude. The $2.5\log_{10}(1+z)$ term accounts for the narrower rest-frame bandwidth of an observed band, for AB magnitudes.
  - The notebook used observed colours, which is fine for its $z \approx 0.1$ galaxies. Our hosts reach $z \approx 1$, so the correction matters: it lowers $\log M_\star$ by a median 0.29 dex.
- **124–130, `logms_bell03_gr`:** Bell et al. (2003) colour–mass-to-light relation in g−r,
  $$\log M_\star = 1.097\,(g-r)_0 - 0.406 - 0.4\,(M_r - 4.65),$$
  as used in `atel_analysis.ipynb`.
  - **Why colour:** redder galaxies have older stars and more mass per unit light, so colour predicts the mass-to-light ratio to about 0.1 dex. Multiplying by the luminosity (from the absolute magnitude $M_r$, with the Sun's $M_{r,\odot} = 4.65$) gives the mass.
  - −0.406 is Bell's zero-point (−0.306, for their "diet Salpeter" initial mass function) shifted by −0.1 dex to a Chabrier IMF, the convention of most modern SED fits, including the host masses we compare with.
  - **Why not the paper's relation:** the paper used Taylor et al. (2011) in g−i with SDSS photometry. Legacy Surveys DR9 has no i band and SDSS does not cover most of our fields, so we use the g−r relation from `atel_analysis.ipynb`.
- **133–136, `logms_taylor11_gi`:** the paper's Taylor et al. (2011) relation, $\log M_\star = 1.15 + 0.70(g-i) - 0.4M_i$. It is used only by `validate`, to reproduce the notebook's numbers with the notebook's own inputs.
- **139–155, `galaxy_dm_samples`:** the Monte Carlo for one galaxy.
  - **Spectroscopic redshift:** all draws use it, and the galaxy is certainly in front ($P_{\rm fg} = 1$).
  - **Photometric redshift only:** draws come from $\mathcal{N}(z_{\rm phot}, \sigma_z)$ truncated to $0.001 \le z \le z_{\rm host}$, as in the paper.
    - **Why truncate:** a galaxy behind the host contributes nothing, so the DM we want is the DM if the galaxy is in front. The lower bound avoids $z \to 0$, where the angular scale and the inferred mass become meaningless.
    - **$P_{\rm fg}$:** $\Phi(\text{upper}) - \Phi(\text{lower})$, the probability that the galaxy really is in front. The paper implicitly sets it to 1. That is reasonable for one sightline at $z = 0.13$, but over many high-redshift hosts it lets galaxies that are probably behind (for example $z_{\rm phot} = 0.5 \pm 0.2$ for a host at 0.3) contribute as if they sat just in front of the host, which inflates the total. Weighting by $P_{\rm fg}$ gives the expected contribution instead. `run_host` keeps both versions.
  - **Mass draws:** Gaussian in $\log M_\star$ with
    $$\sigma = \sqrt{(1.097\,\sigma_g)^2 + (1.497\,\sigma_r)^2 + 0.10^2}.$$
    The coefficients are the derivatives of the Bell relation with respect to $g$ (1.097) and $r$ ($1.097 + 0.4 = 1.497$, because $r$ enters both the colour and the luminosity). For the Taylor relation they are 0.70 and 1.10, as in cell 26. The 0.10 dex term is the systematic scatter.
  - Each draw is mapped through `logmh_from_logms` and `dm_mnfw`. The mass is computed at each draw's own redshift, because the same apparent magnitude means a more luminous, more massive galaxy further away.
  - **Returns:** the 2,000 DM samples (observer frame) and $P_{\rm fg}$.

## Catalogue queries (lines 158–242)

All queries go to NOIRLab's Astro Data Lab and are cached as `ls_cache/<FRB>_<name>.csv`.

- **159–168, `dl_query`:** runs one SQL query, retrying up to 3 times because the service occasionally times out.
- **171–177, `cached`:** returns the cached CSV if present, otherwise queries and saves. Every result therefore depends only on files in the repository, not on the current state of the online catalogues.
- **180–191, `ls_footprint`:** checks whether the field is in LS DR9 (any source with a photo-z within 0.02°). If not, it checks DR10, which adds southern sky from DECam. If neither covers the field, the FRB is skipped, since without the photo-z catalogue no foreground estimate is possible. This is why 49 of our 55 hosts have an estimate.
- **194–207, `ls_candidates`:** the paper's candidate selection, done in SQL. It joins the Tractor photometry catalogue to the LS photo-z catalogue on `ls_id` and keeps sources that:
  - lie within 0.5° of the FRB;
  - are extended (`type` not `PSF` or `DUP`). PSF sources are stars or unresolved quasars, and DUP entries are Gaia sources inside large galaxies that carry no flux of their own;
  - are `brick_primary`. LS is processed in overlapping sky tiles ("bricks"), so a source near a tile edge appears twice; this keeps one copy. The notebook read the bricks directly and did not apply this;
  - satisfy $R_{\rm eff} > 0.018\,\theta$ (the size cut);
  - satisfy $z_{\rm phot} - 2\sigma_z \le z_{\rm host}$: anything consistent with being in front of the host at the 2σ level is kept, so true foreground galaxies are rarely lost to photo-z errors.

  **Why the size cut works:** galaxy size scales with halo size, $R_{\rm eff} \approx 0.018\,r_{200}$ (Somerville et al. 2018). The sightline is inside the halo when $\theta < r_{200}$ in angle, i.e. $\theta < R_{\rm eff}/0.018$, i.e. $R_{\rm eff}/\theta > 0.018$. Both quantities are angles, so the cut needs no redshift, which makes it a cheap first filter. It is approximate; the exact test ($b < r_{200}$ for each draw) happens later in `dm_mnfw`. (The paper writes the ratio inverted; the notebook and this code use $R_{\rm eff}/\theta$.)
- **210–230, `spectroscopic_redshifts`:** all spectroscopic redshifts within 0.5° from DESI DR1 (`zwarn = 0`, galaxies, primary spectrum), SDSS DR17 (`zwarning = 0`, galaxies) and SGA-2020 (HyperLEDA redshifts).
  - **Why:** photo-z errors ($\sigma_z \sim 0.02$–0.1) are the largest uncertainty. A spectroscopic redshift fixes the distance, and therefore the mass and the physical impact parameter, and also says definitively whether the galaxy is in front. The quality flags exclude failed redshift fits.
- **233–242, `sga_wide`:** galaxies of the Siena Galaxy Atlas (SGA-2020) within 30° with $0.0005 < z \le z_{\rm host}$.
  - **Why:** a very nearby galaxy can have a halo spanning degrees on the sky, so it can reach the sightline from far outside the 0.5° search. SGA contains the large, nearby galaxies, with photometry measured on the whole galaxy (the standard Tractor catalogue breaks such galaxies into pieces). The paper did this check with HECATE; SGA is in the same Data Lab service.
  - $z > 0.0005$ removes Local Group objects with unreliable redshift-based distances. The checks confirm that no FRB lies inside the halos of M31 or M33, the only Local Group galaxies that would matter.

## Cleaning the candidate list (lines 245–302)

- **246–247, `separation_arcsec`:** flat-sky angular separation, accurate at the arcsecond scales where it is used.
- **250–265, `attach_specz`:** best redshift for each candidate. It uses the LS compiled spectroscopic redshift if there is one. Otherwise it takes the nearest external spectroscopic redshift within the larger of 1.5″ and $0.2\,R_{\rm eff}$.
  - **Why that tolerance:** 1.5″ matches a fibre to its galaxy without false matches to neighbours. For large galaxies the fibre or catalogue position can sit away from the Tractor centroid, so the tolerance grows with galaxy size.
- **268–278, `identify_host`:** finds the FRB host in the LS catalogue. Among sources within the larger of 5″ and $3R_e$ of the FRB, it takes the one closest in r magnitude to our GALFIT host (if within 1.5 mag), otherwise the nearest.
  - **Why:** the host's own halo and ISM are part of DM$_{\rm host}$, the quantity we are measuring, so it must not be counted as foreground. Matching on brightness as well as position avoids picking a faint neighbour when the FRB is offset from its host's centre.
- **281–287, `remove_host`:** drops the host and its "shreds". A shred is any source within the largest of $2R_e$, $2\times$ the LS size and 2″ of the host, unless its spectroscopic redshift differs from the host's by more than 0.01.
  - **Why:** Tractor often splits a large or clumpy galaxy into several sources (spiral arms, star-forming knots). Each piece would otherwise be treated as a separate galaxy with its own halo right on top of the sightline, which would add a large, spurious DM. A source with a clearly different spectroscopic redshift is a real separate galaxy and is kept.
- **290–293, `mag_err`:** magnitude error from flux and inverse variance, $1.0857/(\text{S/N})$, where $1.0857 = 2.5/\ln 10$. It is capped at 1 mag, and set to 0.3 mag where it is undefined (zero or negative flux), so that a few bad measurements cannot produce absurd mass scatter.
- **296–298, `ab_mag`:** converts SGA fluxes (nanomaggies, zero-point 22.5) to AB magnitudes, dividing by the Milky Way transmission to correct for Galactic dust. LS candidates already come with dust-corrected (`dered_mag_*`) magnitudes.
- **301–302, `sample_summary`:** mean, 16th and 84th percentile of a galaxy's DM samples.

## Per FRB (lines 305–398)

- **306–323, `load_candidates`:** queries the candidates, converts columns to numbers (cached CSVs can contain text), attaches spectroscopic redshifts and removes the host. It then drops galaxies whose spectroscopic redshift puts them behind the host ($z > z_{\rm host} + 0.001$), and any without valid g and r magnitudes, which are needed for the mass.
  - The 0.001 margin keeps galaxies at essentially the host redshift (group members whose velocities put them slightly behind). They are placed at $z_{\rm host}$ in the next step.
- **326–398, `run_host`:** the main loop for one sightline.
  - **330–338:** set up the result row; skip the FRB if it is outside the LS footprint; count candidates and spectroscopic redshifts.
  - **340–343:** three running totals, kept side by side so that the effect of each assumption can be tested:
    - `total`: the main estimate. In each of the 2,000 draws, a photo-z galaxy is included with probability $P_{\rm fg}$, so the mean of `total` is the expected foreground DM and its spread includes the uncertainty about which galaxies are in front;
    - `total_conditional`: the paper's convention, every candidate assumed in front (the sum of each galaxy's mean DM);
    - `total_specz`: spectroscopic galaxies only, the most secure part of the total.
  - **345–370:** for each LS candidate:
    1. Cap a spectroscopic $z$ at $z_{\rm host}$ (see `load_candidates`).
    2. Skip photo-z galaxies that would be brighter than $M_r = -24$ if placed at $\min(z_{\rm phot}, z_{\rm host})$. Even the brightest cluster galaxies are around $-23.5$, so this means the photo-z is wrong: usually a star or an AGN, or a galaxy that is really much closer. Kept, it would be assigned an enormous halo directly on the sightline and dominate the total. `n_photoz_implausible` records how many are skipped.
    3. Draw DM samples. $\sigma_z$ is floored at 0.01 because the catalogue sometimes quotes unrealistically small photo-z errors.
    4. Add the samples to the three totals.
    5. Record one row in `foreground_galaxies.csv` (position, photometry, redshift, $P_{\rm fg}$, mean and 16–84% DM), so each host's total can be traced to its individual galaxies.
  - **372–388:** for each SGA galaxy beyond 0.5°, compute its magnitudes and halo at its nominal mass. If the sightline is inside $r_{200}$, draw samples (0.05 mag errors, which is typical for these bright galaxies, and the HyperLEDA redshift) and add them to all totals.
    - **Why test the nominal mass first:** SGA returns up to tens of thousands of galaxies per field. Almost none reach the sightline, so running the full Monte Carlo for every one would be very slow. A galaxy whose average halo misses could still be reached by a high-mass draw, but such contributions come from the edge of the halo, where the density is lowest.
  - **390–398:** the per-host summary: mean and standard deviation of `total`, the conditional and spec-z-only sums, and the largest single contributor (`top_gal`, `top_dm`), which shows whether a total rests on one galaxy.

## Driver (lines 401–460)

- **405–415, `validate`:** reruns the notebook's six galaxies for FRB 20240114A, with the notebook's SDSS g and i photometry, spectroscopic redshifts, offsets and Taylor+11 masses.
  - **Why:** it checks that the physics chain reproduces the original code with identical inputs, independently of our different photometry and selection.
  - The notebook printed medians of 15.9, 4.6, 0, 104.4, 5.1 and 0 pc cm⁻³. This code gives 16.9, 4.7, 0, 109.7, 5.2 and 0. The difference is entirely the units passed to Colossus in the notebook (see the summary). The end-to-end test below uses different random draws, hence 110.6 there for galaxy 4.
- **418–438, `main`:** runs the validation, then every host with a redshift. It writes:
  - `foreground_dm_per_host.csv`, with total uncertainty $\sigma = \sqrt{(0.5\,\mathrm{DM_{fg}})^2 + \sigma_{\rm MC}^2}$ (the paper's 50% systematic combined with the Monte Carlo spread);
  - `foreground_galaxies.csv`, one row per galaxy;
  - `validation_20240114A.csv`.
- **441–460, `controls`:** runs the same calculation along two random sightlines 0.75° from each host, at the host's redshift, with no host removal and 500 draws per galaxy. The output is `foreground_dm_controls.csv`.
  - **Why:** the Macquart relation gives the mean DM$_{\rm cosmic}$ averaged over all sightlines, and that average already includes the gas in galaxy halos. Subtracting the full foreground DM would therefore remove the typical halo contribution twice. The random sightlines measure that typical contribution at the same redshift and in the same part of the sky (same survey depth and Galactic dust), so the excess over them is what is genuinely special about each FRB sightline. 0.75° puts the control just outside the 0.5° search area, so its candidates do not overlap the FRB's.

## How the outputs are used

`robust_correlations.load` (lines 46–61 of `scripts/robust_correlations.py`) builds two corrected host DMs, in the host's rest frame:

- `dm_host_fg` $= \mathrm{DM_{host}} - \mathrm{DM_{fg}}(1+z_{\rm host})$: the full foreground contribution removed;
- `dm_host_fg_excess` $= \mathrm{DM_{host}} - (\mathrm{DM_{fg}} - \langle \mathrm{DM_{fg}}\rangle_{\rm random})(1+z_{\rm host})$: only the excess over random sightlines removed.

The factor $(1+z_{\rm host})$ converts the observer-frame foreground DM to the host rest frame, the frame in which DM$_{\rm host}$ is defined. The two versions bracket the truth. The first over-corrects (it removes halo gas already in the Macquart mean) and makes DM$_{\rm host}$ negative for a few hosts. The second is the consistent correction but relies on the control sightlines being typical. Both are carried through the correlation tests.

`scripts/foreground_sightline_counts.py` counts the halos intersected per sightline (`SIGHTLINES.md`).

## Summary of the checks

Full numbers are in [CODE_CHECKS.md](CODE_CHECKS.md), produced by `scripts/check_foreground_halo_dm.py`.

**Agreement with the original code.** The physics is the same as `Estimate_DM_halo.ipynb` cell 26: Moster+13 parameters, the Eq. B3 profile normalization, $f_{\rm gas}=0.75$, $\alpha=y_0=2$, $\mu_H=1.3$, $n_e/n_H=1.1667$, zero DM outside $r_{200}$ and the $1/(1+z)$ factor. The candidate selection is the same as `Galaxy_offset_code.ipynb` and the Appendix B text.

- **End-to-end test on the notebook's own inputs.** The six FRB 20240114A galaxies of cell 26 were run through a re-implementation of the notebook and through this code, with the same 2,000 mass draws in both.
  - The notebook re-run reproduces its published medians exactly: 15.9, 4.6, 0, 104.4, 5.1 and 0 pc cm⁻³, a total of 130.0.
  - This code gives 16.9, 4.7, 0, 110.6, 5.2 and 0, a total of 137.3.
- **The only difference is a unit slip in the notebook.** It passes the halo mass to Colossus in $M_\odot$, but Colossus expects $h^{-1}M_\odot$. This makes its concentrations about 3.5% low and its DMs up to 6.5% low. With the unit corrected, the notebook gives exactly this code's numbers (137.3).
- **Component checks.** On an 80-case grid, halo masses agree to 0.0000 dex and the integration to 0.000%.
- **Deliberate differences.** Bell g−r masses (DR9 has no i band), K-corrections, the mean rather than the median per galaxy (means add), $P_{\rm fg}$ weighting, removal of brick duplicates, SGA instead of HECATE, and no optional faint-galaxy cut (that cut is in the notebook but not in the paper). All are listed with their effects in CODE_CHECKS.md.

**Science consistency.** Units, cosmology, impact parameters, the observer-frame $1/(1+z)$ factor and the order of the subtraction are all correct. Issues found:

1. **Gas normalization.** Eq. B3 (paper, notebook and this code) normalizes the gas profile with the dark-matter $\rho_0$. The gas mass inside $r_{200}$ is then only $f_b(c)/f_{\rm dm}(c)\approx0.2$–$0.35$ of $f_{\rm gas}f_bM_{200}$, so halo DMs are about $3.8\times$ lower than the gas-mass-normalized profile of the `frb` package (Milky Way: 14 vs 61 pc cm$^{-3}$). This is not clearly wrong. Along random sightlines the coded halos already supply about 17% of the mean $\rm DM_{cosmic}$, and the gas-mass-normalized version would supply about 64% from LS galaxies alone, which is too much. The coded convention therefore behaves like an effective hot-gas fraction of 0.13–0.26, and $f_{\rm hot}=0.75$ is an upper bound. With the upper bound, DM$_{\rm host}$ vs $q$ is still null, and DM–$|{\rm RM}|$ weakens to $r_S=+0.34$ ($p=0.045$).
2. **Mass error coefficient (fixed).** The r-band term in $\sigma(\log M_*)$ used 0.659 instead of $1.097+0.4=1.497$; it now uses 1.497. The median $\sigma$ rises from 0.126 to 0.147 dex, and per-host foreground DMs change by at most 1.3 pc cm⁻³.
3. **Paper text vs code (the code is right).** The intercept of Eq. B1 is wrong as printed. The 0.018 size ratio is written inverted. The text cites Moster et al. (2010), but cell 26 uses the Moster et al. (2013) parameters.
4. **Coverage.** M31 and M33 are not in SGA, but no FRB lies within their $r_{200}$. The 30° SGA search radius never binds (the largest angular $r_{200}$ is 5.5°).

## In plain English

A fast radio burst travels billions of light years to reach us. On the way, its signal passes through free electrons, and each electron delays the low radio frequencies slightly more than the high ones. The total delay, the dispersion measure (DM), counts how many electrons the burst passed through. We want to know how many of those electrons belong to the galaxy the burst came from, because that tells us about the gas in that galaxy. In our paper the question is whether bursts from edge-on galaxies, which have to cross more of their galaxy's gas disc, pick up more electrons than bursts from face-on galaxies.

To get the host galaxy's share, we subtract what we know about everything else: our own Galaxy, and the average amount of gas between galaxies for a burst at that distance. The problem is that a burst can also pass close to other galaxies that happen to lie in front of its host. Each galaxy is surrounded by a large, thin cloud of hot gas, much bigger than the visible galaxy. If the burst goes through one of these clouds, its extra electrons get wrongly counted as belonging to the host.

This code estimates that extra contribution for every burst. It looks up all the catalogued galaxies near the burst's position on the sky and keeps the ones that are in front of the host and close enough to the line of sight. It estimates how heavy each one is from its brightness and colour, and works out from that how big its gas cloud is and how dense the gas is. It then adds up the electrons the burst would meet passing through each cloud. Because the distances and masses of these galaxies are uncertain, the calculation is repeated thousands of times with slightly different values, which gives an error bar. The same calculation along random directions on the sky tells us how much of this is just the typical amount every burst encounters, which the average we already subtracted includes.

What we found: most bursts do pass through at least one foreground gas cloud. The typical extra contribution is modest, about a quarter of the host's apparent share, and for some bursts it is large. Removing it does not change the paper's conclusion: the host galaxy's share of the electrons still shows no dependence on how tilted the galaxy is. We also checked that this code reproduces the original notebooks; the one difference traces to a units slip in the notebook, and its effect is small.
