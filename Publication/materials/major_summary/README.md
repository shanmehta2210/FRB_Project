# Major summary: the evidence for the three main results

Everything in this folder is rebuilt by `scripts/major_summary.py` from the outputs of the analysis scripts; run `scripts/binned_correlations.py` first. Each figure is saved as PDF (for the paper or slides) and PNG (for quick viewing). Tables have been rendered as figures. Their source CSVs are named below so any number can be traced back.

Conventions used throughout:

- $p_{\rm m}$ is the fraction of 10,000 random control samples of the same size that differ from the control at least as much as the FRB hosts do. A small $p_{\rm m}$ would mean the hosts are unusual.
- $r_S$ is the Spearman rank correlation.
- "BH q" is the Benjamini–Hochberg false-discovery-rate value. It accounts for running many tests at once: with 40 tests, about two will reach $p < 0.05$ by chance alone.
- "MC" means 2,000 Monte Carlo draws through the measurement errors.
- Every Spearman test is paired with an AD test, reported as "AD p". Spearman asks whether $y$ rises or falls with $x$. The AD test splits the hosts at the median of $x$ and asks whether $y$ has the same distribution in the two halves (two-sample Anderson–Darling, permutation $p$). The AD test also catches dependences that are not monotonic, such as a change in spread, but it has less power for a simple trend.
- This AD p is a different test from the AD used for the inclination distributions. There, the host $\cos i$ values are compared with the control sample, and the result is the $p_{\rm m}$ above.
- In the binned-correlation plots, a legend entry such as
  `z 0.01–0.14: r_S = +0.27 (p = 0.31); AD p = 0.52`
  `MC r_S +0.29 [+0.20, +0.38]; p < 0.05 in 1% (AD 0%)`
  reads as follows:
  - The first line gives the Spearman $r_S$ of the plotted pair in that tertile with its two-sided $p$-value, then the AD $p$.
  - The second line is the same $r_S$ recomputed in each Monte Carlo draw, with every value resampled within its error bars and the tertiles re-formed. It gives the median and [16th, 84th percentile] of $r_S$; these are correlation coefficients (so they can be negative), not $p$-values.
  - The percentages are the fractions of draws in which Spearman, and in brackets AD, give $p < 0.05$.
- DM$_{\rm host}$ is in the host rest frame.

| Result | Folder | One-line verdict |
|---|---|---|
| 1. Inclination null hypothesis | `1_inclination_null/` | FRB hosts have the same inclination distribution as field disks ($p_{\rm m} = 0.54$); no burst or host property depends on inclination. |
| 2. Property correlations | `2_correlations_null/` | Within bins of each other, DM$_{\rm host}$ does not track $z$ or $M_*$. The only notable link is with SFR, among $z \gtrsim 0.3$ hosts, and it weakens once errors are included. |
| 3. Late-type hosts | `3_late_type_hosts/` | 50 of 56 hosts are disks or have disk-like profiles; one is early-type. |

---

## 1. Inclination null hypothesis (`1_inclination_null/`)

**Claim.** If FRB detection does not depend on how the host disk is oriented, FRB hosts should share the inclination distribution of field disk galaxies measured the same way: 24,450 HSC galaxies with PSF-convolved Sérsic fits, $r \le 22$, $0.4 < n < 1.5$ and $q > 0.2$. They do. The face-on excess reported by Bhardwaj et al. (2024) is reproduced only when the PSF is ignored.

**Main evidence**

| File | What it shows | Key numbers |
|---|---|---|
| `01_cdf_hosts_vs_control` | CDF of $\cos i$: (a) 56 hosts vs control with the 68% band of random samples of 56; (b) the 23 hosts of Bhardwaj et al. (2024), published vs our remeasurement | Hosts lie inside the band; median $i$ 51.3° vs 54.9° for the control |
| `02_quartile_bins` | Hosts in the four control-quartile bins of $\cos i$ | 14, 19, 10, 13 from face-on to edge-on; 14 expected in each |
| `02b_table_headline_tests` | All headline numbers in one table: medians, quartile counts, $p_{\rm m}$ for 56 and 23 hosts, the measurement-error runs, and statistical power | AD / MW two-sided / MW one-sided $p_{\rm m}$ = 0.54 / 0.46 / 0.23 |
| `03_pvalue_summary_all_tests` | **Every inclination test in one picture.** $p_{\rm m}$ for the adopted test, the Bhardwaj hosts, measurement errors, 11 control variants, 4 values of $q_0$, morphology subsets and redshift bins; the no-PSF fits are shown in grey for contrast | All PSF-convolved tests have $p_{\rm m} > 0.05$; the no-PSF fits and the published values are at $p_{\rm m} \le 0.004$ |
| `04_measurement_error_realizations` | Distribution of $p_{\rm m}$ when every $q$ is redrawn from its error and a 5° systematic is added | Median AD 0.50; no AD or two-sided MW realization below 0.05, and 0.4% of one-sided ones |

**The null result is not caused by low sensitivity.** The headline table gives the power: the AD test detects the published deficit at $p < 0.05$ in 99.95% of samples of 56 and 90.9% of samples of 23. Samples are drawn from the 23 published inclinations.

**Robustness to the choice of control sample and $q_0$**

| File | What it shows | Source |
|---|---|---|
| `05_table_control_variants` | $p_{\rm m}$ for Sérsic-index windows, magnitude limits (control only, or both samples) and colour or star-formation cuts | `appendix/robustness/A2,A3,A5_*.csv` |
| `06a/b/c_control_variants_*` | CDFs for each variant | `appendix/robustness/figA_ctrl_*` |
| `07_table_q0`, `08_q0_cdfs` | $q_0$ = 0.10–0.25 | `appendix/robustness/C1_q0.csv` |
| `10_table_redshift_bins` | Tests in redshift bins, against the full control and against a control matched in photometric redshift | `appendix/redshift_bins/inclination_zbins.csv` |
| `11_control_cosi_histogram` | The control is not uniform in $\cos i$ (39.8% have $i > 60°$ vs 50% for random orientations), which is why we compare with a measured control rather than a uniform distribution | — |

The lowest values in any of these are $p_{\rm m} = 0.055$ (AD, $z < 0.1$ against the matched control, 12 hosts) and 0.09 (one-sided MW, $q_0 = 0.25$). Neither is significant, and neither survives a correction for the number of variants tried.

**Why the earlier study found a face-on excess**

| File | What it shows | Key numbers |
|---|---|---|
| `09_table_psf_fit_variants` | GALFIT and AstroPhot with and without the PSF; the 22 shared hosts, ours vs published | Without the PSF, median $i$ drops from 51.3° to 38.0°, 32 hosts fall in the face-on bin, and $p_{\rm m} \le 10^{-4}$. The published values of the 22 shared hosts give $p_{\rm m} = 0.004$; our values give 0.96. |

**No property depends on inclination**

| File | What it shows | Key numbers |
|---|---|---|
| `12a_dm_host_vs_axis_ratio` | DM$_{\rm host}$ vs $q$, with and without foreground halos; forest plot of $r_S$ for every property | DM$_{\rm host}$: $r_S = -0.06$ ($p = 0.64$) |
| `12_table_correlations_with_axis_ratio` | Outlier-robust statistics for every property vs $q$: AD, bootstrap, Kendall $\tau$, jackknife, Theil–Sen; with foreground halos removed and in subsamples | All 12 DM$_{\rm host}$ variants (with and without foreground halos, corrected and published methods, three subsamples) have $\lvert r_S \rvert \le 0.11$, $p \ge 0.50$ and AD $p \ge 0.08$. The only hit is old-star $A_V$ ($-0.46$, $p = 0.01$; AD $p = 0.003$; $N = 30$), the expected sense for dust: edge-on disks are more attenuated. |
| `13_marnoch_shared_axis_ratios` | Our $q$ vs Marnoch et al. (2026) for the 20 shared hosts | $r_S = 0.86$, median $\lvert\Delta q\rvert = 0.031$, so the measurements agree |
| `14_table_marnoch_combined` | DM$_{\rm host}$ vs $q$: Marnoch alone, ours alone, and the union of 69 hosts | Spearman $-0.43 \to -0.06$, union $-0.21$ ($p = 0.08$). AD $p$: Marnoch 0.032, ours 0.38, union 0.09–0.15. |
| `15_table_correlations_with_q_in_redshift_bins` | Property vs $q$ within redshift bins, Spearman and AD | No consistent trend in any bin. The only AD $p < 0.05$ is old-star $A_V$ (all hosts and $z \ge 0.2$). |
| `16_binned_properties_vs_inclination_grid` | **16 panels:** $z$, $M_*$, SFR and DM$_{\rm host}$ against inclination, within tertiles of each of the four (rows = binning quantity). Error bars on every point; Spearman and AD, also through the errors (MC). | Whole sample: DM$_{\rm host}$ $+0.06$ ($p = 0.64$; AD 0.19), SFR $-0.12$ (0.42; AD 0.46), $M_*$ $+0.11$ (0.43; AD 1.00), $z$ $-0.02$ (0.86; AD 0.53). Through the errors, $p < 0.05$ in at most 1% of draws for Spearman and 6% for AD (5% is chance). |
| `17_binned_properties_vs_inclination_heatmap` | All 64 coefficients of the 16 panels in one colour-coded table, with AD $p$ in every tertile and whole-sample cell | Spearman: 2 of 48 tertiles reach $p < 0.05$. AD: also 2 of 48. Both are chance level, and the lowest BH $q$ is 0.62 (Spearman) and 0.96 (AD). |
| `18_global_tests_both_grids` | Permutation test of each grid as a whole (bottom row: grid B). Inclination is shuffled across hosts 5,000 times, which keeps the dependence between the 48 tertile tests. | Hits 2 vs 2.5 expected under the null ($P = 0.67$ Spearman, 0.70 AD). Family-wise $p$ of the best tertile 0.48 (Spearman) and 0.74 (AD); Fisher combined $p$ 0.73 and 0.82. Grid B is noise. |

The nominal hits in panel 16 are:

- **SFR vs $i$ for the least massive third:** Spearman $r_S = -0.64$, MC $-0.58$, BH $q = 0.62$; AD $p = 0.22$. If real, this is most simply the inclination dependence of dust attenuation leaking into the SED-derived SFRs of dwarfs.
- **$M_*$ vs $i$ for the middle DM third:** $+0.53$, MC $+0.39$, $q = 0.66$; AD $p = 0.034$.
- **SFR vs $i$ for the lowest-SFR third:** AD only, $p = 0.042$; Spearman $p = 0.087$.
- **DM$_{\rm host}$ vs $i$ for the 51 hosts with $M_*$:** AD $p = 0.043$. For all 55 hosts it is AD $p = 0.19$, and through the errors it is significant in 6% of draws, so it does not survive.

---

## 2. Correlations among the other properties under binning (`2_correlations_null/`)

**Claim.** Once each quantity is held fixed by binning, the host properties ($z$, $M_*$, SFR) show only the expected selection trends with each other. DM$_{\rm host}$ does not track redshift or stellar mass. The one exception is a positive DM$_{\rm host}$–SFR link confined to the most distant third of the sample; it should be presented as a lead, not a result.

| File | What it shows | Key numbers |
|---|---|---|
| `01_binned_pairs_grid` | **12 panels:** for each binning quantity (rows: $z$, $M_*$, SFR, DM$_{\rm host}$), the three pairs of the other three, by tertile. Error bars on every point; MC statistics through the errors, re-binning in every draw. | See below |
| `02_binned_pairs_heatmap` | All coefficients of the 12 panels in one colour-coded table, with AD $p$ and the tertiles that survive BH framed | Spearman: 2 of 36 tertiles survive BH, both DM–SFR. AD: 3 of 36 have $p < 0.05$, none survives BH (lowest $q = 0.06$). |
| `03_dm_host_vs_sfr_binned_by_z` | The one notable panel, full size | $z > 0.29$: $r_S = +0.82$ ($p < 0.001$; AD $p = 0.009$). Through the errors $+0.57$ [16–84%: $+0.37$, $+0.72$], $p < 0.05$ in 64% of draws (AD 42%). |
| `04_table_dm_host_tests` | Whole-sample tests of what DM$_{\rm host}$ correlates with: Spearman and AD, each with BH $q$, plus bootstrap, delete-3 jackknife, partial $\tau$ given $z$ and the DM$_{\rm cosmic}$ MC | Only DM–$\lvert$RM$\rvert$ survives Spearman BH ($q = 0.04$; AD $p = 0.036$, $q = 0.19$). Nothing survives AD BH; the lowest is DM–scattering (AD $p = 0.002$, $q = 0.075$). DM–$M_*$: $-0.01$ (AD 0.91). DM–$z$: $+0.08$ (AD 0.92). |
| `05_dm_host_physics_panels` | Six scatter plots: host-halo DM, SFR, $\Sigma_{\rm SFR}$, scattering, RM, and $M_*$ vs $z$ | Host-halo DM is only about 7% of DM$_{\rm host}$ |
| `06_foreground_halo_dm_per_frb` | **Intervening-halo DM for every FRB**, ordered by redshift, using the audited estimate; grey columns are outside the Legacy Survey footprint. (a) Observed DM split into DM$_{\rm MW}$, the diffuse cosmic gas $\langle$DM$_{\rm cosmic}\rangle - \langle$DM$_{\rm fg}\rangle(z)$, this sightline's foreground halos and what is left for the host. (b) Corrected and published foreground DM, with the expectation for a random sightline. (c) DM$_{\rm host}$, standard and with the foreground accounted for. (d) What the audit removed from the published estimate. | 49 FRBs in the footprint. Median foreground DM: 18 pc cm$^{-3}$ (45 published). 3 FRBs are above 100 (13 published). The largest share of an observed DM is 19% (44% published). Foreground is a median 8% of DM$_{\rm host}$. 3 of 49 DM$_{\rm host}$ values are negative (1 with the standard estimate), all within $1\sigma$ of zero. |
| `08_foreground_halo_audit` | **Audit of the foreground estimate and the DM budget** (`foreground_halos/audit/AUDIT.md`). (a) Random sightlines against the halo-model maximum. (b) Halo masses assigned to massive galaxies. (c) Budget closure for each FRB. (d) Distribution of the budget residuals in units of $\sigma$. | Subtracting the full foreground on top of $\langle$DM$_{\rm cosmic}\rangle$ counts the mean halo gas twice; only the excess over random sightlines may be removed. The published estimator gives 2.2 times the physical maximum, mostly because halo masses come from inverting $M_*(M_h)$ (up to $10^{15.5}\,M_\odot$) and every group galaxy gets its own halo. The corrected estimator sits at the bound (1.01). No FRB's budget is more than $1\sigma$ short. DM$_{\rm host}$ vs $q$ stays null and high-$z$ DM–SFR stays at $+0.73$. |
| `09_halo_mass_inversion_test` | **The published halo-mass step, run with the notebook code itself** (`foreground_halos/mass_inversion/MASS_INVERSION.md`). (a) Halo mass at fixed stellar mass: inverted mean relation vs $P(M_h \mid M_*)$, with FRB 20240114A's galaxy 6 and the measured mass of its cluster. (b) DM profile of galaxy 6 for both masses. (c) Galaxy 6 DM over Monte Carlo draws. (d) The FRB 20240114A DM budget, published vs corrected. | The notebook reproduces the paper's 104 pc cm$^{-3}$ for galaxy 6 only through the inverted relation, which gives it $10^{15.33}\,M_\odot$: 23 times its whole cluster ($10^{13.97}$). $P(M_h \mid M_*)$ gives $10^{14.57}$ and 61 pc cm$^{-3}$. DM$_{\rm c,halo}$ falls from 135 to 30–91, the budget DM$_{\rm host}$ rises from 126 to 171–232 (still within $1\sigma$ of H$\alpha$), and the foreground no longer dominates. |
| `07_global_tests_both_grids` | Permutation test of each grid as a whole (top three rows: grid A). DM$_{\rm host}$ is shuffled across hosts for the DM panels; $y$ is shuffled within tertiles for the pairs among $z$, $M_*$ and SFR. | Spearman: 5 hits vs 1.9 expected ($P = 0.048$); the high-$z$ DM–SFR tertile has a family-wise $p$ of 0.008 and Fisher $p = 0.004$. AD: 3 hits ($P = 0.34$), family-wise $p$ 0.056. |

**Reading the 12-panel grid**

- **DM$_{\rm host}$ vs $z$ and vs $M_*$:** no correlation overall or at fixed values of the other quantities. Partial $\lvert\tau\rvert \le 0.06$; through the errors $\lvert r_S \rvert \le 0.05$; whole-sample AD $p \ge 0.90$.
- **$z$ vs $M_*$ and vs SFR:** positive throughout ($r_S \approx +0.3$; 82–92% of error draws significant). This is the flux limit: distant hosts are only found if bright. It is the same for low- and high-DM bursts.
- **$M_*$ vs SFR:** a weak main sequence ($+0.25$), identical in every bin (tertiles agree, $p \ge 0.85$).
- **DM$_{\rm host}$ vs SFR:**
  - Absent at $z < 0.29$.
  - Strong at $z > 0.29$ in the catalogue values ($+0.82$), but it drops to $+0.57$ once the scatter in DM$_{\rm cosmic}$ and the SFR errors are propagated.
  - The middle-mass result ($+0.70$) drops to $+0.42$ (significant in 33% of draws).
  - Over the whole sample, AD finds a difference in the DM$_{\rm host}$ distribution between low- and high-SFR hosts ($p = 0.02$) that Spearman does not ($p = 0.10$). Through the errors it is significant in only 9% of draws.
- **Caveat:** the heatmap shows that every other cell is consistent with chance.

**Are the hits noise? (`07_global_tests_both_grids`)** BH treats the tertile tests as independent, but they share hosts and pairs. The permutation test keeps that dependence and asks whether a grid with no real dependence would look like ours.
- **Grid A:** it would not. The high-$z$ DM–SFR tertile ($p = 0.0001$) appears in fewer than 1% of permutations; this family-wise $p$ is 0.008 over all 36 tertiles. Nothing else survives.
  - DM–SFR at middle mass has a family-wise $p$ of 0.12.
  - The three hits among $z$, $M_*$ and SFR within DM$_{\rm host}$ tertiles are the known flux-limit trends.
  - AD supports this only weakly (family-wise $p$ 0.21 for the DM–SFR tertile).
  - The permutation uses catalogue values, so it rules out a look-elsewhere fluke but not the weakening through the errors.
- **Grid B:** it would. See `1_inclination_null/18_global_tests_both_grids` (the same figure).

---

## 3. The hosts are late-type (`3_late_type_hosts/`)

**Claim.** The FRB hosts are overwhelmingly disk galaxies, which justifies comparing them with a disk control and interpreting $q$ as an inclination.

| File | What it shows | Key numbers |
|---|---|---|
| `01_morphology_overview` | (a) Morphological class of all 56 hosts; (b) distribution of the free Sérsic index, split by whether the host is resolved | 24 reliable Sb–Irr, 12 literature disks, 14 disk-like profiles (together 50 of 56), 5 ambiguous, 1 early-type; 80% have $n < 2.5$, median $n = 1.03$ |
| `02_gini_m20` | Gini–$M_{20}$ diagram with the Lotz et al. (2008) boundaries | 52 of 56 in the Sb–Irr region (24 of the 28 reliable ones) |
| `03_table_disk_evidence` | Host-by-host evidence for the 32 hosts with unreliable or non-Sb–Irr Gini–$M_{20}$: free and statmorph Sérsic indices, resolution, literature morphology | 3 of the 4 non-Sb–Irr hosts are spirals in the literature; FRB 20220509G is the only early-type host |
| `04a–d_gallery_*` | Colour image, $r$-band image, GALFIT model and residual for each host in the table | Spiral arms and bars appear in the residuals of the resolved hosts |
| `05_table_inclination_tests_disk_subsets` | The inclination tests rerun without the early-type host, without the ambiguous hosts, and with only the most secure disks | $p_{\rm m}$ 0.49–0.55 (AD) in every subset |
| `06_table_morphology_by_redshift` | Morphology split by redshift | Disk fraction 0.83–1.00 in every bin |

The control sample also does not depend on the morphology cut: removing it entirely gives AD $p_{\rm m} = 0.74$ (`1_inclination_null/05_table_control_variants`).

---

## Where the underlying numbers live

| Topic | Generated by | Output |
|---|---|---|
| Inclination tests and control variants | `scripts/appendix_robustness.py`, `scripts/plot_results_figures.py`, `scripts/power_bhardwaj_deficit.py` | `appendix/robustness/`, `statistics/`, `figures/` |
| Redshift bins | `scripts/redshift_bins.py` | `appendix/redshift_bins/` |
| PSF fit variants | `scripts/run_nopsf_refits.py`, `scripts/appendix_fit_variants.py` | `appendix/fit_variants/` |
| Correlations with $q$ | `scripts/robust_correlations.py` | `appendix/robust_correlations/` |
| Binned correlations (12 + 16 panels) | `scripts/binned_correlations.py` | `appendix/binned_correlations/` (`RESULTS.md` has every number; `SUMMARY.md` is the walkthrough) |
| Global permutation tests of both grids | `scripts/binned_global_tests.py` | `appendix/binned_correlations/GLOBAL_TESTS.md`, `global_tests*.csv` |
| Per-FRB foreground diagnostic | `scripts/plot_foreground_diagnostic.py` | `foreground_halos/foreground_diagnostic.{pdf,png}` |
| Foreground audit and DM budget | `scripts/foreground_audit.py` | `foreground_halos/audit/` (`AUDIT.md`) |
| DM$_{\rm host}$ physics | `scripts/dm_host_physics.py` | `appendix/dm_host_physics/` |
| Foreground halos | `scripts/foreground_halo_dm.py` | `foreground_halos/` (`WALKTHROUGH.md`, `CODE_CHECKS.md`) |
| Halo-mass inversion test (notebook code) | `scripts/mass_inversion_test.py` | `foreground_halos/mass_inversion/` (`MASS_INVERSION.md`) |
| Products of the published foreground method | (pre-correction runs) | `archive/foreground_published_method/` |
| Morphology | `scripts/disk_evidence.py`, `scripts/disk_subsample_tests.py`, `scripts/host_gallery.py` | `appendix/morphology/`, `appendix/robustness/gini_m20_56.csv` |
| Marnoch comparison | `scripts/marnoch_shared_q.py`, `scripts/marnoch_dm_combined.py` | `marnoch_recheck/` |
