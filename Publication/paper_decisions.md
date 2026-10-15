# Paper decisions and evidence map

Updated: 2026-09-04

## Locked for the working draft

- Target journal: ApJL.
- Main-text structure: Introduction; Methodology; Results and Analysis;
  Discussion; Conclusions. Abstract is deferred.
- The opening cites Bhardwaj, Lee & Ji (2024, *Nature*,
  doi:10.1038/s41586-024-08065-w).
- Production base: 64 hosts with published localizations and usable production
  GALFIT workdirs.
- Current paper sample: 50 hosts, defined by $m_r\leq22$, $q=b/a>0.2$, and
  `confirmed=True`.
- Adopted host geometry: component 1 of the production or documented winning
  re-fit.
- Winning-leg composition: 37 production fits and 13 re-fits.
- Inclination conversion:

  $$
  \cos^2 i=\frac{q^2-q_0^2}{1-q_0^2},\qquad q_0=0.2.
  $$

- Uncertainty protocol: formal fit $q$ error and sky-perturbation half-range
  combined in quadrature on $b/a$; Hubble conversion; then a single $5^\circ$
  photometric floor added on $i$ (not mapped into $q$). The obsolete
  $10\times\sigma_q$ rule is not used. The $5^\circ$ term follows Barnes &
  Sellwood (2003) and the Bhardwaj et al. (2024) method-mismatch scale.
- Do not use `pipeline_summary.json:inclination_mc` for the paper uncertainty:
  that field currently samples only the formal GALFIT covariance and does not
  include the sky or $5^\circ$ terms.
- Primary control: 24,450 HSC disk analogues with
  `goodfits_flag=1`, $0.4<n<1.5$, $r\leq22$, and $q>0.2$.
- HSC comparison: 10,000 samples of 50 galaxies without replacement;
  point-estimate KS and AD tests plus sample-size-matched empirical
  $p$-values.
- Verification is a falsification and human-review layer, not an automatic
  threshold gate.

## Facts that differ from the 2024 paper

- The 2024 study used 23 non-repeating CRAFT/ASKAP and DSA-110 FRBs. The
  current production cohort is not restricted to those discovery surveys.
- The current final 50 includes one known repeater, FRB 20220912A.
- The current host sample is selected through published host localizations,
  not a uniform pipeline PATH $P(O)>0.95$ requirement.
- The current analysis uses GALFIT with a per-field PSF rather than AutoProf
  plus Photutils as the two primary estimators.
- HSC is the primary control. Its catalog uses Lenstronomy, not GALFIT, but
  both analyses fit PSF-convolved single-Sérsic profiles with free projected
  axis ratio.

## Decisions required before submission

1. **Magnitude definition.** The current science-cut code uses
   `pipeline_galfit_results.csv:mag`. It gives 60 hosts with $m_r\leq22$.
   The recommended science field, `mag_final`, gives 59 because FRB 20221116A
   has no reliable final magnitude. That object has $q<0.2$, so both choices
   still give 53 after the joint cuts and the same confirmed 50. The manuscript
   must nevertheless define one field consistently.
2. **Radio-survey selection claim.** The master table's `survey` field mixes
   historical imaging and discovery labels, so it cannot yet support a clean
   count by discovery instrument. Audit discovery provenance before claiming
   sensitivity to “various radio surveys.”
3. **Repeater policy.** Either retain FRB 20220912A and state that no
   repeater-status cut is imposed, or remove it and regenerate every
   $N=50$ control draw and statistical test at $N=49$.
4. **Association threshold.** Decide whether the submitted sample requires a
   numerical PATH threshold. The current 64-host cohort instead requires a
   published localization and uses pipeline AstroPath as a cross-check.
5. **Individual inclination uncertainties.** Generate and freeze the
   asymmetric 16th--84th percentile intervals for Table 1 using the same
   truncated-$q$ Monte Carlo as the CDF.
6. **HSC band matching.** The HSC shapes are fitted in $i$, while the control
   magnitude cut uses catalog $r$. Keep this explicit and report the
   bandpass/resolution sensitivity in the Appendix.
7. **Morphological matching.** Confirm that $0.4<n<1.5$ is the final disk
   definition and decide whether an HSC color or star-forming cut is also
   required.
8. **Final statistical hierarchy.** KS and AD should be the headline
   distribution tests; matched-$N$ empirical values are the fair comparison to
   the CDF envelopes. Other tests belong in the Appendix unless needed.

## Current numerical checks

- Production fits: 64.
- Magnitude cut using current `mag`: 60.
- Joint `mag` and $q>0.2$ cut: 53.
- Confirmed within the joint cut: 50.
- Rejected within the joint cut: FRB 20190711A, FRB 20221101B, and
  FRB 20230930A.
- Imaging across all 64: 51 Legacy Survey, 12 Pan-STARRS1, and 1 HSC.
- Imaging across the final 50: 41 Legacy Survey, 8 Pan-STARRS1, and 1 HSC.
- HSC control pool: 24,450.
- Point-estimate medians: FRB $\cos i=0.610$; HSC $\cos i=0.575$.
- Matched-$N$ values currently on disk: KS $p=0.569$ and AD $p=0.574$.

## Appendix / data-health plan

Keep the main Letter focused on the sample, fitting protocol, uncertainty
model, control sample, and headline tests. Move the following checks to the
Appendix:

- local reduced $\chi^2$ and residual-flux fraction;
- residual Fourier $m=2$ test;
- sky $\pm1\sigma_{\rm sky}$ re-fits;
- free-isophote comparison;
- host $q$ versus PSF ellipticity and host--PSF PA alignment;
- magnitude and resolution dependence;
- independent AstroPhot S\'ersic fits;
- alternate imaging and bandpass checks;
- controlled $n=1$, fixed-sky, and contaminating-star re-fits;
- sensitivity to $q_0$, disk-selection window, and magnitude threshold.

## Repository evidence

- Cohorts: `pipeline_scripts/verification/SCIENCE_CUT_AND_COHORT.md`
- Confirmation: `pipeline_scripts/verification/host_confirmation.csv`
- Per-host decisions: `pipeline_scripts/verification/HOST_TRIAGE_CASES.md`
- Verification methods: `pipeline_scripts/verification/FIT_VERIFICATION_CHECKS.md`
- Pipeline operator manual: `pipeline_scripts/README.md`
- Uncertainty literature: `plots/plots_null/v2/frb_vs_hsc_confirmed50/GALFIT_Q_ERROR_LITERATURE.md`
- Confirmed-50 inputs: `plots/plots_null/v2/frb_vs_hsc_confirmed50/confirmed50_q.csv`
- HSC tests: `plots/plots_null/v2/frb_vs_hsc_confirmed50/TEST_RESULTS.md`
- Historical HSC method comparison (use its 2026-09-04 status warning; some
  counts and neighbor-policy text are stale):
  `plots/plots_null/v2/hsc_audit/GALFIT_PIPELINE_VS_KAWINWANICHAKIJ.md`

## Citation keys to verify

- `2024Natur.634.1065B` — Bhardwaj, Lee & Ji (2024)
- `2021ApJ...921...38K` — Kawinwanichakij et al. (2021)
- `2019AJ....157..168D` — DESI Legacy Imaging Surveys
- `2016arXiv161205560C` — Pan-STARRS1 survey
- `Aihara2019` — HSC-PDR2
- `Aggarwal2021` — PATH
- `1996A&AS..117..393B` — Source Extractor
- `2011ASPC..442..435B` — PSFEx
- `2002AJ....124..266P`, `2010AJ....139.2097P` — GALFIT
- `2019MNRAS.483.4140R` — Statmorph
- `1926ApJ....64..321H` — Hubble inclination relation
- `2007ApJS..172..615H` — GALFIT simulation uncertainties
- `2017ApJ...845..114G` — sky-perturbation uncertainty strategy
