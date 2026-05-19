# Project State & Method History

**Last Updated:** May 19, 2026

This file is the canonical, condensed history of the methods that are still valid in the current repository. Overturned experiments and obsolete paths are intentionally omitted. It is structured so that any AI assistant can read it end-to-end and immediately understand the project architecture, active pipelines, file locations, and technical conventions.

---

## 1. Canonical Repository Layout

### Pipeline Scripts (`pipeline_scripts/`)
The autonomous end-to-end processing chain. A single orchestrator (`master_run.py`) drives all phases over a per-run workdir under `pipeline_scripts/Output/<FRB>_<tag>/`; each phase script is independently runnable.

* `pipeline_scripts/master_run.py` — Orchestrator. Stages inputs into `<output>/.workdir/`, writes per-run YAML copies with CLI overrides applied, runs all phases, and copies the requested deliverables out. After Phase 2, writes workdir **`galfit_config.yaml`** with `mag_zeropoint` = `zp_aper_40px` from `zero_points.json` (overridable via `--galfit-zp`). Tool keywords: `catalog`, `psf`, `photometry`, `astropath`, `galfit`, `all` (includes `sky_fit_audit.json`). Supports CLI overrides for the localization ellipse (`--err-a-arcsec`, `--err-b-arcsec`, `--err-theta-deg`, `--p-u`), `--detect-thresh`, `--pixel-scale`, `--seeing-fwhm`, `--gain`, `--mag-mode`, `--target-snr-min`, `--galfit-zp`. `use_weight_map` defaults to `true` and is automatically flipped to `false` when `--invvar` is omitted.
* `pipeline_scripts/run_all_frbs.py` — Batch driver. Iterates `master_frb_localization.csv`, runs `master_run.py` for every FRB with `<FRB>_flux.fits` (+ optional invvar) in `large_cutouts/`, captures each run's stdout to `<output>/master_run.log`, and prints a final per-FRB summary. By default restricts to rows with `coord_semantics='host'`.
* `pipeline_scripts/SExtractor + PSFEx/` — Phase 1: PSF modeling and source extraction. Two-pass SExtractor + PSFEx (iterative `SAMPLE_MINSN` retry 30 → 10). Hard-coded catalog name `image.cat`. Produces `image.cat`, `image.psf`, `proto_image.fits` (25×25 PSF model), `segmentation_map.fits`, PSFEx diagnostics.
* `pipeline_scripts/photometry + astropath/` — Phase 2: PSF re-photometry as the sole source of all downstream mags; three calibrated ZPs (40 px aperture, PSF, Auto) against PS1 (Vizier `II/349`) with **SkyMapper DR1.1 (`II/358`) fallback** when PS1 returns nothing (typical south of Dec ≈ −30°); AstroPath association via a WSL Conda bridge. Star/galaxy cut uses the uncertainty-aware criterion `(SPREAD_MODEL + 3·SPREADERR_MODEL) < 0.005`; calibrated-magnitude sanity filter `[12, 28]` rejects corrupt-flux artefacts before they poison the prior. AstroPath integration grid step is adaptive so the absolute step is `≤ σ_loc / 5`. All prior knobs live in a single labelled block (`# ASTROPATH PRIOR CONFIGURATION`) reproducing the Aggarwal+2021 "adopted" defaults.
* `pipeline_scripts/galfit_fitting/` — Phase 3. Dynamic cutout generation centred on the AstroPath best host with **containment-based** neighbor handling (`frac ≥ 0.95` → fit; `0.50 ≤ frac < 0.95` → grow ROI then fit; `0 < frac < 0.50` → mask), iterated to convergence; `host_components.csv` ordered with the FRB host as row 0 so GALFIT component 1 = host; `host_sigma` anchored to the empirical sky MAD when the invvar-derived sigma disagrees with sky scatter by more than 2×. Phase 3b builds `galfit.feedme` with **`J)` = `zp_aper_40px`** from Phase-2 `zero_points.json` (written to workdir `galfit_config.yaml` by `master_run.py`), PA = `THETA_IMAGE − 90°`, initial mag = `MAG_40PX + mag_zeropoint`, sky QA + retry, `constraints.txt` (`n ∈ [0.5, 6.0]`, `re ∈ [1.5, 100.0]`), runs `wsl galfit`, and produces the 3-panel `galfit_results.png`. Template defaults: `galfit_fitting/galfit_config.yaml`.

### Active Tool Directories (`tools/`)
* `tools/astrophot/` — AstroPhot scripts, notebooks, logs, and result tables.
* `tools/galfit/` — Legacy GALFIT run directories (`tools/galfit/runs/{FRB}/{no_psf_sigma, with_psf_sigma}`), scripts, fit outputs, and IMR PNG exports.
* `tools/statmorph/` — Statmorph notebook/script plus non-parametric results.
* `tools/photutils/` — Photutils no-PSF ellipse fitting scripts, comparison diagnostics, and CDF plots.
* `tools/simulation/` — Mock galaxy generation, recovery benchmarking, and accuracy evaluation.
* `tools/AstroPath/` — Legacy standalone AstroPath visualization scripts (e.g. `plot_r70_legacy.py`).
* `tools/astropath/` — Expansion-set AstroPath candidate extraction and PATH posterior analysis.
* `tools/Photometry/` — Original prototyping environment (`photometry.ipynb`) and production calibration script (`calibrate_photometry.py`).
* `tools/legacy/` — Legacy Survey fetch/comparison utilities.

### Shared Data/Assets (Root)
* `cropped_host_galaxies/` — Flux and sigma FITS cutouts used in the original 23-target fitting.
* `large_cutouts/` — 10-arcminute Legacy-style flux + inverse-variance cutouts. **Contents are restricted to** `{FRB}_flux.fits` **and** `{FRB}_invvar.fits` **only** (53 FRB pairs on disk as of May 2026); intermediate PSFEx/SExtractor artifacts were removed from this folder.
* `psfs/` — PSF images, including `PSFEx + SExtractor/final_center_psfs/` with 52 normalized 25×25 center-extracted models.
* `plots/` — Generated figures: CDF comparisons, multiband nulls, photutils diagnostics, etc.

### Archive & Reports
* `Archive/csv/` — Archived CSV products from earlier stages (includes e.g. `n1_comparison_results.txt`).
* `Archive/scripts_retired/` — One-off migration scripts kept for reproducibility only (e.g. `integrate_12frbs_legacy.py`).
* `Archive/astropath_diagnostics/` — AstroPath diagnostic PNGs and CSV moved from repo root.
* `Archive/original paper/` — Original paper assets.
* `Reports/00 Galfit (AS) verification/` — AstroPhot vs GALFIT sigma-weighted validation report.
* `Reports/01 PSFEx_Comparison_Report/` — PSFEx + Legacy Tractor comparison report.
* `Reports/02 PhotUtils + Legacy/` — Photutils + Legacy comparison mini-report and Free-n vs Fixed-n tables.
* `Reports/pipeline_galfit_review_handoff.md` — Pipeline GALFIT QA handoff (Tier A/B flags, ZP offset context, per-FRB checklist). Companion CSV: `pipeline_unphysical_fits_review.csv`.

---

## 2. Active Pipeline Architecture

The current production pipeline processes a single FRB target through these sequential stages:

### Stage 1: SExtractor + PSFEx (`pipeline_scripts/SExtractor + PSFEx/`)
* **Input:** `image.fits`, `invvar.fits` (copied from `large_cutouts/`).
* **Process:** 2-pass SExtractor (detection → PSFEx → fit). Constant PSF (Degree 0). Iterative SNR thresholding (30 → 2) to achieve ≥10 accepted stars.
* **Output:** `image.cat` (FITS_LDAC with `NUMBER`, positions, magnitudes including multi-aperture `MAG_APER`), `image.psf`, `proto_image.fits` (25×25 PSF prototype), `segmentation_map.fits`.
* **Script:** `run_psf_pipeline.py`

### Stage 2: Photometry + AstroPath (`pipeline_scripts/photometry + astropath/`)
* **Input:** PSF products from Stage 1 + `image.fits` (+ optional `invvar.fits`).
* **Process:** Second SExtractor pass with `PSF_NAME image.psf` produces `image.psf.cat` (the **single source of truth** for all downstream magnitudes). Three ZPs (40 px aperture, PSF, Auto) calibrated against PS1 (Vizier `II/349`); falls back to SkyMapper DR1.1 (`II/358`) when PS1 returns nothing. Point-source rejection for the host-candidate set uses `(SPREAD_MODEL + 3·SPREADERR_MODEL) < 0.005`; a calibrated-mag sanity filter `[12, 28]` rejects corrupt-flux artefacts. AstroPath association via a WSL Conda bridge with an adaptive integration grid step (`≤ σ_loc / 5` in absolute arcsec). WSL bridge errors are propagated (`sys.exit(1)`) so failures are not silently masked.
* **Output:** `calibrated_photometry_results.csv` (RA/Dec, three calibrated mags, FLUX_RADIUS, SPREAD_MODEL, AstroPath inclusion flag), `zero_points.json` (3 ZPs + N_stars + reference catalog id), `astropath_posteriors.csv`, `astropath_association.png` (display stretch computed on the 1 ′ zoom region), `image.psf.cat`, `image.homo.fits`.
* **Script:** `run_photometry_astropath.py`
* **Config:** `photometry_astropath_config.yaml` + an in-code `# ASTROPATH PRIOR CONFIGURATION` block (Aggarwal+2021 "adopted" defaults: `THETA_PDF='exp'`, `THETA_MAX=6.0`, `THETA_SCALE=1.0`).

### Stage 3: GALFIT Fitting (`pipeline_scripts/galfit_fitting/`)
* **Input:** `image.fits` (+ optional `invvar.fits`), `image.cat`, `segmentation_map.fits`, `proto_image.fits`, FRB RA/Dec. Phase 3a auto-overrides the target with the AstroPath best host (highest `posterior_O` ≥ `--min-astropath-posterior`) when `astropath_posteriors.csv` is present.
* **Process (Cutout Generation — `generate_galfit_cutouts.py`):**
    * Starts from `ROI = host_seg_bbox + --host-pad` (default 20 px). For every other seg island touching the ROI, computes `frac = pixels_in_ROI / total_pixels`:
        * `frac ≥ --contain-thresh` (default 0.95) → **fully contained** → fit as a Sérsic (or mask if `CLASS_STAR ≥ --neighbor-class-star-max`).
        * `--expand-thresh` (0.50) ≤ `frac < 0.95` → **largely filled** → grow ROI to include the source, then fit (or mask if stellar).
        * `0 < frac < 0.50` → **fringe** → mask in-frame pixels only; ROI does **not** grow.
        * `frac == 0` → out of frame, ignored.
    * The expand-and-recategorise loop iterates up to `--max-roi-iterations` (default 6). The FRB host is unconditionally added to the fit set first and never masked.
    * `host_components.csv` is ordered with the FRB host as **row 0** so GALFIT instantiates it as component 1.
    * `host_sigma.fits = 1/√invvar` with `invvar ≤ 0` or non-finite pixels mapped to `σ = 1e30` (also flagged in `host_mask.fits`). The absolute scale is sanity-checked against the empirical sky noise (robust MAD × 1.4826 over unmasked pixels); when the invvar-derived sigma disagrees with the sky scatter by more than 2×, the sigma cutout is multiplied by a single global factor `k = σ_sky / σ_invvar`. This preserves the spatial structure of the invvar map while pinning its absolute scale to the data — handles LS frames delivered with flux and invvar on inconsistent unit conventions. `k` is logged per FRB.
    * Outputs: `host_cutout.fits`, `host_sigma.fits`, `host_mask.fits`, `host_components.csv`, `qa_cutout_mask.png`.
* **Process (GALFIT Execution — `run_galfit_fitting.py`):**
    * Aborts (non-zero exit) if `proto_image.fits` is missing.
    * Resolves **`mag_zeropoint`** from workdir `galfit_config.yaml` → Phase-2 `zero_points.json` `zp_aper_40px` → fallback 22.5. Reported `mag` in `fit.log` is on the `J)` system; initial `mag_init = MAG_40PX + mag_zeropoint` must use the same ZP or the optimizer can start on the wrong flux scale (structural parameters can still converge badly even when `re`/`n`/`b/a` look reasonable).
    * **Sky QA:** seeds global sky from SExtractor `BACKGROUND` on host row 0 in `host_components.csv` (ADU). After pass 1, compares `fit.log` sky to that reference via `scripts/galfit_fitlog_parse.parse_fitlog_sky_level`; if `|Δ| > sky_tolerance_adu` (default **3 ADU**), reruns with `constraints.txt` line `{sky_comp} 1 -3 3` (soft band around the seed). Writes `sky_fit_audit.json`; exits **1** if QA still fails after retry (`sky_check_enabled`, `sky_max_retries` in `galfit_config.yaml`).
    * For each row in `host_components.csv`: emits a Sérsic component seeded with `XC_CUTOUT`, `YC_CUTOUT`, `1/ELONGATION` (axis ratio), and `pa = THETA_IMAGE − 90°` (SExtractor → GALFIT PA convention). Initial mag = `MAG_40PX + mag_zeropoint` (since SExtractor runs with `MAG_ZEROPOINT=0`, this shifts the seed onto GALFIT's `J)` system so the starting model flux is on the same scale as the aperture flux).
    * Per-component `constraints.txt`: `n ∈ [0.5, 6.0]`, `re ∈ [1.5, 100.0]`; optional sky band on retry.
    * Executes `wsl galfit galfit.feedme` (up to two passes when sky QA triggers).
    * Parses `out.fits` (HDU 1 = image, 2 = model, 3 = residual) and renders the 3-panel `galfit_results.png`.
* **GALFIT also produces:** `fit.log`, `galfit.01`.
* **Config:** `pipeline_scripts/galfit_fitting/galfit_config.yaml` (defaults: sky QA toggles, plate scale); per-run copy in workdir with `mag_zeropoint` from `zp_aper_40px`.

---

## 3. Earlier Method Pipeline (Original 23-Target Sample)

### A. Data Preparation
* Verified FRB coordinates and crop regions, then generated per-target cutouts in `cropped_host_galaxies/`.
* Produced/collected PSFs and matched them to each cropped host image.
* Used inverse-variance survey products to derive sigma maps for physically weighted fits.

### B. GALFIT (Sigma-Weighted, Two-Stage)
* Ran sigma-weighted workflows with PSF-inclusive runs in `tools/galfit/runs/`.
* Compiled results in `galfit_sigma_metrics_summary.csv` (root).
* Free-n baseline extracted via `scripts/compile_true_freen_legacy_comparison.py` from intermediate `galfit.0X` files → `galfit_vs_legacy_quick_read.csv`.

### C. AstroPhot (PSF + Sigma)
* Used per-target sigma maps converted to variance. Optional fixed-Sersic mode (`--fixed-n`).
* Canonical output: `tools/astrophot/results/astrophot_psf_sigma_n1_inclination_angles.csv`.

### D. Statmorph Non-Parametric Validation
* Computed CAS, Gini, M20 for all FRB hosts.
* Canonical output: `tools/statmorph/statmorph_nonparam_results.csv`.

### E. Cross-Tool Validation
* AstroPhot vs GALFIT consistency summaries and CDF bias comparisons.
* Statmorph metrics used to check for severe morphology outliers.

---

## 4. Canonical Data Products

### Master Catalogs
* `master_frb_summary.csv` — 55 FRB host candidates, standardized coordinates (decimal and HMS/DMS); legacy wide table including initials, partial GALFIT columns, and notes.
* `master_frb_localization.csv` — **Localization truth table** built from the summary: sky position, HMS/DMS, error ellipse, `z`, `DM`, `status`, `coord_semantics` (host vs burst position), crop box (`xmin`…`ymax` when present), and **`survey`** (`LS`, `PS1`, `LS+PS1`, `unsure`). Generator: `scripts/build_master_frb_localization.py`.
* `master_frb_galfit_from_logs.csv` — **GALFIT parameters from `fit.log`** under `tools/galfit/runs/<FRB>/{no_psf_sigma,with_psf_sigma}/`, including inclinations and parse provenance. Generator: `scripts/build_master_frb_galfit_from_logs.py`; audit: `scripts/audit_master_frb_galfit_csv.py`. Shared parser module: `scripts/galfit_fitlog_parse.py`.
* `galfit_sigma_metrics_summary.csv` — Canonical GALFIT comparison baseline (original 23).
* `galfit_vs_legacy_quick_read.csv` — Free-n GALFIT vs Legacy Survey morphology.
* `galfit_vs_legacy_master.csv` — Full comparison with error metrics.
* `new_16_frbs_galfit_results.csv` — Unified 39-target registry of physical morphology (23 old + 16 new, free-n extraction via second-to-last `fit.log` block parsing with asterisk flag sanitization).

### Pipeline Benchmarking Outputs (root)
* `pipeline_galfit_results.csv` — With-PSF Sérsic parameters parsed from every `pipeline_scripts/Output/<FRB>_all/fit.log` (host = GALFIT component 1 / parser index 0). Columns include **`n_sersic_components`** (rows in `host_components.csv`, sky excluded) and **`compare_ok`** (`True` only when exactly one Sérsic is fitted). **`mag`** retained for reference. Generator: `scripts/compare_pipeline_galfit_vs_master.py`.
* `pipeline_vs_master_galfit_diff.csv` — Join against `master_frb_galfit_from_logs.csv` with `_pipeline`, `_master`, `_delta` for **`chi2nu`, `re`, `n`, `b/a`, `pa`, `inc` only** — **magnitude/flux deltas omitted** because pipeline uses per-field `zp_aper_40px` while legacy master runs use mixed `J)` systems. Carries `n_sersic_components` and `compare_ok`; summary statistics use **`compare_ok=True`** rows only (single-Sérsic deblends). `scripts/analyze_pipeline_vs_master_diff.py` for detailed stats. Three benchmark outliers (`20171020A`, `20220509G`, `20240210A`) excluded from written CSVs by default.
* `pipeline_unphysical_fits_review.csv` — Automated QA flags (heuristic thresholds; many Tier-A cases are false positives after ZP+sky fix — see handoff). Generator: `scripts/flag_pipeline_unphysical_fits.py`. Human workflow: `Reports/pipeline_galfit_review_handoff.md`.
* Per-FRB batch logs: `pipeline_scripts/Output/<FRB>_all/master_run.log` (from `run_all_frbs.py`).

### Cross-Catalog Comparison
* `legacy_vs_galfit_inclination_comparison.csv` — Tractor shape-based inclinations with uncertainty propagation from `shape_e1_ivar`/`shape_e2_ivar`.
* `legacy_vs_galfit_two_inclinations.csv` — Compact collaborator-facing table with `galfit_inc_psf_deg`, `ls_inc_deg`, `ls_inc_err_deg`, `type_ls`, `sersic_n_fit`.
* `photutils_vs_galfit_comparison.csv` — Root-level canonical Photutils comparison.

### AstroPath
* `tools/astropath/results/astropath_expansion_summary.csv` — Expansion-set posterior summary.
* `tools/astropath/results/successful_associations.csv` — High-confidence (P > 0.8) host associations.
* `tools/astropath/data/*_candidates.csv` — Per-FRB candidate tables.

### Plots & Reports
* `plots/plots_astrophot_psf_sigma_n1/CDF_bias_comparison_psf_sigma_n1_fixed.png`
* `plots/plots_multiband_cdf/CDF_bias_multiband_rgb_bands_mc_policy.png`
* `plots/plots_photutils/` — Photutils CDF and delta-inclination diagnostics.
* `Reports/00 Galfit (AS) verification/astrophot_galfit_sigma_validation_2026-03-18.pdf`
* `Reports/01 PSFEx_Comparison_Report/psfex_galfit_report.pdf`
* `Reports/02 PhotUtils + Legacy/Mini Report.pdf`

### PSF Assets
* `psfs/PSFEx + SExtractor/final_center_psfs/` — 52 normalized 25×25 center-extracted PSF models.
* `psfs/PSFEx + SExtractor/final_center_psfs/final_center_psf_manifest.csv` — Processing manifest.

---

## 5. Key Technical Conventions

### Inclination Formula
All analysis streams use the Hubble inclination formula with intrinsic thickness $q_0 = 0.2$:
$$\cos^2 i = \frac{q^2 - q_0^2}{1 - q_0^2}$$

### Legacy Survey Tractor Ellipticity Conversion
* $|e| = \sqrt{e_1^2 + e_2^2}$, then $q = (1 - |e|) / (1 + |e|)$.
* Uncertainty propagated via Monte Carlo draws in $(e_1, e_2)$ space using `shape_e1_ivar`, `shape_e2_ivar`.

### Null catalogs (v1, May 2026)
* **Canonical outputs (repo root):** `LS_catalog_v1_allsky_modelmr.csv`, `SDSS_catalog_v1_allsky_modelmr.csv`. v0 files live in `Archive/csv/null_catalogs_v0/`.
* **Footprint:** joint Legacy DR10 ∩ SDSS DR16 imaging via conservative Dec cut **−30° ≤ Dec ≤ +90°**; random sampling (SDSS: chunked SQL + shuffle; Legacy: TAP `TOP` with shuffle fallback when `ORDER BY RANDOM()` is unsupported).
* **Magnitudes:** compare **model \(r\)** — Legacy `tractor_mag_r` = `22.5 − 2.5 log10(flux_r)` (nanomaggies); SDSS `rmag` = `cmodelMag_r`. Do **not** use Legacy `petroMag_r` (v0 misname).
* **Axis ratio:** `expAB_r` = \(b/a\) from Tractor \(e_1,e_2\) (Legacy) or PhotoObj `expAB_r` (SDSS).
* **CDF modes:** `strict` — require \(q > q_0\) before building \(\cos i\) pools; `inclusive` — finite \(q \in [0,1]\), with \(q \le q_0\) mapped to \(\cos i = 0\) (face-on edge). CLI: `--sample-mode`, `--mag-column`, `--tag`.
* **Sanity:** `scripts/test_null_catalog_sanity.py` asserts ≥10k galaxies per survey after strict cuts.
* **Plots:** `plots/plots_legacy_cdf/v1_null_plots/` (filename tags e.g. `v1_allsky_modelmr_strict`, `_inclusive`, `_bestmodelba`, `_withrex`).

### CDF Plotting Policy
* Monte Carlo error propagation for FRB host curves when `b/a` uncertainties are available.
* Show 68% CI for single null/reference line plots; suppress CI for multi-line overlays.

### Pixel Scale & Seeing (r70 Reference Field)
* Pixel scale: 0.186856 arcsec/pixel (from CD matrix determinant). Default for Legacy cutouts: 0.262 arcsec/pixel.
* r70 seeing: 0.63″ (3.38 pixels).

---

## 6. Chronological Milestones

### March 2026 — Core Buildout & Weighting
* End-to-end FRB host workflow: coordinate validation, cutouts, PSFs, initial parametric fitting.
* Integrated invvar/sigma-based weighting. Standardized sigma-weighted GALFIT outputs.
* Ran AstroPhot with PSF + sigma; fixed `n=1` mode for method alignment.
* Repository consolidation: tools into `tools/`, archives into `Archive/`.
* GALFIT report corrections (display flip, `20220207C` re-initialization with `R_e=7.8`, `b/a=0.8`).
* Verified SDSS CDF-null logic. Legacy Tractor inclination comparison with uncertainty propagation.
* Multi-band null CDF plotting with MC policy controls.
* Photutils no-PSF ellipse fitting (23/23), median isophotal averaging (Jedrzejewski 1987), nclip=3.
* Parameter recovery simulation suite: 140 synthetic galaxies, GALFIT RMSE ≈ 5° (with good init), Photutils RMSE ≈ 8° (robust baseline).
* Free-n Legacy alignment: RMSE ≈ 16.7° vs 19.8° for fixed-n.
* Reports finalized: `00 Galfit verification`, `01 PSFEx Comparison`, `02 PhotUtils + Legacy`.

### April 2026 — Expansion Set & Autonomous Pipelines
* Merged 32 new targets into `master_frb_summary.csv` (total: 55 FRBs).
* Deployed `scripts/fetch_large_cutouts_v2.py`: 10-arcmin FOV, DR10 → DR9 → PS1 priority.
* Achieved 52/55 PSFEx models. Excluded: `20191228A` (unsuitable), `20230718A`/`20230731A` (no coverage).
* Extracted 52 center-localized 25×25 PSFs to `psfs/PSFEx + SExtractor/final_center_psfs/`.
* Built production photometry pipeline (`tools/Photometry/calibrate_photometry.py`): dual ZP, 40px large-aperture recovery, FLAGS monitoring.
* Integrated AstroPath expansion-set workflow: 29/31 successful catalog pulls, high-confidence associations tracked.
* Developed autonomous `pipeline_scripts/` chain: SExtractor+PSFEx → Photometry+AstroPath → GALFIT Fitting.
* Decoupled `run_psf_pipeline.py` to output only `proto_image.fits` + `.psf` + `.cat` + `segmentation_map.fits` (stripped residual/chi/moffat/snap outputs).
* Built segmentation-aware dynamic cutout generator with ≥20% neighbor inclusion / <20% masking logic.
* Implemented autonomous `.feedme` builder with dynamic zero-point triangulation from the photometry pipeline's calibrated magnitudes, Sersic constraints (`n ∈ [0.5, 6]`, `re ∈ [1.5, 100]`), and `wsl galfit` execution with 3-panel QA output.
* Consolidated free-n results for 39 targets via second-to-last `fit.log` block parsing.
* **Pipeline Hardening & r70 Analysis:**
    *   Transitioned to "weight-map agnostic" logic (`use_weight_map` YAML toggle) to handle single-image runs without `invvar.fits`.
    *   Resolved `MAG_AUTO` calibration bias by implementing a dedicated `zp_auto_med` derived from Kron-aperture comparison against PS1.
    *   Automated PSF seeing measurement: The pipeline now computes FWHM from `proto_image.fits` moments and injects it into SExtractor's `SEEING_FWHM` parameter. This ensures the `CLASS_STAR` morphology classifier is correctly calibrated to the local PSF.
    *   Refined candidate selection: Discovered that `FLAGS_MODEL == 0` was dropping valid but complex galaxies; restricted the flag filter to ZP calibration stars only.
    *   Retired a short-lived **`classstar`** directory trial in favor of **PSFEx spread-model** / PSF-aware stellar classification (better than raw `CLASS_STAR` alone for compact galaxies mislabelled in DR10).
    *   Executed r70 analysis: σ = 1" run successfully localized host with P(O) = 0.83; subsequent σ = 0.5" and `classstar` runs ongoing to refine "missing" galaxy detection.

### May 2026 — Master tables, repository hygiene, documentation

* Added **`master_frb_localization.csv`** and **`master_frb_galfit_from_logs.csv`** as reproducible, script-generated counterparts to fields scattered in `master_frb_summary.csv` (see §4). GALFIT columns use robust **`fit.log`** parsing (`scripts/galfit_fitlog_parse.py`): sane-χ² iteration selection and skipping final fixed-`n` refinement when bracketed `n` appears on the last sane iteration.
* **`scripts/compile_galfit_logs.py`** updated to use **`tools/galfit/runs`** and the shared parser.
* **`large_cutouts/`** cleaned to **only** `{FRB}_flux.fits` / `{FRB}_invvar.fits` (53 FRB pairs); removed stray PSFEx defaults, segmentation maps, and per-target diagnostic FITS from that folder.
* Root clutter removed (duplicate SExtractor/photometry defaults, AstroPath diagnostic PNGs/CSV → **`Archive/astropath_diagnostics/`**). **`tmp/`** directory removed—use **`scripts/`** for small utilities only.
* Obsolete **`scripts/`** entries deleted (`_add_inc_to_results.py`, `extract_aper13_photometry.py`, `condense_csvs.py`); **`integrate_12frbs_legacy.py`** moved to **`Archive/scripts_retired/`**.
* **`scripts.md`** rewritten as a complete index of **`scripts/`** plus pipeline/tools cross-references.

### May 2026 — Pipeline orchestration, batch validation, GALFIT hardening

* **Top-level orchestrator `pipeline_scripts/master_run.py`** added. Stages inputs and a per-run YAML (with CLI overrides) into `<output>/.workdir/`, runs all phases, and copies only the requested deliverables out. Tool keywords: `catalog`, `psf`, `photometry`, `astropath`, `galfit`, `all`. Exposes CLI overrides for the AstroPath ellipse (`--err-a-arcsec`, `--err-b-arcsec`, `--err-theta-deg`, `--p-u`), `--detect-thresh`, `--pixel-scale`, `--seeing-fwhm`, `--gain`, `--mag-mode`, `--target-snr-min`, `--galfit-zp`, and `--keep-workdir`.
* **Batch driver `pipeline_scripts/run_all_frbs.py`** added. Iterates `master_frb_localization.csv`, runs `master_run.py` for every FRB with `<FRB>_flux.fits` (+ optional invvar) in `large_cutouts/`, captures stdout to `<output>/master_run.log`, and prints a final per-FRB summary. Restricts by default to `coord_semantics='host'`.
* **Photometric-ZP fallback (Phase 2)**: PS1 (Vizier `II/349`) → **SkyMapper DR1.1 (`II/358`)** when PS1 returns nothing (typical for fields south of Dec ≈ −30°). The catalog actually used is recorded in `zero_points.json` under `reference_catalog`. WSL bridge errors now `sys.exit(1)` so Phase-2 failures are not masked.
* **AstroPath robustness (Phase 2)**:
  * Point-source rejection switched to the uncertainty-aware criterion `(SPREAD_MODEL + 3·SPREADERR_MODEL) < 0.005` so faint / low-SNR galaxies survive.
  * Calibrated-magnitude sanity filter (`[12, 28]`) drops corrupt-flux artefacts (PSF flux ≤ 0 → `mag ≈ 131`) that otherwise blow up `driver_sigma` and swamp real hosts.
  * Integration grid step set **adaptively** so the absolute step in arcsec stays `≤ σ_loc / 5` for every candidate (handles tight localizations with extended candidates).
  * All prior knobs live in one labelled in-code block (`# ASTROPATH PRIOR CONFIGURATION`) reproducing the Aggarwal+2021 "adopted" defaults.
  * `astropath_association.png` stretch is computed on the 1 ′ zoom around the FRB, not the full cutout, so a bright off-FRB source cannot wash the panel to black.
* **GALFIT cutout policy (Phase 3a)** rewritten as a **containment-based** triage: `frac ≥ 0.95` → fit; `0.50 ≤ frac < 0.95` → grow ROI to fully contain, then fit; `0 < frac < 0.50` → mask in-frame pixels only; with iterative ROI expansion. The FRB host is unconditionally added to the fit set (and never masked, even with high `CLASS_STAR`). `host_components.csv` is written with the host as row 0 so GALFIT component 1 = host. The shared `galfit_fitlog_parse.py` parser accepts `sersic_component_index` (default 0) so the comparison tooling reads the host every time.
* **`host_sigma` scale anchoring (Phase 3a)**: when the invvar-derived sigma disagrees with the empirical sky scatter (robust MAD × 1.4826) by more than 2× in either direction, the cutout sigma is multiplied by a single global factor `k = σ_sky / σ_invvar`. This preserves the spatial structure of the invvar map while pinning its absolute scale to the data, fixing a class of Legacy Surveys frames delivered with inconsistent flux / invvar unit conventions (which would otherwise inflate `χ²/ν` by many orders of magnitude). `k` is logged per FRB.
* **GALFIT initial-mag ZP shift (Phase 3b)**: initial magnitude written to the feedme is `MAG_40PX + mag_zeropoint`. SExtractor in this pipeline runs with `MAG_ZEROPOINT=0`, so `MAG_40PX` is raw `−2.5·log10(flux_ADU)`; shifting onto GALFIT's `J)` ZP keeps the starting model flux on the same scale as the aperture flux and avoids the Levenberg-Marquardt descent into "thin-core, peaked-SB" minima.
* **Phase-3b PA convention** finalised as `pa = THETA_IMAGE − 90°` (SExtractor +x → GALFIT +y).
* **Phase-3b hard requirement**: aborts with non-zero exit if `proto_image.fits` is missing rather than silently producing an unconvolved fit.
* **Phase-1 + 2 `detect_thresh`** lowered to **3** (from 10) so faint hosts make it into `image.cat` / `segmentation_map.fits`; PSFEx still filters cleanly via `SAMPLE_MINSN=30`.
* **`use_weight_map`** defaults to `true`. `master_run.py` writes `use_weight_map: <invvar_provided>` into the workdir copy automatically, so the YAML key only needs touching when running a phase script outside the orchestrator.
* **Pipeline-vs-master benchmarking**: `scripts/compare_pipeline_galfit_vs_master.py` parses every pipeline `fit.log` (with-PSF Sérsic, host component) and joins against `master_frb_galfit_from_logs.csv`, writing `pipeline_galfit_results.csv` and `pipeline_vs_master_galfit_diff.csv` to the repo root; `scripts/analyze_pipeline_vs_master_diff.py` produces summary statistics (median, MAD, RMSE, quantiles, Spearman correlations, top-N deviations). **May 19 update:** diff excludes `mag`; adds `n_sersic_components` / `compare_ok`; filters summaries to single-Sérsic fits.

### May 2026 — Pipeline GALFIT QA & magnitude-system alignment

* **Master vs pipeline magnitude offset diagnosed**: apparent ~−2.4 mag median `mag_delta` was not a flux-scaling bug in GALFIT but a **zeropoint mismatch** — pipeline always `J) 22.5`; legacy master (23 hosts) used `J) 25.0` from `scripts/run_galfit_with_sigma.py`. `compare_pipeline_galfit_vs_master.py` now applies a **−2.5 mag correction to `mag_master`** for those feedmes before `mag_delta` (documented in module docstring; no extra CSV column).
* **Unphysical pipeline fits audit**: `scripts/flag_pipeline_unphysical_fits.py` → `pipeline_unphysical_fits_review.csv`. **11 Tier A** degenerate cases (bright `mag` ≈ 4–10, `re` → 100 px, bad sky, multi-Sérsic confusion, or huge χ²/ν); **26 Tier B** suspect (face-on `b/a` floor, `n=6` ceiling, `last_block_first_sersic` parser). **19 FRBs** pass strict automated cuts (`17 ≤ mag ≤ 24`, `re < 50`, `chi2nu ≤ 3`, sane parser, `b/a > 0.15`) pending visual QA.
* **Handoff artifact**: `Reports/pipeline_galfit_review_handoff.md` — tomorrow checklist, Tier tables, failure modes, ZP context, links to `pipeline_scripts/Output/<FRB>_all/` QA files.
* **Sky seed + QA retry (Phase 3b)**: `run_galfit_fitting.py` seeds sky from SExtractor `BACKGROUND`, retries with ±3 ADU constraint when `fit.log` sky drifts; outputs `sky_fit_audit.json`.
* **Null catalogs v1** replace v0 COSMOS / single-brick mismatch; see §5 “Null catalogs (v1)”. Older drivers (`generate_multiband_cdf_null_plot.py`, `generate_sigma_plots.py`, `generate_all_plots.py`) still default to archived v0 paths until a follow-up pass.

### May 19, 2026 — Phase 3b ZP + sky hardening, full batch re-run, pipeline QA review

* **Phase 3b photometric ZP wired to Phase 2:** `master_run.py` calls `write_galfit_config()` after Phase 2, setting `mag_zeropoint` from `zero_points.json` → **`zp_aper_40px`** (40 px aperture calibration vs PS1/SkyMapper). `run_galfit_fitting.py` resolves ZP as: workdir `galfit_config.yaml` → `zero_points.json` → fallback 22.5. Initial feedme magnitude remains `MAG_40PX + mag_zeropoint`. Template: `pipeline_scripts/galfit_fitting/galfit_config.yaml`.
* **Sky QA (unchanged logic, now always active via config):** SExtractor `BACKGROUND` seed; pass-2 retry with ±3 ADU constraint; `sky_fit_audit.json`; non-zero exit if QA fails after retry.
* **Utility `scripts/rerun_pipeline_galfit_phase3b.py`:** Re-run Phase 3b only on existing `Output/*_all/` folders (writes `galfit_config.yaml` from `zp_aper_40px` first). Used to refresh all fits after ZP/sky changes without re-running Phases 1–2.
* **Full end-to-end batch:** `python pipeline_scripts/run_all_frbs.py` on all 53 cutouts in `large_cutouts/` — **39 OK** (host `coord_semantics`), **0 fail**, **14 skip** (signal coordinates; use `--include-signal` to run those), **~22 min** wall time (May 2026). Stdout captured per FRB in `Output/<FRB>_all/master_run.log`. Confirms orchestrator, WSL bridges, σ rescaling (e.g. `20171020A` logs `k ≈ 9.4×10³`), and per-field ZP propagation are intact.
* **Manual GALFIT QA review (human):** Original **11 Tier A** cases from pre-ZP pipeline were re-inspected. Most were **automated false positives** (crude `mag < 15`, `sky > 100 ADU`, or `|Δsky| > 3` by 0.003 ADU). After review, all original Tier A hosts are **accepted for science use** with noted caveats: multi-component stamps use **component 1 = host**; spirals / crowded fields may have high χ²/ν or ugly residuals without invalidating host `re`, `b/a`, `PA`. See updated `Reports/pipeline_galfit_review_handoff.md`.
* **Comparison tooling refresh:** `scripts/compare_pipeline_galfit_vs_master.py` — added `n_sersic_components`, `compare_ok`; **removed `mag` from diff deltas**; summary stats on single-Sérsic subset only (**23 of 36** matched FRBs as of batch re-run). `galfit_fitlog_parse.count_fitted_sersic_components()` added. `flag_pipeline_unphysical_fits.py` reads `sky_fit_audit.json` for `sky_qa_failed`.
* **Regenerated:** `pipeline_galfit_results.csv`, `pipeline_vs_master_galfit_diff.csv` after batch.

**Single-Sérsic comparison snapshot (May 19, post-rerun):** median |Δ| roughly `re` 0.16 px, `n` 0.13, `b/a` 0.01, `inc` 0.82°, `pa` 0.69°. Large `chi2nu` deltas often reflect legacy master χ² blow-ups (invvar scaling), not pipeline structural failure — treat pipeline χ²/ν ≈ O(1) as the meaningful metric when σ anchoring ran.

---

## 7. Technical Runbook: Command-Line Automation Patterns

### SExtractor & PSFEx (Windows → WSL)
* Invoked through Python `subprocess` calls targeting `wsl source-extractor` and `wsl psfex`.
* Typical: `wsl source-extractor <image> -c <config> -CATALOG_NAME <cat>` and `wsl psfex <cat> -c <config>`.
* File paths remain Windows-native; WSL resolves them automatically via `/mnt/c/` when needed.

### AstroPath Engine (Linux WSL + Conda)
* `astropath_pkg` depends on Linux-only libraries (`healpy`, `healpix`).
* Execution from Windows orchestrator: `wsl -e bash -ic 'conda activate frb_project && python <script> <args>'`.
* Convert absolute `C:\` paths to `/mnt/c/...` when passing to WSL commands.

### GALFIT (WSL)
* Executed via `wsl galfit galfit.feedme` from the `galfit_fitting/` directory.
* Outputs `out.fits` (3-HDU: image/model/residual), `fit.log`, `galfit.01`.

### GALFIT Cutout Isolation Logic
* SExtractor `CHECKIMAGE_TYPE SEGMENTATION` defines galaxy boundaries.
* **Containment-based deblending:** neighbors whose seg-island is `≥ contain-thresh` (default 0.95) contained inside the host ROI are fit as additional Sérsic components; neighbors with `expand-thresh ≤ frac < contain-thresh` (default 0.50–0.95) trigger an ROI expansion until they are fully contained, then are fit; "fringe" islands with `frac < 0.50` are masked in `host_mask.fits` only and the ROI is left alone. The FRB host is unconditionally added to the fit set.

---

## 8. Guidance for Future Runs

* Treat paths in this file as authoritative. Verify against the filesystem if in doubt.
* Each `pipeline_scripts/` subdirectory operates on local `image.fits` and `invvar.fits` — do not reference files from other directories.
* Prefer sigma-weighted outputs and fixed-`n` AstroPhot mode when reproducing the professor-facing comparison workflow.
* Keep new outputs under the corresponding `tools/<domain>/results` or `plots/` location. Do not use a repo-root `tmp/` folder (removed); put transient logs under `plots/`, `tools/.../logs`, or user-specific paths outside the tracked tree if needed.
* Functional index of **`scripts/`**: see **`scripts.md`**.
* Report-writing conventions are codified in `.agents/skills/report-writing/SKILL.md`.
* MiKTeX is available at `C:\Users\lenovo\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe` for local PDF compilation.
