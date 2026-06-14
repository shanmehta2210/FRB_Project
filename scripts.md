# Scripts Documentation

This file indexes the **active** analysis and utility code. The canonical end-to-end pipeline lives under `pipeline_scripts/`; this document also lists every Python module in the repository-root **`scripts/`** directory (as of May 18, 2026). Transient scratch files belong under `scripts/` or a user-local path — there is no repo-root `tmp/`.

---

## Autonomous Pipeline (`pipeline_scripts/`)

A single orchestrator (`master_run.py`) chains the three phases together over a per-run workdir. The phase scripts are independently runnable.

- **`pipeline_scripts/master_run.py`** — Orchestrator. Stages inputs, applies CLI YAML overrides, runs all phases, copies deliverables to `Output/<FRB>_<tag>/`. After Phase 2, writes **`galfit_config.yaml`** with `mag_zeropoint` = `zp_aper_40px` from `zero_points.json` (override: `--galfit-zp`). Tool keywords: `catalog`, `psf`, `photometry`, `astropath`, `galfit`, `all` (includes `sky_fit_audit.json`). Optional CLI: localization ellipse, `--detect-thresh`, `--pixel-scale`, `--seeing-fwhm`, `--gain`, `--mag-mode`, `--target-snr-min`, `--keep-workdir`. `use_weight_map` auto from `--invvar`.
- **`pipeline_scripts/run_all_frbs.py`** — Batch driver for `master_run.py` over `large_cutouts/*_flux.fits` + invvar using `master_frb_localization.csv`. Default: `coord_semantics=host` only. Flags: `--frb`, `--skip-existing`, `--include-signal`, `--dry-run`. Captures stdout to `pipeline_scripts/Output/<FRB>_all/master_run.log`.
- **`pipeline_scripts/SExtractor + PSFEx/run_psf_pipeline.py`** — Phase 1. Two-pass SExtractor + PSFEx (iterative `SAMPLE_MINSN` retry 30 → 10). Hard-coded catalog name `image.cat` so downstream phases find deterministic outputs.
- **`pipeline_scripts/photometry + astropath/run_photometry_astropath.py`** — Phase 2. PSF re-photometry (sole source of all downstream mags); three calibrated ZPs (40 px aperture, PSF, Auto) against PS1 (Vizier `II/349`) with **Legacy Survey DR10 Tractor (NOIRLab TAP) co-check** — whichever survey yields more matched stars (≥3); **`zero_points.json`** includes **`field_depth`** (5σ limiting mag at production aperture); default diagnostics **`sep_vs_shape_r.png`** and **`sep_vs_theta_max_phi.png`** via `pipeline_diagnostics.py`. AstroPath association via WSL Conda bridge. Star/galaxy cut uses `(SPREAD_MODEL + 3·SPREADERR_MODEL) < 0.005`; calibrated-magnitude sanity filter `[12, 28]` to reject corrupt-flux artefacts. Integration grid step set adaptively so absolute step ≤ `σ_loc / 5`. All AstroPath priors live in a single labelled block (`# ASTROPATH PRIOR CONFIGURATION`) reproducing the Aggarwal+2021 "adopted" defaults (`THETA_PDF='exp'`, `THETA_MAX=6.0`, `THETA_SCALE=1.0`). Display stretch in `astropath_association.png` is computed on the 1 arcmin zoom around the FRB.
- **`pipeline_scripts/photometry + astropath/field_depth.py`** — Sky σ, 5σ aperture `m_lim` (production 40 px); used by Phase 2 WSL bridge.
- **`pipeline_scripts/photometry + astropath/pipeline_diagnostics.py`** — `sep_vs_shape_r` / `sep_vs_theta_max_phi` plots; shared with M49/R70 standalone scripts.
- **`pipeline_scripts/galfit_fitting/generate_galfit_cutouts.py`** — Phase 3a. AstroPath mode uses `sex_number` from `astropath_posteriors.csv` (Phase-2 **`image.psf.cat`** `NUMBER`; SPREAD re-check), then maps to Phase-1 **`image.cat`** / segmentation ID by sky position (≤1″). Localization mode (`--use-localization-host`) picks the **nearest galaxy** within `--max-host-sep-arcsec` at the CSV RA/Dec using the same SPREAD cut (with coordinate-based SPREAD lookup when Phase-1/2 IDs differ). **Containment-based** neighbor handling: `frac ≥ --contain-thresh` (0.95) → fit; `--expand-thresh` ≤ frac < contain → grow ROI, then fit; `0 < frac < expand-thresh` → mask only. Loop bounded by `--max-roi-iterations`. The FRB host is unconditionally added to the fit set. `host_components.csv` rows are ordered with the host first so GALFIT component 1 = host. `host_sigma` is anchored to the empirical sky MAD when the invvar-derived sigma disagrees with sky scatter by more than 2× (preserves spatial structure, fixes Legacy Surveys frames with inconsistent flux / invvar units).
- **`pipeline_scripts/galfit_fitting/run_galfit_fitting.py`** — Phase 3b. Builds `galfit.feedme` and `constraints.txt` (`n ∈ [0.5, 6.0]`, `re ∈ [1.5, 100.0]`, **`mag ∈ [8, 40]`**). **`J)` ZP** from workdir `galfit_config.yaml` `mag_zeropoint`, else Phase-2 `zero_points.json` **`zp_aper_40px`**, else 22.5. **Sky:** SExtractor `BACKGROUND` seed; GALFIT crash detection; retry with ±`sky_tolerance_adu` only when pass 1 finishes and parsed sky drifts (not on crash); **`sky_fit_audit.json`** includes sky QA plus host **`snr_win`** / optional **`snr_auto`** from `host_components.csv` row 0. Defaults in `galfit_fitting/galfit_config.yaml`; `master_run.py` writes per-run config after Phase 2. `mag_init = MAG_40PX + mag_zeropoint`; `pa = THETA_IMAGE − 90°`.
- **`rerun_pipeline_galfit_phase3b.py`** — Re-run Phase 3b on `pipeline_scripts/Output/*_all/` (`--frb` optional).

Phase configuration files: `pipeline_config.yaml`, `photometry_astropath_config.yaml`, and **`galfit_fitting/galfit_config.yaml`** (`mag_min`/`mag_max`, sky QA, plate scale). Per-run copies staged into `<output>/.workdir/`; `galfit_config.yaml` also written to the workdir after Phase 2 with field-specific `mag_zeropoint`.

---

## Master tables & GALFIT `fit.log` parsing (`scripts/`)

- **`pipeline_asset_paths.py`** — Shared paths to `pipeline_scripts/Output/<FRB>_all/host_cutout.fits` and `host_sigma.fits` (replaces removed `cropped_host_galaxies/`).
- **`build_master_frb_localization.py`** — Builds `master_frb_localization.csv` from `master_frb_summary.csv` (coordinates, error ellipse, `z`, `DM`, `status`, `coord_semantics`, `survey` tag).
- **`merge_confident_hosts_localization.py`** — Merges host positions from `new_confident_hosts.txt` (LaTeX tables) into `master_frb_localization.csv` (`DM_MW`, `DM_exgal`, survey, z, DM, coords); skips CHIME.
- **`build_master_frb_galfit_from_logs.py`** — Builds `master_frb_galfit_from_logs.csv` from `tools/galfit/runs/<FRB>/{no_psf_sigma,with_psf_sigma}/fit.log` (full Sérsic parameters, inclinations, parse strategy, relative log paths).
- **`galfit_fitlog_parse.py`** — Shared parser (dash-line splitting, sane-χ² filtering, fixed-`n` refinement handling). Accepts `sersic_component_index` (default 0) for multi-component fits. **`parse_fitlog_sky_level(log_path)`** — fitted global sky (ADU). **`count_fitted_sersic_components(output_dir)`** — number of Sérsic components fitted (prefers `host_components.csv` row count; sky excluded). Imported by comparison, QA, and pipeline Phase 3b tooling.
- **`compile_galfit_logs.py`** — Writes **`galfit_metrics_summary.csv`** using paths under `tools/galfit/runs/` and `{no_psf_sigma, with_psf_sigma}`.
- **`audit_master_frb_galfit_csv.py`** — Verifies `master_frb_galfit_from_logs.csv` round-trips against fresh parses and checks expansion-16 rows vs `new_16_frbs_galfit_results.csv`.

### Pipeline benchmarking against the master table

- **`compare_pipeline_galfit_vs_master.py`** — Parses pipeline `Output/<FRB>_all/fit.log` (**GALFIT component 1** = `sersic_component_index=0`) and joins `master_frb_galfit_from_logs.csv`. Writes **`pipeline_galfit_results.csv`** (currently **62** production hosts; includes `mag`, **`n_sersic_components`**, **`host_number`**, **`snr_win`** / **`snr_auto`** from `sky_fit_audit.json`, informational **`single_sersic`**) and **`pipeline_vs_master_galfit_diff.csv`** (shape deltas only — **no mag/flux comparison**). Experimental runs belong in `pipeline_scripts/docs/EXCLUDED_RUNS.md`, not `Output/`.
- **`rerun_pipeline_galfit_phase3b.py`** — Re-run Phase 3b on existing `pipeline_scripts/Output/*_all/` (writes `galfit_config.yaml` from `zp_aper_40px`, then `run_galfit_fitting.py`). Optional `--frb`.
- **`analyze_pipeline_vs_master_diff.py`** — Detailed statistics on `pipeline_vs_master_galfit_diff.csv`: mean / median / std, MAD, quantiles, RMSE, robust statistics after excluding `chi2nu` blow-ups, and Spearman correlations across deltas.
- **`flag_pipeline_unphysical_fits.py`** — Heuristic QA on pipeline GALFIT outputs → `pipeline_unphysical_fits_review.csv`. Flags include bright `mag`, `re`/`n` at constraints, face-on `b/a`, high χ²/ν, `bad_sky`, **`sky_qa_failed`** (from `sky_fit_audit.json`), multi-Sérsic blocks. Many flags are **false positives** after per-field ZP + sky QA (see handoff). Workflow: `Reports/pipeline_galfit_review_handoff.md`.

---

## GALFIT: legacy workflow, expansion batch, inclinations (`scripts/`)

- **`setup_galfit_stage1.py`** / **`setup_galfit_stage2.py`** — Stage-1 no-PSF and stage-2 PSF `feedme` generation.
- **`run_single_galfit.py`** — Run GALFIT for one target directory.
- **`run_galfit_with_sigma.py`** — Sigma-weighted batch workflow (fixed `n=1`), aggregates parameters.
- **`run_galfit_16_expansion.py`** — Orchestrates the 16-FRB expansion set (uses **`get_initial_guesses.py`**); writes results rows.
- **`run_final_inclinations.py`** — Final inclination pass over fit logs.
- **`calculate_galfit_inclination.py`** — Parses logs and computes Hubble-style inclination from `b/a`.
- **`compile_true_freen_legacy_comparison.py`** — Free-`n` parameters from `galfit.0X` / logs vs Legacy.
- **`append_old_frbs_galfit_results.py`** — Merges older FRB free-`n` blocks into `new_16_frbs_galfit_results.csv`.

---

## Legacy Survey, SDSS, and cross-catalog comparison (`scripts/`)

- **`compare_galfit_vs_tractor_inclination.py`** — Legacy Tractor inclinations with uncertainty from shape ivars.
- **`query_legacy_tractor_shape_uncertainty_fields.py`** — Schema check for Tractor uncertainty columns.
- **`query_ls_photometry.py`** — TAP queries against `ls_dr10.tractor` for nearby sources.
- **`null_catalog_utils.py`** — Shared null cuts (`prepare_null_sample`, `prepare_null_strict_color_base`, `slice_null_base_by_mag`); morphology: SDSS `filter_sdss_drop_dev_winners` (lnL), Legacy `filter_legacy_spiral_morph` (EXP or \(n\in[0.75,2]\)); CDF constants `SDSS_UR_MAX_CDF=2.3`, `LEGACY_GR_MAX_CDF=0.75`, `LEGACY_CDF_TYPE_EXCLUDE=REX,DEV`; default SDSS shape `expAB_r`. RAM-safe catalog readers with `usecols`.
- **`build_legacy_catalog_csv.py`** — Builds **`LS_catalog_v1_allsky_modelmr.csv`** (Legacy DR10 TAP, `--region joint`, `tractor_mag_r`/`rmag`, shape \(e_1,e_2\); excludes `type=PSF` in query). Default `--top 500000`. v0: `Archive/csv/null_catalogs_v0/`.
- **`build_sdss_null_catalog.py`** — Builds **`SDSS_catalog_v1_allsky_modelmr.csv`** (SDSS DR16 `PhotoObj`, joint Legacy∩SDSS Dec footprint, RA-bin `TOP N` sampling). Default `--top 500000`. Prefer v2 for unbiased full-footprint work.
- **`build_sdss_null_catalog_v2.py`** — Builds **`SDSS_catalog_v2_fullsky_modelmr.csv`** (full SDSS footprint, HTM-hash random strata, `objID` dedupe). Iterates until strict null pool at **`modelMag_r < 20`** ≥ 50k (production cuts: u−r&lt;2.3, lnL exp-wins, `expAB_r`&gt;0.2). Cache: `SDSS_v2_fetch_cache/`. Audit: **`audit_sdss_v2_footprint.py`** → `plots/plots_null/v2_sdss_audit/`.
- **`patch_sdss_profile_winner.py`** — SQL-merge `lnLDeV_r`/`lnLExp_r`; sets `model_winner_is_exp`. Use `--merge-only` with `SDSS_lnl_patch_cache.csv` to re-merge without re-query. **`--footprint full`** for v2 (HTM strata; cache `SDSS_lnl_patch_cache_v2.csv`).
- **`audit_sdss_v2_footprint.py`** — v2 footprint / strict mag20 pool / objID dedupe audit vs optional v1.
- **`validate_sdss_profile_winner.py`** — lnL vs mag-proxy vs `fracDeV_r` audit → `diagnostics/sdss_profile_winner/`.
- **`diagnose_legacy_morphology.py`** — Legacy type / Sérsic-\(n\) pool counts; re-query gate → `diagnostics/legacy_morphology/`.
- **`archive_null_plots_pre_morphology.py`** — Moves `plots/plots_null/` to `Archive/plots_null_pre_morphology_cut/<timestamp>/` before morphology regen.
- **`augment_sdss_v1_colors.py`** — Adds `modelMag_u`, `modelMag_g`, `u_r`, `g_r` to an existing SDSS v1 CSV via sky match to PhotoObj.
- **`plot_sdss_color_cuts.py`** — **SDSS-only** u−r color-cut mag vs b/a (`diagnostics/color_cuts/`); streaming CSV (low RAM).
- **`plot_legacy_color_cuts.py`** — **Legacy-only** \(g-r < 0.75\) mag vs Tractor \(b/a\) (`diagnostics/color_cuts/gr_lt_0p75/`); streaming CSV (low RAM). Default matches null CDF cuts (also \(b/a>0.2\), no REX).
- **`audit_and_plot_null_v1_diagnostics.py`** — Funnel + mag vs b/a + REX only with **`--full`** (loads both v1 catalogs; high RAM). Default: prints safe commands and exits. See **`diagnostics/MEMORY_SAFE_NULL_WORK.md`**.
- **`plot_null_mag_cut_cdfs.py`** — Strict mag-cut CDFs under `mag_cuts/` with SDSS lnL exp-wins + \(u-r<2.3\) and Legacy EXP∪\(n\) + \(g-r<0.75\) (no REX/DEV). Default mag limits 24…15. Do **not** use `--no-clear` on a pre-morphology tree; archive first.
- **`run_sdss_frb_inclination_tests.py`** — Anderson–Darling + Mann–Whitney U for SDSS vs FRB at mag 20/21/22 (strict cos(i), strict i deg, inclusive cos(i)). Writes **`test_results.md`** at repo root.
- **`plot_sdss_quantile_inclination_bins.py`** — 8 equal-count SDSS cos(i) quantile bins (strict, lnL exp-wins, u-r color); `--mag-limit` 20/21/22 → `quantile8_mag{NN}_sdss_strict/`.
- **`enrich_sdss_best_model_ba.py`** — Optional re-merge of `best_model_ba_r` columns into an existing SDSS v1 CSV.
- **`enrich_sdss_null_size.py`** — Adds `n_eff_r` locally; `--query-radii` fetches `expRad_r`/`deVRad_r`/`best_model_re_r` via SDSS SQL.
- **`test_null_catalog_sanity.py`** — Validates null catalogs (`--survey-version v1|v2|both`; v2 checks objID unique + strict mag20 ≥ 50k).
- **`compare_sdss_legacy_null_distributions.py`** — SDSS vs Legacy null CDF overlay (both nulls on one figure). Output **`plots/plots_null/v1_null_plots/`**.
- **`pipeline_null_plot_utils.py`** — Shared CDF/histogram/random-host/sky helpers; temporary per-FRB plate scale (0.25″/pix for 8 northern LS hosts, else 0.262).
- **`plot_pipeline_diagnostics.py`** — Diagnostic suite (default **strict + color** null, tag `v1_allsky_modelmr_strict_color`): inclination FRB vs SDSS bin fractions, mag21 CDFs + combined overlay, Re/n CDFs, random-host test, master deltas. `--section all|hist|cdf|random|master`.
- **`plot_pipeline_null_overlays.py`** — Deprecated (duplicate of `null_cdf_pipeline_*`; use diagnostics).
- **`generate_legacy_cdf_null_plot.py`** — Legacy null CDF + FRB MC (shape-space); v1 CSV defaults; output `plots/plots_null/v1_null_plots/`.
- **`plot_legacy_cdf_inc_mc.py`** — Inclination-space MC vs Legacy null; v1 defaults; output `plots/plots_null/v1_null_plots/`.
- **`generate_galfit_mc_vs_sdss.py`** — GALFIT vs SDSS null (v1 SDSS CSV); output `plots/plots_null/v1_null_plots/`.
- **`compare_legacy_vs_galfit_reff.py`** — `shape_r` vs GALFIT `R_e` in arcsec.
- **`merge_reff_into_legacy_inclination_comparison.py`** — Merges Reff columns into inclination comparison CSVs.
- **`make_two_inclinations_csv.py`** — Compact collaborator table (`galfit_inc_psf_deg`, `ls_inc_deg`, errors).
- **`fetch_legacy_imr_cutouts_for_comparison.py`** — Downloads Legacy image/model/residual cutouts.

---

## Plotting & figure drivers (`scripts/`)

- **`generate_all_plots.py`** — Master plot driver.
- **`generate_sigma_plots.py`** — Sigma-weighted comparison figures.
- **`generate_multiband_cdf_null_plot.py`** — Multi-band null CDF overlays (CI policy). **Still defaults to archived v0 SDSS CSV** — pass `--sdss-csv` explicitly for v1.
- **`generate_top_deviation_imr_panel.py`** — Largest Legacy–GALFIT inclination disagreements, IMR panels.
- **`generate_bin_barchart.py`** — Inclination-bin bar charts.
- **`generate_psf_symmetry_plot.py`** — PSF symmetry from directional profiles.
- **`generate_report_pngs.py`** — FITS → report PNGs.
- **`generate_galfit_mc_vs_sdss.py`** — GALFIT vs SDSS Monte Carlo-style comparison plots.
- **`plot_r70_posterior_cutout.py`** — r70 field cutout with AstroPath overlays.
- **`plot_legacy_cdf_inc_mc.py`** — Combined MC CDF in inclination space vs Legacy null.
- **`generate_sigma_maps.py`** — Rewrites `feedme` to emit sigma maps via GALFIT (batch helper).

---

## PSF construction, PSFEx, and manifests (`scripts/`)

- **`build_catalog_psf.py`** / **`batch_build_catalog_psf.py`** — Empirical PSFs from catalog stars.
- **`downsample_psfs.py`** / **`normalize_psfs.py`** / **`clean_psfs.py`** — PSF resampling and housekeeping.
- **`analyze_psfs.py`** / **`evaluate_moffat_psf.py`** / **`diagnose_psf_stars.py`** — Quality and diagnostics.
- **`extract_center_local_psfs.py`** — Center-local PSF from spatially varying PSFEx model.
- **`psfex_local.py`** — Load/render PSFEx `.psf` models in Python.
- **`populate_psfex_sextractor_input_templates.py`** — Copies templates into per-run PSFEx folder layout under `Archive/psfs/`.
- **`bootstrap_psfex_sextractor_runs_and_download_5arcmin.py`** — Bootstraps PSFEx run trees and fetches 10-arcmin cutouts per manifest (filename legacy).
- **`update_psf_fwhm_summary_with_psfex.py`** — Updates `psf_fwhm_summary.csv` from PSFEx outputs.

---

## Data acquisition, `large_cutouts/`, and invvar sync (`scripts/`)

- **`cutout_download.py`** — Download **one** r-band **10′** cutout (Legacy → PS1 → DES). **Preflight** 32px probes skip tiers with no coverage (~1 min vs ~15 min wasted). `python cutout_download.py FRB [--force]`; `--preflight-only` to check coverage; registry: `large_cutouts/cutout_registry.csv`.
- **`consolidate_new_hosts_logs.py`** — Builds **`pipeline_scripts/new_hosts_master.csv`** + **`new_hosts_master.md`** (46-host cohort: cutout status, pipeline outcomes, batch log appendix). Run after batch runs: `python scripts/consolidate_new_hosts_logs.py`.
- **`consolidate_cutout_logs.py`** — Refreshes `large_cutouts/cutout_registry.csv` from disk + master cohort list.
- **`update_progress.py`** — Refreshes `large_cutouts/cutout_validation.csv`; then run consolidate for master MD/CSV.
- **`refresh_pipeline_status.py`** — Alias: re-runs `consolidate_new_hosts_logs.py` (live scan of `Output/`).
- **`download_missing_hosts.py`** — Batch-fetch cutouts for FRBs in `new_hosts_master.csv` still missing on disk.
- **`cutout_fetch_common.py`** / **`cutout_resample.py`** — shared fetch/resample helpers used by `cutout_download.py`.
- **`audit_cutout_coverage.py`** — Pre-flight probes (Legacy / PS1 / DES) → `large_cutouts/coverage_audit.csv`.
- **`cutout_fetch_common.py`** / **`cutout_resample.py`** — Shared fetch + 0.262″ grid standardization helpers.
- **`fetch_invvar_maps.py`** — Fetches invvar/weight maps into `large_cutouts/` for existing flux cutouts.
- **`fetch_and_sync_invvar_10arcmin.py`** — Syncs invvar with PSFEx manifest under `Archive/psfs/PSFEx + SExtractor`.
- **`sync_and_verify_true_10arcmin_inputs.py`** — Verifies flux/invvar vs manifest expectations.
- **`fetch_frb_cutouts.py`** / **`fetch_specific_frbs.py`** / **`fetch_excluded_frbs.py`** / **`fetch_excluded_panstarrs.py`** — Survey fetch helpers.
- **`fetch_sdss_data.py`** / **`fetch_panstarrs_data.py`** / **`fetch_legacy_survey_data.py`** — Survey-specific queries/cutouts.

Host science cutouts for fitting live under **`pipeline_scripts/Output/<FRB>_all/`** (`host_cutout.fits`, `host_sigma.fits` from Phase 3a). Wide-field inputs remain in **`large_cutouts/`**.

---

## Coordinates, catalogs, consolidation (`scripts/`)

- **`convert_catalog.py`** / **`merge_coords.py`** — Catalog merges and format conversion.
- **`update_coordinates_from_sheet.py`** / **`revert_coordinates.py`** — Coordinate updates from external sheets.
- **`consolidate_and_restructure_csvs.py`** / **`make_nonparam_analysis_readable.py`** — CSV reshaping.
- **`quick_check_master_vs_sheet_coords.py`** — Compares `master_frb_summary.csv` coordinates to the Excel estimate sheet.

---

## Model experiments & validation (`scripts/`)

- **`fit_all_sersic.py`** — Batch Sérsic fitting.
- **`fix_n1_fits.py`** / **`fix_all_n1_fits.py`** / **`compare_n1_fits.py`** — Fixed-`n` vs free-`n` experiments.
- **`check_fits_stats.py`** / **`inspect_discrepancies.py`** — Fit quality and discrepancy review.
- **`compare_sigma_results.py`** — Side-by-side comparison of two GALFIT metric CSVs (e.g. summary vs sigma summary).

---

## Reporting (`scripts/`)

- **`generate_tex_report.py`** — Main LaTeX report build.
- **`build_psfex_latex_report.py`** — PSFEx comparison LaTeX fragments.

---

## Manual inspection (DS9) (`scripts/`)

- **`open_ds9_batch.py`** / **`open_random_frbs_in_ds9.py`** / **`open_remaining_cutouts.py`** — DS9 launch helpers.

---

## Shared helpers (`scripts/`)

- **`get_initial_guesses.py`** — Image moments → structural guesses (imported by the expansion orchestrator and others).

---

## Photutils (`tools/photutils/scripts/`)

- **`run_photutils_ellipse_nopsf.py`** — No-PSF ellipse fitting, median isophotal averaging.
- **`plot_cdf_photutils_vs_galfit_sdss.py`** / **`analyze_photutils_vs_galfit.py`** — CDF and comparison tables. Photutils null overlay defaults: v1 SDSS CSV, `--sample-mode`, `--sdss-q-column` (via shared utils if updated).

---

## Simulation (`tools/simulation/`)

- **`generate_mock_galaxies.py`** / **`run_simulation_fits.py`** / **`evaluate_recovery_accuracy.py`** — Mock grids and recovery metrics.

---

## AstroPath (`tools/astropath/`, `tools/AstroPath/`)

- **`M49 astropath/run_m49_astropath.py`** — DESI M49 field association (4.6′ circular localization). Run via WSL `frb_project` (`run_via_wsl.ps1`). **`run_m49_grid_sweep.py`** — 18-run grid under `outputs/inverse/<pdf>/theta<NN>/`; **`run_m49_pO_prior_sweep.py`** — 54-run P(O) prior sweep. Plots: **`plot_m49_sweep_summary.py`** (objid 0 + P(U)), **`plot_m49_candidate_geometry.py`**, **`plot_path_offset_priors.py`**.
- **`R70 astropath/`** — R70 r-band coadd (FRB 04h17m35.9058s +07d55m51.9812s; **1″** localization, **60″** catalog). See **`R70 astropath/README.md`**.
  - **`run_r70_astropath.py`** — default association; **`plot_fov_overlay()`** writes 1′ cutout with all objid labels (`figures/r70_fov_overlay.png`). AstroPath lazy-imported so overlay regen works on Windows without `healpy`.
  - **`build_r70_photometry_from_psfcat.py`** — Phase-2-style calibration on Windows when WSL Vizier fails.
  - **`measure_r70_field_sensitivity.py`** — 5σ **aperture** depth (\(A=\pi(d/2)^2\) px², not seeing); per-aperture PS1 ZP sweep → `r70_depth_vs_aperture.csv`, `r70_depth_vs_aperture.png`; `r70_field_sensitivity.json`.
  - **`plot_r70_aperture_depth_variance.py`** — multi-panel aperture variance (noise-only vs point-source PSF model); `docs/R70_depth_literature.md`.
  - **`vet_r70_metadata.py`** — FITS vs config gain/pixel-scale report → `inputs/metadata_vet.json`.
  - **`r70_calibration.py`** — shared ZP/sky/aperture helpers.
  - **`run_r70_pO_prior_sweep.py`** — 54-run grid (`outputs/inverse{,1,2}/<pdf>/theta<NN>/`).
  - **`plot_r70_sweep_summary.py`** — sweep lines for **objid 33**, **objid 31**, and **P(U)** (not best-host P(O)).
  - **`plot_r70_candidate_geometry.py`**, **`plot_path_offset_priors.py`**, **`run_via_wsl.ps1`**.
- **`extract_candidates.py`** / **`run_path_analysis.py`** / **`update_notes.py`** / **`open_ds9.py`** — Expansion PATH workflow.
- **`astropath_pkg/scripts/run_expansion_association.py`** — Package-level expansion driver.
- **`plot_r70_legacy.py`** — Legacy standalone r70 visualization.

---

## High-precision photometry (`tools/Photometry/`)

- **`calibrate_photometry.py`** — Production dual-ZP photometry vs PS1, LS validation, 40 px aperture recovery.
- **`photometry.ipynb`** — SExtractor + PSFEx prototyping.

---

## Archived / retired scripts

- **`Archive/scripts_retired/integrate_12frbs_legacy.py`** — One-time merge of 12 expansion FRBs into `galfit_vs_legacy_*` and `master_frb_summary.csv` (kept for reproducibility only).

---

## Notes

- Pipeline GALFIT QA handoff (Tier A/B tables, ZP context, per-FRB checklist): **`Reports/pipeline_galfit_review_handoff.md`**. Machine-readable flags: **`pipeline_unphysical_fits_review.csv`**.
- Canonical Photutils comparison table at repo root: **`photutils_vs_galfit_comparison.csv`**.
- Older Photutils CSV/plot variants: **`Archive/csv/photutils_cleanup_20260326`**, **`Archive/plots/photutils_cleanup_20260326`**.
- Stored text artifact (not executable): **`Archive/csv/n1_comparison_results.txt`**.
- AstroPath diagnostic PNGs/CSVs moved to **`Archive/astropath_diagnostics/`**.

---

## Draft LaTeX

- **`Reports/02 PhotUtils + Legacy/galfitfreevsfix.tex`** — Fixed-`n` vs free-`n` structural comparison draft.
