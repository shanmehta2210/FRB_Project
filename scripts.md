# Scripts Documentation

This file indexes the **active** analysis and utility code. The canonical end-to-end pipeline lives under `pipeline_scripts/`; this document also lists every Python module in the repository-root **`scripts/`** directory (as of May 18, 2026). Transient scratch files belong under `scripts/` or a user-local path — there is no repo-root `tmp/`.

---

## Autonomous Pipeline (`pipeline_scripts/`)

A single orchestrator (`master_run.py`) chains the three phases together over a per-run workdir. The phase scripts are independently runnable.

- **`pipeline_scripts/master_run.py`** — Orchestrator. Stages inputs, applies CLI YAML overrides, runs all phases, copies deliverables to `Output/<FRB>_<tag>/`. After Phase 2, writes **`galfit_config.yaml`** with `mag_zeropoint` = `zp_aper_40px` from `zero_points.json` (override: `--galfit-zp`). Tool keywords: `catalog`, `psf`, `photometry`, `astropath`, `galfit`, `all` (includes `sky_fit_audit.json`). Optional CLI: localization ellipse, `--detect-thresh`, `--pixel-scale`, `--seeing-fwhm`, `--gain`, `--mag-mode`, `--target-snr-min`, `--keep-workdir`. `use_weight_map` auto from `--invvar`.
- **`pipeline_scripts/run_all_frbs.py`** — Batch driver for `master_run.py` over `large_cutouts/*_flux.fits` + invvar using `master_frb_localization.csv`. Default: `coord_semantics=host` only. Flags: `--frb`, `--skip-existing`, `--include-signal`, `--dry-run`. Captures stdout to `pipeline_scripts/Output/<FRB>_all/master_run.log`.
- **`pipeline_scripts/run_all_frbs.py`** — Batch driver. Iterates every FRB in `master_frb_localization.csv` that has flux + invvar in `large_cutouts/`, passes per-FRB ellipses, captures each run's stdout to `<output>/master_run.log`, and prints a final per-FRB summary. By default runs only rows with `coord_semantics='host'`; pass `--include-signal` to add the burst-only rows.
- **`pipeline_scripts/SExtractor + PSFEx/run_psf_pipeline.py`** — Phase 1. Two-pass SExtractor + PSFEx (iterative `SAMPLE_MINSN` retry 30 → 10). Hard-coded catalog name `image.cat` so downstream phases find deterministic outputs.
- **`pipeline_scripts/photometry + astropath/run_photometry_astropath.py`** — Phase 2. PSF re-photometry (sole source of all downstream mags); three calibrated ZPs (40 px aperture, PSF, Auto) against PS1 (Vizier `II/349`) with **SkyMapper DR1.1 (`II/358`) fallback** when PS1 returns nothing; AstroPath association via WSL Conda bridge. Star/galaxy cut uses `(SPREAD_MODEL + 3·SPREADERR_MODEL) < 0.005`; calibrated-magnitude sanity filter `[12, 28]` to reject corrupt-flux artefacts. Integration grid step set adaptively so absolute step ≤ `σ_loc / 5`. All AstroPath priors live in a single labelled block (`# ASTROPATH PRIOR CONFIGURATION`) reproducing the Aggarwal+2021 "adopted" defaults (`THETA_PDF='exp'`, `THETA_MAX=6.0`, `THETA_SCALE=1.0`). Display stretch in `astropath_association.png` is computed on the 1 arcmin zoom around the FRB.
- **`pipeline_scripts/galfit_fitting/generate_galfit_cutouts.py`** — Phase 3a. AstroPath-aware target picker (centres on the candidate with highest `posterior_O` ≥ `--min-astropath-posterior`). **Containment-based** neighbor handling: `frac ≥ --contain-thresh` (0.95) → fit; `--expand-thresh` ≤ frac < contain → grow ROI, then fit; `0 < frac < expand-thresh` → mask only. Loop bounded by `--max-roi-iterations`. The FRB host is unconditionally added to the fit set. `host_components.csv` rows are ordered with the host first so GALFIT component 1 = host. `host_sigma` is anchored to the empirical sky MAD when the invvar-derived sigma disagrees with sky scatter by more than 2× (preserves spatial structure, fixes Legacy Surveys frames with inconsistent flux / invvar units).
- **`pipeline_scripts/galfit_fitting/run_galfit_fitting.py`** — Phase 3b. Builds `galfit.feedme` and `constraints.txt` (`n ∈ [0.5, 6.0]`, `re ∈ [1.5, 100.0]`). **`J)` ZP** from workdir `galfit_config.yaml` `mag_zeropoint`, else Phase-2 `zero_points.json` **`zp_aper_40px`**, else 22.5. **Sky:** SExtractor `BACKGROUND` seed; retry with ±`sky_tolerance_adu` constraint if `|Δ| > 3` ADU; `sky_fit_audit.json`. Defaults in `galfit_fitting/galfit_config.yaml`; `master_run.py` writes per-run config after Phase 2. `mag_init = MAG_40PX + mag_zeropoint`; `pa = THETA_IMAGE − 90°`.
- **`rerun_pipeline_galfit_phase3b.py`** — Re-run Phase 3b on `pipeline_scripts/Output/*_all/` (`--frb` optional).

Phase configuration files: `pipeline_config.yaml`, `photometry_astropath_config.yaml`, and **`galfit_fitting/galfit_config.yaml`** (sky QA defaults, plate scale). Per-run copies staged into `<output>/.workdir/`; `galfit_config.yaml` also written to the workdir after Phase 2 with field-specific `mag_zeropoint`.

---

## Master tables & GALFIT `fit.log` parsing (`scripts/`)

- **`build_master_frb_localization.py`** — Builds `master_frb_localization.csv` from `master_frb_summary.csv` (coordinates, error ellipse, `z`, `DM`, `status`, `coord_semantics`, crop box `xmin`…`ymax`, `survey` tag).
- **`build_master_frb_galfit_from_logs.py`** — Builds `master_frb_galfit_from_logs.csv` from `tools/galfit/runs/<FRB>/{no_psf_sigma,with_psf_sigma}/fit.log` (full Sérsic parameters, inclinations, parse strategy, relative log paths).
- **`galfit_fitlog_parse.py`** — Shared parser (dash-line splitting, sane-χ² filtering, fixed-`n` refinement handling). Accepts `sersic_component_index` (default 0) for multi-component fits. **`parse_fitlog_sky_level(log_path)`** — fitted global sky (ADU). **`count_fitted_sersic_components(output_dir)`** — number of Sérsic components fitted (prefers `host_components.csv` row count; sky excluded). Imported by comparison, QA, and pipeline Phase 3b tooling.
- **`compile_galfit_logs.py`** — Writes **`galfit_metrics_summary.csv`** using paths under `tools/galfit/runs/` and `{no_psf_sigma, with_psf_sigma}`.
- **`audit_master_frb_galfit_csv.py`** — Verifies `master_frb_galfit_from_logs.csv` round-trips against fresh parses and checks expansion-16 rows vs `new_16_frbs_galfit_results.csv`.

### Pipeline benchmarking against the master table

- **`compare_pipeline_galfit_vs_master.py`** — Parses pipeline `Output/<FRB>_all/fit.log` (component 0 = host) and joins `master_frb_galfit_from_logs.csv`. Writes **`pipeline_galfit_results.csv`** (includes `mag` for reference, **`n_sersic_components`**, **`compare_ok`**) and **`pipeline_vs_master_galfit_diff.csv`** (deltas for `chi2nu`, `re`, `n`, `b/a`, `pa`, `inc` only — **no mag/flux comparison**). Summary statistics use **`compare_ok=True`** (single-Sérsic) rows only. Excludes `20171020A`, `20220509G`, `20240210A` by default.
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
- **`null_catalog_utils.py`** — Shared null cuts (`prepare_null_sample`), Hubble \(\cos i\), joint-footprint constants; `--sample-mode strict|inclusive`; `--q-column` for SDSS `expAB_r` vs `best_model_ba_r`.
- **`build_legacy_catalog_csv.py`** — Builds **`LS_catalog_v1_allsky_modelmr.csv`** (Legacy DR10 TAP, `--region joint`, `tractor_mag_r`/`rmag`, shape \(e_1,e_2\); excludes `type=PSF` in query). Default `--top 500000`. v0: `Archive/csv/null_catalogs_v0/`.
- **`build_sdss_null_catalog.py`** — Builds **`SDSS_catalog_v1_allsky_modelmr.csv`** (SDSS DR16 `PhotoObj`, `rmag`=`cmodelMag_r`, `expAB_r`, `best_model_ba_r` from deV vs exp `modelMag_r` winner, RA/Dec, chunked SQL). Default `--top 500000`.
- **`enrich_sdss_best_model_ba.py`** — Optional re-merge of `best_model_ba_r` columns into an existing SDSS v1 CSV (normally handled in `build_sdss_null_catalog.py`).
- **`test_null_catalog_sanity.py`** — Validates v1 catalogs (≥10k strict pool per survey, footprint, column contract).
- **`compare_sdss_legacy_null_distributions.py`** — Primary SDSS vs Legacy null CDF overlay. Defaults: v1 CSVs, `galfit_metrics_summary.csv`, output **`plots/plots_legacy_cdf/v1_null_plots/`**. Flags: `--sample-mode strict|inclusive`, `--mag-column rmag`, `--sdss-q-column` (`expAB_r` or `best_model_ba_r`), `--exclude-types REX` (default), `--include-rex`, `--tag`.
- **`generate_legacy_cdf_null_plot.py`** — Legacy null CDF + FRB MC (shape-space); v1 CSV defaults; output `v1_null_plots/`.
- **`plot_legacy_cdf_inc_mc.py`** — Inclination-space MC vs Legacy null; v1 defaults; output `v1_null_plots/`.
- **`generate_galfit_mc_vs_sdss.py`** — GALFIT vs SDSS null (v1 SDSS CSV; needs `legacy_vs_galfit_two_inclinations.csv`); output `v1_null_plots/`.
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
- **`populate_psfex_sextractor_input_templates.py`** — Copies templates into per-run PSFEx folder layout under `psfs/`.
- **`bootstrap_psfex_sextractor_runs_and_download_5arcmin.py`** — Bootstraps PSFEx run trees and fetches 10-arcmin cutouts per manifest (filename legacy).
- **`update_psf_fwhm_summary_with_psfex.py`** — Updates `psf_fwhm_summary.csv` from PSFEx outputs.

---

## Data acquisition, `large_cutouts/`, and invvar sync (`scripts/`)

- **`fetch_large_cutouts_v2.py`** — Primary 10-arcmin flux + invvar downloader (DR10 → DR9 → PS1).
- **`fetch_invvar_maps.py`** — Fetches invvar/weight maps into `large_cutouts/` for existing flux cutouts.
- **`fetch_and_sync_invvar_10arcmin.py`** — Syncs invvar with PSFEx manifest under `psfs/PSFEx + SExtractor`.
- **`sync_and_verify_true_10arcmin_inputs.py`** — Verifies flux/invvar vs manifest expectations.
- **`fetch_frb_cutouts.py`** / **`fetch_specific_frbs.py`** / **`fetch_excluded_frbs.py`** / **`fetch_excluded_panstarrs.py`** — Survey fetch helpers.
- **`fetch_sdss_data.py`** / **`fetch_panstarrs_data.py`** / **`fetch_legacy_survey_data.py`** — Survey-specific queries/cutouts.
- **`crop_images.py`** / **`update_specific_crops.py`** / **`refine_crop_20191001A.py`** — Crop maintenance.

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
