# Project State

**Updated:** 2026-08-05 · **Tasks:** `tasks.md` · **Scripts:** `scripts.md`

Canonical history for the current repo. Obsolete experiments omitted.

---

## 1. Layout

| Area | Role |
|------|------|
| `pipeline_scripts/` | Production chain: `master_run.py` → selective Phases 1–3b in `Output/<FRB>_<tag>/`; batch via `run_all_frbs.py`. Docs: `pipeline_scripts/README.md` |
| `large_cutouts/` | `{FRB}_flux.fits` + `{FRB}_invvar.fits` only (53 pairs) |
| `PATH/` | Standalone AstroPath validation: `M49 astropath/` (DESI M49), `R70 astropath/` (R70 coadd) — see [`PATH/README.md`](PATH/README.md) |
| `tools/` | Legacy GALFIT, AstroPhot, Statmorph, Photutils, AstroPath, simulation |
| `plots/` | CDFs, nulls, diagnostics |
| `CHIME/` | Isolated CHIME repeater side project (localizations, cutouts, pipeline, inclination CDF vs main repo) — see `CHIME/README.md` |
| `Archive/` | Retired CSVs/scripts; v0 null catalogs |
| `Reports/` | Validation reports; handoff stub → `Archive/reports/pipeline_galfit_review_handoff.md` |

### Pipeline (one line per phase)

1. **SExtractor + PSFEx** — WSL; `DEBLEND_MINCONT=0.005`; verifies `image.cat` / `image.psf` / `proto_image.fits` / `segmentation_map.fits`.
2. **Photometry + AstroPath** — PSF cat = sole mags; ZP PS1 ↔ Legacy TAP co-check; spread-model star cut; `zero_points.json` (`zp_aper` for GALFIT); `astropath_posteriors.csv`. WSL conda bridge; deliverables verified after exit 0.
3. **GALFIT** — **3a:** containment ROI (host row 0 → component 1), σ anchored to sky MAD if invvar scale off; **statmorph** (optional); **3b:** feedme input pre-check, `wsl galfit` (stdin closed, 1 h timeout), `J)=zp_aper`, PA=`THETA_IMAGE−90°`, sky QA, `mag 8–40`. Config: `galfit_fitting/galfit_config.yaml`.

`master_run.py`: stages workdir, CLI overrides (ellipse, seeing, `--galfit-zp`, etc.), **`--outputs` drives phase execution** (not just file collection), copies deliverables, always writes `pipeline_summary.json`. `use_weight_map` follows invvar presence.

---

## 2. Key CSVs

| File | Content |
|------|---------|
| `master_frb_localization.csv` | Host positions, `DM`/`DM_MW`/`DM_exgal`, `z`; `survey` = discovery array (ASKAP, DSA-110, …) from paper table via `merge_confident_hosts_localization.py` (CHIME excluded) |
| `master_frb_galfit_from_logs.csv` | Legacy `tools/galfit/runs` fits |
| `pipeline_galfit_results.csv` | **67** production pipeline hosts (`Output/<FRB>_all/`; 62 original + 5 accepted 2026-07-21); **component 1**; `single_sersic` informational; `mag_final` = ZP-corrected magnitude (LS DR10 / PS1 reference or rescaled pipeline ZP when Phase 2 calibration failed) |
| `pipeline_vs_master_galfit_diff.csv` | Shape deltas only (no mag) |
| `pipeline_unphysical_fits_review.csv` | Heuristic QA (many Tier-A false positives post-ZP) |
| `catalog/LS_catalog_v1_allsky_modelmr.csv`, `catalog/SDSS_catalog_v1_allsky_modelmr.csv` | v1 nulls (model \(r\), joint Dec −30°…+90°) |

Parser: `scripts/galfit_fitlog_parse.py` (`sersic_component_index=0` = host).

**M49 AstroPath:** `outputs/m49_candidates.csv` (224 galaxies); appendix `outputs/top_20_m49_{gaussian,uniform}.csv` (top 20 non-M49 hosts, exponential offset, \(\theta_{\max}=6\)); sensitivity audit `marginal_hosts_ls10_lookup.csv` (uniform sweep, \(P(O|x)>1\%\)).

**CHIME (isolated):** `CHIME/repeater_localizations.csv` (12 repeaters), `CHIME/frb_localizations.csv` (31 total), `CHIME/chime_host_magnitudes.csv`, `CHIME/chime_hosts_inclination.csv` — see `CHIME/README.md`. Not merged into main tables.

**Comparison policy:** All hosts in `pipeline_galfit_results.csv`; no benchmark exclusions; multi-Sérsic hosts kept. Excluded trials: `pipeline_scripts/docs/EXCLUDED_RUNS.md` (e.g. **20220501C** marginal AstroPath). **Weak associations inside the 67** (no published host lit and pipeline \(P(O)\le 0.95\)): `pipeline_scripts/docs/WEAK_ASSOCIATIONS_PRODUCTION67.md` + `production67_lit_astropath_audit.csv`.

---

## 3. Conventions

**Inclination** (all streams): \(\cos^2 i = (q^2 - q_0^2)/(1 - q_0^2)\), \(q_0=0.2\).

**Null catalogs (CDF plots):** **strict only** — \(b/a > q_0=0.2\). **Morphology:** SDSS drop deV winners (`lnLExp_r > lnLDeV_r`), shape **`expAB_r`**; Legacy **`EXP`** or Sérsic \(n\in[0.75,2]\), exclude **`REX`+`DEV`**. **Color:** SDSS \(u-r<2.3\); Legacy \(g-r<0.75\). FRB: mag + \(b/a>0.2\) only. Pre-cut plots: `Archive/plots_null_pre_morphology_cut/`. Methods: [`NULL_CATALOG_DATA_AND_INCLINATION_AUDIT.md`](NULL_CATALOG_DATA_AND_INCLINATION_AUDIT.md) (v1 METHOD file is a stub → Appendix A).

**Plots — nulls:** `plot_null_mag_cut_cdfs.py` → `mag_cuts/` (PNG). SDSS u−r color cuts: `plot_sdss_color_cuts.py` → `diagnostics/color_cuts/` (streaming; SDSS only). Full audit funnel: `audit_and_plot_null_v1_diagnostics.py --full` only (high RAM).

**SDSS vs FRB tests:** `run_sdss_frb_inclination_tests.py` → **`Archive/reports/test_results.md`** (AD + Mann–Whitney at mag 20/21/22; same strict/color cuts as CDF plots).

**Pixel scale:** WCS on `large_cutouts` shows **0.262″/px** (31/39 hosts) and **0.25″/px** (8 northern LS). GALFIT feedme still hardcoded 0.262; diagnostic Re plots use a temporary per-FRB map until **P8** (`tasks.md`).

---

## 4. Milestones (compressed)

| When | Highlights |
|------|------------|
| Mar 2026 | 23-host workflow; sigma GALFIT; AstroPhot/Statmorph/Photutils; simulation RMSE; reports 00–02 |
| Apr 2026 | 55 FRB summary; `large_cutouts` + PSFEx; autonomous `pipeline_scripts/`; containment cutouts; 39-target free-n compile |
| May 2026 | Master localization/GALFIT CSVs; repo hygiene; `master_run` + batch; SkyMapper ZP fallback; pipeline vs master tooling |
| May 19 | ZP/sky/mag constraints; full batch 39 OK; `20240210A` fixed; comparison CSVs refreshed; null v1 plots; `plots_legacy_cdf` → `plots_null` |
| May 19 | **R70 AstroPath** complete: pipeline photometry, 5σ depth, 54-run sweep, FOV 1′ overlay (objid labels), sweep plots for objid 33 / 31 / P(U) |
| May 23 | **Phase 2 diagnostics in production pipeline:** `field_depth` in `zero_points.json`; default `sep_vs_shape_r.png` + `sep_vs_theta_max_phi.png` (`field_depth.py`, `pipeline_diagnostics.py`) |
| May 25 | **46-host cohort** tracking (`new_hosts_master.csv`); cutout ladder Legacy→PS1→DES; Phase 3a **sex_number→seg NUMBER** sky match; ZP doc (AB from PS1/Legacy, not raw nanomaggies); null plots (`legacy_sdss_strict_combined`, SDSS color-cut mag vs b/a); signal-FRB probe — **20220501C** trial excluded from `Output/` (P(O)≈0.67); **62** folders in production `Output/` |
| May 28 | **M49 AstroPath** marginal-host audit: **8** non-M49 objids peak \(>1\%\) in uniform sweep (`marginal_hosts_ls10_lookup.csv`; LS DR10 + SIMBAD via `lookup_ls10_names.py`); sweep plots refreshed |
| Jun 28 | **M49 paper figures + appendix CSVs:** `sweep_P_O_M49_and_P_U_uniform_eellipse.png`, `sep_vs_shape_r.png` (`PAPER_FONT=16`); `top_20_m49_{gaussian,uniform}.csv` via `export_m49_marginal_hosts.py` |
| Jul 11 | **Pipeline robustness + selective execution:** `--outputs` drives phase execution (not just collection); WSL/GALFIT hang hardening (stdin closed, feedme pre-check, 1 h timeout); Phase 1/2 deliverable verification; batch dry-run / tag / list-file fixes |
| Jul 21 | **Production-67 association audit + doc consolidation:** flagged 5 weak associations (`WEAK_ASSOCIATIONS_PRODUCTION67.md`); Verdi+2025 still in prep; promoted several signal→host where published lit exists; merged overlapping null / dust / REX / CHIME docs (stubs + Archive copies); see §7 |
| Aug 2 | **Phase 3a ROI:** replaced containment-fraction neighbor policy with **20 px pad + `re_sep_factor`×Re_neighbor** (default 3; YAML/`--re-sep-factor`). Stars always masked; host remains GALFIT component 1. Shared `sersic_init.effective_re_px` with Phase 3b. |
| Aug 4 | **Docs + production re-run:** README / `PIPELINE_AUDIT` Sci5 / `scripts.md` / CHIME ROI notes updated for Re-separation. Re-run restricted to the **65**-host confirmed-literature cohort (`scripts/build_confirmed_lit_cohort.py` → `production_confirmed_lit_hosts.csv`; drops `20210214G` + the four Verdi-in-prep hosts). New live validator `scripts/watch_production_batch.py` reports `PASS`/`BROKEN` per host during the batch. |
| Aug 4 | **Batch-tag trap found:** `run_all_frbs.py --outputs galfit` writes `Output/<FRB>_galfit/` with the reduced deliverable set and never touches the production `_all` trees, while still printing `OK`. First re-run attempt hit this and was discarded; production refreshes now documented as `--outputs all`. |
| Aug 5 | **Re-separation re-run complete:** 65-host batch → **64 OK, 1 FAIL** (50 min). `20230930A` = tolerated Phase 2 ZP failure (sparse field, `mag_final` falls back to PS1). `20240304B` **removed from production** — Phase 3a found no galaxy within 5″ and no LS/PS1 match within 3″; z = 2.148 non-detection whose 2026-05 fit (`SNR_WIN=0.83`, `re` at floor, `b/a=0.06`) was noise. |
| Aug 5 | **Production tree cleaned to 64.** The 5 weak-association hosts still carried pre-Re-separation fits (`neighbor_policy` `null` / absent) and were being re-parsed into `pipeline_galfit_results.csv`, mixing ROI policies. Moved to `Output/_unconfirmed/` (outside the production glob, nothing deleted). Now **folders = results = cohort = 64**, all `re_separation`; enforced by new `scripts/audit_production_outputs.py`. |
| Aug 5 | **Fit verification suite** (`pipeline_scripts/verification/`, spec + results in `FIT_VERIFICATION_CHECKS.md`). Nine read-only diagnostics over all 64 hosts. All four falsification tests pass on the 53-host science cut: no PSF leakage (q vs `e_psf` r = +0.04, p = 0.77; PA alignment KS p = 0.85), AstroPhot agrees to Δq = −0.006 ± 0.008, sky ±1σ moves q by 0.004 (median), no `m=2` residual signal. Trust tiers **A = 34, B = 13, C = 6**. Exposed: 25/53 hosts have `Re < FWHM` (why the Fourier estimator is usable on only 20), χ²ν tracks SNR (r = +0.41) so it cannot gate, a 6.3σ Δmag–Re trend, and 15/53 with Sérsic `n` pinned at a bound. |

### M49 AstroPath (standalone)

- **Field:** DESI catalog **224** candidates; FRB at 187.4449°, 8.0004°; localization **4.6′** circle; host **objid 0** (M49 / NGC 4472).
- **Models:** `outputs/gaussian/` (`eellipse`) and `outputs/uniform/` (`uniform_eellipse` box-cut); each has full **54-run** prior sweep (`inverse`, `inverse1`, `inverse2` × `uniform`/`core`/`exponential` × \(\theta_{\max}\in\{6,12,18,24,30,36\}\)).
- **Appendix hosts:** `outputs/top_20_m49_gaussian.csv` and `outputs/top_20_m49_uniform.csv` — top 20 non-M49 galaxies by `posterior_O` (inverse P(O), exponential offset, \(\theta_{\max}=6\)); regenerate with `export_m49_marginal_hosts.py`.
- **Sensitivity audit (uniform only):** `marginal_hosts_ls10_lookup.csv` — objids with \(P(O|x)>1\%\) in any sweep cell (LS DR10 `brickname-objid`, SIMBAD name). Top peaks: VCC 1199 (4.6%), NGC 4467 (4.2%), NGVS J122958.95+075800.8 (3.4%).
- **Paper figures:** `figures/sweep_P_O_M49_and_P_U_uniform_eellipse.png` (2×3 box-cut panel); `figures/sep_vs_shape_r.png`; typography `m49_figure_font.PAPER_FONT = 16`.
- **Other plots:** `figures/sweep_P_O_M49_vs_theta_max_{gaussian,uniform}.png` etc.; `python plot_m49_sweep_summary.py [--localization-pdf gaussian|uniform]`.
- **Runbook:** WSL `.\PATH\M49 astropath\run_via_wsl.ps1`; plots native Windows — see `PATH/M49 astropath/README.md`.

### R70 AstroPath (standalone)

- **Field:** `coadded_astrometrically_corrected_rband_r70.fits`; FRB at 64.3996075°, 7.931106°; localization **1″**; catalog **60″**.
- **Photometry:** Phase 1 in `PATH/R70 astropath/pipeline_work/`; calibration via `build_r70_photometry_from_psfcat.py` (Windows Vizier when WSL DNS fails).
- **Depth:** \(m_{\mathrm{lim,5\sigma}} \approx 25.1\) AB (`measure_r70_field_sensitivity.py`); nearest host **objid 33** at 1.82″, \(P(O|x)\approx 0.98\) (default inverse + exp, \(\theta_{\max}=6\)).
- **Sweep:** 54 runs under `outputs/inverse{,1,2}/`; plots in `figures/` — see `PATH/R70 astropath/README.md`.
- **Runbook:** WSL `.\PATH\R70 astropath\run_via_wsl.ps1`; plots native `python plot_r70_sweep_summary.py` etc.

---

## 5. Runbook snippets

```bash
# Full pipeline (default --outputs all)
python pipeline_scripts/master_run.py \
  --image large_cutouts/<FRB>_flux.fits \
  --invvar large_cutouts/<FRB>_invvar.fits \
  --ra <RA> --dec <DEC> --keep-workdir

# AstroPath only (Phases 1 + 2)
python pipeline_scripts/master_run.py \
  --image large_cutouts/<FRB>_flux.fits \
  --invvar large_cutouts/<FRB>_invvar.fits \
  --ra <RA> --dec <DEC> --outputs astropath

# Batch (host coords)
python pipeline_scripts/run_all_frbs.py

# Preview phase selection without running
python pipeline_scripts/master_run.py ... --outputs galfit --dry-run

# Phase 3b only on existing outputs
python scripts/rerun_pipeline_galfit_phase3b.py
```

WSL: `source-extractor`, `psfex`, `conda activate frb_project`, `galfit` via `wsl` from orchestrator scripts. Canonical usage: [`pipeline_scripts/README.md`](pipeline_scripts/README.md).

---

## 6. Guidance

- Paths here are authoritative; verify on disk if stale.
- Science host = **GALFIT component 1** unless a task says otherwise.
- New plots under `plots/`; no repo-root `tmp/`.
- Report style: `.agents/skills/report-writing/SKILL.md`.

---

## 7. Documentation changes (2026-07-21)

### Weak associations (maintained)

| Path | Role |
|------|------|
| `pipeline_scripts/docs/WEAK_ASSOCIATIONS_PRODUCTION67.md` | **Canonical list** of 5 production-67 FRBs without published host lit and without pipeline \(P(O)>0.95\): `20210214G`, `20230913`, `20240104A`, `20240203`, `20250518` |
| `production67_lit_astropath_audit.csv` | Machine-readable lit vs AstroPath audit for all 67 |

Cross-linked from `EXCLUDED_RUNS.md`, `pipeline_scripts/README.md`, and this file §2.

### Doc merges (lossless; old paths kept as stubs or Archive)

| Former / overlapping | Now |
|----------------------|-----|
| `plots/.../NULL_CATALOG_AND_CDF_METHOD.md` (full) | Content → `NULL_CATALOG_DATA_AND_INCLINATION_AUDIT.md` **Appendix A**; live path is a stub |
| `REX_INCLINATION_RESEARCH.md` + `ELLIPTICAL_DISK_PLAN.md` | Merged → `plots/plots_null/v2/ls_audit/REX_AND_ELLIPTICAL_DISK.md`; old names are stubs |
| `EXTINCTION_CORRECTION_PLAN.md` (full) | Archived → `Archive/notes/EXTINCTION_CORRECTION_PLAN.md`; live path stub → `DUST_AND_MEDIAN_BA.md` |
| `Reports/pipeline_galfit_review_handoff.md` (stale N≈39) | Archived → `Archive/reports/…`; live path stub → EXCLUDED_RUNS / WEAK_ASSOCIATIONS |
| `CHIME/chime_phase4_pipeline_qc.md` (full) | Summary appendix in `CHIME/repeater_reported_values_sources.md`; full → `Archive/notes/…`; live stub |
| `CHIME/catalog/README.md` | Slimmed; fixed 12/31 → **16/35**; tiers kept (used by `build_frb_properties.py`); points to `SOURCES_AUDIT.md` |

### Script / doc reference updates

- `scripts/elliptical_disk_model.py`, `SCALED_IS_DEGENERATE_RYDEN.md`, `DUST_AND_MEDIAN_BA.md`, `ls_audit/README.md`, `DES_calibration_and_dust_research.md`, `scripts.md`, `pipeline_scripts/README.md` — point at merged/canonical paths.

### Phase 3a Re-separation ROI (2026-08-02 / docs 2026-08-04)

| Item | Location |
|------|----------|
| Algorithm + CLI | `pipeline_scripts/galfit_fitting/generate_galfit_cutouts.py` (`resolve_neighbor_re_roi`) |
| Shared Re seed | `pipeline_scripts/galfit_fitting/sersic_init.py` (`effective_re_px`) |
| Defaults | `galfit_config.yaml` `cutouts.{host_pad:20, re_sep_factor:3.0, max_roi_iterations:8}` |
| Operator docs | `pipeline_scripts/README.md` §Phase 3a; `docs/PIPELINE_AUDIT.md` Sci5 |
| Tests | `scripts/tests/test_neighbor_re_roi.py` |
| Production refresh | `run_all_frbs.py --list-file production_confirmed_lit_hosts.csv --outputs all --use-localization-host` (64 hosts) |
| Live QC | `scripts/watch_production_batch.py --list-file production_confirmed_lit_hosts.csv --since-epoch <ts>` |
| Post-batch gate | `scripts/audit_production_outputs.py` (folders == results == cohort, 0 unconfirmed) |
