# Project State

**Updated:** 2026-05-25 · **Tasks:** `tasks.md` · **Scripts:** `scripts.md`

Canonical history for the current repo. Obsolete experiments omitted.

---

## 1. Layout

| Area | Role |
|------|------|
| `pipeline_scripts/` | Production chain: `master_run.py` → Phases 1–3 in `Output/<FRB>_<tag>/`; batch via `run_all_frbs.py` |
| `large_cutouts/` | `{FRB}_flux.fits` + `{FRB}_invvar.fits` only (53 pairs) |
| `M49 astropath/` | Standalone M49/DESI AstroPath (4.6′ loc); 54-run prior sweep |
| `R70 astropath/` | Standalone R70 coadd AstroPath (1″ loc, 1′ catalog); 54-run prior sweep |
| `tools/` | Legacy GALFIT, AstroPhot, Statmorph, Photutils, AstroPath, simulation |
| `plots/` | CDFs, nulls, diagnostics |
| `Archive/` | Retired CSVs/scripts; v0 null catalogs |
| `Reports/` | Validation reports; `pipeline_galfit_review_handoff.md` |

### Pipeline (one line per phase)

1. **SExtractor + PSFEx** — 2-pass, `DEBLEND_MINCONT=0.005`, `image.cat` / `image.psf` / `proto_image.fits` / `segmentation_map.fits`.
2. **Photometry + AstroPath** — PSF cat = sole mags; ZP PS1 → SkyMapper fallback; spread-model star cut; `zero_points.json` (`zp_aper_40px` for GALFIT); `astropath_posteriors.csv`.
3. **GALFIT** — **3a:** containment ROI (host row 0 → component 1), σ anchored to sky MAD if invvar scale off; **3b:** `J)=zp_aper_40px`, PA=`THETA_IMAGE−90°`, sky QA, `mag 8–40`, `fit.log` + `galfit_results.png`. Config: `galfit_fitting/galfit_config.yaml`. **Known:** ROI can over-grow (**P7** in `tasks.md`).

`master_run.py`: stages workdir, CLI overrides (ellipse, seeing, `--galfit-zp`, etc.), copies deliverables. `use_weight_map` follows invvar presence.

---

## 2. Key CSVs

| File | Content |
|------|---------|
| `master_frb_localization.csv` | Host positions, `DM`/`DM_MW`/`DM_exgal`, `z`; `survey` = discovery array (ASKAP, DSA-110, …) from paper table via `merge_confident_hosts_localization.py` (CHIME excluded) |
| `master_frb_galfit_from_logs.csv` | Legacy `tools/galfit/runs` fits |
| `pipeline_galfit_results.csv` | **62** production pipeline hosts (`Output/<FRB>_all/`); **component 1**; `single_sersic` informational |
| `pipeline_vs_master_galfit_diff.csv` | Shape deltas only (no mag) |
| `pipeline_unphysical_fits_review.csv` | Heuristic QA (many Tier-A false positives post-ZP) |
| `LS_catalog_v1_allsky_modelmr.csv`, `SDSS_catalog_v1_allsky_modelmr.csv` | v1 nulls (model \(r\), joint Dec −30°…+90°) |

Parser: `scripts/galfit_fitlog_parse.py` (`sersic_component_index=0` = host).

**Comparison policy:** All hosts in `pipeline_galfit_results.csv`; no benchmark exclusions; multi-Sérsic hosts kept. Excluded trials: `pipeline_scripts/docs/EXCLUDED_RUNS.md` (e.g. **20220501C** marginal AstroPath).

---

## 3. Conventions

**Inclination** (all streams): \(\cos^2 i = (q^2 - q_0^2)/(1 - q_0^2)\), \(q_0=0.2\).

**Null catalogs (CDF plots):** **strict only** — \(b/a > q_0=0.2\). **Morphology:** SDSS drop deV winners (`lnLExp_r > lnLDeV_r`), shape **`expAB_r`**; Legacy **`EXP`** or Sérsic \(n\in[0.75,2]\), exclude **`REX`+`DEV`**. **Color:** SDSS \(u-r<2.3\); Legacy \(g-r<0.75\). FRB: mag + \(b/a>0.2\) only. Pre-cut plots: `Archive/plots_null_pre_morphology_cut/`. Methods: `plots/plots_null/.../NULL_CATALOG_AND_CDF_METHOD.md`.

**Plots — nulls:** `plot_null_mag_cut_cdfs.py` → `mag_cuts/` (PNG). SDSS u−r color cuts: `plot_sdss_color_cuts.py` → `diagnostics/color_cuts/` (streaming; SDSS only). Full audit funnel: `audit_and_plot_null_v1_diagnostics.py --full` only (high RAM).

**SDSS vs FRB tests:** `run_sdss_frb_inclination_tests.py` → **`test_results.md`** (AD + Mann–Whitney at mag 20/21/22; same strict/color cuts as CDF plots).

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

### R70 AstroPath (standalone)

- **Field:** `coadded_astrometrically_corrected_rband_r70.fits`; FRB at 64.3996075°, 7.931106°; localization **1″**; catalog **60″**.
- **Photometry:** Phase 1 in `R70 astropath/pipeline_work/`; calibration via `build_r70_photometry_from_psfcat.py` (Windows Vizier when WSL DNS fails).
- **Depth:** \(m_{\mathrm{lim,5\sigma}} \approx 25.1\) AB (`measure_r70_field_sensitivity.py`); nearest host **objid 33** at 1.82″, \(P(O|x)\approx 0.98\) (default inverse + exp, \(\theta_{\max}=6\)).
- **Sweep:** 54 runs under `outputs/inverse{,1,2}/`; plots in `figures/` — see `R70 astropath/README.md`.
- **Runbook:** WSL `.\R70 astropath\run_via_wsl.ps1`; plots native `python plot_r70_sweep_summary.py` etc.

---

## 5. Runbook snippets

```bash
# Single FRB
python pipeline_scripts/master_run.py --frb <FRB> --output pipeline_scripts/Output/<FRB>_all

# Batch (host coords)
python pipeline_scripts/run_all_frbs.py

# Phase 3b only on existing outputs
python scripts/rerun_pipeline_galfit_phase3b.py
```

WSL: `source-extractor`, `psfex`, `conda activate frb_project`, `galfit` via `wsl` from orchestrator scripts.

---

## 6. Guidance

- Paths here are authoritative; verify on disk if stale.
- Science host = **GALFIT component 1** unless a task says otherwise.
- New plots under `plots/`; no repo-root `tmp/`.
- Report style: `.agents/skills/report-writing/SKILL.md`.
