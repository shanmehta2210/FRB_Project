# Tasks

**Updated:** 2026-08-04 · Canonical state: **`progress.md`** · Script index: **`scripts.md`** · Pipeline: **`pipeline_scripts/README.md`**

---

## Active

| ID | Task | Priority |
|----|------|----------|
| — | `flag_pipeline_unphysical_fits.py` with relaxed heuristics | Low |
| — | Point `generate_multiband_cdf_null_plot.py` / sigma plot drivers at v1 null CSVs | Low |
| — | Signal-localized FRBs (`--include-signal`): exploratory only; see `pipeline_scripts/docs/EXCLUDED_RUNS.md` / `WEAK_ASSOCIATIONS_PRODUCTION67.md` | Optional |

**Plots (regenerate):**
```bash
# SDSS u-r color-cut mag vs b/a (low RAM; SDSS only)
python scripts/plot_sdss_color_cuts.py

# Strict null CDFs with u-r / g-r color cuts (mag_cuts/)
python scripts/plot_null_mag_cut_cdfs.py

python scripts/plot_pipeline_diagnostics.py --section all
# Default: strict + color (tag v1_allsky_modelmr_strict_color)
```

**Pipeline refresh:**
```bash
python pipeline_scripts/run_all_frbs.py
python scripts/compare_pipeline_galfit_vs_master.py
```

---

## Done (summary)

### P7 — Phase 3a Re-separation ROI (2026-08-02)
Replaced containment-fraction neighbor expand/mask with `host_pad=20` +
`re_sep_factor × Re_neighbor` (default 3). Host remains GALFIT component 1.
Docs: `pipeline_scripts/README.md`, `docs/PIPELINE_AUDIT.md` Sci5; production
refresh via `run_all_frbs.py --list-file pipeline_galfit_results.csv --include-signal --outputs galfit`.

### Pipeline robustness + selective `--outputs` (2026-07-11)
`--outputs` drives phase execution (dependency map); WSL/GALFIT hang hardening (feedme pre-check, stdin closed, 1 h timeout); Phase 1/2 deliverable verification; batch dry-run / tag / list-file fixes. Docs: `pipeline_scripts/README.md`, `docs/PIPELINE_AUDIT.md` batch 4.

### P8 — plate scale from WCS (2026-06-17)
Per-FRB plate scale from FITS WCS → workdir configs / `pipeline_summary.json` `plate_scale_arcsec_px` + `re_arcsec`. Hardcoded 0.262 / `--pixel-scale` CLI removed.

### Pipeline GALFIT (2026-05-19)
P1–P6: Phase 3b ZP from `zp_aper`; sky QA + `mag 8–40`; comparison CSVs; association flags `pipeline_scripts/docs/WEAK_ASSOCIATIONS_PRODUCTION67.md` (archived handoff: `Archive/reports/pipeline_galfit_review_handoff.md`).

### New-host cohort + Phase 3a fix (2026-05-25)
46-host master table; fetch/coverage audit; AstroPath `sex_number` mapped Phase-2→Phase-1 by coordinates; production `Output/` now **67** hosts (`pipeline_galfit_results.csv`; 62 + 5 GTC-trial fits accepted 2026-07-21); 20220501C signal trial documented in `EXCLUDED_RUNS.md`, not kept in `Output/`.

### Null catalogs v1 (2026-05-18)
T1–T8: LS + SDSS v1 CSVs; `null_catalog_utils.py`; sanity tests. SDSS enriched with `n_eff_r`, `best_model_re_r` (`enrich_sdss_null_size.py`).

### R70 AstroPath (2026-05-19)
Standalone `PATH/R70 astropath/`: pipeline photometry + PSFEx workdir; `build_r70_photometry_from_psfcat.py`; 70 candidates; default host **objid 33**; 54-run P(O) prior sweep; 5σ depth; figures include 1′ FOV overlay (all objids) and sweep lines for **objid 33**, **objid 31**, **P(U)**. Docs: `PATH/R70 astropath/README.md`, `scripts.md` AstroPath section.

### Pipeline diagnostic plots (2026-05-19)
Histograms (inclination 4×22.5° bins; Re/n quartile bins), inclination/Re/n CDFs, random-host inclination test, pipeline vs master Δ plots, sky maps. Driver: `plot_pipeline_diagnostics.py`.

**Plot dirs under `plots/plots_null/`:**
| Dir | Content |
|-----|---------|
| `v1_null_plots/` | SDSS vs LS comparison |
| `v1_null_cdf_inclination/` | Inclination CDFs (`null_cdf_pipeline_*` only) |
| `v1_null_cdf_re/`, `v1_null_cdf_n/` | Size/index CDFs (inclusive) |
| `v1_hist_inclination/` | FRB vs SDSS inclination bin fractions |
| `v1_random_host_inclination/` | Random-draw inclination CDFs |
| `v1_pipeline_vs_master/` | Δi, Δre, Δn Cleveland plots |

Re arcsec in pipeline summaries uses WCS plate scale (P8 done). Older null-plot drivers may still use a temporary per-FRB scale map — check the plot script if regenerating.

### Decisions (locked)
| Item | Choice |
|------|--------|
| Footprint | Joint Legacy∩SDSS Dec −30°…+90°, random sample |
| \(m<21\) cut | Model \(r\) both surveys |
| CDF modes | `strict` only; SDSS \(u-r<2.3\), Legacy \(g-r<0.75\) on nulls |
| Host in pipeline | GALFIT component 1 always |
