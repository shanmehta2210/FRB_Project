# Tasks

**Updated:** 2026-05-25 · Canonical state: **`progress.md`** · Script index: **`scripts.md`**

---

## Active

| ID | Task | Priority |
|----|------|----------|
| **P7** | Phase 3a containment ROI: expansion can grow cutouts too much; extra field / extra Sérsics. Rework padding/expansion. **Host = GALFIT component 1** until then. | Low |
| **P8** | Per-FRB plate scale from cutout WCS → `galfit_config.yaml` / feedme + `plate_scale_arcsec` / `re_arcsec` in `pipeline_galfit_results.csv` (replace hardcoded 0.262 and diagnostic `PLATE_SCALE_BY_FRB` map). | Medium |
| — | `flag_pipeline_unphysical_fits.py` with relaxed heuristics | Low |
| — | Point `generate_multiband_cdf_null_plot.py` / sigma plot drivers at v1 null CSVs | Low |
| — | Signal-localized FRBs (`--include-signal`): exploratory only; see `pipeline_scripts/docs/EXCLUDED_RUNS.md` (20220501C marginal P(O), not in `Output/`) | Optional |

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

### Pipeline GALFIT (2026-05-19)
P1–P6: Phase 3b ZP from `zp_aper_40px`; sky QA + `mag 8–40`; comparison CSVs; handoff `Reports/pipeline_galfit_review_handoff.md`.

### New-host cohort + Phase 3a fix (2026-05-25)
46-host master table; fetch/coverage audit; AstroPath `sex_number` mapped Phase-2→Phase-1 by coordinates; production `Output/` held at **62** hosts (`pipeline_galfit_results.csv`); 20220501C signal trial documented in `EXCLUDED_RUNS.md`, not kept in `Output/`.

### Null catalogs v1 (2026-05-18)
T1–T8: LS + SDSS v1 CSVs; `null_catalog_utils.py`; sanity tests. SDSS enriched with `n_eff_r`, `best_model_re_r` (`enrich_sdss_null_size.py`).

### R70 AstroPath (2026-05-19)
Standalone `R70 astropath/`: pipeline photometry + PSFEx workdir; `build_r70_photometry_from_psfcat.py`; 70 candidates; default host **objid 33**; 54-run P(O) prior sweep; 5σ depth; figures include 1′ FOV overlay (all objids) and sweep lines for **objid 33**, **objid 31**, **P(U)**. Docs: `R70 astropath/README.md`, `scripts.md` AstroPath section.

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

Re arcsec in plots uses temporary per-FRB scale (8 hosts @ 0.25″/pix, else 0.262) until **P8**.

### Decisions (locked)
| Item | Choice |
|------|--------|
| Footprint | Joint Legacy∩SDSS Dec −30°…+90°, random sample |
| \(m<21\) cut | Model \(r\) both surveys |
| CDF modes | `strict` only; SDSS \(u-r<2.3\), Legacy \(g-r<0.75\) on nulls |
| Host in pipeline | GALFIT component 1 always |
