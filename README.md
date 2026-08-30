# FRB host-galaxy inclination pipeline

[![Validate](https://github.com/shanmehta2210/FRB_Project/actions/workflows/validate.yml/badge.svg)](https://github.com/shanmehta2210/FRB_Project/actions/workflows/validate.yml)

Measure **inclination angles** of localized Fast Radio Burst (FRB) host galaxies and compare them to field-galaxy “null” samples.

This repository is the analysis follow-up to [Bhardwaj, Lee & Ji (2024, *Nature*)](https://doi.org/10.1038/s41586-024-08065-w), which reported an observational bias against detecting FRBs in edge-on hosts. The original paper notebooks and Table 1 live under [`Archive/original paper/`](Archive/original%20paper/). Everything else here is the expanded pipeline: host association, GALFIT Sérsic fits, fit verification, and survey-catalog comparisons.

**Not a data dump.** Wide-field cutouts, all-sky null catalogs, and per-host `Output/` trees stay off GitHub (see [What is not in this clone](#what-is-not-in-this-clone)).

---

## Current sample

| Name | Definition | *N* |
|------|------------|----:|
| Production fits | Published host localization + parseable GALFIT in `pipeline_scripts/Output/<FRB>_all/` | **64** |
| Science cut (`in_53`) | Production host with $m_r \le 22$ and $b/a > 0.2$ | **53** |
| Paper sample | Science cut **and** human-confirmed in [`host_confirmation.csv`](pipeline_scripts/verification/host_confirmation.csv) | **50** |

Three science-cut hosts are rejected on morphology / wrong object / sky: `20190711A`, `20221101B`, `20230930A`. See [`SCIENCE_CUT_AND_COHORT.md`](pipeline_scripts/verification/SCIENCE_CUT_AND_COHORT.md).

**Inclination convention** (every stream): Hubble thin-disk inversion with intrinsic axis ratio $q_0 = 0.2$,

$$
\cos^2 i = \frac{q^2 - q_0^2}{1 - q_0^2},\qquad q = b/a.
$$

Use **`mag_final`**, not raw GALFIT `mag`, for any magnitude-dependent cut (`pipeline_galfit_results.csv`). Science host = **GALFIT component 1**.

A first comparison of the confirmed-50 $\cos(i)$ values to an HSC disk analogue ($N = 24{,}450$) is consistent with a typical 50-galaxy draw (matched-*N* KS $p \approx 0.57$). Plots and tests: [`plots/plots_null/v2/frb_vs_hsc_confirmed50/`](plots/plots_null/v2/frb_vs_hsc_confirmed50/).

---

## How the pipeline works

One FRB, one flux + inverse-variance FITS pair, five phases driven by [`pipeline_scripts/master_run.py`](pipeline_scripts/master_run.py):

| Phase | What it does |
|------:|--------------|
| 1 | SExtractor detection + PSFEx PSF (via WSL) |
| 2 | PS1 / Legacy Survey zero-point, then [AstroPath](https://github.com/FRBs/astropath) host association |
| 3a | Host cutout, neighbor mask / joint-fit ROI (Re-separation) |
| Statmorph | Optional CAS + Gini–M$_{20}$ |
| 3b | GALFIT Sérsic (host = component 1) + sky QA → inclination |

`--outputs` selects **which phases run**, not only which files are copied (`all` is the production default). Every run writes `pipeline_summary.json`.

```text
cutout FITS  →  catalog + PSF  →  photometry + AstroPath
                                    ↓
                         host stamp + mask
                                    ↓
                    GALFIT  →  q, i, mag_final
```

Full flags, YAML knobs, and troubleshooting: **[`pipeline_scripts/README.md`](pipeline_scripts/README.md)**. Do not use that file as the GitHub landing page — it is the operator manual.

---

## Repository map

| Path | Role |
|------|------|
| [`pipeline_scripts/`](pipeline_scripts/) | Production chain (`master_run.py`, `run_all_frbs.py`) |
| [`pipeline_scripts/verification/`](pipeline_scripts/verification/) | Fit-quality suite, host triage, re-fits |
| [`scripts/`](scripts/) | Catalogs, CDFs, GALFIT parsers, unit tests |
| [`plots/plots_null/`](plots/plots_null/) | Null-catalog and FRB-vs-survey figures |
| [`tools/`](tools/) | Legacy GALFIT runs, AstroPath package, PSF $b/a$ sims, `zdm` |
| [`GTC data/`](GTC%20data/) | GTC visibility and science-review lists (not the imaging pipeline) |
| [`Archive/`](Archive/) | Retired CSVs, original-paper notebooks, old reports |
| [`Reports/`](Reports/) | Early validation write-ups |

Internal working notes: [`progress.md`](progress.md) (state), [`scripts.md`](scripts.md) (script index), [`tasks.md`](tasks.md).

---

## Key tables

| File | What it is |
|------|------------|
| [`master_frb_localization.csv`](master_frb_localization.csv) | Host coordinates, DM, $z$, survey, literature cite (99 rows) |
| [`production_confirmed_lit_hosts.csv`](production_confirmed_lit_hosts.csv) | 64-host published-localization cohort |
| [`pipeline_galfit_results.csv`](pipeline_galfit_results.csv) | Production GALFIT parameters + `mag_final` |
| [`pipeline_scripts/verification/host_confirmation.csv`](pipeline_scripts/verification/host_confirmation.csv) | Human confirm / reject gate |
| [`repeater_localizations.csv`](repeater_localizations.csv) | 20 localized repeating hosts; citations in [`repeater_localizations_README.md`](repeater_localizations_README.md) |

Rebuild the results table after a batch:

```powershell
python scripts/compare_pipeline_galfit_vs_master.py
python scripts/audit_production_outputs.py
```

---

## Quick start

**Runtime.** Python ≥ 3.10 on Windows; SExtractor, PSFEx, and GALFIT via **WSL**. Phase 2 needs the `frb_project` conda env (AstroPath + `pyvo`).

```bash
conda env create -f pipeline_scripts/environment_frb_project.yml
conda activate frb_project
cd tools/AstroPath/astropath_pkg && pip install -e .
```

Sanity check from PowerShell: `wsl source-extractor -v`, `wsl psfex -v`, `wsl galfit`.

**One host** (RA/Dec from `master_frb_localization.csv`; cutouts are local-only):

```powershell
python pipeline_scripts/master_run.py `
    --image  large_cutouts/20190608B_flux.fits `
    --invvar large_cutouts/20190608B_invvar.fits `
    --ra  334.02040578312 --dec -7.89886810526581 `
    --outputs all
```

**Batch** the 64-host cohort (`--outputs all` is required for production trees):

```powershell
python pipeline_scripts/run_all_frbs.py `
    --list-file production_confirmed_lit_hosts.csv `
    --outputs all `
    --use-localization-host
```

**Tests** (no WSL / no FITS):

```powershell
python -m pytest scripts/tests/ -v
```

CI (`.github/workflows/validate.yml`) compiles the phase scripts and runs that suite on Python 3.11.

---

## Fit verification

Automated checks do **not** decide the paper sample. The suite measures; [`host_confirmation.csv`](pipeline_scripts/verification/host_confirmation.csv) is the human gate.

Start at [`pipeline_scripts/verification/VERIFICATION_README.md`](pipeline_scripts/verification/VERIFICATION_README.md). Physics of the nine diagnostics: [`FIT_VERIFICATION_CHECKS.md`](pipeline_scripts/verification/FIT_VERIFICATION_CHECKS.md).

```powershell
cd pipeline_scripts/verification
python run_verification.py --checks all --jobs 4
```

---

## Null catalogs (field-galaxy comparison)

FRB $\cos(i)$ CDFs are compared to disk-like galaxies in SDSS, Legacy Survey, DES, and HSC after **strict** cuts: $b/a > 0.2$, plus morphology and colour (SDSS $u-r < 2.3$ and $\ln L_{\mathrm{exp}} > \ln L_{\mathrm{deV}}$; Legacy EXP or $n \in [0.75, 2]$, $g-r < 0.75$, no REX/DEV).

| Topic | Where |
|-------|--------|
| How catalogs and $\cos(i)$ are built | [`NULL_CATALOG_DATA_AND_INCLINATION_AUDIT.md`](NULL_CATALOG_DATA_AND_INCLINATION_AUDIT.md) |
| Shared cuts | [`scripts/null_catalog_utils.py`](scripts/null_catalog_utils.py) |
| Mag-cut CDFs | `python scripts/plot_null_mag_cut_cdfs.py` |
| Confirmed-50 vs HSC | `python scripts/plot_frb_vs_hsc_confirmed50.py` |

The all-sky CSVs themselves are **not** in git (`catalog/` is gitignored). Rebuild with `scripts/build_sdss_null_catalog_v2.py`, `build_legacy_catalog_v2_exp.py`, `build_hsc_kawinwanichakij_sample.py`, and siblings listed in [`scripts.md`](scripts.md).

---

## What is not in this clone

| Local / sister path | Why it is absent |
|---------------------|------------------|
| `large_cutouts/` | ~10′ flux + invvar FITS pairs |
| `pipeline_scripts/Output/` | Per-host workdirs (FITS, `fit.log`, PNG QA) |
| `catalog/`, `catalog_downloads/` | All-sky SDSS / LS / DES / HSC tables |
| `CHIME/`, `PATH/` | Split to [AstroPATH-additions](https://github.com/shanmehta2210/AstroPATH-additions) (M49/R70 + CHIME repeaters) |

Excluded or quarantined pipeline trials: [`pipeline_scripts/docs/EXCLUDED_RUNS.md`](pipeline_scripts/docs/EXCLUDED_RUNS.md). Weak associations: [`WEAK_ASSOCIATIONS_PRODUCTION67.md`](pipeline_scripts/docs/WEAK_ASSOCIATIONS_PRODUCTION67.md).

---

## Citation

**Published result**

> Bhardwaj, M., Lee, J. & Ji, K. 2024, Selection bias obfuscates the discovery of fast radio burst sources, *Nature*, **634**, 1065. [doi:10.1038/s41586-024-08065-w](https://doi.org/10.1038/s41586-024-08065-w) · [arXiv:2408.01876](https://arxiv.org/abs/2408.01876)

**This pipeline.** Inclinations and $b/a$ in `pipeline_galfit_results.csv` are from the GALFIT reductions here, not from the 2024 paper. Host coordinates and redshifts should be cited from the papers listed in `master_frb_localization.csv` / [`repeater_localizations_README.md`](repeater_localizations_README.md).

AstroPath priors follow [Aggarwal et al. (2021)](https://doi.org/10.3847/1538-4357/abf3c0).
