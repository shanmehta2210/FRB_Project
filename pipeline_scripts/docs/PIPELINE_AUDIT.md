# Pipeline audit — FRB host galaxy imaging pipeline

> **Last updated**: 2025-06-17 (QoL + Science + Robustness pass)

## Legend

| Tag | Meaning |
|-----|---------|
| **F** | Fixed in an earlier pass (prior to the first refactor) |
| **R** | Correctness risk — fixed in the first refactor |
| **S** | Sloppiness / refactoring opportunity — fixed in the first refactor |
| **P** | Portability / environment — deferred (suggested feature) |
| **Q** | Quality-of-life improvement — implemented in this pass |
| **Sci** | Science / analysis enhancement — implemented in this pass |
| **Rob** | Robustness / testing — implemented in this pass |

---

## 1  Fixes applied in earlier passes (F1–F6)

| ID | File | What was wrong | Resolution |
|----|------|---------------|------------|
| F1 | `field_depth.py` | `measure_sky_sigma` unpacked `sigma_clipped_stats` as `(sigma, _, _)` — used the **mean** instead of the **std** for `sigma_sky`. | Unpacking changed to `(_, _, sigma_sky)`; added `np.isfinite` guard. |
| F2 | `run_psf_pipeline.py` | `CalledProcessError` was caught but the script exited with code 0, creating a false-positive for `master_run`. | Added `sys.exit(1)` in the `except` block. |
| F3 | `run_psf_pipeline.py` | Log message said "> 25 optimal stars" when the code tested `>=`. | Changed message to ">= 25". |
| F4 | `photometry_astropath_config.yaml` | Comment listed `sep_vs_theta_max_phi.png` — file was renamed to `sep_vs_x_max_reff.png`. | Comment updated. |
| F5 | `README.md` | Same stale filename + wrong PSFEx SAMPLE_MINSN retry description. | Updated to the actual 30 → 20 → 10 → 5 → 3 → 2.5 → 2 ladder. |
| F6 | `README.md` | Phase 3a config example showed `use_localization_host: false` — YAML default is `true`. | Corrected. |

---

## 2  Correctness risks — fixed in this refactor (R1–R4)

### R1  `feedme NaN passthrough` (run_galfit_fitting.py)

**Problem:** `build_feedme_and_constraints` read `MAG_40PX` or `MAG_AUTO` and
added the zeropoint without checking for SExtractor's 99.0 sentinel or NaN.
GALFIT would receive an out-of-range or NaN magnitude seed and crash.

**Fix:** The magnitude seed now cascades through `MAG_40PX` → `MAG_AUTO` →
midpoint of `[mag_min, mag_max]`, with `math.isfinite` + abs < 90 guard at
each step. The result is clamped into `[mag_min, mag_max]` so it never
violates the GALFIT constraint block.

### R2  `aperture index assumed ≥ 15` (multiple files)

**Problem:** Aperture handling was split across three places with a hard-coded
index 14 (the 15th slot of a 15-element PHOT_APERTURES ladder). Changing the
ladder width or the production aperture required editing code in multiple
files.

**Fix — new aperture system:**

* `pipeline_shared.py` defines a canonical `DEFAULT_APERTURE_DIAMS_PX` list
  (the historical 15-element ladder) and `resolve_apertures(sex_cfg)`:
  - reads `phot_apertures_px` (list or legacy scalar) + `production_aperture_px`
    from the SExtractor YAML block;
  - defaults to the largest aperture when `production_aperture_px` is null;
  - returns `(aperture_list, production_index, production_diameter)`.
* Phase 1 (`run_psf_pipeline.py`) and Phase 2 (`run_photometry_astropath.py`)
  both call `resolve_apertures()` → the `.param` template uses `{NAPER}` →
  SExtractor's `PHOT_APERTURES` string is generated dynamically.
* Phase 2's WSL bridge script receives `__PROD_APER_INDEX__` /
  `__PROD_APER_DIAM__` via placeholder injection so the zero-point and
  candidate magnitude are always computed at the resolved production aperture.
* Phase 3a receives `--mag-aper-index` from `master_run` and uses it to
  extract `MAG_40PX` from the correct column.
* `field_depth.py` now imports `DEFAULT_APERTURE_DIAMS_PX` from
  `pipeline_shared` and uses `max(apertures)` as the fallback production
  aperture (no longer references a fixed index 14).
* Both YAML configs (`pipeline_config.yaml`, `photometry_astropath_config.yaml`)
  now expose `phot_apertures_px` (list) and `production_aperture_px` (null =
  largest).

### R3  `zp_med is a clipped mean, not a median`

**Problem:** `sigma_clipped_stats` returns `(mean, median, std)`. The
zero-point computation unpacked as `zp_med, _, zp_std = ...`, assigning the
**mean** to `zp_med`. The variable name and docstring promised a median.

**Fix:** Changed unpacking to `_, zp_med, zp_std = ...` so `zp_med` is the
actual sigma-clipped **median**, which is more robust to asymmetric
reference-star outliers. Same fix applied to `zp_p_med` and `zp_auto_med`.

### R4  `standalone defaults differ from YAML`

**Problem:** In-code fallback defaults in `run_photometry_astropath.py` used
`detect_thresh=10` and `target_snr_min=1.0`, while the YAML defaults were
`detect_thresh=3` and `target_snr_min=0.0`.

**Fix:** Changed the in-code fallbacks to `3` and `0.0` respectively, matching
the YAML.

---

## 3  Sloppiness / refactoring — fixed in this refactor (S1–S6)

### S1  `duplicated SExtractor templates`

**Problem:** `TEMPLATE_CONV` (3 × 3 Gaussian kernel) and `TEMPLATE_NNW`
(neural-net weights file, 29 lines of floats) were copy-pasted identically in
Phase 1 and Phase 2.

**Fix:** Both constants now live in `pipeline_shared.py` and are imported by
both phase scripts. The dead copy in `run_photometry_astropath.py` has been
removed. The Phase 1 on-disk `default.nnw` file (formerly tracked in git) is
no longer committed — the template is written at runtime.

### S2  `committed run artifacts`

**Problem:** All three phase directories contained tracked FITS images,
catalogs, feedme files, QA PNGs, XML files, and other outputs from a single
test run. These made `git status` noisy and inflated the repo size.

**Fix:**
* `git rm --cached` for every non-source file across
  `SExtractor + PSFEx/`, `photometry + astropath/`, and `galfit_fitting/`.
* Added comprehensive `.gitignore` rules so future run products are never
  tracked. Only `.py`, `.yaml`, and documentation files remain in those
  directories.

### S3  `Phase 1 cleanup permanently disabled`

**Problem:** The Phase 1 `finally` block was hard-coded to skip removal of
generated SExtractor/PSFEx template files, with a comment "Cleanup Disabled
for Inspection phase!".

**Fix:** Cleanup is now **on by default**. The `--keep-templates` CLI flag
(forwarded by `master_run`) retains the files when manual inspection is needed.

### S4  `print-based logging`

**Problem:** Every phase script used raw `print()` statements with ad-hoc
`[*]` / `[!]` prefixes. There was no timestamp, no severity level, no way to
filter or redirect logs programmatically.

**Fix:** `pipeline_shared.get_logger(name)` provides a `logging.Logger` under
a single `frb_pipeline` root with a unified format:

```
[HH:MM:SS] LEVEL frb_pipeline.phase1 | message
```

All four phase scripts (Phase 1, 2, 3a, 3b) and `master_run.py` now use the
shared logger. `run_all_frbs.py` intentionally keeps its `print()` calls
because its aligned summary table and inline `end=""` progress updates are a
user-facing report, not diagnostic logging.

### S5  `magic numbers`

**Problem:** Science-relevant constants (min PSF stars, convolution-box pad,
sigma-rescale gates, PSF-match radius) were buried in code without
documentation or configuration handles.

**Fix — promoted to YAML:**

| Constant | Config key | File | Default |
|----------|-----------|------|---------|
| Min accepted PSF stars | `psfex.min_accepted_stars` | `pipeline_config.yaml` | 25 |
| GALFIT conv-box pad (px) | `conv_box_pad` | `galfit_config.yaml` | 24 |
| Sigma rescale gate [min, max] | `cutouts.sigma_rescale_min/max` | `galfit_config.yaml` | 0.5, 2.0 |
| No-data sigma sentinel | `cutouts.no_data_sigma` | `galfit_config.yaml` | 1e30 |
| PSF-match radius (arcsec) | `cutouts.psf_match_arcsec` | `galfit_config.yaml` | 0.5 |

### S6  `Optional[Path] type nit in master_run.run_phase`

**Problem:** `cwd` parameter was typed as `str` and its `None` handling was
ambiguous.

**Fix:** Signature changed to `cwd: Path | None = None` (and `cmd` is
now flexibly typed). `str(cwd)` conversion happens inside the function body
only when `cwd` is not `None`.

---

## 4  New infrastructure: `pipeline_shared.py`

A new module `pipeline_scripts/pipeline_shared.py` centralises:

| Export | Purpose |
|--------|---------|
| `TEMPLATE_CONV` / `TEMPLATE_NNW` | SExtractor static template strings (S1) |
| `DEFAULT_APERTURE_DIAMS_PX` | Canonical 15-element aperture ladder |
| `resolve_apertures(sex_cfg)` | YAML → `(list, prod_index, prod_diam)` (R2) |
| `format_phot_apertures(list)` | Render for SExtractor `PHOT_APERTURES` |
| `render_param_template(tmpl, n)` | Substitute `{NAPER}` in `.param` templates |
| `get_logger(name)` | Shared `logging` factory with consistent format (S4) |

Stdlib-only, so it's importable from every Windows-side phase script without
extra dependencies.

---

## 5  Quality-of-life improvements (Q1–Q6) — this pass

### Q1  `--rerun-phase` flag (`master_run.py`)

Skip earlier phases and re-execute from a specific phase (e.g. `--rerun-phase 3b`
to re-fit GALFIT without re-running SExtractor + photometry). Requires an
existing workdir from a prior `--keep-workdir` run. Valid phases:
`1`, `2`, `3a`, `statmorph`, `3b`.

### Q2  Progress ETA (`run_all_frbs.py`)

Running average ETA is printed inline after each FRB completes, using an
exponential moving average of elapsed times.

### Q3  `--parallel N` (`run_all_frbs.py`)

Run up to N FRBs concurrently using `concurrent.futures.ProcessPoolExecutor`.
Default 1 (sequential). When parallel > 1, progress is printed as completion
lines rather than inline updates to avoid interleaving.

### Q4  Auto-refresh `pipeline_galfit_results.csv`

After the batch loop, `run_all_frbs.py` automatically invokes
`scripts/compare_pipeline_galfit_vs_master.py` (if it exists) to keep the
metrics table in sync. Suppress with `--no-auto-refresh`.

### Q5  HTML batch report (`run_all_frbs.py`)

An HTML report (`Output/batch_report.html`) is generated after each batch run,
containing the summary table and embedded base64 thumbnails of the three
diagnostic PNGs for each FRB.

### Q6  `--dry-run` for `master_run.py`

Prints the commands that would be executed for all phases, then exits.

---

## 6  Science / analysis enhancements (Sci1–Sci4) — this pass

### Sci1  Effective radius in arcseconds

`pipeline_summary.json` now includes `re_arcsec` and `re_arcsec_err` alongside
the pixel values. Conversion uses the plate scale from `galfit_config.yaml`
(or `--pixel-scale` CLI override). The scale itself is recorded as
`plate_scale_arcsec_px` in the summary.

### Sci2  Monte Carlo inclination errors

The inclination error is now computed via Monte Carlo (10,000 samples from
`N(b/a, b/a_err)`) in addition to the analytic symmetric propagation. The MC
result reports **asymmetric** 16th/84th percentile CIs, which correctly
capture the non-linear error distribution near `q → q₀` (edge-on). The
summary includes `inclination_mc.{median_deg, err_lo_deg, err_hi_deg, p16_deg,
p84_deg}`. References: Holmberg (1946), Hubble (1926), Padilla & Strauss (2008).

### Sci3  Statmorph integration (Phase 3a.5)

Non-parametric morphology (CAS, Gini-M₂₀) is computed on the host cutout
using `statmorph`, running after Phase 3a (cutouts) and before Phase 3b
(GALFIT). Results are written to `statmorph_results.json` and folded into
`pipeline_summary.json`. **Future note:** statmorph metrics can motivate
GALFIT initial conditions (e.g. high asymmetry → multi-component fit).

### Sci4  PSF quality metrics

PSFEx XML diagnostics (FWHM mean/min/max/std, ellipticity, chi², number of
accepted/rejected stars) are parsed from `psfex.xml` and included in
`pipeline_summary.json` under `psf_quality`.

---

## 7  Robustness / testing (Rob1–Rob4) — this pass

### Rob1  Retry logic for Vizier/TAP queries

External catalog queries (PS1 Vizier, Legacy TAP) in Phase 2 are now wrapped
in a retry function with exponential backoff (3 attempts, base delay 2s). If
a query fails after all retries, the pipeline degrades gracefully (continues
with 0 matches from that catalog). **Future note:** additional survey
backends (SDSS, 2MASS, Gaia) as fallbacks for fields where neither PS1 nor
Legacy has coverage.

### Rob2  Input validation

`master_run.py` validates RA ∈ [0, 360) and Dec ∈ [−90, 90] before execution.
A WCS containment check warns (non-fatal) if the target pixel coordinate
falls outside the image footprint.

### Rob3  Unit tests

`tests/test_pipeline_shared.py`, `tests/test_galfit_fitlog_parse.py`, and
`tests/test_field_depth.py` cover pure-function modules with pytest.

### Rob4  CI validation

`.github/workflows/validate.yml` runs `py_compile` on all pipeline scripts
and `pytest tests/ -v` on push/PR to main.

---

## 8  Remaining items

### P1  WSL path conversion portability

`_to_wsl_path()` in Phase 2 hard-codes `C:` → `/mnt/c/`. On non-Windows or
non-WSL setups, SExtractor / PSFEx / GALFIT invocations would fail. A
`wsl wslpath -u` call or a `GALFIT_CMD` / `SEX_CMD` config key would make the
pipeline OS-agnostic.

### P2  `subprocess.run` binary names

Phase 1 calls `wsl source-extractor`; Phase 2 calls `wsl source-extractor` +
`wsl bash -ic`; Phase 3b calls `wsl galfit`. The binary names and the `wsl`
prefix are hard-coded. A portable solution would expose per-tool command
templates in YAML (e.g. `commands.sextractor: "wsl source-extractor"`).

**Status:** Deferred — earlier attempts to make these configurable introduced
regressions. Both P1 and P2 are left as suggested features.

### Future: Additional survey backends

The retry logic (Rob1) handles transient network failures, but fields outside
PS1/Legacy coverage still get 0 calibration stars. Adding SDSS, 2MASS, or
Gaia as additional reference catalogs (with automatic fallback order) would
improve robustness for these edge cases. This is flagged for future inspection.

### Future: Statmorph-driven GALFIT configuration

Statmorph metrics (Sci3) are intentionally computed before GALFIT. A future
enhancement could use high asymmetry or merger flags to automatically select
multi-component GALFIT models, adjust initial Sérsic index, or flag fits for
manual review.

---

## Resolved (batch 3 — 2026-06-17)

### Env1 — `frb_project` conda env reproducibility   ✅

Added `pipeline_scripts/environment_frb_project.yml` that fully specifies the
WSL-side conda environment. Includes numpy, scipy, astropy, pandas, matplotlib,
astroquery, photutils, healpy, pyvo, etc. with minimum version pins. Install
instructions added to README §2.

### Env2 — Pixel scale from WCS (never hardcoded)   ✅

Removed `pixel_scale` from all three YAML configs. `master_run.py` now calls
`astropy.wcs.utils.proj_plane_pixel_scales` on the input FITS and injects the
result into the workdir YAMLs. Phase scripts compute it as a fallback when run
standalone. The `--pixel-scale` CLI flag has been removed.

### Env3 — `mag_40px` → `mag_aper` rename   ✅

All references to `mag_40px` / `zp_aper_40px` / `MAG_CALIB_APER_40PX` renamed
to `mag_aper` / `zp_aper` / `MAG_CALIB_APER`. Backward-compatible reads
(`zp_data.get("zp_aper") or zp_data.get("zp_aper_40px")`) added in all
consumers. `--mag-mode` choices updated to `mag_aper | mag_psf | mag_auto`.

### Env4 — YAML-configurable AstroPath priors   ✅

Prior knobs (`p_o_method`, `theta_pdf`, `theta_max`, `theta_scale`) moved from
hard-coded constants to the `astropath:` block in
`photometry_astropath_config.yaml`. New candidate prior options exposed:
`inverse_ang` (÷R_eff) and `inverse_ang2` (÷R_eff²) alongside `inverse` and
`identical`. New offset profile `flat` added to AstroPath alongside `exp`,
`uniform`, `core`. AstroPath `priors.py` data model updated to match.

### Env5 — Configurable search radius   ✅

The 1-arcmin candidate search box is now controlled by
`astropath.search_radius_arcsec` (default 60.0). Propagated to:
candidate selection, `astropath_association.png` zoom, P(U) normalisation.

---

## 9  Testing checklist

After any future change to the pipeline scripts, verify:

- [ ] `python -m py_compile pipeline_scripts/<file>.py` for every touched file
- [ ] `python -m pytest tests/ -v` (unit tests)
- [ ] `python -c "from pipeline_shared import resolve_apertures; ..."` smoke test
- [ ] Single-FRB end-to-end: `python pipeline_scripts/master_run.py --image ... --invvar ... --ra ... --dec ... --keep-workdir`
- [ ] Check that `pipeline_scripts/Output/<frb>_all/` contains all expected deliverables
- [ ] Batch: `python pipeline_scripts/run_all_frbs.py --frb <2-3 FRBs> --keep-workdir`
