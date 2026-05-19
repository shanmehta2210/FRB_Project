# Tasks — Null catalogs, pipeline GALFIT, and CDF comparison

**Last updated:** 2026-05-19  
**Context:** (1) SDSS vs Legacy null CDF methodology (T1–T8 below). (2) Autonomous pipeline GALFIT hardening and validation (P1–P6). See `progress.md` for canonical state.

---

### Pipeline GALFIT (completed 2026-05-19)

| ID | Task | Status |
|----|------|--------|
| **P1** | Wire Phase 3b `mag_zeropoint` to Phase-2 `zero_points.json` **`zp_aper_40px`** (`master_run.write_galfit_config`, `run_galfit_fitting.load_mag_zeropoint`, template `galfit_config.yaml`) | `[x]` |
| **P2** | Sky QA on Phase 3b (SExtractor `BACKGROUND` seed, ±3 ADU retry, `sky_fit_audit.json`) — already in code; ensured active via config + batch | `[x]` |
| **P3** | `scripts/rerun_pipeline_galfit_phase3b.py` for Phase-3b-only refresh | `[x]` |
| **P4** | Full batch `run_all_frbs.py` end-to-end (39 host OK, 0 fail, 14 signal skipped); per-FRB logs in `Output/<FRB>_all/master_run.log` | `[x]` |
| **P5** | Manual review of Tier A “unphysical” fits; accept hosts after ZP fix (documented in `Reports/pipeline_galfit_review_handoff.md`) | `[x]` |
| **P6** | Refresh comparison CSVs: `n_sersic_components`, `compare_ok`; exclude **mag** from diff; single-Sérsic summary only | `[x]` |

**Deferred / follow-up (not blocking):**

- Re-run `flag_pipeline_unphysical_fits.py` with relaxed heuristics (negative sky, `mag<15` on bright hosts, χ² when σ OK).
- Re-run Phase 3a only on FRBs where `host_sigma` was not rescaled in an old folder (batch re-run should have fixed this).
- `generate_multiband_cdf_null_plot.py` / sigma plot drivers → v1 null CSV defaults.
- `--include-signal` batch for 14 skipped FRBs if radio-position runs are needed.

**Quick commands:**

```bash
python pipeline_scripts/run_all_frbs.py
python scripts/compare_pipeline_galfit_vs_master.py
python scripts/analyze_pipeline_vs_master_diff.py   # optional deeper stats
python scripts/flag_pipeline_unphysical_fits.py     # heuristic QA only
```

---

### v1 null catalogs (completed 2026-05-18)

| Output | Builder | Notes |
|--------|---------|-------|
| `LS_catalog_v1_allsky_modelmr.csv` | `scripts/build_legacy_catalog_csv.py` | Joint Dec −30°…+90°; `tractor_mag_r` / `rmag`; no misleading `petroMag_r` |
| `SDSS_catalog_v1_allsky_modelmr.csv` | `scripts/build_sdss_null_catalog.py` | DR16 `PhotoObj`; `rmag` = `cmodelMag_r`; RA/Dec exported |
| v0 CSVs | `Archive/csv/null_catalogs_v0/` | COSMOS SDSS + arbitrary Legacy brick — **do not use** for new plots |

Shared cuts / CDF modes: `scripts/null_catalog_utils.py` (`--sample-mode strict|inclusive`, `--mag-column rmag`). Sanity: `scripts/test_null_catalog_sanity.py`. Canonical figures: `plots/plots_legacy_cdf/v1_null_plots/` (filename tags e.g. `v1_allsky_modelmr_strict` / `_inclusive`).

**Deferred (still v0 defaults):** `generate_sigma_plots.py`, `generate_all_plots.py`, `generate_multiband_cdf_null_plot.py` — see Archive README.

---

## Background (read first)

### What null plots do

Scripts build a **background galaxy** cumulative distribution of \(\cos i\) from catalog axis ratios, then compare to FRB host inclinations (GALFIT MC). Shared formula (Hubble, \(q_0 = 0.2\)):

\[
\cos^2 i = \frac{q^2 - q_0^2}{1 - q_0^2}, \quad q \equiv b/a
\]

**Common selection cuts (plot time):**

| Cut | Field | Default |
|-----|-------|---------|
| Brightness | `rmag` / `tractor_mag_r <= 21` (v1; model \(r\)) | `mag_limit=21`, `--mag-column rmag` |
| Axis-ratio floor | `expAB_r > 0.2` (strict mode A) | `q0=0.2`, `--sample-mode strict` |
| Inclusive mode B | finite `expAB_r`; \(q \le q_0\) → \(\cos i = 0\) | `--sample-mode inclusive` |
| Legacy morphology (optional) | `tractor_type` not in exclude list | `REX` excluded in some scripts |

**Key scripts:**

| Script | Role |
|--------|------|
| `scripts/build_legacy_catalog_csv.py` | Builds v1 `LS_catalog_v1_allsky_modelmr.csv` from Legacy DR10 TAP |
| `scripts/build_sdss_null_catalog.py` | Builds v1 `SDSS_catalog_v1_allsky_modelmr.csv` from SDSS DR16 SQL |
| `scripts/null_catalog_utils.py` | Shared mag/q cuts and strict/inclusive CDF pools |
| `scripts/test_null_catalog_sanity.py` | Validates v1 catalogs (≥10k strict pool, footprint) |
| `scripts/convert_catalog.py` | **v0 only** — `SDSS_catalog.csv` from `SDSS_catalogue.txt` (archived) |
| `scripts/compare_sdss_legacy_null_distributions.py` | SDSS vs Legacy null CDF overlay |
| `scripts/generate_legacy_cdf_null_plot.py` | Legacy null + FRB MC (shape-space MC) |
| `scripts/plot_legacy_cdf_inc_mc.py` | Same; default `--exclude-types REX` |
| `scripts/generate_multiband_cdf_null_plot.py` | SDSS r/g/b nulls |
| `scripts/generate_galfit_mc_vs_sdss.py` | GALFIT vs SDSS null |
| `scripts/generate_sigma_plots.py`, `scripts/generate_all_plots.py` | Older SDSS null drivers |
| `tools/photutils/scripts/plot_cdf_photutils_vs_galfit_sdss.py` | Photutils/GALFIT vs SDSS |

### Issue A — Magnitude system mismatch (Caveat 1)

- **SDSS null** uses true **Petrosian** \(r\) (`petroMag_r` ← `rPmag` in `SDSS_catalogue.txt`).
- **Legacy null** uses **`petroMag_r` column that is misnamed**: it is `22.5 - 2.5*log10(flux_r)` from Tractor **`flux_r`** (model flux in \(r\), nanomaggies), per Legacy DR10 docs — **not** Petrosian.
- On-disk check (`SDSS_catalog.csv`, ~19.7k rows): median `petroMag_r - rmag` ≈ **+0.06 mag**; at `petroMag_r ≤ 21` there are **7,983** galaxies vs **8,568** for `rmag ≤ 21`; mean \(\cos i\) shifts by only ~0.01 between those cuts — **small for CDF shape, wrong for strict \(m<21\) interpretation**.
- **FRB hosts** use pipeline photometry on Legacy cutouts (SExtractor + PS1/SkyMapper ZP), not SDSS Petrosian — null should align with **Tractor/model** or host photometry, not Petrosian-only.

### Issue B — Footprint / spatial mismatch (Caveat 2)

- **SDSS catalog** is **`new_SDSS_DR16_cosmos`** — COSMOS only (Dec ≈ 1.2°–3.2° in exported CSV; no RA column in `SDSS_catalog.csv`).
- **Current `LS_catalog.csv`** (85k rows from `SELECT TOP 200000` with **no** `ORDER BY` or sky cut) occupies **RA ≈ 284°–288°, Dec ≈ 44°–46°** — a single Northern Legacy brick, **zero overlap** with COSMOS.
- FRB hosts are **all-sky** on Legacy imaging (`large_cutouts/`). Comparing COSMOS SDSS null to an arbitrary Legacy brick is **not** a same-field or same-survey comparison.
- Legacy DR10 **does** cover COSMOS; the bug is the TAP query design, not missing imaging.

### Archived reference (pre-fix diagnostics)

`Archive/csv/plots_legacy_cdf/null_distribution_stats_sdss_vs_legacy_no_rex.csv` — mean \(\cos i\) SDSS 0.494 vs Legacy 0.484 (similar despite footprint mismatch; do not treat as validation of current methodology).

---

## Task list

Status: `[ ]` todo · `[~]` in progress · `[x]` done

---

### T1 — Rebuild `LS_catalog.csv` with explicit spatial selection

**Status:** `[x]` Done (v1) — 2026-05-18  
**Priority:** High  
**Blocks:** T3, T4, T5

**Problem:** `scripts/build_legacy_catalog_csv.py` `run_query()` uses `SELECT TOP N` with quality cuts only — returns an arbitrary sky patch (currently RA ~285°, Dec ~45°), not COSMOS or FRB-relevant footprint.

**Implementation options (pick one primary; document choice in script help text):**

1. **COSMOS box** (SDSS-comparable null): e.g. `ra BETWEEN 149 AND 151 AND dec BETWEEN 1 AND 3` (tune to COSMOS DR16 footprint).
2. **FRB-host envelope** (science-fair null): query bricks/regions covering all `master_frb_summary.csv` RAs/Decs ± buffer (e.g. 0.5°–1°), or union of brick IDs from host coordinates.
3. **All-sky random sample:** `ORDER BY RANDOM()` or multiple brick queries until `TOP` filled — avoids single-patch bias; heavier TAP load.

**Requirements:**

- Add CLI flags, e.g. `--region cosmos|frb_envelope|all` and/or `--ra-min`, `--ra-max`, `--dec-min`, `--dec-max`.
- Add `ORDER BY` when using `TOP` (document that unordered `TOP` is non-reproducible).
- Log row count, sky bounds, and type histogram after build (existing summary prints are a good template).
- **Implemented:** `--region joint` (default) Dec −30°…+90°; TAP `ORDER BY RANDOM()` with post-fetch shuffle fallback; output `LS_catalog_v1_allsky_modelmr.csv` (500k rows; strict pool 28k+).

**Acceptance criteria:**

- Catalog coordinates overlap the intended region (verify min/max RA/Dec in a one-liner or unit test).
- `compare_sdss_legacy_null_distributions.py` compares samples from **comparable** sky when using COSMOS mode.

**Files:** `scripts/build_legacy_catalog_csv.py`, `LS_catalog_v1_allsky_modelmr.csv` (output)

---

### T2 — Fix magnitude column naming and semantics in Legacy catalog

**Status:** `[x]` Done (v1) — 2026-05-18  
**Priority:** High  
**Blocks:** T3, T4

**Problem:** `build_legacy_catalog_csv.py` writes Tractor model \(r\) mag into column `petroMag_r`, implying Petrosian parity with SDSS.

**Implementation:**

- Rename output column to `tractor_mag_r` or `model_mag_r` (keep `rmag` as alias if needed for backward compat).
- Optionally keep `petroMag_r` as deprecated duplicate with comment in CSV header doc only if breaking scripts is unacceptable — prefer a clean rename + script updates (T3).
- Document in script docstring: `flux_r` → AB mag via `22.5 - 2.5*log10(flux_r)` nanomaggy convention (Legacy DR10).
- Update `progress.md` §SDSS Null Construction to describe Legacy magnitude correctly.

**Acceptance criteria:**

- No Tractor magnitude labeled `petroMag_r` without explicit deprecation warning.
- `LS_catalog.csv` column names match physical meaning.

**Files:** `scripts/build_legacy_catalog_csv.py`, `progress.md`

---

### T3 — Unify magnitude cuts across null plotting scripts

**Status:** `[x]` Done (v1 plot drivers) — 2026-05-18  
**Priority:** High  
**Depends on:** T1, T2 (for Legacy); can partially do SDSS-only first

**Problem:** Null scripts cut on `petroMag_r <= mag_limit` for both surveys — apples-to-oranges for Legacy.

**Implementation:**

- **Legacy null scripts:** cut on `tractor_mag_r` / `model_mag_r` (post-T2 name) or `rmag` if identical to flux-derived \(r\).
- **SDSS null scripts:** add CLI flag `--mag-column {petro,model}` defaulting to:
  - `model` (`rmag`) when comparing to Legacy/FRB pipeline, or
  - `petro` when explicitly reproducing historical SDSS Petrosian null.
- Update cut logic in:
  - `scripts/compare_sdss_legacy_null_distributions.py`
  - `scripts/generate_legacy_cdf_null_plot.py`
  - `scripts/plot_legacy_cdf_inc_mc.py`
  - `scripts/generate_multiband_cdf_null_plot.py` (SDSS: already uses `rmag`/`gmag` for band limits; document that `expAB_r` is always r-band shape)
  - `scripts/generate_galfit_mc_vs_sdss.py`
  - `scripts/generate_sigma_plots.py`, `scripts/generate_all_plots.py` — **deferred** (v0 paths)
  - `tools/photutils/scripts/plot_cdf_photutils_vs_galfit_sdss.py` — **updated**
  - `tools/astrophot/scripts/plot_cdf_bias_astrophot_vs_galfit_psf.py` — **updated**
- Use same default `mag_limit=21` but document that SDSS `model` ≈ Legacy Tractor for fair comparison.

**Acceptance criteria:**

- SDSS vs Legacy overlay uses **equivalent flux definitions** at the brightness cut (recommend: both model \(r\)).
- Plot labels state which magnitude is used (e.g. “Tractor \(m_r\)” vs “Petrosian \(m_r\)”).
- Re-generated null CDFs archived under `plots/plots_legacy_cdf/` with new tag (e.g. `_modelmr_cosmos`).

**Files:** scripts listed above; plot output dirs

---

### T4 — Re-run SDSS vs Legacy null comparison with fixed catalogs

**Status:** `[x]` Done — 2026-05-18 (`plots/plots_legacy_cdf/v1_null_plots/`)  
**Priority:** Medium  
**Depends on:** T1, T2, T3

**Actions:**

1. Rebuild `LS_catalog.csv` (T1, COSMOS or FRB envelope per project decision).
2. Run `scripts/compare_sdss_legacy_null_distributions.py` with aligned mag columns and document `--tag`.
3. Update stats CSVs in `plots/plots_legacy_cdf/` or `Archive/csv/plots_legacy_cdf/`.
4. Note in commit/message whether KS/Wasserstein shifts vs archived `null_distribution_stats_sdss_vs_legacy_no_rex.csv`.

**Acceptance criteria:**

- `null_overlay_*.png/pdf` and `null_distribution_stats_*.csv` reflect fixed methodology.
- README or `progress.md` points to new canonical null comparison tag.

---

### T5 — Add optional FRB-footprint Legacy null (recommended science path)

**Priority:** Medium  
**Depends on:** T1, T2, T3

**Rationale:** FRB hosts are measured on Legacy cutouts worldwide; the fairest null may be Legacy galaxies in similar footprints, not COSMOS-only SDSS.

**Implementation:**

- Extend T1 with `--region frb_envelope` using `master_frb_summary.csv` (and/or `master_frb_localization.csv`) coordinates.
- New script or flag in `generate_legacy_cdf_null_plot.py` to use `LS_catalog_frb_envelope.csv` as default for FRB bias plots.
- Keep COSMOS-matched SDSS null as optional reference, not the only null.

**Acceptance criteria:**

- At least one plot pipeline uses FRB-footprint Legacy null by default for host comparison.
- Documented in `tasks.md` / `progress.md` which null is “primary” for the paper.

---

### T6 — SDSS catalog: export RA and document COSMOS provenance

**Status:** `[x]` Done (v1 SQL builder exports `RA_ICRS`, `DE_ICRS`) — 2026-05-18  
**Priority:** Low  
**Depends on:** none

**Problem:** `SDSS_catalog.csv` has `DE_ICRS` but no RA; COSMOS provenance only in `SDSS_catalogue.txt` header string.

**Implementation:**

- Update `scripts/convert_catalog.py` to preserve RA if present in source (inspect `SDSS_catalogue.txt` column list from first CSV row after quote-split).
- Add `scripts/README` note or comment in `convert_catalog.py`: source = `new_SDSS_DR16_cosmos`, SDSS DR16 COSMOS field.
- If RA unavailable in export, document that spatial matching to Legacy must use external table or re-query SDSS CAS.

**Acceptance criteria:**

- Future spatial cross-matches SDSS↔Legacy are possible without re-downloading full parent file.

**Files:** `scripts/convert_catalog.py`, `SDSS_catalog.csv`

---

### T7 — TAP query reproducibility and tests

**Status:** `[x]` Done — `scripts/test_null_catalog_sanity.py`; TAP `RANDOM()` falls back to shuffle  
**Priority:** Low  
**Depends on:** T1

**Implementation:**

- Add `ORDER BY objid` (or `random()`) to Legacy TAP query when using `TOP`.
- Optional: small `scripts/test_null_catalog_sanity.py` that asserts:
  - `expAB_r` ∈ [0, 1] for all rows
  - Legacy catalog RA/Dec within expected box for selected region
  - After cuts, null pool size ≥ FRB sample size (from `galfit_sigma_metrics_summary.csv` / `master_frb_summary.csv`)

**Acceptance criteria:**

- Re-running `build_legacy_catalog_csv.py` with same flags yields same row count and coordinate bounds.

---

### T8 — Documentation sweep

**Status:** `[x]` Done — 2026-05-18 (`tasks.md`, `progress.md`, `scripts.md`); pipeline section added 2026-05-19  
**Priority:** Low  
**Depends on:** T1–T5 (after implementation)

**Update:**

- `progress.md` — SDSS Null Construction, Legacy catalog, caveats resolved; May 19 pipeline batch + comparison policy.
- `scripts.md` — null-catalog builders, pipeline orchestrator, comparison/QA scripts.
- This file — null tasks T1–T8 `[x]`; pipeline tasks P1–P6 `[x]`.

---

## Decisions needed (human)

Record choice here before T1/T5 implementation:

| Decision | Options | Chosen |
|----------|---------|--------|
| Primary Legacy null footprint | COSMOS-only / FRB envelope / all-sky random | **Joint Legacy∩SDSS Dec (−30°…+90°), random sample** |
| Primary magnitude for \(m<21\) cut | SDSS+Legacy both model \(r\) / SDSS Petrosian + Legacy model (document mismatch) | **Both model \(r\)** (`rmag`, `tractor_mag_r`) |
| Break CSV column names? | Rename `petroMag_r` → `tractor_mag_r` (breaking) vs keep alias | **`tractor_mag_r` + `rmag`; no Legacy `petroMag_r`** |
| CDF sample modes | strict \(q>q_0\) only / + inclusive edge | **Both** (`--sample-mode strict\|inclusive`) |

---

## Quick verification commands (for implementers)

```bash
python scripts/test_null_catalog_sanity.py

python scripts/build_legacy_catalog_csv.py --region joint --top 500000
python scripts/build_sdss_null_catalog.py --top 500000

python scripts/compare_sdss_legacy_null_distributions.py \
  --galfit-csv galfit_metrics_summary.csv \
  --sample-mode strict --tag v1_allsky_modelmr_strict
```

---

## Related constants (do not change without project-wide review)

- \(q_0 = 0.2\) — `scripts/galfit_fitlog_parse.py` `Q0_DEFAULT`, all null scripts
- `mag_limit = 21.0` — default CLI
- Legacy type exclusion: `REX` in `plot_legacy_cdf_inc_mc.py`, `compare_sdss_legacy_null_distributions.py` (`--exclude-types`)
