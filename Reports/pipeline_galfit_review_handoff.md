# Pipeline GALFIT review handoff

**Created:** 2026-05-18 · **Updated:** 2026-05-19  
**Purpose:** Working notes for pipeline GALFIT QA and comparison against the legacy master catalog.

**Status (May 19):** Full batch re-run complete (39 host OK). Per-field `zp_aper_40px` + sky QA active. Manual review concluded most automated Tier-A flags were false positives; hosts accepted for inclination work when `galfit_results.png` / component 1 look reasonable. Multi-component stamps: use component 1 only; set `compare_ok=False` in diff CSV.

---

## Quick start

```bash
python pipeline_scripts/run_all_frbs.py                    # full batch
python scripts/compare_pipeline_galfit_vs_master.py        # no mag deltas; n_sersic + compare_ok
python scripts/analyze_pipeline_vs_master_diff.py          # optional
python scripts/flag_pipeline_unphysical_fits.py            # heuristic only
python scripts/rerun_pipeline_galfit_phase3b.py --frb NAME # Phase 3b only
```

| Artifact | Path |
|----------|------|
| All pipeline fit parameters | `pipeline_galfit_results.csv` |
| Pipeline vs master (shape only) | `pipeline_vs_master_galfit_diff.csv` |
| Heuristic flags | `pipeline_unphysical_fits_review.csv` |
| Per-FRB outputs | `pipeline_scripts/Output/<FRB>_all/` |
| Sky QA audit (per run) | `pipeline_scripts/Output/<FRB>_all/sky_fit_audit.json` |
| Per-FRB run log | `pipeline_scripts/Output/<FRB>_all/master_run.log` |

---

## Sky seed + QA retry (Phase 3b)

`run_galfit_fitting.py` now:

1. Seeds GALFIT sky from SExtractor **`BACKGROUND`** on host row 0 in `host_components.csv` (ADU).
2. After pass 1, compares `fit.log` sky to that reference.
3. If `|Δ| > 3` ADU (default), reruns with `constraints.txt` line `{sky_comp} 1 -3 3`.
4. Writes `sky_fit_audit.json` (`sky_ref_adu`, `sky_pass1_adu`, `sky_pass2_adu`, `passed`).

Re-run Phase 3b only:

```bash
python pipeline_scripts/galfit_fitting/run_galfit_fitting.py --dir pipeline_scripts/Output/<FRB>_all
```

Config overrides in workdir `galfit_config.yaml`: `sky_check_enabled`, `sky_tolerance_adu`, `sky_max_retries`.

---

## Scope

- **36 FRBs** in `pipeline_galfit_results.csv` (from `pipeline_scripts/Output/*_all/fit.log`).
- **3 additional** pipeline folders **not** in that CSV (benchmark exclusions): `20171020A`, `20220509G`, `20240210A`.
- **39** total pipeline output folders under `pipeline_scripts/Output/`.

---

## What “physical” looks like (J)=22.5)

At `J) 22.5000` in `galfit.feedme` (nanomaggy / Legacy AB convention):

| Parameter | Typical OK range | Red flag |
|-----------|------------------|----------|
| `mag` | **17 – 24** | **&lt; 15** (too bright) or **&gt; 25** |
| `re` | ~1.5 – 15 px | **≥ 99 px** (pinned at constraint max **100**) |
| `n` | 0.5 – 4 | **≥ 5.95** (at ceiling **6.0**) |
| `b/a` | 0.2 – 0.9 | **≤ 0.12** → inclination forced to **90°** |
| `chi2nu` | ~0.4 – 2 | **&gt; 5** (often bad model / wrong components) |
| `parse_strategy` | `last_sane_single_sersic` | `last_block_first_sersic` → verify correct log block |

**Not trustworthy:** converged `mag` ≈ **4 – 10** with **`re` → 100** — model is absorbing diffuse light / sky / neighbors, not a normal host at \(m_r \sim 20\).

---

## Tier A — review first (7 FRBs)

After **2026-05-18** Phase 3b refresh (`zp_aper_40px` + sky QA). Regenerate with:

```bash
python scripts/rerun_pipeline_galfit_phase3b.py
python scripts/compare_pipeline_galfit_vs_master.py
python scripts/flag_pipeline_unphysical_fits.py
```

| FRB | mag | re | χ²/ν | Flags | Output folder |
|-----|-----|-----|------|-------|-------------|
| **20211127I** | 14.95 | 55.8 | 0.49 | Still bright; sky QA failed (ref≈115 ADU); large Re | `pipeline_scripts/Output/20211127I_all` |
| **20220319D** | 15.68 | 43.1 | 1.54 | n=6, 2 Sérsic in last block | `pipeline_scripts/Output/20220319D_all` |
| **20220307B** | 19.86 | 5.6 | 0.97 | Sky QA failed (ref≈5 ADU) | `pipeline_scripts/Output/20220307B_all` |
| **20220912A** | 19.52 | 5.2 | 0.94 | Sky ≈ −16 in fit.log | `pipeline_scripts/Output/20220912A_all` |
| **20211212A** | 16.47 | 12.9 | **9.62** | High χ²/ν; `last_block_first_sersic` | `pipeline_scripts/Output/20211212A_all` |
| **20220725A** | 18.14 | 6.6 | **6.78** | High χ²/ν | `pipeline_scripts/Output/20220725A_all` |
| **20171020A** | ~17* | ~26* | huge* | **Not in results CSV**; multi-Sérsic | `pipeline_scripts/Output/20171020A_all` |

\*From last `fit.log` block only (excluded from `pipeline_galfit_results.csv`).

### Resolved by ZP + sky QA (no longer Tier A)

Previously degenerate bright-mag / Re→100 cases now in normal mag range:

| FRB | Was (mag, re) | Now (mag, re) |
|-----|----------------|---------------|
| **20210807D** | 5.73, 100 | **17.48**, 21.9 |
| **20220207C** | 8.26, 23.3 | **18.14**, 5.6 (Tier B only if flagged) |
| **20220825A** | 10.35, 5.5 | **19.52**, 5.2 (ok) |

`20240210A` remains excluded from results CSV (Tier B in review file; sky/multi-Sérsic).

---

## Tier B — suspect (14 FRBs in CSV + 2 excluded folders)

May be usable for inclination depending on science case; inspect before trusting.

### Face-on / b/a floor (i → 90°)

| FRB | mag | b/a | inc |
|-----|-----|-----|-----|
| 20200430A | 21.87 | 0.05 | 90.0 |
| 20200906A | 20.03 | 0.12 | 90.0 |
| 20221106A | 21.57 | 0.12 | 90.0 |
| 20240310A | 20.55 | 0.08 | 90.0 |

### n pinned at 6.0

`20190611B`, `20190711A`, `20220914A` (also `20220319D` in Tier A).

### Parser used `last_block_first_sersic` (not `last_sane_single_sersic`)

`20191001A`, `20210410D`, `20220105A`, `20220319D`, `20220914A`, `20221106A`, `20240310A` — confirm the parsed block is the **host** (GALFIT component 1).

### Other Tier B (mostly bright SExtractor seeds or large uncertainties)

`20190608B`, `20190714A`, `20210320C`, `20211203C`, `20220310F`, `20220920A`, `20221012A`, `20230526A`, `20230708A`, `20230902A`, `20231226A`, `20240201A`, `20240208A`, `20240304A`, `20240318A`, `20220509G` (not in results CSV; χ²/ν ≈ 9.6).

Full detail: see `pipeline_unphysical_fits_review.csv` columns `flags` and `notes`.

---

## Strict “OK for now” subset (19 FRBs)

Passed automated cuts: `17 ≤ mag ≤ 24`, `re < 50`, `chi2nu ≤ 3`, `last_sane_single_sersic`, `b/a > 0.15`.

```
20180924B  20190102C  20190608B  20190611B  20190711A
20190714A  20210320C  20211203C  20220310F  20220920A
20221012A  20230526A  20230708A  20230902A  20231226A
20240201A  20240208A  20240304A  20240318A
```

Still spot-check QA images; some may have large parameter uncertainties.

---

## Per-FRB checklist (use for every Tier A / B case)

Under `pipeline_scripts/Output/<FRB>_all/`:

| # | File | What to check |
|---|------|----------------|
| 1 | `qa_cutout_mask.png` | Correct host? Neighbors fitted vs masked? |
| 2 | `host_components.csv` | Row 0 = host; `MAG_40PX`, `FLUX_RADIUS`, `ELONGATION` |
| 3 | `galfit.feedme` | `J) 22.5`; initial `3) mag` = `MAG_40PX + 22.5`; component count |
| 4 | `fit.log` | Last sane block: mag, re, n, b/a, sky, χ²/ν |
| 5 | `galfit_results.png` | Model vs residual structure |
| 6 | `host_cutout.fits` / `host_sigma.fits` | Flux scale; sigma anchoring (`k` in master_run log) |

---

## Likely failure modes to test

1. **SExtractor seed:** `MAG_ZEROPOINT=0` → `mag_init = MAG_40PX + 22.5`. Very negative raw mags → wrong initial flux (see `run_galfit_fitting.py` comments on thin-core degeneracy).
2. **Multi-Sérsic confusion:** Host must be **row 0** in `host_components.csv` / **first** `sersic` in `fit.log` (parser `sersic_component_index=0`).
3. **Re constraint:** `re` max **100 px** in `constraints.txt` — fits hit ceiling when modeling extended glare.
4. **Sky:** Check `sky_fit_audit.json`; if `retried: true` or `passed: false`, inspect whether SExtractor `BACKGROUND` matches the cutout (global vs local background).
5. **Sigma / invvar:** Bad weights → odd sky and χ²; check `host_sigma.fits` vs empirical sky MAD.
6. **ROI / masking:** Fringe neighbors not fully masked; containment logic in `generate_galfit_cutouts.py`.
7. **Wrong log block parsed:** `last_block_first_sersic` vs `last_sane_single_sersic` in `galfit_fitlog_parse.py`.

---

## Pipeline vs master magnitude offset (context)

When comparing to `master_frb_galfit_from_logs.csv`:

- **Pipeline:** always `J) 22.5` in feedme.
- **Legacy master (23 hosts):** `J) 25.0` from `scripts/run_galfit_with_sigma.py` → reported mag **+2.5** vs pipeline for same flux.
- **Legacy expansion (16 hosts):** already `J) 22.5`.

`scripts/compare_pipeline_galfit_vs_master.py` **subtracts 2.5 mag from `mag_master`** for hosts whose master feedme has `J) 25.0` before computing `mag_delta` (documented in script docstring; no extra CSV column).

After correction, median `mag_delta` ≈ **0.06 mag** (not ~−2.4).

---

## Related scripts

| Script | Role |
|--------|------|
| `scripts/flag_pipeline_unphysical_fits.py` | Builds `pipeline_unphysical_fits_review.csv` |
| `scripts/compare_pipeline_galfit_vs_master.py` | Builds `pipeline_vs_master_galfit_diff.csv` |
| `scripts/analyze_pipeline_vs_master_diff.py` | Summary stats on deltas |
| `scripts/galfit_fitlog_parse.py` | Shared fit.log parser |
| `pipeline_scripts/galfit_fitting/run_galfit_fitting.py` | Feedme: `mag = MAG_40PX + mag_zeropoint` |
| `pipeline_scripts/galfit_fitting/generate_galfit_cutouts.py` | Cutout, mask, `host_components.csv` |

---

## Suggested workflow for tomorrow

1. Open this file + `pipeline_unphysical_fits_review.csv`.
2. Work through **Tier A** table (11 FRBs) using the per-FRB checklist.
3. For each fix, note: seed mag, mask, sigma, parser block, re-run `master_run.py` for that FRB.
4. Re-run `flag_pipeline_unphysical_fits.py` and `compare_pipeline_galfit_vs_master.py`.
5. For inclination science, prefer the **strict OK** list (19) until Tier A/B are resolved or excluded.

---

## Open tasks (see also `tasks.md`)

- Rebuild `LS_catalog.csv` with correct sky footprint for null CDFs.
- Harmonize Legacy vs SDSS magnitude column naming (`petroMag_r` vs Tractor `flux_r`).
