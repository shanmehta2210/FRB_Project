# Pipeline GALFIT review handoff

**Created:** 2026-05-18 · **Updated:** 2026-05-19  
**Purpose:** Working notes for pipeline GALFIT QA and comparison against the legacy master catalog.

**Status (May 19):** Full batch complete — **39 host FRBs** with `fit.log` in `pipeline_scripts/Output/*_all/`. Per-field `zp_aper_40px`, sky QA, and **`mag ∈ [8, 40]`** constraints active. **`20240210A`** re-run OK (7 Sérsic; host component 1 stable). For inclination science: always use **GALFIT component 1** (`host_components.csv` row 0; parser `sersic_component_index=0`). Multi-Sérsic stamps are acceptable; `single_sersic` in CSVs is informational only.

---

## Quick start

```bash
python pipeline_scripts/run_all_frbs.py                    # full batch
python scripts/compare_pipeline_galfit_vs_master.py        # host component 1; all 39 FRBs
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

`run_galfit_fitting.py`:

1. Seeds GALFIT sky from SExtractor **`BACKGROUND`** on host row 0 in `host_components.csv` (ADU).
2. Runs pass 1 (sky free). Detects **GALFIT crashes** in WSL output; does **not** treat a crash as sky drift.
3. If pass 1 finishes and `|sky_fit − sky_ref| > 3` ADU, reruns with `constraints.txt` line `{sky_comp} 1 -3 3`.
4. Writes `sky_fit_audit.json` (`sky_ref_adu`, `sky_pass1_adu`, `sky_pass2_adu`, `galfit_pass1_ok`, `failure_reason`, `passed`).

Re-run Phase 3b only:

```bash
python pipeline_scripts/galfit_fitting/run_galfit_fitting.py --dir pipeline_scripts/Output/<FRB>_all
```

Config in workdir `galfit_config.yaml`: `mag_min`/`mag_max`, `sky_check_enabled`, `sky_tolerance_adu`, `sky_max_retries`.

---

## Scope

- **39 FRBs** in `pipeline_galfit_results.csv` and `pipeline_vs_master_galfit_diff.csv` (all host pipeline outputs under `pipeline_scripts/Output/*_all/`).
- **No benchmark exclusions** — includes `20171020A`, `20220509G`, `20240210A`.
- **23** single-Sérsic, **16** multi-Sérsic (host metrics always from component 1).

---

## Host identification policy

| Rule | Detail |
|------|--------|
| Catalog order | Row 0 in `host_components.csv` = FRB host |
| GALFIT | Component 1 = first `sersic` in `fit.log` |
| Parser | `parse_fitlog_file(..., sersic_component_index=0)` |
| Comparison CSVs | Structural columns = host only; `n_sersic_components` counts all fitted galaxies |

---

## What “physical” looks like

At per-field `J)` from `zp_aper_40px` in `galfit.feedme`:

| Parameter | Typical OK range | Red flag |
|-----------|------------------|----------|
| `mag` | **8 – 40** (constraint band) | Pinned at 8 or 40 |
| `re` | ~1.5 – 50 px | **≥ 99 px** (constraint max **100**) |
| `n` | 0.5 – 4 | **≥ 5.95** (ceiling **6.0**) |
| `b/a` | 0.2 – 0.9 | **≤ 0.12** → inclination → **90°** |
| `chi2nu` | ~0.4 – 15 (multi-Sérsic higher) | **> 50** with bad residuals |

**Note:** High χ²/ν on multi-Sérsic / large-ROI stamps is often dominated by undetected field galaxies and extra fitted components, not necessarily a bad host `b/a`.

---

## Per-FRB checklist

Under `pipeline_scripts/Output/<FRB>_all/`:

| # | File | What to check |
|---|------|----------------|
| 1 | `qa_cutout_mask.png` | Correct host? ROI too large? (see **P7** in `tasks.md`) |
| 2 | `host_components.csv` | Row 0 = host |
| 3 | `galfit.feedme` | `J)` = field ZP; `3) mag` seeds; component count |
| 4 | `constraints.txt` | `mag 8.0 to 40.0` per Sérsic |
| 5 | `fit.log` | **First** `sersic` line = host; sky; χ²/ν |
| 6 | `sky_fit_audit.json` | `galfit_pass1_ok`, `passed`, `failure_reason` |
| 7 | `galfit_results.png` | Model vs residuals on host |

---

## Likely failure modes

1. **ROI too large (P7):** containment expansion captures excess field → extra Sérsics, high χ²/ν.
2. **Multi-Sérsic:** Host still component 1; use host row for inclination.
3. **GALFIT crash:** `failure_reason: galfit_crash_pass1` — check `fit.log`, mag constraints, component count.
4. **Sky QA:** `sky_qa_failed` only when parsed sky drifts > 3 ADU after a successful fit.
5. **Sigma / invvar:** check `host_sigma` rescaling in `master_run.log`.
6. **Parser block:** `last_block_first_sersic` is normal for multi-Sérsic; host is still line 0 in that block.

---

## Pipeline vs master magnitude (context)

Diff table **excludes mag**. Pipeline uses per-field `zp_aper_40px`; legacy master mixed `J)` systems.

---

## Related scripts

| Script | Role |
|--------|------|
| `scripts/compare_pipeline_galfit_vs_master.py` | Results + diff CSVs (all 39 FRBs, host component 1) |
| `scripts/flag_pipeline_unphysical_fits.py` | Heuristic QA |
| `scripts/galfit_fitlog_parse.py` | Shared fit.log parser |
| `pipeline_scripts/galfit_fitting/run_galfit_fitting.py` | Phase 3b |
| `pipeline_scripts/galfit_fitting/generate_galfit_cutouts.py` | Phase 3a |

---

## Open tasks

See `tasks.md` — **P7** (low priority): rework Phase 3a ROI extension / padding.
