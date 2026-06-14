# New-host cohort — master log (46 FRBs)

Last consolidated: **2026-05-24 19:55:14 UTC** (`python scripts/consolidate_new_hosts_logs.py`)

Machine-readable table: `pipeline_scripts/new_hosts_master.csv`

This file replaces scattered `new_hosts_*.txt`, `new_hosts_pipeline_status.csv`,
`new_hosts_pipeline_batch.log`, and `large_cutouts/PROGRESS.md` for the 46-host cohort.
Cutout inventory for all on-disk cutouts (including non-cohort): `large_cutouts/cutout_registry.csv`.

## Survey / fetch policy

- r-band 10′ cutouts (2290 px @ 0.262″/px): **Legacy → PS1 → DES** (no SkyMapper).
- One FRB at a time: `python scripts/cutout_download.py <FRB>`
- Batch pipeline: `python pipeline_scripts/run_all_frbs.py --use-localization-host --list-file pipeline_scripts/new_hosts_master.csv`
  (use column `frb`; or filter CSV to `paired_cutout_on_disk == True`).
- Production `Output/` is the **62** hosts in `pipeline_galfit_results.csv` only; marginal/experimental runs → `pipeline_scripts/docs/EXCLUDED_RUNS.md`.

## Summary

| Metric | Count |
|--------|------:|
| Cohort FRBs | 46 |
| Paired flux+invvar on disk | 41 |
| No survey coverage (manual cutout) | 5 |
| Pipeline complete | 23 |
| Pipeline partial | 14 |
| Pipeline host_missing | 8 |

## Cohort table (key columns)

| FRB | Dec | Cutout | Pipeline | P(O) max | Nearest src ″ | Batch | Notes |
|-----|-----|--------|----------|---------:|--------------:|-------|-------|
| 20121102A | +33.1 | ok | host_missing | 0.000 | 7.3 | not_in_batch |  |
| 20180301A | +4.7 | ok | partial | — | — | not_in_batch | on disk |
| 20190523A | +72.5 | ok | complete | 0.985 | — | not_in_batch |  |
| 20190614D | +73.7 | ok | host_missing | — | 11.5 | not_in_batch |  |
| 20201123A | -50.8 | no_coverage | partial | — | — | skipped_no_cutout | No Legacy/PS1/DES coverage at host position; manual cutout r… |
| 20201124A | +26.1 | ok | partial | — | — | not_in_batch |  |
| 20220204A | +69.7 | ok | host_missing | 0.727 | 6.1 | not_in_batch |  |
| 20220222C | -28.0 | ok | partial | — | 19.2 | not_in_batch | on disk |
| 20220224C | -22.9 | ok | partial | — | 10.1 | not_in_batch | on disk |
| 20220418A | +70.1 | ok | complete | 0.937 | — | not_in_batch |  |
| 20220506D | +72.8 | ok | host_missing | — | 12.2 | not_in_batch |  |
| 20220717A | -19.3 | ok | complete | 0.992 | — | not_in_batch | on disk |
| 20220726A | +69.9 | ok | host_missing | — | 19.1 | not_in_batch |  |
| 20221029A | +72.5 | ok | partial | — | 20.0 | not_in_batch |  |
| 20221101B | +70.7 | ok | complete | 0.991 | — | not_in_batch |  |
| 20221113A | +70.3 | ok | host_missing | — | 17.0 | not_in_batch |  |
| 20221116A | +72.7 | ok | partial | — | — | not_in_batch |  |
| 20221219A | +71.6 | ok | complete | 0.953 | — | not_in_batch |  |
| 20230124A | +71.0 | ok | complete | 0.994 | — | not_in_batch |  |
| 20230125D | -31.5 | no_coverage | partial | — | — | skipped_no_cutout | No Legacy/PS1/DES coverage at host position; manual cutout r… |
| 20230307A | +71.7 | ok | complete | 0.996 | — | not_in_batch |  |
| 20230501A | +70.9 | ok | partial | — | — | not_in_batch |  |
| 20230521B | +71.1 | ok | partial | — | 22.5 | not_in_batch |  |
| 20230613A | -27.1 | ok | complete | 0.997 | — | not_in_batch | on disk |
| 20230626A | +71.1 | ok | complete | 0.998 | — | not_in_batch |  |
| 20230628A | +72.3 | ok | complete | 0.916 | — | not_in_batch |  |
| 20230712A | +72.6 | ok | complete | 0.983 | — | not_in_batch |  |
| 20230718A | -40.5 | no_coverage | partial | — | — | skipped_no_cutout | No Legacy/PS1/DES coverage at host position; manual cutout r… |
| 20230731A | -56.8 | no_coverage | partial | — | — | skipped_no_cutout | No Legacy/PS1/DES coverage at host position; manual cutout r… |
| 20230814B | +73.0 | ok | host_missing | — | 7.3 | not_in_batch |  |
| 20230907D | +8.7 | ok | complete | 0.992 | — | not_in_batch | on disk |
| 20230913 | +70.8 | ok | partial | — | — | not_in_batch |  |
| 20230930A | +41.4 | no_coverage | missing | — | — | skipped_no_cutout | No Legacy/PS1/DES coverage at host position; manual cutout r… |
| 20231020B | -37.8 | ok | complete | 0.746 | — | not_in_batch | on disk |
| 20231120A | +73.3 | ok | complete | 0.999 | — | not_in_batch |  |
| 20231123B | +70.8 | ok | complete | 0.989 | — | not_in_batch |  |
| 20231220A | +73.7 | ok | complete | 0.995 | — | not_in_batch |  |
| 20240104A | +72.8 | ok | complete | 0.000 | 6.9 | not_in_batch |  |
| 20240119A | +71.6 | ok | complete | 0.319 | — | not_in_batch |  |
| 20240123A | +71.9 | ok | host_missing | — | 5.3 | not_in_batch |  |
| 20240203 | +73.9 | ok | partial | — | — | not_in_batch |  |
| 20240213A | +74.1 | ok | complete | 0.839 | — | not_in_batch |  |
| 20240215A | +70.2 | ok | complete | 0.988 | — | not_in_batch |  |
| 20240229A | +70.7 | ok | complete | 0.900 | — | not_in_batch |  |
| 20240304B | +11.8 | ok | complete | 0.000 | 11.1 | not_in_batch |  |
| 20250518 | +71.3 | ok | complete | 0.896 | — | not_in_batch |  |

## List tags (from former sidecar files)

- **high_north_ps1** — Dec ≳ +70°; PS1-only footprint for calibration.
- **former_skymapper** — originally SkyMapper cutouts; re-download with Legacy ladder.
- **no_coverage** — no Legacy/PS1/DES at host; see notes column.

## No-cutout FRBs (detail)

- **20201123A** — No Legacy/PS1/DES coverage at host position; manual cutout required.
- **20230125D** — No Legacy/PS1/DES coverage at host position; manual cutout required.
- **20230718A** — No Legacy/PS1/DES coverage at host position; manual cutout required.
- **20230731A** — No Legacy/PS1/DES coverage at host position; manual cutout required.
- **20230930A** — No Legacy/PS1/DES coverage at host position; manual cutout required.

## Full batch run log (preserved)

Source: former `new_hosts_pipeline_batch.log` (41-host batch, `--use-localization-host`).

```
[batch] Reconstructed summary from new_hosts_master.csv (original new_hosts_pipeline_batch.log was consolidated then removed).


[batch] 0 FRB(s) with batch metadata in master CSV.
```
