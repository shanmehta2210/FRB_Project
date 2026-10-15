# New-host cohort — master log (46 FRBs)

Last consolidated: **2026-10-10 23:59:54 UTC** (`python scripts/consolidate_new_hosts_logs.py`)

Machine-readable table: `pipeline_scripts/new_hosts_master.csv`

This file replaces scattered `new_hosts_*.txt`, `new_hosts_pipeline_status.csv`,
`new_hosts_pipeline_batch.log`, and `large_cutouts/PROGRESS.md` for the 46-host cohort.
Cutout inventory for all on-disk cutouts (including non-cohort): `large_cutouts/cutout_registry.csv`.

## Survey / fetch policy

- r-band 10′ cutouts (2290 px @ 0.262″/px): **Legacy → PS1 → DES** (no SkyMapper).
- One FRB at a time: `python scripts/cutout_download.py <FRB>`
- Batch pipeline: `python pipeline_scripts/run_all_frbs.py --use-localization-host --list-file pipeline_scripts/new_hosts_master.csv`
  (use column `frb`; or filter CSV to `paired_cutout_on_disk == True`).

## Summary

| Metric | Count |
|--------|------:|
| Cohort FRBs | 47 |
| Paired flux+invvar on disk | 42 |
| No survey coverage (manual cutout) | 5 |
| Pipeline complete | 24 |
| Pipeline partial | 1 |
| Pipeline host_missing | 0 |

## Cohort table (key columns)

| FRB | Dec | Cutout | Pipeline | P(O) max | Nearest src ″ | Batch | Notes |
|-----|-----|--------|----------|---------:|--------------:|-------|-------|
| 20121102A | +33.1 | ok | missing | — | — | not_in_batch |  |
| 20180301A | +4.7 | ok | missing | — | — | not_in_batch |  |
| 20190523A | +72.5 | ok | complete | 0.985 | — | not_in_batch |  |
| 20190614D | +73.7 | ok | missing | — | — | not_in_batch |  |
| 20201123A | -50.8 | no_coverage | missing | — | — | skipped_no_cutout |  |
| 20201124A | +26.1 | ok | missing | — | — | not_in_batch |  |
| 20220204A | +69.7 | ok | missing | — | — | not_in_batch |  |
| 20220222C | -28.0 | ok | missing | — | — | not_in_batch |  |
| 20220224C | -22.9 | ok | missing | — | — | not_in_batch |  |
| 20220418A | +70.1 | ok | complete | 0.936 | — | not_in_batch |  |
| 20220506D | +72.8 | ok | missing | — | — | not_in_batch |  |
| 20220717A | -19.3 | ok | complete | 0.992 | — | not_in_batch |  |
| 20220726A | +69.9 | ok | missing | — | — | not_in_batch |  |
| 20221029A | +72.5 | ok | missing | — | — | not_in_batch |  |
| 20221101B | +70.7 | ok | complete | 0.991 | — | not_in_batch |  |
| 20221113A | +70.3 | ok | missing | — | — | not_in_batch |  |
| 20221116A | +72.7 | ok | complete | 0.973 | — | not_in_batch |  |
| 20221219A | +71.6 | ok | complete | 0.953 | — | not_in_batch |  |
| 20230124A | +71.0 | ok | complete | 0.994 | — | not_in_batch |  |
| 20230125D | -31.5 | no_coverage | missing | — | — | skipped_no_cutout |  |
| 20230307A | +71.7 | ok | complete | 0.996 | — | not_in_batch |  |
| 20230501A | +70.9 | ok | missing | — | — | not_in_batch |  |
| 20230521B | +71.1 | ok | missing | — | — | not_in_batch |  |
| 20230613A | -27.1 | ok | complete | 0.997 | — | not_in_batch |  |
| 20230626A | +71.1 | ok | complete | 0.998 | — | not_in_batch |  |
| 20230628A | +72.3 | ok | complete | 0.916 | — | not_in_batch |  |
| 20230712A | +72.6 | ok | complete | 0.983 | — | not_in_batch |  |
| 20230718A | -40.5 | no_coverage | missing | — | — | skipped_no_cutout |  |
| 20230731A | -56.8 | no_coverage | missing | — | — | skipped_no_cutout |  |
| 20230814B | +73.0 | ok | missing | — | — | not_in_batch |  |
| 20230907D | +8.7 | ok | complete | 0.992 | — | not_in_batch |  |
| 20230913G | +70.8 | ok | complete | 1.000 | — | not_in_batch |  |
| 20230930A | +41.4 | no_coverage | partial | — | — | skipped_no_cutout |  |
| 20231020B | -37.8 | ok | complete | 0.745 | — | not_in_batch |  |
| 20231120A | +73.3 | ok | complete | 0.999 | — | not_in_batch |  |
| 20231123B | +70.8 | ok | complete | 0.989 | — | not_in_batch |  |
| 20231220A | +73.7 | ok | complete | 0.995 | — | not_in_batch |  |
| 20240104A | +72.8 | ok | missing | — | — | not_in_batch |  |
| 20240119A | +71.6 | ok | complete | 0.317 | — | not_in_batch |  |
| 20240123A | +71.9 | ok | missing | — | — | not_in_batch |  |
| 20240203D | +73.9 | ok | complete | 1.000 | — | not_in_batch |  |
| 20240213A | +74.1 | ok | complete | 0.840 | — | not_in_batch |  |
| 20240215A | +70.2 | ok | complete | 0.988 | — | not_in_batch |  |
| 20240224A | +73.5 | ok | complete | 0.998 | — | not_in_batch |  |
| 20240229A | +70.7 | ok | complete | 0.899 | — | not_in_batch |  |
| 20240304B | +11.8 | ok | missing | — | — | not_in_batch |  |
| 20250518 | +71.3 | ok | missing | — | — | not_in_batch |  |

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

## Batch run log

When `new_hosts_pipeline_batch.log` exists, consolidate copies it to `new_hosts_batch_log_archive.txt` and embeds the full text below. Historical 41-host batch (May 2026, `--use-localization-host`): 41 FRBs, 2112 s total, 0 failed — per-FRB timing and ZP details were in the deleted sidecar log; re-run batch or check git history if you need the verbatim file.

```
[batch] Reconstructed summary from new_hosts_master.csv (original new_hosts_pipeline_batch.log was consolidated then removed).


[batch] 0 FRB(s) with batch metadata in master CSV.
```
