# GTC pipeline trial cohort (13 FRBs)

Last consolidated: **2026-06-24 16:56:11 UTC** (`python "GTC data/pipeline_trial/consolidate_trial_cohort.py"`)

Machine-readable table: `GTC data/pipeline_trial/cohort_manifest.csv`

## Purpose

Archival imaging trial runs for GTC-target FRBs that still need host imaging.
**Update 2026-07-21:** the 5 trial fits that completed (`20210214G`,
`20221116A`, `20230913`, `20230930A`, `20240203`) were **accepted into the
production `pipeline_galfit_results.csv` cohort (now 67 hosts)**; the 8
no-`fit.log` trial folders were deleted from `pipeline_scripts/Output/`
(see `pipeline_scripts/docs/EXCLUDED_RUNS.md`).
Mark bad fits in `cohort_manifest.csv` (`fit_disposition=bad_fit`) then re-run
consolidate; removable paths are listed in `excluded_bad_fits.csv`.

## Selection (frozen)

1. In `master_frb_localization_needs_imaging.csv` (no prior pipeline GALFIT fit)
2. Not in `visibility/exclude_pipeline_fitted_frbs.csv` (production 62)
3. Passes >=1 rigorous GTC night, 2026-06-24 .. 2026-07-24
4. Paired `large_cutouts/{FRB}_flux.fits` + `_invvar.fits` on disk

## Run pipeline batch

```bash
python pipeline_scripts/run_all_frbs.py \
  --list-file "GTC data/pipeline_trial/frb_list.txt" \
  --include-signal \
  --skip-existing
```

(`--include-signal` required: 4/13 have `coord_semantics=signal`.)

## Summary

| Metric | Count |
|--------|------:|
| Cohort FRBs | 13 |
| Paired cutouts | 13 |
| Pipeline complete (fit.log + galfit_results.png) | 5 |
| Fit disposition pending review | 13 |
| Marked bad_fit / exclude | 0 |

## Cohort table

| FRB | Dec | Semantics | Cutout | GTC nights | Pipeline | Disposition | Preview |
|-----|-----|-----------|--------|------------|----------|-------------|---------|
| 20210117A | -16.15 | host | nan | 23 | host_missing | pending | 20210117A_1arcmin.png |
| 20210214G | -5.83 | signal | nan | 17 | complete | pending | 20210214G_1arcmin.png |
| 20210809C | +1.33 | signal | nan | 31 | host_missing | pending | 20210809C_1arcmin.png |
| 20220204A | +69.72 | host | legacy | 31 | partial | pending | 20220204A_1arcmin.png |
| 20220506D | +72.83 | host | ps1 | 31 | host_missing | pending | 20220506D_1arcmin.png |
| 20221116A | +72.65 | host | ps1 | 8 | complete | pending | 20221116A_1arcmin.png |
| 20230501A | +70.92 | host | ps1 | 31 | host_missing | pending | 20230501A_1arcmin.png |
| 20230521A | -2.39 | signal | nan | 26 | host_missing | pending | 20230521A_1arcmin.png |
| 20230521B | +71.14 | host | nan | 31 | host_missing | pending | 20230521B_1arcmin.png |
| 20230814B | +73.03 | host | ps1 | 31 | host_missing | pending | 20230814B_1arcmin.png |
| 20230913 | +70.79 | host | ps1 | 31 | complete | pending | 20230913_1arcmin.png |
| 20230930A | +41.42 | host | ps1 | 31 | complete | pending | 20230930A_1arcmin.png |
| 20240203 | +73.90 | host | ps1 | 31 | complete | pending | 20240203_1arcmin.png |
