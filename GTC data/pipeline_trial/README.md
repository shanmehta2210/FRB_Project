# GTC pipeline trial cohort

Trial imaging-pipeline runs for **13 FRBs** selected for GTC follow-up. Kept
separate from the production 62-host catalog so bad or marginal fits can be
removed without polluting `pipeline_galfit_results.csv`.

## Files

| File | Purpose |
|------|---------|
| `cohort_manifest.csv` | Master record: selection criteria, cutout/GTC metadata, pipeline status, fit QA disposition |
| `cohort_manifest.md` | Human-readable summary (auto-generated) |
| `frb_list.txt` | One FRB per line — pass to `run_all_frbs.py --list-file` |
| `excluded_bad_fits.csv` | FRBs marked `bad_fit` or `exclude` in the manifest (safe-removal list) |
| `consolidate_trial_cohort.py` | Refresh manifest from disk after pipeline runs or QA edits |

## Workflow

1. Run pipeline trials:
   ```bash
   python pipeline_scripts/run_all_frbs.py \
     --list-file "GTC data/pipeline_trial/frb_list.txt" \
     --include-signal \
     --skip-existing
   ```
2. Refresh records:
   ```bash
   python "GTC data/pipeline_trial/consolidate_trial_cohort.py"
   ```
3. After visual QA, edit `cohort_manifest.csv`:
   - `fit_disposition`: `pending` → `keep` | `bad_fit` | `exclude`
   - `fit_disposition_reason`: short note (required for bad_fit/exclude)
4. Re-run consolidate; check `excluded_bad_fits.csv` for paths to remove from `Output/`.

**Rule:** Trial outputs stay out of production until `exclude_from_production` is
manually cleared and the FRB is promoted through the normal catalog workflow.
