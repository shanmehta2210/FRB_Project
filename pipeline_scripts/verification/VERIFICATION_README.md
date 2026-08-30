# Verification documentation map

Everything under `pipeline_scripts/verification/` that decides whether a GALFIT
host geometry is trustworthy for the inclination sample. **Start here**, then
follow the doc that matches the job.

## Which doc for which job

| Job | Doc |
|---|---|
| Physics / implementation of the 9 checks | [`FIT_VERIFICATION_CHECKS.md`](FIT_VERIFICATION_CHECKS.md) (+ HTML twin) |
| Science cut (`in_53`) vs full 64 | [`SCIENCE_CUT_AND_COHORT.md`](SCIENCE_CUT_AND_COHORT.md) |
| Confirm / reject hosts; CSV conventions | [`HOST_CONFIRMATION_WORKFLOW.md`](HOST_CONFIRMATION_WORKFLOW.md) |
| Long-form per-FRB case notes | [`HOST_TRIAGE_CASES.md`](HOST_TRIAGE_CASES.md) |
| Independent sky (B / E / F) | [`SKY_PROTOCOL.md`](SKY_PROTOCOL.md) |
| Re-fits, reject grid, fixed-n / fixed-sky / masked-star PSF | [`REFIT_AND_REJECT_GRID.md`](REFIT_AND_REJECT_GRID.md) |
| Deeper-survey / g-band inspection stamps | `fetch_alt_imaging.py` → `alt_imaging/` |
| Manual GALFIT sandbox | [`SANDBOX.md`](SANDBOX.md) |
| Visual panels, stretches, pptx | [`VISUAL_PANELS.md`](VISUAL_PANELS.md) |
| PSFEx / GALFIT / verification PSF normalization | [`PSF_NORMALIZATION.md`](PSF_NORMALIZATION.md) |
| Re-fits staging layout (short) | [`Re-fits/README.md`](Re-fits/README.md) |

**Checks ≠ triage.** The suite measures; `host_confirmation.csv` is the human
gate for the paper. Aggregate trust tiers (A/B/C) are automatic flags, not the
same thing as `confirmed=True`.

## Tree (scripts + products)

```
verification/
├── VERIFICATION_README.md          ← this index
├── FIT_VERIFICATION_CHECKS.md      check bible
├── HOST_TRIAGE_CASES.md            per-FRB narratives
├── HOST_CONFIRMATION_WORKFLOW.md
├── SCIENCE_CUT_AND_COHORT.md
├── SKY_PROTOCOL.md
├── REFIT_AND_REJECT_GRID.md
├── SANDBOX.md
├── VISUAL_PANELS.md
├── PSF_NORMALIZATION.md
├── host_confirmation.csv           machine gate (frb, confirmed, notes)
├── vercommon.py                    shared IO / geometry / cohort
├── run_verification.py             production suite orchestrator
├── run_refit.py                    stage + GALFIT + full suite into Re-fits/
├── fetch_alt_imaging.py            HSC/DES/LS hips2fits 1′ stamps (mag 21–22 + spirals)
├── run_reject_grid.py              n1 / sky / n1_sky grid for rejects
├── run_sandbox.py                  hand-edited feedme GALFIT
├── sky_protocol.py                 independent sky B/E/F
├── aggregate.py                    metrics CSV, flags, plots, contact sheet
├── build_confirmed_pptx.py         one panel per slide deck
├── checks/                         chi2, rff, fourier, psf, mag, isophote,
│                                   sky, astrophot, visual
├── Re-fits/<FRB>/                  staged re-fits + sandboxes
└── outputs/
    ├── per_host/<FRB>/             check JSON + panel.png
    ├── panels/                     published panels (+ confirmed-leg copies)
    ├── tables/                     metrics, flags, population_summary
    ├── plots/
    └── confirmed_fit_panels.pptx
```

Production GALFIT workdirs remain in `pipeline_scripts/Output/<FRB>_all/` —
verification never overwrites them.

## One-command cheatsheet

```bash
# Full suite (all 64), then aggregate
python run_verification.py --checks all --jobs 4

# Visual panels only (after stretch / header changes)
python run_verification.py --checks visual --force --jobs 4 --no-aggregate

# Independent sky for one host
python sky_protocol.py 20181112A

# Re-fit: n=1 + protocol sky, full verification
python run_refit.py 20181112A --sky-from-protocol --fix-n 1 --label n1_sky

# Mag 21–22 deeper r/i + Fourier-spiral g-band inspection stamps
python fetch_alt_imaging.py

# Masked-star PSF companion (host Sérsic kept)
python run_refit.py 20221101B --add-psf-at-masked-star --psf-dmag 2 --label psf

# Reject grid (production copy + three legs)
python run_reject_grid.py

# Manual sandbox
python run_sandbox.py --init                 # seed 20220509G + 20230930A
python run_sandbox.py 20220509G              # run feedme → panel.png

# Confirmed-sample PowerPoint
python build_confirmed_pptx.py
```

See [`FIT_VERIFICATION_CHECKS.md` §2](FIT_VERIFICATION_CHECKS.md) for suite CLI
details and [`REFIT_AND_REJECT_GRID.md`](REFIT_AND_REJECT_GRID.md) for re-fit
flags.
