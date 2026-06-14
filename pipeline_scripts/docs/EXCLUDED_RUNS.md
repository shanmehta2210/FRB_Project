# Excluded / experimental pipeline runs

Runs documented here are **not** kept under `pipeline_scripts/Output/`.
The production tree holds exactly the **62** hosts in `pipeline_galfit_results.csv`
(one folder per FRB: `<FRB>_all/` with a parseable `fit.log`).

---

## 20220501C (signal localization; marginal AstroPath)

| Item | Value |
|------|--------|
| `coord_semantics` | `signal` (burst position, not published host) |
| Mode tested | `--include-signal --use-astropath-host` |
| Best candidate | objid=14, `sex_number=371` (Phase 2 `image.psf.cat`) |
| Segmentation ID | maps to `NUMBER=349` in Phase 1 `image.cat` (sky match) |
| `posterior_O` / `posterior_U` | **0.670 / 0.330** |
| Sep from burst | ~0.55 arcsec |
| GALFIT (trial) | Completed once (2026-05-25); host mag ~21.8 AB, `chi2nu` ~0.97 |

**Why excluded:** Clears the pipeline floor (`min_astropath_posterior=0.05`) but is
not a secure host association (~33% posterior on “unassociated”). Morphology from
that trial describes a plausible nearby galaxy, not a confirmed host.

**Technical note (fixed 2026-05-25):** AstroPath records `sex_number` from
`image.psf.cat` (Phase 2). Phase 3a must map to Phase 1 `image.cat` / segmentation
`NUMBER` by sky position — direct equality fails when the two passes assign
different IDs to the same source.

**Other signal FRBs (ellipse errors in `master_frb_localization.csv`):**
`20191228A`, `20210117A`, `20210407E`, `20210912A`, `20220918A` — AstroPath
P(O) ≪ 0.05 or no galaxy within 5″; `20230718A`, `20230731A` — no cutout.
None are in the production `Output/` set.

**Do not** fold 20220501C into confident-host or inclination CDF samples without
explicit review. Burst localization confidence in `new_confident_hosts.txt` (1.0)
is **not** AstroPath P(O).

---

## 20221101B (localization vs SExtractor radius)

Real host ~5.25″ from CSV position; default 5″ nearest-galaxy cut can miss it.
Re-run with `--use-astropath-host` or widen `--max-host-sep-arcsec` if revisiting.
See repo-root `New fitnotes.txt`.
