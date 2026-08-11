# Excluded / experimental pipeline runs

Runs documented here are **not** kept under `pipeline_scripts/Output/`.
Production `Output/` holds **64** hosts — 67 as of 2026-07-21, + 3
published-host re-runs kept 2026-07-30 (`20181112A`, `20220501C`, `20220918A`),
− `20240304B` deleted 2026-08-05 (host undetected, below), − the 5 weak
associations quarantined to `Output/_unconfirmed/` the same day (see
[`WEAK_ASSOCIATIONS_PRODUCTION67.md`](WEAK_ASSOCIATIONS_PRODUCTION67.md)).
`pipeline_galfit_results.csv` and `production_confirmed_lit_hosts.csv` are both
64 rows; `scripts/audit_production_outputs.py` asserts folders == results ==
cohort. Hosts with a failed Phase 2 zero-point carry a corrected `mag_final`
(reference-survey or rescaled-ZP magnitude — see
`pipeline_scripts/reference_photometry.py`).

**Weak associations** (no published host lit and/or pipeline \(P(O)\le 0.95\),
or no usable `mag_final`): see
[`WEAK_ASSOCIATIONS_PRODUCTION67.md`](WEAK_ASSOCIATIONS_PRODUCTION67.md).
As of 2026-08-05 those folders are moved to `Output/_unconfirmed/`, not deleted —
restore and re-run if the host is later published.

---

## Published-host re-run batch (2026-07-30)

Tried **20** FRBs that (a) have `coord_semantics=host`, (b) have a
**published** localization cite (`has_published_localization_ref=yes`),
(c) are **not** Verdi+2025 in-prep / no-publication rows, (d) are **not**
ATel-only (`20230814B` skipped), and (e) have paired cutouts — using
`--use-localization-host` under the `mag_final` framework.

| Outcome | FRBs |
|---------|------|
| **Kept** (fit + `mag_final`) | `20181112A`, `20220501C`, `20220918A` |
| Fit but **no** `mag_final` → deleted | `20180301A`, `20201124A` |
| No `fit.log` / Phase 3 fail → deleted | the other 15 (see `_rerun_eval_20260730.csv`) |

List / log / eval: `_rerun_published_hosts_20260730.txt`,
`_rerun_batch_20260730.log`, `_rerun_eval_20260730.csv` in this `docs/` folder.

---

## GTC Jun-2026 trial folders removed (no `fit.log`)

Deleted from `pipeline_scripts/Output/` (2026-07-11) because host fitting never
produced a usable `fit.log` (Phase 2/3a failures or blank hosts). Also purged
from `pipeline_galfit_results.csv` / diff refresh:

`20121102A`, `20210117A`, `20210809C`, `20220204A`, `20220506D`, `20230501A`,
`20230521A`, `20230521B`, `20230814B`

(Several of these were re-tried 2026-07-30 with the mag framework; none of that
retry set produced a keepable fit except where listed above.)

See `GTC data/pipeline_trial/cohort_manifest.csv` for per-FRB failure notes.

Also removed leftover selective-`--outputs` shells for `20121102A`
(`*_astropath`, `*_catalog_psf`, `*_galfit`, `*_statmorph`) that only contained
staged `.workdir` inputs and never produced a proper `<FRB>_all/` tree.
`20121102A` remains `host_missing` in pipeline trials. Inputs remain at
`large_cutouts/20121102A_*.fits`.

---

## 20220501C (historical note; now in production)

| Item | Value |
|------|--------|
| `coord_semantics` | `host` (Shannon et al. 2025, PASA; published host localization) |
| 2026-05 trial | `--include-signal --use-astropath-host` → \(P(O)\approx0.67\) |
| 2026-07-30 | Re-run with `--use-localization-host` → **kept** (`mag_final` OK, \(b/a=0.44\)) |

Earlier exclusion was for an AstroPath-only signal trial. The host is now on
disk in `Output/20220501C_all/`. Pipeline AstroPath \(P(O)\) is still ~0.67;
use the published host position for science, not the AstroPath pick alone.

**Technical note (fixed 2026-05-25; SPREAD lookup hardened 2025-06):** AstroPath records `sex_number` from
`image.psf.cat` (Phase 2). Phase 3a must map to Phase 1 `image.cat` / segmentation
`NUMBER` by sky position — direct equality fails when the two passes assign
different IDs to the same source. All SPREAD lookups (`_spread_for_catalog_index`,
neighbor masking in the ROI loop) now sky-match; never index Phase 2 SPREAD by
Phase 1 `NUMBER` alone.

---

## 20240304B — host undetected (removed 2026-08-05)

`Output/20240304B_all/` deleted; inputs kept at `large_cutouts/20240304B_*.fits`.

| Item | Value |
|------|-------|
| Localization | MeerKAT, `host` semantics, **z = 2.148** (Caleb et al. 2025, arXiv:2508.01648) |
| Phase 3a (2026-08-05) | **exit 1** — *"No galaxy (SPREAD+3·SPREADERR ≥ 0.005) within 5.0″ of the target position"* |
| Reference photometry | `no_match` — no LS DR10 **and** no PS1 source within 3.0″ |
| Prior fit (2026-05-23) | `SNR_WIN=0.83`, `re=1.5` (at constraint floor), `b/a=0.06`, `inc=90±45°`, `mag=23.0±1.07` |

The old fit was a Sérsic fit to noise. It predates the SPREAD sky-match hardening
and the Phase-2→Phase-1 segmentation ID mapping fix, so the May run most likely
locked onto a spurious or stellar detection that the current star/galaxy cut
rejects. At z = 2.148 — the highest-redshift host in the sample — a
non-detection in this imaging is the expected outcome, and the Phase 3a refusal
is correct behaviour, not a pipeline regression.

Revisit only with deeper imaging; do **not** widen `--max-host-sep-arcsec` to
force an association here.

---

## 20221101B (localization vs SExtractor radius)

Real host ~5.25″ from CSV position; default 5″ nearest-galaxy cut can miss it.
Re-run with `--use-astropath-host` or widen `--max-host-sep-arcsec` if revisiting.
See `Archive/notes/New fitnotes.txt`.
