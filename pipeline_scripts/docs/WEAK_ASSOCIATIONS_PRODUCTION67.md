# Weak host associations inside the production 67

**Maintained note — last updated 2026-08-04.**

> **Cohort:** production `Output/` holds **69** hosts (`20240304B` removed
> 2026-08-05, host undetected at z = 2.148 — see [`EXCLUDED_RUNS.md`](EXCLUDED_RUNS.md));
> the five in §B are excluded
> from the **64**-host confirmed-literature cohort used for production re-runs and
> science samples. Regenerate that list with
> `python scripts/build_confirmed_lit_cohort.py` →
> [`production_confirmed_lit_hosts.csv`](../../production_confirmed_lit_hosts.csv).
> Because the §B five are not re-run, their `Output/<FRB>_all/` trees can lag the
> current pipeline — check `cutout_meta.json` `neighbor_policy` before reuse.

These FRBs appear in `pipeline_galfit_results.csv` (production `Output/<FRB>_all/`)
but **must not** be treated as secure hosts in inclination CDFs or host-property
samples without an explicit override. Two disjoint reasons are tracked below.

Machine-readable audit of all 67: [`production67_lit_astropath_audit.csv`](../../production67_lit_astropath_audit.csv)
(root). Related exclusions (never in `Output/`): [`EXCLUDED_RUNS.md`](EXCLUDED_RUNS.md).

---

## A. No usable host magnitude (3 FRBs)

**3 / 67** hosts have **no survey r-band magnitude** (LS DR10 footprint miss +
PS1 DR1 too shallow / no match within 3″). Under the `mag_final` policy
(`pipeline_scripts/reference_photometry.py`) they are `mag_final_source =
unavailable` with empty `mag_final` — the raw GALFIT `mag` column is **not**
trusted for analysis (failed Phase 2 ZP and/or `mag_err > 1`).

| FRB | GALFIT `mag`±`mag_err` | Why no survey mag | Notes |
|-----|------------------------|-------------------|-------|
| **20221116A** | 17.25 ± 2.23 | High-dec PS1 field; nearest PS1 ~8.5″, LS DR9/10 empty | Sharma+2024 published host; ZP broken (`zp_aper_std` ≫ 1); `mag_err > 1` |
| **20230913** | 14.04 ± 0.11 | Host below PS1 catalog depth; LS empty | Verdi+2025 (in prep.); Phase 2 failed |
| **20240104A** | 17.84 ± 0.02 | Bright fitted source has **no** PS1 entry within 10″; LS empty | Verdi+2025 (in prep.); fitted object CLASS_STAR≈0.98, 6.9″ from DSA loc — may be a foreground star |

**Policy:** do not use these three for any magnitude-dependent science. Structural
parameters (re, n, b/a, PA, inclination) from `fit.log` remain on disk but the
host association / identity of **20240104A** and **20230913** is also weak
(§B). Prefer dropping all three from analysis samples that need a host mag.

---

## B. Weak published / AstroPath association (5 FRBs)

These do **not** meet the association bar used for science:

- **published** host-localization literature, **or**
- pipeline AstroPath \(P(O\mid x) > 0.95\).

| FRB | `coord_semantics` | Why no published lit | Pipeline AstroPath \(P(O)\) | Notes |
|-----|-------------------|----------------------|----------------------------|-------|
| **20210214G** | `signal` | Multi-beam / no secure host paper (Shannon+2025: voltages disabled) | **0.714** | GTC-trial accept 2026-07-21; not a published host |
| **20230913** | `host` | Cite is **Verdi+2025 (in prep.)** only | *none* (Phase 2 failed; loc-mode fit) | Also in §A (no `mag_final`) |
| **20240104A** | `host` | Verdi+2025 (in prep.) only | ~0 (\(P(U)=1\)) | Also in §A; TeX \(P_{\rm host}=0.99\) (in prep.) |
| **20240203** | `host` | Verdi+2025 (in prep.); not in Hussaini+2025 | *none* (Phase 2 failed; loc-mode fit) | TeX \(P_{\rm host}=0.996\) (in prep.); has PS1 `mag_final` |
| **20250518** | `host` | Verdi+2025 (in prep.) only | **0.896** | TeX \(P_{\rm host}=0.99\) (in prep.); pipeline best \(P(O)&lt;0.95\) |

Also among the five GTC accepts: **20221116A** and **20230930A** *do* have published
cites (Sharma+2024 / etc.) for association — but **20221116A** is still barred
from mag-dependent analysis by §A.

---

## Policy

1. **Do not** fold §A or §B into confident-host or strict inclination CDF samples
   by default. §A is mandatory for any cut that needs a host magnitude.
2. If Verdi+2025 (or another paper) **publishes** any of them with a usable host, move
   the FRB out of §B, fill `master_frb_localization.csv` `repeater_source`, and
   refresh `production67_lit_astropath_audit.csv`.
3. Re-running AstroPath with better cutouts / priors may rescue a case only if
   \(P(O)>0.95\) **and** the fitted host matches the intended galaxy — re-audit both.
4. **Quarantined (2026-08-05).** The §B five now live in
   `pipeline_scripts/Output/_unconfirmed/<FRB>_all/`, outside the production
   glob, so `compare_pipeline_galfit_vs_master.py` cannot pick them up and no
   stale fit can leak into `pipeline_galfit_results.csv`. Nothing was deleted —
   move a folder back up one level if its host is later published, then re-run
   it before use. Their fits predate the Re-separation ROI
   (`neighbor_policy` is `null` or absent), so they are **not** comparable to
   the current 64.
5. Science magnitudes: use `mag_final` only (never raw GALFIT `mag`). Prefer survey
   when `mag_err > 1` or the pipeline ZP is untrusted.

---

## Changelog

| Date | Change |
|------|--------|
| 2026-07-21 | Initial list from lit vs pipeline AstroPath audit of production 67. |
| 2026-07-21 | Added §A: 3 hosts with no survey / no `mag_final` (`20221116A`, `20230913`, `20240104A`). Mag policy: prefer survey when GALFIT `mag_err > 1`. |
| 2026-07-30 | Published-host re-run: kept `20181112A`, `20220501C`, `20220918A` in `Output/` (production **70**). `20220918A` has `mag_final` but \(b/a=0.10\) (not inclination-reliable). |
| 2026-08-04 | §B formalised as a machine-readable cohort: `scripts/build_confirmed_lit_cohort.py` → `production_confirmed_lit_hosts.csv` (**65** = 70 − §B five). Production batches run that list only. |
| 2026-08-05 | `20240304B` removed from `Output/` (Phase 3a found no galaxy within 5″; z = 2.148 non-detection). Production **69**, cohort **64**. |
| 2026-10-11 | Verdi et al. 2026 (arXiv:2609.30616) published `20230913G`, `20240104A`, `20240203D`, `20240224A` (astropath $P>0.9$, Keck/LRIS $z$). `20230913` → `20230913G` and `20240203` → `20240203D` promoted out of §B on fresh PS1 DR2 cutouts centred on the PS1 host (the earlier PS1 cutouts carried a mirrored WCS, fixed in `scripts/cutout_resample.py`); `20240224A` added. Phase 2 now calibrates on 86–233 PS1 stars, so none of the three is in §A. Production **67**, cohort **67**. `20240104A` stays quarantined (z = 1.33; fitted object may be a star). Triage: `20230913G`, `20240203D` rejected (unresolved, $R_e$ floor); `20240224A` confirmed (quiescent-host caveat). |
| 2026-10-11 | Re-ran the four pre-existing PS1 production hosts (`20220717A`, `20221101B`, `20221116A`, `20230930A`) on reprojected cutouts. Phase 2 now calibrates on 313-993 PS1 stars (ZP scatter 0.07-0.23 mag; previously 3-7 stars, scatter up to 52 mag). `mag_final` now comes from the pipeline ZP for `20221101B` (18.13), `20230930A` (16.10), `20220717A` (21.22); `20221116A` leaves §A (PS1 DR1 Kron 20.89 at 5.0″). Shapes essentially unchanged (b/a 0.07→0.05, 0.75→0.74, 0.16→0.14, 0.84→0.78); all four stay rejected in triage. AstroPath P(O) fell for `20220717A` (0.99→0.17) and `20221116A` (0.97→0.26) because the CSV host coordinates are 3-decimal (1.7″ / 5.0″ from the galaxy) and the 1″ default ellipse now applies to a correct WCS — a coordinate-precision issue, not a host-identity change. `Re-fits/20221101B` legs predate this re-run. |
| 2026-10-13 | Added two Gordon et al. (2025) CRAFT hosts identified from the Marnoch+2026 overlap: `20240525A` (LS, $q=0.365$, Tier A, confirmed) and `20240615B` (LS, $q=0.94$, bright face-on spiral, Tier A, confirmed). Production is now 69 hosts; the science sample becomes 54 once statistics are refreshed (55 if `20240203D` recovers). `20191228A` was attempted on its existing LS cutout and rejected (host undetected; $m_r\approx 22.3$). DMs for the two new FRBs are not yet published (Gordon+2025 defer radio properties); `DM`, `DM_MW`, `DM_exgal` left blank in the master list. |
| 2026-10-18 | Repeater policy: the science sample is restricted to non-repeaters. `20220912A` (Ravi+2023 repeater; clean fit) set to `confirmed=False` in `host_confirmation.csv`; it stays in production and the master list. The only other repeater in production, `20190711A`, was already rejected; `20180916B` (rejected) is no longer a candidate on policy grounds. Science sample 56 → **55**. CHIME repeater hosts will not be added. |
| 2026-10-18 | Literature batch of 14 non-repeaters added to master (now 118 rows): nine Dial+2026 CRAFT HTR2 bursts (`20240807A`, `20241014B`, `20250106B`, `20250511A`, `20250607C`, `20251024B`, `20260324E`, `20260331E`, `20260608D`; FRB position + $z$ published, host = our AstroPATH pick, $P\ge0.93$), MeerTRAP `20231210F`/`20231010A`, CHIME Outriggers `20250316A` (NGC 4141), DSA-110 `20220208A`. 12 new production fits (`pipeline_galfit_results.csv` now 84 rows). Triage: 11 confirmed, `20231010A` fails cuts, `20220208A` host undetected in PS1, `20250520A` has no LS/DES imaging. Science sample 55 → **66** pending the statistics refresh. Details in `HOST_TRIAGE_CASES.md`. |
