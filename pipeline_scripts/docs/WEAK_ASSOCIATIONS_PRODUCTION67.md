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
