# Audit verification — NULL_CATALOG_DATA_AND_INCLINATION_AUDIT.md

**Date:** 2026-05-28  
**Catalog:** `catalog/SDSS_catalog_v1_allsky_modelmr.csv`  
**Production CDF mag limit:** `modelMag_r < 21`

## Automated checks

| Check | Result |
|-------|--------|
| `python scripts/audit_sdss_exp_winner_pool.py --mag-limit 21` | **PASS** — CDF pool N=27,409; **0** lnL violations |
| `python scripts/test_null_catalog_sanity.py` | **PASS** — SDSS strict+color+morph @ m<21: N=27,409 (min 10,000) |

## Row counts vs audit document

| Metric | Audit doc | Verified (2026-05-28) | Match |
|--------|-----------|------------------------|-------|
| Catalog rows | ~457k | **457,190** | Yes |
| Finite lnL pairs | ~87% | **86.7%** (396,327) | Yes |
| CDF pool @ m<21 | (§6 ~10k at mag20; production uses 21) | **27,409** | Expected (mag 21) |
| lnL violations in pool | 0 | **0** | Yes |
| Median cos(i) @ production pool | §6 mag20 ~0.61 | **0.579** @ m<21 | Consistent (mag 21 slightly lower) |

## Section-by-section

| Section | Verdict | Notes |
|---------|---------|-------|
| §1 Reference line | **PASS** | Matches `plot_inclination_cdf_overlay` + `cdf_envelope` |
| §2 SDSS download / SQL | **PASS** | Matches `build_sdss_null_catalog.py`; `mode=1` not in SQL (doc correct) |
| §2 Profile winner | **PASS** | `sdss_exp_wins_lnl_mask` strict `lnLExp_r > lnLDeV_r` |
| §3 Hubble cos(i) | **PASS** | `Q0=0.2` in `null_catalog_utils.py` |
| §4 CDF pipeline order | **PASS** | `prepare_null_strict_color_base` → `slice_null_base_by_mag` (mag last) |
| §5 Stage table | **PARTIAL** | See ordering caveat below |
| §6 mag<20 interpretation | **PASS** | Qualitative story correct; use mag21 N=27,409 for production |
| §7 Sampling checks | **PASS** | Dec footprint, TOP bias, no isotropic draw — all accurate |
| §8 Claims table | **PASS** | lnL applied; uniform diagonal not expected after b/a selection |
| §11 Column contract | **PASS** | `modelMag_r`, `expAB_r`, `u-r<2.3`, lnL |

## §5 ordering caveat (documentation only)

`scripts/decompose_cosi_cdf_bias.py` applies **`modelMag_r < 20` before `u-r`**, so its intermediate N/medians differ from **production order**:

| Stage | Production order (this folder) | Decompose script (audit §5 @ mag20) |
|-------|--------------------------------|-------------------------------------|
| After catalog | N=457,190; med cos i=0.258 | Same |
| After mag cut | (applied last in production) | N=49,157 @ m<20 first |
| After u-r | N=183,294; med=0.228 | N=14,498 @ m<20 path |
| After lnL | N=114,899; med=0.271 | N=10,490 |
| After strict b/a | N=71,424; med=0.494 | N=9,581 |
| Final + mag | **N=27,409** @ m<21 | — |

Production-pool N at m<21 equals `prepare_null_strict_color_base` + `slice_null_base_by_mag` (**27,409**).

## Shape distribution spot checks

| Claim (§5.1) | Verified |
|--------------|----------|
| ~17% at expAB_r = 0.05 | **16.5%** |
| ~38% at expAB_r ≤ 0.2 | **38.2%** |

## Overall

**The audit document reflects the true code and implementation.** Update §5 intermediate table labels if you want numbers that match production cut order; use `diagnostics/sdss_cut_evolution/stage_summary.csv` for production-order statistics.

## Evolution vs production overlay (same pool, different curve)

At `modelMag_r < 21`, the **galaxy table is identical** (N = 27,409; same median cos(i) ≈ 0.579).

| Plot type | What is drawn |
|-----------|----------------|
| `sdss_cut_evolution/cdf_03_mag_lt_21.png` | **Full-pool empirical CDF** (all N in pool) |
| `mag_cuts/.../sdss_strict/null_cdf_inclination.png` | **`cdf_envelope`**: mean of 10⁴ CDFs from **N_FRB = 41** random subsets of that pool |

So the green curve in the FRB overlay is **not** the ECDF of 27,409 galaxies. See `README.md` and `cdf_final_vs_production_overlay_mag21.png`.

**Reproduce:**

```bash
python scripts/audit_sdss_exp_winner_pool.py --mag-limit 21
python scripts/test_null_catalog_sanity.py
python scripts/plot_sdss_null_cut_evolution_cdfs.py
```
