# Jimin field galaxy counts — verified (DR16)

**Full write-up (counts + VizieR COSMOS cross-check + CDFs): see [`ANALYSIS.md`](ANALYSIS.md).**  
**Production SDSS v2 vs Jimin cuts: see [`PROD_VS_JIMIN.md`](PROD_VS_JIMIN.md).**

**Field:** RA ∈ [148, 152], Dec ∈ [0, 4] → **16 deg²** (fully covered: 9 runs, 123 fields, Dec spans 0→4).  
**Release:** SDSS DR16 `PhotoObj` + `Photoz` via astroquery.  
**Script:** `scripts/_jimin_field_count_verify.py` → `field_count_verify.csv`.  
**Re-run date:** 2026-07-24 (re-counted from scratch; matches prior funnel).

## Mag cut clarification

Historical Jimin SQL: `p.r BETWEEN 12 AND 21`.

| Cut | N (type_r=3, mode=1, clean=1) |
|-----|-------------------------------:|
| `p.r BETWEEN 12 AND 21` | **47,351** |
| `p.r < 21` | 47,350 |
| `p.r <= 21` | 47,351 |
| `p.petroMag_r BETWEEN 12 AND 21` | 44,599 |

So “mag < 21” vs inclusive ≤21 / BETWEEN 12–21 is **the same to ±1 object** here.  
`p.r` ≡ **modelMag_r** (not Petrosian); probe: `max|r − modelMag_r| = 0`.

## Does ~3000 galaxies / deg² to r≈21 make sense?

Yes. Yasuda et al. (2001), [astro-ph/0105545](https://arxiv.org/abs/astro-ph/0105545), Table 2: differential galaxy counts in **r\*** (Petrosian, extinction-corrected), Northern equatorial stripe. Summing \(N_{r^*}\) from 12.0–21.0:

| Source | Cumulative density (12 ≲ r ≲ 21) |
|--------|----------------------------------:|
| Yasuda Table 2 (sum of \(N_{r^*}\)) | **≈ 3370 deg⁻²** |
| This field, `type_r=3` + **model** `r` ∈ [12,21] + mode/clean | **2959 deg⁻²** (88% of Yasuda) |
| This field, same + **petroMag** ∈ [12,21] | **2787 deg⁻²** (83% of Yasuda) |
| Expected in 16 deg² @ Yasuda density | ≈ 54k |
| Observed (model, pre–extra cuts) | **47,351** |

Slightly below Yasuda is expected: we impose `clean=1`, use observed (not dereddened) mags, and a small 16 deg² patch has cosmic variance. Order of magnitude and fraction-of-literature are solid — **not** an undercount of the field.

(Do **not** divide Yasuda’s “900,000 galaxies” by ~440 deg² ≈ 2045/deg²: their bins use different areas, so that ratio is not the cumulative density.)

Without a mag cut, `type_r=3` alone is **9004 deg⁻²** (144k in the box) — includes objects fainter than 21 down to the detection limit; Yasuda’s 20.5–21.0 bin alone is already ~1216 deg⁻² differential, so totals well above the <21 sample are normal.

## Full cut funnel (modelMag / historical)

| Stage | N | / deg² | Dropped vs previous | What it removes |
|-------|--:|-------:|--------------------:|-----------------|
| All `PhotoObj` in box | 653,604 | 40,850 | — | Everything (multi-epoch detections, stars, junk) |
| + `mode=1` (primary) | 386,868 | 24,179 | 266,736 (41%) | Duplicate detections of same object |
| + `clean=1` | 250,395 | 15,650 | 136,473 (35%) | Photometric flags / unclean photometry |
| + `type_r=3` (galaxy in r) | 144,070 | 9,004 | 106,325 (42%) | Stars / non-galaxies (`type=3` gives 149,555 — nearly same) |
| + **`r` ∈ [12, 21]** | **47,351** | **2,959** | **96,719 (67%)** | **Faint (and rare bright) galaxies — the mag limit** |
| + `lnLStar_r < -10` | 33,045 | 2,065 | 14,306 (30%) | Star-like / PSF-preferring objects |
| + Photoz `nnAvgZ > 0` | 32,978 | 2,061 | 67 (0.2%) | No usable photo-z row / non-positive z |
| + `score > 0.8` | **28,043** | **1,753** | 4,935 (15%) | Low photo-z “score” (Jimin BASE, no morph) |
| + `lnLDeV < lnLExp` (**V2**) | **17,657** | **1,104** | 10,386 (37%) | deV-preferred (bulge-like) profiles |
| + `fracDeV_r = 0` (**V1**) | **9,120** | **570** | 8,537 (48% of V2) | Any cModel bulge light (`fracDeV>0`) |
| V2 + `expAB_r > 0.2` | 16,551 | 1,034 | 1,106 from V2 | Edge-on / floor b/a ≤ q0 |
| V1 + `expAB_r > 0.2` | 8,445 | 528 | 675 from V1 | Same |

### Extra checks

- **`fracDeV_r = 0` alone** (no explicit lnL): **9,120** = V1. So every pure-cModel-exp object already has lnLExp winning; the AND lnL clause is redundant once `fracDeV=0` is required.
- Mag cut is the **largest single drop** after star/galaxy separation (67% of type_r=3 galaxies are outside 12–21, almost all faint).
- Morphology (lnL then fracDeV) removes another **~⅔ of the BASE** (28k → 9k for V1).

## Bottom line

1. Field density to **r ≈ 21** is **~3000 galaxies deg⁻²** → **~47k** in 16 deg² — consistent with Yasuda+2001 (~3370 deg⁻²).  
2. V1 ending at **9,120** is **not** a TOP/truncation bug; it is what remains after star reject + Photoz/score + **fracDeV=0**.  
3. V2 full is **17,657**; SkyServer `TOP 10000` would truncate only V2.  
4. Mag cut uses **modelMag (`p.r`)**, not Petrosian; swapping to petro **lowers** N slightly.
