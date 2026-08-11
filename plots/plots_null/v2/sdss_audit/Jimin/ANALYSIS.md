# Jimin SDSS DR16 — field counts, cuts, and VizieR COSMOS cross-check

Verified 2026-07-24. Scripts: `scripts/_jimin_field_count_verify.py`, `scripts/analyze_jimin_vizier_cosmos.py`, `scripts/build_jimin_sdss_replication.py`.

---

## 1. Field and mag cut

| | |
|--|--|
| Footprint (Jimin SQL) | RA ∈ [148, 152], Dec ∈ [0, 4] = **16 deg²** |
| Coverage | Full (9 runs, 123 fields) |
| Release | SDSS DR16 |
| Historical mag SQL | `p.r BETWEEN 12 AND 21` |

`p.r` ≡ **modelMag_r** (not Petrosian). Live probe: `max|r − modelMag_r| = 0`.

| Mag expression | N (`type_r=3`, mode=1, clean=1) |
|----------------|--------------------------------:|
| `r BETWEEN 12 AND 21` | **47,351** |
| `r < 21` | 47,350 |
| `r ≤ 21` | 47,351 |
| `petroMag_r BETWEEN 12 AND 21` | 44,599 |

“mag < 21” vs inclusive ≤21 is the same to ±1 object in this box.

---

## 2. Do ~3000 galaxies / deg² to r≈21 make sense?

Yes. Yasuda et al. (2001) [astro-ph/0105545](https://arxiv.org/abs/astro-ph/0105545) Table 2: sum of differential \(N_{r^*}\) (Petrosian) from 12→21 ≈ **3370 deg⁻²**.

| Source | Density (12 ≲ r ≲ 21) |
|--------|----------------------:|
| Yasuda Table 2 | ≈ 3370 deg⁻² |
| This field, model `r` ∈ [12,21] | **2959 deg⁻²** (88%) |
| This field, petro ∈ [12,21] | 2787 deg⁻² (83%) |
| Expected @ Yasuda × 16 deg² | ≈ 54k |
| Observed pre–extra cuts | **47,351** |

Slightly below Yasuda is expected (`clean=1`, observed not dereddened mags, small-patch variance). Field is **not** undercounted.

Do not use 900k / ~440 deg² ≈ 2045/deg² from the Yasuda abstract — bin areas differ, so that ratio is not cumulative density.

---

## 3. Cut funnel (modelMag, historical Jimin)

| Stage | N | /deg² | Drop | Removes |
|-------|--:|------:|-----:|---------|
| All PhotoObj in box | 653,604 | 40,850 | — | Multi-epoch, stars, junk |
| + `mode=1` | 386,868 | 24,179 | 41% | Duplicates |
| + `clean=1` | 250,395 | 15,650 | 35% | Unclean photometry |
| + `type_r=3` | 144,070 | 9,004 | 42% | Non-galaxies |
| + **`r` ∈ [12,21]** | **47,351** | **2,959** | **67%** | Faint (mag limit) |
| + `lnLStar_r < −10` | 33,045 | 2,065 | 30% | Star-like |
| + Photoz `nnAvgZ > 0` | 32,978 | 2,061 | 0.2% | No usable photo-z |
| + `score > 0.8` | **28,043** | **1,753** | 15% | Low photo-z score |
| + lnLExp wins (**V2**) | **17,657** | **1,104** | 37% | deV-preferred |
| + `fracDeV_r = 0` (**V1**) | **9,120** | **570** | 48% of V2 | Any cModel bulge |
| V2 + `expAB_r > 0.2` | 16,551 | 1,034 | 7% | ba ≤ q0 |
| V1 + `expAB_r > 0.2` | 8,445 | 528 | 7% | ba ≤ q0 |

Notes:
- V1 = 9,120 with or without TOP — not truncation; full COUNT matches.
- `fracDeV=0` alone also = 9,120 (lnL AND is redundant once fracDeV is exact 0).
- V2 without TOP = 17,657; SkyServer `TOP 10000` would truncate V2 only.
- Cos(i) uses **`expAB_r`**, Hubble \(q_0=0.2\).

CSV: `field_count_verify.csv`, `count_funnel.csv`, `mag_cut_compare.csv`.

---

## 4. VizieR file `catalog/SDSS_DR16_cosmos.txt` — is it the same as Jimin?

**No. Same survey (DR16), different query, footprint, cuts, and axis-ratio column.**

### VizieR query (from file header)

Source: CDS VizieR `V/154/sdss16` (2023-11-20).

| Cut | VizieR COSMOS txt | Jimin SkyServer |
|-----|-------------------|-----------------|
| Geometry | Cone: (150.1255, 2.2108), **r = 1°** (~π deg²) | Box: RA 148–152, Dec 0–4 (**16 deg²**) |
| mode / clean | `mode==1`, `clean==1` | Same |
| Galaxy class | `class=3`, `rc==3`, `rs==0` | `type_r=3` |
| Mag | **`rPmag < 22`** (Petrosian) | **`p.r` ∈ [12,21]** (model) |
| Photo-z | `<zph> ≥ 0` | `nnAvgZ > 0` **and** `score > 0.8` |
| Star lnL | — | `lnLStar_r < −10` |
| Morph | — | V1: fracDeV=0 + lnL; V2: lnL |
| Axis ratio in file | **`rdVell` = deVAB_r** | **`expAB_r`** (selected; deV also fetched) |
| N | **19,739** | V1 9,120 / V2 17,657 (full box) |

Converted CSV: `catalog/SDSS_DR16_cosmos.csv`.

### Cross-match (1.5″)

| Check | Result |
|-------|--------|
| Jimin V1 objects inside 1° cone | 1,208 |
| Jimin V2 objects inside 1° cone | 2,299 |
| Those found in VizieR file | **100%** (both) |
| VizieR rows matching any V2 | 2,299 (11.6% of VizieR) |
| median \|rdVell − deVAB_jimin\| | **0.000224** (same column) |
| median \|rdVell − expAB\| | **0.058** (different column) |

So: Jimin spiral-selected objects in the COSMOS cone are a **subset** of the VizieR dump. The VizieR dump is **not** the Jimin selection — it is brighter/fainter-wider (petro&lt;22), has no lnL/fracDeV/score/lnLStar cuts, and ships **deV** ellipticity, not exp.

Density: 19,739 / π ≈ **6280 deg⁻²** to petroMag&lt;22 (deeper than mag=21), vs ~2960 deg⁻² to model r≤21 in the big box — consistent with going ~1 mag deeper.

---

## 5. Hubble CDFs from the VizieR file

Using `rdVell` (deVAB_r), ba &gt; 0.2, Hubble \(q_0=0.2\). Outputs: `vizier_cosmos/`.

| Sample | N (ba&gt;0.2) | median cos(i) |
|--------|------------:|--------------:|
| VizieR all (rPmag&lt;22), deVAB | 14,682 | 0.476 |
| VizieR rPmag≤21, deVAB | 8,369 | 0.540 |
| VizieR rmag≤21, deVAB | 8,828 | 0.533 |
| Jimin V1 in cone, **expAB** | 1,112 | 0.523 |
| Jimin V2 in cone, **expAB** | 2,146 | 0.541 |
| Matched ∩ V2: deVAB | 1,998 | 0.500 |
| Matched ∩ V2: expAB | 2,146 | 0.541 |

Plots:
- `vizier_cosmos/plots/cdf_*_deVAB.png`
- `vizier_cosmos/plots/cdf_overlay_vizier_mag_variants.png`
- `vizier_cosmos/plots/cdf_overlay_vizier_vs_jimin.png` (deVAB vs expAB — not apples-to-apples)
- `vizier_cosmos/plots/cdf_matched_deVAB_vs_expAB.png` (same galaxies; column difference only)

**Do not treat VizieR CDFs as a Jimin replication** unless you restrict to the matched subset and switch to expAB (requires joining SkyServer columns).

---

## 6. Bottom line

1. Jimin field counts to mag≈21 are literature-consistent (~3000 deg⁻²).
2. V1 N=9,120 is real after cuts, not a fetch bug.
3. `SDSS_DR16_cosmos.txt` is a **related but different** VizieR COSMOS cone dump (petro&lt;22, deVAB, no morph cuts).
4. All Jimin V1/V2 objects in that cone appear in the txt; most txt objects are **not** in Jimin V1/V2.
5. For inclination work matching Jimin, keep using **expAB_r** from the SkyServer catalogs under `catalog/v1_*.csv`, `catalog/v2_*.csv`.
