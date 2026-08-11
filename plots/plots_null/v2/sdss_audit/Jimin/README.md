# Jimin / advisor SDSS DR16 query replication

Reproduction of the historical SkyServer SQL selection used by the advisor + Jimin,
plus Hubble cos(i) CDFs. Release: **DR16**. Mag cut: **modelMag_r (p.r)** (`p.r BETWEEN 12 AND 21`).


**Mag cut confirmation:** SQL uses `p.r BETWEEN 12 AND 21`. In SDSS `PhotoObj`, shorthand `r` ≡ `modelMag_r` (exact match in DR16 probe: `max|r−modelMag_r|=0`). It is **not** `petroMag_r`. See sibling `petroMag/` for the Petrosian variant.

## Morphology cuts (two versions)

| Version | Morph SQL | Meaning |
|---------|-----------|---------|
| **V1** (historical) | `fracDeV_r = 0` **AND** `lnLDeV_r < lnLExp_r` | Pure cModel exponential **and** pure-exp likelihood wins |
| **V2** (weaker) | `lnLDeV_r < lnLExp_r` only | Pure-exp likelihood wins (SDSS `modelMag` rule) |

Shared base: RA 148–152, Dec 0–4, `p.r BETWEEN 12 AND 21`, `mode=1`, `clean=1`, `type_r=3`,
`lnLStar_r < -10`, Photoz `nnAvgZ > 0`, `score > 0.8`, no TOP (full COUNT).

After fetch: keep **`expAB_r` > 0.2**, Hubble cos(i) with `q0=0.2`
(`cosi_array_from_df(..., q_col='expAB_r')`). Not `deVAB_r`.

## Why V1 is ~9k (not a fetch bug)

`COUNT(*)` **without TOP** on DR16 for the exact V1 WHERE with **model** `p.r` returns
**9,120** — identical to the fetched catalog. Historical `TOP 10000` never limited V1.

Ablated funnel (`scripts/_jimin_count_audit.py`, DR16, no TOP, **model** mag):

| Stage | N |
|-------|--:|
| Box + mode/clean + `type_r=3` | 144,070 |
| + `r` in [12, 21] | 47,351 |
| + `lnLStar_r < -10` | 33,045 |
| + Photoz `nnAvgZ > 0` | 32,978 |
| + `score > 0.8` (base, no morph) | 28,043 |
| + lnLExp wins (**V2 full**) | **17,657** |
| + `fracDeV_r = 0` (**V1 full**) | **9,120** |

Petrosian swap (`petroMag_r` in [12, 21], same other cuts): V2 **17,138**, V1 **8,668**
— slightly *fewer*, not more. Overlap (V2): both mags in range 17,096; model-only 561;
petro-only 42.

## What `fracDeV` and `lnL*` mean (SDSS photometry)

Per [SDSS DR16 magnitudes](https://www.sdss4.org/dr16/algorithms/magnitudes/):

1. Fit pure exponential → `lnLExp_r`
2. Fit pure de Vaucouleurs → `lnLDeV_r`
3. Linear mix of those two (ellipses fixed) → `fracDeV_r` = deV weight in cModel

So `lnLDeV < lnLExp` only asks which *single* profile wins. `fracDeV = 0` requires the
**composite** to put zero bulge light — the stronger pure-disk cut.

## Results

| Version | N (SQL) | N (ba > 0.2) | median cos(i) | frac of V2 with fracDeV=0 |
|---------|--------:|----------------:|--------------:|--------------------------:|
| V1 fracDev0+lnL | 9,120 | 8,445 | 0.5458 | — |
| V2 lnL only | 17,657 | 16,551 | 0.5484 | 0.5165 |

TOP truncated: V1=False, V2=False (top_n=None).

## Files

- **`ANALYSIS.md`** — full field-count verification + VizieR COSMOS cross-check
- `FIELD_COUNTS.md`, `field_count_verify.csv`, `count_funnel.csv`
- `sql/v1_fracDev0_and_lnL.sql`, `sql/v2_lnL_exp.sql`
- `catalog/v1_*.csv`, `catalog/v2_*.csv`, `catalog/SDSS_DR16_cosmos.txt` (+ `.csv`)
- `plots/cdf_*.png`, `summary.csv`
- `vizier_cosmos/` — CDFs from VizieR file (deVAB) + crossmatch summary

```bash
python scripts/build_jimin_sdss_replication.py
python scripts/build_jimin_sdss_replication.py --mag-cut petro
python scripts/analyze_jimin_vizier_cosmos.py
```
