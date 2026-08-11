# Production SDSS v2 vs Jimin criteria

Live DR16 counts in RA [148,152] × Dec [0,4] (16 deg²).  
Script: `scripts/_jimin_vs_prod_funnel.py` → `prod_vs_jimin_funnel.csv`.

**Note:** The production HTM-random CSV (`SDSS_catalog_v2_fullsky_modelmr.csv`) has **0** rows in this box (sparse full-sky sample). Field contrasts below use complete SkyServer `COUNT(*)`, not that CSV.

---

## Side-by-side criteria

| Cut | Production SDSS v2 | Jimin (advisor SQL) |
|-----|--------------------|---------------------|
| Footprint | Full SDSS sky (HTM random) | Fixed 4°×4° box |
| `mode=1`, `clean=1` | Yes (at fetch) | Yes |
| Galaxy type | **`type = 3`** | **`type_r = 3`** |
| Mag | Post-hoc **`modelMag_r ≤ limit`** (usually 21) | SQL **`p.r` ∈ [12, 21]** (= modelMag) |
| Star reject | Type only | **`lnLStar_r < −10`** |
| Photo-z | No | **INNER JOIN Photoz**, `nnAvgZ > 0` |
| Photo-z quality | No | **`score > 0.8`** |
| Disk morph | **`lnLExp_r > lnLDeV_r`** (CDF) | Same for V2; V1 also **`fracDeV_r = 0`** |
| Color / late-type | **`u − r < 2.3`** (CDF) | **No** |
| Axis ratio | **`expAB_r > 0.2`**, Hubble | Same |
| `fracDeV` | Stored; **not** used as hard cut | V1: exact **0** |

Shared scientifically: clean primary galaxies, model *r*, exponential-likelihood winner, `expAB_r`, Hubble \(q_0=0.2\).

---

## What Jimin adds that we don’t

### 1. `lnLStar_r < −10`
Rejects objects whose light is better fit by a PSF (star-like).

| Pool | N in box |
|------|--------:|
| Base (type_r=3, mag 12–21) | 47,351 |
| After lnLStar | 33,045 |
| **Drop** | **−14,306 (−30%)** |

On production-style base + lnL only: 26,130 → 20,732 (**−21%**).

**Value:** Real at the faint end — `type=3` leaks compact/star contaminants. Cost: can also remove compact galaxies. For an inclination null, moderately useful purity cut; not morphology of disks.

### 2. Photoz `nnAvgZ > 0`
| After lnLStar | 33,045 |
| After nnAvgZ | 32,978 |
| **Drop** | **−67 (−0.2%)** |

**Value:** Negligible. Almost every object already has a Photoz row. Does not improve shape selection.

### 3. `score > 0.8`
| After Photoz | 32,978 |
| After score | 28,043 |
| **Drop** | **−4,935 (−15%)** |

On prod+lnL: 26,130 → 21,963 with Photoz+score (**−16%**).

**Value:** Keeps higher-confidence photometric redshifts. For a **pure geometric** inclination null (no distance weighting), this is weakly motivated — you pay 15% for redshift-quality that cos(i) does not use. Useful if the null is meant to match a photo-z-selected host population.

### 4. `fracDeV_r = 0` (Jimin V1 only)
| Jimin V2 (lnL) | 17,657 |
| Jimin V1 (+fracDeV=0) | 9,120 |
| **Drop** | **−8,537 (−48% of V2)** |

On prod+lnL: 26,130 → 13,373 (**−49%**).

**Value:** Strongest “pure exponential cModel” cut. Scientifically meaningful if you want zero bulge light in the composite fit — stricter than lnL alone (which only asks which *single* profile wins). Cost: halves the disk pool; many exp-winning galaxies still have `fracDeV > 0`. For FRB-host-like late-type disks, **lnL (+ color) is usually enough**; fracDeV=0 is optional purity, not required for Hubble-from-expAB.

### 5. `type_r = 3` vs `type = 3`
47,504 → 47,351 (**−0.3%**). Effectively irrelevant.

---

## What we add that Jimin doesn’t

### `u − r < 2.3` (Strateva-style late-type)
| Jimin V2 | 17,657 |
| + u−r &lt; 2.3 | 12,177 |
| **Drop** | **−5,480 (−31%)** |

Production path in box: lnL 26,130 → +u−r 17,352 (**−34%** of lnL pool).

**Value:** High for disk-inclination science. Separates blue/late-type from red/early-type. Jimin’s morph cuts do not replace this — red exp-winners remain without a color cut.

---

## Mag 12–21 field: how hard Jimin cuts relative to a shared base

Shared starting point ≈ production fetch + Jimin mag window:

**`mode=1`, `clean=1`, `type=3`, `modelMag_r` ∈ [12, 21]** → **N = 47,504**

| Stage | N | Remaining | Cut this step |
|-------|--:|----------:|--------------:|
| Prod-style base | 47,504 | 100% | — |
| → `type_r=3` | 47,351 | 99.7% | −0.3% |
| → + `lnLStar < −10` | 33,045 | 69.6% | −30% |
| → + Photoz `nnAvgZ>0` | 32,978 | 69.4% | −0.2% |
| → + `score>0.8` | 28,043 | 59.0% | −15% |
| → + lnLExp wins (**Jimin V2**) | **17,657** | **37.2%** | −37% |
| → + `fracDeV=0` (**Jimin V1**) | **9,120** | **19.2%** | −48% of V2 |

So relative to “all clean galaxies in the field at 12≤r≤21”:

- Jimin **V2** keeps **~37%**
- Jimin **V1** keeps **~19%**

Production CDF path on the same base (no Jimin extras):

| Stage | N | Remaining |
|-------|--:|----------:|
| Base | 47,504 | 100% |
| + lnLExp wins | 26,130 | 55.0% |
| + u−r &lt; 2.3 | 17,352 | 36.5% |
| + expAB &gt; 0.2 | **16,406** | **34.5%** |

Strict end states in this box (ba &gt; 0.2):

| Selection | N |
|-----------|--:|
| Production (lnL + u−r + ba) | **16,406** |
| Jimin V2 (Jimin extras + lnL + ba) | **16,551** |
| Jimin V1 (+fracDeV=0 + ba) | **8,445** |

N(V2) ≈ N(prod) here, but the **galaxies differ**: Jimin drops stars/photo-z-failures and keeps some red exp-winners; production drops red galaxies and keeps objects failing score/lnLStar.

---

## How valuable are Jimin’s extras? (summary)

| Addition | Fraction removed (typical) | Value for inclination null |
|----------|---------------------------:|----------------------------|
| `type_r` vs `type` | ~0% | None |
| `lnLStar < −10` | ~20–30% | Moderate (star purity) |
| `nnAvgZ > 0` | ~0% | None |
| `score > 0.8` | ~15% | Low–moderate (only if photo-z selection matters) |
| `fracDeV = 0` (V1) | ~50% of lnL pool | High purity, **harsh**; optional |
| *(ours)* `u−r < 2.3` | ~30% of lnL pool | **High** for late-type disks |

**Practical takeaway:** Jimin’s distinctive stack is mostly **star-lnL + photo-z score + (optional) fracDeV=0**. The morph cut we already share (lnL) is the main disk selector. For FRB-host-like nulls, production’s **u−r** is more on-target than Jimin’s **score**; **fracDeV=0** is the only Jimin-only cut that strongly changes the physical sample (and the CDF size).
