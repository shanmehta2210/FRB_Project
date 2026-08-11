# Do Jimin's `lnLStar` / `score` cuts bias cos(i)? — and why is production more face-on?

Empirical test on the **same field** (RA[148,152]×Dec[0,4], model r ≤ 21), one DR16 fetch of all 47,351 galaxies (`type_r=3`, mode, clean) with every diagnostic column, then sliced in pandas. Cos(i): Hubble, `expAB_r`, q0=0.2, ba>0.2.

Scripts: `scripts/analyze_jimin_cut_bias.py`. Data: `cut_bias/`.

---

## TL;DR

1. **`type=3` does NOT guarantee "not a star."** It is a coarse concentration cut. `lnLStar` is stricter and still removes ~20%.
2. **Jimin's cuts do NOT skew face-on.** `lnLStar` mildly skews **edge-on**; `score`, `nnAvgZ`, `fracDeV=0` are shape-neutral. All Jimin cuts together move median cos(i) by only **−0.007**.
3. **The reason production looks more face-on is OUR `u−r < 2.3` color cut** (+0.025 in median cos(i)), which removes dust-reddened **edge-on** disks. This is a real inclination-dependent selection effect in *our* pipeline, not Jimin's.

---

## 1. Does `type=3 + mode + clean` already exclude stars?

No — only coarsely. Per [SDSS classification](https://www.sdss.org/dr17/algorithms/classify/), `type` is set by a single concentration test:

```
extended (type=3) if  psfMag − cmodelMag > 0.145
```

- It is a **binary morphology proxy**, valid at ~95% to r≈21, and known to misclassify compact galaxies and bright-nucleus Seyferts.
- `lnLStar_r` is the **actual PSF-fit log-likelihood** (probability the object is consistent with a point source). `lnLStar_r < −10` demands the object be *strongly* inconsistent with a PSF — a much harder "extended" requirement than `type=3`.

**They are not equivalent:** among `type_r=3` galaxies (lnL-exp base), **4,780 / 24,222 (20%)** still fail `lnLStar_r < −10`. So Jimin's `lnLStar` is a genuine extra purity cut, mostly catching **compact, marginally-resolved** objects that `type=3` let through.

---

## 2. Direction of bias: cos(i) of KEPT vs REMOVED (base = lnL-exp winners)

| Cut | N kept | median cos(i) kept | N removed | median cos(i) removed | shift (removed−kept) |
|-----|-------:|-------------------:|----------:|----------------------:|---------------------:|
| `lnLStar_r < −10` | 19,442 | 0.5491 | 4,780 | **0.5639** | **+0.015** |
| `score > 0.8` | 20,406 | 0.5519 | 3,816 | 0.5553 | +0.003 |
| `nnAvgZ > 0` | 24,142 | 0.5523 | 80 | 0.5630 | +0.011 (n=80) |
| `fracDeV_r = 0` | 12,142 | 0.5532 | 12,080 | 0.5514 | −0.002 |
| **`u−r < 2.3` (ours)** | 16,347 | **0.5773** | 7,875 | **0.5029** | **−0.074** |

Reading:
- **`lnLStar`** removes objects that are *slightly more face-on* (removed median 0.564 > kept 0.549). So it pulls the surviving sample **toward edge-on**, not face-on. Mechanism: the roundest galaxies look most PSF-like, so they fail `lnLStar` most (roundest b/a bin fails at **32%** vs ~20% typical). But the effect is small.
- **`score`, `nnAvgZ`, `fracDeV=0`** are essentially **shape-neutral** (|shift| ≤ 0.003, and fracDeV=0 halves N with no cos(i) change).
- **Our `u−r < 2.3`** is the big one: removed (red) objects are **edge-on** (0.503), kept (blue) are **face-on** (0.577). It strongly biases the sample face-on.

---

## 3. Why our color cut biases face-on: dust

Edge-on disks have longer dust sightlines → redder `u−r`. In this field:

| b/a bin | median u−r | frac u−r > 2.3 |
|---------|-----------:|---------------:|
| 0.2–0.3 (edge-on) | 2.14 | 42% |
| 0.5–0.6 | 1.91 | 32% |
| 0.9–1.0 (face-on) | 1.82 | 27% |

Spearman(b/a, u−r) = **−0.13** (p≈1e−92): edge-on ⇒ redder. Red (u−r>2.3) median b/a = 0.50; blue median b/a = 0.58. So cutting `u−r > 2.3` preferentially deletes **edge-on** galaxies → face-on-biased CDF. This is the classic inclination–extinction selection effect.

---

## 4. Same-field CDFs: production vs Jimin decomposed

All in the identical field, model r ≤ 21, expAB_r, ba>0.2:

| Selection | N | median cos(i) |
|-----------|--:|--------------:|
| Pool (gal + mag + lnL) | 24,222 | 0.5523 |
| **Production** (+ `u−r < 2.3`) | 16,347 | **0.5773** |
| **Jimin V2** (+ lnLStar + photoz + score) | 16,551 | **0.5484** |
| Jimin V1 (+ fracDeV=0) | 8,445 | 0.5458 |

**Gap production − Jimin V2 = +0.029**, decomposed from the common pool (0.5523):

- Production: `u−r` cut → **+0.025** (face-on)
- Jimin V2: lnLStar+score+photoz → **−0.004** (edge-on)
- Sum: 0.025 − (−0.004) = **0.029** ✓ (fully accounted)

Incremental Jimin stack (each step's median shift):

| Stage | N | median | Δ |
|-------|--:|-------:|--:|
| gal+mag+lnL | 24,222 | 0.5523 | — |
| +lnLStar<−10 | 19,442 | 0.5491 | −0.003 |
| +nnAvgZ>0 | 19,423 | 0.5491 | 0.000 |
| +score>0.8 (V2) | 16,551 | 0.5484 | −0.001 |
| +fracDeV=0 (V1) | 8,445 | 0.5458 | −0.003 |

Total Jimin effect on median: **−0.0065**. Negligible and toward edge-on.

Plot: `cut_bias/cdf_prod_vs_jimin_samefield.png`.

---

## 5. Answers to the questions

- **Do `score`/`lnLStar` skew our morphology toward face-on?** No. `score` is shape-neutral; `lnLStar` skews *slightly edge-on* (removes ~20% of galaxies, biased to the round/compact ones, shifting median cos(i) by −0.003). Neither creates a face-on excess.
- **Doesn't `type=3` already remove stars?** Only coarsely (`psfMag−cmodelMag>0.145`). `lnLStar<−10` is a stronger PSF-inconsistency cut and legitimately removes ~20% more compact interlopers. Reasonable purity cut, minor shape cost.
- **Why is production (mag<21) more face-on than Jimin?** Because of **our `u−r < 2.3` color cut**, which removes dust-reddened edge-on disks (+0.025 median cos(i)). Jimin has no color cut, so it keeps those edge-ons and sits lower. The difference is **our** selection effect, not a Jimin bias.

**Implication:** if we want the least inclination-biased null, the color cut is the item to scrutinize (it is inherently inclination-dependent via dust), far more than Jimin's `score`/`lnLStar`. A dust-aware or color-free late-type selection would reduce our face-on bias.
