# Deep analysis: expAB_r, fixed-q₀ clip, and joint-strict diagnostic plots

**Date:** 2026-06-17  
**Catalog:** `catalog/SDSS_catalog_v2_fullsky_modelmr.csv` (N = 1.9M HTM-sampled PhotoObj rows)  
**Plots analysed (cos i q₀ sensitivity):**

| # | PNG | CSV |
|---|-----|-----|
| A | `cosi_mag_bin_fixed_q0_clip.png` | `cosi_per_mag_bin_fixed_q0_clip.csv` |
| B | `cosi_mag_bin_joint_strict.png` | `cosi_per_mag_bin_joint_strict.csv` |

**Related expAB_r plots (no Hubble):** see [`EXPAB_R_BA_PLOTS.md`](EXPAB_R_BA_PLOTS.md) — `ba_mag_joint_panel*.png` (full / ur+lnL / strict).

Reference production cos(i) curve: `cosi_per_mag_bin.csv` / `cosi_mag_joint_panel.png`.

---

## Executive summary

**The faint-end behaviour is real in the catalog columns, but it is not “all galaxies have b/a ≈ 0”.**

1. **Raw expAB_r (full catalog):** At mag 22.0–22.5 the **median** `expAB_r` ≈ **0.10** (~30% at SDSS floor 0.05). See [`EXPAB_R_BA_PLOTS.md`](EXPAB_R_BA_PLOTS.md).

2. **Plot A (fixed q₀ = 0.2 clip):** At mag ≳ 22.5 the **median cos(i) can hit exactly 0** when **more than half** of the pool receives `cos(i) = 0` via the rule `q ≤ 0.2 → edge-on`. This does **not** mean there are no galaxies with `q > 0.2`; it means the **median** is dominated by the clip pileup. At mag 22.5–23.0 in the `q_min = 0.05` pool, **~68%** have `q ≤ 0.2` (hence `cos(i) = 0`) while **~32%** still have `q > 0.2`.

3. **Plot B (joint strict):** The production line (`q_min = q₀ = 0.2`) **never** assigns `cos(i) = 0` from clipping — every galaxy has `q > 0.2` by construction. Faint bins still show falling median cos(i) (e.g. **0.35** at mag 22.5–23.0, **0.29** at 23.0–23.5) because the **measured** `expAB_r` values in the surviving sample are modest (median **0.40** and **0.35**), not because the sample is empty.

**Bottom line:** The strange faint-end shapes are a combination of (a) SDSS `expAB_r` quantization/floor at 0.05, (b) low-S/N profile fitting pushing faint galaxies toward small axis ratios, and (c) for plot 2 only, an intentional **edge-on assignment** for `q ≤ 0.2`. They are **not** evidence that the pipeline is secretly applying the Hubble formula to plot 1.

---

## 1. What each plot is doing

### Plot 1 — `ba_mag_joint_panel.png`

| Item | Definition |
|------|------------|
| Sample | **Entire v2 CSV** (1.9M rows) |
| Cuts | None (only finite `modelMag_r`, `expAB_r ∈ [0, 1]`) |
| y-axis | **Median of raw `expAB_r`** per 0.5 mag bin |
| Transform | **None** — no Hubble, no q₀, no u-r, no lnL |

### Plot 2 — `cosi_mag_bin_fixed_q0_clip.png`

| Item | Definition |
|------|------------|
| Sample | Production pre-pool per line: **u-r**, **lnL exp-wins**, **`expAB_r > q_min`** |
| `q_min` | Varies: **0.05, 0.10, 0.15, 0.20** (pool rebuilt each line) |
| Hubble q₀ | **Fixed at 0.20** for all lines |
| Edge rule | If `q ≤ 0.2`, set **`cos(i) = 0`** (via `hubble_cosi_from_ba`) |
| y-axis | Median cos(i) per mag bin |

### Plot 3 — `cosi_mag_bin_joint_strict.png`

| Item | Definition |
|------|------------|
| Sample | Same cuts as plot 2 |
| `q_min` | Varies: **0.05, 0.10, 0.15, 0.20** |
| Hubble q₀ | **Equals `q_min`** each line |
| Edge rule | **None** — pool already requires `q > q₀` |
| y-axis | Median cos(i) per mag bin |

At **`q_min = q₀ = 0.20`**, plots 2 and 3 **coincide** and match `cosi_per_mag_bin.csv`.

---

## 2. Plot 1 — raw `expAB_r` vs magnitude

### 2.1 What the curve shows

| Mag bin | N (full catalog) | Median `expAB_r` | Mean `expAB_r` | frac `= 0.05` | frac `≤ 0.2` | frac `> 0.2` |
|---------|------------------|------------------|----------------|---------------|--------------|--------------|
| 18.0–18.5 | 17,452 | 0.625 | 0.585 | 0.021 | 0.077 | 0.923 |
| 20.0–20.5 | 133,794 | 0.581 | 0.563 | 0.021 | 0.077 | 0.923 |
| 21.0–21.5 | 322,661 | 0.427 | 0.432 | 0.087 | 0.220 | 0.780 |
| 21.5–22.0 | 440,619 | 0.278 | 0.323 | 0.201 | 0.411 | 0.589 |
| **22.0–22.5** | **416,618** | **0.101** | **0.226** | **0.292** | **0.619** | **0.381** |
| 22.5–23.0 | 154,137 | 0.095 | 0.174 | 0.230 | 0.751 | 0.249 |
| 23.0–23.5 | 11,972 | 0.097 | 0.164 | 0.212 | 0.769 | 0.231 |

**Global catalog:** median `expAB_r` = **0.345**; **15.4%** exactly **0.05**; **36.2%** ≤ 0.2.

### 2.2 Answering: “Are all b/a values basically zero at mag 22–22.5?”

**No.** Three separate points:

1. **The median is ~0.10, not ~0.** Most galaxies report small but non-zero `expAB_r`.

2. **~30% are exactly 0.05** — the documented SDSS PhotoObj storage floor (`SDSS_BA_FLOOR_MIN = 0.05` in `null_catalog_utils.py`). In mag 22.0–22.5, **124,348 / 416,618** rows have `expAB_r = 0.05` exactly. This is a **catalog quantization limit**, not a physical “zero axis ratio”.

3. **~38% still have `expAB_r > 0.2`** even in the raw bin. Bright-face-on and edge-on galaxies coexist; the **median** is pulled down by the floor pileup and the `≤ 0.2` tail.

### 2.3 Why mean ≫ median in faint bins

At mag 22.0–22.5: mean = **0.226** vs median = **0.101**.

The distribution is **bimodal / heavy-tailed**: a large floor-dominated component near 0.05–0.15 and a substantial tail with `expAB_r` ~ 0.3–0.7. The mean is sensitive to the upper tail; the median tracks the dominant low-`q` pileup.

Top reported values in mag 22.0–22.5 (by count):

| `expAB_r` | Count |
|-----------|-------|
| **0.050** | **124,348** |
| 0.100 | 4,972 |
| 0.099 | 3,852 |
| 0.150 | 3,348 |
| 0.350 | 2,827 |

### 2.4 Physical / pipeline interpretation

Fainter galaxies in SDSS DR16 PhotoObj:

- Have **lower S/N** → exponential profile fits are less stable → **`expAB_r` biased small**.
- Hit the **0.05 storage floor** when the fitted axis ratio would be smaller (values are truncated/quantized in the database).
- Include more **true edge-on** and **bulge-dominated** systems at fixed apparent magnitude.

Plot 1 therefore shows **what SDSS reports**, including measurement limits. It does **not** by itself prove a physical edge-on excess at fixed luminosity — but it **does** show that any inclination analysis using raw faint `expAB_r` inherits these systematics.

### 2.5 Noisy bins at mag ≳ 24

Bins above mag ~23 have **N < 2,000** in the full catalog (e.g. 994 at 23.5–24.0). Medians jump (e.g. 0.072 at 23.5–24, then 0.286 at 26.0–26.5 with N = 71) — **small-number noise**, not a physical turnover.

---

## 3. Plot 2 — fixed Hubble q₀ = 0.2 with edge-on clip

### 3.1 Mechanism

For each galaxy:

```
if expAB_r ≤ 0.2:
    cos(i) = 0
else:
    cos(i) = sqrt((q² − 0.2²) / (1 − 0.2²))
```

Galaxies with **`q_min < q ≤ 0.2`** are **kept in the pool** (if they pass u-r and lnL) but forced **edge-on**. Lower `q_min` → more galaxies in the pileup band → **lower median cos(i)**.

### 3.2 Key per-bin numbers (`q_min = 0.05` line)

| Mag bin | Pool N | Median cos(i) | Mean cos(i) | Interpretation |
|---------|--------|---------------|-------------|----------------|
| 21.5–22.0 | 122,871 | 0.364 | 0.366 | ~majority with cos(i) > 0 |
| 22.0–22.5 | 121,352 | 0.166 | 0.246 | ~43% clipped to 0 |
| **22.5–23.0** | **64,404** | **0.000** | **0.127** | **>50% clipped → median = 0** |
| **23.0–23.5** | **6,772** | **0.000** | **0.097** | same |
| 23.5–24.0 | 352 | 0.000 | 0.177 | tiny N |

**Mag 22.0–22.5 detail (`q_min = 0.05` pool):**

| Quantity | Value |
|----------|-------|
| N | 121,352 |
| Median `expAB_r` | 0.258 |
| Fraction with `q ≤ 0.2` | **42.7%** |
| Fraction with `cos(i) = 0` | **42.7%** |
| Median cos(i) | **0.166** |

**Mag 22.5–23.0 detail (`q_min = 0.05` pool):**

| Quantity | Value |
|----------|-------|
| N | 64,404 |
| Median `expAB_r` | 0.117 |
| Fraction with `q ≤ 0.2` | **67.9%** |
| Fraction with `cos(i) = 0` | **67.9%** |
| Median cos(i) | **0.000** (exactly) |

### 3.3 Answering: “At high mag, is there no q > 0.2?”

**No — that is not what the clip plot shows.**

- In the **raw** catalog at mag 22.5–23.0, **24.9%** of galaxies have `expAB_r > 0.2` (38,353 / 154,137).
- In the **clip pool** at the same bin, **~32%** still have `q > 0.2` (100% − 67.9%).
- The **median cos(i) = 0** occurs because **>50%** are assigned **exactly 0** by the clip rule, not because `q > 0.2` galaxies are absent.

The **mean** cos(i) stays **> 0** (0.127 at 22.5–23.0) — the face-on tail remains in the mean.

### 3.4 Effect of varying `q_min` (at mag 22.0–22.5)

| `q_min` | Median cos(i) | Why |
|---------|---------------|-----|
| 0.05 | 0.166 | Largest pileup band (0.05, 0.2] |
| 0.10 | 0.317 | Fewer clipped |
| 0.15 | 0.365 | Fewer clipped |
| 0.20 | 0.410 | No clip zone (all q > 0.2) — matches plot 3 |

This is the intended **edge-on pileup** diagnostic: holding Holmberg q₀ = 0.2 fixed while admitting thinner disks.

### 3.5 Spurious spikes at mag ≳ 25.5

Lines with `q_min ≤ 0.15` show median cos(i) → **~0.99** at mag 25.5–26.0 with **N ≈ 40–65** per bin. These are **not robust** — the production pool has only **43–46** galaxies there. Do not over-interpret.

---

## 4. Plot 3 — joint strict (`q_min = q₀`)

### 4.1 Mechanism

Each line uses **`expAB_r > q_min`** and Hubble with the **same** `q₀ = q_min`. No sub-threshold pileup: every galaxy in the pool has `q > q₀`, so cos(i) is always the positive Hubble solution.

### 4.2 Production line (`q_min = q₀ = 0.20`) — matches `cosi_per_mag_bin.csv`

| Mag bin | N | Median cos(i) | Median `expAB_r` |
|---------|---|---------------|------------------|
| 21.0–21.5 | 86,360 | 0.504 | 0.533 |
| 21.5–22.0 | 97,030 | 0.450 | 0.484 |
| 22.0–22.5 | 69,564 | 0.410 | 0.449 |
| 22.5–23.0 | 20,694 | 0.354 | 0.400 |
| 23.0–23.5 | 1,833 | 0.294 | 0.350 |
| 23.5–24.0 | 129 | 0.469 | 0.501 |

**Every galaxy has `expAB_r > 0.2` by construction** (verified: frac `> 0.2` = 1.000 at all bins).

### 4.3 Answering: “Does the strict plot lose all q > 0.2 at high mag?”

**No.** At mag 23.0–23.5 the strict production pool still has **N = 1,833** with median `expAB_r` = **0.35** and median cos(i) = **0.29**. The sample shrinks because **most** faint raw galaxies fail `expAB_r > 0.2`, but the survivors are **not** all edge-on.

**Selection effect:** At mag 22.0–22.5, only **38%** of raw catalog galaxies have `q > 0.2`; after u-r + lnL + strict cut, **69,564** remain — the strict pool is **biased toward the high-`q` tail** of the faint population.

### 4.4 Comparing `q_min` lines at faint bins

Lower `q_min` (thinner intrinsic disk in Hubble):

- **Admits more galaxies** (larger N).
- Maps each `q` to a **higher** cos(i) for the same axis ratio (smaller q₀).
- At mag 22.0–22.5: median cos(i) ranges from **0.253** (`q_min = 0.05`) to **0.410** (`q_min = 0.20`).

Unlike plot 2, there is **no cos(i) = 0 floor** from clipping — faint medians stay **> 0** until bins become tiny.

### 4.5 Agreement with plot 2 at `q_min = q₀ = 0.20`

Per-bin medians are **identical** (max difference 0.0). When the pool cut equals Hubble q₀, clip and strict are the same operation.

---

## 5. Cross-plot comparison

### 5.1 Same faint bin, three views (mag 22.0–22.5)

| View | Median y | What drives it |
|------|----------|----------------|
| **Plot 1** raw `expAB_r` | **0.101** | SDSS floor + faint measurement bias; **no** inclination math |
| **Plot 2** clip, `q_min=0.05` | **0.166** cos(i) | 43% forced to cos(i)=0 at q₀=0.2 |
| **Plot 2** clip, `q_min=0.20` | **0.410** cos(i) | Same as strict production |
| **Plot 3** strict, `q_min=0.05` | **0.253** cos(i) | Low q₀ Hubble; no zero pileup |
| **Plot 3** strict, `q_min=0.20` | **0.410** cos(i) | Production null mapping |

### 5.2 Why plot 1 “looks like” plot 3 but is not the same

Both fall at faint magnitudes, but:

- Plot 1 tracks **reported axis ratio** (bounded below by 0.05).
- Plot 3 tracks **Hubble cos(i)** on a **selected** subsample with `q > q₀`.

A median `expAB_r` ≈ 0.10 does **not** equal median cos(i) ≈ 0.10. For q₀ = 0.2, `expAB_r` = 0.10 implies cos(i) = **0** (below q₀). The production strict median cos(i) at the same raw population bin is **0.41** because galaxies with `q ≤ 0.2` are **removed**, not clipped.

### 5.3 Role of morphology cuts (mag 22.0–22.5)

| Stage | N | Median `expAB_r` | frac `≤ 0.2` |
|-------|---|------------------|--------------|
| Raw catalog | 416,618 | 0.101 | 0.619 |
| After u-r | 274,731 | 0.108 | 0.609 |
| After lnL exp-wins | 166,397 | 0.128 | 0.582 |
| Strict `q > 0.2` (production) | 69,564 | 0.449 | 0.000 |

lnL and u-r change the medians **slightly** at this bin; the **strict q cut** is the largest change (removes ~58% of lnL survivors).

---

## 6. Implications for the null CDF audit

| Finding | Implication |
|---------|-------------|
| Raw faint `expAB_r` piles up at **0.05** | Any analysis on unfiltered v2 `expAB_r` is dominated by catalog systematics below ~0.15. |
| **Production strict pool** removes `q ≤ 0.2` | Median cos(i) at mag ≤ 21 ≈ **0.58** is a property of **selected** discs, not raw PhotoObj. |
| **Clip mode** (plot 2) | Shows sensitivity to treating thin disks as edge-on at fixed q₀ = 0.2; medians can hit **0** when clip fraction > 50%. |
| **Joint strict** (plot 3) | Shows sensitivity to **both** sample definition and q₀; no artificial cos(i)=0 floor. |
| Mag ≳ 24 bins | Too few galaxies for stable medians — treat as indicative only. |

### Recommended reading order

1. **Plot 1** — is the mag trend in the **catalog column**?
2. **Plot 3** (`q_min = 0.2`) — production null behaviour (`cosi_per_mag_bin.csv`).
3. **Plot 2 vs 3** — how much does **clip vs exclude** matter for q₀ = 0.2 sensitivity?

---

## 7. Reproduction

```bash
# Regenerate all formal audit outputs including the three plots
python scripts/audit_sdss_v2_cosi_formal.py

# Quick per-bin floor statistics (example)
python -c "
import sys, numpy as np, pandas as pd
sys.path.insert(0, 'scripts')
from null_catalog_utils import read_sdss_null_catalog
from pipeline_null_plot_utils import DEFAULT_SDSS_V2
df = read_sdss_null_catalog(DEFAULT_SDSS_V2)
mag = pd.to_numeric(df['modelMag_r'], errors='coerce')
ba = pd.to_numeric(df['expAB_r'], errors='coerce')
m = (mag>22)&(mag<=22.5)&ba.notna()
b = ba[m].to_numpy()
print('N', len(b), 'med', np.median(b), 'frac_0.05', (b==0.05).mean())
"
```

---

## 8. References in this repo

- `NULL_CATALOG_DATA_AND_INCLINATION_AUDIT.md` §5.1 — global floor and `expAB_r ≤ 0.2` fractions (v1; same PhotoObj fields).
- `scripts/null_catalog_utils.py` — `SDSS_BA_FLOOR_MIN = 0.05`, `hubble_cosi_from_ba()`, pool builders.
- `FORMAL_COSI_AUDIT.md` §3.3, §3.7–3.9 — plot inventory and high-level interpretation.
