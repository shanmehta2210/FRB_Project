# SDSS v2: cos(i) vs magnitude — full pipeline audit

**Question:** Why does a brighter `modelMag_r` cut shift the null CDF toward face-on (higher cos i), when orientations are expected to be isotropic (mean cos i ≈ 0.5)?

**Short answer:** No implementation bug was found. The **full production pool before any magnitude limit** has mean cos(i) ≈ **0.501** and median ≈ **0.495**. The shift appears **only when a cumulative magnitude cut removes the faint tail**, where SDSS `expAB_r` is systematically **lower** (more edge-on in projection). That is a **composition effect** in the data (and/or SDSS shape measurement), not a wrong formula or wrong axis-ratio column.

---

## 1. Pipeline trace (v2)

| Step | Script / function | What it does | Inclination risk? |
|------|-------------------|--------------|-------------------|
| SQL fetch | `build_sdss_null_catalog_v2.py` | HTM-stratified random `PhotoObj`; `type=3`, `clean=1`, `mode=1`; **no Dec clip, no mag cut** | None — random on sky, not on orientation |
| Columns stored | `build_catalog()` | `modelMag_r`, `expAB_r`, `lnLExp_r`, `lnLDeV_r`, …; `best_model_ba_r` computed but **not used** in CDF pools | OK if CDF uses `expAB_r` only |
| lnL patch | (in SQL) | 100% lnL coverage in v2 build | Rows without lnL dropped at lnL cut only |
| u-r cut | `filter_sdss_ur` | `u-r < 2.3` | No direct i dependence |
| Morphology | `filter_sdss_drop_dev_winners` | `lnLExp_r > lnLDeV_r`; then **`expAB_r` only** | Correct exp profile |
| Shape cut | `apply_strict_q_cut` | **`expAB_r > 0.2`** (not `best_model_ba_r`) | Truncates very edge-on; symmetric in i |
| Magnitude | `slice_null_base_by_mag` | **`modelMag_r <= limit`** | **Changes galaxy mix** (see below) |
| cos(i) | `hubble_cosi_from_ba(expAB_r, q0=0.2)` | Standard Hubble formula | Correct direction: high b/a → high cos i |
| CDF plot | `plot_sdss_null_cut_evolution_cdfs.py` | Full-pool ECDF of cos(i); MC overlay optional | Plotting does not bias the pool |

Cut order (production `prepare_null_strict_color_base`): u-r → lnL → strict b/a.  
Cut order (evolution diagnostic): strict b/a → u-r → lnL.  
**Audit:** both orders give **identical** N and mean cos(i) at mag 21/22.

---

## 2. Smoking gun: mean cos(i) by pipeline stage

From `sdss_cut_evolution/stage_summary.csv` (v2):

| Stage | N | mean cos(i) | median cos(i) | median expAB_r |
|-------|---|-------------|---------------|----------------|
| After strict b/a only | 1,212,324 | 0.497 | 0.489 | 0.519 |
| + u-r | 658,630 | 0.487 | 0.477 | 0.508 |
| + lnL exp-wins | 434,962 | **0.501** | **0.495** | 0.525 |
| **+ modelMag_r ≤ 21** | 159,144 | **0.567** | **0.584** | 0.606 |

**Before the magnitude cut, mean cos(i) is already ≈ 0.5.**  
The face-on shift is introduced entirely by **`modelMag_r <= 21`**, not by lnL, colour, or expAB handling.

---

## 3. Per magnitude bin (production pool, no cumulative cut)

From `cosi_mag_bias_audit.csv`:

| modelMag_r bin | N | mean cos(i) | median expAB_r |
|----------------|---|-------------|----------------|
| 21.0 – 21.5 | 86,360 | **0.505** | 0.533 |
| 21.5 – 22.0 | 97,030 | 0.463 | 0.484 |
| 22.0 – 22.5 | 69,564 | 0.429 | 0.449 |
| 20.0 – 20.5 | 40,823 | 0.573 | 0.613 |
| 19.5 – 20.0 | 24,714 | 0.590 | 0.633 |

**Fainter bins are more edge-on in projection** (lower `expAB_r` → lower cos i).  
**Brighter bins sit near or above cos i ≈ 0.5–0.6.**

A cumulative cut `mag <= 21` **drops the faint edge-on bins** (21.5+) and keeps the brighter face-on-rich bins → mean cos(i) rises.

| Cumulative cut | mean cos(i) |
|----------------|-------------|
| mag ≤ 23 | **0.501** |
| mag ≤ 22 | 0.522 |
| mag ≤ 21 | 0.567 |
| mag ≤ 20 | 0.589 |

Including more faint galaxies **restores** mean cos(i) toward 0.5.

---

## 4. Code checks performed

- **`expAB_r` used everywhere** in production pools and CDFs; `best_model_ba_r @ mag21` gives different statistics (not used in production).
- **Cut order** A vs B: identical results.
- **Mag column:** `modelMag_r` only (not `cmodelMag_r` / `rmag`).
- **Hubble formula:** `cos²(i) = (b/a² − q₀²)/(1 − q₀²)`; values clipped to [0, 1].
- **HTM sampling:** no orientation or magnitude pre-selection in SQL.
- **Plotting:** ECDF is a monotone sort of pool cos(i); MC envelope only affects overlay plots, not pool statistics.

---

## 5. Interpretation (for discussion with supervisor)

1. **Isotropic orientations do not imply mean cos(i) = 0.5 after observing b/a cuts.**  
   Even for truly random orientations, selecting on **projected** `expAB_r > 0.2` skews cos(i). Here the full pool after colour + lnL still lands near 0.5, so that is not the dominant effect.

2. **The mag-dependent shift is a sample composition effect:**  
   In this SDSS exp-winning late-type sample, **`expAB_r` increases with brightness** (median b/a rises from ~0.45 at 22–22.5 to ~0.61 at 20–20.5). A brighter magnitude limit removes faint, edge-on-rich bins.

3. **Possible astrophysical / instrumental causes (not code bugs):**
   - SDSS `expAB_r` measurement bias vs SNR (edge-on disks harder to fit at low flux).
   - Mixture of true disk thickness + inclination with a flux-limited sample.
   - lnL exp-winner may correlate with light concentration at fixed apparent mag (less tested here; effect is small before mag cut).

4. **Reference line on CDF plots:** The dashed “Uniform” diagonal (F = cos i) assumes **uniform cos(i) in the final sample**. That is **violated** by any cut that changes the mag–b/a mix, even if intrinsic orientations are isotropic.

---

## 6. Regenerate audit

```bash
python scripts/audit_sdss_v2_cosi_mag_bias.py
```

Outputs: `plots/plots_null/v2/sdss_audit/cosi_mag_bias_audit.csv`

---

## 7. Formal diagnostics (extended)

The parent audit above uses **production** pools. For deeper mag-bin work see **`formal/`**:

| Document | Content |
|----------|---------|
| [`formal/FORMAL_COSI_AUDIT.md`](formal/FORMAL_COSI_AUDIT.md) | Hypothesis tests, permutation, sky N(m) |
| [`formal/EXPAB_R_BA_PLOTS.md`](formal/EXPAB_R_BA_PLOTS.md) | Raw PhotoObj `expAB_r` (3 cut levels); SDSS provenance; median b/a vs median cos(i) |
| [`formal/THREE_PLOT_DEEP_ANALYSIS.md`](formal/THREE_PLOT_DEEP_ANALYSIS.md) | q_min vs q₀ clip vs joint-strict cos(i) overlays; SDSS 0.05 floor |

```bash
python scripts/audit_sdss_v2_cosi_formal.py
```

**Key findings (formal):**

- Median **raw** `expAB_r` on the full 1.9M catalog falls to ~0.10 at mag 22–22.5 (SDSS **0.05 storage floor** + faint S/N), not because all b/a ≈ 0.
- Median **cos(i)** on the strict pool has the **same bin-to-bin shape** as median **expAB_r** because Hubble(q₀=0.2) is **monotonic** → `median(cos i) = Hubble(median b/a)` per bin (values differ by ~0.02–0.06). See `formal/ba_cosi_strict_overlay.png`.
- **Clip** vs **joint strict** q₀ tests: [`formal/cosi_mag_bin_fixed_q0_clip.png`](formal/cosi_mag_bin_fixed_q0_clip.png), [`formal/cosi_mag_bin_joint_strict.png`](formal/cosi_mag_bin_joint_strict.png).
