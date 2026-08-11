# Formal cos(i) audit — SDSS v2 null pool

All artifacts for this study live in `plots/plots_null/v2/sdss_audit/formal/` (decoupled from the parent v2 audit plots).

**Script:** `python scripts/audit_sdss_v2_cosi_formal.py`  
**Subsample stability:** `python scripts/audit_sdss_v2_cosi_subsample_stability.py`

---

## 1. Question

Why does a cumulative `modelMag_r` cut shift the null CDF toward face-on (higher cos i), when the full production pool before any mag limit has mean/median cos(i) ≈ 0.5?

We test three hypotheses:

1. **Code bug** — wrong axis ratio column or Hubble formula.
2. **Small-sample noise** — finite null draws mimic a face-on shift.
3. **Composition** — fainter galaxies in the pool have systematically lower `expAB_r` (more edge-on in projection); removing them via a mag cut changes the median cos(i).

---

## 2. Pool definition

| Stage | Cuts | N |
|-------|------|---|
| Raw v2 catalog | HTM-random `PhotoObj` galaxies | 1,900,000 |
| Production pool | u-r < 2.3 → lnLExp > lnLDeV → expAB_r > 0.2 | **434,962** |
| mag ≤ 21 (strict null) | + modelMag_r ≤ 21 | **159,144** |
| mag ≤ 22 | + modelMag_r ≤ 22 | 342,534 |

cos(i) from `expAB_r` via Hubble formula with q₀ = 0.2 (same as null CDF pipeline).

---

## 3. Results summary

### 3.1 Isotropy and mag-cut tests (`cosi_hypothesis_tests.csv`)

| Stage | N | Median cos(i) | 95% bootstrap CI | Permutation p |
|-------|---|---------------|------------------|---------------|
| Pre-mag (full production pool) | 434,962 | **0.495** | [0.494, 0.496] | — |
| mag ≤ 19 | 16,951 | 0.606 | [0.602, 0.613] | ≈ 0.001 |
| mag ≤ 20 | 55,679 | 0.613 | [0.610, 0.616] | ≈ 0.001 |
| mag ≤ 21 | 159,144 | **0.584** | [0.583, 0.586] | ≈ 0.001 |
| mag ≤ 22 | 342,534 | 0.523 | [0.523, 0.523] | ≈ 0.001 |

KS test on pre-mag cos(i) vs Uniform(0,1) rejects (p ≈ 0) because N is huge; the **effect size** is tiny (median 0.495 vs 0.5). Permutation: shuffle mag labels, re-apply mag cut; observed median at mag≤21 is far above the null — **mag and cos(i) are coupled in the data**.

### 3.2 Correlations (`cosi_mag_correlation.csv`)

| Quantity | Value |
|----------|-------|
| Spearman ρ(mag, cos i) | **−0.259** |
| Spearman ρ(expAB_r, cos i) | **1.000** (by construction: cos i = f(expAB_r)) |
| Partial ρ(mag, cos i \| expAB_r) | **+0.020** |

Fainter galaxies tend to lower cos(i). After conditioning on `expAB_r`, mag adds almost nothing — the mag-cut effect works **through b/a**, not an independent mag→orientation channel.

### 3.3 Per 0.5 mag bin (`cosi_per_mag_bin.csv`, `cosi_mag_joint_panel.png`)

| Mag bin | frac pool | Median cos(i) |
|---------|-----------|---------------|
| 18.0–19.0 | ~2.7% | **0.62** |
| 19.5–20.0 | ~5.7% | **0.61** |
| 20.0–20.5 | ~9.4% | **0.59** |
| 20.5–21.0 | ~14.4% | **0.56** |
| 21.0–21.5 | ~19.9% | **0.50** |
| 21.5–22.0 | ~22.3% | **0.45** |

Bright bins are face-on; the faint tail (mag > 21) is edge-on. A cumulative cut at mag 21 **drops ~42% of the pool** (everything with mag > 21) and raises the median cos(i).

### 3.4 Mixture and simulation nulls

**Mixture (`cosi_mixture_decomposition.csv`):**

| Scenario | Median cos(i) at mag≤21 |
|----------|-------------------------|
| Observed | 0.584 |
| Uniform reweight across mag bins | 0.593 |
| Weighted bin medians | 0.592 |

Reweighting does **not** remove the face-on shift — the bright bins themselves are face-on.

**Simulation (`cosi_simulation_null.csv`):**

| Null | median full | median mag≤21 |
|------|-------------|---------------|
| Isotropic cos(i) ~ U(0,1) | 0.498 | 0.498 |
| Shuffled cos(i) vs mag | 0.495 | 0.494 |
| **Observed** | **0.495** | **0.584** |

If orientation were isotropic and independent of mag, the mag cut would **not** move the median. Only the observed (mag, expAB_r) joint distribution produces the shift.

### 3.5 Subsample stability (`cosi_subsample_stability.csv`, `.png`)

Pool: strict production cuts + mag ≤ 21 (N = 159,144).  
Method: 500 random subsamples per size n (log-spaced 2,000 … 159,144); mean of subsample medians.

| n | Mean median cos(i) | 16–84% band |
|---|-------------------|-------------|
| 2,000 | 0.585 | 0.577 – 0.592 |
| 10,000 | 0.584 | 0.581 – 0.587 |
| 159,144 | 0.584 | (exact) |

Even n = 2,000 gives mean median ≈ **0.585**, not 0.5. The mag≤21 face-on value is a **property of the selected pool**, not subsample size. The band narrows with n (sampling noise only).

### 3.6 Sky magnitude distribution validation (`mag_counts_sky_validation.csv`)

**Why we compare N(m).** The cos(i)–mag coupling could be an artifact if our null pool’s magnitude histogram were unrepresentative of the on-sky galaxy population — e.g. a catalog build bug, a non-uniform HTM draw, or cuts that silently re-weight magnitudes in a way that mimics composition-driven cos(i) shifts. Before interpreting the mag-cut face-on shift, we check that the **empirical dN/dm in our sample matches what SDSS reports for the real sky** in the same magnitude system (`modelMag_r`).

Three normalized dN/dm curves (0.5 mag bins, 15.5–21.5):

| Reference | Role |
|-----------|------|
| **Full v2 catalog** (1.9M HTM galaxies, no morphology cuts) | Same pipeline and `modelMag_r` as the null pool; internal “what the sky looks like in this catalog.” |
| **Strict null pool** (434k after u-r / lnL / expAB cuts) | Production pool whose cos(i) we audit. |
| **Literature sky counts** (`reference/sdss_modelmag_r_counts.csv`) | [Yasuda et al. 2001](https://doi.org/10.1086/322093) commissioning counts re-mapped from Petrosian r* to **modelMag_r** using Δ(r*−modelMag_r)=+0.055 mag from [Stoughton et al. 2002](https://doi.org/10.1086/324741) / [Abazajian et al. 2009](https://doi.org/10.1088/0067-0049/182/2/543) SDSS photometry algorithms (PSF-convolved profile fits). No published all-sky modelMag count table exists; this is the closest literature match to our magnitude definition. |

| Mag range | pool / raw catalog | pool / modelMag sky ref |
|-----------|-------------------|-------------------------|
| 17–18 | 0.84 – 1.04 | 0.44 – 1.04 |
| 18–20.5 | **0.97 – 1.09** | **0.90 – 1.02** |
| 21.0–21.5 | ~0.95 | ~1.04 |

Over 18–20.5 the strict-pool shape tracks both the **full catalog** and the **modelMag-adjusted SDSS sky counts** within ~10%. Global χ² tests fail (morphology-selected exp sample vs all galaxies; pool extends past Yasuda’s faint limit), but the **shape agreement** supports the conclusion that our mag distribution is not pathological — the cos(i) shift at mag≤21 is driven by (mag, expAB_r) composition, not a bogus N(m).

Petrosian-only Yasuda counts remain in `reference/yasuda2001_r_counts.csv` for historical comparison.

### 3.7 Median expAB_r per mag bin (`ba_per_mag_bin*.csv`, `ba_mag_joint_panel*.png`)

Median **raw PhotoObj `expAB_r`** per 0.5 mag bin — **no Hubble mapping**. Three cut variants:

| Variant | Files | Cuts |
|---------|-------|------|
| Full catalog | `ba_per_mag_bin.csv`, `ba_mag_joint_panel.png` | None (finite b/a only) |
| Morphology | `ba_per_mag_bin_ur_lnl.csv`, `ba_mag_joint_panel_ur_lnl.png` | u-r + lnL exp-wins; all b/a |
| Strict | `ba_per_mag_bin_strict.csv`, `ba_mag_joint_panel_strict.png` | + `expAB_r > 0.2` |

**Strict pool vs cos(i):** `median(cos i) = Hubble(median expAB_r)` per bin (monotonic map; cos i ~0.02–0.06 lower). Overlay: `ba_cosi_strict_overlay.png`.

**Docs:** [`EXPAB_R_BA_PLOTS.md`](EXPAB_R_BA_PLOTS.md), [`THREE_PLOT_DEEP_ANALYSIS.md`](THREE_PLOT_DEEP_ANALYSIS.md) §2.

### 3.8 cos(i) with fixed Hubble q₀ = 0.2 and edge-on clip (`cosi_per_mag_bin_fixed_q0_clip.csv`, `cosi_mag_bin_fixed_q0_clip.png`)

Four curves vary **`q_min`**; Hubble **q₀ = 0.2** fixed; `q ≤ 0.2 → cos(i) = 0`. Median cos(i) hits **0** when >50% of a bin is clipped — **not** because `q > 0.2` galaxies are absent. See [`THREE_PLOT_DEEP_ANALYSIS.md`](THREE_PLOT_DEEP_ANALYSIS.md) §3.

### 3.9 cos(i) with joint strict q_min = q₀ (`cosi_per_mag_bin_joint_strict.csv`, `cosi_mag_bin_joint_strict.png`)

`q_min = q₀` varied jointly; no clip pileup. Production line matches §3.3. Faint bins still contain galaxies with `q > 0.2` (e.g. N ≈ 1,800 at mag 23–23.5). See [`THREE_PLOT_DEEP_ANALYSIS.md`](THREE_PLOT_DEEP_ANALYSIS.md) §4.

### 3.10 Cut survival vs mag (`cut_survival_vs_mag.csv`, `.png`)

From the 1.9M catalog, fraction passing each cut vs mag bin:

- **u-r:** ~33% at mag 18 → ~53% at mag 21.5 (faint galaxies more likely to pass u-r — less quenched).
- **lnL exp-wins:** ~25% at mag 18 → ~33% at mag 21.5.
- **expAB_r > 0.2:** combined survival falls at mag > 21.5 (~22% at 21.5–22).

Morphology cuts reshape N(m) but do not obviously invert the cos(i)–mag trend; the joint (mag, expAB_r) pattern in the production pool drives the CDF shift.

---

## 4. Analysis and interpretation

### 4.1 No evidence for a pipeline bug

Before any mag cut, median cos(i) = **0.495** with a tight bootstrap CI. Cut order (strict b/a vs u-r vs lnL) does not matter (see parent `cosi_mag_bias_audit.csv`). cos(i) is computed from `expAB_r` after lnL exp-wins — consistent with production null CDFs.

**Conclusion:** the Hubble mapping and column choice are not the source of the face-on shift at mag ≤ 21.

### 4.2 The mag cut effect is real and statistically formalized

Permutation p ≈ 0.001 for every cumulative limit tested. Δmedian(mag≤21 vs full) = **+0.089**. The effect is:

- **Largest** at mag ≤ 20 (median 0.613) because the faintest edge-on tail is still excluded.
- **Moderate** at mag ≤ 21 (0.584) — production null cut.
- **Small** at mag ≤ 22 (0.523) — most of the faint edge-on population is back in.

This matches the per-bin table: excluding mag > 21 removes bins with median cos(i) ≈ 0.45–0.50.

### 4.3 Mechanism: mag–b/a composition, not orientation physics

Three lines agree:

1. **Partial Spearman** — mag does not predict cos(i) once expAB_r is fixed.
2. **Simulation null** — isotropic cos(i) stays at 0.50 under mag cut.
3. **Per-bin medians** — monotonic drop in cos(i) from mag ~19 to mag ~22.

Interpretation: in this SDSS `expAB_r` sample, **fainter galaxies measure smaller b/a** (more edge-on in projection). Possible contributors (not mutually exclusive):

- SDSS profile fitting noise at low S/N (faint → noisier axis ratio).
- Real dust/inclination selection ([Shao et al. 2007](https://doi.org/10.1086/511131)).
- Mix of morphologies at the faint end (exp-winner fraction vs mag).

This is **not** the same as claiming disks are physically edge-on preferentially at fainter magnitudes in an isotropic universe; it is the **measured** (mag, expAB_r) joint in our catalog.

### 4.4 Subsample stability rules out small-N artifacts

For null CDF work at mag ≤ 21 (N ≈ 159k), even random subsamples of a few thousand galaxies reproduce median cos(i) ≈ 0.58. FRB comparisons or MC draws with similar pool size will inherit this **composition**, not an artifact of drawing too few galaxies.

### 4.5 Implications for null CDF interpretation

| Observation | Implication |
|-------------|-------------|
| Full pool cos(i) ≈ 0.5 | Null axis-ratio pipeline is consistent with isotropic projection **before** mag cut. |
| mag ≤ 21 median ≈ 0.58 | CDFs at mag 20/21/22 will sit **face-on of uniform** if drawn from this pool — **expected** given composition, not necessarily a bug. |
| mag ≤ 22 median ≈ 0.52 | Loosening the cut toward the full pool median. |
| dN/dm ~ sky over 18–20.5 | Mag distribution matches full catalog and modelMag-adjusted SDSS counts — **not** a driver of the cos(i) skew. |

When comparing FRB hosts to nulls, match **the same mag cut and cuts** and interpret CDF offset in light of the **0.58 reference median**, not 0.5.

### 4.6 Caveats

- Literature sky counts use Yasuda et al. 2001 Petrosian r* re-mapped to modelMag_r; offset is typ. 0.05–0.08 mag for mixed morphologies, 0.055 mag adopted here ([Stoughton et al. 2002](https://doi.org/10.1086/324741); [Abazajian et al. 2009](https://doi.org/10.1088/0067-0049/182/2/543)).
- Strict pool is exp-selected; bright-end ratios vs all-galaxy counts differ by design.
- Very faint bins (mag > 23) have small N and noisy cos(i).
- Spearman(expAB, cos i) = 1 is tautological; use partial correlation and simulations for mechanism.
- KS p ≈ 0 on isotropy is a large-N effect; rely on median and simulation nulls for practical conclusions.
- **Median summary statistic:** For strict pools, `median(cos i) = Hubble(median expAB_r)` per mag bin (monotonic map). Use **median**, not mean, when comparing mag-bin tracks to uniform cos(i) = 0.5; see §3.7 and `ba_cosi_strict_overlay.png`.

---

## 5. File inventory

| File | Description |
|------|-------------|
| `cosi_hypothesis_tests.csv` | Isotropy, bootstrap CIs, permutation tests |
| `cosi_mag_correlation.csv` | Spearman, partial Spearman, OLS |
| `cosi_per_mag_bin.csv` | Median cos(i) and pool fraction per 0.5 mag |
| `cosi_mixture_decomposition.csv` | Reweighting scenarios |
| `cosi_simulation_null.csv` | Isotropic and shuffled nulls |
| `cosi_subsample_stability.csv` / `.png` | Bootstrap median vs subsample size |
| `mag_counts_sky_validation.csv` / `mag_counts_comparison.png` | N(m) sky validation (pool vs catalog vs literature) |
| `cut_survival_vs_mag.csv` / `.png` | Cut efficiency vs mag |
| `cosi_mag_joint_panel.png` | cos(i) and N(m) vs mag |
| `ba_per_mag_bin.csv` / `ba_mag_joint_panel.png` | Median raw expAB_r (full v2 catalog) |
| `ba_per_mag_bin_ur_lnl.csv` / `ba_mag_joint_panel_ur_lnl.png` | + u-r, lnL exp-wins; all b/a |
| `ba_per_mag_bin_strict.csv` / `ba_mag_joint_panel_strict.png` | + expAB_r > 0.2 |
| `ba_cosi_strict_overlay.png` / `ba_cosi_strict_comparison.csv` | Strict pool: median b/a vs median cos(i) |
| `EXPAB_R_BA_PLOTS.md` | SDSS expAB_r provenance, cut definitions, median monotonicity |
| `cosi_per_mag_bin_fixed_q0_clip.csv` / `cosi_mag_bin_fixed_q0_clip.png` | Vary q_min, fixed Hubble q₀=0.2, edge-on clip |
| `cosi_per_mag_bin_joint_strict.csv` / `cosi_mag_bin_joint_strict.png` | Joint q_min = q₀ strict overlay |
| `THREE_PLOT_DEEP_ANALYSIS.md` | Detailed interpretation of plots 1–3 (floor, clip, selection) |
| `reference/sdss_modelmag_r_counts.csv` | modelMag_r sky counts (Yasuda + Stoughton offset) |
| `reference/yasuda2001_r_counts.csv` | Petrosian r* counts (Yasuda 2001 Table 2) |

---

## 6. References

- [Yasuda et al. 2001, AJ 122, 1104](https://doi.org/10.1086/322093) — SDSS commissioning galaxy number counts (photometric catalog verification; Petrosian r*)
- [Stoughton et al. 2002, AJ 123, 485](https://doi.org/10.1086/324741) — SDSS Early Data Release; defines PSF-convolved **modelMag** photometry
- [Abazajian et al. 2009, ApJS 182, 543](https://doi.org/10.1088/0067-0049/182/2/543) — SDSS DR7; Petrosian vs model/cmodel magnitude offsets
- [Smith et al. 2002, AJ 123, 1381](https://doi.org/10.1086/420800) — SDSS spatial homogeneity of galaxy counts (extinction / calibration check)
- [Shao et al. 2007, ApJ 654, 898](https://doi.org/10.1086/511131) — inclination-dependent LF / dust in SDSS spirals
- [Driver et al. 2007, MNRAS 379, 1022](https://doi.org/10.1111/j.1365-2966.2007.11865.x) — LF vs inclination
- [Blanton et al. 2003, AJ 125, 2348](https://doi.org/10.1086/320405) — SDSS luminosity function
- Parent audit: `../COSI_MAG_BIAS_AUDIT.md`
