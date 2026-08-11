# Null catalog data provenance and cos(i) CDF audit

**Date:** 2026-05-28  
**Purpose:** Document **exactly** how SDSS/Legacy null catalogs are built, how cos(i) is derived, and why the inclination CDF **does not** cross 0.5 at cos(i)=0.5 even when the lnL exponential-winner cut is applied correctly.

**Quick conclusion:** The dashed “Uniform” line (CDF = cos i) is the reference for **isotropic random viewing angles**. This pipeline **does not** draw random orientations. It uses **SDSS PhotoObj projected axis ratios** (`expAB_r`), applies **magnitude, colour, profile, and b/a cuts**, then maps b/a → cos(i) with the Hubble formula. Bright, late-type, exponential-winning galaxies in SDSS are **face-on biased in projected shape**. That explains median cos(i) ≈ **0.61** at mag&lt;20 with strict cuts—not a failure of the lnL filter.

---

## 1. Reference line on the plots

| Element | Definition |
|--------|------------|
| **X-axis** | cos(i) from Hubble inversion of catalog b/a (not measured inclination). |
| **Y-axis** | Empirical CDF of that cos(i) sample. |
| **Black dashed “Uniform”** | F(c) = c — i.e. **isotropic orientations** with P(cos i) uniform on [0, 1]. |
| **Green curve** | SDSS null pool after all cuts; MC band = 10⁴ bootstrap draws of **N = N_FRB** galaxies without replacement from that pool (`cdf_envelope` in `scripts/pipeline_null_plot_utils.py`, seed 42). |

**Important:** Matching the uniform diagonal requires the **final sample** to have uniformly distributed cos(i). Our sample is filtered on **observed b/a**, magnitude, and colour. That generically breaks uniformity even if lnL exp-wins are perfect.

---

## 2. Data download — exact sources

### 2.1 SDSS DR16 (`catalog/SDSS_catalog_v1_allsky_modelmr.csv`)

| Item | Value |
|------|--------|
| **Service** | SDSS SkyServer SQL via `astroquery.sdss.SDSS.query_sql` |
| **Table** | `PhotoObj` (alias `p`) |
| **Build script** | `scripts/build_sdss_null_catalog.py` |
| **Patch lnL** | `scripts/patch_sdss_profile_winner.py` + `scripts/merge_lnl_patch_into_sdss.py` |
| **Cache** | `catalog/SDSS_lnl_patch_cache.csv` (optional re-merge) |

**Per-RA chunk SQL** (12 bins, `TOP` 50k–80k per bin, `type=3` GALAXY, `clean=1`):

```sql
SELECT TOP {N}
    p.ra, p.dec,
    p.cmodelMag_r, p.petroMag_r,
    p.modelMag_r, p.modelMag_u, p.modelMag_g,
    p.deVMag_r, p.expMag_r,
    p.lnLDeV_r, p.lnLExp_r,
    p.deVAB_r, p.expAB_r,
    p.fracDeV_r, p.expRad_r, p.deVRad_r,
    p.type
FROM PhotoObj AS p
WHERE p.ra >= {ra_min} AND p.ra < {ra_max}
  AND p.dec >= -30 AND p.dec <= 90
  AND p.type = 3
  AND p.clean = 1
  AND p.cmodelMag_r > 0 AND p.cmodelMag_r < 90
  AND p.deVAB_r > 0 AND p.deVAB_r <= 1
  AND p.expAB_r > 0 AND p.expAB_r <= 1
```

**Post-query (build script):**

- Deduplicate `(ra, dec, cmodelMag_r)`; target ~500k rows (shuffle seed 42).
- `rmag` ← `cmodelMag_r` (composite; **not** used for CDF mag cut).
- `u_r` ← `modelMag_u - modelMag_r`, `g_r` ← `modelMag_g - modelMag_r`.
- CDF magnitude: **`modelMag_r`**; CDF axis ratio: **`expAB_r`** only after morphology.

**Sampling bias (SDSS):** `TOP N` per RA slice is **not** a volume-limited sample. It favours **bright** objects in each slice. Combined with `dec ∈ [-30°, 90°]`, the parent sample is **not** all-sky uniform in depth or orientation.

**Current file on disk:** ~457k rows (some rows lost during an earlier bad merge); ~**87%** have finite `lnLExp_r`/`lnLDeV_r` after patch.

**Profile winner (production):**

```text
keep  iff  lnLExp_r > lnLDeV_r   (strict; ties → drop)
       and both lnL finite
```

Implemented in `sdss_exp_wins_lnl_mask()` / `filter_sdss_drop_dev_winners()` in `scripts/null_catalog_utils.py`. CDF pools: **0 violations** checked by `scripts/audit_sdss_exp_winner_pool.py`.

---

### 2.2 Legacy Survey DR10 (`catalog/LS_catalog_v1_allsky_modelmr.csv`)

| Item | Value |
|------|--------|
| **Service** | NOIRLab TAP `https://datalab.noirlab.edu/tap` |
| **Table** | `ls_dr10.tractor` |
| **Build script** | `scripts/build_legacy_catalog_csv.py` |

**TAP query (default `TOP 500000`, `brick_primary=1`, joint Dec):**

```sql
SELECT TOP {N}
    objid, ra, dec, type, brick_primary,
    flux_g, flux_r, flux_i, flux_z,
    sersic, shape_r,
    shape_e1, shape_e2,
    shape_e1_ivar, shape_e2_ivar
FROM ls_dr10.tractor
WHERE brick_primary = 1
  AND type <> 'PSF'
  AND flux_r > 0
  AND shape_e1/2 and ivars finite, |e| < 1
  AND dec >= -30 AND dec <= 90   -- if --region joint
ORDER BY RANDOM()   -- or ORDER BY objid + shuffle
```

**Derived:**

- `tractor_mag_r` = 22.5 − 2.5 log₁₀(flux_r) (= `rmag` alias).
- `expAB_r` = (1 − |e|) / (1 + |e|) from `shape_e1`, `shape_e2`.
- `gmag` from flux_g; `g-r` = gmag − rmag.

**Legacy CDF cuts:** `g-r < 0.75`; exclude `REX`, `DEV`; keep `EXP` or `0.75 ≤ rdVrad ≤ 2`; strict `expAB_r > 0.2`. **No lnL** (Tractor only).

---

### 2.3 FRB hosts (`pipeline_galfit_results.csv`)

| Item | Value |
|------|--------|
| **Inclination** | GALFIT `inc` (degrees) → cos(i) for CDF |
| **b/a** | GALFIT `b_a`; strict pool: **`b_a > 0.2`** |
| **Magnitude** | GALFIT `mag`; cut `mag ≤` mag limit (20, 21, …) |
| **Colour** | **No** host colour cut on FRB side |

---

## 3. cos(i) mapping (Hubble thin disk)

**Code:** `hubble_cosi_from_ba()` in `scripts/null_catalog_utils.py`, `Q0 = 0.2`.

\[
\cos^2 i = \frac{(b/a)^2 - q_0^2}{1 - q_0^2}, \quad q_0 = 0.2
\]

- If \((b/a)^2 < q_0^2\) or invalid → **cos(i) = 0** (i = 90°).
- If \((b/a) > 1\) → cos(i) = 1.

**This is an invertible map from projected b/a to inclination** assuming fixed intrinsic thickness—not a draw from an isotropic orientation model.

---

## 4. CDF pipeline (plotting)

**Driver:** `scripts/plot_null_mag_cut_cdfs.py`

**SDSS null pool order** (`prepare_null_strict_color_base`):

1. Load CSV (`read_sdss_null_catalog`, minimal `usecols` incl. lnL).
2. `u-r < 2.3` (`SDSS_UR_MAX_CDF`).
3. **`lnLExp_r > lnLDeV_r`** (finite lnL only).
4. **`expAB_r > 0.2`** (strict).
5. Slice `modelMag_r < mag_limit` (e.g. 20).

**Plot:**

- `cosi = cosi_array_from_df(pool, q_col="expAB_r")`.
- `cdf_envelope(cosi, n_sample=N_FRB, n_draws=10000)` → green curve + 68% band.
- FRB: `frb_hosts_for_cdf` → GALFIT inc → cos(i).

---

## 5. Stage-by-stage decomposition (SDSS, mag &lt; 20)

Reproduced with `scripts/decompose_cosi_cdf_bias.py` → `cosi_cdf_audit/stage_summary_mag20.csv`.

| Stage | N | Median cos(i) | CDF(cos i = 0.5) | Median expAB_r |
|-------|---|---------------|------------------|----------------|
| Full catalog (no cuts) | 457,190 | **0.26** | 0.70 | 0.32 |
| `modelMag_r < 20` | 49,157 | **0.55** | 0.45 | 0.57 |
| + `u-r < 2.3` | 14,498 | 0.57 | 0.42 | 0.59 |
| + lnL exp-wins | 10,490 | 0.58 | 0.40 | 0.60 |
| + strict `b/a > 0.2` | **9,581** | **0.61** | **0.35** | 0.63 |

**Monte Carlo checks:**

| Model | Median cos(i) | CDF(0.5) |
|--------|---------------|----------|
| Isotropic cos(i) ~ U(0,1) | 0.50 | 0.50 |
| Isotropic cos(i), then strict b/a&gt;0.2 via Hubble | 0.50 | 0.50 |

So **strict b/a alone does not** move the median off 0.5 for a true isotropic population. The shift comes from **which galaxies enter the sample** (especially **magnitude** and **SDSS expAB_r distribution**).

### 5.1 Why the full catalog looks “edge-on”

- ~**17%** of rows have `expAB_r == 0.05` (PhotoObj storage floor).
- ~**38%** have `expAB_r ≤ 0.2` → Hubble maps many to **cos(i) ≈ 0**.
- Median cos(i) over **all** galaxies ≈ 0.26 (edge-on pile-up in projected shape).

**v2 mag-bin detail:** [`plots/plots_null/v2/sdss_audit/formal/THREE_PLOT_DEEP_ANALYSIS.md`](plots/plots_null/v2/sdss_audit/formal/THREE_PLOT_DEEP_ANALYSIS.md) §2; raw `expAB_r` plots in [`formal/EXPAB_R_BA_PLOTS.md`](plots/plots_null/v2/sdss_audit/formal/EXPAB_R_BA_PLOTS.md).

### 5.2 Why mag &lt; 20 looks “face-on”

- Bright subsample median **expAB_r ≈ 0.57** → median cos(i) ≈ **0.55** (close to uniform **before** strict b/a).
- Adding **strict b/a &gt; 0.2** removes the remaining low-b/a tail → median → **~0.61** (matches your plot).

**The lnL exp-winner step** changes median cos(i) only slightly (0.577 → 0.610 after strict at mag 20 path); it is **not** the main cause of the offset from the diagonal.

---

## 6. Interpreting your mag&lt;20 plot (N=9581, median cos i ≈ 0.65 on CDF)

The plot shows **CDF(cos i) ≈ 0.35 at cos i = 0.5** (35% below 0.5) → median cos i ≈ **0.61–0.65**. That matches the table above.

**This is not evidence that deV profiles remain in the pool** (audit: **0** lnL violations in CDF sample).

It **is** evidence that:

1. **Projected-shape selection** (bright + late-type + high b/a) favours face-on systems.
2. The **uniform line is the wrong null** for “what we should see after our cuts” unless we explicitly want isotropic **orientations** independent of b/a.
3. FRB hosts track the same biased null (red on green) — comparison is internally consistent.

---

## 7. “All-sky” and sampling checks

| Check | Finding |
|-------|---------|
| **Dec footprint** | Both catalogs: Dec ∈ [-30°, 90°] (Legacy∩SDSS overlap), not full sphere. |
| **RA coverage** | 0–360° in bins; not uniform mass/luminosity per bin (`TOP` per slice). |
| **Magnitude** | CDF uses `modelMag_r < limit` — surface-brightness / size bias toward face-on discs. |
| **Colour** | `u-r < 2.3` removes red sequence; changes mix, not orientation draw. |
| **lnL patch** | ~13% rows lack lnL → excluded from CDF (not silently kept as deV). |
| **Hubble q₀** | Fixed 0.2; FRB GALFIT uses same convention for strict b/a. |

There is **no** step that draws isotropic cos(i); therefore **no step** should be expected to produce a uniform cos(i) CDF after b/a-based selection.

---

## 8. What would be “surprising” vs “expected”

| Claim | Verdict |
|-------|--------|
| “lnL cut not applied” | **False** — filter is strict; pool audit passes. |
| “CDF should hit 0.5 at cos i = 0.5 for our null” | **Only if** orientations are isotropic **before** b/a selection. They are not. |
| “Median cos i ≈ 0.65 is impossible” | **Possible** when selecting high-b/a, bright, exponential-winning discs in SDSS. |
| “Plot unchanged from pre-lnL era” | Pool **N** dropped (~133k → ~10k at mag&lt;21); shape similar because **mag + b/a** dominate. |

---

## 9. Scripts for reproduction

```bash
# Verify every CDF galaxy has lnLExp > lnLDeV
python scripts/audit_sdss_exp_winner_pool.py --mag-limit 20

# Stage-by-stage median cos(i) table
python scripts/decompose_cosi_cdf_bias.py

# Regenerate mag-cut CDFs
python scripts/plot_null_mag_cut_cdfs.py --mag-limits 20 21 22
```

Outputs: `cosi_cdf_audit/stage_summary_mag20.csv`, `cosi_cdf_audit/cdf_curves_mag20.csv`.

### 9.1 v2 formal mag-bin diagnostics

Full write-up: [`plots/plots_null/v2/sdss_audit/formal/FORMAL_COSI_AUDIT.md`](plots/plots_null/v2/sdss_audit/formal/FORMAL_COSI_AUDIT.md).

```bash
python scripts/audit_sdss_v2_cosi_formal.py
```

| Output family | Purpose |
|---------------|---------|
| `ba_mag_joint_panel*.png` | Median **raw** `expAB_r` (no Hubble); full / ur+lnL / strict cuts |
| `ba_cosi_strict_overlay.png` | Strict pool: median b/a vs median cos(i) |
| `cosi_mag_bin_fixed_q0_clip.png` | Vary `q_min`, fixed q₀=0.2, edge-on clip |
| `cosi_mag_bin_joint_strict.png` | `q_min = q₀` varied jointly |

Supporting docs: `formal/EXPAB_R_BA_PLOTS.md`, `formal/THREE_PLOT_DEEP_ANALYSIS.md`.

---

## 10. Recommendations (if you want CDF nearer “uniform”)

1. **Change the reference curve** — e.g. bootstrap from a model with isotropic cos(i) → implied b/a → apply same cuts (see MC row in §5).
2. **Do not cut on observed b/a before CDF** — use inclusive cos(i) for shape comparison only (mode B; not used in current production CDFs because of i=90° pile-up).
3. **Volume-limited sample** — replace `TOP` SQL with a proper flux-limited selection.
4. **Restore full 500k SDSS catalog** — `python scripts/build_sdss_null_catalog.py` then re-patch lnL for complete coverage.

---

## 11. File / column contract (CDF)

| Survey | Mag cut | Shape for cos(i) | Morphology |
|--------|---------|------------------|------------|
| SDSS | `modelMag_r` | **`expAB_r`** only | `lnLExp_r > lnLDeV_r`, `u-r < 2.3` |
| Legacy | `tractor_mag_r` / `rmag` | `expAB_r` | no REX/DEV; EXP or n∈[0.75,2]; `g-r < 0.75` |
| FRB | GALFIT `mag` | GALFIT `inc` | `b_a > 0.2` only |

---

*This document supersedes informal assumptions that the black “Uniform” diagonal is the expected outcome for the morphology-cut null. It is the expected outcome for **isotropic cos(i)** only.*

---

## Appendix A — v1 cut order, funnel sizes, column quick reference

Merged from the former `plots/plots_null/v1_null_cdf_inclination/diagnostics/NULL_CATALOG_AND_CDF_METHOD.md` (now a stub). Machine-readable funnel: `plots/plots_null/v1_null_cdf_inclination/diagnostics/cut_funnel.csv` (`python scripts/audit_and_plot_null_v1_diagnostics.py --full`). Pre-morphology-cut plots: `Archive/plots_null_pre_morphology_cut/`.

### A.1 Catalog files (v1 on disk)

| File | Build script | Rows |
|------|----------------|------|
| `catalog/SDSS_catalog_v1_allsky_modelmr.csv` | `scripts/build_sdss_null_catalog.py` | ~501,885 |
| `catalog/LS_catalog_v1_allsky_modelmr.csv` | `scripts/build_legacy_catalog_csv.py` | 500,000 |

Footprint: RA 0–360°, Dec **−30° to +90°**. Cuts applied at plot time in `scripts/null_catalog_utils.py`.

### A.2 Cut order (CDF pools)

1. Magnitude (`mag <= limit`)
2. Colour — SDSS \(u-r < 2.3\); Legacy \(g-r < 0.75\)
3. SDSS: drop deV winners (`lnLExp_r > lnLDeV_r`); shape **`expAB_r` only**
4. Legacy: drop `REX`, `DEV`; keep `EXP` or \(0.75 \le n \le 2\); shape `expAB_r`
5. Strict: \(b/a > 0.2\)

Hubble \(q_0 = 0.2\). FRB hosts: `pipeline_galfit_results.csv`; mag + GALFIT \(b/a > 0.2\); no host colour cut.

### A.3 Funnel — m &lt; 21 (post morphology)

| Survey | After mag | After color | After morphology | After \(b/a>0.2\) | **Final** |
|--------|-----------|-------------|------------------|-------------------|-----------|
| Legacy | 43,113 | 5,980 (g−r) | 4,577 (EXP∪n) | 4,374 | **4,374** |
| SDSS | 156,343 | 29,844 (u−r + lnL exp) | — | 27,474 | **27,474** |

SDSS lnL patch matches ~77% of catalog rows; unmatched rows drop at exp-winner. Legacy CDF pool at m&lt;21 is &lt; 5,000 after morphology.

### A.4 Column quick reference

| Column | Survey | CDF use |
|--------|--------|---------|
| `modelMag_r` | SDSS | Mag cut |
| `expAB_r` | SDSS / Legacy | Shape / Hubble \(i\) |
| `lnLExp_r`, `lnLDeV_r`, `model_winner_is_exp` | SDSS | Profile winner |
| `tractor_mag_r` | Legacy | Mag cut |
| `rdVrad`, `tractor_type` | Legacy | Sérsic \(n\) / EXP·REX·DEV |

**Changelog (v1 method):** 2026-05-27 late-type morphology cuts; 2026-07-21 folded into this canonical audit.
