# COSMOS HST vs SDSS b/a audit — catalog & methodology decisions

**Date:** 2026-07-02  
**Science goal:** Replicate the SDSS v2 null audit style (median `expAB_r` vs `modelMag_r`, mag-bin composition, optional Re slices) **separately** for HST and SDSS over the COSMOS field. If HST (PSF-safe) shows a **flatter** mag–b/a trend than SDSS (PSF-limited), attribute the SDSS-only effect to **PSF and/or pipeline**; if trends match, PSF is not the dominant explanation.

**Status:** **Complete (Zurich HST + SDSS no-colour).** Farmer download deferred; see [ZURICH_CATALOG_NOTES.md](ZURICH_CATALOG_NOTES.md).

---

## 0. Production mode for this run (2026-07)

**HST source:** COSMOS Zurich morphology v1.0 (`cosmos_morph_zurich_1.0.tbl`), GIM2D parameters.  
**Colour:** **None on both surveys** (ablation vs v2 full-sky strict null with u−r < 2.3).

The symmetric u−r policy in §0 below remains documented for v2 parity but is **not** used in the current Zurich production run.

---

## 0a. Symmetric colour-cut policy (v2 parity; not current production)

**Decision:** **Both surveys get the same colour cut: `u − r < 2.3`.**

| Survey | Colour definition | Threshold |
|--------|-------------------|-----------|
| **SDSS DR17** | `u_r = modelMag_u − modelMag_r` | `< 2.3` |
| **HST (Farmer)** | `u_r = CFHT_u_MAG − HSC_r_MAG` | `< 2.3` |

**Why apply colour to both (not neither):**
- v2 SDSS null audit uses `u-r < 2.3` (`SDSS_UR_MAX_CDF`) to select late-type / spiral-like galaxies.
- Dropping colour on SDSS only would compare a **colour-selected SDSS pool** to an **all-colour HST pool** — confounds PSF with stellar-population mix.
- Dropping colour on **both** is valid but loses parity with v2; we keep the cut on **both** instead.

**Why `CFHT_u − HSC_r` on HST (not `CFHT_u − ACS_F814W`):**
- SDSS cut is **u − r** in AB mags, not u − I.
- Farmer provides **`CFHT_u_MAG`** and **`HSC_r_MAG`** on the same sources — closest **u/r analogue** in COSMOS2020.
- **Shape and mag bins stay on ACS:** `b_a` and `ACS_F814W_MAG` unchanged; colour is a **sample selection** cut only (same role as SDSS `u-r` before reading `expAB_r`).

**HST colour validity:** require `CFHT_u_VALID == 1` and `HSC_r_VALID == 1` (Farmer flags untrustworthy flux/mag).

**Rejected:** `CFHT_u − ACS_F814W < 2.3` — that is **u − I**, not u − r; would not match SDSS semantics.  
**Rejected:** asymmetric cuts (colour on SDSS only) — breaks apples-to-apples.

**Sensitivity variant (not production):** no colour on either survey (`--no-color-cut`) for ablation.

---

## 0. Analysis template (what we are copying from v2 SDSS audit)

| v2 SDSS audit product | COSMOS analogue |
|----------------------|-----------------|
| `ba_mag_joint_panel_strict.png` — median raw b/a vs 0.5 mag bin | Same for HST and SDSS, side-by-side |
| `cosi_mag_joint_panel.png` — optional Hubble cos(i) overlay | **Deferred** for pass 1 (see §6) |
| `cut_funnel_v2.csv` — per-cut survival | SDSS only (HST has no lnL) |
| `cosi_mag_correlation.csv` — Spearman(mag, b/a) | Both surveys |
| `re_slice_mag_ba_panel.png` — b/a vs mag in Re bins | Both, using survey-native size column |
| Formal permutation / isotropic nulls | Phase 2 (after pass-1 trends) |

**Primary readout:** median **projected b/a vs magnitude** on each survey, **no cross-match required** for pass 1 (same sky volume, independent pipelines).

---

## 1. Footprint decision

### Chosen: **ACS/WFC F814W contiguous mosaic footprint**

| Bound | Value (J2000) | Source |
|-------|---------------|--------|
| RA | 149.43° – 150.80° | Scoville et al. 2007 ACS mosaic |
| Dec | 1.57° – 2.89° | same |
| Area | ~1.64 deg² | contiguous HST ACS I-band |

**Reasons:**
- HST shape columns in our chosen catalog (`ACS_*`) are measured on **ACS F814W** mosaics; outside this footprint ACS structural params are missing or unreliable.
- SDSS fully covers this box (Dec ≈ +2°).
- Using the full 2 deg² COSMOS2020 extent would include regions **without** ACS structural measurements → HST sample would be a mix of measured and NaN shapes, biasing medians.

**Rejected:** Full 2 deg² COSMOS box — asymmetric HST coverage.  
**Rejected:** Cone around (150.119°, 2.206°) only — arbitrary; ACS rectangle is the published HST limit.

**SQL box for SDSS** (CasJobs / `astroquery`):

```sql
WHERE p.ra BETWEEN 149.43 AND 150.80
  AND p.dec BETWEEN 1.57 AND 2.89
```

---

## 2. HST catalog decision

### Chosen (production): **COSMOS Zurich morphology v1.0** (`cosmos_morph_zurich_1.0.tbl`)

| Item | Value |
|------|--------|
| **File** | `plots/plots_null/v2/sdss_audit/COSMOS/data/cosmos_morph_zurich_1.0.tbl` |
| **Rows** | ~131,532 (full catalog); ~34,232 after strict ACS cuts + GIM2D |
| **Cite** | Sargent et al. 2007; COSMOS morphology collaboration |

### HST magnitude column

| Column | `ACS_MAG_AUTO` |
|--------|----------------|
| System | AB |
| Band | HST/ACS I (AUTO photometry) |

### HST axis ratio column

| Definition | `b_a = 1 - ELL_GIM2D` |
|------------|----------------------|
| Fitting | GIM2D on ACS I-band |
| Depth | Reliable to I_AB ≲ 22.5 |

### HST size column

| Column | `Re_arcsec = R_0P5_GIM2D` |
|--------|---------------------------|

### HST quality cuts (production pool)

| # | Cut | Column | Reason |
|---|-----|--------|--------|
| 1 | Footprint | `ra`, `dec` in ACS box | §1 |
| 2 | Galaxy classifier | `ACS_MU_CLASS == 1` | Drop stars |
| 3 | Star/galaxy | `STELLARITY == 0` | Galaxies only |
| 4 | Clean flag | `JUNKFLAG == 0` | Drop spurious |
| 5 | GIM2D fit | `FLUX_GIM2D > 0`, finite `ELL_GIM2D` | Successful morphology |
| 6 | Strict b/a | `b_a > 0.2` | Match SDSS strict null |
| 7 | Finite size | `Re_arcsec > 0` | Re-slice diagnostics |

**No colour cut** on HST (matched to SDSS `--no-color-cut` for this audit).

### Deferred: **COSMOS2020 Farmer** (`COSMOS2020_FARMER_R1_v2.2_p3.fits.gz`)

| Item | Value |
|------|--------|
| **URL** | https://irsa.ipac.caltech.edu/data/COSMOS/tables/cosmos2020/COSMOS2020_FARMER_R1_v2.2_p3.fits.gz |
| **Size** | ~2.0 GB compressed |
| **Rows** | ~964,506 (Farmer subset) |
| **Cite** | Weaver et al. 2022, ApJS 258, 11 + Farmer package (Zalesky et al. 2023) |

### Why Farmer (not alternatives)

| Option | Verdict | Reason |
|--------|---------|--------|
| **COSMOS2020 Classic** | Rejected for HST shape | Aperture photometry; **no** `ACS_A_WORLD` / `ACS_B_WORLD` profile axes on ACS |
| **ACS I-band catalog 2007** (Leauthaud; AstroPath feather) | Rejected as primary | ~1.1M sources, SExtractor `kron_radius` + `mu_class` only; **no** clean b/a column; older astrometry; already a subset of ACS. Useful as **sanity check** only |
| **Zurich morphology** (GIM2D) | **Chosen (production)** | `ELL_GIM2D = 1 − b/a`; reliable I ≲ 22.5; PSF-safe ACS reference |
| **COSMOS2020 Farmer** | Deferred | Download issues; superseded by user-supplied Zurich catalog |

### HST magnitude column

| Column | `ACS_F814W_MAG` |
|--------|-----------------|
| System | AB |
| Band | HST/ACS F814W (λ_c ≈ 8333 Å) |
| Provenance | Farmer/Tractor fit on ACS mosaic |

**Reasons:**
- Shape (`ACS_A/B_WORLD`) and mag are from the **same ACS F814W** fit → internally consistent (critical after Re-slice lesson: don’t mix mag from one band with shape from another).
- F814W is the canonical deep HST band in COSMOS (ACS treasury).

**Rejected:** `HSC_i_MAG` — ground-based, PSF ≈ 0.6″, defeats PSF-isolation goal.  
**Rejected:** Detection-band `izYJHKs` chi-mean mag — different image than ACS structural fit.

**Caveat (document, don’t fix in pass 1):** SDSS uses **`modelMag_r`** (r-band, λ ≈ 6230 Å). HST uses **I_F814W**. We compare **trend shapes** (slope of median b/a vs mag), not absolute mag offsets. Pass 1 is **not** a flux-limited matched sample across surveys.

### HST axis ratio column

| Definition | `b_a = ACS_B_WORLD / ACS_A_WORLD` |
|------------|-----------------------------------|
| Units | Both in degrees on sky → ratio dimensionless |
| Convention | B = semi-minor, A = semi-major → **b/a ∈ (0, 1]** |

**Reasons:**
- Direct projected axis ratio, analogous to SDSS `expAB_r`.
- From Farmer profile fit on **0.03″/pixel** ACS data (effective PSF FWHM ≈ **0.09–0.12″** in I-band) vs SDSS **~1.4″** FWHM at **0.396″/pix**.

**Rejected:** `ACS_ELONGATION` from raw ACS2007 (a/b not b/a) — would need inversion; Farmer axes are cleaner.  
**Rejected:** `ELL_GIM2D` from Zurich — different fitting code, depth limit.

### HST size column (Re-slice diagnostics)

| Column | `Re_arcsec = ACS_A_WORLD × 3600` |
|--------|----------------------------------|
| Meaning | Semi-major axis of ACS F814W profile (proxy for effective radius) |

**Reasons:** Same Farmer ACS fit as b/a; enables Re bins [1,2), [2,3), [3,∞)″ like v2 audit.

**Caveat:** This is **semi-major axis of the fitted profile**, not Petrosian or exp `expRad_r`. Use **only within HST** or **only within SDSS** slices, never compare absolute Re across surveys.

### HST quality cuts (production pool)

Applied in order:

| # | Cut | Column | Reason |
|---|-----|--------|--------|
| 1 | Footprint | `ra`, `dec` in ACS box | §1 |
| 2 | Finite ACS photometry | `ACS_F814W_MAG` finite, `ACS_F814W_FLUX > 0` | Need mag |
| 3 | Finite shape | `ACS_A_WORLD > 0`, `ACS_B_WORLD > 0` | Need b/a |
| 4 | Galaxy classifier | `ACS_MU_CLASS == 1` | Drop stars (same spirit as SDSS `type=3`) |
| 5 | Fake detections | `ACS_MU_CLASS != 3` | IRSA recommendation (CR/spurious) |
| 6 | Farmer model OK | `MODEL_FLAG == 0` | Drop non-converged / >0.6″ drifted fits |
| 7 | **Colour (symmetric)** | `CFHT_u_MAG − HSC_r_MAG < 2.3` | Match SDSS `u-r < 2.3` (§0) |
| 8 | Colour valid | `CFHT_u_VALID == 1`, `HSC_r_VALID == 1` | Trustworthy u and r in Farmer |
| 9 | Strict b/a | `b_a > 0.2` | Match SDSS strict null (`Q0 = 0.2`) |
| 10 | Clean imaging mask | `FLAG_COMBINED == 0` | Inside clean UltraVISTA+HSC+Suprime region per COSMOS2020 |

**Optional sensitivity cut (not in production):** `ACS_F814W_MAGERR < 0.1` — run as variant if scatter is huge.

**Explicitly NOT applied to HST:**
- No lnL exp-wins (N/A on Farmer ACS).

---

## 3. SDSS catalog decision

### Chosen: **SDSS DR17 `PhotoObj` via CasJobs or `astroquery.sdss.SDSS.query_sql`**

| Item | Value |
|------|--------|
| Table | `PhotoObj` with **`mode = 1`** (primary; matches v2) |
| Release | **DR17** (`data_release=17` in astroquery; CasJobs context `DR17`) |

**Reasons:**
- User-specified DR17 (final SDSS-IV imaging release; same schema as DR16 for PhotoObj shape columns).
- Same production cuts as v2 null audit (`null_catalog_utils.py`), including **`u-r < 2.3`**.
- CasJobs handles ~10⁵ rows in COSMOS box without SkyServer timeout.
- `astroquery.sdss` defaults to DR17.

**Note:** COSMOS was observed early in SDSS; DR16 vs DR17 differences in this box are expected to be negligible (recalibration only), but we use **DR17** consistently for this project.

**Rejected:** DR16 — superseded by user choice.  
**Rejected:** Archived v0 `SDSS_catalog.csv` (COSMOS) — **no RA**, known provenance issues.  
**Rejected:** SAS `photoObj` FITS field files — unnecessary complexity for first pass.

### SDSS columns (mirror v2)

| Role | Column |
|------|--------|
| ID | `objID` |
| Position | `ra`, `dec` |
| Magnitude (mag bins) | **`modelMag_r`** |
| Axis ratio | **`expAB_r`** |
| Size | **`expRad_r`** (arcsec, exponential) |
| Colour | `modelMag_u`, `modelMag_g` → `u_r` |
| Morphology | `lnLExp_r`, `lnLDeV_r` |
| Quality | `type`, `clean`, `mode` |

### SDSS production cuts (identical to v2 strict null)

From `null_catalog_utils.py` / `prepare_null_strict_color_base`:

| # | Cut | Value |
|---|-----|-------|
| 1 | Footprint | ACS box (§1) |
| 2 | Object type | `type = 3` (GALAXY) |
| 3 | Quality | `clean = 1`, `mode = 1` |
| 4 | Colour | `u-r < 2.3` (`modelMag_u − modelMag_r`; `SDSS_UR_MAX_CDF`) | §0 — symmetric with HST |
| 5 | Profile winner | `lnLExp_r > lnLDeV_r` (finite lnL) |
| 6 | Strict b/a | `expAB_r > 0.2` |
| 7 | Finite | `modelMag_r`, `expAB_r`, `expRad_r` > 0 |

**Reason:** Any SDSS trend difference vs HST must not come from **different cuts** than the v2 audit we are comparing against.

---

## 4. Analysis design: separate vs cross-matched

### Pass 1 (primary): **Separate ensemble medians**

- Build HST pool and SDSS pool independently in the same footprint.
- Plot median b/a vs mag for each.
- Compare **functional form** (slope, faint-end floor, bin-to-bin Δmedian).

**Reason:** Science question is whether **SDSS as a pipeline** introduces mag-dependent b/a bias relative to a PSF-safe reference population in the **same volume**. That is a **distribution-level** test; forcing 1:1 cross-matches introduces match-radius bias, different depth limits, and blended-source ambiguity.

### Pass 2 (optional): **Matched sample**

- Cross-match Farmer `ra,dec` ↔ SDSS `ra,dec`, radius **0.5″** (tune with histogram).
- Require both pass production cuts.
- Plot Δ(b/a) = SDSS `expAB_r` − HST `b_a` vs mag.

**Use when:** Pass 1 shows different trends and we need object-level confirmation.

---

## 5. cos(i) / Hubble mapping — deferred in pass 1

| Decision | **Pass 1: raw b/a only** (like `ba_mag_joint_panel*.png`) |
|----------|-----------------------------------------------------------|

**Reasons:**
- HST Farmer b/a is already a direct axis ratio; applying Hubble with `q0=0.2` assumes thin-disk geometry not encoded in Farmer fits.
- SDSS v2 audit showed much of the faint-end pathology is in **raw `expAB_r`**, not cos(i) (floor at 0.05, Re confounds).
- Keeps HST and SDSS comparison on the **same observable**: projected b/a.

**Pass 2:** Add `cos(i)` from Hubble for SDSS only if comparing to FRB hosts; add HST inclination only if we run GALFIT on cutouts.

---

## 6. Expected outcomes and interpretation rules

| Observation | Interpretation |
|-------------|----------------|
| SDSS median b/a drops with fainter `modelMag_r`; HST **flat** to mag ~21–22 | Supports **PSF and/or SDSS PhotoObj pipeline** bias |
| Both slopes similar | PSF unlikely dominant; composition, selection, or both pipelines share systematics |
| HST faint-end pile-up at small b/a | HST pipeline/SNR limit (less likely — ACS depth to ~28 AB) |
| SDSS hits 0.05 floor; HST does not | Classic SDSS quantization artifact; strong PSF/pipeline smoking gun |
| Different N(m) between surveys | Expected (different depth, cuts); compare trends not raw counts |

**PSF scale reminder:**

| Survey | Pixel scale | PSF FWHM (typical) |
|--------|-------------|-------------------|
| ACS F814W | 0.03″ (drizzled) | ~0.1″ |
| SDSS r | 0.396″ | ~1.4″ |

Galaxies with `expRad_r` ≲ 1–2″ are **resolved in HST** but **barely resolved in SDSS** — Re-slice plots are diagnostic, not PSF-proof on their own.

---

## 7. Data acquisition plan (when implemented)

### HST
1. Download `COSMOS2020_FARMER_R1_v2.2_p3.fits.gz` from IRSA.
2. Read with `astropy.table.Table`; select columns in §2.
3. Apply footprint + cuts; write `plots/plots_null/v2/sdss_audit/COSMOS/cosmos_hst_farmer_acs_strict.csv`.

### SDSS
1. CasJobs / `query_sql` box query with columns in §3.
2. Apply same cuts in Python (`null_catalog_utils`) or SQL.
3. Write `plots/plots_null/v2/sdss_audit/COSMOS/cosmos_sdss_dr17_strict.csv`.

### Plots → `plots/plots_null/v2/sdss_audit/COSMOS/plots/` (mirror `plots/plots_null/v2/sdss_audit/`)

---

## 8. File inventory

| File | Description |
|------|-------------|
| `plots/plots_null/v2/sdss_audit/COSMOS/CATALOG_DECISIONS.md` | This document |
| `plots/plots_null/v2/sdss_audit/COSMOS/ZURICH_CATALOG_NOTES.md` | Zurich column map and cut funnel |
| `plots/plots_null/v2/sdss_audit/COSMOS/data/cosmos_morph_zurich_1.0.tbl` | Raw Zurich morphology table |
| `plots/plots_null/v2/sdss_audit/COSMOS/cosmos_hst_zurich_strict.csv` | HST production pool |
| `plots/plots_null/v2/sdss_audit/COSMOS/cosmos_sdss_dr17_nocolor_strict.csv` | SDSS production pool (no colour) |
| `plots/plots_null/v2/sdss_audit/COSMOS/cosmos_sdss_dr17_strict.csv` | SDSS with u−r cut (legacy) |
| `plots/plots_null/v2/sdss_audit/COSMOS/plots/hst_ba_mag_joint_panel.png` | HST median b/a vs mag |
| `plots/plots_null/v2/sdss_audit/COSMOS/plots/sdss_ba_mag_joint_panel.png` | SDSS median b/a vs mag |
| `plots/plots_null/v2/sdss_audit/COSMOS/plots/hst_sdss_ba_mag_overlay.png` | Both on one axes |
| `plots/plots_null/v2/sdss_audit/COSMOS/plots/hst_sdss_ba_mag_overlay_mag22p5.png` | Overlay at GIM2D reliable depth |
| `plots/plots_null/v2/sdss_audit/COSMOS/plots/cosi_cdf_mag*.png` | Field-level cos(i) CDF comparison |
| `plots/plots_null/v2/sdss_audit/COSMOS/plots/matched_delta_ba_vs_mag.png` | Pass-2 matched Δ(b/a) |
| `scripts/build_cosmos_sdss_catalog.py` | SDSS DR17 box query + strict cuts |
| `scripts/build_cosmos_hst_zurich_catalog.py` | HST Zurich GIM2D pool |
| `scripts/plot_cosmos_ba_mag_audit.py` | Pass-1 b/a–mag panels + Re slices |
| `scripts/plot_cosmos_null_cdfs.py` | cos(i) CDF overlays |
| `scripts/plot_cosmos_matched_ba.py` | Pass-2 cross-match diagnostics |

---

## 9. Open items

1. **Pass-2 cross-match** radius (0.5″ vs 1.0″).
2. **Images** — not needed for pass 1; GALFIT on COSMOS cutouts is phase 3.

**Resolved:**
- ~~DR16 vs DR17~~ → **DR17**
- ~~HST colour cut~~ → **`CFHT_u − HSC_r < 2.3` on HST; `modelMag_u − modelMag_r < 2.3` on SDSS** (symmetric)

---

## 10. References

- Weaver et al. 2022, ApJS 258, 11 — COSMOS2020
- Zalesky et al. 2023, ApJS — The Farmer
- Leauthaud et al. 2007 — ACS COSMOS catalog
- Scoville et al. 2007, ApJS 172, 38 — ACS mosaic
- Stoughton et al. 2002; Abazajian et al. 2009 — SDSS modelMag / PhotoObj
- Repo: `NULL_CATALOG_DATA_AND_INCLINATION_AUDIT.md`, `plots/plots_null/v2/sdss_audit/`
