# Null catalog and inclination CDF methodology (v1)

This document describes how the v1 null comparison catalogs are built, which cuts remove which rows, and how inclination CDF plots are produced. Machine-readable funnel counts live in [`cut_funnel.csv`](cut_funnel.csv) (regenerate with `python scripts/audit_and_plot_null_v1_diagnostics.py --full`).

Pre-morphology-cut plots are archived under `Archive/plots_null_pre_morphology_cut/`.

---

## 1. Catalog files and row counts

| File | Build script | Rows on disk |
|------|----------------|--------------|
| [`SDSS_catalog_v1_allsky_modelmr.csv`](../../../SDSS_catalog_v1_allsky_modelmr.csv) | [`scripts/build_sdss_null_catalog.py`](../../../scripts/build_sdss_null_catalog.py) | **501,885** |
| [`LS_catalog_v1_allsky_modelmr.csv`](../../../LS_catalog_v1_allsky_modelmr.csv) | [`scripts/build_legacy_catalog_csv.py`](../../../scripts/build_legacy_catalog_csv.py) | **500,000** |

Both builders target ~500k objects. Neither applies null CDF cuts at build time; cuts are applied in [`scripts/null_catalog_utils.py`](../../../scripts/null_catalog_utils.py) when plotting.

**Footprint (both surveys):** right ascension 0–360°, declination **−30° to +90°** (`JOINT_DEC_MIN` / `JOINT_DEC_MAX`).

---

## 2. SDSS catalog provenance

**Source:** SDSS DR16 `PhotoObj` via `astroquery.sdss.SDSS.query_sql`, chunked in RA.

**Profile winner (r-band):** `model_winner_is_exp` from **`lnLExp_r > lnLDeV_r`** (strict; ties count as deV). Patched into the v1 CSV with [`scripts/patch_sdss_profile_winner.py`](../../../scripts/patch_sdss_profile_winner.py). Validation: [`diagnostics/sdss_profile_winner/`](sdss_profile_winner/).

**Null CDF defaults (SDSS):** magnitude **`modelMag_r`**; axis ratio **`expAB_r` only** (after dropping deV winners). Do not use `best_model_ba_r` in CDF pools.

---

## 3. Legacy Survey catalog provenance

**Source:** Legacy Survey DR10 `ls_dr10.tractor` via NOIRLab TAP.

**Null CDF defaults (Legacy):** magnitude **`tractor_mag_r`**; axis ratio **`expAB_r`**; exclude **`REX` and `DEV`**; keep **`tractor_type == EXP`** or **Sérsic `rdVrad` ∈ [0.75, 2]**.

**Diagnostics** (mag vs b/a, REX fraction): still use the mag+color pool with **REX visible** (REX-only type exclusion on hexbins, not the full CDF morphology cut).

Morphology diagnostic: [`diagnostics/legacy_morphology/`](legacy_morphology/).

---

## 4. Null sample cuts

Implemented in [`scripts/null_catalog_utils.py`](../../../scripts/null_catalog_utils.py). CDF order:

1. Magnitude cut (`mag <= limit`)
2. **SDSS:** \(u-r < 2.3\); **Legacy:** \(g-r < 0.75\)
3. **SDSS:** drop deV profile winners (`lnLExp_r > lnLDeV_r`)
4. **Legacy:** drop `REX`, `DEV`; keep EXP or \(0.75 \le n \le 2\)
5. **Strict:** \(b/a > 0.2\)

### 4.3c Late-type morphology cuts (CDF pools)

| Survey | Rule | Axis ratio for Hubble \(i\) |
|--------|------|-----------------------------|
| **SDSS** | Drop all deV-winning galaxies; then **`expAB_r` only** | `expAB_r` |
| **Legacy** | **`EXP`** or **`rdVrad` ∈ [0.75, 2]**; no **`DEV`** in CDF pools | `expAB_r` |

Hubble formula with \(q_0 = 0.2\) (see §4.4 in prior docs).

---

## 5. Cut funnel — final pool sizes (post morphology)

From [`cut_funnel.csv`](cut_funnel.csv) (SDSS: `modelMag_r` + `expAB_r`; Legacy: `rmag` + `expAB_r`).

### m < 21 (reference)

| Survey | After mag | After color | After morphology | After \(b/a>0.2\) | **Final** |
|--------|-----------|-------------|------------------|-------------------|-----------|
| Legacy | 43,113 | 5,980 (g−r) | 4,577 (EXP∪n) | 4,374 | **4,374** |
| SDSS | 156,343 | 29,844 (u−r + lnL exp) | — | 27,474 | **27,474** |

**Notes:**

- SDSS lnL patch matches ~77% of catalog rows on disk; unmatched rows lack `lnL` and are dropped at the exp-winner stage.
- Legacy CDF pool at m&lt;21 is **&lt; 5,000** after morphology; see `legacy_morphology_summary.md` for re-query guidance.

---

## 6. CDF plotting mechanism

**Scripts:**

- Mag cuts: [`scripts/plot_null_mag_cut_cdfs.py`](../../../scripts/plot_null_mag_cut_cdfs.py) → `mag_cuts/mag{15–24}/legacy_strict/` and `sdss_strict/`.
- Pipeline: [`scripts/plot_pipeline_diagnostics.py`](../../../scripts/plot_pipeline_diagnostics.py).

**FRB hosts:** [`pipeline_galfit_results.csv`](../../../pipeline_galfit_results.csv); mag + GALFIT \(b/a > 0.2\); no host color cut.

---

## 7. Column quick reference

| Column | Survey | CDF use |
|--------|--------|---------|
| `modelMag_r` | SDSS | Mag cut |
| `expAB_r` | SDSS | Shape / Hubble \(i\) |
| `lnLExp_r`, `lnLDeV_r` | SDSS | Profile winner |
| `model_winner_is_exp` | SDSS | Filter (1 = keep) |
| `tractor_mag_r` | Legacy | Mag cut |
| `expAB_r` | Legacy | Shape |
| `rdVrad` | Legacy | Sérsic \(n\) for morphology |
| `tractor_type` | Legacy | EXP / REX / DEV / … |

---

## 8. Changelog

- **2026-05-27:** Late-type morphology cuts — SDSS lnL exp-wins + `expAB_r`; Legacy EXP∪n, no REX/DEV in CDF pools. Pre-cut plots archived under `Archive/plots_null_pre_morphology_cut/`.
