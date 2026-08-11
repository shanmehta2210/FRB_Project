# DES calibration + dust extinction research notes

Sources: Tarsitano et al. 2018 (MNRAS, [arXiv:1807.10767](https://arxiv.org/abs/1807.10767)); DES Y1 morph FITS on NCSA; literature on inclination-dependent dust.

Working sample: `catalog/DES_y1_morph_sample_500k.csv` (fixed 2026-07-14: `mag_r` no longer double-applies `MAG_CAL`).

---

## 1. Verdict: how to use Tarsitano calibration

**Use the main parametric columns as-is. They are already calibrated.**

| Quantity | Use this column | `*_CAL_*` means |
|----------|-----------------|-----------------|
| Magnitude | `MAG_SERSIC_X` | `MAG_CAL_X` = η already applied |
| Size | `RE_X` | `RE_CAL_X` = η already applied |
| Sérsic n | `N_SERSIC_X` | `N_SERSIC_CAL_X` = η already applied |
| Ellipticity | `ELLIPTICITY_SERSIC_X` | `ELLIPTICITY_SERSIC_CAL_X` = η already applied |
| Axis ratio | `ba = 1 − ELLIPTICITY_SERSIC_X` | — |

From Appendix B (paper):

> MAG_SERSIC_X — GALFIT value… **already includes** the calibration listed in MAG_CAL_X.  
> RE_X … **already calibrated**. The correction is reported in RE_CAL_X.  
> N_SERSIC_X … **calibrated**; correction in N_SERSIC_CAL_X.  
> ELLIPTICITY_SERSIC_X … **corrected**; calibration accessible through ELLIPTICITY_SERSIC_CAL_X.  
> Absence of calibration → correction value set to **99**.  
> Science-ready: `FIT_STATUS_X = 1` (≡ validated+calibrated pack). Cleaner: also `SN_X > 30`. Reliable morph typically for `MAG_AUTO ≲ 21–21.5`.

**Wrong uses we tried earlier**

1. `ba = 1 − ELLIPTICITY_SERSIC_CAL` — CAL is η, not calibrated ε → ba often >1 / nonsense.
2. `mag = MAG_SERSIC + MAG_CAL` — **double-counts** (builder bug; fixed). Shift was small (~0.02–0.03 mag median) but wrong in principle.
3. Ignoring calibration entirely for ellipticity — unnecessary: the ellipticity column is already the calibrated one.

**Correct derived products for this project**

```text
ba_r  = 1 - ELLIPTICITY_SERSIC_R          # already calibrated
mag_r = MAG_SERSIC_R                      # already calibrated
n_r   = N_SERSIC_R                        # already calibrated
# optional diagnostics:
eps_uncal = ELLIPTICITY_SERSIC_R - ELLIPTICITY_SERSIC_CAL_R
mag_uncal = MAG_SERSIC_R - MAG_CAL_R
```

Builder: `scripts/build_des_y1_morph_sample.py` now matches Appendix B.

---

## 2. Exact calibration methodology (Tarsitano+2018 §4.4)

### 2.1 Goal

GALFIT on real DES coadds is biased by PSF, noise, and parameter covariances (esp. faint / low S/N). They quantify bias with **UFig-BCC** image simulations (~10⁷ objects, DES-Y1-like), then apply cell-wise corrections.

### 2.2 4D calibration grid

Nodes (hypervolume “cells”):

| Axis | Nodes |
|------|--------|
| mag | [14.5, 23.5] step 1 |
| size r | [0.5, 16.5] px step 2 |
| n | {0.2, 2, 4, 10} |
| ε | {[0,0.3), [0.3,0.6), [0.6,1]} |

In each cell, model (truth) and fit (measured) distributions have medians \(\hat\mu^i\), \(\mu^i\) and scatters \(\hat\sigma^i\), \(\sigma^i\) for \(i \in \{mag, r, n, \varepsilon\}\).

### 2.3 Correction definition (Eq. 8)

\[
\eta^{i} = \hat\mu^{i} - \mu^{i}
\qquad\text{(truth median − measured median)}
\]

So for an object in that cell:

\[
p_{\mathrm{cal}} = p_{\mathrm{meas}} + \eta^{i}
\]

Catalog stores **both** \(p_{\mathrm{cal}}\) (main column) and \(\eta^{i}\) (`*_CAL_*`).

Cell reliability weight (Eq. 9):

\[
w = \sqrt{\sum_i (\hat\sigma_i / \hat m_i)^2}
\]

Large \(w\) → less trustworthy cell correction (maps use squares/pentagons).

### 2.4 What the simulations show

- Magnitude usually well recovered (\(\eta^{mag}\sim 10^{-3}\) for ~99% of converged fits).
- **GALFIT tends to recover larger sizes and ellipticities** (and often larger n) at low S/N → typical \(\eta^{\varepsilon}<0\) (subtract ellipticity when calibrating).
- Cutting **S/N < 30** removes most size/mag outliers.
- Calibration maps (Fig. 9 + App. A for g,r): strongest size/n corrections at faint mag; many cells need little correction.

### 2.5 Empirical check on our 500k sample

Valid CAL (|η|<10): N=346 328 with 15 < mag < 23 and ε∈[0,1].

| Cut | N | median ba (cal) | median ba (uncal) | median cos i (q₀=0.2, cal) |
|-----|---:|----------------:|------------------:|---------------------------:|
| all | 346328 | 0.680 | 0.483 | 0.663 |
| SN≥30 | 103132 | 0.676 | 0.593 | 0.659 |
| mag≤21 | 58797 | 0.609 | 0.573 | 0.587 |
| mag≤21.5 | 95814 | 0.607 | 0.558 | 0.585 |
| SN≥30 & mag≤21.5 | 78284 | 0.652 | 0.604 | 0.633 |
| 0.4<n<1.5 | 152413 | 0.690 | 0.489 | 0.674 |
| 0.4<n<1.5 & SN≥30 | 46995 | 0.632 | 0.556 | 0.612 |
| 0.4<n<1.5 & mag≤21 | 23435 | 0.616 | 0.579 | 0.594 |
| 0.4<n<1.5 & mag≤21.5 & SN≥30 | 35932 | 0.614 | 0.571 | 0.592 |

η^ε by S/N (median): SN[0,20) → −0.16; [20,30) → −0.10; [30,50) → −0.057; [50,100) → −0.037; ≥100 → −0.019.

**Interpretation:** uncalibrated ba (~0.48) looks closer to LS EXP (~0.45) only because GALFIT **overestimates** ellipticity before correction. The paper’s science product is the **calibrated** ba (~0.61–0.69 depending on cuts). Do not use uncalibrated values to “match” LS.

### 2.6 Why DES calibrated ba still ≫ LS EXP

Calibration does **not** erase the DES–LS gap. Remaining drivers (from prior audit):

1. **LS Tractor REX gate:** free ellipticity must beat round REX; face-ons → REX (ba≡1) and are absent from `type=EXP` → EXP sample edge-on biased.
2. **Population:** free-n DES Sérsic sample mixes disks + early-types (early-types peak ba~0.7–0.8).
3. **Depth / S/N:** paper itself says full morph reliability ~mag ≲ 21–21.5; our faint-heavy sample inflates scatter and residual bias even after η.

Science-ready DES EXP-analogue for comparison: `0.4 < n < 1.5`, `SN_R > 30`, `MAG_SERSIC_R ≤ 21.5` → median ba ≈ **0.614**, median cos i ≈ **0.59**.

---

## 3. Dust extinction and inclination / axis-ratio samples

### 3.1 Effect (why median cos i ~ 0.6 can look “natural”)

For dusty disks, edge-on sightlines have higher optical depth:

- Edge-ons appear **fainter** (and redder) than face-ons of the same intrinsic luminosity.
- In a **magnitude-limited** catalog, edge-ons drop out preferentially.
- Observed ba / cos i distributions become **face-on weighted** → median cos i rises above the thin-disk isotropic value (~0.5 for ba, or ~0.5 for cos i if flat).

User note matches literature: edge-on spirals are shifted faintward; the surviving sample looks rounder / more face-on. Median cos i ~ 0.6 is a common symptom of this selection, not proof of a wrong inclination formula.

Intrinsic thickness (q₀ ~ 0.1–0.25) also truncates the most edge-on ba, raising the median further.

### 3.2 Classic extinction laws (catalog / TF corrections)

| Form | Typical use |
|------|-------------|
| \(\Delta m \propto \log(a/b)\) or \(\propto \sec i\) | Older RC3 / Holmberg-style catalog corrections |
| \(\Delta M_r = 1.27(\log q)^2\) | Unterborn & Ryden 2008 (SDSS exponential disks); better than pure \(\propto\log q\) |
| \(\tau \propto \cos i\) (optically thick slab) | Shao et al. 2007: recover face-on LF + band-dependent γ |
| \(\Delta m(i)\) from radiative transfer (Tuffs/Popescu) | Driver et al. 2007: separate bulge vs disc attenuation |
| Empirical \(\eta_g \approx 1.2\) mag at b/a=0.1 | Yip et al. 2010 (composite spectra) |

Face-on magnitude (Unterborn & Ryden):

\[
M_r^{f} = M_r - 1.27(\log q)^2
\]

Then rebuild luminosity functions / absolute-mag cuts on \(M_r^{f}\).

Driver et al. (MGC): fit attenuation vs \(1-\cos i\) empirically, then use Tuffs et al. RT models to also remove residual **face-on** opacity (\(\tau_B^f \approx 3.8\pm0.7\) centrally — optically thick centers). Correct disc and bulge separately (bulge light seen through dusty disc).

### 3.3 How literature usually “corrects” for dust in shape / inclination work

Practical recipes used in papers (often combined):

1. **Correct magnitudes to face-on, then re-cut**  
   Apply \(\Delta M(q)\) (or RT table), redefine the magnitude-limited sample on face-on mag. Removes the first-order selection bias against edge-ons (Unterborn & Ryden; TF distance work).

2. **Volume-limited samples**  
   Absolute-mag + redshift box deep enough that dimming does not eject edge-ons from the volume (Yip et al.; many SDSS shape papers). Still need care near the faint absolute limit.

3. **Model the selection function**  
   Fit inclination-binned LFs; infer optical depth and the face-on LF jointly (Shao et al. 2007). Explicitly predicts under-representation of small-ba objects.

4. **Wavelength choice**  
   NIR / mid-IR selection strongly reduces dust bias (2MASS, WISE). Optical morph catalogs remain dust-sensitive.

5. **Do not confuse dust selection with measurement bias**  
   Dust changes **who enters the sample**. PSF/GALFIT bias changes **measured ba** of those who enter. Tarsitano η corrects (2); Unterborn-style \(\Delta M(q)\) corrects (1). Both can raise median cos i; they are different physics.

6. **Shape-distribution modeling with dust**  
   Padilla & Strauss 2008 and related SDSS shape papers explicitly note dust under-represents edge-ons in mag-limited samples when inferring intrinsic thickness / ellipticity.

### 3.4 Implication for our null / inclination CDFs

- A mag-cut CDF of cos i with median ~0.6 for disk-like galaxies is **compatible with dusty, magnitude-limited selection**, especially if no face-on mag correction is applied.
- LS EXP median ba ~0.45 / cos i ~0.4–0.5 is **not** the dust-free isotropic disk expectation either: Tractor REX filtering removes face-ons, pushing the opposite way.
- Fair survey comparison needs: (a) same morph class, (b) science-ready S/N and mag, (c) either volume limits or explicit \(\Delta M(q)\) re-selection, (d) awareness that LS EXP ≠ “all exponential disks.”

Suggested next test (not run here): apply Unterborn \(\Delta m_r = 1.27(\log ba)^2\) to DES/LS apparent mags, keep objects with \(m^{f} \le m_{\lim}\), recompute median ba and cos i CDFs.

**Result + recommended baseline:** [`../DUST_AND_MEDIAN_BA.md`](../DUST_AND_MEDIAN_BA.md).  
**Full staged plan (archived):** [`Archive/notes/EXTINCTION_CORRECTION_PLAN.md`](../../../../Archive/notes/EXTINCTION_CORRECTION_PLAN.md) (stub: [`EXTINCTION_CORRECTION_PLAN.md`](EXTINCTION_CORRECTION_PLAN.md)).

---

## 4. Why `ba_vs_mag_scatter` looks discontinuous at faint mag

Plot: `plots/plots_null/v2/des_audit/ba_vs_mag_scatter.png` (full 500k, **no** science-ready cuts). The blocky / stepped look at mag ≳ 22 is **real in the catalog product**, not a plotting bug and not evidence we applied η wrong.

### 4.1 What the data do (500k sample, calibrated columns)

| mag bin | N | median ba | median SN_R | median η^ε | note |
|---------|--:|----------:|------------:|-----------:|------|
| 21.00–21.25 | 16135 | 0.595 | 44.9 | −0.036 | still near science-ready |
| 21.25–21.50 | 20882 | 0.610 | 37.0 | −0.043 | paper edge ≈ 21.5 |
| 21.50–21.75 | 24338 | 0.665 | 30.1 | −0.082 | ba jump; SN~30 |
| 22.00–22.25 | 37646 | 0.653 | 19.9 | −0.097 | almost all SN<30 |
| 22.50–22.75 | 52905 | 0.725 | 13.4 | **−0.239** | discrete η jump |
| 23.25–23.50 | 50528 | 0.705 | 7.6 | −0.341 | deep in unreliable zone |
| 24.00–24.25 | 5962 | 0.695 | 4.8 | −0.508 | |
| 24.50–24.75 | 2805 | **0.156** | 4.4 | **0.000** | cliff: Re~0.02 px junk |

Two visual features:

1. **Stepped density / ba jumps near ~21.5 and ~22.5**  
   - Crossing the paper’s reliability edge (`MAG ≲ 21.5`, `SN > 30`).  
   - **Piecewise 4D calibration grid** uses mag nodes every **1 mag** ([14.5, 23.5]). Objects in adjacent cells get different constant η^ε → discrete jumps in calibrated ba (η^ε median −0.11 → −0.24 across 22.25–22.75).  
   - Low-S/N GALFIT also biases n/Re/ε before η is applied; η cannot fully erase that.

2. **Hard cliff near mag ≈ 24.5**  
   - Median ba collapses 0.69 → 0.15; Re collapses to ~0.02 px; η^ε → 0 (outside / no usable calibration cell).  
   - These are pathological faint fits, not a physical edge-on population. Horizontal stripes at ba~0.05 are GALFIT hitting near-zero axis-ratio extremes.

Also: ~3.3% of the full sample has ba ≤ 0.05 (fit-bound debris), almost all at faint mag.

### 4.2 How we should use the catalog for this project

The discontinuity is a **warning that the audit plot includes objects the authors say not to trust**. For inclination / null CDFs:

- Prefer **`SN_R > 30` and `MAG_SERSIC_R ≤ 21.5`** (paper science-ready; App. B also recommends `FIT_STATUS_R = 1` when available).  
- Under those cuts, median ba ≈ **0.65** (full types) / ≈ **0.61** for `0.4 < n < 1.5` — smooth, no cliff.  
- Do **not** interpret mag≳23 structure in the uncut scatter as astrophysics.

---

## 5. Reconfirmation: are we using Y1 morph correctly, and is this the right catalog?

### 5.1 Column usage — confirmed correct

| Our product | Source | Status |
|-------------|--------|--------|
| `ba_r` | `1 − ELLIPTICITY_SERSIC_R` | Correct: main column **already calibrated** (App. B) |
| `mag_r` | `MAG_SERSIC_R` | Correct: **already includes** `MAG_CAL_R` (fixed double-add) |
| `n_r`, `re_r` | `N_SERSIC_R`, `RE_R` | Correct: already calibrated |
| `*_CAL_*` | stored η for diagnostics | Correct: **not** drop-in calibrated values |

Ellipticity convention: paper defines ε = 1 − (b/a) from GALFIT axis ratio → `ba = 1 − ε` is the intended axis ratio.

### 5.2 Is Tarsitano DES Y1 morph the catalog we should use?

**Yes, for this application** (compare DES structural axis-ratio / Sérsic-disk pools to LS Tractor EXP and HSC Kawinwanichakij Sérsic fits).

| DES product | What it is | Use for ba / inclination null? |
|-------------|------------|--------------------------------|
| **Tarsitano+2018 Y1 morph** (what we use) | Public GALFIT single-Sérsic + ZEST+; n, Re, mag, **axis ratio**, PA; sim-calibrated | **Yes** — only large public DES catalog built for structural morphology |
| DES Y1/Y3 **metacalibration** shear catalogs | Weak-lensing shapes (ε for cosmology) | **No** — not Sérsic structural ba; different ε definition and selection |
| DES **MOF / SOF / fitvd** (Y3/Y6 Gold photometry) | Multi-epoch photometry / bulge–disk flux models for cosmology | **Not a drop-in** morph catalog; no Tarsitano-style public Sérsic ba release for Y3/Y6 |
| DESDM SExtractor only | Detection + Kron/FLUX_RADIUS | Too crude for inclination CDFs |

There is **no later DES data release that supersedes Tarsitano with a comparable all-sky Sérsic morph catalog**. Y3/Y6 advanced lensing and photometry pipelines; they did not republish a GALFIT morph catalog of this type. So for DES structural ba vs LS/HSC, Y1 morph is the right public choice.

### 5.3 Caveats (honest, not blockers)

1. **Depth / area**: Y1 (~1800 deg², shallower) vs full LS / HSC — compare at matched mag and S/N, not raw full samples.  
2. **REX-like selection absent**: DES keeps free-ε Sérsics; LS `type=EXP` excludes face-on REX winners → DES ba medians higher even when morph is correct.  
3. **Our builder uses `FIT_AVAILABLE_R==1`**, not the fuller App. B `FIT_STATUS_R==1` pack (outliers / overlap / S/G). Prefer adding `FIT_STATUS` + `SN>30` + mag≤21.5 for production CDFs when that column is ingested.  
4. **Uncut audit plots will look discontinuous** past mag~21.5 by construction of the catalog + calibration grid; that does not invalidate the calibrated columns inside the science-ready regime.

### 5.4 Absolute short verdict

- Calibration approach: **correct** (use main columns; do not add `*_CAL_*` again).  
- Catalog choice: **correct** for DES structural morphology / inclination comparison.  
- Faint-mag discontinuity in the full-sample scatter: **expected artifact** of leaving the recommended mag/S/N regime + discrete η cells — restrict science samples accordingly.

---

## 6. Action items taken

1. Documented Tarsitano η methodology and Appendix B column semantics.
2. Confirmed `ELLIPTICITY_SERSIC_*` / `MAG_SERSIC_*` are the calibrated products; `*_CAL_*` are applied η.
3. Fixed DES builder + CSVs: `mag_r = MAG_SERSIC_R` (no `+ MAG_CAL`).
4. Reported calibrated vs uncalibrated ba/cos i under paper-like cuts.
5. Summarized standard dust-correction approaches for inclination-selected catalogs.
6. Diagnosed faint-mag discontinuity in DES `ba_vs_mag_scatter` (calibration grid + reliability cliff + junk beyond ~24.5).
7. Reconfirmed Y1 morph is the appropriate DES catalog for this application vs metacal/MOF.

---

## 7. Key references

- Tarsitano et al. 2018, MNRAS, morph catalog + UFig-BCC calibration ([arXiv:1807.10767](https://arxiv.org/abs/1807.10767)); App. B column dictionary; science-ready cuts §5 / conclusions.
- DES Y1 morph release: https://des.ncsa.illinois.edu/releases/y1a1/gold/morphology
- Gatti et al. 2021 (DES Y3 shear) — metacalibration; not a morph replacement ([arXiv:2011.03408](https://arxiv.org/abs/2011.03408)).
- Unterborn & Ryden 2008, ApJ 687, 976 — \(\Delta M_r = 1.27(\log q)^2\).
- Shao et al. 2007, ApJ 659, 1159 — inclination-dependent LFs / optically thick disks.
- Driver et al. 2007, MNRAS 379, 1022 — MGC attenuation–inclination + Tuffs RT.
- Yip et al. 2010, ApJ 709, 780 — inclination-dependent composite spectra.
- Padilla & Strauss 2008, MNRAS 388, 1321 — SDSS shapes + dust selection note.
- Holmberg 1958; Disney et al. 1993 — classical opacity / selection-bias debates.
