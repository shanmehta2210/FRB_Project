# Plan: inclination-dependent dust / extinction correction

**Result summary + chosen baseline:** [`../DUST_AND_MEDIAN_BA.md`](../DUST_AND_MEDIAN_BA.md).

Goal: correct (or control for) the bias that **edge-on disks are dimmed by internal dust**, so they preferentially drop out of magnitude-limited samples and the surviving cos(i) / b/a distributions look too face-on.

This plan is for **selection / luminosity correction**, not for morph-fit calibration (Tarsitano η, HSC `corrected_*`). Those fix measurement bias; dust fixes **who enters the sample**.

Scope for our stacks: LS EXP (optionally scaled), DES Y1 morph EXP-analogue, HSC Kawinwanichakij EXP-analogue; mag cuts ≤20/21/22; Hubble cos(i) with q₀=0.2.

---

## 1. Problem statement (what we are correcting)

For a dusty thin disk, line-of-sight optical depth rises toward edge-on:

- Observed magnitude: \(m_{\rm obs} = m_{\rm face} + A(q)\) with \(A(q)\ge 0\), \(A\to 0\) as \(q\to 1\).
- In a cut \(m_{\rm obs} \le m_{\lim}\), edge-ons need to be intrinsically brighter to pass → under-represented at fixed depth.
- Observed median cos(i) rises (often ~0.55–0.65 for late-types). That can be **expected**, not a bug in Hubble’s formula.

**What dust does *not* explain:** LS Tractor `type=EXP` edge-on bias from the REX gate (face-ons → REX). That is a **pipeline selection** effect; keep LS ba≤0.8 scaling separate from dust.

**Circular dependency:** \(q\) (or cos i) is both (a) the science variable and (b) the argument of \(A(q)\). Any correction that uses \(q\) to brighten edge-ons will change who sits in the mag cut and therefore the cos(i) CDF. That is intended for selection bias, but it must be documented and sensitivity-tested.

---

## 2. Competing methods (document all; pick later)

### Method A — Empirical face-on magnitude + re-cut (Unterborn & Ryden 2008)

**Law (r-band, bright exponential disks, SDSS):**

\[
\Delta M_r = 1.27\,(\log_{10} q)^2,\qquad
M_r^{f} = M_r - \Delta M_r
\]

(equivalently for apparent mag in the same band: \(m^{f} = m - \Delta m(q)\)).

Linear \(\Delta M \propto \log q\) was a worse fit; near face-on (\(q\gtrsim 0.5\)) dimming is weak.

**How to apply to our catalogs (no redshift needed):**

1. Restrict to disk-like pool (EXP / EXP-analogue; already done).
2. Compute \(\Delta m_r(q) = 1.27\,(\log_{10} q)^2\) with \(q=b/a\).
3. \(m^{f} = m_r - \Delta m_r\).
4. Rebuild pools with \(m^{f} \le 20,21,22\) (instead of \(m\le\ldots\)).
5. Recompute cos(i) CDFs / medians / KS vs uniform.

**Pros:** Simple; designed exactly for mag-limited shape samples; published coefficient for exponential disks in *r*.  
**Cons:** Calibrated on bright (\(M_r\lesssim -19\)) SDSS exponentials; may over-correct dwarfs / faint high-z; uses observed \(q\) (noise + thickness).  
**Best first method for us.**

Variants to A/B against the same pipeline:

| Variant | Formula | Source |
|---------|---------|--------|
| A1 (baseline) | \(1.27(\log q)^2\) | Unterborn & Ryden 2008 |
| A2 | \(\gamma_r \log(a/b)\) with \(\gamma_r\sim 0.5\)–\(1.4\) | TF / Shao-style linear laws |
| A3 | Masters-like quadratic \(c_1\log q + c_2(\log q)^2\) | Masters et al. 2003 (NIR); UR08 also quote a mixed fit |
| A4 | Cap: \(\Delta m=0\) for \(q>0.5\) (UR08 qualitative) | optional |

---

### Method B — Classic TF / catalog extinction \(A=\gamma\log(a/b)\)

\[
A_\lambda = \gamma_\lambda\,\log_{10}(a/b) = -\gamma_\lambda\,\log_{10} q
\]

Tully et al. (1998): \(\gamma\) **luminosity-dependent** (brighter disks dustier); often proxied by HI linewidth \(W\) when \(M\) is circular.

Typical optical \(\gamma_r\sim 1\)–\(1.4\) for luminous spirals; near zero for faint dwarfs. Edge-on vs face-on relative attenuation \(\sim 0.85\gamma\) over the observed \(q\) range (Devour & Bell 2016).

**Pros:** Ubiquitous; easy.  
**Cons:** Linear law systematically mis-fits vs \((\log q)^2\) (UR08; Masters NIR); luminosity dependence needs a mass/luminosity proxy we may lack for pure photometry samples.  
**Status:** Implement as sensitivity suite (A2), not sole baseline.

---

### Method C — Optically thick slab / LF joint fit (Shao et al. 2007)

Model: optical depth \(\tau \propto \cos i\) (or related); likelihood recovers both face-on LF and amplitude \(\gamma\) from inclination-binned LFs. Finds optically thick disks; \(M^*\) dims by ~1.2 mag (*u*) to ~0.5 mag (*z*) vs face-on.

**Pros:** Proper statistical control of selection in a flux+redshift sample; band-dependent extinction curve (\(\tau\propto\lambda^{-n}\), \(n\approx 0.96\)).  
**Cons:** Needs redshift + volume machinery (NYU-VAGC style); not drop-in for LS/DES photometric mag cuts alone.  
**Status:** **Phase-2** if we build spectroscopic overlap subsets (SDSS/DES/HSC with z). Use as physics prior / sanity check, not first code path.

---

### Method D — Radiative-transfer attenuation–inclination (Driver / Tuffs / Popescu)

Empirical \(\Delta M\) vs \(1-\cos i\), then Tuffs et al. RT models to also remove residual **face-on** opacity (Driver et al. 2007: central \(\tau_B^f\sim 3.8\)). Separate **disc vs bulge** attenuation (bulge light seen through dusty disc).

**Pros:** Physically motivated; corrects absolute (not only relative-to-face-on) attenuation.  
**Cons:** Needs B/T or bulge–disc decompositions; heavy assumptions; overkill for null cos(i) CDFs.  
**Status:** Optional later if we care about absolute luminosities / SFR; **not** required for first-order selection bias on cos(i).

---

### Method E — Spectral / continuum empirical laws (Yip et al. 2010)

Composite spectra vs axis ratio; e.g. \(\eta_g\sim 1.2\) mag continuum extinction for highly inclined (\(b/a\sim 0.1\)) vs face-on in a volume-limited SDSS box.

**Pros:** Direct spectral constraint.  
**Cons:** Central fiber aperture; volume-limited calibration sample; convert \(\eta(q)\) to a full \(A(q)\) curve for mag re-cuts needs an assumed functional form.  
**Status:** Use to **cross-check** amplitude of A1/A2 at the edge-on end, not as primary re-cut law.

---

### Method F — Bypass dust with selection design (no \(\Delta m\) formula)

| Approach | Idea | Fit for us? |
|----------|------|-------------|
| F1 Volume-limited | Absolute mag + z box so edge-ons cannot fall out (Yip; many SDSS shape papers) | Only where reliable photo-z / spec-z exist (HSC photoz; DES/LS harder) |
| F2 NIR / MIR selection | Select on WISE/2MASS (Devour & Bell 2016); optical used only for ba | Strong methodologically; needs cross-match |
| F3 Inclination-binned completeness model | Forward-model \(P({\rm detect}|q,m)\) | Most rigorous; most work |

Devour & Bell stress: **measuring** \(\gamma\) from optically selected samples is itself biased unless selection is inclination-independent (they use WISE + per-bin redshift caps). Lesson for us: prefer applying a **published** \(A(q)\) to re-cut, or NIR-selected subsets, over refitting \(\gamma\) on our optical mag-limited pools.

---

## 3. Recommended staged plan

### Stage 0 — Preconditions (do not skip)

1. Science-ready morph pools only (DES: `SN>30`, mag≲21.5 / `FIT_STATUS` when ingested; HSC: `goodfits`; LS: EXP + optional scaled ba≤0.8 for REX, **separate** from dust).
2. Freeze EXP-analogue definition (`0.4<n<1.5` DES/HSC).
3. Record **pre-correction** CDF medians (already in `v2/comparisons/summary.csv`).
4. Never mix “dust correction” with “LS ba scaling” in the same sentence in plots — two different physics knobs.

### Stage 1 — Primary implementation (Method A1)

**Deliverable:** script e.g. `scripts/apply_inclination_extinction.py`

For each survey CSV:

```text
q      = ba
dm     = 1.27 * (log10(q))**2          # A1; clamp q in (0,1]
m_face = m_r - dm
keep   = (q > q0) & (m_face <= m_lim)  # DES/HSC
# LS: also optional (q <= 0.8) for REX scaling AFTER or BEFORE dust — test both orders
cosi   = Hubble(q, q0)
```

Outputs under `plots/plots_null/v2/extinction/`:

| Artifact | Content |
|----------|---------|
| `summary_before_after.csv` | N, median cos(i), median ba per survey × mag × {raw, A1} |
| `cdfs/mag{20,21,22}_*.png` | CDFs before/after overlay per survey |
| `compare/mag{20,21,22}.png` | LS / DES / HSC after A1 |
| `delta_m_vs_ba.png` | \(\Delta m(q)\) curve + histogram of applied \(\Delta m\) |
| `funnel.csv` | how many objects enter/leave each mag cut after correction |

**Success metrics (not pass/fail — report):**

- Change in median cos(i) (expect **decrease** toward ~0.5 if dust selection dominated face-on bias).
- Change in N at each mag cut (expect **more** edge-ons retained → larger N or redistribution).
- DES/HSC should move closer to each other if dust was a shared bias; LS may still differ due to REX.

### Stage 2 — Competing-law suite (same pipeline, Methods A2–A4 + B)

Grid of \(\Delta m(q)\) laws; same re-cut machinery. Table of median cos(i) vs law. Identify which laws move medians by ≲0.02 vs ≳0.05 (robust vs sensitive).

Suggested coefficient set for A2: \(\gamma_r \in \{0.5, 1.0, 1.27, 1.4\}\).

### Stage 3 — Amplitude sanity checks (Methods E, literature)

- At \(q=0.2\): A1 gives \(\Delta m = 1.27(\log 0.2)^2 \approx 0.89\) mag.  
- At \(q=0.1\): A1 ≈ 1.27 mag — comparable to Yip \(\eta_g\sim 1.2\) and Shao edge-on vs face-on ~0.65 in *r* for LF \(M^*\) (different statistic).  
- Flag if any law implies \(\Delta m > 2\) mag inside our \(q\) range (unphysical for this application).

### Stage 4 — Optional / heavier (Methods C, D, F)

Only if Stage 1–2 show dust is important **and** residual DES–HSC–LS tension remains:

1. **F1:** HSC photo-z volume box + absolute mag cut; repeat CDFs.  
2. **F2:** WISE cross-match for LS/DES subset; select on \(W1\), measure optical cos(i).  
3. **C:** Spectroscopic SDSS overlap for a Shao-style check.  
4. **D:** Skip unless we need dust-free absolute magnitudes.

### Stage 5 — Report integration

Append results to `DES_calibration_and_dust_research.md` §3 and a short `v2/extinction/README.md`. Keep dust results visually separate from LS scaled CDFs in `v2/comparisons/`.

---

## 4. Survey-specific notes

| Survey | Mag column | ba | Dust method notes |
|--------|------------|-----|-------------------|
| **DES** | `mag_r` = calibrated `MAG_SERSIC_R` | `1−ε` calibrated | Best optical target for A1; apply only in science-ready SN/mag regime |
| **HSC** | prefer `corrected_mag` when `calib_flag=1` for Stage 1+ (morph bias ≠ dust) | `fitted_q` (floor 0.1) | No `corrected_q`; ba floor already truncates extreme edge-ons |
| **LS** | `modelMag_r` | from e1,e2 | Apply A1 on full EXP; report **with and without** ba≤0.8 REX scaling as a 2×2 matrix (dust × REX) |

**Order of operations (LS):** preferred default = morph/REX cuts first, then dust re-cut on remaining objects. Also run dust-first then REX cap as sensitivity.

**Galactic extinction:** catalogs are generally already MW-extinction corrected (DES GOLD / HSC cmodel / LS). Do **not** double-apply Schlegel maps unless a column is explicitly uncorrected.

---

## 5. What we will *not* do in Stage 1

- Refit \(\gamma\) from our mag-limited optical samples (biased estimator — Devour & Bell).
- Apply Tuffs RT / bulge–disc without decompositions.
- Treat LS REX scaling as a dust correction.
- Claim “dust-free intrinsic cos(i)” — A1 only removes **differential** edge-on vs face-on selection to first order.
- Use uncalibrated DES ellipticity or double-add `MAG_CAL`.

---

## 6. Concrete work packages

| # | Task | Depends on | Output |
|---|------|------------|--------|
| WP1 | Implement `delta_m(q, law=...)` + `m_face` helpers in `null_catalog_utils` | — | shared API |
| WP2 | A1 re-cut CDFs for DES, HSC, LS at mag 20/21/22 | WP1 | `v2/extinction/` plots + summary |
| WP3 | Law grid A2–A4 | WP2 | sensitivity table |
| WP4 | LS 2×2 (dust × ba≤0.8) | WP2 | matrix of medians |
| WP5 | Decision memo: is dust enough to explain DES/HSC face-on medians vs LS? | WP2–4 | short md verdict |
| WP6 (optional) | HSC volume-limited F1 or WISE F2 | WP5 | follow-up |

---

## 7. Decision criteria after Stage 1–2

| Observation | Interpretation | Next |
|-------------|----------------|------|
| After A1, DES & HSC median cos(i) → ~0.5±0.05 | Dust selection was a major driver of face-on bias | Prefer A1 (or nearest grid law) for science CDFs |
| Medians barely move (≲0.02) | Dust selection sub-dominant at these depths / pools | Keep raw mag cuts; document null result |
| DES moves, HSC does not (or vice versa) | Depth / n-mix / fit floor differences | Per-survey laws or F1 volume cut |
| LS still far from DES/HSC after dust + REX scaling | Pipeline / population mismatch dominates | Do not force dust to reconcile surveys |

---

## 8. Key references

- Unterborn & Ryden 2008, ApJ 687, 976 — \(\Delta M_r=1.27(\log q)^2\); face-on-corrected shape sample.
- Tully et al. 1998, AJ 115, 2264 — \(A=\gamma\log(a/b)\); luminosity-dependent \(\gamma\).
- Shao et al. 2007, ApJ 659, 1159 — inclination-binned LFs; optically thick disks.
- Masters et al. 2003, 2010 — NIR / Galaxy Zoo; quadratic / bilinear laws; luminosity trends.
- Driver et al. 2007, MNRAS 379, 1022 — attenuation–inclination + Tuffs RT; bulge vs disc.
- Yip et al. 2010, ApJ 709, 780 — spectral extinction vs axis ratio.
- Devour & Bell 2016, MNRAS 459, 2054 — WISE-selected \(\gamma(M,{\rm sSFR})\); **selection bias warning**.
- Padilla & Strauss 2008 — dust under-represents edge-ons in mag-limited shape samples.
- Disney et al. 1993 — classical opacity / selection-effect caution.

---

## 9. Bottom line

**First implementation:** Method **A1** (Unterborn & Ryden) — subtract \(1.27(\log q)^2\) from *r*-band mag, re-apply mag cuts, redo cos(i) CDFs for DES/HSC/LS.

**In parallel as competitors (same code path):** linear \(\gamma\log(a/b)\) grid and optional \(q>0.5\) zero-dimming cap.

**Defer:** Shao LF likelihood, Tuffs RT, refitting \(\gamma\) on our optical samples.

**Keep separate:** LS ba≤0.8 / cos(i) scaling (REX), DES/HSC morph calibration, and this dust selection correction.
