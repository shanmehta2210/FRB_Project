# Kawinwanichakij+2021 HSC Lenstronomy vs our GALFIT pipeline

> **Status (2026-09-04):** Kawinwanichakij et al. use Lenstronomy, not
> GALFIT. Both analyses use PSF-convolved single-Sérsic models, so the
> comparison is one of model class rather than identical software. Some
> cohort counts and the containment-based neighbor description later in this
> historical audit predate the current 64-host, Re-separation pipeline. For
> paper-facing definitions use `pipeline_scripts/verification/SCIENCE_CUT_AND_COHORT.md`,
> `pipeline_scripts/README.md`, and
> `plots/plots_null/v2/frb_vs_hsc_confirmed50/README.md`.

Comparison for using HSC as an inclination **null**.  
HSC: Kawinwanichakij et al. 2021, ApJ 921, 38 ([arXiv](https://arxiv.org/abs/2109.09778)); catalog from IPMU.  
Ours: `pipeline_scripts/` → `pipeline_galfit_results.csv`.

Skew flags (impact on **axis ratio / cos(i)** unless noted):

| Flag | Meaning |
|------|---------|
| **LOW** | Unlikely to move median cos(i) by ≳0.02 |
| **MED** | Can matter at the few-percent level or in tails |
| **HIGH** | Can systematically shift the cos(i) CDF / selection |

---

## Executive summary

| Area | Similarity | Skew for cos(i) null |
|------|------------|---------------------|
| Single free-n, PSF-convolved Sérsic model | **High** | — |
| Local empirical PSF | **High** (both local) | LOW |
| Neighbor mask + simultaneous fit | **High** (same idea, different rules) | MED |
| Free sky | **Partial** (they free gradients; we freeze) | MED |
| Band (HSC *i* vs our *r*) | Different | **LOW–MED** for *q* at fixed z; see §8 |
| Pre-fit QC / depth | HSC much deeper / stricter | Selection, not shape algorithm |
| Morphological *k*-correction | They correct **Re**; we do none | **LOW for *q***; HIGH for size science |
| *n* ceiling | 8 vs 6 | LOW for disk null |
| Survey / seeing / pixel scale | Different | MED (resolution → roundness at faint end) |

Bottom line: the two pipelines use the **same class of fit** (single Sérsic,
free *n*, free *q*, local PSF, neighbor deblend), but different optimizers.
Main null risks are **survey resolution + faint-end rounding**,
**sky-gradient policy**, and **how we select disk analogues** (*n* window),
not the *i* vs *r* band choice for nearby hosts.

---

## 1. Side-by-side methodology

### 1.1 Data / imaging

| Item | Kawinwanichakij HSC | Our pipeline |
|------|--------------------|--------------|
| Survey | HSC-SSP (Wide + Deep/UltraDeep) | Legacy → PS1 → DES ladder |
| Fit band | **i** (median seeing ~0.6″) | **r** (seeing typically ~0.9–1.2″ on LS/PS1) |
| Pixel scale | HSC ~0.168″/px | Standardized **0.262″/px** |
| Typical depth | Fits to *i*_cmodel ≲ 24.5 | Hosts: median *m_r* ≈ 20; field depth varies |

**Skew: MED (resolution).** Poorer seeing + coarser pixels tend to round faint/compact galaxies (inflate *b/a*). HSC is sharper → better for faint disks. At FRB host brightness (≲21–22), both are usually well-resolved; effect grows toward the faint end.

### 1.2 Pre-fit selection / QC

| Cut | HSC | Ours |
|-----|-----|------|
| Extended | `i_extendedness_value = 1` | SPREAD_MODEL star cut; host must be galaxy |
| Brightness | `i_cmodel_mag < 24.5`, err ≤ 0.1 | No hard mag ceiling for fit; mag seed in [8, 40] |
| Depth / exposures | `inputcount` ≥4 (*g/r*), ≥6 (*i/z/y*) | Whatever the cutout survey provides |
| PSF vs cmodel | `PSFMAG − CMODELMAG > 0.2` | Implicit via SPREAD / CLASS_STAR |
| Sample role | Millions of field galaxies | ~70 FRB hosts |

**Skew: HIGH for sample definition, LOW for fit physics.** HSC’s cuts define a deep, clean parent; ours define a small host list. For null CDFs we already subsample (EXP-analogue *n*, mag limits). Do **not** treat HSC’s 24.5 limit as matching an FRB cut of 21 unless you apply the same mag cut on both.

### 1.3 Postage stamps & neighbors

| Item | HSC | Ours |
|------|-----|------|
| Stamp | Cutout around target | Host seg bbox + 15 px; ≤512 px side |
| Faint / far neighbors | **Masked** (weight → 0) | Mask if containment frac &lt; 0.88 |
| Bright / overlapping | **Simultaneous Sérsic** | Fit if frac ≥ 0.995 (or expand ROI if ≥ 0.88) |
| Stars | Masked (not fit) | Masked (`CLASS_STAR` / SPREAD) |
| Caps | (paper: simultaneous when needed) | ≤25 components; extended hosts → host-only |

**Skew: MED.** Same design philosophy. Differences in thresholds can change *b/a* when a bright companion sits on the major axis (under-masking → artificial elongation; over-masking → rounder / wrong *Re*). Production: 26/70 hosts multi-component. For null galaxies HSC already applied their policy; we cannot re-fit them.

**To look more like HSC:** document containment thresholds against their “close/bright → simultaneous” rule; audit FRB hosts with neighbors for *b/a* sensitivity (mask-all vs fit).

### 1.4 PSF

| Item | HSC | Ours |
|------|-----|------|
| Type | Local empirical HSC pipeline PSF at target coords | Local **PSFEx** from the 10′ field (`proto_image.fits`, 25×25) |
| Spatial variation | Per-object HSC model | PSFVAR_DEGREES = 0 (constant across stamp) |

**Skew: LOW–MED.** Both convolve with a local empirical PSF — correct approach. Our PSF is constant on the stamp (usually fine for ≤512 px). Wrong/too-broad PSF → rounder *b/a* (**MED** at low S/N).

**To look more like HSC:** keep local PSF (already); optional: allow mild PSF variation if stamps grow.

### 1.5 Sky background

| Item | HSC | Ours |
|------|-----|------|
| Sky amplitude | Free in GALFIT | Free (seed from SExtractor BACKGROUND) |
| Sky gradients *dx, dy* | **Free** | **Fixed to 0** |
| Fallback | Aggressive mask + retry if gradient fit fails | If \|sky_fit − sky_ref\| &gt; 3 ADU → retry with ±3 ADU soft constraint |

**Skew: MED for *Re* / outer isophotes; LOW–MED for *q*.** Sky errors mainly bias size and *n*; *b/a* can shift if residual gradients align with the major axis. HSC’s free gradients are safer on deep mosaics with residual sky / ICL.

**Highest-value change toward HSC:** free `dsky/dx`, `dsky/dy` (with a fallback that freezes them if the fit blows up) — already their design.

### 1.6 GALFIT model & constraints

| Parameter | HSC | Ours |
|-----------|-----|------|
| Model | Single 2D Sérsic (target) | Single 2D Sérsic (host = comp 1) |
| Free | *x, y, mag, Re, n, q, PA* | Same |
| *n* range | **0.5 – 8** | **0.5 – 6** |
| *Re* | Free (physical priors in paper QC) | **1.5 – 100** px (≤300 if extended) |
| *q = b/a* | Free | Free (no hard bounds in constraints) |
| Mag | Free | Free, constrained ~8–40 |
| Multi-component science model | No bulge+disk for science table | No bulge+disk; neighbors only |

**Skew: LOW for disk null.** Ceiling *n*=6 vs 8 rarely matters once we keep EXP-analogues (`0.4 < n < 1.5`). Hitting *n*=0.5 or *n*=6 floor/ceiling can indicate failure (**MED** for those objects — flag them).

**To look more like HSC:** raise *n* max to 8 (cosmetic for disks); keep free *q*.

### 1.7 Wavelength / *k*-corrections

| Item | HSC | Ours |
|------|-----|------|
| Fit band | Observed *i* | Observed *r* |
| Size correction | Empirical shift of *Re* → rest-frame **5000 Å** | None |
| Axis ratio *q* | Catalog science uses **fitted** *q* (observer-frame *i*); we use `fitted_q` → `ba` | Observer-frame *r* *b/a* |

**Skew for inclination: LOW–MED.** Kawinwanichakij’s *k*-correction targets **size–mass**, not *q*. Axis ratio is only weakly band-dependent between adjacent optical bands for the same galaxy (color gradients change *Re* and *n* more than projected *q*; see GAMA MegaMorph / Vulcani-type results). Dust can make optical *q* slightly wavelength-dependent in edge-ons (**MED** for dusty disks).

---

## 2. Does *i*-band (HSC) vs *r*-band (us) matter?

**For cos(i) from *b/a*: usually small — not a showstopper.**

Reasons it is OK for a null:

1. *r* and *i* are adjacent; for typical FRB hosts (*z* ≲ 0.5) both probe similar stellar light.
2. Inclination uses **projected axis ratio**, which is far less wavelength-sensitive than *Re* or *n*.
3. HSC’s published *k*-correction is for **Re**, not for *q*; our null already uses `fitted_q` / `ba` in the observer-frame fit band.

When it is *not* negligible (**MED**):

- Strong color gradients + dusty edge-ons: bluer bands can look thicker / rounder.
- High-*z* HSC galaxies where observed *i* ≈ rest-UV vs our local *r* ≈ rest-*r* — population mismatch, not a per-object band bug. Mitigate with photo-*z* / mass cuts on the null if needed.

**Practical rule:** matching **mag limit + disk selection (*n*) + *q* &gt; *q*_0** matters more than forcing the same filter name. Prefer comparing at the same **observer-frame** depth in each survey’s fit band (*i* for HSC, *r* for FRBs), or convert FRB *r* ↔ HSC *i* with a mean color if you need a joint cut.

---

## 3. Depth: can we use FRBs to 21.5 or 22?

**Yes for the HSC null — HSC is deeper than that.**

- HSC fit selection goes to *i* ≈ **24.5**; our EXP-analogue HSC table still has tens of thousands with `fitted_mag` / `imag` ≤ 22 (median `fitted_mag` ≈ 22.8 in the EXP-analogue file).
- FRB `mag_final` (*r*, *N*=67 finite), with `b_a > 0.2`:

| Mag cut | *N* FRB hosts |
|--------:|-------------:|
| ≤ 21.0 | 42 |
| ≤ 21.5 | 44 |
| ≤ 22.0 | **52** |
| ≤ 22.5 | 53 |

So moving the joint cut from 21 → **21.5 or 22** is limited by **FRB sample size and our fit quality at the faint end**, not by HSC catalog depth.

**Caveats (MED–HIGH for our fits, not HSC):**

- Our imaging (LS/PS1) is shallower and softer-seeing than HSC; *b/a* errors and PSF-rounding grow faintward (see `tools/psf_ba_mag_sim`).
- Prefer mag cuts where FRB GALFIT χ² / constraint hits remain clean; audit *n* stuck at 0.5/6 and extreme *q*.

---

## 4. Difference register (what to change)

### Already similar (keep)

- Single free-*n* Sérsic for the science object  
- Free *q*, *PA*, *Re*, mag, centroid  
- Local empirical PSF convolution  
- Hybrid neighbor **mask + simultaneous Sérsic**  
- Free sky pedestal  

### Decisions (2026-08)

- **Sky gradients stay fixed at 0.** Shallower LS/PS1/DES cutouts do not support a robust free sky plane the way deep HSC mosaics do; free gradients risk absorbing galaxy light. Amplitude-only free sky (+ existing ±3 ADU QA) is the right choice for us.
- **Neighbor policy → move toward HSC.** Prefer their split: faint/far → **mask**; bright/close overlapping → **simultaneous Sérsic**. Planned change to our containment rules in `generate_galfit_cutouts.py` (not yet implemented).

### Change to move **toward** HSC (priority order)

| Change | Effort | Expected skew reduction | Status |
|--------|--------|-------------------------|--------|
| 1. Align neighbor policy with “far/faint → mask; close/bright → simultaneous” | Med | MED in crowded fields | **Planned** |
| 2. Faint-end QA: flag/reject PSF-dominated *q* (compare to sims) | Med | HIGH for mag≳21.5 | Optional |
| 3. Raise *n* max 6 → 8 | Trivial | LOW | Optional |
| 4. Optional: report *i*-band or *r*−*i* for hosts when multi-band exists | Med | LOW for *q* | Optional |
| ~~Free sky gradients~~ | — | — | **Rejected** (depth) |

### Do **not** need for cos(i) null parity

- Morphological *k*-correction of *Re* to 5000 Å (size science only)  
- Identical pre-fit `inputcount` / cmodel&lt;24.5 (apply **matched mag cuts** instead)  
- Fitting on HSC *i* cutouts for every FRB (nice-to-have, not required if *q* is stable)  
- Free sky *dx/dy* (see decision above)  

### Null-catalog side (already in repo)

- EXP analogue: `0.4 < fitted_sersic < 1.5` + `goodfits_flag=1` → `catalog/HSC_kawinwanichakij_exp_analogue.csv`  
- Use `ba` (= `fitted_q`), Hubble *q*_0 = 0.2, same as FRB strict pool  

---

## 5. Residual systematics checklist (when claiming HSC = null)

1. **Resolution mismatch** — HSC sharper than LS/PS1; expect HSC less face-on-biased at faint mag.  
2. **Band + dust** — small *r*/*i* *q* difference possible for edge-on dusty disks.  
3. **Selection** — HSC photo-*z* / mass / quiescent flags exist; FRB hosts are a special population — match late-type via *n* (and optionally color), not raw mag alone.  
4. **Neighbor policy** — irreducible difference unless we re-reduce HSC; treat as random error unless crowded hosts dominate.  
5. **Sky gradients** — intentionally left fixed (shallower data).  
6. **Neighbor policy** — main remaining algorithmic parity lever on our side.

---

## 6. Key paths

| Role | Path |
|------|------|
| Our GALFIT runner | `pipeline_scripts/galfit_fitting/run_galfit_fitting.py` |
| Neighbor cutouts | `pipeline_scripts/galfit_fitting/generate_galfit_cutouts.py` |
| PSFEx | `pipeline_scripts/SExtractor + PSFEx/run_psf_pipeline.py` |
| Pipeline README | `pipeline_scripts/README.md` |
| FRB results | `pipeline_galfit_results.csv` |
| HSC sample builder | `scripts/build_hsc_kawinwanichakij_sample.py` |
| HSC EXP analogue | `scripts/extract_exp_analogue_des_hsc.py` → `catalog/HSC_kawinwanichakij_exp_analogue.csv` |
| Paper | Kawinwanichakij et al. 2021, ApJ 921, 38 |

---

## 7. One-line verdict

Pipelines are **the same species of fit** (local-PSF single-Sérsic GALFIT with neighbor deblend). For an inclination null, **i vs r is a second-order issue**; prioritize **matched mag cuts**, **disk *n* selection**, **faint-end PSF bias control**, and (on our side) **free sky gradients**. HSC depth comfortably supports FRB cuts at **21.5–22**.
