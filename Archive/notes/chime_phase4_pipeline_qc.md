# Phase-4 repeater pipeline QC (2026-07-02)

Manual review after `run_all_frbs.py --frb 20190208A 20190417A 20240114A 20251229A` with `--use-localization-host`.

| FRB | Pipeline GALFIT | Use for inclination? | Magnitude source |
|-----|-----------------|----------------------|------------------|
| **20240114A** | OK — seg #1134 @ 0.40″ from CSV | **Yes** | GALFIT m_r≈21.79 (no dedicated lit r in hand; acceptable) |
| **20251229A** | OK — seg #1129 @ 0.35″ from CSV; Phase 2 AstroPath failed | **Yes** (inc/b/a only) | **Literature** r=18.6 (PS1, ATel #17709); pipeline m=11.12 is uncalibrated |
| **20190417A** | **FIXED** (host override) — seg #1457 @ 0.03″ | **Yes** (after override rerun) | GALFIT m_r≈21.24; lit r=21.47 (Gemini) |
| **20190208A** | **REJECT** — true host absent from cutout | **No** | **Literature** r=27.32±0.16 (GTC O4, Hewitt et al. 2024) |

---

## FRB 20190417A — wrong host / artifact (user flag confirmed)

**Symptoms:** `galfit_results.png` shows saturated ring-like artifacts; χ²/ν = **28.2**; `MAG_AUTO = −6` (non-physical); `SNR_WIN ≈ 0.006`; statmorph flag = 2.

**What happened:**
- CSV target = EVN burst position (294.774566, +59.326897).
- Phase 3a picked **SExtractor #1442** via seg-map at target pixel, but #1442 centroid is **3.69″** from target and **4.17″** from the literature host.
- QA plot: yellow **Target** (+) sits near empty sky; green **Host** circle is on a bright extended artifact to the west with most pixels masked as bad/fringe.

**Literature host (correct):** PRS **20190417A-S1** / optical dwarf (Ibik et al. 2024; PRECISE 2025)
- **Coords:** 19:39:05.82, +59:19:36.7 (J2000) → 294.77425, +59.32686
- **Offset from EVN burst:** ~0.6″ (same galaxy)
- **z** = 0.12817; **m_r** = 21.47 (Gemini); R50 ≈ 0.3″
- **M\*** ≈ 7.6×10⁷ M⊙ (Prospector, Moroianu et al. 2026); **SFR** ≈ 0.19 M⊙ yr⁻¹; **12+log(O/H)** ≈ 8.4 (~0.2 Z⊙)
- **PRS:** 190±40 μJy @ 1.4 GHz, <23 pc; variable RM +3958…+5061 rad m⁻²

**Action:** Fixed via CHIME-only `host_overrides.json` + SPREAD patch + rerun (`rerun_20190417A.py`). Default run picked artifact #1442; override targets psf #1487 → image.cat #1457 (mr≈21.45 in Phase 2 photometry). GALFIT after fix: χ²/ν≈0.99, m_r≈21.24, inc≈67°.

---

## FRB 20190417A — host override (2026-07-02)

**Root cause:** Ibik host centroid falls on empty seg pixel; EVN burst pixel lands on saturated artifact #1442. True dwarf host is **image.cat #1457** (Phase-2 psf #1487, m_r≈21.45) but fails the SPREAD galaxy cut (compact R50≈0.3″).

**Fix (CHIME-only, no changes to `pipeline_scripts/`):**
1. `CHIME/pipeline_scripts/host_overrides.json` — pipeline centre at psf #1487 coords; patch SPREAD for #1487
2. `python CHIME/pipeline_scripts/rerun_20190417A.py` — patch workdir + `master_run --rerun-phase 3a`
3. `run_all_frbs.py` auto-applies override for this FRB after Phase 2

**Result:** host #1457, χ²/ν=0.985, m=21.24, inc=67.1°, b/a=0.43.

---

## FRB 20190208A — associated host **missing** from Legacy cutout

**Symptoms:** Phase 3a: *"no galaxy within 5.0″ — falling back to AstroPath host"* → fit spurious m_r≈20.1 object (P(O)=0.32, 8.7″ from literature position).

**Why:** True host **GTC source O4** is **r = 27.32 ± 0.16** (Hewitt et al. 2024 ApJL 977 L4) — ~8 mag fainter than Legacy 5σ depth (m_lim ≈ 19.4). The host is **not detected** in `large_cutouts/` Legacy DR10 imaging.

**Literature host properties (O4):**
| Property | Value | Source |
|----------|-------|--------|
| Name | GTC source O4 | Hewitt et al. 2024 Table 2 |
| **r (AB)** | **27.32 ± 0.16** (27.17±0.16 in table; abstract uses 27.32) | GTC/OSIRIS |
| Offset from EVN burst | **0.10″** | same |
| PATH P(O\|x) | **0.9995** | Hewitt et al. 2024 |
| z | **not obtained** (too faint for GMOS spectroscopy) | same |
| z_max (DM) | ~0.83 conservative | DM 580 pc cm⁻³ |
| Luminosity | ≲ 10⁸ L⊙ (at z_max); possibly ~10⁶.8 L⊙ at z≈0.19 if O8 interloper | Figure 5 |
| Morphology | extremely faint dwarf; lowest-luminosity FRB host known | same |
| RM | few×10¹ rad m⁻² (modest) | same |
| PRS | none compact (EVN+VLA) | same |
| EVN position | 18:54:11.27, +46:55:21.67 (±260 mas, 2σ) | PRECISE |

**Action:** Exclude from GALFIT/inclination sample. Flag in `chime_host_magnitudes.csv`. Requires GTC/deep imaging cutout for any morphological analysis.

---

## FRB 20251229A — GALFIT OK; photometry failed

- Phase 2 AstroPath WSL bridge failed → **no ZP** → pipeline **m = 11.12 is not trustworthy**.
- Use **literature:** SDSS J204123.23+160126.54, **r = 18.6** (PS1-WISE, ATel #17709); **z = 0.1275 ± 0.0002** (GTC, ATel #17856); M* ~ 10¹⁰ M⊙; SFR 0.7±0.2 M⊙ yr⁻¹; 12+log(O/H) ~ 8.4.
- GALFIT inc ≈ 63.6°, b/a ≈ 0.48 at seg #1129 (0.35″ from CSV host) — **usable for inclination** if combined with literature mag.

---

## FRB 20240114A — acceptable

- Seg #1134 @ 0.40″ from CSV; χ²/ν = 0.76; inc ≈ 55.5°, b/a ≈ 0.59.
- No change needed.
