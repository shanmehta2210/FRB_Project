# Root `repeater_localizations.csv` — catalog + citation guide

Centralized catalog of **every localized repeating FRB host** recognized in this repository (20 rows).

**Regenerate CSV:**
```bash
python scripts/build_root_repeater_localizations.py
```

This README is the **citation map** for values in that CSV. Use the papers listed under each FRB when writing the manuscript. Do not cite this README as a primary source — cite the papers below.

---

## Contents

| Group | Count | Origin |
|-------|------:|--------|
| CHIME-discovered | 16 | `CHIME/repeater_localizations.csv` |
| Non-CHIME discovery | 4 | `master_frb_localization.csv` + literature (`20190520B`) |
| **Total** | **20** | |

Overlaps that appear in both CHIME and master catalogs (`20190711A`, `20201124A`, `20220912A`) are listed **once**; CSV coords/DM follow the CHIME build.

---

## How to read each entry

For every FRB:

| Role | Meaning |
|------|---------|
| **Repeater confirmation** | Paper/telegram that establishes the source repeats |
| **Host / localization** | Paper used for host R.A./Dec. (and burst localization if distinct) |
| **Catalog properties** | Exact papers for *z*, DM, *m*<sub>r</sub> as stored in the CSV |
| **Inclination / *b*/*a*** | This work’s GALFIT (not a literature citation), if present |

DOIs / ATel links below were checked against `CHIME/SOURCES_AUDIT.md`, `CHIME/catalog/host_association_papers.md`, `Archive/notes/repeater_audit_notes.md`, `master_frb_localization.csv`, and (for Phase 4) `CHIME/papers/phase4_repeater_literature.md`.

---

## Non-CHIME discovery (4)

### FRB 20121102A (Arecibo)

| Role | Citation | DOI / ID |
|------|----------|----------|
| **Repeater confirmation** | Spitler et al. 2016, Nature **531**, 202 (first confirmed repeater) | [10.1038/nature17168](https://doi.org/10.1038/nature17168) |
| **Host / localization (VLBI + PRS)** | Chatterjee et al. 2017, Nature **541**, 58 | [10.1038/nature20797](https://doi.org/10.1038/nature20797) |
| **Host *z* / optical ID** | Tendulkar et al. 2017, ApJL **834**, L7 | [10.3847/2041-8213/834/2/L7](https://doi.org/10.3847/2041-8213/834/2/L7) |
| **Catalog *m*<sub>r</sub> = 23.73** | Gordon et al. 2023, ApJ **952**, 122, Table 1 (GMOS-N *r*; cites Chatterjee/Tendulkar) | [10.3847/1538-4357/ace5aa](https://doi.org/10.3847/1538-4357/ace5aa) |
| **Catalog *z*, DM, DM_MW** | From `master_frb_localization.csv` (Chatterjee/Tendulkar lineage; DM_MW NE2001 in master build) | — |
| Inclination | — (no usable GALFIT in this work) | — |

### FRB 20180301A (Parkes)

| Role | Citation | DOI / ID |
|------|----------|----------|
| **Repeater confirmation + host association** | Bhandari et al. 2022, AJ **163**, 69 (marks source as repeater; host properties) | [10.3847/1538-3881/ac3aec](https://doi.org/10.3847/1538-3881/ac3aec) |
| **Catalog *m*<sub>r</sub> = 21.21** | Gordon et al. 2023 Table 1 (NOT *r*; cites Bhandari et al. 2022) | [10.3847/1538-4357/ace5aa](https://doi.org/10.3847/1538-4357/ace5aa) |
| **Catalog coords, *z*, DM, DM_MW** | `master_frb_localization.csv` row (survey RealFAST; values from Bhandari et al. 2022 lineage) | — |
| Inclination | — | — |

### FRB 20190520B (FAST)

| Role | Citation | DOI / ID |
|------|----------|----------|
| **Repeater confirmation + localization + PRS + DM + *z*** | Niu et al. 2022, Nature **606**, 873 | [10.1038/s41586-022-04755-5](https://doi.org/10.1038/s41586-022-04755-5) |
| **Catalog *m*<sub>r</sub> = 22.16** | Gordon et al. 2023 Table 1 (SOAR *r*; cites Niu et al. 2022) | [10.3847/1538-4357/ace5aa](https://doi.org/10.3847/1538-4357/ace5aa) |
| **Catalog DM = 1204.7, DM_MW ≈ 113** | Niu et al. 2022 (DM); MW contribution as used in Niu et al. / follow-up DM budget papers | — |
| **Catalog host coords** | Gordon et al. 2023 Table 1 / Niu et al. 2022 (16:02:04.27, −11:17:17.3) | — |
| Inclination | — (not in master CSV; no GALFIT here) | — |
| Note | Recognized in `CHIME/SOURCES_AUDIT.md` as excluded from the CHIME-only catalog; included in this **root** file | — |

### FRB 20230814B (DSA-110; aka FRB 20230814A)

| Role | Citation | DOI / ID |
|------|----------|----------|
| **Repeater confirmation + localization** | Ravi, ATel **#16191** (2023 Aug 16) — DSA-110 discovery + interferometric localization of the *repeating* source (2 bursts; DSA name FRB 20230814A) | [ATel 16191](https://www.astronomerstelegram.org/?read=16191) |
| **Host *z* = 0.553, DM = 696.40 (refereed)** | Li et al. 2025, ApJ **989**, 77, Table 1 (entry FRB 20230814A, ref. 14) — DM<sub>IGM</sub>–*z* host compilation | [10.3847/1538-4357/adeb72](https://doi.org/10.3847/1538-4357/adeb72) |
| **Catalog coords, DM_MW** | `master_frb_localization.csv` (survey DSA-110; ATel #16191 as `repeater_source`) | — |
| **Catalog *m*<sub>r</sub>** | blank — no published *r* adopted in this catalog | — |
| Inclination | — (`pipeline_galfit_results.csv` lists `missing_fit_log`) | — |
| Caveat | Repeater **status** rests on the DSA-110 telegram (no dedicated refereed repeater paper found as of 2026-07); host *z*/DM **are** in a refereed compilation (Li et al. 2025). Not in Law et al. 2024 — that DSA catalog explicitly excludes repeaters | — |

---

## CHIME-discovered (16)

### FRB 20180814A

| Role | Citation | DOI / ID |
|------|----------|----------|
| **Repeater confirmation** | CHIME/FRB Collaboration 2019, Nature **566**, 235 — *A second source of repeating fast radio bursts* (180814.J0422+73 = the 2nd repeater ever discovered) | [10.1038/s41586-018-0864-x](https://doi.org/10.1038/s41586-018-0864-x) |
| **Host coords, *z*, DM, *m*<sub>r</sub>** | Michilli et al. 2023, ApJ **950**, 134 (`tab:r2galaxies` / host sections; Pan-STARRS *r* = 17.15) | [10.3847/1538-4357/accf89](https://doi.org/10.3847/1538-4357/accf89) |
| Inclination / *b*/*a* | This work — CHIME GALFIT `20180814A_all` | — |

### FRB 20180916B

| Role | Citation | DOI / ID |
|------|----------|----------|
| **Repeater confirmation** | CHIME/FRB Collaboration 2019 (discovery as repeater) | Nature **566**, 235 — [10.1038/s41586-019-0893-0](https://doi.org/10.1038/s41586-019-0893-0) |
| **Host / EVN localization, *z*, DM** | Marcote et al. 2020, Nature **577**, 190 | [10.1038/s41586-020-2300-2](https://doi.org/10.1038/s41586-020-2300-2) |
| **Catalog *m*<sub>r</sub> = 16.17** | Gordon et al. 2023 Table 1 (SDSS *r*; cites Marcote) | [10.3847/1538-4357/ace5aa](https://doi.org/10.3847/1538-4357/ace5aa) |
| Inclination / *b*/*a* | This work — GALFIT `20180916B_all` | — |

### FRB 20181030A

| Role | Citation | DOI / ID |
|------|----------|----------|
| **Repeater confirmation** | CHIME/FRB Collaboration 2019, ApJL **885**, L24 (eight new repeating FRBs; includes 181030.J1054+73 / 20181030A) | [10.3847/2041-8213/ab4a80](https://doi.org/10.3847/2041-8213/ab4a80) |
| **Host (NGC 3252), *z*, DM, *m*<sub>r</sub>** | Bhardwaj et al. 2021, ApJL **919**, L24 | [10.3847/2041-8213/ac223b](https://doi.org/10.3847/2041-8213/ac223b) |
| Inclination / *b*/*a* | This work — GALFIT `20181030A_all` | — |

### FRB 20190110C

| Role | Citation | DOI / ID |
|------|----------|----------|
| **Repeater confirmation + host (*z*, DM, *m*<sub>r</sub>, coords)** | Ibik et al. 2024, ApJ **961**, 99 (RN3-host; PATH-secure) | [10.3847/1538-4357/ad0893](https://doi.org/10.3847/1538-4357/ad0893) |
| Inclination / *b*/*a* | This work — GALFIT `20190110C_all` | — |

### FRB 20190208A

| Role | Citation | DOI / ID |
|------|----------|----------|
| **Repeater confirmation** | Fonseca et al. 2020, ApJL **891**, L6 (source #1 = 190208.J1855+46) | [10.3847/2041-8213/ab7208](https://doi.org/10.3847/2041-8213/ab7208) · [arXiv:2001.03595](https://arxiv.org/abs/2001.03595) |
| **Host / EVN localization + *m*<sub>r</sub> = 27.32** | Hewitt et al. 2024, ApJL **977**, L4 (GTC *r*; PATH 99.95%; DSA name/CHIME source, EVN localization) | [10.3847/2041-8213/ad8ce1](https://doi.org/10.3847/2041-8213/ad8ce1) · [arXiv:2410.17044](https://arxiv.org/abs/2410.17044) |
| **Catalog DM** | Michilli et al. 2023 (consistent with EVN); DM_MW = NE2001 in CHIME build | [10.3847/1538-4357/accf89](https://doi.org/10.3847/1538-4357/accf89) |
| **Catalog *z*** | blank — no spectroscopic *z* (Hewitt et al. 2024) | — |
| Inclination | **Not used** — true host absent from Legacy cutout; pipeline GALFIT is spurious | — |

### FRB 20190303A

| Role | Citation | DOI / ID |
|------|----------|----------|
| **Repeater confirmation** | Fonseca et al. 2020, ApJL **891**, L6 (source 190303.J1353+48, DM = 222.4) | [10.3847/2041-8213/ab7208](https://doi.org/10.3847/2041-8213/ab7208) · [arXiv:2001.03595](https://arxiv.org/abs/2001.03595) |
| **Host coords, *z*, DM (via Outriggers burst 20231204A)** | CHIME/FRB Collaboration 2025, ApJS **280**, 6 (KKO Outriggers catalog) | [10.3847/1538-4365/addbda](https://doi.org/10.3847/1538-4365/addbda) · [arXiv:2502.11217](https://arxiv.org/abs/2502.11217) |
| **Catalog *m*<sub>r</sub> = 15.50** | Michilli et al. 2023 (SDSS Petrosian *r*, brighter merger member) | [10.3847/1538-4357/accf89](https://doi.org/10.3847/1538-4357/accf89) |
| Inclination / *b*/*a* | This work — GALFIT `20190303A_all` | — |

### FRB 20190417A

| Role | Citation | DOI / ID |
|------|----------|----------|
| **Repeater confirmation** | Fonseca et al. 2020, ApJL **891**, L6 (source = 190417.J1939+59) | [10.3847/2041-8213/ab7208](https://doi.org/10.3847/2041-8213/ab7208) · [arXiv:2001.03595](https://arxiv.org/abs/2001.03595) |
| **Host *z*, host centroid, *m*<sub>r</sub> = 21.47 (Gemini)** | Ibik et al. 2024, ApJ **961**, 99 (PRS 20190417A-S1 host) | [10.3847/1538-4357/ad0893](https://doi.org/10.3847/1538-4357/ad0893) |
| **Mas localization + PRS + RM + DM** | Kirsten et al. / CHIME–PRECISE 2025, ApJL | [10.3847/2041-8213/ae28c7](https://doi.org/10.3847/2041-8213/ae28c7) · [arXiv:2509.05174](https://arxiv.org/abs/2509.05174) |
| Inclination / *b*/*a* | This work — GALFIT after CHIME host override (seg #1457) | — |

### FRB 20190711A

| Role | Citation | DOI / ID |
|------|----------|----------|
| **Repeater confirmation (ASKAP repeats)** | Kumar et al. 2021, MNRAS **500**, 2525 | [10.1093/mnras/staa3436](https://doi.org/10.1093/mnras/staa3436) |
| **Host association / *z* (ASKAP/CRAFT)** | Heintz et al. 2020, ApJ **903**, 152 | [10.3847/1538-4357/abb6fb](https://doi.org/10.3847/1538-4357/abb6fb) |
| **Catalog *m*<sub>r</sub> = 23.54, DM cross-check** | Gordon et al. 2023 Table 1 (GMOS-S *r*) | [10.3847/1538-4357/ace5aa](https://doi.org/10.3847/1538-4357/ace5aa) |
| Note | Also detected/listed as a CHIME repeater; localization facility is ASKAP/CRAFT | — |
| Inclination / *b*/*a* | This work — GALFIT `20190711A_all` | — |

### FRB 20191106C

| Role | Citation | DOI / ID |
|------|----------|----------|
| **Repeater confirmation (RN3)** | Ibik et al. 2024 (initially marginal PATH) | [10.3847/1538-4357/ad0893](https://doi.org/10.3847/1538-4357/ad0893) |
| **Secure host upgrade (Outriggers burst 20231128A): coords, *z*, DM** | CHIME/FRB Collaboration 2025, ApJS **280**, 6 | [10.3847/1538-4365/addbda](https://doi.org/10.3847/1538-4365/addbda) |
| **Catalog *m*<sub>r</sub> = 17.306** | Ibik et al. 2024 (DESI *r*) | [10.3847/1538-4357/ad0893](https://doi.org/10.3847/1538-4357/ad0893) |
| Inclination / *b*/*a* | This work — GALFIT `20191106C_all` | — |

### FRB 20200120E

| Role | Citation | DOI / ID |
|------|----------|----------|
| **Repeater confirmation + M81 association** | Bhardwaj et al. 2021, ApJL **911**, L34 | [10.3847/2041-8213/ac0532](https://doi.org/10.3847/2041-8213/ac0532) |
| **EVN localization (GCl 01), DM** | Kirsten et al. 2022, Nature **602**, 585 | [10.1038/s41586-021-04354-w](https://doi.org/10.1038/s41586-021-04354-w) |
| **Catalog *m*<sub>r</sub> = 17.09** | *V*-band proxy (Perelmuter & Racine 1995 / Harris M81-G01); see `chime_host_magnitudes.csv` | — |
| Inclination / *b*/*a* | This work — GALFIT `20200120E_all` | — |

### FRB 20200223B

| Role | Citation | DOI / ID |
|------|----------|----------|
| **Repeater confirmation + host (*z*, DM, *m*<sub>r</sub>, coords)** | Ibik et al. 2024, ApJ **961**, 99 | [10.3847/1538-4357/ad0893](https://doi.org/10.3847/1538-4357/ad0893) |
| Inclination / *b*/*a* | This work — GALFIT `20200223B_all` | — |

### FRB 20201124A

| Role | Citation | DOI / ID |
|------|----------|----------|
| **Repeater confirmation** | Lanman et al. 2022, ApJ **927**, 59 (CHIME repeating activity of 20201124A); Fong et al. 2021 also treat it as a repeater host | [10.3847/1538-4357/ac4bc7](https://doi.org/10.3847/1538-4357/ac4bc7) |
| **Host association, *z*, photometry (*m*<sub>r</sub> = 17.904)** | Fong et al. 2021, ApJL **919**, L23 | [10.3847/2041-8213/ac242b](https://doi.org/10.3847/2041-8213/ac242b) |
| **Mas FRB position (pair with Fong host)** | Nimmo et al. 2022, ApJL **938**, L26 | [10.3847/2041-8213/ac540f](https://doi.org/10.3847/2041-8213/ac540f) |
| **Catalog DM cross-check** | Gordon et al. 2023 Table 1 | [10.3847/1538-4357/ace5aa](https://doi.org/10.3847/1538-4357/ace5aa) |
| Inclination / *b*/*a* | This work — GALFIT `20201124A_all` | — |

### FRB 20220912A

| Role | Citation | DOI / ID |
|------|----------|----------|
| **Discovery (CHIME)** | McKinven & CHIME/FRB Collaboration, ATel **#15679** (2022 Oct 15) — first announcement | [ATel 15679](https://www.astronomerstelegram.org/?read=15679) |
| **Repeater confirmation (refereed)** | Ravi et al. 2023, ApJL **949**, L3 (treats it as "the repeating FRB source FRB 20220912A"); further refereed characterization: Hewitt et al. 2024, MNRAS **529**, 1814; Sheikh et al. 2024, MNRAS **527**, 10425 | [10.3847/2041-8213/acc4b6](https://doi.org/10.3847/2041-8213/acc4b6) · [10.1093/mnras/stae632](https://doi.org/10.1093/mnras/stae632) · [10.1093/mnras/stad3630](https://doi.org/10.1093/mnras/stad3630) |
| **Host / DSA-110 localization, *z* = 0.0771, DM, *m*<sub>r</sub> (PS1 *r* = 19.65)** | Ravi et al. 2023, ApJL **949**, L3 (host PSO J347.2702+48.7066; arXiv:2211.09049) | [10.3847/2041-8213/acc4b6](https://doi.org/10.3847/2041-8213/acc4b6) |
| Inclination / *b*/*a* | This work — GALFIT `20220912A_all` | — |

### FRB 20240114A

| Role | Citation | DOI / ID |
|------|----------|----------|
| **Repeater confirmation / CHIME discovery** | Shin & CHIME/FRB, ATel **#16420** | [ATel 16420](https://www.astronomerstelegram.org/?read=16420) |
| **Host / MeerKAT+EVN localization, *z*, DM** | Snelders et al. 2025, ApJL | [10.3847/2041-8213/ae0b68](https://doi.org/10.3847/2041-8213/ae0b68) · [arXiv:2506.11915](https://arxiv.org/abs/2506.11915) |
| **Supporting host *z* / dwarf context** | Bhardwaj et al. 2024, ApJL **971**, L51 | [10.3847/2041-8213/ad64d1](https://doi.org/10.3847/2041-8213/ad64d1) |
| **Catalog *m*<sub>r</sub> = 21.79** | This work — GALFIT (no separate literature *r* adopted) | — |
| Inclination / *b*/*a* | This work — GALFIT `20240114A_all` | — |

### FRB 20240209A

| Role | Citation | DOI / ID |
|------|----------|----------|
| **Repeater confirmation + localization + host association + *m*<sub>r</sub> = 16.79** | Shah et al. 2025, ApJL **979**, L21 | [10.3847/2041-8213/ad9ddc](https://doi.org/10.3847/2041-8213/ad9ddc) |
| **Host *z* = 0.1384 + stellar population** | Eftekhari et al. 2025, ApJL **979**, L22 | [10.3847/2041-8213/ad9de2](https://doi.org/10.3847/2041-8213/ad9de2) |
| **Catalog DM** | Shah et al. 2025 (burst properties / NE2001 DM_MW in CHIME build) | — |
| Inclination / *b*/*a* | This work — GALFIT `20240209A_all` | — |

### FRB 20251229A ⚠ preliminary (ATel-only)

| Role | Citation | DOI / ID |
|------|----------|----------|
| **Repeater confirmation + DM** | CHIME/FRB, ATel **#17574** | [ATel 17574](https://www.astronomerstelegram.org/?read=17574) |
| **Outriggers localization + host association + *m*<sub>r</sub> = 18.6** | Kahinga & CHIME/FRB, ATel **#17709** | [ATel 17709](https://www.astronomerstelegram.org/?read=17709) |
| **Host *z*, SFR, *M*★, metallicity** | Bhardwaj et al., ATel **#17856** (+ correction **#17863**) | [ATel 17856](https://www.astronomerstelegram.org/?read=17856) · [ATel 17863](https://www.astronomerstelegram.org/?read=17863) |
| Inclination / *b*/*a* | This work — GALFIT `20251229A_all` (use ATel mag, not pipeline mag) | — |
| Caveat | **No arXiv or refereed journal paper** as of 2026-07. Cite ATels only; flag as preliminary. | — |

---

## Column groups in the CSV

1. **Localization** — same schema as `CHIME/repeater_localizations.csv` / `master_frb_localization.csv`
2. **Provenance** — `discovery_facility`, `localization_facility`, `catalog_origin`
3. **Photometry** — `mag_r`, `mag_r_err`, `mag_band`, `mag_source`
4. **Morphology** — `inc_deg`, `inc_err_deg`, `b_a` (this work’s GALFIT where available)
5. **Notes** — caveats

---

## Related files (not replaced)

| File | Role |
|------|------|
| `CHIME/repeater_localizations.csv` | CHIME-only subset (16) for CHIME pipeline batch |
| `master_frb_localization.csv` | Full ASKAP/MeerKAT/DSA host sample |
| `CHIME/repeater_reported_values_sources.md` | Per-quantity provenance for the CHIME 16 |
| `CHIME/SOURCES_AUDIT.md` | Full CHIME catalog audit |
| `Archive/notes/repeater_audit_notes.md` | Master-catalog repeater audit (6 `repeater=yes` rows) |
| `CHIME/catalog/host_association_papers.md` | Paper tier list for CHIME hosts |
| `CHIME/papers/phase4_repeater_literature.md` | Phase-4 download / DOI guide |

---

## Citation hygiene

- Prefer the **primary** paper in each role above; use Gordon et al. 2023 only as a **compilation cross-check** or for *m*<sub>r</sub> when that is the adopted catalog value.
- Inclinations and *b*/*a* from this pipeline are **this work**, not literature.
- Truly telegram-only: `20251229A` (all values). For `20230814B`, the *repeater status* is telegram-only (ATel #16191), but host *z*/DM are in a refereed compilation (Li et al. 2025, ApJ **989**, 77). For `20220912A`, discovery was a telegram (ATel #15679) but repetition is confirmed in refereed papers (Ravi et al. 2023; Hewitt et al. 2024; Sheikh et al. 2024).
- **Fonseca et al. 2020** = *Nine New Repeating FRB Sources from CHIME/FRB*, **ApJL 891, L6**, DOI [10.3847/2041-8213/ab7208](https://doi.org/10.3847/2041-8213/ab7208), arXiv:2001.03595. (Do **not** confuse with `ab75f5`/vol. 908 — that was a prior typo, now fixed for 20190208A / 20190417A.)
- Author spelling in print: **Michilli** (ApJ 950, 134); repo folder name may say “Michili”.
