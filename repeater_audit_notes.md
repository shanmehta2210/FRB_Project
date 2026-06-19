# FRB repeater audit — working notes

Companion to `master_frb_localization.csv` columns `repeater` and `repeater_source`.

**Scope:** `coord_semantics=host` only (84 rows). Signal-localized rows stay blank unless noted.

**Survey scope:** ASKAP / MeerKAT / DSA / Legacy-oriented hosts — **not CHIME-primary sample**. FRB **20180916B** (CHIME repeater) is excluded and not in this table.

**Status (2026-06-11):** 81/84 host rows have verified `repeater` + `repeater_source`. Three DSA hosts lack a published repeater classification (see §Pending).

Apply updates: `python scripts/apply_repeater_audit.py`

---

## CHIME internal check

| Check | Result |
|-------|--------|
| `survey` column contains `CHIME` | **None** (0/84 hosts) |
| FRB 20180916B in table | **Absent** (excluded by design) |
| CHIME-discovered sources present via other localizers | **3** (see below) |

**CHIME-associated hosts in this sample** (localized by ASKAP/DSA/MeerKAT — not CHIME survey column):

| FRB | repeater | Notes |
|-----|----------|-------|
| 20201124A | yes | CHIME discovery; ASKAP host — Fong et al. 2021 |
| 20220912A | yes | CHIME discovery; DSA host — Ravi et al. 2023 |
| 20230814B | yes | DSA repeater ATel #16191 (source name FRB 20230814A); CHIME VOEvent coincidence |

**Conclusion:** No host-semantics row is a CHIME-survey localization. CHIME-only repeater 20180916B is not in the sample.

---

## Reference set (user calibration image)

Heintz et al. 2020, ApJ 903, 152 (doi:10.3847/1538-4357/abb6fb) applies to several reference matches; primary discovery papers are preferred where they state repeater status explicitly.

| FRB | repeater | `repeater_source` in CSV |
|-----|----------|--------------------------|
| 20180924B | no | Bannister et al. 2019 (CRAFT) |
| 20190102C | no | Heintz et al. 2020 |
| 20190523A | no | Ravi et al. 2019 (DSA-10) |
| 20190608B | no | Chittidi et al. 2021 |
| 20190611B | no | Gordon et al. 2023 Table 1 |
| 20190614D | no | Law et al. 2020 (RealFAST) |
| **20190711A** | **yes** | Kumar et al. 2021 |
| 20190714A | no | Heintz et al. 2020 |
| 20191001A | no | Bhandari et al. 2020 |
| 20200430A | no | Gordon et al. 2023 Table 1 |

---

## Pending (no published repeater citation)

Host association from internal pipeline (`new_confident_hosts.txt` cites **Verdi+2025, in prep.**); no peer-reviewed repeater statement found as of 2026-06-11:

| FRB | survey | z |
|-----|--------|---|
| 20230913 | DSA-110 | 0.3024 |
| 20240104A | DSA-110 | 1.33 |
| 20250518 | DSA-110 | 0.6392 |

---

## Literature registry (host rows)

Retrieval date: **2026-06-11**. Repeater column: **Y** = observed repeater; **N** = apparent non-repeater / not observed to repeat in cited work.

### Repeaters (6)

| FRB | rep | Source | DOI / ID |
|-----|-----|--------|----------|
| 20121102A | Y | Chatterjee et al. 2017, Nature 541, 58 | 10.1038/nature20797 |
| 20180301A | Y | Bhandari et al. 2022, AJ 163, 69 | 10.3847/1538-3881/ac3aec |
| 20190711A | Y | Kumar et al. 2021, MNRAS 500, 2525 | 10.1093/mnras/staa3436 |
| 20201124A | Y | Fong et al. 2021, ApJL 919, L23 | 10.3847/2041-8213/ac242b |
| 20220912A | Y | Ravi et al. 2023, ApJL 949, L3 | 10.3847/2041-8213/acc4b6 |
| 20230814B | Y | Ravi, ATel #16191 (2023 Aug 16); DSA name FRB 20230814A | https://www.astronomerstelegram.org/?read=16191 |

### Gordon et al. 2023 — Table 1 repeater column (CRAFT / MeerTRAP hosts)

| FRB | rep | DOI |
|-----|-----|-----|
| 20180924B | N | 10.3847/1538-4357/ace5aa |
| 20190102C | N | (also Heintz 2020) |
| 20190608B | N | (also Chittidi 2021) |
| 20190611B | N | |
| 20190714A | N | |
| 20191001A | N | |
| 20200430A | N | |
| 20200906A | N | |
| 20210320C | N | |
| 20210410D | N | (also Caleb 2023) |
| 20210807D | N | |
| 20211127I | N | |
| 20211203C | N | |
| 20211212A | N | |
| 20220105A | N | |

### Law et al. 2024 — DSA first catalog (11 FRBs; not observed to repeat)

ApJ 967, 29 — doi:10.3847/1538-4357/ad3736

20220207C, 20220307B, 20220310F, 20220319D, 20220418A, 20220506D, 20220509G, 20220825A, 20220914A, 20220920A, 20221012A — all **N**.

### Sharma et al. 2024 — DSA extended localized sample

Nature 635, 61 — doi:10.1038/s41586-024-08074-9

Used for additional DSA hosts in Nature/Connor tables — all **N** (single-burst localized population study).

### Connor et al. 2024 — DSA baryon / extended FRB sample

arXiv:2409.16952 — doi:10.48550/arXiv.2409.16952

20230521B — **N** (single-burst localization; highest-z entry in that work).

### Hussaini et al. 2025 — localized FRB foreground-LSS catalog

ApJL 993, L27 — doi:10.3847/2041-8213/ae0a49

20231220A, 20240119A, 20240123A, 20240203, 20240213A, 20240215A, 20240229A — **N** (tabulated localized bursts; no repeat reported in that paper).

### Gordon et al. 2025 — CRAFT spatial mapping (33 non-repeaters in sample)

ApJ 993, 119 — doi:10.3847/1538-4357/ae0298

20220725A, 20221106A, 20230526A, 20230902A, 20231226A, 20240201A, 20240208A, 20240210A, 20240304A, 20240310A, 20240318A — **N**.

### Other primary sources

| FRB | rep | Source | DOI |
|-----|-----|--------|-----|
| 20171020A | N | Li et al. 2023, PASA 40, 29 | (PASA; URL in CSV) |
| 20180924B | N | Bannister et al. 2019, Science 365, 565 | 10.1126/science.aaw5903 |
| 20190523A | N | Ravi et al. 2019, Nature 572, 352 | 10.1038/s41586-019-1389-7 |
| 20190614D | N | Law et al. 2020, ApJ 899, 161 | 10.3847/1538-4357/aba4ac |
| 20201123A | N | Rajwade et al. 2022, MNRAS 514, 1961 | 10.1093/mnras/stac1450 |
| 20220610A | N | Ryder et al. 2023, Science 382, 294 | 10.1126/science.adf2678 |
| 20220717A | N | Rajwade et al. 2024, MNRAS 532, 3881 | 10.1093/mnras/stae1652 |
| 20230708A | N | Muller et al. 2026, ApJ 1001, 118 | 10.3847/1538-4357/ae5060 |
| 20230930A | N | Anna-Thomas et al. 2025, ApJ 993, 221 | 10.3847/1538-4357/ae1014 |
| 20240304B | N | Caleb et al. 2025, arXiv:2508.01648 | 10.48550/arXiv.2508.01648 |
| 20220222C–20231020B | N | Pastor-Marazuela et al. 2026, MNRAS 545, f2144 | 10.1093/mnras/staf2144 |

Pastor-Marazuela et al. abstract: *15 apparently non-repeating FRBs* (arXiv:2507.05982).

### Sherman et al. 2023 / 2024 (DSA polarimetry; non-repeating sample)

| Paper | DOI | Note |
|-------|-----|------|
| Sherman et al. 2023, ApJL 957, L8 | 10.3847/2041-8213/ad0380 | 10 as yet nonrepeating DSA FRBs |
| Sherman et al. 2024, ApJ 964, 131 | 10.3847/1538-4357/ad275e | 25 as yet nonrepeating DSA FRBs |

---

## Progress log

| Date | Action |
|------|--------|
| 2026-06-11 | Filled 81/84 host rows; CHIME check passed; 3 DSA hosts pending published repeater lit |
| Earlier | Batches 1–2; reference-set Heintz 2020 notes |
