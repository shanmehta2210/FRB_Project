# FRB repeater audit — working notes

Companion to `master_frb_localization.csv` columns `repeater` and `repeater_source`.

Scope: fill **host-localized** rows only (`coord_semantics=host`); leave signal-localized rows blank unless noted.

**Survey scope:** ASKAP / MeerKAT / DSA / Legacy-oriented hosts only — **not CHIME**. FRB 20180916B (CHIME repeater) is **excluded** and will not be added to this sample.

**Primary citation for reference-set matches:** Heintz et al. 2020, ApJ 903, 152 (doi:10.3847/1538-4357/abb6fb)

---

## User reference set — matches in `master_frb_localization.csv`

All use `repeater_source` = Heintz et al. 2020 (doi:10.3847/1538-4357/abb6fb).

| FRB | `coord_semantics` | repeater | Notes |
|-----|-------------------|----------|-------|
| 20180924B | host | no | |
| 20181112A | signal | no | signal-localized; filled for reference only |
| 20190102C | host | no | |
| 20190523A | host | no | |
| 20190608B | host | no | |
| 20190611B | host | no | |
| 20190614D | host | no | |
| **20190711A** | host | **yes** | only repeater in reference set |
| 20190714A | host | no | |
| 20191001A | host | no | |
| 20200430A | host | no | |

**Not in sample (by design):** FRB **20180916B** — CHIME repeater; outside survey scope for this project.

---

## Progress log

| Batch | FRBs | Status |
|-------|------|--------|
| 1 | 20171020A–20190714A (first 5 hosts in CSV order) | filled |
| 2 | 20191001A, 20200906A, 20210320C, 20210410D, 20210807D | filled, pending user review |

Reference-set rows above supersede earlier per-FRB citations where they overlap (20180924B, 20190102C, 20190608B, 20190714A, 20191001A now Heintz et al. 2020).
