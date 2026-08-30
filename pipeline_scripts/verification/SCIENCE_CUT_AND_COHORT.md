# Science cut and cohorts

## Definitions

| name | definition | where |
|---|---|---|
| **all64** | Every FRB with a production GALFIT workdir | `vercommon.cohort("all64")` |
| **science cut / in_53** | `mag ≤ 22` and `b/a > 0.2` | `vercommon.MAG_CUT`, `BA_CUT`; column `in_53` |
| **confirmed** | Human triage pass in `host_confirmation.csv` | not the same as `in_53` |

Cuts use production catalog `mag` and `b_a` from `pipeline_galfit_results.csv`
(merged in `vercommon.cohort`).

## When the cut applies

| Stage | Cohort |
|---|---|
| Run verification suite | **all64** by default (`--cohort all64`) |
| Aggregate plots / population tests | Filter to **in_53** for science-facing figures |
| Visual triage walk | **all64** — do not skip faint hosts during first pass ([`HOST_TRIAGE_CASES.md`](HOST_TRIAGE_CASES.md)) |
| Paper inclination sample | **in_53 ∩ confirmed** |
| Reject grid (`run_reject_grid.py`) | Historically **in_53 ∩ confirmed=False**; list can lag CSV as hosts flip |

## Interaction with confirmation

- A host can be `confirmed=True` **outside** the cut (recorded, not in the 53).
- A host can be `in_53` and `confirmed=False` (rejected for morphology / sky / etc.).
- Trust tiers A/B/C from `aggregate.py` are automatic; they do **not** override
  the CSV.

## Current counting convention

After triage:

- Science cut size: **53**
- Confirmed in cut (paper sample): **50**
- Rejected in cut (review pile): **3** — `20190711A` (wrong object / star), `20221101B` (morphology uncertain; MW star contaminates host), `20230930A` (structured diffuse light / bad sky)

Total CSV rows: **64**. Recompute with:

```python
import pandas as pd
import vercommon as vc
hc = pd.read_csv("host_confirmation.csv")
hc["confirmed"] = hc["confirmed"].astype(str).str.lower().eq("true")
m = hc.merge(vc.cohort("all64")[["frb", "mag", "b_a", "in_53"]], on="frb")
print((m.confirmed & m.in_53).sum(), "/", m.in_53.sum())
```

## Selection-bias note

Cutting on measured `b/a` and mag can itself bias inclination distributions.
See [`FIT_VERIFICATION_CHECKS.md` §1 / §9](FIT_VERIFICATION_CHECKS.md). Triage
documents *fit quality*, not a re-derivation of the survey selection function.
