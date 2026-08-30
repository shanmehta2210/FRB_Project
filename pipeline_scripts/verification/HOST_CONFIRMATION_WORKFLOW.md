# Host confirmation workflow

## What `confirmed` means

From [`HOST_TRIAGE_CASES.md`](HOST_TRIAGE_CASES.md):

> `confirmed = True` means residuals / discrepancies are judged **physically
> reasonable** and the host’s \(q\) (etc.) may go into the paper directly.

It is a **human** decision. It is not:

- passing every automated flag,
- trust tier A,
- or “χ²/ν ≈ 1”.

## CSV schema — `host_confirmation.csv`

| column | type | meaning |
|---|---|---|
| `frb` | string | FRB name (e.g. `20171020A`) |
| `confirmed` | bool-like | see tokens below |
| `notes` | string | short reason + panel/leg pointers |

Accepted truthy/falsey tokens (`aggregate._load_confirmation`):  
`True`/`False`, `true`/`false`, `yes`/`no`, `1`/`0`. Blank = unreviewed.

One row per fitted host (64).

## Recommended notes convention

Keep notes short; put the long argument in `HOST_TRIAGE_CASES.md`.

When confirmation uses a **re-fit leg** (not production), cite leg + panel:

```text
confirmed on n=1 fixed re-fit; panel outputs/panels/20181112A_n1.png
confirmed on n1+sky; unresolved; panel outputs/panels/20230526A_n1_sky.png
confirmed on sandbox host+PSF; panel outputs/panels/20220509G_psf.png
REJECTED - whole lotta diffuse light / structured background; bad sky
```

Patterns the tooling understands for alternate panels:

- `outputs/panels/<FRB>_n1.png`
- `outputs/panels/<FRB>_sky.png`
- `outputs/panels/<FRB>_n1_sky.png`
- `outputs/panels/<FRB>_psf.png` (sandbox host+PSF archive)

## End-to-end loop

1. **Measure** — `python run_verification.py --checks all --jobs 4`
2. **Look** — `outputs/panels/<FRB>.png` (synced from `per_host/<FRB>/panel.png`)
3. **Write case** — add/update section in `HOST_TRIAGE_CASES.md`
4. **Gate** — set `confirmed` + `notes` in `host_confirmation.csv`
5. **If reject needs rescue** — [`REFIT_AND_REJECT_GRID.md`](REFIT_AND_REJECT_GRID.md) / [`SANDBOX.md`](SANDBOX.md)
6. **Publish winning panel** — copy leg panel to `outputs/panels/<FRB>_<leg>.png` and point `notes` at it
7. **Aggregate** — `python run_verification.py --aggregate-only` (merges confirmation into metrics)

Production `Output/<FRB>_all/` is never overwritten by confirmation.

## Relationship of artifacts

| artifact | role |
|---|---|
| `host_confirmation.csv` | machine-readable gate |
| `HOST_TRIAGE_CASES.md` | narrative evidence |
| `outputs/panels/<FRB>.png` | production (or last synced) visual |
| `outputs/panels/<FRB>_<leg>.png` | paper-truth panel when production was superseded |
| `outputs/tables/fit_verification_metrics.csv` | includes `confirmed` after aggregate |
| `outputs/confirmed_fit_panels.pptx` | deck of confirmed panels ([`VISUAL_PANELS.md`](VISUAL_PANELS.md)) |

## Confirmed-on-refit hosts (convention)

If production \(n\) or sky is pathological but a constrained leg is physical:

- leave production suite products as historical,
- confirm the **leg**,
- publish `outputs/panels/<FRB>_<leg>.png`,
- note unresolved / PSF-trap caveats when \(R_e/\mathrm{FWHM}\) is tiny.

`q` is **never** held fixed in standard legs — only \(n\) and/or sky
([`REFIT_AND_REJECT_GRID.md`](REFIT_AND_REJECT_GRID.md)).
