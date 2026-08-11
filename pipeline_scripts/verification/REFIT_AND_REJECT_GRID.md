# Re-fits and reject grid

Staging lives under `Re-fits/<FRB>/`. Production `Output/<FRB>_all/` is
**never** overwritten. Short layout note: [`Re-fits/README.md`](Re-fits/README.md).

## Design rules

- **`q` is never fixed** in standard legs — only host \(n\) and/or sky.
- A re-fit is not “new `out.fits` only”: the **full verification suite** is
  re-run into the leg directory (RFF, isophotes, Fourier, sky±, AstroPhot,
  visual, …).
- Sky± legs inside a parent that freezes \(n\) also freeze \(n\)
  (`host_n_held_fixed` / `refit_meta.json`).
- AstroPhot locks \(n\) when GALFIT does (avoids runaway \(n\sim 8\)).

## `run_refit.py`

Stages sidecars from production, edits `galfit.feedme`, runs GALFIT, runs
`run_verification.run_host_dir` into the same workdir.

```bash
python run_refit.py 20181112A --fix-n 1 --label n1
python run_refit.py 20181112A --sky-from-protocol --label sky
python run_refit.py 20181112A --sky-from-protocol --fix-n 1 --label n1_sky
python run_refit.py 20181112A --sky-adu 6.3e-5 --fix-n 1 --label n1_sky
```

| flag | effect |
|---|---|
| `--fix-n N` | host Sérsic index fixed at `N` (first Sérsic only) |
| `--sky-adu X` | sky component fixed at `X` ADU |
| `--sky-from-protocol` | sky = `Re-fits/<FRB>/sky_protocol.json` consensus |
| `--label` | subdir name under `Re-fits/<FRB>/` |
| `--checks` | subset of suite (default all) |

### Staging / feedme edits

Copied sidecars (non-exhaustive): `host_cutout.fits`, `host_sigma.fits`,
`host_mask.fits`, `proto_image.fits`, `constraints.txt`, `image.cat`,
`psfex.xml`, `cutout_meta.json`, photometry/ZP JSONs, …

Feedme edits (`_edit_feedme_fixed`):

- Reseed free Sérsic params from production best-fit (unless disabled).
- Sky line → fixed (`flag 0`) or free reseeds.
- Host `n` → fixed when requested; strip `n` rows from constraints.
- When sky is fixed, drop sky constraint rows for that component.

Writes `refit_meta.json` (what was held) and `refit_summary.json` (fit +
verification status).

### Panel titles

Refit panels can show `[n1]`, `[sky]`, or `[n1+sky]` via `refit_meta.json`
(see `checks/visual.py`).

## Standard legs

| label | \(n\) | sky | typical use |
|---|---|---|---|
| `n1` | fixed 1 | free | tame bullshit \(n\) / Re–n degeneracy |
| `sky` | free | protocol fixed | stop free-sky runaway |
| `n1_sky` | fixed 1 | protocol fixed | both pathologies |

Directory:

```
Re-fits/<FRB>/
  sky_protocol.json
  n1/   sky/   n1_sky/     # full suite each
  panel_production.png     # byte copy of outputs/panels/<FRB>.png
  panel_n1.png
  panel_sky.png
  panel_n1_sky.png
  reject_grid_summary.json
  sandbox/                 # optional; see SANDBOX.md
```

## `run_reject_grid.py`

For science-cut rejects: builds the four panels above.

1. **Byte-copy** production panel → `panel_production.png` (never regenerate).
2. Run `sky_protocol` if needed.
3. Run legs `n1`, `sky`, `n1_sky` via `run_refit`.
4. Copy each leg’s `panel.png` up to `panel_<leg>.png`.

Cohort selection: prefers `host_confirmation.csv` ∩ `in_53` ∩ `confirmed=False`;
falls back to a hardcoded list that can go **stale** as hosts are confirmed —
recompute from CSV when in doubt.

## Promoting a winning leg

1. Visually pick the leg (and check suite JSONs under that leg).
2. Set `host_confirmation.csv` → `confirmed=True` with notes citing the leg.
3. Copy `Re-fits/<FRB>/<leg>/panel.png` → `outputs/panels/<FRB>_<leg>.png`.
4. Update [`HOST_TRIAGE_CASES.md`](HOST_TRIAGE_CASES.md).

## What re-fits do *not* do

- Do not change the science cut definition.
- Do not auto-update `IN53_REJECTS` hardcodes after confirmation flips.
- Do not replace sandbox hand-tuning ([`SANDBOX.md`](SANDBOX.md)) for
  multi-component weirdos.
