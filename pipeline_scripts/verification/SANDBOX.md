# Manual GALFIT sandbox

Script: `run_sandbox.py`.  
Workdirs: `Re-fits/<FRB>/sandbox/`.

Use this when reject-grid legs are not enough — hand-edited feedme, extra PSF
components, frozen centroids, iterative `galfit.NN` archives. **Does not** run
the verification suite by itself (you can run checks into the sandbox later).

Default weird hosts seeded by `--init`: `20220509G`, `20230930A`.

## When to use sandbox vs `run_refit`

| | `run_refit` | `run_sandbox` |
|---|---|---|
| Feedme | auto-edited from production | you edit by hand |
| Suite | full verification every time | GALFIT + quick 3-panel only |
| Archives | overwrites leg products | keeps `galfit.NN` / named files |
| Best for | controlled n/sky legs | multi-component rescue |

## Workflow

```bash
# Once: copy cutout/sigma/mask/PSF/constraints + seed feedme from production
python run_sandbox.py --init
python run_sandbox.py 20220509G --init --reseed-feedme   # overwrite feedme

# Edit Re-fits/20220509G/sandbox/galfit.feedme

# Run GALFIT → out.fits, fit.log, panel.png
python run_sandbox.py 20220509G

# Rebuild panel only
python run_sandbox.py 20220509G --panel-only
```

### Archive hygiene

- Runner refreshes only `fit.log`, `out.fits`, `panel.png`.
- **Does not delete** `galfit.NN` or renamed archives (`cent_fix`, `psf`, …).
- Rename good restart files yourself, e.g. `galfit.01` → `cent_fix` or `psf`.
- Optional: copy `out.fits` → `out_psf.fits` when freezing a milestone.

### Full visual suite on a sandbox fit

After a good `out.fits` / `fit.log`:

1. Ensure sidecars exist (`image.cat`, `psfex.xml`, …) — copy from production.
2. `python -c` / small script calling `run_verification.run_host_dir(frb, sandbox, sandbox, CHECK_ORDER, force=True)`.
3. Publish `panel.png` → `panel_psf.png` and `outputs/panels/<FRB>_psf.png` if confirmed.

`20220509G` was confirmed this way (host Sérsic + PSF for NE star).

## Quick panel stretch (sandbox `make_panel`)

| panel | stretch |
|---|---|
| data / model | grayscale **asinh** over data **1–99%** percentiles (soft=10) |
| residual | linear `RdBu_r`, ±98th percentile of \|resid\| (**not** /σ) |

The full verification visual uses /σ for the residual
([`VISUAL_PANELS.md`](VISUAL_PANELS.md)).

## Lessons from the two weird hosts

### `20220509G` — confirmed (sandbox `psf`)

- Dual nucleus / star NE of host; single Sérsic left a catastrophe residual.
- Fix: free-centroid Sérsic + free-centroid **PSF** component at the bright peak.
- χ²/ν dropped ~160 → ~3; published as `outputs/panels/20220509G_psf.png`.
- Caveat: merger-like / unresolved structure still in notes.

### `20230930A` — rejected

- Huge diffuse / structured background + spiderweb pattern noise; ZP failed.
- Free sky runs away (physically expected); \(q\) not inclination-usable.
- Protocol / n=1 experiments do not fix structured light under the host.
- CSV: `REJECTED - whole lotta diffuse light / structured background; bad sky`.

Literature pointers for structured sky vs GALFIT’s planar sky: see triage
discussion and Peng TOP10 / TFAQ (sky #1 systematic); destriping helps lines,
not nebulosity-on-host.
