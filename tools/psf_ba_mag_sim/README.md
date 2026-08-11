# PSF-over-Magnitude b/a Overcorrection Simulation

Deterministic GalSim forward model + dual GALFIT fits (with / without PSF) on a rigid
parameter grid to isolate PSF-deconvolution b/a bias as a pure fitting artifact.

## Environment (WSL only)

GALFIT runs in WSL (`/usr/local/bin/galfit`). GalSim is installed via conda-forge into
the `frb_project` env (Python 3.14):

```bash
~/miniforge3/bin/conda install -n frb_project -c conda-forge galsim pyyaml
```

All scripts use `~/miniforge3/envs/frb_project/bin/python`.

**Note:** GALFIT aborts on `/mnt/c` paths. `run_galfit_grid.py` stages each fit under
`~/.psf_ba_mag_sim/<galaxy_id>/<mode>/` (WSL-native) and archives logs to `outputs/fits/`.

## Quick start

From WSL, in this directory:

```bash
PY=~/miniforge3/envs/frb_project/bin/python

# 1) Generate mocks (135 grid points at n_realizations=1)
$PY generate_mocks.py

# 2) Dual GALFIT fits (270 total: psf + nopsf per galaxy)
$PY run_galfit_grid.py

# 3) Merge, aggregate, plot
$PY analyze.py
```

Smoke test (2 mocks, 1 galaxy fit):

```bash
$PY generate_mocks.py --limit 2
$PY run_galfit_grid.py --limit 1
$PY analyze.py
```

## Configuration (`config.yaml`)

| Block | Key parameters |
|-------|----------------|
| `grid` | `intrinsic_ba`, `intrinsic_re_arcsec`, mag sweep 17–24 step 0.5 |
| `physics` | SDSS-like pixel scale 0.396″/px, 1″ Moffat PSF, sky 21.2 mag/arcsec² |
| `realizations` | `n_realizations: 1` **or** `n_realizations_by_mag` thresholds |
| `galfit` | `modes: [psf, nopsf]`, xy bounds, mag limits |

### Scaling up realizations

Single count for all magnitudes:

```yaml
realizations:
  n_realizations: 1
```

More fits at faint magnitudes (disable `n_realizations` when using this):

```yaml
realizations:
  n_realizations_by_mag:
    default: 1
    ">=21.0": 5
    ">=23.0": 20
```

## Grid size

3 b/a × 3 Re × 15 magnitudes × N realizations = **135 × N** mock galaxies.
Each galaxy gets **2** GALFIT runs (PSF / no-PSF) unless `--modes` is narrowed.

## Outputs

```
outputs/
  mocks/<galaxy_id>/mock.fits, sigma.fits
  mocks/psf.fits
  fits/<galaxy_id>/{psf,nopsf}/fit.log, out.fits, galfit.feedme
  catalogs/truth_catalog.csv
  catalogs/fit_results.csv
  catalogs/merged.csv, aggregated.csv, summary.csv
  plots/ba_vs_mag_combined.png   # 3x3 grid: true b/a vs GALFIT+PSF vs GALFIT no-PSF
```

## GALFIT feedme policy

- **Object 1 (sersic):** seeded at exact truth; n locked to 1; mag, Re, b/a, PA free;
  position bounded ±1 px only (no mag constraint — that breaks GALFIT on these stamps).
- **Object 2 (sky):** seeded at **0** on sky-subtracted stamps; **left free**.
- **C) none** — GALFIT auto-sigma (custom sigma.fits breaks magnitude).
- **Fit A:** `D) psf.fits` — **Fit B:** `D) none`.

### Known issues (v1 run)

The first full run used auto-tuned ZP (sky ~26 e-/pix), custom sigma maps, and a mag
constraint in `constraints.txt`. That produced numerically unstable fits: magnitude
snapped to ~25 while b/a looked fine, and 68/270 fits crashed. See `audit_fits.py` and
`outputs/reports/FIT_AUDIT.md`. **Regenerate mocks and re-fit after the v2 config.**

## Expected result

For Re = 0.5″ and 1.0″, recovered b/a (PSF run) should detach from truth and dive
toward ~0.2 as mag crosses ~20 (SNR@Re ~ 1). No-PSF run should show the opposite:
rounder / inflated b/a from noise, not PSF overcorrection.
