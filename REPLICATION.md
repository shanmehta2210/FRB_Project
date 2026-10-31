# Replicating the FRB host-inclination analysis

Everything needed to regenerate the paper-level results (N = 56 confirmed, non-repeating hosts with published DM) from a fresh clone. Large imaging (FITS) is not in git; it is either re-downloadable or only needed to re-run the GALFIT stage, whose outputs (fit logs, tables) are tracked.

## 1. Clone

```bash
git clone --recurse-submodules https://github.com/shanmehta2210/FRB_Project.git
cd FRB_Project
```

## 2. Python environments

| File | Use |
|---|---|
| `env/requirements_venv.txt` | main analysis venv (Python 3.13.12) |
| `env/requirements_venv_astrophot.txt` | AstroPhot verification venv (Python 3.10.11) |
| `env/requirements_system_python.txt` | system Python snapshot (reference only) |
| `pipeline_scripts/environment_frb_project.yml` | conda alternative |

```bash
python -m venv venv && venv/Scripts/activate      # Windows; source venv/bin/activate on Linux/macOS
pip install -r env/requirements_venv.txt
```

Key pins: numpy 2.4.2, scipy 1.15.3, pandas 2.2.3, astropy 7.2.0, matplotlib 3.10.3, mwprop 2.0.5 (NE2001), pyvo 1.8, astrophot 0.16.13.

External binaries (only to re-run image fitting): GALFIT, SExtractor, PSFEx.

## 3. Restore the HSC control catalogue

The HSC-PDR2 (Kawinwanichakij et al. 2021) control is stored compressed in `data/catalog_snapshots/`:

```bash
python scripts/restore_catalogs.py     # writes catalog/*.csv and checks SHA256
```

## 4. Local AstroPath modifications

The `FRBs/astropath` submodule is pinned at `4ab7d74b`; our changes are kept as a patch:

```bash
git submodule update --init                     # if not cloned with --recurse-submodules
git -C tools/AstroPath/astropath_pkg apply ../local_changes/astropath_pkg.patch
mkdir -p tools/AstroPath/astropath_pkg/scripts
cp tools/AstroPath/local_changes/scripts/*.py tools/AstroPath/astropath_pkg/scripts/
pip install -e tools/AstroPath/astropath_pkg
```

## 5. Re-run the statistics (no imaging needed)

Inputs are tracked: `master_frb_localization.csv`, `pipeline_scripts/verification/host_confirmation.csv`, GALFIT results under `pipeline_scripts/Output/`, literature host properties under `Publication/materials/host_properties/`.

```bash
python scripts/plot_frb_vs_hsc_confirmed50.py   # builds confirmed50_q.csv (N = 56) + CDF figures
python scripts/run_frb_hsc_confirmed50_tests.py # KS, AD, Kuiper, CvM, MWU, binomial edge-on
python scripts/run_frb_hsc_error_mc.py          # Protocol B error Monte Carlo
python scripts/run_frb_hsc_nature23_tests.py    # 4/23 Bhardwaj et al. comparison
python scripts/run_host_property_correlations.py# DM_host, RM, tau, M*, SFR, ... vs inclination
python scripts/recheck_marnoch2026.py           # Marnoch et al. 2026 reproduction
python scripts/appendix_controls_tables.py
```

Outputs land in `plots/plots_null/v2/` and `Publication/materials/`. The summary of all results is `Publication/materials/PROGRESS_REPORT_2026-10.md`.

`recheck_marnoch2026.py` parses `papers/marnoch2026.txt` (text of arXiv:2610.13062) if present, otherwise uses the tracked `Publication/materials/marnoch_recheck/marnoch2026_table1.csv`.

## 6. Re-running the imaging stage (optional)

1. Cutouts: `large_cutouts/cutout_registry.csv` lists every host's position and source survey. Re-fetch with `python scripts/cutout_download.py <FRB>` (Legacy, then PS1, then DES).
2. PSFs and fits: `pipeline_scripts/master_run.py` / `pipeline_scripts/run_all_frbs.py` (see `pipeline_scripts/README.md`).

## 7. Literature PDFs

Not tracked (copyright). `Publication/materials/host_properties/PAPERS_MANIFEST.csv` lists every file used, with arXiv IDs where known, so they can be re-downloaded from arXiv.

## Not in git

| Item | Size | Recovery |
|---|---|---|
| `large_cutouts/**/*.fits` | 8.2 GB | re-download via registry |
| `pipeline_scripts/Output/**/*.fits, *.cat, *.png` | ~4.3 GB | re-run pipeline; text outputs are tracked |
| `catalog/` (DES, SDSS, LS audit catalogues) | 2.6 GB | `scripts/build_*`/`fetch_*` scripts; HSC is restored from snapshot |
| `catalog_downloads/` | 1.5 GB | re-query |
| paper PDFs | 0.5 GB | `PAPERS_MANIFEST.csv` |
| `CHIME/`, `PATH/` | 4.4 GB | separate AstroPATH-additions repo / external backup |
