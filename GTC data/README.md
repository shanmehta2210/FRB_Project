# GTC data

Observation planning and GTC-specific inputs for this project. **Not** part of the
SExtractor → AstroPath → GALFIT imaging pipeline in `pipeline_scripts/`.

## Layout

| Path | Purpose |
|------|---------|
| `environment.yml` | Conda env `gtc_visibility` (astropy, astroplan, pandas) |
| `gtc_portal_coordinates.txt` | Legacy portal copy-paste list (99 FRBs); superseded by visibility batch |
| `visibility/` | Rigorous visibility filtering scripts and reports |
| `production_62/` | Tiered GTC science review for the 62 production fitted hosts |
| `pipeline_trial/` | 13-FRB archival pipeline trial cohort (manifest + bad-fit tracking) |

## Environment (WSL)

```bash
conda env create -f "GTC data/environment.yml"   # first time
conda activate gtc_visibility
conda env update -f "GTC data/environment.yml" --prune
```

The imaging pipeline uses `frb_project` (`pipeline_scripts/environment_frb_project.yml`).
`astroplan` lives only in `gtc_visibility`.
