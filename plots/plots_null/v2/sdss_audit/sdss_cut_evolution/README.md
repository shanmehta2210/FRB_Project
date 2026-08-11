# SDSS v2 null cos(i) cut evolution

Catalog: `catalog/SDSS_catalog_v2_fullsky_modelmr.csv` (HTM random, full footprint).

## Directory layout

| Path | Content |
|------|---------|
| `cdf_00_…` – `cdf_03_…` | Cut evolution @ `modelMag_r <= 21` |
| `final_by_mag/` | Full-pool empirical CDF @ mag 20, 21, 22 |
| `final_by_mag/cdf_final_empirical_panel_mag20_21_22.png` | **3-panel** ECDF + dashed uniform |
| `mc_sensitivity/` | ECDF vs MC null at several subsample sizes |
| `final_by_mag/final_by_mag_summary.csv` | Median cos(i) audit |

## Production pools (strict + u-r + lnL exp-wins)

| mag limit | N pool | median cos(i) | x @ CDF=0.5 |
|-----------|--------|---------------|-------------|
| 20 | 55,679 | 0.613 | 0.615 |
| 21 | 159,144 | 0.584 | 0.585 |
| 22 | 342,534 | 0.523 | 0.525 |

## Regenerate

```bash
python scripts/plot_sdss_null_cut_evolution_cdfs.py \
  --sdss-csv catalog/SDSS_catalog_v2_fullsky_modelmr.csv \
  --out-dir plots/plots_null/v2/sdss_audit/sdss_cut_evolution \
  --final-mag-limits 20 21 22
```
