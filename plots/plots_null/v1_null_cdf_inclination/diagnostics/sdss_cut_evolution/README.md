# SDSS null cos(i) cut evolution

## Directory layout

| Path | Content |
|------|---------|
| `cdf_00_…` – `cdf_03_…` | Cut evolution @ `modelMag_r <= 21` (strict b/a first) |
| `final_by_mag/` | **Full-pool empirical CDF** @ mag 20, 21, 22 |
| `mc_sensitivity/` | Same pools: ECDF vs MC null at several subsample sizes |
| `mc_sensitivity/mc_audit.csv` | Numeric audit (see below) |
| `MC_AUDIT.txt` | How production `cdf_envelope` works |

## Production overlay (`mag_cuts/…/sdss_strict/`)

**Pool:** `catalog/SDSS_catalog_v1_allsky_modelmr.csv` → `u-r < 2.3`, lnL exp-wins, `expAB_r > 0.2`, `modelMag_r <= limit`.

**Green curve:** `cdf_envelope(pool, n_sample=N_FRB, n_draws=10000)` in `pipeline_null_plot_utils.py`.

- **10,000** independent draws.
- Each draw: pick **N_FRB** galaxies (41 @ mag21) **without replacement**, sort cos(i), build step CDF, evaluate on 100-point grid.
- Plotted curve = **mean** over 10,000 draws (68% band = 16–84% of draws).

This is **not** the ECDF of all ~27k galaxies.

## Two readings (easy to confuse)

| Quantity | mag21 typical | Meaning |
|----------|---------------|---------|
| **cdf_at_cosi_0.5** (y at x=0.5) | ~0.38 | Fraction of sample with cos i **<** 0.5 |
| **x_at_cdf_0.5** | ~0.58 | **Median** cos i (where CDF crosses 0.5) |

Evolution `final_by_mag` plots use **full-pool ECDF** (all N in pool).  
Overlay uses **MC with n=41**. Those differ in **smoothing**, not in median (~0.58): subsampling 41 from 27k does **not** push y@0.5 to ~0.65.

See `mc_sensitivity/mc_sensitivity_mag21.png`: curves converge as `n_sample` increases.

## Regenerate

```bash
python scripts/plot_sdss_null_cut_evolution_cdfs.py
```
