# Color cuts for null-catalog diagnostics (SDSS only)

Uses `u_r` from the augmented v1 SDSS catalog.
Density panel: **matplotlib hexbin** (`gridsize=50`, `mincnt=1`, counts ≥1 only).
See **`HEXBIN_AND_SDSS_BANDS_AUDIT.md`** for bin geometry and horizontal bands.

| u-r cut | Role | N in cut |
|---------|------|----------|
| < 3.5 | Broad blue / star-forming population | 326,568 |
| < 2.2 | Near Strateva et al. (2001) u*-r* ~ 2.22 separator | 198,447 |
| < 1.5 | Strongly blue tail | 104,921 |

Legacy g-r proxy plots: **not generated** (deferred).

## Legacy Survey (`plot_legacy_color_cuts.py`)

Tractor $g-r = \mathrm{gmag} - \mathrm{rmag}$; axis ratio = Tractor $b/a$ from $e_1$, $e_2$ (column `expAB_r`).

| g-r cut | Role | N in cut |
|---------|------|----------|
| < 0.75 | Diagnostic pool, $b/a > 0.2$, exclude REX (matches null CDF) (no mag limit) | 26,971 |

Output: `color_cuts/gr_lt_0p75/mag_vs_ba_legacy.png`
