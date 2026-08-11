# Legacy Survey v2 EXP — scaled b/a audit

Catalog: `catalog/LS_catalog_v2_fullsky_exp.csv` (Tractor `type=EXP`).

Face-on cap matches CDF `strict_scaled`: keep **b/a ≤ 0.8**
(including b/a < q0=0.2; strict q0 cut is CDF-only). Plot the stretched
axis ratio

$$(b/a)' = (b/a) / 0.8$$

so the former face-on edge at b/a=0.8 maps to 1.

| Metric | Value |
|--------|-------|
| Rows (finite mag + b/a ≤ 0.8) | 1,997,809 |
| RA span | 0.000 … 360.000 |
| Dec span | -89.557 … 84.767 |
| median modelMag_r | 22.643 |
| median (b/a)' | 0.5557 |
| mean (b/a)' | 0.5553 |

## Plots

| File | Content |
|------|---------|
| `mag_histogram_ba.png` | Counts per 0.5 mag bin |
| `median_ba_vs_mag.png` | Median (b/a)' per mag bin |
| `ba_vs_mag_scatter.png` | (b/a)' vs mag hexbin (y→1 at former ba=0.8) |

## Regenerate

```bash
python scripts/audit_ls_v2_mag_ba.py --scaled
```
