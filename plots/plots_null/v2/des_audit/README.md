# DES Y1 morph (Tarsitano+2018) — mag / b/a audit

Catalog: `catalog/DES_y1_morph_sample_500k.csv`

| Metric | Value |
|--------|-------|
| Rows (finite mag + b/a) | 500,000 |
| RA span | 0.000 … 360.000 |
| Dec span | -66.731 … 2.510 |
| median mag | 22.577 |
| median b/a | 0.6853 |
| mean b/a | 0.6337 |

Selection: full sample  |  no color / q / n cuts

## Axis ratio

$b/a=1-\varepsilon$ (already calibrated; App. B)

## Plots

| File | Content |
|------|---------|
| `mag_histogram_ba.png` | Galaxy counts per 0.5 mag bin (bright open bin `<15`; display through mag ≤ 26) |
| `median_ba_vs_mag.png` | Median b/a per mag bin (same display cut) |
| `ba_vs_mag_scatter.png` | b/a vs mag (log hexbin + thin scatter overlay; mag ≤ 26) |

CSV companions retain the full magnitude range: `mag_histogram_ba.csv`, `median_ba_vs_mag.csv`.

## Regenerate

```bash
python scripts/audit_des_hsc_mag_ba.py --survey des
```
