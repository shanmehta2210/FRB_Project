# HSC Kawinwanichakij+2021 — mag / b/a audit

Catalog: `catalog/HSC_kawinwanichakij_sample_500k.csv`

| Metric | Value |
|--------|-------|
| Rows (finite mag + b/a) | 499,956 |
| RA span | 33.652 … 353.905 |
| Dec span | -5.990 … 56.812 |
| median mag | 23.201 |
| median b/a | 0.6170 |
| mean b/a | 0.6143 |

Selection: full sample  |  no color / q / n cuts

## Axis ratio

$b/a=$ fitted_q

## Plots

| File | Content |
|------|---------|
| `mag_histogram_ba.png` | Galaxy counts per 0.5 mag bin (bright open bin `<15`; display through mag ≤ 26) |
| `median_ba_vs_mag.png` | Median b/a per mag bin (same display cut) |
| `ba_vs_mag_scatter.png` | b/a vs mag (log hexbin + thin scatter overlay; mag ≤ 26) |

CSV companions retain the full magnitude range: `mag_histogram_ba.csv`, `median_ba_vs_mag.csv`.

## Regenerate

```bash
python scripts/audit_des_hsc_mag_ba.py --survey hsc
```
