# HSC Kawinwanichakij+2021 — mag / b/a audit

Catalog: `catalog/HSC_kawinwanichakij_exp_analogue.csv`

| Metric | Value |
|--------|-------|
| Rows (finite mag + b/a) | 140,020 |
| RA span | 33.654 … 353.902 |
| Dec span | -5.989 … 56.812 |
| median mag | 22.755 |
| median b/a | 0.5790 |
| mean b/a | 0.5874 |

Selection: EXP analogue  0.4<n<1.5  |  no color / q cuts

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
python scripts/audit_des_hsc_mag_ba.py --survey hsc --exp-analogue
```
