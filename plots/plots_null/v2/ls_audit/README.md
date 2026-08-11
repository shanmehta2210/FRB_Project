# Legacy Survey v2 EXP catalog audit

Catalog: `catalog/LS_catalog_v2_fullsky_exp.csv` (Tractor `type=EXP` only, full LS footprint).

| Metric | Value |
|--------|-------|
| Rows (finite mag + b/a) | 2,000,000 |
| RA span | 0.000 … 360.000 |
| Dec span | -89.557 … 84.767 |
| median modelMag_r | 22.642 |
| median b/a | 0.4447 |
| mean b/a | 0.4447 |

## Axis ratio

All galaxies are already EXP, so shape is taken from Tractor ellipticity:

$$|e| = \sqrt{e_1^2 + e_2^2},\qquad b/a = (1 - |e|) / (1 + |e|)$$

No Hubble transform; no color or `b/a > q_0` cuts in this audit.

## Plots

| File | Content |
|------|---------|
| `mag_histogram_ba.png` | Galaxy counts per 0.5 mag bin (`modelMag_r`; bright open bin `<15`; display through mag ≤ 26) |
| `median_ba_vs_mag.png` | Median b/a per mag bin (same display cut) |
| `ba_vs_mag_scatter.png` | b/a vs mag (log hexbin + thin scatter overlay; mag ≤ 26) |

CSV companions retain the full magnitude range: `mag_histogram_ba.csv`, `median_ba_vs_mag.csv`.

## Mag-cut cos(i) CDFs (`cdfs/`)

Hubble \(\cos(i)\) from \(b/a(e_1,e_2)\) with \(q_0=0.2\) and \(b/a > q_0\). No color cut.

| File | Cut | N | median cos(i) |
|------|-----|---|---------------|
| `mag20.png` | ≤ 20 | 25,274 | 0.496 |
| `mag21.png` | ≤ 21 | 138,605 | 0.455 |
| `mag22.png` | ≤ 22 | 524,666 | 0.425 |
| `overlay.png` | all three | — | — |

```bash
python scripts/plot_ls_v2_mag_cut_cdfs.py
python scripts/plot_ls_v2_mag_cut_cdfs.py --mode strict_scaled
```

## Physical elliptical-disk scaling

Ad-hoc REX stretch lives in [`scaled/`](scaled/) (`b/a≤0.8`, `(b/a)/0.8`).

Physically motivated tracks (Ryden 2004 / Padilla & Strauss 2008):

| Dir | Method |
|-----|--------|
| [`scaled_ryden/`](scaled_ryden/) | Ryden shape fit + Unterborn \((\log q)^2\) edge-on re-add |
| [`scaled_padilla/`](scaled_padilla/) | Joint Padilla shape + \(E_0\) (photometric \(\psi(\theta)\)) |

Plan + REX research: [`REX_AND_ELLIPTICAL_DISK.md`](REX_AND_ELLIPTICAL_DISK.md). Compare CDFs: `scaled_ryden_padilla_cdf_compare.png`.

```bash
python scripts/fit_ls_scaled_elliptical.py --mode both
```

## Regenerate

```bash
python scripts/audit_ls_v2_mag_ba.py
python scripts/audit_ls_v2_mag_ba.py --scaled
```
