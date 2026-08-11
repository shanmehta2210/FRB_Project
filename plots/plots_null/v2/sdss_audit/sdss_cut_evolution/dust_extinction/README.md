# SDSS v2 — Unterborn A1 dust correction (mag ≤ 21)

Production null pool: `u−r < 2.3`, `lnLExp > lnLDeV`, `expAB_r > 0.2`.

## Method

Unterborn & Ryden 2008 A1:

$$\Delta m = 1.27\,(\log_{10} q)^2,\qquad m^f = m_{\rm obs} - \Delta m(q)$$

- **raw:** keep `modelMag_r ≤ 21`
- **A1:** keep `m^f ≤ 21` (re-admits edge-ons dust-faded past the observed limit)

## Results

| Mode | N | median cos(i) | median b/a |
|------|--:|--------------:|-----------:|
| raw m_obs≤21 | 159,144 | 0.5843 | 0.6064 |
| A1 m^f≤21 | 182,937 | 0.5423 | 0.5677 |
| A1-added only | 23,793 | 0.2670 | 0.3293 |

Δ median cos(i) (A1 − raw) = **-0.0420**

## Files

- `cdf_raw_vs_a1_mag21.png`
- `ba_hist_raw_vs_a1_added.png`
- `summary_mag21.csv`

```bash
python scripts/plot_sdss_v2_dust_a1_mag21.py
```
