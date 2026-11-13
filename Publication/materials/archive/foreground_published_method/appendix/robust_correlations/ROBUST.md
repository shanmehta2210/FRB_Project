# Outlier-robust correlations with axis ratio $q$

Bootstrap: 10,000 resamples (seed 42). Jackknife: leave one host out. Theil-Sen slope in property units per unit $q$, 95% interval. $P(r_S \le -0.44)$: bootstrap fraction at or below the Marnoch et al. value. AD $p$: two-sample Anderson-Darling test of the property below vs above the median $q$ (5,000 permutations).

| property | $N$ | $r_S$ (Fisher 95%) | $r_S$ $p$ | AD $p$ | $r_S$ bootstrap 95% | Kendall $\tau$ ($p$) | $\tau$ bootstrap 95% | jackknife $r_S$ range | most influential host ($r_S$ without) | Theil-Sen slope (95%) | $P(r_S\le-0.44)$ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| DM$_{\rm host}$ (rest frame) | 55 | -0.06 (-0.32, +0.20) | 0.638 | 0.370 | (-0.32, +0.20) | -0.04 (0.63) | (-0.22, +0.13) | -0.12 to -0.04 | 20230907D (-0.12) | -38 (-239, 159) | 0.003 |
| DM$_{\rm host}$, foreground halos removed | 49 | -0.04 (-0.32, +0.24) | 0.764 | 0.675 | (-0.33, +0.25) | -0.03 (0.74) | (-0.23, +0.17) | -0.10 to +0.00 | 20230124A (-0.10) | -51 (-288, 165) | 0.003 |
| DM$_{\rm host}$, foreground excess over random sightlines removed | 49 | -0.08 (-0.35, +0.21) | 0.605 | 0.652 | (-0.35, +0.21) | -0.06 (0.56) | (-0.25, +0.13) | -0.13 to -0.03 | 20230907D (-0.13) | -67.6 (-278, 147) | 0.005 |
| DM$_{\rm host}$, foreground removed (all photo-$z$ galaxies in front) | 49 | -0.07 (-0.35, +0.21) | 0.623 | 0.902 | (-0.35, +0.22) | -0.04 (0.67) | (-0.24, +0.16) | -0.12 to -0.02 | 20231120A (-0.02) | -42.6 (-292, 178) | 0.005 |
| DM$_{\rm host}$, spectroscopic foreground only | 49 | -0.06 (-0.34, +0.23) | 0.684 | 0.408 | (-0.34, +0.23) | -0.04 (0.72) | (-0.23, +0.16) | -0.12 to -0.02 | 20230124A (-0.12) | -32.6 (-243, 147) | 0.004 |
| $\log|{\rm RM}|$ | 41 | -0.02 (-0.32, +0.29) | 0.912 | 0.346 | (-0.31, +0.30) | -0.00 (0.99) | (-0.19, +0.21) | -0.09 to +0.04 | 20230124A (-0.09) | -0.00544 (-1.46, 1.12) | 0.001 |
| $\log\tau_{\rm scatt}$ | 29 | -0.17 (-0.50, +0.21) | 0.385 | 0.643 | (-0.55, +0.27) | -0.11 (0.41) | (-0.40, +0.21) | -0.28 to -0.07 | 20240312D (-0.28) | -1.03 (-2.92, 1.48) | 0.092 |
| Linear polarization | 24 | -0.33 (-0.65, +0.08) | 0.114 | 0.500 | (-0.69, +0.12) | -0.24 (0.10) | (-0.53, +0.08) | -0.45 to -0.25 | 20220914A (-0.45) | -0.467 (-1.08, 0.188) | 0.299 |
| $A_V$, old stars | 30 | -0.46 (-0.70, -0.12) | 0.010 | 0.003 | (-0.67, -0.15) | -0.28 (0.03) | (-0.46, -0.06) | -0.53 to -0.41 | 20220207C (-0.53) | -0.429 (-1.06, -0.0207) | 0.547 |
| DM$_{\rm host}$ (rest frame), m_r < 21 | 45 | -0.09 (-0.38, +0.21) | 0.543 | 0.080 | (-0.38, +0.19) | -0.07 (0.52) | (-0.26, +0.12) | -0.16 to -0.06 | 20230907D (-0.16) | -54 (-273, 168) | 0.008 |
| DM$_{\rm host}$ (rest frame), n < 2 | 46 | -0.04 (-0.32, +0.26) | 0.806 | 0.585 | (-0.32, +0.26) | -0.02 (0.85) | (-0.21, +0.17) | -0.10 to +0.00 | 20230124A (-0.10) | -10.5 (-177, 168) | 0.002 |
| DM$_{\rm host}$ (rest frame), bulge fits removed | 52 | -0.07 (-0.33, +0.21) | 0.637 | 0.252 | (-0.34, +0.20) | -0.04 (0.64) | (-0.22, +0.13) | -0.12 to -0.04 | 20230907D (-0.12) | -38.4 (-254, 168) | 0.003 |
| DM$_{\rm host}$, foreground halos removed, m_r < 21 | 39 | -0.04 (-0.35, +0.28) | 0.814 | 0.467 | (-0.36, +0.28) | -0.03 (0.76) | (-0.26, +0.18) | -0.11 to +0.02 | 20230124A (-0.11) | -41.5 (-256, 173) | 0.007 |
| DM$_{\rm host}$, foreground halos removed, n < 2 | 42 | -0.00 (-0.31, +0.30) | 0.989 | 0.483 | (-0.31, +0.30) | -0.00 (0.99) | (-0.21, +0.21) | -0.06 to +0.05 | 20230124A (-0.06) | -0.495 (-274, 204) | 0.002 |
| DM$_{\rm host}$, foreground halos removed, bulge fits removed | 46 | -0.04 (-0.33, +0.25) | 0.787 | 0.575 | (-0.33, +0.26) | -0.03 (0.77) | (-0.23, +0.17) | -0.10 to +0.01 | 20230124A (-0.10) | -46.9 (-290, 167) | 0.004 |

## Foreground-halo DM over the sample (rest frame of the host)

LS footprint: 49 hosts; outside: 6. Median 58, 84th percentile 195, maximum 500 pc cm$^{-3}$; 27 hosts above 50 and 15 above 100. Median fraction of DM$_{\rm host}$: 0.24. Hosts negative after correction: 6. Spearman of DM$_{\rm fg}$ with $q$: $r_S = +0.14$ ($p = 0.33$; AD $p = 0.76$). Median expected (random-sightline) term: 44 pc cm$^{-3}$.
