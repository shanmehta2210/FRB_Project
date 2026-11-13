# Outlier-robust correlations with axis ratio $q$

Bootstrap: 10,000 resamples (seed 42). Jackknife: leave one host out. Theil-Sen slope in property units per unit $q$, 95% interval. $P(r_S \le -0.44)$: bootstrap fraction at or below the Marnoch et al. value. AD $p$: two-sample Anderson-Darling test of the property below vs above the median $q$ (5,000 permutations).

| property | $N$ | $r_S$ (Fisher 95%) | $r_S$ $p$ | AD $p$ | $r_S$ bootstrap 95% | Kendall $\tau$ ($p$) | $\tau$ bootstrap 95% | jackknife $r_S$ range | most influential host ($r_S$ without) | Theil-Sen slope (95%) | $P(r_S\le-0.44)$ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| DM$_{\rm host}$ (rest frame) | 55 | -0.06 (-0.32, +0.20) | 0.638 | 0.370 | (-0.32, +0.20) | -0.04 (0.63) | (-0.22, +0.13) | -0.12 to -0.04 | 20230907D (-0.12) | -38 (-239, 159) | 0.003 |
| DM$_{\rm host}$, foreground halos removed | 49 | -0.09 (-0.36, +0.20) | 0.560 | 0.349 | (-0.37, +0.20) | -0.06 (0.51) | (-0.26, +0.13) | -0.14 to -0.04 | 20230124A (-0.14) | -42.6 (-275, 136) | 0.005 |
| DM$_{\rm host}$, foreground excess over random sightlines removed | 49 | -0.09 (-0.37, +0.19) | 0.520 | 0.348 | (-0.36, +0.19) | -0.06 (0.57) | (-0.24, +0.14) | -0.15 to -0.05 | 20230124A (-0.15) | -53.9 (-256, 144) | 0.007 |
| DM$_{\rm host}$, foreground removed (all photo-$z$ galaxies in front) | 49 | -0.07 (-0.35, +0.21) | 0.615 | 0.348 | (-0.35, +0.21) | -0.06 (0.56) | (-0.25, +0.14) | -0.13 to -0.03 | 20230124A (-0.13) | -45.4 (-288, 135) | 0.005 |
| DM$_{\rm host}$, spectroscopic foreground only | 49 | -0.06 (-0.33, +0.23) | 0.699 | 0.461 | (-0.34, +0.23) | -0.04 (0.70) | (-0.23, +0.15) | -0.12 to -0.02 | 20230124A (-0.12) | -20.6 (-219, 172) | 0.003 |
| DM$_{\rm host}$, foreground excess removed (published method) | 49 | -0.09 (-0.36, +0.20) | 0.534 | 0.492 | (-0.36, +0.18) | -0.06 (0.55) | (-0.25, +0.13) | -0.15 to -0.05 | 20230124A (-0.15) | -71 (-279, 154) | 0.005 |
| $\log|{\rm RM}|$ | 41 | -0.02 (-0.32, +0.29) | 0.912 | 0.335 | (-0.32, +0.29) | -0.00 (0.99) | (-0.20, +0.21) | -0.09 to +0.04 | 20230124A (-0.09) | -0.00544 (-1.46, 1.12) | 0.001 |
| $\log\tau_{\rm scatt}$ | 29 | -0.17 (-0.50, +0.21) | 0.385 | 0.648 | (-0.55, +0.27) | -0.11 (0.41) | (-0.40, +0.21) | -0.28 to -0.07 | 20240312D (-0.28) | -1.03 (-2.92, 1.48) | 0.092 |
| Linear polarization | 24 | -0.33 (-0.65, +0.08) | 0.114 | 0.492 | (-0.70, +0.12) | -0.24 (0.10) | (-0.53, +0.08) | -0.45 to -0.25 | 20220914A (-0.45) | -0.467 (-1.08, 0.188) | 0.306 |
| $A_V$, old stars | 30 | -0.46 (-0.70, -0.12) | 0.010 | 0.004 | (-0.67, -0.15) | -0.28 (0.03) | (-0.46, -0.06) | -0.53 to -0.41 | 20220207C (-0.53) | -0.429 (-1.06, -0.0207) | 0.543 |
| DM$_{\rm host}$ (rest frame), m_r < 21 | 45 | -0.09 (-0.38, +0.21) | 0.543 | 0.086 | (-0.38, +0.20) | -0.07 (0.52) | (-0.26, +0.13) | -0.16 to -0.06 | 20230907D (-0.16) | -54 (-273, 168) | 0.009 |
| DM$_{\rm host}$ (rest frame), n < 2 | 46 | -0.04 (-0.32, +0.26) | 0.806 | 0.585 | (-0.33, +0.25) | -0.02 (0.85) | (-0.22, +0.17) | -0.10 to +0.00 | 20230124A (-0.10) | -10.5 (-177, 168) | 0.003 |
| DM$_{\rm host}$ (rest frame), bulge fits removed | 52 | -0.07 (-0.33, +0.21) | 0.637 | 0.261 | (-0.34, +0.20) | -0.04 (0.64) | (-0.23, +0.13) | -0.12 to -0.04 | 20230907D (-0.12) | -38.4 (-254, 168) | 0.003 |
| DM$_{\rm host}$, foreground excess over random sightlines removed, m_r < 21 | 39 | -0.11 (-0.41, +0.21) | 0.506 | 0.156 | (-0.43, +0.20) | -0.06 (0.57) | (-0.28, +0.14) | -0.19 to -0.06 | 20230124A (-0.19) | -51 (-263, 171) | 0.021 |
| DM$_{\rm host}$, foreground excess over random sightlines removed, n < 2 | 42 | -0.07 (-0.37, +0.24) | 0.648 | 0.310 | (-0.38, +0.24) | -0.04 (0.72) | (-0.25, +0.17) | -0.14 to -0.02 | 20230124A (-0.14) | -38.2 (-225, 159) | 0.009 |
| DM$_{\rm host}$, foreground excess over random sightlines removed, bulge fits removed | 46 | -0.09 (-0.37, +0.21) | 0.566 | 0.268 | (-0.37, +0.20) | -0.05 (0.63) | (-0.24, +0.15) | -0.15 to -0.05 | 20230124A (-0.15) | -44.1 (-256, 157) | 0.007 |

## Foreground-halo DM over the sample (rest frame of the host)

LS footprint: 49 hosts; outside: 6. Median 23 (published method: 58), 84th percentile 99, maximum 194 pc cm$^{-3}$; 15 hosts above 50 and 8 above 100. Median fraction of DM$_{\rm host}$: 0.08. Hosts negative after removing the excess: 1; after the full subtraction: 2. Spearman of DM$_{\rm fg}$ with $q$: $r_S = +0.15$ ($p = 0.32$; AD $p = 0.46$). Median expected (random-sightline) term: 25 pc cm$^{-3}$.
