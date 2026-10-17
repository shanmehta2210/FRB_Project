# Appendix controls: intrinsic thickness $q_0$ and the HSC Sérsic-index window

Source: `scripts/appendix_controls_tables.py`, which reuses the loaders of `scripts/plot_paper_inclination_figures.py` (figures `figA_q0`, `figA_nrange`). Point estimates only (winning-leg $q$, no error smearing). HSC control: Kawinwanichakij et al. parent sample, `goodfits_flag = 1`, $r \le 22$. One-sided $p$ tests the direction claimed in 2024 (FRB hosts more edge-on, i.e. smaller $\cos i$). Matched-$N$ KS $p$ is the fraction of 10,000 $N$-sized HSC subsamples with $D \ge D_{\rm obs}$.

## 1. Intrinsic axis ratio $q_0$

$\cos^2 i = (q^2 - q_0^2)/(1 - q_0^2)$ applied to both samples with the same $q_0$; the $q > q_0$ cut is re-applied to both, so the FRB sample size changes with $q_0$ (hosts with $q \le 0.2$ that pass the magnitude and verification cuts re-enter at lower $q_0$).

| $q_0$ | $N_{\rm FRB}$ | $N_{\rm HSC}$ | hosts added / removed vs 0.2 | median $\cos i$ FRB / HSC | $f(i>60^\circ)$ FRB / HSC | KS $D$ | KS $p$ (2-sided) | KS $p$ (edge-on) | KS $p$ matched-$N$ | MWU $p$ (2-sided) | MWU $p$ (edge-on) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.10 | 58 | 24,976 | + 20240201A (0.18) | 0.640 / 0.585 | 0.293 / 0.377 | 0.119 | 0.356 | 0.180 | 0.356 | 0.322 | 0.161 |
| 0.15 | 58 | 24,870 | + 20240201A (0.18) | 0.634 / 0.579 | 0.310 / 0.390 | 0.117 | 0.378 | 0.191 | 0.370 | 0.348 | 0.174 |
| 0.20 | 57 | 24,450 | — | 0.642 / 0.575 | 0.298 / 0.398 | 0.119 | 0.365 | 0.184 | 0.358 | 0.334 | 0.167 |
| 0.25 | 53 | 23,699 | − 20220207C (0.23), 20231120A (0.25), 20231210F (0.22), 20240213A (0.24) | 0.653 / 0.572 | 0.302 / 0.402 | 0.150 | 0.164 | 0.082 | 0.154 | 0.115 | 0.057 |

## 2. Sérsic-index window of the HSC control

FRB sample fixed (57, $q_0 = 0.2$). The HSC pool is re-selected with different $n$ windows (the adopted window is $0.4 < n < 1.5$).

| $n$ window | $N_{\rm HSC}$ | median $n$ | median $\cos i$ HSC | $f(i>60^\circ)$ HSC | KS $D$ | KS $p$ (2-sided) | KS $p$ (edge-on) | KS $p$ matched-$N$ | MWU $p$ (2-sided) | MWU $p$ (edge-on) |
|---|---|---|---|---|---|---|---|---|---|---|
| $0.4 < n < 1.5$ | 24,450 | 0.92 | 0.575 | 0.398 | 0.119 | 0.365 | 0.184 | 0.357 | 0.336 | 0.168 |
| $0.5 < n < 1.2$ | 18,200 | 0.85 | 0.572 | 0.405 | 0.124 | 0.322 | 0.162 | 0.325 | 0.300 | 0.150 |
| $0.8 < n < 1.2$ | 10,652 | 0.98 | 0.581 | 0.393 | 0.112 | 0.442 | 0.224 | 0.433 | 0.406 | 0.203 |
| $0.4 < n < 2.0$ | 31,212 | 1.05 | 0.582 | 0.386 | 0.109 | 0.471 | 0.239 | 0.455 | 0.418 | 0.209 |

FRB sample: median $\cos i$ = 0.642, $f(i>60^\circ)$ = 0.298.
