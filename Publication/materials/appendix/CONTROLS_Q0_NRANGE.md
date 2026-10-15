# Appendix controls: intrinsic thickness $q_0$ and the HSC Sérsic-index window

Source: `scripts/appendix_controls_tables.py`, which reuses the loaders of `scripts/plot_paper_inclination_figures.py` (figures `figA_q0`, `figA_nrange`). Point estimates only (winning-leg $q$, no error smearing). HSC control: Kawinwanichakij et al. parent sample, `goodfits_flag = 1`, $r \le 22$. One-sided $p$ tests the direction claimed in 2024 (FRB hosts more edge-on, i.e. smaller $\cos i$). Matched-$N$ KS $p$ is the fraction of 10,000 $N$-sized HSC subsamples with $D \ge D_{\rm obs}$.

## 1. Intrinsic axis ratio $q_0$

$\cos^2 i = (q^2 - q_0^2)/(1 - q_0^2)$ applied to both samples with the same $q_0$; the $q > q_0$ cut is re-applied to both, so the FRB sample size changes with $q_0$ (hosts with $q \le 0.2$ that pass the magnitude and verification cuts re-enter at lower $q_0$).

| $q_0$ | $N_{\rm FRB}$ | $N_{\rm HSC}$ | hosts added / removed vs 0.2 | median $\cos i$ FRB / HSC | $f(i>60^\circ)$ FRB / HSC | KS $D$ | KS $p$ (2-sided) | KS $p$ (edge-on) | KS $p$ matched-$N$ | MWU $p$ (2-sided) | MWU $p$ (edge-on) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.10 | 52 | 24,976 | + 20240201A (0.18) | 0.625 / 0.585 | 0.288 / 0.377 | 0.114 | 0.478 | 0.243 | 0.469 | 0.535 | 0.268 |
| 0.15 | 52 | 24,870 | + 20240201A (0.18) | 0.619 / 0.579 | 0.308 / 0.390 | 0.112 | 0.502 | 0.255 | 0.490 | 0.569 | 0.284 |
| 0.20 | 51 | 24,450 | — | 0.610 / 0.575 | 0.294 / 0.398 | 0.115 | 0.476 | 0.241 | 0.462 | 0.538 | 0.269 |
| 0.25 | 48 | 23,699 | − 20220207C (0.23), 20231120A (0.25), 20240213A (0.24) | 0.614 / 0.572 | 0.312 / 0.402 | 0.138 | 0.293 | 0.147 | 0.283 | 0.312 | 0.156 |

## 2. Sérsic-index window of the HSC control

FRB sample fixed (51, $q_0 = 0.2$). The HSC pool is re-selected with different $n$ windows (the adopted window is $0.4 < n < 1.5$).

| $n$ window | $N_{\rm HSC}$ | median $n$ | median $\cos i$ HSC | $f(i>60^\circ)$ HSC | KS $D$ | KS $p$ (2-sided) | KS $p$ (edge-on) | KS $p$ matched-$N$ | MWU $p$ (2-sided) | MWU $p$ (edge-on) |
|---|---|---|---|---|---|---|---|---|---|---|
| $0.4 < n < 1.5$ | 24,450 | 0.92 | 0.575 | 0.398 | 0.115 | 0.476 | 0.241 | 0.465 | 0.543 | 0.271 |
| $0.5 < n < 1.2$ | 18,200 | 0.85 | 0.572 | 0.405 | 0.119 | 0.428 | 0.216 | 0.415 | 0.493 | 0.247 |
| $0.8 < n < 1.2$ | 10,652 | 0.98 | 0.581 | 0.393 | 0.108 | 0.558 | 0.286 | 0.561 | 0.630 | 0.315 |
| $0.4 < n < 2.0$ | 31,212 | 1.05 | 0.582 | 0.386 | 0.105 | 0.589 | 0.302 | 0.580 | 0.650 | 0.325 |

FRB sample: median $\cos i$ = 0.610, $f(i>60^\circ)$ = 0.294.
