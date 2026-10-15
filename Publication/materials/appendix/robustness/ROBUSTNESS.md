# Appendix robustness checks (51 hosts vs HSC)

Source: `scripts/appendix_robustness.py` (draws per null: 10,000; error-MC realizations: 2,000; matched-N draws use the main-text streams, seed 42 for AD and 43 for MW / $\chi^2$, so the adopted row equals the main-text run; error-MC seed 20241025). Hosts: `plots/plots_null/v2/frb_vs_hsc_confirmed50/confirmed50_q.csv` (51, winning-leg $q$). Control: Kawinwanichakij et al. (2021) HSC-SSP parent, `goodfits_flag = 1`, $r \le 22$, $0.4 < n < 1.5$, $q > 0.2$ unless stated. Inclinations from $\cos^2 i = (q^2 - q_0^2)/(1 - q_0^2)$, $q_0 = 0.2$.

$p_{\rm m}$ = fraction of random control subsamples of the FRB sample size whose statistic against the full pool is at least as extreme as the FRB value (two-sided AD statistic; two-sided MW $|U - nm/2|$; one-sided MW $U$ in the edge-on-deficit direction, i.e. FRB $\cos i$ larger; quartile $\chi^2$ with bins at the pool quartiles). Quartile counts are listed face-on to edge-on.

## A1. The control sample is not uniform in $\cos i$

`figA_control_hist.pdf` (final version: `scripts/appendix_figures.py`): HSC disk pool ($N = 24,450$, $r \le 22$, $0.4 < n < 1.5$, $q > 0.2$). Median $\cos i = 0.575$ ($i = 54.9^\circ$); fraction with $i > 60^\circ$: 39.8% (random orientation: $60^\circ$ and 50%). KS against uniform in $\cos i$: $D = 0.158$, $p = 0.0e+00$. The pool is depleted at both ends: at $\cos i \to 0$ because the $q > q_0$ cut, the magnitude limit and dust remove the most inclined disks, and because the Hubble formula maps $q$ just above $q_0$ onto a narrow range of $\cos i$; at $\cos i \to 1$ because measured axis ratios rarely reach $q = 1$ (intrinsic disk ellipticity, non-axisymmetric structure, PSF and noise). Tests must therefore be made against this pool, not against a uniform distribution, and the quartile bins of Figure 2 are defined on it.

## A2. Sérsic-index window of the control

FRB sample fixed (51). Pool re-selected; quartile bins recomputed from each pool.

| sample | $N_{\rm FRB}$ | $N_{\rm pool}$ | median $\cos i$ FRB / pool | quartile counts | AD $p_{\rm m}$ | MW 2-sided $p_{\rm m}$ | MW 1-sided $p_{\rm m}$ | $\chi^2$ $p_{\rm m}$ |
|---|---|---|---|---|---|---|---|---|
| $0.4 < n < 1.5$ (adopted) (median $n$ = 0.92) | 51 | 24,450 | 0.610 / 0.575 | 11,19,10,11 | 0.61 | 0.56 | 0.27 | 0.25 |
| $0.4 < n < 2.0$ (median $n$ = 1.05) | 51 | 31,212 | 0.610 / 0.582 | 11,18,10,12 | 0.71 | 0.65 | 0.33 | 0.40 |
| $0.4 < n < 2.5$ (median $n$ = 1.16) | 51 | 36,700 | 0.610 / 0.587 | 11,18,10,12 | 0.79 | 0.76 | 0.38 | 0.41 |
| $n < 2.5$ (median $n$ = 1.14) | 51 | 37,670 | 0.610 / 0.586 | 11,18,10,12 | 0.79 | 0.72 | 0.37 | 0.40 |
| all $n$ (no morphology cut) (median $n$ = 1.56) | 51 | 52,682 | 0.610 / 0.617 | 10,14,15,12 | 0.77 | 0.66 | 0.66 | 0.79 |

## A3. Magnitude limit of the control

| sample | $N_{\rm FRB}$ | $N_{\rm pool}$ | median $\cos i$ FRB / pool | quartile counts | AD $p_{\rm m}$ | MW 2-sided $p_{\rm m}$ | MW 1-sided $p_{\rm m}$ | $\chi^2$ $p_{\rm m}$ |
|---|---|---|---|---|---|---|---|---|
| control $r \le 22$, all 51 hosts | 51 | 24,450 | 0.610 / 0.575 | 11,19,10,11 | 0.61 | 0.56 | 0.27 | 0.25 |
| control $r \le 21$, all 51 hosts | 51 | 7,423 | 0.610 / 0.592 | 10,19,10,12 | 0.71 | 0.80 | 0.40 | 0.24 |
| both samples $r \le 21$ | 41 | 7,423 | 0.642 / 0.592 | 9,16,7,9 | 0.62 | 0.56 | 0.28 | 0.21 |
| control $r \le 20$, all 51 hosts | 51 | 1,426 | 0.610 / 0.598 | 10,19,10,12 | 0.79 | 1.00 | 0.50 | 0.22 |
| both samples $r \le 20$ | 29 | 1,426 | 0.610 / 0.598 | 7,10,6,6 | 0.81 | 0.74 | 0.36 | 0.70 |

## A4. Magnitude-matched control draws

Each null draw takes, from every $r$ bin, the same number of control galaxies as there are FRB hosts in that bin (without replacement); $p_{\rm m}$ is the fraction of such magnitude-matched draws that differ from the full pool at least as much as the FRB hosts. The FRB statistic itself is unchanged (51 hosts vs the full pool); only the null changes.

| $r$ bin | FRB hosts | control galaxies | hosts |
|---|---|---|---|
|  r <= 18 | 10 | 1 | 20171020A, 20190608B, 20210807D, 20211127I, 20211212A, 20220319D, 20220509G, 20220725A, 20231120A, 20240210A |
| 18 < r <= 19 | 11 | 90 | 20191001A, 20220207C, 20220920A, 20221012A, 20221106A, 20230124A, 20230307A, 20231123B, 20231226A, 20240213A, 20240318A |
| 19 < r <= 20 | 7 | 1,335 | 20210320C, 20220307B, 20220825A, 20220912A, 20220914A, 20230628A, 20230907D |
| 20 < r <= 21 | 13 | 5,997 | 20180924B, 20190714A, 20200906A, 20211203C, 20220310F, 20220418A, 20230526A, 20230613A, 20230626A, 20231220A, 20240215A, 20240229A, 20240304A |
| 21 < r <= 22 | 10 | 17,027 | 20181112A, 20190102C, 20190523A, 20210410D, 20220501C, 20221219A, 20230902A, 20231020B, 20240119A, 20240208A |

The HSC pool is essentially empty brighter than $r = 18$ (bright galaxies saturate / are excluded from the Kawinwanichakij et al. fits), so a bin-by-bin match of the ten $r \le 18$ hosts is impossible. Two runs are made: the bright end merged into a single $r \le 19$ bin (21 hosts drawn from 91 control galaxies, so those draws are strongly repeated), and hosts with $r > 19$ only.

| sample | $N_{\rm FRB}$ | $N_{\rm pool}$ | median $\cos i$ FRB / pool | quartile counts | AD $p_{\rm m}$ | MW 2-sided $p_{\rm m}$ | MW 1-sided $p_{\rm m}$ | $\chi^2$ $p_{\rm m}$ |
|---|---|---|---|---|---|---|---|---|
| all 51, bright end merged into one $r \le 19$ bin | 51 | 24,450 | 0.610 / 0.575 | 11,19,10,11 | 0.65 | 0.55 | 0.16 | 0.26 |
| hosts with $r > 19$ only | 30 | 24,450 | 0.610 / 0.575 | 5,14,5,6 | 0.42 | 0.56 | 0.34 | 0.05 |

## A5. Colour / star-formation cut on the control

HSC $g - i$ is observed-frame (no $k$-correction); `quiescent_flag` is the Kawinwanichakij et al. (2021) UVJ-based flag. The FRB sample is unchanged.

| sample | $N_{\rm FRB}$ | $N_{\rm pool}$ | median $\cos i$ FRB / pool | quartile counts | AD $p_{\rm m}$ | MW 2-sided $p_{\rm m}$ | MW 1-sided $p_{\rm m}$ | $\chi^2$ $p_{\rm m}$ |
|---|---|---|---|---|---|---|---|---|
| star-forming only (`quiescent_flag = 0`) | 51 | 23,331 | 0.610 / 0.577 | 11,18,10,12 | 0.62 | 0.58 | 0.29 | 0.40 |
| blue, $g - i < 1.2$ | 51 | 11,720 | 0.610 / 0.583 | 11,18,10,12 | 0.73 | 0.63 | 0.32 | 0.40 |
| red, $g - i \ge 1.2$ | 51 | 12,730 | 0.610 / 0.569 | 11,19,10,11 | 0.47 | 0.48 | 0.24 | 0.25 |

## B. Gini–$M_{20}$ (statmorph on the GALFIT cutouts, PSF and masks)

`figB_gini_m20.pdf` (final version: `scripts/appendix_figures.py`). Lotz et al. (2008) boundaries: mergers $G > -0.14\,M_{20} + 0.33$; E/S0/Sa $G > 0.14\,M_{20} + 0.80$ (below the merger line); Sb–Irr otherwise.

All 51: late 47, early 3, merger 1. Reliable only (`flag = 0`, $\langle S/N\rangle_{\rm pix} \ge 2.5$, $R_e/{\rm FWHM} \ge 1.0$; $N = 30$): late 26, early 3, merger 1. Marked unreliable: 21 (flag $\ne 0$: 7; low S/N: 4; marginally resolved: 16).

Hosts outside the late-type region:

| FRB | class | reliable | $G$ | $M_{20}$ |
|---|---|---|---|---|
| 20190608B | early | True | 0.544 | -1.84 |
| 20220319D | early | True | 0.587 | -2.40 |
| 20220509G | merger | True | 0.539 | -1.03 |
| 20231120A | early | True | 0.534 | -2.01 |

Hosts described as quiescent in the literature:

| FRB | class | reliable | $G$ | $M_{20}$ | $F(G,M_{20})$ | note |
|---|---|---|---|---|---|---|
| 20190523A | late | False | 0.349 | -1.57 | -1.15 | flag; Re/FWHM<1.0 |
| 20210807D | late | True | 0.456 | -1.80 | -0.46 |  |
| 20220509G | merger | True | 0.539 | -1.03 | -0.58 |  |
| 20221012A | late | True | 0.466 | -1.69 | -0.49 |  |

Per-host values: `gini_m20_51.csv`.

## C1. Intrinsic thickness $q_0$

$\cos^2 i = (q^2 - q_0^2)/(1 - q_0^2)$ and the cut $q > q_0$ applied to both samples. Hosts with $q \le 0.2$ that pass the magnitude and verification cuts re-enter at lower $q_0$.

| $q_0$ | hosts added | hosts removed |
|---|---|---|
| 0.10 | 20240201A ($q=0.18$) | — |
| 0.15 | 20240201A ($q=0.18$) | — |
| 0.20 | — | — |
| 0.25 | — | 20220207C ($q=0.23$), 20231120A ($q=0.25$), 20240213A ($q=0.24$) |

| $q_0$ | $N_{\rm FRB}$ | $N_{\rm pool}$ | median $\cos i$ FRB / pool | quartile counts | AD $p_{\rm m}$ | MW 2-sided $p_{\rm m}$ | MW 1-sided $p_{\rm m}$ | $\chi^2$ $p_{\rm m}$ |
|---|---|---|---|---|---|---|---|---|
| 0.10 | 52 | 24,976 | 0.625 / 0.585 | 11,19,10,12 | 0.63 | 0.54 | 0.27 | 0.30 |
| 0.15 | 52 | 24,870 | 0.619 / 0.579 | 11,19,10,12 | 0.65 | 0.56 | 0.29 | 0.30 |
| 0.20 | 51 | 24,450 | 0.610 / 0.575 | 11,19,10,11 | 0.61 | 0.55 | 0.27 | 0.25 |
| 0.25 | 48 | 23,699 | 0.614 / 0.572 | 11,18,10,9 | 0.31 | 0.30 | 0.15 | 0.27 |

## C2. Measurement uncertainties (Monte Carlo)

2000 realizations. Each host: $q \sim \mathcal{N}(q, \sigma_q)$ truncated to $q_0 < q \le 1$ (redrawn, not clipped; median $\sigma_q = 0.011$), converted to $i$, plus $\mathcal{N}(0, 5^\circ)$ reflected at $0^\circ$ and $90^\circ$. The same $5^\circ$ reflected scatter is applied to every control galaxy in every realization, so both samples are treated alike. $p_{\rm m}$ in each realization is read against a null of 10,000 draws of 51 from scattered control pools.

| quantity | 16% | median | 84% | fraction of realizations with $p_{\rm m} < 0.05$ |
|---|---|---|---|---|
| AD $p_{\rm m}$ | 0.34 | 0.57 | 0.80 | 0.0% |
| MW two-sided $p_{\rm m}$ | 0.32 | 0.51 | 0.75 | 0.0% |
| MW one-sided $p_{\rm m}$ | 0.16 | 0.25 | 0.38 | 0.1% |
| quartile $\chi^2$ $p_{\rm m}$ | 0.20 | 0.51 | 0.83 | 2.6% |
| median $\cos i$ | 0.596 | 0.620 | 0.642 | |
| $N(i > 60^\circ)$ | 15 | 17 | 18 | |

Largest fraction of realizations in which any test reaches $p_{\rm m} < 0.05$: **2.6%**. Median quartile counts (face-on to edge-on): 13, 16, 12, 11.
