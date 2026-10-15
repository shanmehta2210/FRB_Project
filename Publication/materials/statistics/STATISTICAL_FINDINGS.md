# Inclination tests: confirmed sample and the 2024 sample of 23

Tracking note for the HSC comparison. Point estimates are in `TEST_RESULTS.md`. Protocol B error draws: `scripts/run_frb_hsc_error_mc.py` (2,000 draws). The 2024-host subset: `scripts/run_frb_hsc_nature23_tests.py`, raw output `nature23_tests.json`.

## Current sample, after the 20211203C refit

The corrected host passes $m_r < 22$ and $b/a > 0.2$, so the confirmed in-cut sample is 51. Its point values are $q = 0.92$, $\cos i = 0.917$, $i = 23.6^\circ$. It is not one of the $i > 60^\circ$ hosts. All 23 of the 2024 hosts are now in this sample.

### 51 hosts, point estimates

Median $\cos i = 0.610$ (HSC 0.575). Mean $\cos i = 0.585$ (HSC 0.567).

| Test | Pool $p$ | Matched-$N$ $p$ |
|---|---:|---:|
| Anderson–Darling | $\ge 0.25$ | 0.606 |
| KS, face-on, $D^- = 0.115$ | 0.241 | 0.237 |
| Mann–Whitney, face-on | 0.271 | 0.274 |
| Quartile $\chi^2 = 4.14$, counts 11, 10, 19, 11 | 0.247 | 0.254 |
| $i > 60^\circ$, 15/51 vs HSC 0.398 | 0.083 | 0.081 |

$P(K \ge 15 \mid n = 51, p = 4/23) = 0.024$. Against the HSC-scaled 2024 rate 0.138, $p = 0.0030$.

### 51 hosts, Protocol B (2,000 draws; 16/50/84)

Median $\cos i$: 0.594, 0.618, 0.642. Anderson–Darling stays at the $p \ge 0.25$ cap. One-sided KS $p$: 0.10, 0.24, 0.42. One-sided Mann–Whitney $p$: 0.16, 0.27, 0.40. Quartile $\chi^2$ $p$: 0.22, 0.52, 0.80. $N(i > 60^\circ)$: 15, 17, 18. $P(K \ge k)$ vs $4/23$: 0.0013, 0.0035, 0.020. Vs the scaled rate: $7 \times 10^{-5}$, $2.5 \times 10^{-4}$, 0.0024.

### Original 23, point estimates

Median $\cos i = 0.610$. Mean $\cos i = 0.579$.

| Test | Pool $p$ | Matched-$N$ $p$ |
|---|---:|---:|
| Anderson–Darling | $\ge 0.25$ | 0.930 |
| KS, face-on, $D^- = 0.106$ | 0.554 | 0.557 |
| Mann–Whitney, face-on | 0.417 | 0.423 |
| Quartile $\chi^2 = 1.52$, counts 6, 4, 8, 5 (expected 5.75) | 0.677 | 0.733 |
| $i > 60^\circ$, 7/23 vs HSC 0.398 | 0.243 | 0.244 |

$P(K \ge 7 \mid n = 23, p = 4/23) = 0.090$. Against the scaled rate, $p = 0.032$.

### Original 23, Protocol B (2,000 draws; 16/50/84)

Median $\cos i$: 0.563, 0.597, 0.629. Anderson–Darling stays capped. One-sided KS $p$: 0.28, 0.47, 0.66. One-sided Mann–Whitney $p$: 0.29, 0.41, 0.54. Quartile $\chi^2$ $p$: 0.41, 0.68, 0.92. $N(i > 60^\circ)$: 7, 8, 9. $P(K \ge k)$ vs $4/23$: 0.011, 0.035, 0.090. Vs HSC, $P(K \le k)$: 0.24, 0.40, 0.56.

The sections below are the previous frozen sample of 50, which excluded 20211203C.

## Samples and quantity

Working sample at the time of the sections below: 50 confirmed hosts with $m_r < 22$ and $b/a > 0.2$. No host has $m_r = 22$ (closest 21.98 and 21.97), so $m_r < 22$ and $m_r \le 22$ are the same cut. Table 1 lists 56: the magnitude sample after one host with no reliable magnitude and three morphology failures. Six of those 56 have $b/a \le 0.2$, where the Hubble formula is undefined ($q_0 = 0.2$). Those six are not in these tests.

Control: Kawinwanichakij et al. (2021) HSC-SSP disks, `goodfits_flag = 1`, photometric $r \le 22$, strict $0.4 < n < 1.5$, $b/a > 0.2$, $N = 24{,}450$. Shapes are $i$-band; the magnitude cut is catalog $r$. HSC is the control because the shapes are PSF-convolved single free-$n$ Sérsic fits (Lenstronomy). FRB shapes are GALFIT. Same model class, different code.

$$
\cos^2 i = \frac{(b/a)^2 - 0.2^2}{1 - 0.2^2}.
$$

$i = 0^\circ$ is face-on ($\cos i = 1$). $i = 90^\circ$ is edge-on ($\cos i = 0$). A uniform-in-$\cos i$ population has median $i = 60^\circ$. The HSC pool does not: median $\cos i = 0.575$ ($i = 54.9^\circ$), and 39.8% have $i > 60^\circ$.

FRB values are winning-leg point estimates unless a section says Protocol B. Protocol B draws $q \sim \mathcal{N}(q, \sigma_q)$ with $\sigma_q = \sqrt{q_{\rm err}^2 + \sigma_{q,{\rm sky}}^2}$, clips to $0.2 < q \le 1$, converts with Hubble, then adds $\mathcal{N}(0, 5^\circ)$ on $i$ and clips to $[0^\circ, 90^\circ]$. Median $\sigma_q$ on the 50 is 0.011. HSC is never smeared.

## How to read the probabilities

Pool $p$: the 50 (or 22) against all 24,450 HSC galaxies. Matched-$N$ $p$: fraction of 10,000 draws of the same size, without replacement, whose statistic against the pool is at least as large as the FRB statistic. Prefer matched-$N$. It is the 2024 procedure (they drew 23 from SDSS, 10,000 times).

The alternative has one direction: fewer edge-on hosts, meaning larger $\cos i$. Anderson–Darling has no one-sided form here. Kolmogorov–Smirnov (`alternative='less'`: FRB CDF below HSC) and Mann–Whitney (`alternative='greater'`) are one-sided. Mann–Whitney face-on $p = 0.35$ and edge-on $p = 0.65$ on the 50. They are not the same number.

SciPy Anderson–Darling probabilities are capped at 0.25 when the samples agree. A printed 0.25 means $p \ge 0.25$. Quote $p$, not the statistic (it is negative when the samples agree).

## Confirmed 50, point estimates

| | FRB | HSC |
|---|---:|---:|
| $N$ | 50 | 24,450 |
| median $\cos i$ | 0.610 | 0.575 |
| mean $\cos i$ | 0.579 | 0.567 |

FRB minus HSC: median $+0.035$, mean $+0.012$. Slightly face-on.

| Test | Statistic | Pool $p$ | Matched-$N$ $p$ |
|---|---:|---:|---:|
| Anderson–Darling | $-0.418$ (do not quote) | $\ge 0.25$ | 0.574 |
| KS, face-on | $D^- = 0.107$ | 0.294 | 0.286 |
| KS, opposite | $D^+ = 0.073$ | 0.556 | — |
| KS, two-sided | $D = 0.107$ | 0.573 | 0.569 |
| Mann–Whitney, face-on | $U = 630{,}824$ | 0.348 | 0.352 |
| Mann–Whitney, opposite | | 0.652 | — |
| Mann–Whitney, two-sided | rank-biserial $r = -0.032$ | 0.695 | — |
| Quartile $\chi^2$ | 4.56 | 0.207 | 0.220 |

$D^- + D^+ = 0.181$ equals the Kuiper $V$. Rank-biserial $r > 0$ would mean FRB ranks toward smaller $\cos i$ (more edge-on). $r = -0.032$ is a slight face-on rank shift. $U/(50 \times 24{,}450) = 0.516$, against 0.5 for identical ranks.

HSC quartile edges, edge-on to face-on. Each bin is 25% of HSC, so 12.5 of 50 if the samples match. These are not the 2024 degree boundaries.

| Bin | $\cos i$ | $i$ | FRB count | Expected |
|---|---|---|---:|---:|
| 1 | 0.020–0.389 | $88.8^\circ$–$67.1^\circ$ | 11 | 12.5 |
| 2 | 0.389–0.575 | $67.1^\circ$–$54.9^\circ$ | 10 | 12.5 |
| 3 | 0.575–0.755 | $54.9^\circ$–$41.0^\circ$ | 19 | 12.5 |
| 4 | 0.755–0.999 | $41.0^\circ$–$2.6^\circ$ | 10 | 12.5 |

Counts 11, 10, 19, 10. The edge-on bin is full. The excess is in bin 3. Bins 1+2 are the HSC edge-on half, $i > 54.9^\circ$: 21 hosts, against 25 expected ($\chi^2 = 1.28$, $p = 0.258$ for the two-bin split).

Other splits, low $\cos i$ to high, all with expected counts $\ge 10$:

| Bins | Expected | Counts | $\chi^2$ | $p$ |
|---:|---:|---|---:|---:|
| 2 | 25 | 21, 29 | 1.28 | 0.258 |
| 3 | 16.7 | 13, 21, 16 | 1.96 | 0.375 |
| 4 | 12.5 | 11, 10, 19, 10 | 4.56 | 0.207 |
| 5 | 10 | 8, 8, 12, 15, 7 | 4.60 | 0.331 |

### Fixed angle versus the quartiles

The 21 is $i > 54.9^\circ$ (HSC median). The 15 is $i > 60^\circ$. In 2024 those cuts coincided because $60^\circ$ was the SDSS median. Here $i > 60^\circ$ is the most edge-on 39.8% of HSC, not the most edge-on half. All 11 hosts in bin 1 have $i > 67.1^\circ$. Four more of the 15 sit in bin 2 between $60^\circ$ and $67.1^\circ$. The other 6 in bin 2 sit between $54.9^\circ$ and $60^\circ$ and are in the 21 only.

| Cut | FRB | FRB fraction | HSC fraction | Expected | One-sided pool $p$ | Matched-$N$ $p$ |
|---|---:|---:|---:|---:|---:|---:|
| $i > 45^\circ$ | 34/50 | 0.680 | 0.681 | 34.0 | 0.547 | — |
| $i > 60^\circ$ | 15/50 | 0.300 | 0.398 | 19.9 | 0.100 | 0.098 |
| $i > 75^\circ$ | 4/50 | 0.080 | 0.102 | 5.1 | 0.414 | — |

The $i > 60^\circ$ alternative is fewer edge-on hosts than HSC, so the probability is $P(K \le 15)$.

### 2024 numbers these are compared with

2024 sample: 23 hosts, SDSS control, $m_r < 21$, parametric inclinations (conclusion unchanged with non-parametric isophotes). Quartiles of SDSS $\cos i$:

| Bin | $i$ | Hosts |
|---|---|---:|
| 1 | $0^\circ$–$42^\circ$ | 11 |
| 2 | $42^\circ$–$60^\circ$ | 8 |
| 3 | $60^\circ$–$76^\circ$ | 3 |
| 4 | $76^\circ$–$90^\circ$ | 1 |

Bins 3+4 are $i \ge 60^\circ$: 4/23 = 0.174 on the parametric angles. The second angle in their table gives 2/23. Published bins use the parametric set. Text: 48% have $i < 42^\circ$ (11/23 = 47.8%); one host has $i > 75^\circ$. Anderson–Darling $p = 0.0003$. Mann–Whitney $p = 0.00012$. They read 11 hosts in bin 1 as an all-sky rate low by about a factor of two (44 expected, 23 observed), Poisson, ignoring survey differences.

Same observed count of 15 with $i > 60^\circ$, three nulls:

| Null | Rate | Expected out of 50 | Probability |
|---|---:|---:|---:|
| HSC fraction at $i > 60^\circ$ | 0.398 | 19.9 | $P(K \le 15) = 0.100$ |
| 2024 absolute rate, 4/23 | 0.174 | 8.7 | $P(K \ge 15) = 0.020$ |
| 2024 rate scaled to HSC, $(0.174/0.5)\times 0.398$ | 0.138 | 6.9 | $P(K \ge 15) = 0.0024$ |

The most face-on HSC quartile ($i < 41.0^\circ$) contains 10/50 against 12.5 expected. There is no 2024-style pile-up in that bin.

## Confirmed 50, Protocol B errors

2,000 draws. Entries are the 16th, 50th, and 84th percentiles. HSC fixed.

| Quantity | Point estimate | 16%, median, 84% |
|---|---:|---|
| median $\cos i$ | 0.610 | 0.591, 0.613, 0.636 |
| mean $\cos i$ | 0.579 | 0.564, 0.577, 0.588 |
| Anderson–Darling $p$ | $\ge 0.25$ | $\ge 0.25$ across the central 68% |
| one-sided KS $p$ | 0.29 | 0.13, 0.28, 0.47 |
| one-sided Mann–Whitney $p$ | 0.35 | 0.22, 0.34, 0.48 |
| quartile $\chi^2$ $p$ | 0.21 | 0.18, 0.51, 0.83 |
| $N(i > 60^\circ)$ | 15 | 15, 17, 18 |
| quartile counts | 11, 10, 19, 10 | median 11, 11, 16, 12 |

Fraction of draws with $p < 0.05$: Anderson–Darling 0.003, KS 0.018, Mann–Whitney 0.000, $\chi^2$ 0.044. Fraction with $N(i > 60^\circ) \ge 15$: 0.90. Fifteen is the low end of the error range.

| Null for the varying count $k$ | 16%, median, 84% of the tail probability |
|---|---|
| $P(K \le k)$ at the HSC rate 0.398 | 0.100, 0.25, 0.35 |
| $P(K \ge k)$ at 4/23 | 0.0013, 0.0035, 0.020 |
| $P(K \ge k)$ at the scaled rate 0.138 | $7\times 10^{-5}$, $2.5\times 10^{-4}$, 0.0024 |

The upper end of the 2024 rows is the point-estimate result, recovered when a draw lands on $k = 15$.

## 2024 hosts, remeasured

Source list: `Archive/original paper/frb_sample.txt` (23 names). Tests use this pipeline’s winning-leg $\cos i$ for those hosts that sit in the confirmed 50, against the same HSC pool. Matched draws are of 22.

FRB 20211203C is confirmed but is not in the 50. The apparent magnitude fails $m_r < 22$: GALFIT `mag` and `mag_final` are both 22.17, and the Legacy Surveys DR10 catalog magnitude is 22.23 at 0.044 arcsec from the fit. $b/a = 0.44$, Hubble $i = 66.4^\circ$. The verification note calls it faint and barely resolved ($R_e/\mathrm{FWHM} \approx 0.69$, SNR $\lesssim 1$), with a possible sky problem ($\Delta q_{\rm sky} \approx 0.06$).

The pipeline did not fit the galaxy used by Gordon et al. (2023) and the 2024 paper. Gordon's position, 13:38:15.00, $-31$:22:48.2, is 0.42 arcsec from AstroPath candidate 143 (aperture magnitude 20.10) and 9.2 arcsec from candidate 136, the source that was fit. Candidate 136 has Legacy Surveys DR10 model magnitude 22.23. The stored coordinate is labeled host and sits on candidate 136, so AstroPath's prior was already there: posterior 0.55 for 136, about 0 for 143, and 0.45 for an unseen host. The 19.64 in the 2024 table is Gordon's Galactic-extinction-corrected SOAR $r$ for their galaxy, which passes $m_r < 21$. The "AP agrees" note is AstroPhot agreeing with GALFIT on candidate 136, not a match to the Gordon host. The tests below use the other 22 hosts.

| | These 22 | HSC |
|---|---:|---:|
| median $\cos i$ | 0.604 | 0.575 |
| mean $\cos i$ | 0.564 | 0.567 |

Median still slightly face-on ($+0.029$). Mean is $0.003$ below HSC. Rank-biserial $r = +0.015$: a slight edge-on rank shift, the opposite sign from the 50, and the same size as noise.

| Test | Statistic | Pool $p$ | Matched-$N$ $p$ |
|---|---:|---:|---:|
| Anderson–Darling | $-0.797$ (do not quote) | $\ge 0.25$ | 0.855 |
| KS, face-on | $D^- = 0.095$ | 0.634 | 0.634 |
| Mann–Whitney, face-on | $U = 264{,}936$, $r = +0.015$ | 0.548 | 0.542 |
| Quartile $\chi^2$ | 2.00 | 0.572 | 0.617 |
| $i > 60^\circ$ vs HSC | 7/22 = 0.318 vs 0.398 | 0.296 | 0.295 |

Quartile counts, edge-on to face-on: 6, 4, 8, 4. Expected 5.5 per bin. The edge-on bin is 6 against 5.5. The most face-on bin is 4 against 5.5.

HSC expectation at $i > 60^\circ$ is $22 \times 0.398 = 8.8$. Observed 7.

| Null | Rate | Expected out of 22 | Probability |
|---|---:|---:|---:|
| HSC | 0.398 | 8.8 | $P(K \le 7) = 0.296$ |
| 2024 absolute rate, 4/23 | 0.174 | 3.8 | $P(K \ge 7) = 0.073$ |
| Scaled onto HSC | 0.138 | 3.0 | $P(K \ge 7) = 0.025$ |

Protocol B, 2,000 draws of these 22. HSC fixed.

| Quantity | Point estimate | 16%, median, 84% |
|---|---:|---|
| median $\cos i$ | 0.604 | 0.553, 0.583, 0.614 |
| Anderson–Darling $p$ | $\ge 0.25$ | $\ge 0.25$ across the central 68% |
| one-sided KS $p$ | 0.63 | 0.34, 0.56, 0.76 |
| one-sided Mann–Whitney $p$ | 0.55 | 0.40, 0.54, 0.67 |
| quartile $\chi^2$ $p$ | 0.57 | 0.28, 0.65, 0.82 |
| $N(i > 60^\circ)$ | 7 | 7, 8, 9 |
| quartile counts | 6, 4, 8, 4 | median 5, 6, 7, 4 |

Fraction of draws with $p < 0.05$: Anderson–Darling 0.0005, KS 0.003, Mann–Whitney 0. $P(K \le k)$ against HSC: 0.30, 0.46, 0.63. $P(K \ge k)$ against 4/23: 0.008, 0.027, 0.073.

The remeasured 22 do not reproduce the 2024 Anderson–Darling or Mann–Whitney result. They sit on the HSC distribution. The edge-on count is 7 rather than 4, and the probability of that under the published 4/23 rate is 0.073 (0.025 after scaling to the HSC base rate). With $N = 22$ that exclusion is weaker than the same comparison on the 50.

### Host-by-host

$i$ is this pipeline’s Hubble inclination from the winning-leg $q$. $i_{2024}$ is the first (parametric) angle in the 2024 table. $\sigma_q$ is Protocol B on $q$ only, before the $5^\circ$ floor.

| FRB | $q$ | $\cos i$ | $i$ (deg) | $i_{2024}$ (deg) | $\sigma_q$ |
|---|---:|---:|---:|---:|---:|
| 20171020A | 0.400 | 0.354 | 69.3 | 69 | 0.0006 |
| 20180924B | 0.680 | 0.663 | 48.4 | 34 | 0.0112 |
| 20190102C | 0.347 | 0.290 | 73.2 | 50 | 0.1311 |
| 20190608B | 0.700 | 0.685 | 46.8 | 37 | 0.0009 |
| 20190714A | 0.400 | 0.354 | 69.3 | 51 | 0.0103 |
| 20191001A | 0.540 | 0.512 | 59.2 | 53 | 0.0006 |
| 20200906A | 0.340 | 0.281 | 73.7 | 63 | 0.0006 |
| 20210320C | 0.700 | 0.685 | 46.8 | 28 | 0.0200 |
| 20210410D | 0.390 | 0.342 | 70.0 | 41 | 0.0867 |
| 20210807D | 0.710 | 0.696 | 45.9 | 41 | 0.0075 |
| 20211127I | 0.973 | 0.972 | 13.5 | 21 | 0.0044 |
| 20211212A | 0.830 | 0.822 | 34.7 | 33 | 0.0002 |
| 20220207C | 0.230 | 0.116 | 83.3 | 80 | 0.0022 |
| 20220307B | 0.630 | 0.610 | 52.4 | 44 | 0.0335 |
| 20220310F | 0.690 | 0.674 | 47.6 | 39 | 0.0302 |
| 20220319D | 0.800 | 0.791 | 37.8 | 25 | 0.0151 |
| 20220509G | 0.544 | 0.516 | 58.9 | 60 | 0.0013 |
| 20220825A | 0.620 | 0.599 | 53.2 | 46 | 0.0505 |
| 20220912A | 0.630 | 0.610 | 52.4 | 43 | 0.0424 |
| 20220914A | 0.510 | 0.479 | 61.4 | 45 | 0.0108 |
| 20220920A | 0.860 | 0.854 | 31.4 | 27 | 0.0003 |
| 20221012A | 0.530 | 0.501 | 59.9 | 45 | 0.0003 |

Seven of these 22 have $i > 60^\circ$: 20171020A, 20190102C, 20190714A, 20200906A, 20210410D, 20220207C, 20220914A. In the 2024 parametric bins the edge-on four were 20171020A, 20200906A, 20220207C, and 20220509G. 20220509G is now $58.9^\circ$.

## Appendix tests on the 50

Same point estimates. Redundant with the four main tests. Cramér–von Mises pool $p$ is unusable when $N_{\rm HSC} \gg 50$ (printed $5.7 \times 10^{-6}$); the matched value is the valid one.

| Test | Statistic | Pool $p$ | Matched-$N$ $p$ |
|---|---:|---:|---:|
| KS, two-sided | $D = 0.107$ | 0.573 | 0.569 |
| Cramér–von Mises | $T = 2.141$ | unusable | 0.503 |
| Kuiper | $V = 0.181$ | — | 0.346 |
| Mann–Whitney, two-sided | $r = -0.032$ | 0.695 | — |
| Brunner–Munzel | $W = -0.420$ | 0.676 | — |
| Epps–Singleton | $W = 4.32$ | 0.364 | — |
| Welch $t$ | $t = 0.384$ | 0.702 | — |
| Wasserstein-1 | $W_1 = 0.0334$ | — | 0.498 |

## Limits

Point-estimate $p$-values do not include $\sigma_q$ or the $5^\circ$ floor. The Protocol B sections do, on the FRB side only. Matched-$N$ for the error draws was not run; those percentiles are pool probabilities of each noisy realization. The $5^\circ$ term dominates $\sigma_q = 0.011$ once both are on $i$. Quartile edges in degrees belong to this HSC pool, not to the 2024 edges at $42^\circ$, $60^\circ$, and $76^\circ$.
