# Results-section figures

Built by `scripts/plot_results_figures.py` (seed 20240924; 10,000 control draws). Inputs: winning-leg GALFIT $q$ (`statistics/confirmed50_q.csv`), AstroPhot $q$ from the verification run (`host_properties/seed_56.csv` and `seed_51.csv`, column `q_astrophot`), HSC exponential-disk control pool ($r \le 22$, $q > 0.2$; $N = 24,450$). $\cos^2 i = (q^2 - q_0^2)/(1 - q_0^2)$, $q_0 = 0.2$.

## fig1_cdf

(a) 56 hosts vs control with the 68% band of samples of 56. (b) The 23 hosts of Bhardwaj et al. (2024): published inclinations (grey), GALFIT remeasurement (red), AstroPhot remeasurement (orange), band for samples of 23. The two remeasurements lie on top of each other and inside the band; the published curve does not.

## AstroPhot − GALFIT inclination differences

| sample | $N$ | median | mean | RMS | $1.4826\,$MAD | $|\Delta i| > 10^\circ$ |
|---|---|---|---|---|---|---|
| AstroPhot − GALFIT, the 23 of Bhardwaj+2024 | 23 | +0.61 | +1.19 | 2.98 | 0.73 | 1 |
| AstroPhot − GALFIT, all 56 | 55 | +0.59 | +1.67 | 4.96 | 0.75 | 3 |
| GALFIT (this work) − published 2024 | 23 | +9.43 | +10.04 | 12.94 | 8.17 | 10 |
| AstroPhot (this work) − published 2024 | 23 | +10.40 | +11.22 | 13.94 | 8.68 | 12 |

Hosts of the 23 with $|i_{\rm AstroPhot} - i_{\rm GALFIT}| > 5^\circ$:

| FRB | published | GALFIT | AstroPhot | $R_e$/FWHM |
|---|---|---|---|---|
| 20210410D | 41 | 70.0 | 75.2 | 0.92 |
| 20220509G | 60 | 58.9 | 71.0 | 4.58 |

## fig2_bins

Bins are the quartiles of the control sample's $\cos i$ (edge-on to face-on). Expected 14 per bin.

| bin | $i$ range | FRB count | 68% range of control samples of 56 |
|---|---|---|---|
| 1 | 67°–90° | 13 | 11–17 |
| 2 | 55°–67° | 10 | 11–17 |
| 3 | 41°–55° | 19 | 11–17 |
| 4 | 0°–41° | 14 | 11–17 |

## LaTeX

Figure environments with the updated panel-(b) caption: `results_figures.tex`.
