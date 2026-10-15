# Appendix: non-parametric morphology (statmorph) of the 51 hosts

Source: `scripts/statmorph_51_full.py` (statmorph v0.5 run on the adopted GALFIT cutout of each host, `host_cutout.fits`, with `host_sigma.fits` as weight map, the PSFEx prototype as PSF and the GALFIT star mask; segmentation = 1.5σ detection on a Gaussian-smoothed image, central segment kept). Table generator: `scripts/appendix_statmorph_section.py`. Per-host numbers: `statmorph_51.csv`; LaTeX: `statmorph_51_table.tex`.

## A. Independent inclination check

statmorph fits its own single Sérsic model (different optimiser, PSF convolution and masking from GALFIT); AstroPhot is the third, fully independent fitter already used in the paper. All three are converted to inclination with $\cos^2 i = (q^2 - q_0^2)/(1 - q_0^2)$, $q_0 = 0.2$.

| pair | $N$ | median $\Delta i$ | RMS | $1.4826\,$MAD | $|\Delta i| > 10^\circ$ |
|---|---|---|---|---|---|
| statmorph Sérsic − GALFIT | 51 | +0.16° | 4.11° | 0.85° | 3 |
| AstroPhot − GALFIT | 50 | +0.54° | 5.12° | 0.80° | 3 |
| AstroPhot − statmorph | 50 | +0.35° | 6.32° | 1.51° | 3 |
| statmorph − GALFIT, flag = flag_sersic = 0 | 41 | +0.20° | 3.60° | 0.75° | 1 |

Per-host standard deviation of the three Sérsic inclinations: median 0.69°, mean 1.78°.

The RMS of 4–5° between independent Sérsic codes on the same images is what motivates the $5^\circ$ inclination floor; the core scatter is below 1° and the RMS is set entirely by the marginally resolved tail ($R_e/{\rm FWHM} < 1$):

| FRB | $R_e$/FWHM | $i$ GALFIT | $i$ statmorph | $i$ AstroPhot | note |
|---|---|---|---|---|---|
| 20190523A | 0.55 | 45 | 44 | 64 | faint ($m_r = 21.7$), $R_e$/FWHM = 0.55 |
| 20210410D | 0.92 | 70 | 90 | 75 | GALFIT $n$ at the 6 cap; statmorph $q = 0.18$ |
| 20220501C | 0.64 | 59 | 48 | 84 | statmorph `flag = 1`; faint ($m_r = 21.8$) |
| 20220509G | 4.58 | 59 | 71 | 71 | statmorph `flag_sersic = 2`; see B |

Moment-based ellipticities (statmorph `ellipticity_asymmetry`, isophote fits at $1R_e$) are not PSF-corrected and sit $\sim 8$–$11^\circ$ low systematically; they are not used.

## B. Disk classification

Adopted metric: the Gini–$M_{20}$ bulge statistic $F(G, M_{20})$ of Rodriguez-Gomez et al. (2019), $F < 0$ = disk-dominated side. Supporting metrics: the Lotz et al. (2008) $G$–$M_{20}$ classes, concentration $C$ and the Sérsic index.

| metric | disk side | of |
|---|---|---|
| $F(G,M_{20}) < 0$ | 48 | 51 |
| Lotz+2008 late (Sb–Irr) | 47 | 51 (early 3, merger 1) |
| $C < 3.5$ | 49 | 51 |
| statmorph $n < 2.5$ | 44 | 51 |
| GALFIT $n < 2.5$ | 41 | 51 |

Hosts not on the disk side of $F$, with the literature morphology:

| FRB | $F$ | Lotz | $C$ | $n$ GALFIT / statmorph | literature |
|---|---|---|---|---|---|
| 20190608B | +0.00 | early | 2.99 | 3.1 / 2.4 | barred spiral, FRB on a spiral arm (Chittidi+2021; Mannings+2021 HST) |
| 20220319D | +0.61 | early | 3.99 | 6.0 / 16.2 | IRAS 02044+7048, barred spiral SBa, FRB on the eastern arm (Ravi+2025; Law+2024) |
| 20231120A | +0.08 | early | 3.54 | 1.9 / 1.8 | 2MFGC 07440: 2MASS Flat Galaxy Catalogue, edge-on by selection ($a/b \ge 3$; Mitronova+2004); spiral arms (Sharma+2024); edge-on disk (Leung+2025) |
| 20220509G | -0.58 | merger ($A = 0.60$) | 1.85 | 0.5 / 5.2 | quiescent early-type in Abell 2311 (Sharma+2023; Connor+2023); literature contradicts a disk |

All three $F \ge 0$ hosts are spirals in the literature (two barred spirals, one catalogued flat edge-on galaxy whose concentration is raised by projection). The one host whose disk nature is contradicted rather than supported is FRB 20220509G: statmorph puts it on the disk side ($F = -0.58$) but with an unconstrained Sérsic shape (GALFIT $n$ at the 0.5 floor, statmorph $n = 5.2$, `flag_sersic = 2`, $A = 0.60$), and Sharma+2023 / Connor+2023 classify it as a quiescent early-type galaxy. It was part of the 2024 sample ($i = 60^\circ$).

## C. Per-host table

| FRB | $m_r$ | $R_e$/FWHM | $q$ GALFIT | $i$ GALFIT | $q$ sm | $i$ sm | $i$ AP | $n$ GALFIT | $n$ sm | $G$ | $M_{20}$ | $F$ | $C$ | $A$ | Lotz | flags |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 20171020A | 15.10 | 11.48 | 0.40 | 69 | 0.40 | 70 | 69 | 0.54 | 0.49 | 0.44 | -1.25 | -0.89 | 2.27 | 0.05 | late | 0/0 |
| 20180924B | 20.24 | 0.92 | 0.68 | 48 | 0.67 | 49 | 50 | 3.91 | 3.32 | 0.53 | -1.74 | -0.12 | 2.85 | 0.02 | late | 0/0 |
| 20181112A | 21.98 | 0.71 | 0.86 | 31 | 0.86 | 31 | 32 | 1.02 | 0.80 | 0.51 | -1.59 | -0.32 | 2.30 | 0.01 | late | 1/0 |
| 20190102C | 21.67 | 0.96 | 0.35 | 73 | 0.42 | 68 | 70 | 0.50 | 0.35 | 0.39 | -1.55 | -0.94 | 2.18 | -0.15 | late | 1/0 |
| 20190523A | 21.72 | 0.55 | 0.72 | 45 | 0.73 | 44 | 64 | 0.62 | 0.18 | 0.35 | -1.57 | -1.15 | 2.23 | 0.01 | late | 1/0 |
| 20190608B | 17.46 | 1.95 | 0.70 | 47 | 0.69 | 48 | 47 | 3.13 | 2.43 | 0.54 | -1.84 | +0.00 | 2.99 | 0.03 | early | 0/0 |
| 20190714A | 20.24 | 1.04 | 0.40 | 69 | 0.41 | 69 | 71 | 0.50 | 0.44 | 0.42 | -1.50 | -0.82 | 2.31 | -0.02 | late | 0/0 |
| 20191001A | 18.41 | 1.56 | 0.54 | 59 | 0.54 | 59 | 60 | 0.51 | 0.51 | 0.43 | -1.55 | -0.76 | 2.29 | 0.11 | late | 0/0 |
| 20200906A | 20.13 | 1.56 | 0.34 | 74 | 0.34 | 74 | 75 | 0.65 | 0.62 | 0.41 | -1.42 | -0.92 | 2.43 | 0.10 | late | 0/0 |
| 20210320C | 19.80 | 0.96 | 0.70 | 47 | 0.71 | 46 | 48 | 1.13 | 0.95 | 0.44 | -1.63 | -0.65 | 2.44 | -0.04 | late | 0/0 |
| 20210410D | 21.12 | 0.92 | 0.39 | 70 | 0.18 | 90 | 75 | 6.00 | 3.98 | 0.53 | -1.65 | -0.22 | 2.52 | -0.01 | late | 0/0 |
| 20210807D | 15.93 | 3.86 | 0.71 | 46 | 0.68 | 48 | 45 | 3.58 | 1.45 | 0.46 | -1.80 | -0.46 | 2.63 | -0.08 | late | 0/0 |
| 20211127I | 14.25 | 5.87 | 0.97 | 14 | 0.96 | 17 | 15 | 3.14 | 2.17 | 0.45 | -1.87 | -0.42 | 2.75 | 0.06 | late | 0/0 |
| 20211203C | 20.11 | 0.48 | 0.92 | 24 | 0.92 | 24 | 25 | 0.56 | 0.51 | 0.45 | -1.60 | -0.61 | 2.30 | -0.00 | late | 0/0 |
| 20211212A | 16.46 | 3.37 | 0.83 | 35 | 0.82 | 35 | 35 | 1.09 | 0.94 | 0.41 | -1.69 | -0.75 | 2.49 | 0.17 | late | 0/0 |
| 20220207C | 18.12 | 5.86 | 0.23 | 83 | 0.22 | 85 | 84 | 0.50 | 0.34 | 0.41 | -1.29 | -1.04 | 2.15 | -0.18 | late | 0/0 |
| 20220307B | 19.87 | 1.41 | 0.63 | 52 | 0.62 | 54 | 53 | 0.58 | 0.44 | 0.39 | -1.44 | -1.04 | 2.21 | -0.21 | late | 0/0 |
| 20220310F | 20.83 | 1.14 | 0.69 | 48 | 0.69 | 48 | 49 | 1.50 | 1.19 | 0.43 | -1.65 | -0.68 | 2.42 | -0.00 | late | 1/0 |
| 20220319D | 16.08 | 6.66 | 0.80 | 38 | 0.82 | 36 | 37 | 6.00 | 16.24 | 0.59 | -2.40 | +0.61 | 3.99 | -0.05 | early | 0/2 |
| 20220418A | 20.88 | 1.33 | 0.76 | 42 | 0.76 | 42 | 42 | 0.50 | 0.37 | 0.36 | -1.37 | -1.22 | 2.07 | -0.03 | late | 0/0 |
| 20220501C | 21.80 | 0.64 | 0.55 | 59 | 0.68 | 48 | 84 | 1.31 | 1.46 | 0.43 | -1.43 | -0.82 | 2.20 | -0.29 | late | 1/0 |
| 20220509G | 16.66 | 4.58 | 0.54 | 59 | 0.38 | 71 | 71 | 0.50 | 5.24 | 0.54 | -1.03 | -0.58 | 1.85 | 0.60 | merger | 0/2 |
| 20220725A | 17.76 | 1.76 | 0.60 | 55 | 0.60 | 55 | 55 | 0.96 | 0.96 | 0.48 | -1.66 | -0.45 | 2.49 | 0.05 | late | 0/0 |
| 20220825A | 19.98 | 1.36 | 0.62 | 53 | 0.59 | 55 | 54 | 0.68 | 0.43 | 0.40 | -1.26 | -1.08 | 2.20 | -0.12 | late | 0/0 |
| 20220912A | 19.52 | 1.31 | 0.63 | 52 | 0.62 | 53 | 52 | 0.86 | 0.68 | 0.45 | -1.57 | -0.67 | 2.31 | -0.15 | late | 0/0 |
| 20220914A | 19.65 | 1.25 | 0.51 | 61 | 0.50 | 62 | 62 | 2.55 | 1.74 | 0.45 | -1.66 | -0.61 | 2.40 | -0.01 | late | 0/0 |
| 20220920A | 18.13 | 1.90 | 0.86 | 31 | 0.86 | 31 | 32 | 1.12 | 1.04 | 0.47 | -1.69 | -0.44 | 2.50 | 0.03 | late | 0/0 |
| 20221012A | 18.96 | 1.47 | 0.53 | 60 | 0.54 | 60 | 61 | 2.01 | 1.64 | 0.47 | -1.69 | -0.49 | 2.53 | -0.00 | late | 0/0 |
| 20221106A | 18.09 | 2.46 | 0.67 | 49 | 0.67 | 49 | 49 | 1.79 | 1.42 | 0.50 | -1.77 | -0.26 | 2.71 | 0.01 | late | 0/0 |
| 20221219A | 21.97 | 1.25 | 0.50 | 62 | 0.43 | 67 | 69 | 0.50 | 0.01 | 0.39 | -1.45 | -1.03 | 2.09 | -0.21 | late | 1/1 |
| 20230124A | 18.58 | 1.00 | 0.90 | 26 | 0.91 | 26 | 27 | 1.93 | 1.66 | 0.46 | -1.64 | -0.54 | 2.46 | 0.02 | late | 0/0 |
| 20230307A | 18.18 | 6.58 | 0.80 | 38 | 0.79 | 39 | 37 | 6.00 | 2.78 | 0.49 | -1.81 | -0.26 | 2.90 | 0.03 | late | 0/0 |
| 20230526A | 20.09 | 0.39 | 0.77 | 40 | 0.77 | 41 | 45 | 1.85 | 2.46 | 0.46 | -1.60 | -0.58 | 2.40 | -0.00 | late | 0/0 |
| 20230613A | 20.31 | 1.11 | 0.84 | 34 | 0.84 | 34 | 34 | 0.62 | 0.59 | 0.42 | -1.60 | -0.76 | 2.23 | -0.02 | late | 0/0 |
| 20230626A | 20.00 | 0.53 | 0.73 | 44 | 0.73 | 45 | 45 | 6.00 | 6.18 | 0.48 | -1.63 | -0.47 | 2.50 | -0.02 | late | 0/0 |
| 20230628A | 19.56 | 1.18 | 0.73 | 44 | 0.73 | 45 | 45 | 1.04 | 0.92 | 0.43 | -1.63 | -0.68 | 2.33 | -0.02 | late | 0/0 |
| 20230902A | 21.54 | 0.42 | 0.63 | 52 | 0.64 | 51 | 54 | 0.88 | 0.80 | 0.47 | -1.55 | -0.56 | 2.32 | -0.04 | late | 1/0 |
| 20230907D | 19.41 | 1.22 | 0.85 | 33 | 0.85 | 32 | 32 | 1.58 | 1.29 | 0.46 | -1.72 | -0.49 | 2.59 | -0.02 | late | 0/0 |
| 20231020B | 21.75 | 0.59 | 0.59 | 55 | 0.60 | 55 | 57 | 0.67 | 0.64 | 0.46 | -1.59 | -0.59 | 2.29 | -0.03 | late | 0/0 |
| 20231120A | 15.47 | 7.87 | 0.25 | 81 | 0.25 | 81 | 81 | 1.91 | 1.84 | 0.53 | -2.01 | +0.08 | 3.54 | 0.08 | early | 0/0 |
| 20231123B | 18.59 | 3.15 | 0.38 | 71 | 0.38 | 71 | 71 | 0.70 | 0.65 | 0.41 | -1.59 | -0.80 | 2.51 | 0.05 | late | 0/0 |
| 20231220A | 20.22 | 1.06 | 0.26 | 81 | 0.20 | 90 | — | 0.57 | 0.51 | 0.44 | -1.52 | -0.73 | 2.37 | -0.04 | late | 0/0 |
| 20231226A | 18.54 | 2.10 | 0.76 | 42 | 0.75 | 43 | 42 | 0.65 | 0.58 | 0.41 | -1.50 | -0.88 | 2.22 | -0.02 | late | 0/0 |
| 20240119A | 21.44 | 1.36 | 0.75 | 42 | 0.79 | 39 | 42 | 0.50 | 0.03 | 0.35 | -1.23 | -1.39 | 1.87 | -0.06 | late | 0/0 |
| 20240208A | 21.80 | 0.65 | 0.40 | 69 | 0.42 | 68 | 69 | 1.15 | 1.02 | 0.48 | -1.64 | -0.45 | 2.51 | 0.01 | late | 0/0 |
| 20240210A | 14.50 | 7.26 | 0.59 | 55 | 0.59 | 55 | 55 | 1.00 | 0.97 | 0.48 | -1.54 | -0.50 | 2.73 | 0.31 | late | 0/0 |
| 20240213A | 18.35 | 2.52 | 0.24 | 82 | 0.24 | 82 | 83 | 0.51 | 0.51 | 0.45 | -1.54 | -0.66 | 2.52 | 0.11 | late | 0/0 |
| 20240215A | 20.06 | 0.81 | 0.48 | 64 | 0.48 | 64 | 64 | 2.03 | 1.46 | 0.45 | -1.53 | -0.69 | 2.36 | 0.04 | late | 0/0 |
| 20240229A | 20.84 | 0.56 | 0.63 | 52 | 0.58 | 56 | 54 | 6.00 | 5.94 | 0.48 | -1.68 | -0.42 | 2.47 | -0.04 | late | 0/2 |
| 20240304A | 20.81 | 0.93 | 0.66 | 50 | 0.67 | 49 | 51 | 1.24 | 0.95 | 0.44 | -1.52 | -0.73 | 2.35 | -0.05 | late | 0/0 |
| 20240318A | 18.47 | 2.53 | 0.43 | 67 | 0.42 | 68 | 68 | 0.70 | 0.64 | 0.43 | -1.62 | -0.72 | 2.40 | 0.03 | late | 0/0 |
