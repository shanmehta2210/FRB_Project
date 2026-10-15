# statmorph vs GALFIT, 51-host working sample

Inclinations from cos^2 i = (q^2 - q0^2)/(1 - q0^2), q0 = 0.2, q = 1 - ellipticity.

- Sersic-fit ellipticity, all hosts: N=51, mean=+0.81, median=+0.16, RMS=4.11, std=4.07, 1.4826*MAD=0.85, |d|>10 deg: 3
- Sersic-fit ellipticity, flag=0 & flag_sersic=0 (41): N=41, mean=+0.98, median=+0.20, RMS=3.60, std=3.51, 1.4826*MAD=0.75, |d|>10 deg: 1
- Moment (asymmetry-centre) ellipticity, all: N=51, mean=-11.66, median=-10.82, RMS=17.53, std=13.22, 1.4826*MAD=12.55, |d|>10 deg: 31
- Moment ellipticity, flag=0 (44): N=44, mean=-10.92, median=-10.55, RMS=15.10, std=10.55, 1.4826*MAD=11.40, |d|>10 deg: 25

## Classification counts

- Lotz+2008 G-M20: {'late': 47, 'early': 3, 'merger': 1}
- F(G,M20) bulge statistic < 0 (disk side): 48 / 51
- statmorph Sersic n < 2.5: 44; GALFIT n < 2.5: 41
- Concentration C < 3.5 (Conselice 2003 late-type side): 49

## Cross-code inclination scatter (same images, q0 = 0.2)

| pair | N | median | RMS | 1.4826 MAD | outliers > 10 deg |
|---|---|---|---|---|---|
| statmorph Sersic - GALFIT | 51 | +0.16 | 4.11 | 0.85 | 20210410D (+20.0), 20220509G (+12.0), 20220501C (-10.6) |
| AstroPhot - GALFIT | 50 | +0.54 | 5.12 | 0.80 | 20220501C (+25.1), 20190523A (+19.7), 20220509G (+12.1) |
| AstroPhot - statmorph | 50 | +0.35 | 6.32 | 1.51 | |
| isophote (1 Re, PSF-uncorrected) - GALFIT | 43 | -7.35 | 14.5 | 11.0 | systematic: moments are rounded by the PSF |

Moment-based (PSF-uncorrected) ellipticities are biased low by ~11 deg and are not a fair comparator; the three Sersic fits are.

## Hosts not on the disk side of F(G, M20) (Rodriguez-Gomez et al. 2019; F < 0 = disk)

| FRB | F | Lotz 2008 | C | n (GALFIT / statmorph) | literature |
|---|---|---|---|---|---|
| 20190608B | 0.00 | early (on the line) | 2.99 | 3.1 / 2.4 | barred spiral, FRB on arm (Chittidi+2021 HST+VLT; Mannings+2021) |
| 20231120A | +0.08 | early | 3.54 | 1.9 / 1.8 | 2MFGC 07440, 2MASS Flat Galaxy Catalogue (a/b >= 3 selection, Mitronova+2004); spiral arms (Sharma+2024); edge-on disk (Leung+2025) |
| 20220319D | +0.61 | early | 3.99 | 6.0 (cap) / 16 | IRAS 02044+7048, barred spiral SBa, FRB on eastern arm (Ravi+2025); face-on barred spiral (Law+2024) |

Lotz merger class: 20220509G (G = 0.54, M20 = -1.03, A = 0.60). F = -0.58 (disk side), but the literature calls this host a quiescent early-type elliptical (Sharma+2023; Connor+2023, Abell 2311 member). GALFIT n sits at the 0.5 floor and statmorph n = 5.2 with flag_sersic = 2, so the Sersic shape is not constrained. This is the one host where the disk assumption is contradicted by the literature, not supported by it.

Flagged inclination outliers: 20210410D (GALFIT q = 0.39 with n at the 6 cap, statmorph q = 0.18, AstroPhot i = 75; Re/FWHM = 0.92), 20220501C (Re/FWHM = 0.64, statmorph flag = 1, three codes give 48/59/84 deg), 20190523A (AstroPhot 64 vs 45; Re/FWHM = 0.55).
