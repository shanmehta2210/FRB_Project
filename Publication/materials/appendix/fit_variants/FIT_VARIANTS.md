# Inclination tests with alternative axis-ratio measurements

Control: HSC disk pool ($N = 24,450$). Matched-draw $p_{\rm m}$ from 10,000 draws (main-text seeds).

| measurement | $N$ | median $i$ | quartile counts (face-on to edge-on) | AD $p_{\rm m}$ | MW 2-sided $p_{\rm m}$ | MW 1-sided $p_{\rm m}$ | $\chi^2$ $p_{\rm m}$ |
|---|---|---|---|---|---|---|---|
| GALFIT, PSF-convolved (adopted) | 56 | 51.3° | 14,19,10,13 | 0.54 | 0.46 | 0.23 | 0.4 |
| AstroPhot, PSF-convolved | 56 | 53.2° | 13,18,8,17 | 0.44 | 0.99 | 0.5 | 0.24 |
| GALFIT, no PSF | 56 | 38.0° | 32,13,6,5 | 0.0001 | 0.0001 | 0.0001 | 0.0001 |
| AstroPhot, no PSF | 56 | 38.1° | 31,14,6,5 | 0.0001 | 0.0001 | 0.0001 | 0.0001 |
| Shared with Bhardwaj et al. (2024): this work | 22 | 52.8° | 5,7,4,6 | 0.96 | 0.85 | 0.43 | 0.87 |
| Shared with Bhardwaj et al. (2024): their AutoProf values | 22 | 42.5° | 11,7,2,2 | 0.0044 | 0.0036 | 0.0026 | 0.013 |

## Agreement with the published inclinations of Bhardwaj et al. (2024) (22 shared hosts)

| our measurement | median $i_{\rm ours} - i_{\rm B24}$ | RMS difference |
|---|---|---|
| GALFIT, PSF-convolved | +9.2° | 13.1° |
| GALFIT, no PSF | -0.1° | 3.4° |
| AstroPhot, no PSF | +0.0° | 3.4° |
