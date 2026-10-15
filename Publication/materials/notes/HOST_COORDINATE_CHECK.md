# Host-coordinate check, first two papers

Compared the galaxy actually fit to the host or FRB coordinate in the paper. A match means the published position lands on the source that was fit, or on no other catalog source. Separation is in arcseconds.

## Gordon et al. 2023 (17 hosts in this pipeline)

| FRB | Separation | Our LS $r$ | Their magnitude | Reading |
|---|---:|---:|---:|---|
| 20211203C | 9.15 | 22.23 | 19.64 | Wrong galaxy. A source at 0.42 arcsec from their position, aperture mag 20.10, was not the one fit. Refit below. |
| 20190711A | 6.61 | 18.33 | 23.54 | Wrong object. Already rejected: the fit is on a brighter source; their host is $r = 23.54$. |
| 20210320C | 0.07 | 19.77 | 19.47 | Same galaxy. The Gordon et al. 2025 host coordinate lands on the source that was fit. |
| 20210807D | 4.78 | 17.79 | 17.17 | Same galaxy. Half-light radius is about 24 arcsec, so 5 arcsec is inside it. |
| 20191001A | 3.17 | 18.44 | 18.36 | Same galaxy. Magnitudes agree. |
| 20190608B | 2.92 | 17.56 | 17.41 | Same galaxy. Magnitudes agree. |
| 20190611B | 2.08 | 22.94 | 22.15 | Same nearest source. |
| ten others | $\le 1.13$ | | | Position and magnitude agree. |

The ten are 20200906A, 20210410D, 20180924B, 20190714A, 20200430A, 20220105A, 20190102C, 20181112A, 20211127I, and 20211212A.

## Law et al. 2024, Table 2 (10 hosts in this pipeline)

All ten land within 1.12 arcsec of the published host coordinate: 20220207C, 20220307B, 20220310F, 20220319D, 20220418A, 20220509G, 20220825A, 20220914A, 20220920A, 20221012A. 20220506D is in that table and not in this pipeline.

Two large galaxies are brighter in their Pan-STARRS magnitude than in Legacy Surveys (20220207C, 16.15 vs 18.71; 20220319D, 13.25 vs 16.20). The coordinates still match.

## Remaining papers

Sharma et al. 2024 host coordinates match the fitted galaxy for every overlap except 20221116A. The matches are 20230124A, 20230628A, 20231123B, 20230307A, 20230626A, 20230712A, 20231120A, and 20221219A (0.90 arcsec). 20221116A is a different object: the Sharma host is 8.9 arcsec from the source that was fit, on the other side of the burst, and nothing in the Pan-STARRS catalog sits on the published coordinate. That fit is already rejected (empty stamp, no reference magnitude).

Gordon et al. 2025 host coordinates match, all within 1.3 arcsec: 20220725A, 20221106A, 20230526A, 20230708A, 20230902A, 20231226A, 20240201A, 20240210A, 20240304A, 20240310A, 20240318A, and the corrected 20211203C (0.07 arcsec). 20220918A is the exception. The cutout is centered 0.5 arcsec from their host, but the only detection is 3.5 arcsec away and the half-light radius is 0.4 arcsec. 20240208A is not in that table.

Single-host and small-sample papers, separation of the fitted source from the published host or burst:

| FRB | Separation | Coordinate used | Reading |
|---|---:|---|---|
| 20220717A | 0.08 | Rajwade et al. 2024 PATH host | Same galaxy. |
| 20220912A | 0.08 | Ravi et al. 2023 PSO host | Same galaxy. |
| 20171020A | 0.73 | Li et al. 2023, ESO 601-G036 | Same galaxy. Half-light radius is 6.4 arcsec. |
| 20230930A | 0.00 | Anna-Thomas et al. 2025 host | Same galaxy. |
| 20230613A | 0.72 | Pastor-Marazuela et al. 2026 burst | Same galaxy. |
| 20230907D | 0.93 | Pastor-Marazuela et al. 2026 burst | Same galaxy. |
| 20231020B | 0.17 | Pastor-Marazuela et al. 2026 burst | Same galaxy. |
| 20220501C | 0.57 | Shannon et al. 2025 burst | Same galaxy. |
| 20240208A | 1.31 | Shannon et al. 2025 burst | Same nearest galaxy. Half-light radius is 0.7 arcsec. |
| 20230708A | 0.21 | Muller et al. 2026 / Gordon host | Same galaxy. |
| 20240215A | 0.53 | Hussaini et al. 2025 burst | Same galaxy. |
| 20240213A | 2.7 | Hussaini et al. 2025 burst | Nearest detection. Half-light radius is 1.7 arcsec. |
| 20240229A | 2.2 | Hussaini et al. 2025 burst | Nearest detection. Half-light radius is 0.9 arcsec. |
| 20240119A | 4.4 | Hussaini et al. 2025 burst | Nearest detection, and no other source lies on the burst. Half-light radius is 1.0 arcsec, so the burst sits well off the centroid. |

20240304B has no pipeline fit. Caleb et al. 2025 report no ground-based host to $R = 26.2$. Verdi et al. (in prep.) have no published coordinate: 20230913, 20240104A, 20240203, 20250518. 20210214G is a signal row, not a literature host. Ravi et al. 2019 for 20190523A was not re-extracted in this pass.

## 20211203C refit

The old fit is in `pipeline_scripts/Output/20211203C_all_wronghost`. The new run centers on 13:38:15.00, $-31$:22:48.2. AstroPath selects source 143 with posterior 0.86. Legacy Surveys DR10 model $r = 20.06$. GALFIT: $m_r = 20.11$, $b/a = 0.92$, $n = 0.56$, $R_e = 0.48$ arcsec, $i = 23.6^\circ$, $\chi^2_\nu = 0.73$. PSFEx measures the seeing as 0.93 arcsec. The 1.83 arcsec value is the width of the collapsed PSF stamp in the field-depth helper, not the image seeing. $R_e = 0.48$ arcsec is about half that PSF, so the source is only marginally resolved and $i = 23.6^\circ$ is not a clean conflict with the 2024 value of $7^\circ$. Sky shifts $q$ by 0.0006. AstroPhot gives $q = 0.91$. GALFIT minus the LS magnitude is 0.04 mag. The source now passes $m_r < 22$ and $b/a > 0.2$. It is not yet in the frozen list of 50.
