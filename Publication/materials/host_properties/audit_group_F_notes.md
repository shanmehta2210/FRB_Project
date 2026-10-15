# Audit group F: cross-check against Bhardwaj, Lee & Ji (2024) Table 1

Scope: the 23 hosts of Bhardwaj, Lee & Ji 2024 (Nature 634, 1065; arXiv:2408.01876v2). Inputs checked: `Archive/original paper/frb_sample.txt`, `papers/bhardwaj2024_2408.01876.pdf` (the file `papers/bhardwaj2024.pdf` is the same 24-page arXiv v2 build), `audit_seed_picks.csv`, `host_properties_wide.csv`, and the papers Bhardwaj cites: Gordon+2023 (ref. [5]), Law+2024 arXiv:2307.03344 (ref. [15]), Lee-Waddell+2023 ([47]), Ravi+2023 ([48]), plus Bhandari+2022, Sharma+2024, Connor+2025, Ravi+2023b (20220319D), Mahony+2018. Output rows are in `audit_group_F.csv` (159 rows: inclination, $z$, SFR, $\log M_*$, $m_r$, DM$_{\rm host}$ for all 23, plus one transcription row).

## 1. Transcription check (`frb_sample.txt` vs PDF Table 1)

All 23 rows were compared cell by cell (both inclination values, bin, SFR, log sSFR, mass, $z$, DM$_{\rm host}$, $m_r$, reference number).

- One error: row 6, **20220310F**, `log(sSFR)` is `-9.43` in the file; the PDF has $-9.34$. The PDF value is the self-consistent one ($\log(4.25/0.93\times10^{10}) = -9.34$); $-9.43$ is the value of row 1 (20200906A), so this is a copy slip.
- Everything else matches, including the two inclination values per host, the bin labels and the reference numbers.
- Bold rows in Table 1 (Pan-STARRS DR2 images used; all others DESI DR10) are not recoverable from the text file; from the PDF font they are rows 3, 7, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23 (20211127I, 20220914A, 20171020A, 20190714A, 20220207C, 20220509G, 20220825A, 20220912A, 20220307B, 20220319D, 20210807D, 20210320C).

## 2. Per-host comparison with `audit_seed_picks.csv`

### Where Bhardwaj's numbers come from

- CRAFT/MeerKAT hosts (ref. [5], 12 hosts): SFR and mass are Gordon+2023 Table 3 (Prospector, SFR$_{0-100\,\rm Myr}$); $z$ is the Gordon Table 3 *Prospector-derived* redshift (hence 0.3437, 0.2909, 0.2796 instead of the spectroscopic 0.3439, 0.2912, 0.2797); $m_r$ is Gordon Table 1.
- DSA-110 hosts (ref. [15], 9 hosts): $z$, SFR$_{100\,\rm Myr}$ and $\log M_*$ are exactly Law+2024 Table 3 (the arXiv:2307.03344 catalogue). Our picks use the Sharma+2024 (Nature 635, 61) re-fits of the same hosts. The $m_r$ values are **not** from Law+2024 Table 2 nor from Sharma+2024 (see section 4).
- 20171020A (ref. [47]): Lee-Waddell+2023 ($\log M_* = 8.64$, SFR $=0.09$, $z = 0.00867$, $m_r = 14.97$ from their Table 1 WISE$\times$SCOSPZ column).
- 20220912A (ref. [48]): Ravi+2023 ($z = 0.0771$, $\log M_* = 10.0$, $m_r = 19.65$ PS1). The SFR $=2.5$ is not in Ravi+2023.

### Differences and resolutions

| FRB | quantity | Bhardwaj | our pick | resolution |
|---|---|---|---|---|
| 20200906A | mass | $1.33\times10^{10}$ ($\log 10.12$) | 10.37 (Gordon+2023) | **Bhardwaj is wrong for the reference they cite**: Gordon+2023 Table 3 gives $\log M_* = 10.37$ ($2.34\times10^{10}$). $1.33\pm0.37\times10^{10}$ is the CIGALE mass in Bhandari+2022 Table 1. Their log sSFR $-9.43$ mixes the Bhandari mass with the Gordon SFR (with Gordon's mass it would be $-9.68$). Keep 10.37. |
| 20220912A | SFR | 2.5 | blank | Ravi+2023 give only SFR $\gtrsim 0.1$ (H$\alpha$, central arcsec$^2$). 2.5 M$_\odot$ yr$^{-1}$ is the radio-implied SFR of Hewitt+2024 Sec. 5 (hypothetical, if the nearby radio source were star formation). Untraceable to [48]; keep blank / H$\alpha$ limit. |
| 20220307B | $z$ | 0.2481 | 0.2507 (Connor+2025) | Law+2024 (0.248123, pPXF) and Sharma+2024 (0.2481) agree; Connor+2025's compilation table has 0.2507 with no comment. $\Delta z = 0.0026$ is far beyond spectroscopic error. Recommend 0.2481 (host papers; Connor+2025 is a compilation). Marked `conflict`. |
| 9 DSA hosts | SFR, $\log M_*$ | Law+2024 Table 3 | Sharma+2024 ED Table 2 | Both refereed; Sharma+2024 is the later fit by the same group with more photometry. Beyond quoted errors for 20220207C (1.16 vs 0.71), 20220310F (4.25 vs 0.97), 20220319D (0.21 vs 0.42), 20220825A (1.62 vs 7.91), 20220914A (0.72 vs 0.46), 20220920A (4.10 vs 1.62); within errors for 20220307B, 20220509G, 20221012A. Masses differ by $\le 0.09$ dex except 20220307B (10.04 vs 10.14). Marked `conflict`, Sharma+2024 recommended (= current picks). |
| 9 DSA hosts | $m_r$ | 18.24, 19.93, 20.72, 15.39, 17.73, 20.03, 19.77, 18.16, 19.14 | Sharma+2024 catalogue mags | Bhardwaj's values match neither Law+2024 Table 2 (16.15, 20.45, 21.20, 13.25, 16.51, 19.36, 20.05, 19.46, 19.67) nor Sharma+2024 ED Table 1; offsets range from $-1.30$ (20220920A) to $+2.14$ mag (20220319D). Marked `corrected` to the Sharma+2024 value. |
| 20220914A, 20220920A, 20221012A | $m_r$ pick source | | 20.05, 19.46, 19.67 "sharma2024" | These are the Law+2024 Table 2 PS1 values; Sharma+2024 ED Table 1 lists 20.20, 19.43, 19.64 (BASS $r$). The pick values/source labels need correcting. |
| 20190102C | $m_r$ | 20.77 | 20.77 (src bhardwaj2024) | Gordon+2023 Table 1 lists 20.77 with filter **VLT/FORS2 $I$**, not $r$ (Table 5: no $r$ entry; $g = 21.87$, $I = 20.77$, $z = 20.54$). Both Bhardwaj and our pick carry an $I$-band magnitude in the $r$ column. Use our GALFIT $m_r$ or flag. |
| 20211127I | $m_r$ | 14.97 | 14.96 | Gordon+2023 Table 1 says 14.97, Table 5 says $14.96\pm0.01$ (same SOAR measurement). Not a real difference. |
| 20171020A | $m_r$ | 14.97 | 14.97 (src bhardwaj2024) | Source is Lee-Waddell+2023 Table 1 (WISE$\times$SCOSPZ $r$, SuperCOSMOS based); their Table 6 Pan-STARRS1 $r = 15.141$. Re-label the pick source as leewaddell2023; consider the PS1 value. The coincidence with 20211127I (also 14.97) is real, not a transcription slip. |
| 12 CRAFT hosts | $z$ | Gordon Table 3 Prospector $z$ | spectroscopic $z$ | differ by $\le 0.0003$; confirmed. |

Everything else (CRAFT SFRs and masses, 20171020A, 20220912A mass/$z$/$m_r$, all 23 inclination pairs) is confirmed against the cited paper and equals the current pick.

## 3. DM$_{\rm host}$ convention

### What the paper says (exact text)

The only description is Table 1 footnote d (arXiv v2 p. 6; the footnote marker is typeset on the "Host $m_r$" column header, but the text clearly refers to DM$_{\rm host}$):

> "d: These values are estimated after subtracting the fiducial contributions of the Milky Way ISM[49], and intergalactic medium (using the median value of the Macquart relation[50]) from the FRB DMs, which are reported in their respective discovery papers."

with [49] = Cordes & Lazio 2002 (NE2001) and [50] = Macquart et al. 2020. The main text adds only (p. 3):

> "In Figure 2b, we show a box-plot illustrating the range of DM values contributed by FRB hosts (DMhost) in each of the four bins. The DMhost for our sample FRB hosts are listed in Table 1."

Nothing in the Methods addresses DM$_{\rm host}$. Hence the paper does **not** state: which cosmology or $F$ enters the Macquart median, whether a Milky Way halo term is subtracted, or whether the value is observer-frame or rest-frame ($\times(1+z)$).

### What the numbers imply

For each host I recomputed DM$_{\rm host}$ with the discovery-paper DM and NE2001 (Law+2024 Table 1 for DSA; Bannister+2019, Macquart+2020, Heintz+2020, Bhandari+2020b/2022, James+2022, Caleb+2023, Mahony+2018 for the others), the median of the Macquart+2020 $p(\Delta)$ with Planck 2018 ($F = 0.2$ and $F = 0.32$), no halo term, in both frames; YMW16 was evaluated at the host positions with the ATNF web calculator (20220509G from Connor+2023).

- **Rest frame.** Observer-frame values miss by an amount that grows with $z$ (e.g. 20200906A: 309 vs 223 obs / 305 rest; 20210410D: 468 vs 410 / 464; 20220914A: 538 vs 487 / 542). With rest-frame scaling 13/23 hosts agree with our `dm_host_restframe_calc` to within 10 pc cm$^{-3}$ and 18/23 within 18. The column is therefore $({\rm DM} - {\rm DM}_{\rm MW} - \widetilde{\rm DM}_{\rm cosmic})(1+z)$ even though the footnote does not say so.
- **No halo term.** Adding 30-50 pc cm$^{-3}$ would move every value down by 30-70 pc cm$^{-3}$; the data are consistent with none.
- **Median amplitude.** For the 18 hosts not affected by the DM$_{\rm MW}$ issue below, the rms residual against our $F = 0.2$ median is 12 pc cm$^{-3}$ (mean $+5.6$); $F = 0.25$-$0.32$ gives rms 7 (mean $\approx 0$). Their median is thus slightly lower than ours at $z \gtrsim 0.28$ (NE2001 + $F = 0.32$ reproduces 20180924B, 20190102C, 20210320C, 20211203C, 20220310F to $\le 4$), but this is degenerate with the YMW16 explanation for four of those five.
- **DM$_{\rm MW}$ is not uniformly NE2001.** Six hosts are reproduced only (or far better) by YMW16: 20171020A (needs 25.5; NE2001 38, YMW16 24.7), 20180924B (needs 26.8; 40.5 / 27.6), 20210807D (needs 91.3; 121.2 / 93.7), 20211203C (needs 47.0; 63.4 / 48.4), 20220307B (needs 191.8; 135.7 / 187.0), 20220825A (needs 92.3; 79.7 / 87.0). For 20210807D and 20220307B no choice of $F$ or frame helps; only YMW16 works. Conversely 20210410D (needs 56.3; NE2001 56.2, YMW16 42.2), 20220509G, 20220920A and 20190714A are reproduced by NE2001 and not by YMW16. The simplest reading is that the Galactic term was taken from whichever model/value each discovery paper tabulated or quoted first, not from a single NE2001 run.
- **20220319D and 20220310F** are reproduced by neither model: 20220319D (DM 110.95) needs DM$_{\rm MW} \approx 87$-$89$, which is the H$\alpha$-emission-measure DM$_{\rm ISM} \approx 89$ that Ravi+2023b derive for this sightline (NE2001 133, YMW16 188 both exceed the total DM and give negative values); 20220310F needs DM$_{\rm MW} \approx 23$ with $F = 0.2$ (NE2001 45, YMW16 40) or NE2001 with the $F = 0.32$ median (19 vs 22).

### Reconciliation table (pc cm$^{-3}$)

DM, NE2001 = discovery-paper values; YMW16 = ATNF calculator at the host position; $\widetilde{\rm DM}_{\rm cosmic}$ = median of Macquart+2020 at $z_B$ ($F = 0.2$, Planck 2018); "ours" = `dm_host_restframe_calc` from `host_properties_wide.csv` (uses our picked DM and DM$_{\rm MW}$, which differ from the discovery-paper inputs by at most 4 pc cm$^{-3}$, the largest being 20210410D with DM 575 vs 578.78, so "ours" and a recomputation from the discovery inputs agree to $\le 5$); "DM$_{\rm MW}$ needed" = the Galactic DM that reproduces Bhardwaj's value in the rest frame with the $F = 0.2$ median.

| FRB | $z_B$ | DM | NE2001 | YMW16 | $\widetilde{\rm DM}_{\rm cosmic}$ (F=0.2) | Bhardwaj | ours (rest, NE2001) | B $-$ ours | obs-frame NE2001 | rest NE2001 F=0.32 | rest YMW16 F=0.2 | DM$_{\rm MW}$ needed | closest recipe |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 20200906A | 0.3688 | 577.80 | 36.0 | 37.9 | 319.1 | 309 | 304.9 | +4.1 | 222.7 | 325.2 | 302.3 | 33.0 | NE2001, F=0.2, rest |
| 20210410D | 0.1415 | 578.78 | 56.2 | 42.2 | 112.5 | 468 | 463.8 | +4.2 | 410.1 | 472.4 | 484.1 | 56.3 | NE2001, F=0.2, rest |
| 20211127I | 0.0469 | 234.83 | 42.5 | 31.5 | 35.1 | 171 | 164.0 | +7.0 | 157.2 | 165.2 | 176.1 | 36.4 | either (NE2001 +6, YMW16 $-5$) |
| 20211203C | 0.3437 | 636.20 | 63.4 | 48.4 | 295.2 | 395 | 373.3 | +21.7 | 277.6 | 391.1 | 393.1 | 47.0 | YMW16, F=0.2 (or NE2001, F=0.32) |
| 20211212A | 0.0707 | 206.00 | 27.1 | 27.5 | 53.9 | 127 | 134.0 | $-7.0$ | 125.0 | 135.2 | 133.5 | 33.5 | NE2001/YMW16 equal, $-7$ |
| 20220310F | 0.4780 | 462.24 | 45.4 | 39.5 | 424.7 | 22 | $-13.2$ | +35.2 | $-7.9$ | 19.3 | $-2.9$ | 22.7 | NE2001, F=0.32 (DM$_{\rm MW}$ 23 otherwise) |
| 20220914A | 0.1139 | 631.28 | 55.2 | 51.1 | 89.2 | 538 | 542.1 | $-4.1$ | 486.9 | 545.3 | 546.9 | 59.1 | NE2001, F=0.2, rest |
| 20220920A | 0.1582 | 314.99 | 40.3 | 33.4 | 126.8 | 170 | 171.1 | $-1.1$ | 147.9 | 176.4 | 179.4 | 41.4 | NE2001, F=0.2, rest |
| 20221012A | 0.2847 | 441.08 | 54.4 | 50.5 | 240.1 | 190 | 190.5 | $-0.5$ | 146.6 | 201.7 | 193.3 | 53.1 | NE2001, F=0.2, rest |
| 20171020A | 0.0087 | 114.10 | 38.0 | 24.7 | 6.3 | 83 | 70.4 | +12.6 | 69.8 | 71.3 | 83.8 | 25.5 | YMW16 only |
| 20190714A | 0.2365 | 504.13 | 38.5 | 31.2 | 196.0 | 335 | 333.5 | +1.5 | 269.6 | 343.3 | 342.5 | 37.2 | NE2001, F=0.2, rest |
| 20220207C | 0.0430 | 262.38 | 79.3 | 83.3 | 32.1 | 149 | 157.4 | $-8.4$ | 151.0 | 158.0 | 153.3 | 87.4 | YMW16 ($-4$) / NE2001 ($-8$) |
| 20220509G | 0.0894 | 269.53 | 55.2 | 46.0 | 69.0 | 154 | 158.3 | $-4.3$ | 145.3 | 160.3 | 168.3 | 59.2 | NE2001, F=0.2, rest |
| 20220825A | 0.2414 | 651.24 | 79.7 | 87.0 | 200.4 | 445 | 460.6 | $-15.6$ | 371.1 | 470.9 | 451.7 | 92.3 | YMW16 ($-7$) |
| 20220912A | 0.0770 | 219.46 | 125.2 | 122.2 | 59.0 | 34 | 37.9 | $-3.9$ | 35.3 | 39.5 | 41.3 | 128.9 | NE2001, F=0.2, rest |
| 20220307B | 0.2481 | 499.27 | 135.7 | 186.9 | 206.5 | 126 | 193.3 | $-67.3$ | 157.0 | 206.7 | 132.0 | 191.8 | YMW16 only |
| 20180924B | 0.3212 | 361.42 | 40.5 | 27.6 | 274.1 | 80 | 62.3 | +17.7 | 46.9 | 78.2 | 78.9 | 26.8 | YMW16, F=0.2 (or NE2001, F=0.32) |
| 20190102C | 0.2909 | 363.60 | 57.3 | 43.3 | 245.8 | 94 | 79.3 | +14.7 | 60.5 | 91.9 | 96.2 | 45.0 | NE2001, F=0.32 (or YMW16, F=0.2) |
| 20190608B | 0.1178 | 338.70 | 37.2 | 26.6 | 92.5 | 238 | 234.8 | +3.2 | 209.0 | 236.8 | 245.5 | 33.3 | NE2001, F=0.2, rest |
| 20191001A | 0.2342 | 506.92 | 44.2 | 31.1 | 193.9 | 340 | 332.2 | +7.8 | 268.8 | 341.5 | 347.9 | 37.5 | NE2001, F=0.2, rest |
| 20220319D | 0.0112 | 110.98 | 133.3 | 211.0 | 8.1 | 16 | $-30.7$ | +46.7 | $-30.5$ | $-30.1$ | $-109.4$ | 87.0 | none; needs DM$_{\rm MW}\approx 87$-89 (Ravi+2023b H$\alpha$ EM) |
| 20210807D | 0.1293 | 251.90 | 121.2 | 93.7 | 102.1 | 66 | 32.5 | +33.5 | 28.6 | 36.0 | 63.4 | 91.3 | YMW16 only |
| 20210320C | 0.2796 | 384.80 | 42.2 | 30.4 | 235.4 | 148 | 137.4 | +10.6 | 107.2 | 150.2 | 152.3 | 33.8 | NE2001, F=0.32 (or YMW16) |

The five hosts singled out in the task:

- 20180924B (80 vs our 62): YMW16 (27.6 vs NE2001 40.5) with $F = 0.2$ gives 78.9; alternatively NE2001 with $F = 0.32$ gives 78.2. Bhardwaj's $z$ (0.3212) is the same as ours.
- 20210807D (66 vs 33): only YMW16 (93.7 vs NE2001 121.2) works, giving 63.4; $z$ identical. No $F$ or frame choice recovers 66 with NE2001.
- 20220307B (126 vs 193): only YMW16 (187.0 vs NE2001 135.7) works, giving 132; with our NE2001 pick (135.7) the rest-frame value is 193 and the observer-frame value 157. $z_B = 0.2481$ (ours 0.2507 is the questionable Connor+2025 value; it changes the result by $< 1$).
- 20220310F (22 vs $-13$): needs DM$_{\rm MW} \approx 23$ with our median (NE2001 45.4, YMW16 39.5 both too large) or the same NE2001 with a lower median (our $F = 0.32$ run gives 19). $z$ identical (0.4780 vs 0.479).
- 20220319D (16 vs $-31$): needs DM$_{\rm MW} \approx 87$-$89$; this matches the H$\alpha$-EM estimate of DM$_{\rm ISM}\approx 89$ in Ravi+2023b Sec. 3 (their pulsar-based alternatives, 46.7 and 65.3, give 57 and 38). NE2001 (133) and YMW16 (188-211) both exceed the total DM.

Verdict for the CSV: `confirmed` where $|{\rm B} - {\rm ours}| \le 10$ (13 hosts), `conflict` otherwise (10 hosts: 20171020A, 20180924B, 20190102C, 20210320C, 20210807D, 20211203C, 20220307B, 20220310F, 20220319D, 20220825A). `verified_value` is our homogeneous NE2001 / $F = 0.2$ / rest-frame number; the notes column gives the recipe that reproduces Bhardwaj.

## 4. $m_r$ column

The Methods never say where $m_r$ comes from. The only references to magnitudes are the selection statement ("hosts with an apparent r-band magnitude (mr) < 21 AB"), the image sources (Table 1 footnote a: Pan-STARRS DR2 $r$ for bold rows, DESI DR10 $r$ otherwise; Methods: "deepest available archival r-band images of the host (DESI or Pan-STARRS if DESI does not cover the FRB field-of-view)"), and the Extended Data Fig. 4 legend, which labels the SDSS cut as "rPmag < 21 AB" (SDSS Petrosian $r$). Table 1 footnote d is attached to the "Host $m_r$" header but describes DM$_{\rm host}$ (the DM column has no footnote), so no footnote actually documents $m_r$.

Empirically:

- The 12 CRAFT/MeerKAT values are Gordon+2023 Table 1 host magnitudes verbatim (mixed surveys: DECaLS, DES, Pan-STARRS, SOAR $r$; and the FORS2 $I$ band for 20190102C).
- 20171020A is the Lee-Waddell+2023 Table 1 catalogue value (14.97, WISE$\times$SCOSPZ), 20220912A the PS1 value quoted by Ravi+2023 (19.65).
- The 9 DSA values match no published catalogue magnitude (Law+2024 PS1/Legacy Table 2, Sharma+2024 ED Table 1). The two large, nearby galaxies are much fainter in Bhardwaj (20220207C 18.24 vs 16.15; 20220319D 15.39 vs 13.25), which is what a point-source or small-aperture catalogue magnitude (e.g. PS1 `rPSFMag`, or a Tractor-deblended piece) would give; 20220920A is 1.3 mag brighter. These look like the authors' own catalogue look-ups or aperture measurements on the PS1 DR2 / DESI DR10 images, not the literature host magnitudes. They should not be used as literature $m_r$; the CSV replaces them with the Sharma+2024 ED Table 1 values. (I could not query the PS1/Legacy catalogues from this machine to pin down the exact aperture.)

## 5. Other observations

- The PDF Table 1 column headers carry footnote letters a-d but "d" is on $m_r$ and there is no letter on DM$_{\rm host}$; this is a typesetting slip in the paper.
- Bhardwaj's reference [15] is the arXiv version of Law+2024 (2307.03344); the values above were checked against the January 2024 draft in `papers/law2024.pdf` (the ApJ 967, 29 version was not re-checked).
- For every DSA host the Bhardwaj $z$ is the Law+2024 pPXF value rounded to 4 decimals; our Connor+2025 picks differ at the 4th decimal except 20220307B (see section 2).
- Our `dm_mw` pick for 20211212A is 27 (Shannon+2025/James+2022 NE2001); Connor+2025 and the repo master file give 38.8. Bhardwaj's 127 is reproduced to $-7$ with 27 and to $+6$ with 38.8, so it does not discriminate.
- Temporary scripts for this audit were kept outside the repository; nothing other than `audit_group_F.csv` and this file was written.
