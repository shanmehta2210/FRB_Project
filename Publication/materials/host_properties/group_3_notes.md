# Group 3 host-property sweep: notes

Hosts: 20230526A, 20230613A, 20230902A, 20230907D, 20231020B, 20231226A, 20240210A, 20240304A, 20240318A, 20240208A.
Output: `group_3.csv` (397 rows), `group_3.bib`. PDFs and extracted text in `papers/`.

## Papers read (all fetched as PDF from arXiv and text-extracted with `pdfplumber`)

| key | reference | used for |
|---|---|---|
| `shannon2025` | Shannon et al. 2025, PASA 42, e036 (arXiv:2408.02083) | CRAFT/ICS catalogue: DM, DM_MW (NE2001), z, S/N, width, fluence, astrometry, PATH (Table 3/4), FORS2/HAWK-I photometry (Table 5), per-host notes (Sect 5.2) |
| `scott2025` | Scott et al. 2025, PASA 42, e133 (arXiv:2505.17497) | High-time-resolution: DM_S/N, DM_struct, widths, scattering (tau_obs, alpha, tau_1GHz), polarisation fractions, RM, RM_MW, scintillation |
| `gordon2025` | Gordon et al. 2025, ApJ 993, 119 (arXiv:2506.06453) | Host coordinates, z, Sersic re/n (R and K), K-band M/L stellar masses, offsets (Table 4), coflash model preferences, new 20240318A spectroscopy (Table C2) |
| `glowacki2025` | Glowacki et al. 2025, PASA 42, e157 (arXiv:2506.23403) | Rest-frame tau and RM (Table 1), SF flag / m_R / log SFR / [SII] ratio (Table 2), GALFIT Reff, offset, inclination, b/a (Table 3; values credited to Marnoch et al. in prep) |
| `muller2026` | Muller et al. 2026, ApJ 1001, 118 (arXiv:2506.20774; "Muller et al. 2025" in other papers) | FURBY X-shooter: z, m_R (Galactic-extinction corrected), Halpha SFR, metallicity, line fluxes, Balmer-absorption flag |
| `leung2025` | Leung et al. 2025, ApJL 991, L25 (arXiv:2507.16816) | New SED stellar masses and rest-frame DM_host^Macquart for 20230526A, 20231226A, 20240210A |
| `athukoralalage2026` | Athukoralalage et al. 2026, arXiv:2608.23868 (preprint, no journal yet) | Rest-frame DM_host (NE2025 + YT20, Connor 2025 cosmic prior) and rest-frame |RM|_host for the same three hosts |
| `hussaini2025` | Hussaini et al. 2025, ApJL 993, L27 (arXiv:2506.04186) | DM_exgal with NE2001 (mwprop, 30 kpc) for 20231226A, 20240210A, 20240304A |
| `roxburgh2025` | Roxburgh et al. 2025, PASA 42, e164 (arXiv:2507.06594) | MeerKAT HI for 20240210A: M_HI, W50, R_HI, M200, radio SFR, asymmetries, environment |
| `pastormarazuela2026` | Pastor-Marazuela et al. 2026, MNRAS 545, staf2144 (arXiv:2507.05982) | MeerTRAP hosts 20230613A, 20230907D, 20231020B: DM, NE2001, positions, z_spec/z_phot, PATH, DECam+WISE photometry, CIGALE M*/SFR/Z/age, foreground structures |

Papers checked but not usable for these hosts: Baptista et al. 2024 (sample ends 2022; none of our bursts), Gordon et al. 2023 (none), Eftekhari et al. 2025 (20240209A only), Sharma et al. 2024 / Law et al. 2024 / Connor et al. 2023 (DSA, none), Caleb et al. 2023 (20210410D only), the DSA-110 host paper arXiv:2609.30616 (none of ours), Rajwade et al. 2024 (earlier MeerTRAP bursts; the three 2023 MeerTRAP hosts are only in Pastor-Marazuela 2026). "Marnoch et al. in prep" (GALFIT morphology of CRAFT hosts) is not yet on arXiv as of 2026-10-07; its values are only available second-hand through Glowacki 2025 Table 3. No ATels were found carrying host properties for any of these bursts; the ASKAP hosts were first published directly in Shannon 2025 / Gordon 2025 and the MeerTRAP hosts in Pastor-Marazuela 2026.

## Global issues to keep in mind

- **Stellar masses disagree by >1 dex for the three low-z CRAFT hosts.** Gordon 2025 Table 3 uses a K-band $M/L = 0.6$ on the HAWK-I Ks magnitudes from Shannon 2025 Table 5 (10.73, 10.34, 10.33 for 20230526A, 20231226A, 20240210A). Leung 2025 (Table 1, flag b, new SED fits, $\le 0.16$ dex) gives 9.31, 9.59, 9.67, and Athukoralalage 2026 adopts the Leung values. For 20230526A the K-band value is driven by a suspicious Ks = 17.1 (R - Ks ~ 4 mag) in Shannon Table 5. Both are recorded; the Leung SED masses look more trustworthy for 20230526A and 20231226A; for 20240210A (r = 14.9 Seyfert spiral at z = 0.024) neither is obviously wrong.
- **Halpha SFRs:** Glowacki 2025 Table 2 lists log SFR values credited to Muller et al. that are 0.47-0.6 dex below the published Muller 2026 Table 3 SFRs (e.g. 20230526A: $10^{-0.98} = 0.105$ vs 0.359). Same pattern for every FURBY host in that table, so Glowacki used an earlier draft. Published Muller 2026 values should take precedence.
- **No published sSFR for any Group 3 host.** Glowacki Table 2 has no sSFR entries for them; nobody has combined Muller SFRs with Leung/Gordon masses. `log_sSFR` is therefore absent from the CSV.
- **No SED-based A_V, ages, colours, or stellar metallicities** exist for the seven CRAFT hosts; Prospector fits (Gordon 2023 style) have not been published for any of them.
- **DM_MW (NE2001) values in Shannon 2025 Table 1** look unreliable for several hosts: 145 for 20231226A (b ~ +47 deg; Hussaini's NE2001 implies 38.07), 98 for 20240208A, 30 for 20240304A (Hussaini implies 73.76). seed_51.csv already uses independent NE2001 values; Shannon's are recorded as published.
- **Inclination / axis ratio:** the only published b/a and inclinations are in Glowacki 2025 Table 3 (GALFIT, Marnoch et al. in prep), for 20230902A, 20231226A, 20240210A, 20240304A. None for 20230526A, 20240318A, 20240208A or the MeerTRAP hosts. Gordon 2025 fits b/a and PA but does not tabulate them.
- **Offsets:** Gordon 2025 (probability-map, Table 4) and Glowacki/Marnoch (GALFIT centre, Table 3) differ substantially for 20240210A (4.38 vs 2.05 kpc) and 20240304A (8.73 vs 0.86 kpc). Both recorded.

## Per-host notes

### FRB 20230526A
- Read: Shannon 2025, Scott 2025, Gordon 2025, Glowacki 2025, Muller 2026, Leung 2025, Athukoralalage 2026.
- z = 0.1570 (X-shooter, Shannon; Muller re-measures 0.1570). Earlier DES photo-z 0.25 +/- 0.10.
- DM: Shannon 361.4(2). **Scott 2025 Table 1 lists DM_S/N = 316.2 +/- 0.2 and DM_struct 316.1** -- a digit transposition; Athukoralalage uses 361.4.
- Shannon Sect 5.2.6 says the R-band image is saturated at the nucleus so only an upper limit is possible, but Table 5 lists R = 21.08(1); Muller gives m_R = 21.03 +/- 0.01 (extinction corrected). HAWK-I Ks = 17.1(2) gives R - Ks ~ 4, implausible, and inflates the Gordon K-band mass (10.73) vs Leung SED mass (9.31).
- Scattering: tau_obs = 1.16 ms at 1272 MHz, tau_1GHz = 2.75 ms, components 'a' (ambiguous); scintillation bandwidth 2.6 MHz.
- RM = 613 +/- 2 (observed); rest-frame |RM|_host 807 +/- 7 (Athukoralalage), log 2.907 (Glowacki).
- SFR 0.359 +/- 0.003 Msun/yr (Halpha, Muller); 12+log(O/H) = 8.5 +/- 0.1 (O3N2). No sSFR published.
- Morphology: Gordon re = 0.75" (2.05 kpc, n = 0.32) in R; offsets 1.76 kpc / 0.86 re. No GALFIT b/a or inclination in Glowacki Table 3.
- DM_host (rest): 183 (Leung), 189 +42/-57 (Athukoralalage).

### FRB 20230613A (MeerTRAP)
- Read: Pastor-Marazuela 2026 only (Sect 4.1.5, Tables 1, D1, D2, D3).
- z_spec = 0.3923 (Keck/DEIMOS), z_phot 0.42 +/- 0.03, z_Macquart 0.48 +0.14/-0.26.
- DM 483.51 +/- 0.01 (Table 1) but Table D1 prints 483.82; NE2001 = 30.
- Host G1 (DESI DR10) r = 20.132 +/- 0.004, 0.73" from FRB, R_hl 1.39", P(O|x) = 99.9% with P(U) = 7%.
- CIGALE: log M* = 10.13 +/- 0.18, SFR = 4.5 (0.3 dex), log(Z/Zsun) < -0.13, age 2.0 +/- 0.8 Gyr. DECam griz + WISE W1/W2.
- Foreground: Abell 4038 (z ~ 0.03) 5.5' away could add 100-200 pc cm^-3 (M500 printed as 5.78e10 Msun -- almost certainly 1e14).
- Nothing on inclination, morphology, or Sersic fits. Offset 3.4 +/- 3.2 kpc.

### FRB 20230902A
- Read: Shannon, Scott, Gordon, Glowacki, Muller.
- z = 0.3619 everywhere. DM 440.1(1); NE2001 34.
- Photometry: R 21.522(6), Ks 20.36(4) (Shannon); m_R 21.491 +/- 0.006 (Muller, extinction corrected).
- Mass: only Gordon K-band 10.25. SFR 0.72 +/- 0.01 (Muller); Glowacki's log SFR -0.61 is the earlier-draft value.
- GALFIT (Glowacki Table 3): Reff 2.14 kpc, b/a = 0.39 +/- 0.01, i = 69.67 +/- 0.77 deg, offset 2.63 +/- 2.38 kpc (1.23 Reff). Gordon: re 0.46" (2.37 kpc), n = 0.99, offset 4.44 kpc (1.87 re).
- Scattering tau_1GHz 0.072 ms (multiple components); RM 164.8; L/I 0.91.
- Gordon Sect 5.2 notes very weak coflash model preference (Profile marginal).

### FRB 20230907D (MeerTRAP)
- Read: Pastor-Marazuela 2026 (Sect 4.1.8, Tables 1, D1, D3).
- z_spec = 0.4638 (Table 1/D1; text 0.464 +/- 0.015; D3 0.4637), Keck/LRIS 2024 Jan 5. z_phot 0.44 +0.03/-0.02.
- DM 1030.79 +/- 0.04 (D1 prints 1031.11); NE2001 29. Sightline passes 2.7 kpc from the LSB galaxy UGC 7596 (z = 0.0019, Virgo southern filament): Virgo ~230 and UGC 7596 ~50 pc cm^-3 budgeted; z_Macquart after subtraction 0.84 +0.18/-0.44.
- Host r = 19.83 +/- 0.09 (authors' pipeline) vs 19.50 +/- 0.01 (DESI catalogue); two fainter DESI sources G2/G3 at 1.39" and 3.18". PATH with P(U) = 40%: 94.2% (3 candidates) or 98.5% (G1 only); Table D1 prints 0.9585. All three recorded.
- CIGALE (griz + W1-W4): log M* = 10.88 +/- 0.2, SFR 14.91 (0.33 dex), log Z < -0.11, age 2.3 +/- 1.1 Gyr. Offset 0.92" = 5.7 +/- 4.2 kpc; R_hl 1.39".

### FRB 20231020B (MeerTRAP)
- Read: Pastor-Marazuela 2026 (Sect 4.1.10, Tables 1, D1, D3).
- z_spec = 0.4775 (DEIMOS, 730 s; [OIII], Hbeta, [OII]); z_phot 0.46 +0.07/-0.08; z_Macquart 1.04 +0.25/-0.42 -- large DM excess attributed to a foreground photometric cluster (z ~ 0.33, 127" away, up to ~300 pc cm^-3 observer frame).
- DM 952.2 +/- 0.3 (D1 prints 953.84); NE2001 34.
- Host r = 21.81 +/- 0.02 (D1/D3) or 21.79 +/- 0.01 (text); 0.18" from FRB; R_hl 0.88"; P(O|x) 99.8% with P(U) = 35%.
- CIGALE: log M* = 10.75 +/- 0.3, SFR 13.58 with a 1 dex error (poorly constrained), log Z < 0.09, age 2.2 +/- 1.4 Gyr. Offset 1.2 +/- 3.3 kpc.

### FRB 20231226A
- Read: Shannon, Scott, Gordon, Glowacki, Muller, Leung, Athukoralalage, Hussaini.
- Host WISEA J102127.29+060634.5; FRB on a spiral arm of a grand-design spiral (Gordon Sect 4.1 cites it as the grand-design example).
- **Redshift conflict:** 0.1569 (Shannon, Leung, Athukoralalage, Hussaini), 0.157 (Gordon), 0.1570 (Muller) vs **0.1539 in all Glowacki 2025 tables** (unexplained; affects their Reff/offset in kpc by ~2%).
- DM 329.9(1); DM_struct 328.73. Shannon NE2001 = 145 is inconsistent with Hussaini's DM_exgal (implies 38.07).
- Photometry R 18.990(6), Ks 18.08(2); Muller m_R 18.942.
- Mass: Gordon K-band 10.34 vs Leung SED 9.59 (+Athukoralalage). SFR 0.359 +/- 0.005 (Halpha); 12+log(O/H) 8.5 (N2); strong Balmer absorption flagged (post-starburst-like spectrum).
- GALFIT (Glowacki): Reff 5.83 kpc, b/a 0.725 +/- 0.003, i = 44.64 +/- 0.22 deg, offset 5.91 +/- 1.41 kpc (1.01 Reff). Gordon: re 1.83" (5.02 kpc, n 0.79), offset 4.84 kpc (0.97 re).
- Scattering poorly constrained (tau_obs 0.1 +/- 0.07 ms, alpha -1 +/- 3); not used by Glowacki. RM 428.4 +/- 0.3; |RM|_host 556 +/- 7.
- DM_host (rest): 152 (Leung), 152 +42/-55 (Athukoralalage).

### FRB 20240210A
- Read: Shannon, Scott, Gordon, Glowacki, Leung, Athukoralalage, Hussaini, Roxburgh.
- Host WISEA J003506.47-281619.1, r = 14.9 Seyfert 1 spiral, z = 0.023686 from SSRS (archival). Glowacki flags SF = Y despite the Seyfert classification.
- FRB detected near a primary-beam null, so Shannon's fluence (26 Jy ms) is likely underestimated and Scott flags the polarisation calibration as possibly wrong.
- Gordon host centre (00:35:06.52 -28:16:19.3) is ~9" from the FRB; offsets 9.00" = 4.38 kpc = 1.18 re (two-component disk+bulge fit, re 7.64"/8.17"). Glowacki/Marnoch GALFIT gives Reff 3.90 kpc, offset 2.05 +/- 0.41 kpc, b/a 0.6051 +/- 0.0005, i = 54.35 +/- 0.04 deg. Gordon Table 5: optical Residual model preferred over GC Spherical (dBIC 3.44), consistent with the spiral-arm location.
- Masses: Gordon K-band 10.33 vs Leung SED 9.67.
- No optical SFR published; Roxburgh gives radio-continuum SFR 1.19 +0.66/-0.43 Msun/yr, log M_HI = 9.61 +/- 0.05, W50 = 166 km/s, R_HI 16.8 kpc, log M200 = 11.38; HI asymmetries all low, host "settled"; lies in an HI-rich group (neighbour 1.5' SW at z = 0.02360, >= 4 galaxies within 10 Mpc).
- DM_host (rest): 206 (Leung), 188 +17/-22 (Athukoralalage); DM_exgal 255.10 (Hussaini).
- Scintillation bandwidth 2.7 MHz at 863.5 MHz.

### FRB 20240304A
- Read: Shannon, Scott, Gordon, Glowacki, Hussaini.
- Shannon Table 1 lists z as pending but Sect 5.2.14 gives X-shooter z = 0.2423 (Halpha, Hbeta, [OII], [OIII], [NII], [SII]); Gordon 0.242; Hussaini 0.2423.
- DM 652.6(5); DM_struct 653.4 +5.2/-3.5. Shannon NE2001 = 30 vs Hussaini-implied 73.76.
- Scattering is odd: alpha = +3.5 +/- 1.3 (positive) with tau_obs 2.51 ms at 877 MHz, tau_1GHz 4 +/- 0.5 ms. Scott Table 3 caption says no RM could be fit, yet the table lists RM = 489.7 +/- 0.8; Glowacki carries no RM for it. Recorded with a caution note.
- Photometry R 20.723(7), Ks 20.01(5). Mass only Gordon K-band 9.99. No SFR, metallicity, or Halpha fluxes published (not in the FURBY release).
- GALFIT (Glowacki): Reff 4.13 kpc, b/a 0.638 +/- 0.005, i = 51.77 +/- 0.42 deg, offset 0.86 +/- 1.57 kpc (0.21 Reff) -- strongly disagrees with Gordon's 8.73 +2.91/-3.02 kpc (2.37 re). Gordon also notes the coflash model preference differs between optical and NIR.

### FRB 20240318A
- Read: Shannon, Scott, Gordon, Glowacki.
- Shannon: host WISEA J100134.45+373659.7 (<~1.2" away), SDSS g = 19.5, photo-z 0.12 +/- 0.04; not included in Shannon's PATH table or photometry table. Gordon 2025 obtained a new Keck/LRIS spectrum (2025 Jan 2): z = 0.1120 with [OII], Hbeta, [OIII], [NII], Halpha, [SII] emission and CaII H&K absorption (Table C2).
- Imaging is Keck/LRIS R + MMT/MMIRS K (not VLT). Gordon: re 2.40" (4.95 kpc, n 0.85), K-band mass 10.31, offset 2.63 kpc (0.53 re).
- Scattering tau_1GHz 0.128 ms (multiple components); RM -48.5 +/- 0.3; scintillation bandwidth 4.1 MHz; PA microstructure shown in Scott Fig 8.
- No published m_r (only SDSS g = 19.5), no SFR, no inclination/b/a (absent from Glowacki Table 3), no DM_host.

### FRB 20240208A
- Read: Shannon, Scott, Glowacki. **Not in Gordon 2025** (seed_51.csv cites Gordon 2025 as host citation for this burst; that is wrong -- Shannon 2025 Sect 5.2.12 is the only host-identification text).
- Shannon: coincident (<~4") with SDSS J103654.96-005712.2, g = 22.8, photo-z 0.4 +/- 0.1; no optical follow-up at that time. Localisation is the poorest of the group (sigma 0.81 x 1.36"), 14 antennas.
- Shannon Table 1 prints the same UTC (2024-02-08 20:00:54.246) for 20240201A and 20240208A; the 20240201A row is the typo.
- **Glowacki 2025 lists z = 0.2385** (spectroscopic, provenance not stated; presumably CRAFT follow-up to appear in Marnoch et al. in prep) and SF = Y. This is inconsistent with the SDSS photo-z but is the only redshift in print.
- DM 260.2(3); Shannon NE2001 = 98 (suspiciously high). Scattering tau_1GHz 1.0 +/- 0.45 ms (single component); RM -73.7 +/- 1.4; L/I 0.94.
- No host photometry beyond SDSS g, no mass, SFR, size, offset, PATH, or morphology published.
