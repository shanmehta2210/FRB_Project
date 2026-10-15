# Group 2 host-property sweep: notes

Hosts: 20210320C, 20210410D, 20210807D, 20211127I, 20211203C, 20211212A, 20220501C, 20220725A, 20221106A
(all CRAFT/ASKAP ICS localisations except 20210410D, which is MeerKAT/MeerTRAP).

Output: `group_2.csv` (1133 rows, 20 sources), `group_2.bib`. PDFs in `papers/`.

## Papers read (source_short -> paper)

| source_short | paper | used for |
|---|---|---|
| gordon2023 | Gordon+2023 ApJ 952 122 | Prospector SED (logM*, SFR 0-100 Myr, t_m, Z, A_V), multiband photometry (Table 5), z, galaxy type |
| gordon2025 | Gordon+2025 ApJ 993 119 | GALFIT Sersic n, R_e, offsets (I/R and K band), normalised offsets, GC-model comparisons, K-band M/L masses for 20220725A/20221106A |
| shannon2025 | Shannon+2025 PASA 42 e036 | DM, S/N, fluence, boxcar width, NE2001 DM_MW, PATH (P(O), P(O|x), theta, phi), FORS2/HAWK-I photometry (Table 5, not extinction corrected), host IDs for 2022 hosts |
| caleb2023 | Caleb+2023 MNRAS 524 2064 | 20210410D everything (MeerKAT burst, SOAR/VISTA photometry, PATH, DM budget, H-alpha EM DM_host) |
| glowacki2023 | Glowacki+2023 ApJ 949 25 | 20211127I WALLABY HI (M_HI, W50, asymmetries, 2MASS inclination 19.8 deg), SFR indicators |
| roxburgh2025 | Roxburgh+2025 PASA (MeerKAT HI) | 20211127I and 20211212A HI masses, W50, asymmetries, radio SFR, environment |
| bernalescortes2025 | Bernales-Cortes+2025 A&A | MUSE H-alpha direct DM_host (rest frame) for 20210320C and 20210410D, host halo masses |
| james2022 | James+2022 MNRAS 516 4862 | early DM/z/PATH table (20210320C, 20210807D, 20211127I, 20211203C, 20211212A) |
| baptista2024 | Baptista+2024 ApJ 965 57 | DM, NE2001, z for 20211203C, 20220501C, 20220725A, 20221106A |
| khrykin2024 | Khrykin+2024 ApJ 973 151 | FLIMFLAM DM_IGM / foreground-halo / host-halo decomposition for 20211127I and 20211212A |
| sammons2023 | Sammons+2023 MNRAS 525 5653 | 20210320C two-screen scattering, scintillation bandwidth |
| bhardwaj2024 | Bhardwaj+2024 Nature 634 1065 | Table 1 inclinations (two methods), compiled M*, SFR, sSFR, m_r, DM_host |
| leung2025 | Leung+2025 ApJL | rest-frame DM_host^Macquart with YT20 halo, Chabrier-rescaled masses (20210807D, 20211127I, 20211212A, 20220725A) |
| li2025 | Li+2025 ApJ | BAGPIPES re-fits (M*, SFR, sSFR), PetroFit b/a inclinations with q0=0.2, projected offsets and areas |
| connor2025 | Connor+2025 Nat. Astron. | compiled DM, NE2001, z (Extended Data) |
| hussaini2025 | Hussaini+2025 ApJL 993 L27 | compiled DM, DM_exgal, coords |
| glowacki2025 | Glowacki+2025 PASA 42 e157 | Table 1 burst props (rest-frame tau, |RM|, pol fractions), Table 2 host props (incl. H-alpha EW, linear Z), **Table 3 GALFIT R_e, offset, inclination, b/a (Marnoch et al. in prep)** |
| scott2025 | Scott+2025 PASA 42 e133 | high-time-resolution DM (S/N- and structure-maximising), widths, tau_obs/alpha/tau_1GHz, L/I, V/I, P/I, RM, RM_MW, scintillation |
| muller2026 | Muller+2026 ApJ 1001 118 | FURBY X-Shooter: line fluxes, H-alpha SFR, m_R, BPT class for 20220725A and 20221106A |
| athukoralalage2026 | Athukoralalage+2026 arXiv:2608.23868 | rest-frame DM_host (NE2025 + YT20), |RM|_host for 20210410D, 20210807D, 20211127I, 20211212A, 20220725A |

Not used / not fetched: Bhandari+2023 ApJ 948 67 is the FRB 20210117A dwarf-host paper and does not cover 20211127I; Marnoch+2023 does not cover any Group 2 host; Marnoch et al. (in prep) holds the GALFIT b/a values tabulated in Glowacki+2025 Table 3; Pastor-Marazuela+2026 only mentions 20210410D in passing (MeerTRAP integration time) with no host numbers. The "Glowacki 2024" paper named in the brief appears to be a mis-citation of Glowacki+2023 (WALLABY, 20211127I); the 20211212A HI data are in Roxburgh+2025.

## Conventions and systematic differences to watch

- **DM_host frames.** Bhardwaj+2024 Table 1 "DM_Host" is the observer-frame residual $DM - DM_{\rm NE2001} - \langle DM_{\rm cosmic}\rangle$ with no MW-halo term (e.g. 148 for 20210320C). Bernales-Cortes+2025, Leung+2025 and Athukoralalage+2026 quote rest-frame values (multiplied by $1+z$) and subtract a MW halo term (40, YT20, YT20 respectively). Caleb+2023 quote $DM_{\rm host}/(1+z)=361.4$ for 20210410D. These are **not** interchangeable; the CSV `method` strings encode frame and halo model.
- **Stellar masses.** Gordon+2023 (Prospector, Kroupa IMF) is the common root; Bhardwaj+2024, Gordon+2025, Glowacki+2025, Bernales-Cortes+2025 copy it. Leung+2025 subtract 0.03 dex for Chabrier. Li+2025 BAGPIPES re-fits differ by up to 0.5 dex (20210807D: 10.41 vs 10.97; 20211127I: 10.03 vs 9.48). Gordon+2025 derive K-band M/L masses for 20220725A (11.26) and 20221106A (11.27) that are 0.5 dex above Leung+2025 / Li+2025 (10.71, 10.70 / 10.78).
- **SFR.** Gordon+2023 SFR is integrated over 0-100 Myr. Li+2025 BAGPIPES SFRs are systematically different (often higher). Glowacki+2023/Roxburgh+2025 radio and WISE SFRs for 20211127I are ~2-3 Msun/yr vs Prospector 35.8 (see below).
- **log sSFR.** Gordon+2023 do not tabulate sSFR; the values in Bhardwaj+2024 (and `derived_from_table3` rows here) are log10(SFR_med / M*_med). Gordon+2023 text quotes per-draw medians for two hosts that differ (20210807D: -11.41 vs -11.17; 20211127I: -8.29 vs -7.93).
- **Inclinations.** Four independent published sets: Bhardwaj+2024 (parametric and isophote, r band), Li+2025 (PetroFit b/a, q0=0.2), Glowacki+2025 Table 3 (GALFIT b/a, Marnoch in prep; conversion not stated), Glowacki+2023 (2MASS, 20211127I only). They disagree by up to 20 deg for 20210320C (27-28 vs 48.13) and 20220725A (58.58, no Bhardwaj value).
- **Photometry.** Gordon+2023 Table 5 is Galactic-extinction corrected; Shannon+2025 Table 5 FORS2/HAWK-I integrated magnitudes are **not** corrected; Muller+2026 FORS2 R is corrected. Differences of 0.1-0.5 mag between the two FORS2 sets are mostly extinction.
- **Redshift precision.** Gordon+2025 Table 1 rounds z to 3 dp. James+2022/Leung+2025/Roxburgh+2025 use z=0.0715 for 20211212A vs 0.0707 in Gordon+2023/Shannon+2025/Scott+2025; Khrykin+2024 use 0.0713.
- **PATH.** Shannon+2025 Table 3/4 posteriors assume P(U)=0 and are given to 4 dp; James+2022 values are earlier and slightly lower (20210807D: 0.957 vs 1.0000).

## Per-host notes

### FRB 20210320C

Host at 13:37:50.13 -16:07:21.2 (Gordon+2025), z=0.2797 (X-shooter; 0.2796 Prospector). DM 384.8 (Shannon; 384.59 structure-maximising in Scott+2025), NE2001 42 (Connor+2025 tabulate 39.2). Prospector: logM*=10.37, SFR=3.51, t_m=4.56 Gyr, log Z*=-0.82, A_V,young=1.26; star-forming. SOAR r=19.47. GALFIT (Gordon+2025, I band): n=1.33, R_e=0.89"=3.83 kpc, offset 0.83"=3.57 kpc, 0.93 R_e. Glowacki+2025 Table 3 gives a different R_e (4.87 kpc) and offset (2.22 kpc), i=48.13 deg, b/a=0.684.
Inclination conflict: Bhardwaj+2024 28/27 deg, Li+2025 27.62 deg, Glowacki+2025 48.13 deg.
DM_host: Bernales-Cortes+2025 direct MUSE 113+/-4 (rest, ISM local 79 + halo 34), Macquart 78, max 224; Bhardwaj+2024 148 (observer frame). Sammons+2023: two-screen scattering, tau=0.247 ms at 824 MHz, scintillation bandwidth 0.91 MHz, screen within ~3 Mpc of source. Scott+2025: tau_obs=0.193 ms (828 MHz), alpha=-4.4, RM=288.8 rad/m2 (edge-beam, calibration flagged), L/I=0.86.
Oddity: James+2022 and Bernales-Cortes+2025 give DM 384.8/385 while Connor/Hussaini/Scott use 384.59.

### FRB 20210410D

MeerKAT/MeerTRAP localisation (Caleb+2023). Host J214420.69-791904.8, z=0.1415 (GMOS). DM 578.78+/-2 (S/N-max), 571.16 (structure-max), 572.65 (scattering-corrected); NE2001 56.2, YMW16 42.2; halo 40 (PZ19). Prospector (Caleb Table 1 / Gordon Table 3 differ slightly: logM* 9.46 vs 9.47, log Z* -1.08 vs -1.04): SFR 0.03, t_m 6.78 Gyr, transitioning. SOAR g/r/i/z + VISTA Y/J (r=20.65). PATH P(O|x)=0.9957 (second candidate 0.0043 at 4.69"). Offset 0.78"=2.9 kpc, phi=0.83".
DM_host: Caleb Macquart budget 361.4 observer frame (rest ~412), H-alpha EM 207.5 (source frame); Bernales-Cortes+2025 direct 39+/-3 (ISM local 22.1 + halo 17), Macquart 409, max 481; Athukoralalage+2026 406 (+44/-59) rest; Bhardwaj+2024 468 observer frame. Note the 10x spread between direct H-alpha (39) and budget (~410) methods: Caleb attribute the excess to dense circumburst/ISM gas; strong scattering tau=29.4 ms at 1 GHz supports a host origin.
RM: Caleb+2023 had no full-pol data (only Stokes I); the 3D RM synthesis peak at -77929 rad/m2 is a 6-sigma imaging artefact candidate, not a measurement. Glowacki+2025 Table 1 nonetheless tabulates log(|RM_ex|(1+z)^2)=1.705 for this FRB; origin unclear.
No Gordon+2025/Glowacki+2025 morphology (not a CRAFT host). Inclination: Bhardwaj 41/41, Li 42.48 deg. Caleb SFR lower error printed as "0.061" in the PDF (taken as 0.01).

### FRB 20210807D

Host at 19:56:52.81 -00:45:44.5, z=0.1293 (James+2022 give 0.12927). DM 251.9+/-0.2 (Shannon; Connor/Hussaini/Li/Athukoralalage use 251.3); NE2001 121. Localised with the 10 s hardware-correlator visibilities (no voltage dump), so no high-time-resolution burst properties in Scott+2025 beyond the Glowacki+2025 rest-frame tau (log 1.682). Prospector: logM*=10.97, SFR=0.63, t_m=8.36 Gyr, quiescent (the only quiescent host in Gordon+2023), log sSFR -11.17 (table) / -11.41 (text). Pan-STARRS r=17.17. GALFIT: two-component fits in both I and K (n=0.52/2.92 I; 0.31/2.05 K; R_e 11.2+3.3 kpc I), offset 4.81"=11.21 kpc (I), 4.58"=10.69 kpc (K), ~1.0-1.1 R_e. Glowacki+2025 Table 3: R_e=15.13 kpc, offset 8.92 kpc, i=46.17, b/a=0.707.
DM_host: Leung+2025 -59 (rest, YT20 halo 73; negative allowed), Athukoralalage+2026 20 (+29/-16) flagged "consistent with zero", Bhardwaj+2024 66 (observer frame). Masses: Li+2025 BAGPIPES 10.41 vs Prospector 10.97.
Inclination: Bhardwaj 41/41, Li 44.5, Glowacki+2025 46.17 deg.

### FRB 20211127I

Host WALLABY J131913-185018 at 13:19:13.97 -18:50:16.7, z=0.0469 (6dFGS 0.04695). DM 234.83+/-0.08 (Shannon/Scott; Connor/Hussaini 234.97); NE2001 42-43; Khrykin ISM+halo 82.5. Burst: S/N 37.9, fluence 35 Jy ms, RM=-67 rad/m2 (no pol calibrator), L/I=0.244 (low), V/I=0.129, scintillation bandwidth 2.88 MHz, tau_obs 0.025+/-0.02 ms.
Stellar: Prospector logM*=9.48, SFR=35.83 (0-100 Myr), t_m=3.85 Gyr, log Z*=-0.53, log Z_gas=+0.29, A_V,young 1.22, star-forming; SOAR r=14.97 (Table 1) / 14.96 (Table 5). SFR conflict: Prospector 35.8 vs Glowacki+2023 radio 2.91, WISE W3 2.57, NUV 1.22, Prospector last-30-Myr 0.45; Roxburgh+2025 radio 2.82; Li+2025 1.33 (quality flag 4, bad fit). Gordon+2023 note the Prospector SFR for this host is extreme; Glowacki+2023 Table 2 lists 0.45 for the 30 Myr Prospector SFR. Treat 35.8 with caution.
HI: WALLABY M_HI=6.5e9 (log 9.81), W50=150, A_2D=0.23, A_flux=1.12 (possible low-S/N asymmetry); MeerKAT (Roxburgh) log M_HI=9.99, W50=158, A_2D=0.091, "clearly undisturbed". Glowacki+2023 inclination 19.8 deg from 2MASS axis ratio; Roxburgh+2025 use i<10 deg cut note.
Morphology: GALFIT two-component in I (R_e 8.57"/0.93"), single in K (R_e 3.89"=3.62 kpc, n=2.18); flocculent spiral; offset 2.38"=2.21 kpc (I), 0.28 R_e. Glowacki+2025 Table 3: R_e=10.11 kpc, offset 2.10 kpc, i=18.03, b/a=0.953.
DM_host: Leung+2025 109 rest, Athukoralalage+2026 108 (+23/-27) rest, Bhardwaj+2024 171 observer; Khrykin+2024 host halo only 24.1, DM_IGM 61.6, foreground halos 22.6. |RM|_host=70+/-6 (Athukoralalage).
Inclination: Bhardwaj 21/26, Li 13.99, Glowacki+2023 19.8, Glowacki+2025 18.03 deg.

### FRB 20211203C

Host at 13:38:15.03 -31:22:48.5, z=0.3439 (0.3437 Prospector; Baptista 0.344). DM 636.2+/-0.4 (Shannon; Connor/Li 635.0); NE2001 63.4. Burst: scattered (tau_obs 1.66 ms at 891 MHz, alpha=-9.7+/-2.4 — unphysically steep, single-component fit), RM=34.3, L/I=0.57. Prospector: logM*=9.76, SFR=15.91, t_m=2.47 Gyr, log Z*=0.00, star-forming; SOAR r=19.64. GALFIT (R band): n=0.60, R_e=0.46"=2.28 kpc, offset 0.86"=4.26 kpc, 1.87 R_e (the largest normalised offset in the group). Glowacki+2025 Table 3: R_e=2.41 kpc, offset 3.01 kpc, i=25.72, b/a=0.905.
James+2022 had no redshift yet. Li+2025 flag PS1 Kron fluxes as inconsistent (quality 3). Inclination: Bhardwaj 7/11, Li 13.37, Glowacki+2025 25.72 deg. DM_host only from Bhardwaj+2024 (395 observer frame); no rest-frame value published.

### FRB 20211212A

Host SDSS J102924.22+012139.2 at 10:29:24.23 +01:21:38.9, z=0.0707 (James/Leung/Roxburgh 0.0715; Khrykin 0.0713). DM 206.0+/-5 (Shannon; Scott HTR 200+2/-1; Connor/Hussaini 209.0); NE2001 27 (James/Shannon) vs 39 (Leung/Connor) — different NE2001 evaluations, both tabulated. Burst: tau_obs 1.8 ms at 1491 MHz, RM=21+/-7 (corner beam, calibration flagged), fluence 131 Jy ms. Prospector: logM*=10.28, SFR=0.73, t_m=5.83 Gyr, star-forming, log sSFR -10.42; SOAR r=16.44; SDSS u, 2MASS JHK, WISE. GALFIT two-component in I and K (n=0.21/2.03 I), offset 1.67"=2.27 kpc (I), 0.43 R_e. Glowacki+2025 Table 3: R_e=4.41 kpc, offset 2.16 kpc, i=35.58, b/a=0.822.
HI (Roxburgh+2025, MeerKAT): log M_HI=9.65, W50=87, R_HI=25.9 kpc, undisturbed, two faint HI neighbours ~2 Mpc away; M200 not derived because inclination too low (<10 deg implied by their note). Radio continuum 0.31 mJy, radio SFR 1.92.
DM_host: Leung+2025 81 rest, Athukoralalage+2026 84 (+24/-30) rest, Bhardwaj+2024 127 observer; Khrykin host halo 35.8, DM_IGM 171.2+/-106.
Inclination: Bhardwaj 33/35, Li 33.26, Glowacki+2025 35.58, Roxburgh "<10" (inferred) — the Roxburgh note conflicts with all optical values.

### FRB 20220501C

Host at 23:29:31.02 -32:29:26.5 found only in FORS2 I band (no catalogued galaxy); z=0.381 (X-shooter H-alpha, [OIII]). DM 449.5+/-0.2 (Scott HTR 449.6), NE2001 31; fluence 32 Jy ms; tau_obs 0.35 ms, RM=35.5, L/I=0.68, scintillation 9.2 MHz. PATH P(O|x)=1.0 with prior P(O)=1.0, theta 0.17", phi 0.9". FORS2 g=21.49, I=20.47 (uncorrected). No Prospector fit in Gordon+2023 and no Gordon+2025 morphology (host too faint). Glowacki+2025 Table 3 (Marnoch in prep) gives R_e=2.50 kpc, offset 2.64 kpc, offset/R_e=1.06, i=36.65+/-4.83, b/a=0.81+/-0.05 — the only published size/inclination.
**No published stellar mass, SFR, sSFR, DM_host or r-band magnitude** (only g and I). Li+2025 list DM_exc but no SED fit.

### FRB 20220725A

Host WISEA J233315.68-355925.0 at 23:33:15.69 -35:59:24.9, z=0.1926 (SOAR H-alpha; X-shooter in Muller). DM 290.4+/-0.3 (Scott HTR 290.1); NE2001 31 (Shannon) / 30.7 (Baptista) / 29 (Leung). Burst: fluence 72 Jy ms, tau_obs 2.29 ms at 1149 MHz, alpha=-1.94 (Glowacki+2025 prints +1.94), RM=-26.3 (corner beam, flagged; RM_MW -190.7+/-49.8), circular pol 0.13. Spiral (FORS2; GALFIT n=1.16 in R). FORS2 g=18.53, R=17.84, I=17.23, HAWK-I K=16.26 (uncorrected); Muller R=17.806 (corrected). GALFIT R band R_e=1.47"=4.76 kpc, K band 1.95"=6.31 kpc; offset 0.63"=2.03 kpc (R), 0.43 R_e. Glowacki+2025 Table 3: R_e=5.80 kpc, offset 1.53 kpc, i=58.58, b/a=0.548 (lowest b/a in the group).
Masses: Gordon+2025 K-band M/L 11.26; Leung+2025 NED-LVS 10.71+/-0.3; Li+2025 BAGPIPES 10.70+/-0.10 (SFR 10.73, log sSFR -9.67). Muller+2026: H-alpha SFR 3.97 but host is AGN on BPT (metallicity excluded); line fluxes tabulated. Not in Gordon+2023 and not in Bhardwaj+2024 (no inclination from that paper). DM_host: Leung+2025 76 rest, Athukoralalage+2026 86 (+45/-49) rest; |RM|_host 52+/-7.
Glowacki+2025 Table 2 sparse row "17.81 0.22 0.43 1.20 0.04": last pair taken as [SII] ratio; first pair recorded as `log_SFR` with an UNCERTAIN flag (does not match log of Muller SFR).

### FRB 20221106A

Host WISEA J034649.07-253411.7 at 03:46:49.10 -25:34:10.9, z=0.2044 (SOAR+X-shooter; Muller 0.2043; Scott Table 1 typo 2.044). DM 343.8+/-0.8 (Shannon; Baptista 344.0; Scott HTR 343.2); NE2001 35. Burst: fluence 80 Jy ms, tau_obs 0.182 ms at 1650 MHz, RM=444+/-1 (RM_MW 34.7), L/I=0.862, highest rest-frame scattering in the group (Glowacki log tau=2.774), scintillation 2+/-1 MHz. PATH P(O|x)=0.9708 (second candidate at 3.24", mag 21.07). FORS2 R=18.35, HAWK-I K=16.40 (uncorrected); Muller R=18.322 (corrected), strong Balmer absorption, AGN on BPT, H-alpha SFR 0.535. GALFIT R: n=2.22, R_e=1.69"=5.75 kpc; K: n=2.65, 1.62"=5.51 kpc; offset 1.33"=4.52 kpc (R), 0.79 R_e. Glowacki+2025 Table 3: R_e=6.72 kpc, offset 3.90 kpc, i=50.90, b/a=0.65.
Masses: Gordon+2025 K-band 11.27; Li+2025 10.78 (SFR 2.37, log sSFR -10.39). Not in Gordon+2023, Bhardwaj+2024, Leung+2025 or Athukoralalage+2026. **No published DM_host.** Inclination: Li 39.69, Glowacki+2025 50.90 deg.

## Residual issues

- Gordon+2023 does not list u-r colours in Table 3 for these hosts in the extracted text; `u_r_color` therefore absent.
- Bhardwaj+2024 Table 1 marks Pan-STARRS-imaged hosts in boldface, which is lost in text extraction; the image source per host is flagged as unknown in the CSV notes.
- Scott+2025 Table 1 prints the 20211127I structure-maximising DM upper error as "+14"; recorded as +0.14 with a note.
- Glowacki+2025 Table 3 does not state the band of the GALFIT fits nor the b/a -> inclination conversion.
