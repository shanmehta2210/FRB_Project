# Audit notes: Group A (CRAFT/ASKAP earlier sample, 14 hosts)

Hosts: 20171020A, 20180924B, 20181112A, 20190102C, 20190608B, 20190714A, 20191001A, 20200906A, 20210320C, 20210410D, 20210807D, 20211127I, 20211203C, 20211212A.

Companion files: `audit_group_A.csv` (296 rows: 97 confirmed, 2 corrected, 35 conflict, 155 new, 7 unavailable) and `audit_group_A.bib` (`yamanaka2024`, `shannon2018`; everything else is already in `host_properties.bib`).

## Papers opened (text extracts re-read for this audit)

Primary: Shannon et al. 2025 (PASA 42 e036; Tables 1, 2, 3, 5), Gordon et al. 2023 (ApJ 954 80; Tables 1, 3, 5, Sec. 5.1), Gordon et al. 2025 (arXiv:2506.06453; Tables 1, 3, 4), Glowacki et al. 2025 (arXiv:2506.23403; Tables 1, 2, 3), Li et al. 2025 (ApJ 989 77; Table 2), Bhardwaj et al. 2024 (Nature; Table 1), Lee-Waddell et al. 2023 (Tables 1, 6, 7, Sec. 3.3, 3.6), Mahony et al. 2018 (Tables 1, 3), Bannister et al. 2019 (Table 1, S1.14), Prochaska et al. 2019 (Table 1, S3), Macquart et al. 2020 (Methods), Bhandari et al. 2020a (Table 1, 3), Bhandari et al. 2020b (Table 1), Heintz et al. 2020 (Tables 1, 3), Day et al. 2020 (Tables 3, 4), Cho et al. 2020, Mannings et al. 2021 (Table 3), Chittidi et al. 2021 (Table 3, Sec. 3.3-3.4), Bhandari et al. 2022 (Tables 1, A3, A4), James et al. 2022 (Table 2), Caleb et al. 2023 (Table 1, Sec. 4), Glowacki et al. 2023 (Table 2, Sec. 3.1), Sammons et al. 2023, Yamanaka et al. 2024 (Table 3), Khrykin et al. 2024, Scott et al. 2025 (Tables 2, 3), Connor et al. 2025 (compilation), Leung et al. 2025 (Table 1), Roxburgh et al. 2025 (Tables 1, 3), Bernales-Cortes et al. 2025 (Table 4), Athukoralalage et al. 2026 (compilation).

External checks: SDSS DR18 SkyServer (20190608B spectroscopic z; 20211212A no spectrum within 12 arcsec), Crossref for DOIs.

Conventions: `method` strings follow SWEEP_PROTOCOL; Prospector rows = Gordon 2023 setup (non-parametric 8-bin continuity SFH, Kroupa 2001 IMF, Kriek & Conroy dust, dynesty; SFR = mean over last 100 Myr; AGN template on for 20190608B, 20190714A, 20200906A, 20210807D, 20211127I). Glowacki 2025 inclinations reproduce exactly from their $b/a$ with $\cos^2 i = (q^2-q_0^2)/(1-q_0^2)$, $q_0 = 0.2$. Li 2025 uses the same formula (their eq. 2) on PetroFit single-Sersic fits to DESI LS imaging. Bhardwaj 2024 Table 1 lists two inclinations (parametric fit, isophote fit); the parametric one is treated as the pick.

## Global findings

1. **DM / DM_MW**: all 13 ASKAP hosts agree with Shannon 2025 Table 1 to the quoted precision. Differences from discovery papers are at the $\lesssim 1$ pc cm$^{-3}$ level and reflect S/N- vs structure-maximising definitions (recorded in notes). 20210410D (MeerTRAP) is not in Shannon 2025: the seed DM 575 is a rounded compilation value; **corrected to 578.78 $\pm$ 2.0** (Caleb 2023 Table 1). Connor 2025 DM_MW values differ from the primary papers for 20171020A (36.5 vs 38), 20210320C (39.2 vs 42) and 20211212A (38.8 vs 27); primary values kept.
2. **Redshifts**: all confirmed. Note Bhardwaj 2024 used the Prospector-fit redshifts of Gordon 2023 (0.3212, 0.2909, 0.2342, 0.2796, 0.3437) rather than the spectroscopic ones; the table should use the spectroscopic values (Shannon 2025). **20211212A is a genuine conflict**: 0.0707 (Shannon 2025, Gordon 2023 SOAR spectrum) vs 0.0715 (James 2022 "Deller in prep", propagated to Leung 2025 and Roxburgh 2025) and 0.0713 (Khrykin 2024). Recommend 0.0707.
3. **Masses/SFRs**: every Gordon 2023 Table 3 value re-read and confirmed (12 hosts). 20171020A is ProSpect (Chabrier, skewed-normal SFH, "current" SFR, not 0-100 Myr) from Lee-Waddell 2023 Sec. 3.6. Bhardwaj 2024 Table 1 lists $M_* = 1.33 \times 10^{10}$ for 20200906A (the Bhandari 2022 CIGALE value) while citing Gordon 2023 ($2.34 \times 10^{10}$); the seed (10.37) is correct. 20211127I Prospector SFR 35.83 is a flagged outlier (radio/WISE/NUV give 1-3 M$_\odot$/yr).
4. **m_r**: Gordon 2023 Table 1 confirmed for 11 hosts. **20190102C "m_r = 20.77" is VLT/FORS2 I-band**, not r (Gordon 2023 Table 1 Filter column); the only r-band literature value is LCOGT 22.06 $\pm$ 0.51 (Bhandari 2020a). 20171020A seed 14.97 is a SuperCOSMOS photographic R; PS1 r = 15.141 (Lee-Waddell 2023 Table 6) recommended. 20211203C SOAR r = 19.64 is 0.4-0.6 mag brighter than FORS2 R (20.26), PATH photometry (20.29) and LS DR10 (20.06) for the same galaxy (zero-point/aperture issue, flagged).
5. **Inclinations**: for 10 of 13 hosts with multiple photometric determinations, Glowacki 2025 (GALFIT, FORS2/HST) is 5-20 deg higher than Bhardwaj 2024 and Li 2025 (which agree with each other to $\lesssim 3$ deg except 20211127I/20211203C). Only two kinematic values exist: 20190608B 37 $\pm$ 3 (Chittidi 2021, KCWI) agrees with Bhardwaj/Li and not Glowacki (59); 20191001A 61 $\pm$ 5 (Yamanaka 2024, ALMA CO) agrees with Glowacki (62.6) and not Bhardwaj/Li (53). So neither photometric set is systematically "right".
6. **Offsets**: 20191001A seed 11.1 kpc (Bhandari 2022) is **superseded**: Shannon 2025 moved the FRB position by ~2.2 arcsec with CELEBI; Gordon 2025 gives 3.22 (+1.56/-1.45) kpc and Glowacki 2025 3.93 $\pm$ 1.57 kpc. Mannings 2021, Heintz 2020, Li 2025 and Yamanaka 2024 all use the old position.
7. 2025-2026 updates checked: Gordon 2025 (offsets, Re: all 11 ASKAP hosts except 20181112A/20171020A/20210410D), Glowacki 2025 ($b/a$, $i$, offsets, RM compilation), Leung 2025 (Chabrier-rescaled masses, no new measurements), Bernales-Cortes 2025 (DM_host for 20210320C), Athukoralalage 2026 (compilation only; 20210807D z 0.12927, DM 251.3), Roxburgh 2025 (HI masses, radio SFRs for 20211127I, 20211212A).

## Per-host notes

### FRB 20171020A (ESO 601-G036; fly's-eye, no offset)
- z 0.00867 (Lee-Waddell 2023 Table 1; Mahony 2018 HI 0.008672). DM 114.1 (Shannon 2018 via Lee-Waddell/Mahony; not in Shannon 2025). DM_MW NE2001 38, YMW16 26 (Mahony 2018 Table 1); Connor 2025 lists 36.5.
- $\log M_*$ 8.64 (+0.03/-0.15), SFR 0.09 $\pm$ 0.01: ProSpect, BC03, Chabrier, skewed-normal SFH, evolving metallicity (Lee-Waddell 2023 Sec. 3.6; GALEX+PS1+VISTA+WISE). Not a 0-100 Myr SFR. Halpha SFR 0.13 (Mahony 2018 Table 3). Mahony CIGALE-free mass 9e8.
- m_r: seed 14.97 = WISExSCOSPZ photographic R (Lee-Waddell Table 1). PS1 Petrosian r = 15.141 $\pm$ 0.002 (Table 6, SF11-corrected) is the proper CCD value.
- Inclination: Bhardwaj 69/65; Li 2025 68.75. No Glowacki/Gordon 2025 entries. No published Re or offset (fly's-eye 0.38 deg$^2$ region).
- Extras: 12+log(O/H) 8.3 $\pm$ 0.2, Sc, star-forming low-metallicity (Mahony 2018).

### FRB 20180924B (DES J214425.25-405400.81)
- z 0.3214 (Shannon 2025; Bannister 2019 0.3214 $\pm$ 0.0002). Gordon 2023 / Heintz 2020 / Bhardwaj 2024 print 0.3212. DM 362.4 $\pm$ 0.2 (Shannon 2025; Bannister 361.42, Day 2020 HTR 362.16). DM_MW 41 (40.5 Bannister).
- Prospector: $\log M_*$ 10.39 $\pm$ 0.02, SFR 0.62 (+0.32/-0.24), AGN N, $t_m$ 5.63 Gyr, log Z -0.14, A_V 1.10. Halpha 0.88 $\pm$ 0.26 (Heintz 2020). Classified star-forming by Gordon 2023 although Bannister 2019 call it quiescent/LINER-like.
- m_r 20.33 $\pm$ 0.01 (DECaLS, Gordon 2023). Offset 3.4 $\pm$ 0.8 (Bhandari 2022; Heintz 3.4 $\pm$ 0.5; Mannings 3.37 $\pm$ 0.56; Gordon 2025 3.47 (+0.81/-0.83); Glowacki 3.86 $\pm$ 0.94). Re: Gordon 2025 2.22 kpc (HST); Glowacki 7.74 is discrepant.
- Inclination conflict: Bhardwaj 34/36, Li 34.07 vs Glowacki 49.34 $\pm$ 0.66 ($b/a$ 0.669) and Bannister 2019 GALFIT $b/a$ = 0.62 $\pm$ 0.05 ($\Rightarrow$ ~53 deg). HST shows an early-type barred spiral.
- RM 22 $\pm$ 2 (Day 2020; Bannister 14 $\pm$ 1), tau 0.68 ms.

### FRB 20181112A (DES J214923.66-525815.28)
- z 0.4755 (Prochaska 2019 0.4755 $\pm$ 0.00015). DM 589.0 $\pm$ 0.3 (Shannon 2025; Prochaska 589.27; Cho 2020 589.265). DM_MW 40.
- Prospector: 9.87 $\pm$ 0.07, SFR 1.54 (+0.99/-0.65), AGN N. Halpha 0.37 (Heintz). m_r 21.68 $\pm$ 0.05 (DES). Offset 1.7 $\pm$ 19.2 (Heintz 2020; poor 3.8x2.4 arcsec localisation; Glowacki 2.37 $\pm$ 13.64). Not in Gordon 2025.
- Inclination: only Glowacki 2025 27.41 $\pm$ 3.10 ($b/a$ 0.89); not in Bhardwaj 2024, no Li 2025 value. Re 5.51 $\pm$ 0.23 (Glowacki; Heintz 7.2 $\pm$ 1.7). RM 10.9 $\pm$ 0.9.

### FRB 20190102C
- z 0.2912 (Bhandari 2020a 0.2913 $\pm$ 0.0002; Bhardwaj used Prospector 0.2909). DM 364.5 $\pm$ 0.3 (Day 2020 HTR 364.545). DM_MW 57.
- Prospector: 9.69 (+0.09/-0.11), SFR 0.40 (+0.31/-0.11), AGN N. Halpha 0.86 (Heintz).
- **m_r corrected (filter)**: 20.77 $\pm$ 0.05 is VLT/FORS2 I (Gordon 2023 Table 1, Filter column). Only r-band value: LCOGT 22.06 $\pm$ 0.51 (Bhandari 2020a Table 3). Shannon 2025 Table 5: FORS2 g 22.59, I 21.10 (uncorrected). Table must either relabel or substitute.
- Offset 2.3 $\pm$ 4.2 (Bhandari 2022; Heintz 2.0 $\pm$ 2.2; Mannings 2.26; Gordon 2025 5.06 (+3.08/-2.51) median from MC over the 0.8x0.9 arcsec ellipse; Glowacki 1.35 $\pm$ 2.67). All consistent.
- Inclination conflict: Bhardwaj 50/50, Li 45.66 vs Glowacki 62.78 $\pm$ 1.04 ($b/a$ 0.49). Re 4.27 (Gordon 2025 HST), 5.00 (Mannings), 4.71 (Glowacki).
- RM -105 $\pm$ 1 (pulse 2; Day 2020), tau 0.041 ms.

### FRB 20190608B (SDSS J221604.90-075355.9)
- z 0.1178 (SDSS DR18 plate 719/52203/452: 0.1177805). DM 339.5 $\pm$ 0.5 (Day 2020 HTR 340.05). DM_MW 37.
- Prospector: 10.56 $\pm$ 0.02, SFR 7.03 (+1.43/-1.15), AGN Y (LINER/broad Halpha nucleus). Halpha 0.69 (Heintz), 1.2 (Bhandari 2020a). CIGALE 10.06 (Heintz; the 10.4 in Bhandari 2020a was a database error per Heintz App. A).
- m_r 17.41 (DECaLS). Offset 6.5 $\pm$ 0.8 (Bhandari 2022; Heintz 6.6; Mannings 6.52; Gordon 2025 6.14 $\pm$ 0.79; Glowacki 6.38).
- **Inclination seed corrected 36 -> 37**: seed took the isophote (second) value; parametric is 37. Kinematic: Chittidi 2021 i_gas = 37 $\pm$ 3 (KCWI Hbeta, qubefit; 3DBarolo 37), i_stellar = 26 $\pm$ 3 (isophote $b/a$ = 0.90). Li 36.81. Glowacki 59.21 $\pm$ 0.17 ($b/a$ 0.54; single Sersic on a barred spiral) is discrepant. Recommend 37 $\pm$ 3 (kinematic).
- Re: Gordon 2025 3.87/1.23 kpc (two components), Mannings disk 7.37, Heintz 2.8, Glowacki 4.33. RM 353 $\pm$ 2; tau 3.3 ms.

### FRB 20190714A
- z 0.2365 (Keck/LRIS, Heintz 2020). DM 504.7 $\pm$ 0.3 (Li 2025 lists 504.13). DM_MW 39.
- Prospector: 10.22 $\pm$ 0.04, SFR 1.89 (+1.22/-0.72), AGN Y. Halpha 0.65 (Heintz). m_r 20.34 $\pm$ 0.03 (PS1; Bhandari 2022 Kron 20.32).
- Offset 2.7 $\pm$ 1.8 (Bhandari 2022; Heintz 1.9 $\pm$ 1.1; Mannings 1.97; Gordon 2025 3.46 (+1.79/-1.64); Glowacki 2.91).
- Inclination conflict: Bhardwaj 51/50, Li 51.10 vs Glowacki 72.00 $\pm$ 0.12 ($b/a$ 0.363). Re 3.54 (Gordon 2025), 3.85 (Mannings), 4.02 (Glowacki). No RM.

### FRB 20191001A
- z 0.234 (Bhandari 2020b 0.2340 $\pm$ 0.0001; Bhardwaj used 0.2342). DM 506.92 $\pm$ 0.04 (detection; structure-maximised 507.90, used by Li 2025). DM_MW NE2001 44 / YMW16 31.
- Prospector: 10.73 (+0.07/-0.08), SFR 18.28 (+17.24/-8.95), AGN N. Halpha 8.06 (Heintz). m_r 18.36 (DECaLS).
- **Offset conflict / superseded**: seed 11.1 $\pm$ 0.8 (Bhandari 2022, 2020 ASKAP position 21:33:24.37 -54:44:51.86). Shannon 2025 Table 2 CELEBI position 21:33:24.41 -54:44:53.9 ($\pm$0.46, 0.43 arcsec), ~2.2 arcsec south. Gordon 2025 Table 4: 3.22 (+1.56/-1.45) kpc (0.85 arcsec; 0.80 Re); Glowacki 3.93 $\pm$ 1.57. Recommend 3.22. Yamanaka 2024's "FRB 11 kpc north of centre" uses the old position.
- Inclination conflict: Bhardwaj 53/56, Li 53.47 vs Glowacki 62.60 $\pm$ 0.12 ($b/a$ 0.493) and **kinematic 61 $\pm$ 5** (Yamanaka 2024 Table 3, ALMA CO(2-1), 3DBarolo; PA 134 $\pm$ 6, V_rot 205 $\pm$ 26). Recommend the kinematic value.
- Re 4.03 (Gordon 2025), 6.23 (Mannings), 6.25 (Glowacki), CO half-light 6.9 $\pm$ 0.7. Gas mass log 10.36 (Yamanaka). RM 55.5 $\pm$ 0.9; tau 3.3 ms at 824 MHz.

### FRB 20200906A (DES J033358.99-140459.2)
- z 0.3688 (Keck/DEIMOS, Bhandari 2022). DM 577.8 $\pm$ 0.2. DM_MW NE2001 36 / YMW16 38.
- Prospector: 10.37 $\pm$ 0.05, SFR 4.93 (+3.46/-2.34), AGN Y. Bhardwaj 2024 Table 1 prints $1.33 \times 10^{10}$ (Bhandari 2022 CIGALE 10.12) under a Gordon citation; Gordon is $2.34 \times 10^{10}$. Seed 10.37 correct. Halpha 0.48 (Bhandari 2022). m_r 19.95 (DES).
- Offset 5.9 $\pm$ 2.0 (Bhandari 2022 position); Gordon 2025 5.01 (+2.41/-2.33) with the CELEBI position 03:33:58.93 -14:04:58.8 (moved ~2.3 arcsec but offset unchanged within errors); Glowacki 6.45 $\pm$ 2.45.
- Inclination conflict: Bhardwaj 63/59, Li 63.00 vs Glowacki 78.99 $\pm$ 0.13 ($b/a$ 0.274). Re 7.65 (Gordon 2025), 7.58 (Bhandari 2022), 9.01 (Glowacki). No RM.

### FRB 20210320C
- z 0.2797 (SOAR; Bhardwaj used 0.2796). DM 384.8 $\pm$ 0.3 (Scott 2025 HTR 384.6). DM_MW 42 (Connor 2025 39.2).
- Prospector: 10.37 (+0.05/-0.06), SFR 3.51 (+2.44/-1.45), AGN N. m_r 19.47 $\pm$ 0.02 (SOAR). Offset 3.57 (+1.78/-1.76) (Gordon 2025; Glowacki 2.22 $\pm$ 1.61; Li 3.22).
- Inclination: seed Li 27.62 (Table 2, not Table 1), Bhardwaj 28/27 vs Glowacki 48.13 $\pm$ 0.23 ($b/a$ 0.684). Re 3.83 (Gordon 2025), 4.87 (Glowacki).
- RM 288.8 $\pm$ 0.2 (Scott 2025; flagged off-axis), tau 0.247 ms at 824 MHz (Sammons 2023), DM_host 113 $\pm$ 4 rest-frame (Bernales-Cortes 2025 Table 4).

### FRB 20210410D (J214420.69-791904.8; MeerTRAP)
- Primary is Caleb 2023, not Connor 2025. z 0.1415 (GMOS). **DM corrected 575 -> 578.78 $\pm$ 2.0** (S/N-maximising; structure 571.16 $\pm$ 0.97; scattering-corrected 572.65). DM_MW NE2001 56.2 / YMW16 42.2.
- Prospector: 9.47 $\pm$ 0.05 (Caleb prints the same fit as 9.46 (+0.05/-0.06)), SFR 0.03 (+0.03/-0.01), AGN N; "transitioning" class. m_r 20.65 $\pm$ 0.03 (SOAR; Caleb Table 1 same). Offset 2.9 kpc (Caleb 2023, no error; Li copies). No Re (not in Gordon 2025; Glowacki row blank).
- Inclination: Li 42.48, Bhardwaj 41/41; no Glowacki value. Consistent. tau 29.4 ms at 1 GHz; no RM (Stokes I only).

### FRB 20210807D
- z 0.1293 (Keck/LRIS; Athukoralalage 0.12927). DM 251.9 $\pm$ 0.2 (compilations 251.3). DM_MW 121.
- Prospector: 10.97 $\pm$ 0.02, SFR 0.63 (+0.18/-0.17), AGN Y; only quiescent host in Gordon 2023 (borderline). Li Bagpipes 10.41 / 11.33 is inconsistent (0.56 dex lower mass, 18x higher SFR). m_r 17.17 (PS1).
- Offset 11.21 (+1.20/-1.18) (Gordon 2025 I; K 10.69; Glowacki 8.92 $\pm$ 0.96; Li 11.89). Re 11.18/3.31 kpc (two comp.; Glowacki single 15.13).
- Inclination: Li 44.50, Glowacki 46.17 $\pm$ 0.32 ($b/a$ 0.707), Bhardwaj 41/41. Consistent within 5 deg. No RM.

### FRB 20211127I (WALLABY J131913-185018)
- z 0.0469 (6dFGS 0.04695). DM 234.83 $\pm$ 0.08 (compilations 234.97). DM_MW 43 (42.5 James).
- Prospector: 9.48 (+0.06/-0.02), AGN Y. **SFR 35.83 (+1.02/-1.46) flagged**: log sSFR -7.93 is the sample maximum; radio 1.4 GHz 2.91 (+1.70/-1.07) (Glowacki 2023), MeerKAT 2.82 (Roxburgh 2025), WISE W3 2.57, NUV 1.22, Prospector last-30-Myr 0.45. Keep for the uniform column but flag; Li Bagpipes 10.03/1.33 is a bad fit (quality 4).
- m_r 14.97 $\pm$ 0.01 (SOAR; Table 5 prints 14.96). Offset 2.21 (+0.46/-0.47) (Gordon 2025; Glowacki 2.10; Li 1.59). Re 7.97/0.86 (Gordon 2025 I; Glowacki 10.11).
- Inclination: Li 13.99, Glowacki 18.03 $\pm$ 0.31 ($b/a$ 0.953), 2MASS 19.8 (Glowacki 2023), Bhardwaj 21/26. All near face-on; no conflict beyond photometric scatter. Roxburgh 2025 applied a V_rot correction (i >= 10 deg).
- RM -67 $\pm$ 1, tau 0.025 ms (Scott 2025); log M_HI 9.99 (Roxburgh).

### FRB 20211203C (host at 13:38:15.0 -31:22:48.3)
- Host identity verified: Gordon 2023 Table 1 FRB position 13:38:15.00 -31:22:48.2 is the host PATH candidate 204.5626 -31.3801 (Shannon 2025 Table 3, P(O|x) = 1.000, 0.57 arcsec) and the Gordon 2025 Table 1 host centroid 13:38:15.03 -31:22:48.5. All published values recorded refer to this galaxy (HOST_COORDINATE_CHECK.md refit target).
- z 0.3439 (SOAR; Bhardwaj/Li use Prospector 0.3437). DM 636.2 $\pm$ 0.4 (Scott structure 635.16; compilations 635). DM_MW 63.
- Prospector: 9.76 (+0.07/-0.09), SFR 15.91 (+2.82/-2.98), AGN N, rising SFH. m_r 19.64 $\pm$ 0.03 (SOAR, Gordon 2023) confirmed as printed, but **flagged**: FORS2 R 20.258 $\pm$ 0.004 and g 20.842 (Shannon 2025 Table 5, uncorrected, A_r ~ 0.1), PATH photometry 20.29, LS DR10 r 20.06 all put the galaxy ~0.4-0.6 mag fainter; SOAR g 20.32 vs FORS2 g 20.84 shows the same offset, so it is a zero-point/aperture issue for the same source.
- Offset 4.26 (+2.34/-1.78) (Gordon 2025 R; 0.86 arcsec; 1.87 Re; Glowacki 3.01; Li 5.04). Re 2.28 kpc = 0.46 arcsec (barely resolved; Glowacki 2.41).
- Inclination conflict: Li 13.37, Bhardwaj 7/11 vs Glowacki 25.72 $\pm$ 0.43 ($b/a$ 0.905); repo refit on the correct host gives $b/a$ 0.92, i 23.6, agreeing with Glowacki. Recommend Glowacki with a "barely resolved" caveat. Bhardwaj 7 cannot be traced to coordinates.
- RM 34.3 $\pm$ 1.2, tau 1.66 ms at 891 MHz (Scott 2025).

### FRB 20211212A (SDSS J102924.22+012139.2)
- **z conflict**: 0.0707 (Shannon 2025 Tables 1-2; Gordon 2023 SOAR 2021-12-08) vs 0.0715 (James 2022 Table 2 "Deller in prep"; Leung 2025; Roxburgh 2025) and 0.0713 (Khrykin 2024). No SDSS spectrum (DR18 checked). Recommend 0.0707.
- DM 206 $\pm$ 5 (1631.5 MHz; Scott 2025 HTR 200 +2/-1; compilations 209). DM_MW 27 (James 27.1; Connor 38.8 and Leung 39 differ).
- Prospector: 10.28 (+0.05/-0.06), SFR 0.73 (+0.62/-0.39), AGN N. Radio SFR 1.92 (Roxburgh). m_r 16.44 (SOAR; SDSS petro 16.62 corrected). Offset 2.27 (+0.62/-0.59) (Gordon 2025; Glowacki 2.16; Li 2.68). Re 5.26/2.58 (Glowacki 4.41).
- Inclination: Li 33.26, Glowacki 35.58 $\pm$ 0.16 ($b/a$ 0.822), Bhardwaj 33/35: consistent. Roxburgh 2025 Table 3 note b treats this host as too face-on (< 10 deg) for a V_rot correction, without giving its inclination source; qualitatively inconsistent with all three photometric values.
- RM 21 $\pm$ 7 (flagged off-axis), tau 1.8 ms at 1491 MHz (Scott 2025); log M_HI 9.65 (Roxburgh).

## Unavailable quantities
- 20171020A: offset (fly's-eye), Re.
- 20210410D: Re.
- RM: 20190714A, 20200906A, 20210807D, 20210410D (none published).
