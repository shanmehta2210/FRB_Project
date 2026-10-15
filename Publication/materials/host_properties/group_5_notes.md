# Group 5 sweep notes (12 DSA-110 hosts)

Hosts: 20221219A, 20230124A, 20230307A, 20230626A, 20230628A, 20231120A, 20231123B, 20231220A, 20240119A, 20240213A, 20240215A, 20240229A.

All 12 are DSA-110 non-repeaters. Seven (2022-12 to 2023-11 discoveries) are in Sharma et al. (2024, Nature 635, 61); the other five (20231220A, 20240119A, 20240213A, 20240215A, 20240229A) were first published in Connor et al. (2025, Nat. Astron., arXiv:2409.16952) and re-tabulated in Hussaini et al. (2025, ApJL 993, L27). Verdi et al. (2026, arXiv:2609.30616; the "in prep." paper mentioned in the brief now has an arXiv posting) gives DM_host, RM, RM_host, B_parallel and scattering for all 12. Faber et al. (2024, arXiv:2405.14182) is a dedicated paper on 20221219A. Leung et al. (2025, ApJL 991, L25) and Athukoralalage et al. (2026, arXiv:2608.23868) are compilations that add an inclination estimate for 20231120A and independent DM_host values.

## Papers read

| source_short | paper | arXiv | what it gave |
|---|---|---|---|
| sharma2024 | Sharma et al. 2024, Nature 635, 61 | 2409.16964 (LaTeX source + Nature SI MOESM1) | coords, PATH, z, m_r, E(B-V), full Prospector SED outputs, multi-band photometry, spectroscopy log, offset d, R_e, D_L, rest-frame colours, class, BPT/morphology remarks |
| connor2025 | Connor et al. 2025, Nat. Astron. | 2409.16952 | DM_obs, DM_ex, DM_MW (NE2001), z for all 12; coords, PATH, m_r, E(B-V) for the 5 new hosts; sightline notes for 20221219A and 20231220A |
| hussaini2025 | Hussaini et al. 2025, ApJL 993, L27 | 2506.04186 | FRB coords (deg), z, DM_obs, DM_exgal (observer frame, NE2001 ISM subtracted) |
| verdi2026 | Verdi et al. 2026 (DSA-110 host environments) | 2609.30616 | z, DM, rest-frame DM_host, RM, RM_host, B_par,host, tau at 1 GHz, Galactic latitude (Table A2); intervening-galaxy notes |
| faber2024 | Faber et al. 2024 (FRB 20221219A sightline) | 2405.14182 (preprint only; no journal DOI found) | burst properties, PATH, z, logM*, SFR, r, DM budget, |RM| limit, environment |
| leung2025 | Leung et al. 2025, ApJL 991, L25 | 2507.16816 | 20231120A inclination (InclinationZoo), 20230628A logM* and DM_host^Macquart |
| athukoralalage2026 | Athukoralalage et al. 2026 | 2608.23868 | rest-frame DM_host (NE2025 + YT20) and compiled logM* for 20230124A, 20230628A, 20231120A |
| moustakas2023 / bilicki2014 / mitronova2004 | SGA-2020, 2MPZ, 2MFGC catalogues (via SIMBAD/NED) | - | host identifiers, b/a, PA, grz mags and photo-z for the bright 20231120A host only |

Law et al. (2024, ApJ 967, 29) and Sherman et al. (2024, ApJ 964, 131) were downloaded and searched; neither contains any Group 5 host (their samples end before 2022-12). Ravi et al. (2025) was not found as a separate host paper (the DSA-110 sample paper cited in Hussaini is Sharma/Connor); nothing was fetched under that name. Connor et al. (2025) "Supplementary Table 1" (spectroscopy setups for the 9 new hosts) is not in the arXiv PDF, so the instrument/date of the spectra for 20231220A, 20240119A, 20240213A, 20240215A, 20240229A is not recorded.

## Conventions used in group_5.csv

- Sharma values are cited to the published Nature Supplementary Information table numbers (Supp. Table 1 basic properties, Supp. Table 2 spectroscopy, Supp. Table 3 photometry, Supp. Table 6 derived properties); the arXiv source labels the same tables Extended Data Tables 1-2 and Supplementary Tables.
- Sharma SFR is given twice (20 Myr and 100 Myr averages). `sSFR_Gyr` is the tabulated specific SFR (2-decimal, Gyr^-1); `log_sSFR` is derived here two ways (from the tabulated sSFR, and from SFR_100Myr/M*); they differ by up to ~0.3 dex because of rounding and because the tabulated sSFR is a separate posterior quantity.
- `offset_kpc` from Sharma is the projected galactocentric offset "d" with asymmetric errors that propagate the 90% localization ellipse; the errors are large for the low-z hosts.
- `Re_kpc` from Sharma has no quoted uncertainty.
- Connor `DM_exgal` = DM_obs - DM_NE2001 (halo not subtracted); Connor adopt a fixed 30 pc cm^-3 halo in their MCMC. Hussaini `DM_exgal` is the same definition; `DM_MW` for Hussaini rows is derived here as DM_obs - DM_exgal.
- Verdi `DM_host` and `RM_host` are rest frame; tau is observer frame scaled to 1 GHz with alpha = 4; upper limits carry `method=upper_limit`.
- Leung `DM_host` (DM_host^Macquart) and Athukoralalage `DM_host` are rest frame, using NE2001+YT20 and NE2025+YT20 respectively for the Milky Way.
- `frb_ra_deg` / `frb_dec_deg` rows carry the radio localization; `host_ra_deg` / `host_dec_deg` the PATH host centroid.

## Per-host notes

### FRB 20221219A ("Nihari")
- Sharma: z = 0.5530 +/- 0.0002 (Keck-I/LRIS 2023-06-13), P_host 0.99 (0.996 in Fig. 1), m_r = 22.61 +/- 0.20 (BASS), logM* = 10.21 (+0.03/-0.04), SFR_100 = 1.78, class SF, d = 5.92 (+14.96/-4.03) kpc, R_e = 5.61 kpc, A_V = 0.99, logZ = -0.44, t_m = 3.59 Gyr. Only PS1 r,i,z,y and SDSS u,g,r photometry (faint host).
- Faber: z = 0.554, logM* = 10.20 (+0.03/-0.04; abstract has the errors swapped), SFR 1.78 (same as Sharma SFR_100), r = 22.6 +/- 0.2, P_host 0.996, tau = 19.2 +/- 2.7 ms at 1.4 GHz, S/N 8.9, intrinsic FWHM 0.12 ms, L/I >= 0.42-0.50, |RM| < 345 rad m^-2, DM budget MW 44 +/- 10, IGM+ICM 608, IGH1 10, IGH2 14, DM_host = 62 (+48/-43) rest frame.
- Verdi: DM_host = 255.87 (+129.91/-136.14) rest frame (no foreground-halo/cluster correction, hence much larger than Faber's 62), tau(1 GHz) = 73.758 +/- 10.372 ms, no RM, b = 33.5 deg.
- Environment: foreground cluster J171039.6+713427 (z ~ 0.14) and two foreground galaxies at b = 43 and 36 kpc; scattering attributed to an intervening CGM cloudlet, not the host. Exclude from any host-ISM scattering analysis.
- Connor/Hussaini: DM_obs 706.70, DM_MW 44.4, DM_ex 662.3.
- No published inclination or axis ratio. Faber calls the host "Milky Way-like".

### FRB 20230124A ("Leonidas")
- Sharma: z = 0.0939 +/- 0.0002 (Keck-II/DEIMOS 2023-05-14), P_host 0.99 (0.999), m_r = 19.10 +/- 0.20, logM* = 9.46, SFR_100 = 0.75, SFR_20 = 0.09, sSFR 0.16 Gyr^-1, class SF, d = 1.14 (+3.73/-0.81) kpc, R_e = 1.53 kpc, A_V = 0.00, f_AGN = 0.28 (large AGN fraction in the SED model). Note the tabulated logM* and several other parameters have +/-0.00 errors, i.e. the posterior is unrealistically narrow.
- Connor: DM_obs 590.60, DM_MW 38.5, DM_ex 552.1. Verdi tabulates DM = 590.00 (conflict of 0.6 pc cm^-3; both recorded).
- Verdi: DM_host = 477.74 (+48.87/-42.85) rest frame (very large), RM = 1084.47 +/- 12.93, RM_host = 1293.90 +/- 21.03 rad m^-2, B_par = 3.44 uG, tau < 7.305 ms.
- Athukoralalage: DM_host = 477 (+30/-41), logM* 9.46.
- No inclination/axis ratio published.

### FRB 20230307A ("Phineas")
- Sharma: z = 0.2706 +/- 0.0001 (DEIMOS 2023-05-14), P_host 0.97 (0.968), m_r = 19.48 +/- 0.19, logM* = 10.76, SFR_100 = 0.46, SFR_20 = 0.30, sSFR 0.01 Gyr^-1, class transition, d = 7.93 (+1.39/-1.27) kpc (localization is unusually tight, 0.2"), R_e = 3.98 kpc, t_m = 7.07 Gyr, f_AGN = 0.24. BPT: Seyfert; tau-linear declining SFH. Photometry: GALEX NUV, PS1 grizy, SDSS griz, BASS g r, MzLS z, WISE W1-W3.
- Connor: DM_obs 608.90, DM_MW 37.6, DM_ex 571.3. Sharma Fig. 1 label says DM = 580; Verdi tabulates DM = 610.15. All recorded.
- Verdi: DM_host = 392.09 (+89.58/-87.35), RM = -476.67 +/- 0.03, RM_host = -761.65 +/- 12.14, B_par = -2.71 uG, tau = 0.390 (+0.018/-0.012) ms at 1 GHz.
- No inclination/axis ratio published.

### FRB 20230626A ("Fatima")
- Sharma: z = 0.3270 +/- 0.0001 (P200/DBSP 2023-08-26), P_host 0.99 (0.999), m_r = 19.97 +/- 0.19, logM* = 10.44 +/- 0.04, SFR_100 = 0.98, SFR_20 = 1.04, class transition, d = 3.71 (+14.44/-2.61) kpc, R_e = 1.96 kpc, f_AGN = 0.68. E(B-V) is tabulated as 0.806, which is not credible for this high-latitude (b = 40 deg) field; it looks like a copy error (the same 0.806 appears for several unrelated hosts). Photometry: PS1 grizy, BASS g r, MzLS z, WIRC J (2023-07-31), WISE W1 W2.
- Connor: DM_obs 451.20, DM_MW 39.2, DM_ex 412.0.
- Verdi: DM_host = 157.87 (+81.10/-85.15), RM = -576.71 +/- 12.42, RM_host = -1015.47 +/- 31.09, B_par = -10.54 uG, tau < 0.414 ms.
- No inclination/axis ratio published.

### FRB 20230628A ("Tildis")
- Sharma: z = 0.1270 +/- 0.0001 (DEIMOS 2023-07-13), P_host 0.95 (0.953), m_r = 19.94 +/- 0.20, logM* = 9.29 +/- 0.03, SFR_100 = 0.14, SFR_20 = 0.14, sSFR 0.04 Gyr^-1, class SF, d = 3.85 (+2.94/-2.55) kpc, R_e = 2.43 kpc, A_V = 0.22. Photometry: GALEX FUV/NUV, PS1 grizy, BASS g r.
- Connor: DM_obs 345.15, DM_MW 39.1, DM_ex 306.0; z listed as 0.1265 (vs 0.1270 in Sharma).
- Verdi: DM_host = 191.21 (+51.49/-51.07), RM = -17.39 +/- 2.24, RM_host = -6.18 +/- 9.19, B_par ~ 0, tau < 0.591 ms.
- Leung: logM* = 9.26 ("published" value; slightly different from Sharma's 9.29), DM_host^Macquart = 224. Included in Leung's face-on sample, i.e. Leung judged it to have i < 70 deg (no number given). Athukoralalage: logM* 9.26, DM_host = 193 (+35/-48).
- No numerical inclination/axis ratio published.

### FRB 20231120A ("Mizu")
- Redshift conflict. Sharma SI / arXiv source: z = 0.0368 +/- 0.0001 from their own Keck-I/LRIS spectrum (2023-12-17), with D_L = 163.69 Mpc consistent with that z; Verdi and Athukoralalage adopt 0.0368. Connor Table 1 and Hussaini Table 1 list z = 0.0700 (no comment). NED/2MPZ has a photometric z = 0.0664 for the host. The seed file uses 0.07. The two values give very different DM_host (Verdi 322 at 0.0368 vs implied ~290 at 0.07) and a factor ~3.5 in luminosity distance; the Sharma spectroscopic value is the only primary measurement and should be adopted, but the ApJL text should flag the discrepancy.
- Host is 2MFGC 07440 = LEDA 2754635 = 2MASX J09355633+7317059: a 2MASS Flat Galaxy Catalog object (J-band a/b >= 3). SGA-2020 ellipse gives b/a = 0.295 (0.644 x 0.190 arcmin), PA = 31 deg. Leung et al. (2025) quote i ~ 81 deg (InclinationZoo on SDSS/PS1) and exclude it as edge-on. Sharma note spiral arms are visible. This is the one Group 5 host with a published inclination.
- Sharma: P_host 0.99 (0.987), m_r = 15.41 +/- 0.01 (BASS), logM* = 10.40, SFR_100 = 0.40, SFR_20 = 0.70, sSFR 0.01 Gyr^-1, class transition, d = 0.87 (+1.36/-0.60) kpc, R_e = 5.48 kpc, A_V,young 0.54, A_V,old 0.72, logZ = -1.12, t_m = 7.54 Gyr, M_r^0 = -20.41. Full PS1/BASS/MzLS/2MASS/WISE photometry.
- Connor: DM_obs 438.90, DM_MW 43.8, DM_ex 395.1. Verdi: DM_host = 322.28 (+43.33/-39.88) rest frame, tau = 21.904 (+2.010/-2.523) ms at 1 GHz, no RM; intervening galaxy at b = 34 kpc. Athukoralalage: DM_host = 333 (+21/-26).
- Large DM_host and large scattering for an edge-on disk: consistent with a long path through the host ISM.

### FRB 20231123B ("Squanto")
- Sharma: z = 0.2621 +/- 0.0001 (LRIS 2024-02-14), P_host 0.90 (0.908; lowest in the group, host-normalised offset is large), m_r = 18.53 +/- 0.20, logM* = 11.04, SFR_100 = 4.85, SFR_20 = 2.78, sSFR 0.03 Gyr^-1, class transition, d = 14.40 (+6.78/-4.77) kpc, R_e = 14.57 kpc (largest in the sample), M_r^0 = -21.85, (u-r)^0 = 3.11. Spiral arms visible. Photometry PS1 grizy, BASS g r, MzLS z, WISE W1-W4.
- Connor: DM_obs 396.70, DM_MW 40.2, DM_ex 356.5; z 0.2625 (Sharma 0.2621).
- Verdi: DM_host = 106.51 (+63.07/-65.40), RM = 379.37 +/- 18.10, RM_host = 603.99 +/- 32.87, B_par = 7.02 uG, tau < 1.666 ms.
- No inclination/axis ratio published.

### FRB 20231220A ("Gemechu")
- Connor Table 2: FRB 08:15:38.09 +73:39:35.70 (+/-0.7", 0.5"), host 08:15:38.47 +73:39:34.76, P_host 0.99, z = 0.3355 (spec), m_r = 20.49 +/- 0.04 (PS1 r), E(B-V) tabulated 0.806 (same suspicious placeholder value). DM_obs 491.20, DM_MW 49.9, DM_ex 441.3. Sightline passes NGC 2523 (b ~ 70 kpc) and within ~30 kpc of UGC 4279.
- Hussaini: DM_exgal 441.32. Verdi: DM_host = 185.52 (+87.89/-92.25), RM = 37.78 +/- 1.91, RM_host = 89.31 +/- 9.39, B_par = 1.07 uG, tau = 0.641 (+0.073/-0.077) ms; intervening galaxy at b = 78 kpc.
- No stellar mass, SFR, offset, R_e or morphology published anywhere found.

### FRB 20240119A ("Nikhil")
- Connor Table 2: FRB 14:57:52.12 +71:36:42.33 (+/-1.3", 0.7"), host 14:57:53.01 +71:36:40.99, P_host 0.98, z = 0.3760, m_r = 21.20 +/- 0.02 (BASS), E(B-V) 0.024. Connor Table 1 instead lists z = 0.3700 with a photometric-redshift dagger; Hussaini copies 0.3700; Verdi adopts 0.3760. Treat the redshift as uncertain (possibly photometric) until a spectrum is published.
- DM_obs 483.10, DM_MW 37.9, DM_ex 445.2. Verdi: DM_host = 155.79 (+86.06/-89.85), no RM, tau < 3.132 ms.
- No stellar mass, SFR, offset, R_e or morphology published.

### FRB 20240213A ("Annie")
- Connor Table 2: FRB 11:04:40.39 +74:04:31.40 (+/-1.0", 0.6"), host 11:04:40.27 +74:04:28.86, P_host 0.99, z = 0.1185 (spec), m_r = 19.33 +/- 0.04 (BASS), E(B-V) 0.074. DM_obs 357.40, DM_MW 40.1, DM_ex 317.3.
- Verdi: DM_host = 209.52 (+50.42/-49.44), RM = -233.91 +/- 69.61, RM_host = -267.92 +/- 87.81, B_par = -1.79 uG, tau = 6.533 (+0.773/-0.614) ms at 1 GHz (strongly scattered).
- No stellar mass, SFR, offset, R_e or morphology published.

### FRB 20240215A ("Bubble")
- Connor Table 2: FRB 17:53:45.90 +70:13:56.50 (+/-0.8", 0.5"), host 17:53:45.97 +70:13:56.18, P_host 0.99, z = 0.2100 (spec), m_r = 20.36 +/- 0.02 (BASS), E(B-V) 0.037. DM_obs 549.50, DM_MW 48.0, DM_ex 501.5.
- Verdi: DM_host = 352.81 (+75.54/-72.95), RM = 374.32 +/- 65.21, RM_host = 527.60 +/- 97.05, B_par = 2.07 uG, tau = 6.775 (+0.987/-1.202) ms (strongly scattered).
- No stellar mass, SFR, offset, R_e or morphology published.

### FRB 20240229A ("Casey")
- Connor Table 2: FRB 11:19:56.05 +70:40:34.40 (+/-0.8", 0.6"), host 11:19:56.48 +70:40:34.67, P_host 0.99, z = 0.2870, m_r = 21.12 +/- 0.06 (BASS), E(B-V) 0.018. Connor Table 1 flags z = 0.2870 as photometric (dagger); Hussaini and Verdi list it without comment. Treat as photo-z.
- DM_obs 491.15, DM_MW 37.9, DM_ex 453.2. Verdi: DM_host = 234.82 (+85.29/-88.10), RM = -6.15 +/- 0.06, RM_host = 10.16 +/- 8.32, B_par 0.16 (errors tabulated as -0.01/-0.41, sic), tau < 0.196 ms.
- No stellar mass, SFR, offset, R_e or morphology published.

## Cross-paper conflicts (summary)

| FRB | quantity | values |
|---|---|---|
| 20231120A | z | 0.0368 (Sharma SI, Verdi, Athukoralalage) vs 0.0700 (Connor, Hussaini) vs 0.0664 photo-z (2MPZ) |
| 20240119A | z | 0.3760 (Connor Table 2, Verdi) vs 0.3700 photometric (Connor Table 1, Hussaini) |
| 20240229A | z type | photometric per Connor Table 1; not flagged elsewhere |
| 20230124A | DM_obs | 590.60 (Connor/Hussaini) vs 590.00 (Verdi) |
| 20230307A | DM_obs / DM_exgal | 608.90 / 571.3 (Connor) vs 610.15 (Verdi) vs DM_exgal label 580 (Sharma Fig. 1) |
| 20230628A | logM* | 9.29 +/- 0.03 (Sharma) vs 9.26 (Leung, Athukoralalage) |
| 20221219A | DM_host | 62 (+48/-43) with foreground cluster/halo budget (Faber) vs 255.87 (+130/-136) without (Verdi) |
| 20221219A | logM* errors | +0.03/-0.04 (Faber Table 2, Sharma) vs +0.04/-0.03 (Faber abstract) |
| 20230626A, 20231220A | E(B-V) | 0.806 tabulated; implausible, likely placeholder |

## Inclination / axis-ratio coverage

Only 20231120A has anything published: i ~ 81 deg (Leung et al. 2025, InclinationZoo) and catalogue b/a = 0.295 (SGA-2020), plus its 2MFGC flat-galaxy selection. Sharma et al. note visible spiral arms for 20231120A and 20231123B only. No other Group 5 host has a published inclination, axis ratio, Sersic index or position angle.
