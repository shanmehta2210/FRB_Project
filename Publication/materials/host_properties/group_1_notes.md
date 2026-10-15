# Group 1 literature sweep: notes

Hosts: 20171020A, 20180924B, 20190102C, 20190608B, 20190714A, 20191001A, 20190523A, 20181112A, 20200906A.
Output: `group_1.csv` (905 rows, 18 sources), `group_1.bib`, PDFs in `papers/`.

## Papers read (all fetched successfully from arXiv)

| key | reference | arXiv | used for |
|---|---|---|---|
| mahony2018 | Mahony et al. 2018, ApJL 867, L10 | 1810.04354 | 20171020A discovery/host search |
| leewaddell2023 | Lee-Waddell et al. 2023, PASA 40, e029 | 2305.17960 | 20171020A host revisited (the "Li et al. 2023 PASA 40 29" in `seed_51.csv` is this paper) |
| bannister2019 | Bannister et al. 2019, Science 365, 565 | 1906.11476 | 20180924B |
| prochaska2019 | Prochaska et al. 2019, Science 366, 231 | 1909.11681 | 20181112A |
| ravi2019 | Ravi et al. 2019, Nature 572, 352 | 1907.01542 | 20190523A |
| cho2020 | Cho et al. 2020, ApJL 891, L38 | 2002.12539 | 20181112A HTR/RM |
| macquart2020 | Macquart et al. 2020, Nature 581, 391 | 2005.13161 | DM, z, offsets (5 hosts) |
| day2020 | Day et al. 2020, MNRAS 497, 3335 | 2005.13162 | scattering, RM, HTR widths |
| bhandari2020a | Bhandari et al. 2020, ApJL 895, L37 | 2005.13160 | CIGALE M*, SFR, photometry (4 hosts) |
| bhandari2020b | Bhandari et al. 2020, ApJL 901, L20 | 2008.12488 | 20191001A |
| chittidi2021 | Chittidi et al. 2021, ApJ 922, 173 | 2005.13158 | 20190608B kinematics / inclination |
| simha2020 | Simha et al. 2020, ApJ 901, 134 | 2005.13157 | 20190608B DM budget |
| heintz2020 | Heintz et al. 2020, ApJ 903, 152 | 2009.10747 | compilation (8 hosts) |
| mannings2021 | Mannings et al. 2021, ApJ 917, 75 | 2012.11617 | HST morphology/offsets (5 hosts) |
| bhandari2022 | Bhandari et al. 2022, AJ 163, 69 | 2108.01282 | compilation + 20200906A discovery |
| gordon2023 | Gordon et al. 2023, ApJ 954, 80 | 2302.05465 | Prospector SED fits (7 hosts). The protocol cites "ApJ 952 122"; the published record is ApJ 954, 80 |
| shannon2025 | Shannon et al. 2025, PASA 42, e036 | 2408.02083 | CRAFT ICS catalogue (DM, S/N, fluence, PATH) |
| bhardwaj2024 | Bhardwaj et al. 2024, Nature (2408.01876) | 2408.01876 | inclinations and adopted host values in the repo's own paper |

Not fetched: none. Gordon et al. 2025, Sharma et al. 2024, Law et al. 2024, Hussaini et al. 2025, Connor et al. 2023 are not relevant to these nine hosts (no DSA/CHIME/MeerKAT hosts in this group except 20190523A, which is in none of them).

## Systematic conflicts across compilations

- **Stellar mass.** Three model families are in play: CIGALE (bannister2019, bhandari2020a, heintz2020, bhandari2022), Prospector (gordon2023; ravi2019 for 20190523A) and ProSpect (leewaddell2023). Gordon+23 Prospector masses are systematically 0.1-0.5 dex higher than the CIGALE values for the same host. `bhardwaj2024` Table 1 adopts Gordon+23 values for the CRAFT hosts.
- **SFR.** CIGALE SFRs in heintz2020/bhandari2022 are long-timescale SED SFRs; gordon2023 gives 0-100 Myr Prospector SFRs, often 2-5x higher (e.g. 20190608B: 0.69 vs 7.03; 20191001A: 8.06 vs 18.28). Discovery papers quote H-alpha or UV SFRs. All are recorded with the indicator in `method`.
- **DM_MW.** NE2001 is the default in all CRAFT papers; YMW16 recorded where quoted. shannon2025 rounds NE2001 values to integers.
- **DM_obs.** Early papers quote the detection DM; day2020 and shannon2025 give refined HTR/structure-optimised DMs (differences 0.5-1.5 pc cm-3).
- **Offsets.** Pre-HST offsets (macquart2020, bhandari2020a) differ from HST F160W offsets (mannings2021) by up to 1 kpc; heintz2020 adopted earlier positions. shannon2025 re-derived FRB positions and separations for the PATH table.
- **DM_host.** Only chittidi2021 (20190608B), bannister2019 (20180924B), mahony2018 (20171020A) and ravi2019 (upper limit) give something that is a measurement or model constraint. `bhardwaj2024` "DM_host" values are observer-frame residuals after subtracting NE2001 and the median Macquart relation, not independent measurements. prochaska2019 assumes 50/(1+z).

## Per-host notes

### FRB 20171020A (host ESO 601-G036)
- Papers: mahony2018, leewaddell2023, bhandari2022 (Table A3), bhardwaj2024.
- z = 0.00867 (HIPASS/6dF; HI). DM = 114.1, NE2001 38 (YMW16 26). Fly's Eye detection; localisation region is tens of arcmin so **no offset or PATH posterior in the usual sense exists**. leewaddell2023 give two association probabilities: 42% from a `p(m_r)` prior alone and 98% from a joint spatial+redshift prior.
- M*: mahony2018 log M* = 8.95 (photometric M/L); leewaddell2023 ProSpect 8.50 (+0.11/-0.08) and 8.64 (alternative setup). bhandari2022 carries the Mahony value. bhardwaj2024 adopt 8.6.
- SFR: mahony2018 0.13 (H-alpha); leewaddell2023 0.09 (ProSpect). sSFR -9.69 (leewaddell2023 derived), -9.84 (bhandari2022).
- Metallicity disagrees: 12+log(O/H) = 8.3 (mahony2018, O3N2 integrated) vs 8.45 (leewaddell2023, HII complex A). Two different apertures.
- Gas: log M_HI 9.36 (mahony2018 HIPASS), 9.2 for the galaxy alone / 9.32 incl. tail (leewaddell2023 ATCA); interacting with ESO 601-G037. W50 = 126 km/s.
- Inclination: bhardwaj2024 69 (parametric) / 65 (isophote). Mahony's HIPASS rotation gives v_rot ~ 80 km/s (not inclination corrected).
- Odd: seed_51 cites "Li et al. 2023 PASA 40 29"; no such paper. It is Lee-Waddell et al. 2023, PASA 40, e029.

### FRB 20180924B (host DES J214425.25-405400.81)
- Papers: bannister2019, macquart2020, day2020, bhandari2020a, heintz2020, mannings2021, bhandari2022, gordon2023, shannon2025, bhardwaj2024.
- z: 0.3214 (bannister2019, macquart2020, shannon2025) vs 0.3212 (heintz2020, gordon2023, bhardwaj2024). DM 361.42 (discovery) / 362.16 (day2020 HTR) / 362.4 (bhandari2020a, shannon2025).
- M*: 10.35 (bannister2019 CIGALE), 10.12 (heintz2020/bhandari2022/mannings2021 CIGALE), 10.39 (gordon2023 Prospector). SFR: 2.9 (bannister2019 H-alpha), <2 (bannister2019 SED), 0.88 (heintz2020 CIGALE), 0.62 (gordon2023). sSFR <-10.0, -10.18, -10.6.
- RM: 14 +/- 1 (bannister2019) vs 22 +/- 2 (day2020 HTR); both recorded. Scattering 0.68 ms at 1.27 GHz (day2020).
- Offset: 3.9 kpc (bannister2019), 3.5 (macquart2020/bhandari2020a), 3.4 (heintz2020/bhandari2022), 3.37 (mannings2021 HST).
- Morphology: bannister2019 GALFIT b/a = 0.62 +/- 0.05, PA -24, n ~ 2; mannings2021 find spiral arms + bar in HST residuals. bhardwaj2024 i = 34/36 deg. Note bannister2019 b/a = 0.62 would imply i ~ 52 deg for a thin disk; this is from ground-based DES r with poorly resolved galaxy.
- Age: bannister2019 t_age 5.55 Gyr (CIGALE), heintz2020 0.38 Gyr (CIGALE "age" column, different definition), gordon2023 t_m 5.63 Gyr.
- DM_host: bannister2019 Table S9 host halo DM 51 +/- 32 (top-hat 1 R_vir model); bhardwaj2024 residual 80.

### FRB 20190102C (host DES J212939.76-792832.5)
- Papers: macquart2020, day2020, bhandari2020a, heintz2020, mannings2021, bhandari2022, gordon2023, shannon2025, bhardwaj2024.
- z: 0.291/0.2913 (macquart2020, bhandari2020a), 0.2912 (heintz2020, gordon2023, shannon2025), 0.2909 (bhardwaj2024). DM 363.6 (macquart2020) vs 364.5 (day2020, bhandari2020a, shannon2025).
- RM: macquart2020 and bhandari2020a tabulate +110 (unsigned, preliminary); day2020 -105 +/- 1 (pulse 2), -128 +/- 7 (pulse 1). Sign conflict; day2020 is the dedicated polarimetry paper.
- M*: 9.50 (bhandari2020a), 9.53 (heintz2020/mannings2021), 9.67 (bhandari2022), 9.69 (gordon2023). SFR 1.5 (bhandari2020a H-alpha), 0.86 (heintz2020 CIGALE), 0.40 (gordon2023). sSFR -9.33 / -9.74 / -10.09.
- m_r: bhandari2020a 22.06 (DES r) vs bhardwaj2024 20.77. Large discrepancy; bhardwaj2024 Table 1 likely quotes a deeper VLT/FORS2 or model magnitude. Both recorded; flag for the paper.
- Offset: 1.5 (+3.4/-1.5) kpc (macquart2020/bhandari2020a), 2.0 +/- 2.2 (heintz2020), 2.26 +/- 2.22 (mannings2021 HST), 2.3 (bhandari2022).
- Re: 5.3 kpc (bhandari2020a GALFIT), 4.4 (heintz2020), 4.87 isophotal / 5.0 GALFIT (mannings2021). Mannings: no clear symmetric structure in residuals. bhardwaj2024 i = 50/50 deg.

### FRB 20190608B (host SDSS J221604.90-075356.0)
- Papers: macquart2020, day2020, bhandari2020a, chittidi2021, simha2020, heintz2020, mannings2021, bhandari2022, gordon2023, shannon2025, bhardwaj2024.
- z = 0.1178 (SDSS; 0.11778 in bhandari2020a). DM 338.7 (macquart2020) / 339.5 (bhandari2020a, shannon2025) / 339.79 (chittidi2021) / 340.05 (day2020). Most scattered CRAFT burst: tau = 3.3 ms at 1.27 GHz; RM 353 +/- 2.
- **Only host in this group with a published kinematic inclination**: chittidi2021 i_gas = 37 +/- 3 deg (KCWI H-beta, thin disk 39.2, thick disk 35.7) and i_stellar = 26 +/- 3 deg from isophote ellipticity e = 0.10 (b/a = 0.90). bhardwaj2024 give 37/36 deg, consistent with the kinematic value and not with the stellar-isophote value. v_circ = 141 +/- 8 km/s (3DBarolo 165); log M_halo,dyn 11.96.
- M*: 10.4 (bhandari2020a, chittidi2021), 10.06 (heintz2020/bhandari2022/mannings2021), 10.56 (gordon2023). SFR: 1.2 (bhandari2020a H-alpha), 0.69 (heintz2020 CIGALE), 7.03 (gordon2023 0-100 Myr). Order-of-magnitude SFR spread; the galaxy hosts an AGN/LINER nucleus with broad H-alpha (bhandari2020a) so H-alpha SFRs are contaminated.
- m_r: 17.65 and 17.38 (bhandari2020a, two apertures), 17.55 (heintz2020), 17.41 (gordon2023/bhardwaj2024); shannon2025 PATH uses 17.15.
- DM_host: chittidi2021 149 +/- 45 (H-alpha EM for ISM, 94 +/- 38, plus halo model); simha2020 ISM 47-117 plus halo 15-41; bhardwaj2024 residual 238.
- Offset: 7.0 (macquart2020), 6.8 (bhandari2020a), 6.6 (heintz2020), 6.52 (mannings2021 HST), 6.5 (bhandari2022). FRB sits in the NW spiral arm on the kinematic major axis.
- Re: heintz2020 2.8 kpc; mannings2021 isophotal 2.65 vs GALFIT disk 7.37 (two-component fit; bulge+bar make single values unreliable, bhandari2020a say GALFIT Re "not a true representation").

### FRB 20190714A (host DES J121555.12-130116.6)
- Papers: heintz2020, mannings2021, bhandari2022, gordon2023, shannon2025, bhardwaj2024. No dedicated discovery paper; first host characterisation is heintz2020.
- z = 0.2365. DM 504.7 (shannon2025 only; no earlier DM tabulation in the papers read). NE2001 39.
- M*: 10.17 (heintz2020/mannings2021), 10.15 (bhandari2022), 10.22 (gordon2023). SFR 0.65 (CIGALE) vs 1.89 (gordon2023). sSFR -10.34 / -9.94.
- m_r: heintz2020 20.69 vs bhandari2022 20.32 / gordon2023 20.34. Both recorded.
- Offset: 1.9 +/- 1.1 (heintz2020), 1.97 +/- 1.18 (mannings2021 HST), 2.7 +/- 1.8 (bhandari2022). Spiral arms seen only in HST residuals.
- Highest 12+log(O/H) in group (9.03). shannon2025 note ~2/3 of extragalactic DM is attributable to foreground halos (Simha+23). bhardwaj2024 i = 51/50 deg; residual DM_host 335 is an overestimate given the foreground excess.

### FRB 20191001A (host DES J213324.44-544454.65)
- Papers: bhandari2020b, heintz2020, mannings2021, bhandari2022, gordon2023, shannon2025, bhardwaj2024.
- z = 0.234 (0.2342 in bhardwaj2024). DM 506.92; NE2001 44, YMW16 31. Detected at 824 MHz (tau 3.3 ms, RM 55.5). Most luminous host in gordon2023.
- **Two candidate hosts at a common redshift.** bhandari2022 PATH gives 0.60 / 0.40; adopted the closer one. shannon2025 re-derived position gives P(O|x) = 0.9995 for the adopted host and 0.00 for the other, and reduces the angular separation from 2.74-2.86 arcsec to 1.18 arcsec. This changes the physical offset from ~11 kpc (heintz2020, bhandari2022, mannings2021 10.49) to ~4.4 kpc if the shannon2025 position is adopted (not published as kpc; only arcsec recorded).
- M*: 10.7 (bhandari2020b), 10.67 (heintz2020/bhandari2022/mannings2021), 10.73 (gordon2023). SFR 8.06 (CIGALE) vs 18.28 (gordon2023); bhandari2020b note persistent 5.5 GHz radio emission consistent with SF.
- Morphology: spiral, FRB in outskirts / on an arm (mannings2021). bhardwaj2024 i = 53/56 deg.

### FRB 20190523A (host PSO J207.0643+72.4708)
- Papers: ravi2019, heintz2020, bhandari2022, gordon2023 (exclusion only). DSA-10 burst; not in CRAFT compilations.
- z = 0.660. DM 760.8; NE2001 37; MW halo assumed 50-80; extragalactic 644-674. Scattering 1.4 ms at 1 GHz. Fluence >280 Jy ms.
- **Association is insecure.** ravi2019 P_chance < 0.10; heintz2020 P_chance 0.073 (Sample C), with a second galaxy S2 at z = 0.363; bhandari2022 PATH 0.82; gordon2023 excluded it (P(O|x) < 0.9) so there is no Prospector fit.
- M*: 11.07 (ravi2019 Prospector, delayed-tau), 10.79 (heintz2020 CIGALE), 11.67 (bhandari2022 CIGALE). Nearly 0.9 dex spread. SFR: <1.3 (ravi2019), <0.09 (heintz2020 upper limit), 0.09 (bhandari2022 tabulates the limit as a value; sSFR -12.74 derived from it). Quiescent / massive.
- Age: ravi2019 6.6 Gyr; heintz2020 0.69 Gyr (CIGALE column). A_V 0.3 (ravi2019).
- Offset: heintz2020 3.79 arcsec, 27 +/- 23 kpc (large because the DSA-10 localisation ellipse is 3x8 arcsec).
- No DM_host measurement; ravi2019 quote <~150 as the possible host column (from scattering argument). Missing: log_sSFR (none published directly; bhandari2022 derived value recorded), m_r from Pan-STARRS 22.1 / LRIS 22.01 / 21.91.

### FRB 20181112A (host DES J214923.66-525815.28)
- Papers: prochaska2019, cho2020, macquart2020, bhandari2020a, heintz2020, bhandari2022 (PATH only), gordon2023, shannon2025.
- z = 0.4755. DM 589.27 (589.265 cho2020; 589.0 bhandari2020a/shannon2025). NE2001 40.2 (bhandari2020a) / 42 (prochaska2019 DM budget) / 40 (shannon2025). RM 10.9 (prochaska2019), 10.5 pulse 1 / 25 pulse 3 (cho2020). Width <40 us.
- Sightline passes 29 kpc from a foreground galaxy FG-181112 (z = 0.3674, log M* 10.69); heintz2020 therefore put it in Sample C (P_chance 0.026), bhandari2022 PATH 0.67 and excluded it; shannon2025 PATH 0.93 (foreground P(O|x) 0.07). gordon2023 kept it (their own criterion) and fitted it.
- M*: 9.41 (prochaska2019), 9.42 (bhandari2020a), 9.6 (heintz2020), 9.87 (gordon2023). SFR 0.6 (prochaska2019/bhandari2020a H-alpha), 0.4 (prochaska2019 SED), 0.37 (heintz2020), 1.54 (gordon2023). sSFR -9.62 / -9.68.
- m_r: 21.73 (prochaska2019 DES), 21.61 (bhandari2020a), 21.68 (heintz2020/gordon2023).
- Offset poorly constrained: 3.1 (+15.7/-3.1) kpc (macquart2020/bhandari2020a), 1.7 +/- 19.2 (heintz2020), 0.28 arcsec. Not in mannings2021 (no HST). Re 0.64 arcsec, n = 1.9 (prochaska2019 FORS2 GALFIT); 8 +/- 3 kpc (bhandari2020a) vs 7.2 +/- 1.7 (heintz2020).
- Not in bhardwaj2024 Table 1; no inclination published. DM_host only as the prochaska2019 assumption (50/(1+z) = 32).

### FRB 20200906A (host DES J035343.61-140538.57)
- Papers: bhandari2022 (discovery of host), gordon2023, shannon2025, bhardwaj2024.
- z = 0.3688. DM 577.8; NE2001 36 (YMW16 38). Fluence 59 (+25/-10) bhandari2022 vs 53 +/- 3 shannon2025. Width 6.0 ms.
- M*: 10.12 (bhandari2022 CIGALE), 10.37 (gordon2023 Prospector). bhardwaj2024 adopt 10.12 (the Bhandari value) here while using Gordon values for other hosts, but take SFR 4.93 from gordon2023; their sSFR -9.43 is therefore a mixed-source number and differs from gordon2023's own -9.68. Both recorded.
- SFR 0.48 (bhandari2022 CIGALE) vs 4.93 (gordon2023 0-100 Myr). H-alpha flux 6.49e-16 (DEIMOS).
- Offset 5.9 +/- 2.0 kpc (bhandari2022); PATH 1.00. Re 7.58 kpc (GALFIT). No HST data; no morphology beyond "star-forming". bhardwaj2024 i = 63/59 deg.
- m_r 20.15 (bhandari2022) vs 19.95 (gordon2023/bhardwaj2024).

## Quantities not found for any host
- Published axis ratio only for 20180924B (bannister2019 b/a = 0.62) and 20190608B (chittidi2021 b/a = 0.90). No GALFIT axis ratios are tabulated in mannings2021 (only Re and offsets).
- No independent inclination outside chittidi2021 (20190608B) and bhardwaj2024.
- Environment class (cluster/group/field) is given only for 20171020A (group / interacting pair) and 20190523A (possible unconfirmed group).
