# Audit group D notes (round 2)

Hosts: 20230626A, 20230628A, 20231120A, 20231123B (Sharma et al. 2024); 20231220A, 20240119A, 20240213A, 20240215A, 20240229A (Connor et al. 2025 new sources); 20190523A (DSA-10, Ravi et al. 2019); 20220912A (repeater, Ravi et al. 2023 / Hewitt et al. 2024).

Audit date: 2026-10-07. Files written: `audit_group_D.csv` (202 rows), `audit_group_D.bib` (2 new keys: `makarov2014`, `bhusare2026`). New papers/extracts placed in `papers/`: `connor2025_natastron_ESM_MOESM1.pdf/.txt` (Nature Astronomy supplementary file), `li2025_2507.01270.txt` (text of the arXiv PDF; the published IOP HTML `li2025_2507.01270.html` was already present), `frb20220912A_radio2026_2607.05950.pdf/.txt`.

Primary sources opened: Sharma et al. 2024 Nature SI (`sharma2024_nature_SI_MOESM1.txt`, Supplementary Tables 1, 2, 3, 6) and arXiv source; Connor et al. 2025 arXiv v2 (Tables 1, 2) plus the published ESM Supplementary Table 1; Hussaini et al. 2025 Table 1; Verdi et al. 2026 Table A2; Li et al. 2025 Table 1 (arXiv PDF and IOP HTML, identical numbers); Ravi et al. 2019 Table 1 / Methods; Ravi et al. 2023 Table 2; Hewitt et al. 2024 Sec. 4-5; Leung et al. 2025 Sec. 2; Bhandari et al. 2022 Tables A3/A4; Heintz et al. 2020 Tables 1, 3, 6.

## Verdict summary

| FRB | z | DM_obs / DM_MW | log M* | SFR | m_r | offset | inclination | Re |
|---|---|---|---|---|---|---|---|---|
| 20230626A | confirmed 0.327 (Sharma, DBSP) | confirmed 451.2 / 39.2 | confirmed 10.44 | confirmed 0.98 (100 Myr) | confirmed 19.97 BASS r | confirmed 3.71 | new 27.5 (Li+2025) | confirmed 1.96 |
| 20230628A | confirmed 0.1265 (=Sharma 0.127) | confirmed 345.15 / 39.1 | confirmed 9.29 | **corrected**: 0.14 is also the 100 Myr value (+0.04/-0.03); seed said 20 Myr | confirmed 19.94 | confirmed 3.85 | new 30.7 (Li+2025) | confirmed 2.43 |
| 20231120A | **conflict resolved: 0.0368** (Sharma LRIS) | confirmed 438.9 / 43.8 | confirmed 10.40 | confirmed 0.40 (100 Myr) | confirmed 15.41 | confirmed 0.87 | **conflict**: 81 (Leung, approximate) vs 71.7 (Li+2025) | confirmed 5.48 |
| 20231123B | **conflict (minor)**: 0.2625 (Connor) vs 0.2621+/-0.0001 (Sharma); recommend 0.2621 | confirmed 396.7 / 40.2 | confirmed 11.04 | confirmed 4.85 (100 Myr) | confirmed 18.53 | confirmed 14.40 | new 65.2 (Li+2025) | confirmed 14.57 |
| 20231220A | confirmed 0.3355 (DBSP spectrum per ESM) | confirmed 491.2 / 49.9 | **new 10.35+/-0.10** (Li+2025 Bagpipes) | **new 5.88+/-3.06** (Bagpipes) | confirmed 20.49 PS1 r | unavailable | new 47.9 | unavailable |
| 20240119A | confirmed 0.3760 spec (Table 1 "0.3700 photometric" is stale) | confirmed 483.1 / 37.9 | **new 9.56+/-0.14** | **new 2.51+/-1.54** | confirmed 21.20 BASS r | unavailable | new 31.7 | unavailable |
| 20240213A | confirmed 0.1185 (DBSP) | confirmed 357.4 / 40.1 | **new 10.13+/-0.11** | **new 4.94+/-2.65** | confirmed 19.33 | unavailable | new 70.1 | unavailable |
| 20240215A | confirmed 0.2100 (DBSP) | confirmed 549.5 / 48.0 | **new 9.34+/-0.10** (quality 1, no NIR) | **new 2.60+/-1.39** | confirmed 20.36 | unavailable | new 41.0 | unavailable |
| 20240229A | confirmed 0.2870; **corrected** z_type photo -> spec (LRIS 2024 May 3 per ESM); source Connor not Hussaini | confirmed 491.15 / 37.9 | **new 10.08+/-0.20** | **new 1.80+/-2.08** | confirmed 21.12 | unavailable | new 23.8 | unavailable |
| 20190523A | confirmed 0.660+/-0.002 | confirmed 760.8 / 37 | confirmed 11.07 (Ravi Prospector; Li Bagpipes 10.76) | confirmed upper limit <1.3 (no usable later SED SFR) | confirmed 21.91 PS1 Kron | confirmed 27.2 | unavailable | confirmed 3.28 |
| 20220912A | confirmed 0.0771 | confirmed 219.46 / 125.24 (repo NE2001) | confirmed 10.0 (Ravi Prospector; Li Bagpipes 9.60) | **new SED SFR 0.92+/-0.56** (Li+2025 Bagpipes); H-alpha >0.1 and radio 2.5 confirmed; new local radio 1.3 (Bhusare+2026) | confirmed 19.65 PS1 | confirmed 0.8 (VLBI) | confirmed 51 (Ravi); also 43/46 (Bhardwaj+2024), 42.8 (Li+2025) | confirmed 2.2 |

## Task 3: FRB 20231120A redshift

**Recommendation: adopt $z = 0.0368 \pm 0.0001$ (spectroscopic, Keck-I/LRIS). Confidence: high.**

(a) Which value each DSA-110 paper uses:

- Sharma et al. 2024 (Nature 635, 61): $z = 0.0368 \pm 0.0001$, Supplementary Table 1 (host row 20231120A / Mizu), from a Keck-I/LRIS spectrum taken 2023 Dec 17 UT (B400/3400 + R400/8500, 1.0 arcsec slit, 0.80 arcsec seeing, Supplementary Table 2), redshift from a pPXF fit. The same SI tabulates $D_L = 163.69$ Mpc and $M_r^0 = -20.41$ for $m_r = 15.41$, both of which are only consistent with $z = 0.0368$ (at $z = 0.0700$, $D_L \approx 320$ Mpc and $M_r \approx -21.9$). The host was also classified as a spiral with visible arms in archival imaging; at $z=0.0368$ its $R_e = 5.5$ kpc and $\log M_* = 10.4$ are normal for a 2MFGC flat spiral.
- Connor et al. 2025 (Nature Astronomy; arXiv 2409.16952 v2 of 24 Jun 2025 and the public `frbsample_connor0924.csv`): Table 1 lists 0.0700 with no dagger (i.e. labelled spectroscopic) and no error. The published ESM Supplementary Table 1 lists spectroscopic observations only for the nine *new* hosts (20240213A, 20240215A, 20220831A, 20231220A, 20240229A, 20240119A, 20230814A, 20240123A, 20230521B); 20231120A is not among them, so Connor et al. did not take an independent spectrum - the value is a transcription/placeholder, numerically close to the 2MPZ photometric redshift (0.066). No discussion of a discrepancy anywhere in the paper or ESM.
- Hussaini et al. 2025 (ApJL 993, L27): Table 1 lists 0.0700 and cites "Sharma et al. (2024)" as the reference, although Sharma gives 0.0368. Clearly copied from the Connor compilation. No discussion.
- Li et al. 2025 (ApJ, refereed): Table 1 lists 0.0368 with reference 36 = Sharma et al. 2024.
- Athukoralalage et al. 2026 (arXiv 2608.23868): Table 1 lists 0.0368.
- Verdi et al. 2026 (arXiv 2609.30616, the most recent DSA-110 paper): Table A2 lists 0.0368 and uses it in the DM budget (DM_host = 322 pc cm$^{-3}$). No comment on the Connor value.

No DSA-110 paper discusses the discrepancy.

(b) Independent spectroscopic redshift: none exists.

- NED (queried 2026-10-07, `datasearch?objname=2MFGC+07440&search_type=Redshifts`): exactly one redshift entry, $cz = 19902$ km/s, $z = 0.066387$, refcode 2014ApJS..210....9B (Bilicki et al. 2014, 2MPZ photometric redshift catalogue), "no unc. reported", qualifier "Uncertain origin"; NED main-table flag `PUN` (photometric, uncertain). This is the 0.0664 value quoted in the protocol; it is not a spectrum.
- SDSS DR18 SkyServer: `SpecObj` within 0.5 arcmin of (143.98453, +73.28497): no rows (field is outside the SDSS spectroscopic footprint).
- DESI DR1 and EDR (`desi_dr1.zpix`, `desi_edr.zpix` via NOIRLab Data Lab TAP, 36 x 11 arcsec box around the host): no rows.
- LAMOST via VizieR cone (30 arcsec): V/164 (DR5), V/156 (DR7), V/165 (DR8/9 LRS): no rows (V/153 timed out; Dec +73 is north of the LAMOST footprint anyway).
- 2MRS (VizieR J/ApJS/199/26) and 6dFGS (VII/259): no rows (6dF is southern; 2MRS is limited to $K_s \le 11.75$ and this galaxy has $K_s = 12.4$).
- HyperLEDA (VII/237): PGC 2754635 present at 0.4 arcsec, no velocity in the VizieR table.

So the Keck/LRIS spectrum of Sharma et al. 2024 is the only spectroscopic redshift of this galaxy.

(c) Most recent refereed paper: Li et al. 2025 (ApJ) adopts 0.0368; Hussaini et al. 2025 (ApJL, later in 2025) lists 0.0700 but is a compilation that mis-cites Sharma. Under the protocol rule (compilations do not override the measurement paper) the value is 0.0368. The two 2026 arXiv papers (Athukoralalage, Verdi) also use 0.0368.

Consequences: all Sharma-derived quantities ($\log M_*$, SFR, offset 0.87 kpc, $R_e$ 5.48 kpc, $D_L$) are at $z=0.0368$ and remain valid. If 0.0700 were used, offsets and sizes would scale by $\times 1.9$ and masses by $\times 3.6$. The repo already uses 0.0368.

## Task 2: gap hosts - searches and results

Search tools and queries (all run 2026-10-07; raw logs in the session):

1. ADS: the `https://ui.adsabs.harvard.edu/search/q=%22FRB%2020240119A%22` style UI pages are JavaScript-rendered and the ADS API needs a token, so ADS could not be queried programmatically from this machine; the arXiv full-text search and Semantic Scholar were used instead (same literature coverage for 2024-2026 preprints).
2. arXiv search (`https://arxiv.org/search/?query=%22FRB+2024xxxxA%22&searchtype=all`, title/abstract/comments full-text): "FRB 20231220A", "FRB 20240119A", "FRB 20240213A", "FRB 20240215A", "FRB 20240229A", "FRB 20190523A", "DSA-110 host galaxies": 0 results each (these FRBs only appear inside tables, which arXiv does not index). "FRB 20220912A": 30 results, all burst-physics papers except Ravi+2023 (2211.09049), Hewitt+2024 (2312.14490), Sheikh+2024 (2312.07756) and the new **Bhusare et al. 2026 (2607.05950)** radio-continuum paper (see 20220912A below). Also `export.arxiv.org/api/query?search_query=all:"FRB 2024xxxxA"`: 0 hits for all five gap hosts.
3. Semantic Scholar API (`api.semanticscholar.org/graph/v1/paper/search?query=...`): query "FRB 20220912A host galaxy" returned 474 results; the first 30 contain no new host-SED paper for any group-D host other than the ones above (relevant 2025-2026 compilations seen: Glowacki+2025 PASA 2506.23403, Bernales-Cortes+2025 A&A 2501.14063, 2606.10538 PATH update, 2606.25214 redshift-dependent DM_host, none of which add primary SED values for DSA-110 hosts). Per-FRB queries "FRB 20231220A", "FRB 20240119A", "FRB 20240213A", "FRB 20240215A", "FRB 20240229A" (short forms "FRB 2024xxxx" too) were rate-limited (HTTP 429) on every attempt (4-6 retries each, two sessions). Two broader queries did complete: "FRB 20190523A host galaxy" (473 results; top 30 contain no 20190523A-specific paper after Bhandari+2022, only compilations such as 2511.01195 and 2606.25214) and "DSA-110 FRB host galaxies" (16 results; the only DSA-110 host papers are Sharma+2023/2024, Law+2024, Ravi+2023, Connor+2023, Faber+2024, Verdi+2026 - no new DSA-110 host catalogue with SEDs for the 2024 hosts).
4. Targeted follow-up of the compilations named in the protocol: Gordon et al. 2025 (CRAFT only), Athukoralalage et al. 2026 (compiled masses; for 20230628A/20231120A only, copied from Sharma), Leung et al. 2025 (compiled masses; no DSA new hosts), Glowacki et al. 2025 (CRAFT), Verdi et al. 2026 (DM/RM/scattering only, no SED), Hussaini et al. 2025 (z/DM only). The Connor et al. 2025 Nature Astronomy ESM (downloaded) contains only the spectroscopy log (Supp. Table 1) and spectra figures; no photometry or SED fits.
5. Round-1 `host_properties_long.csv` already contained Li et al. 2025 rows for CRAFT hosts but not for the DSA hosts; the DSA rows of Li Table 1 were extracted now.

Result: the only published stellar masses / SFRs for the five Connor et al. 2025 hosts are the **Bagpipes** fits of **Li et al. 2025 (ApJ, arXiv 2507.01270, Table 1)**. Method: Bagpipes (Carnall+2018), BC03 SPS, delayed-tau SFH, Calzetti dust, Cloudy nebular emission; photometry = GALEX + SDSS/DESI-LS + Pan-STARRS Kron + 2MASS/WISE where available; IMF not stated (Bagpipes/BC03 default is Kroupa). SFR is the model's current SFR (not a 100 Myr average). Li also tabulates a PetroFit Sersic axis ratio converted to inclination with $q_0 = 0.2$, and $A = r_e^2/\cos i$. Values entered with `verdict=new`:

| FRB | log M* | SFR (Msun/yr) | log sSFR | i (deg) | A (kpc^2) | quality |
|---|---|---|---|---|---|---|
| 20231220A | 10.35 +/- 0.10 | 5.88 +/- 3.06 | -9.57 +/- 0.26 | 47.87 | 261.75 +/- 1.05 | 0 |
| 20240119A | 9.56 +/- 0.14 | 2.51 +/- 1.54 | -9.14 +/- 0.32 | 31.73 | 444.51 +/- 5.91 | 0 |
| 20240213A | 10.13 +/- 0.11 | 4.94 +/- 2.65 | -9.44 +/- 0.29 | 70.12 | 93.64 +/- 0.15 | 0 |
| 20240215A | 9.34 +/- 0.10 | 2.60 +/- 1.39 | -8.91 +/- 0.28 | 40.97 | 121.18 +/- 0.34 | 1 (no NIR) |
| 20240229A | 10.08 +/- 0.20 | 1.80 +/- 2.08 | -9.79 +/- 1.55 | 23.77 | 196.57 +/- 1.49 | 0 |

No DSA-110 Prospector values, no $r$-band magnitudes beyond Connor Table 2, no offsets and no $R_e$ exist for these five hosts. `re_kpc` could be back-derived from Li's $A$ and $i$ (noted in the CSV, not entered as a value).

Spectroscopy for the five (Connor ESM Supp. Table 1): 20231220A, 20240213A, 20240215A: P200/DBSP 2024 Apr 4 UT; 20240119A: Keck-I/LRIS 2024 Mar 31 UT; 20240229A: Keck-I/LRIS 2024 May 3 UT. All five redshifts are therefore spectroscopic; the daggers on 20240119A (0.3700) and 20240229A (0.2870) in Connor Table 1 are stale. For 20240119A Table 2 gives 0.3760, which Verdi+2026 and Li+2025 adopt; recommend 0.3760.

### 20220912A

- SED SFR: none in Ravi+2023 (their Prospector fit constrained the mass; SFR reported only as H-alpha $> 0.1$ Msun/yr in the central arcsec$^2$). Bhardwaj+2024 Table 1 lists "SFR_0-100Myr = 2.5" citing Ravi+2023, but that number is the Hewitt+2024 radio SFR. **Li et al. 2025 Bagpipes: SFR = 0.92 +/- 0.56 Msun/yr, log M* = 9.60 +/- 0.19 (quality 0)** - the only published SED SFR; entered as `new`.
- New 2026 paper found: Bhusare et al. 2026 (arXiv 2607.05950; ApJ per Semantic Scholar) - VLA 4-26 GHz continuum source at the FRB position, radio SFR 1.3 +/- 0.1 Msun/yr (alpha = -0.73; 0.9 +/- 0.1 for alpha = -0.5), $\Sigma_{SFR} \approx 13$ Msun/yr/kpc$^2$; it is a local star-forming knot, not a host-integrated SFR. Recorded for completeness.
- Inclination: Ravi+2023 51 deg (half-light isophote, convention unstated); Bhardwaj+2024 43 +/- 1 (parametric) and 46 +/- 1 (isophote); Li+2025 42.78 (Sersic, $q_0 = 0.2$). The three $q_0=0.2$ values agree.
- DM_MW: no paper tabulates NE2001; repo value 125.24 kept. Li+2025 DM_exc = 99.03 implies 120.4 - different NE2001 evaluation; flagged.

### 20190523A

- SFR: Ravi+2019 upper limit $< 1.3$ Msun/yr ([OII]); Heintz+2020 $< 0.09$ (H-beta). Li+2025 Bagpipes SFR = 30.84 +/- 48.40 (quality 1) is unconstrained and not usable. **No later SED SFR exists** - tabulate as $< 1.3$.
- Mass: Ravi+2019 Prospector 11.07 +/- 0.06; Heintz+2020 CIGALE 10.79; Bhandari+2022 CIGALE (Kron) 11.67; Li+2025 Bagpipes 10.76 +/- 0.20. Keep Ravi (pick).
- Host association remains insecure (PATH $P(O|x) = 0.82$ in Bhandari+2022; excluded by Gordon+2023). No inclination/axis ratio anywhere.

## Per-host notes (Task 1)

**20230626A (Fatima).** All Sharma Supp. Table 6 values confirmed (Class T, log M* 10.44, SFR100 0.98 +0.33/-0.24, SFR20 1.04, tm 5.41, log Z -0.72, A_V 0.29, (u-r) 2.39, d 3.71 +14.44/-2.61, Re 1.96). m_r 19.97 +/- 0.19 BASS r (Supp. Table 1); E(B-V) printed as 0.806 is a typo (SFD ~0.03 at b = 40 deg; the same 0.806 appears for 20231220A in Connor Table 2). z 0.327 from P200/DBSP 2023 Aug 26. Li+2025 Bagpipes: 10.62 +/- 0.13, 0.70 +/- 2.44, i = 27.51. RM -576.71 (Verdi, arXiv only).

**20230628A (Tildis).** Correction: Supp. Table 6 has both SFR20 = 0.14 +0.02/-0.01 and SFR100 = 0.14 +0.04/-0.03; round 1 only recorded the 20 Myr column and the seed method says `prospector_sed_20Myr`. The pick value is unchanged (0.14) but should be labelled 100 Myr with +0.04/-0.03. Mass 9.29 +/- 0.03 (Leung/Athukoralalage compile 9.26). z: Connor 0.1265 vs Sharma 0.127 (3 d.p.) - same measurement. Li+2025: 9.43 +/- 0.13, SFR 1.40 +/- 0.55, i = 30.70.

**20231120A (Mizu).** Redshift: see Task 3. Host = 2MFGC 07440 / PGC 2754635 (flat-galaxy-catalogue member, J-band a/b >= 3). Inclination conflict: Leung+2025 quote "i ~ 81 deg (Sharma et al. 2024)" in the text of Sec. 2 while Sharma tabulate no inclination (Leung used it only to exclude the host from their face-on sample; method unstated); Li+2025 Table 1 gives 71.74 (PetroFit Sersic b/a, $q_0=0.2$, i.e. b/a ~ 0.37); HyperLEDA logR25 = 0.41 (b/a = 0.39 -> 70 deg); SGA-2020 b/a = 0.295 (-> 77 deg). Recommend Li+2025 71.7 deg as the literature value (only refereed, documented measurement); the host is edge-on in every convention. Sharma SED values confirmed (log M* 10.40, SFR100 0.40, SFR20 0.70, tm 7.54, log Z -1.12, A_V 0.54/0.72, (u-r) 2.69, d 0.87, Re 5.48). m_r 15.41 +/- 0.01 (Supp. Table 1; Supp. Table 3 lists BASS r 15.41 +/- 0.21, PS1 r 15.39 +/- 0.01). Scattering 21.9 ms at 1 GHz (Verdi) with an intervening galaxy at 34 kpc. No RM.

**20231123B (Squanto).** z: Connor/Hussaini/Verdi 0.2625 vs Sharma 0.2621 +/- 0.0001 (Li copies 0.2621). Minor conflict; recommend Sharma 0.2621 (measurement paper with error; D_L 1346.36 Mpc matches). All SED values confirmed (log M* 11.04 +/- 0.01, SFR100 4.85 +0.22/-0.21, SFR20 2.78, tm 4.86, log Z +0.19, A_V 0.37/0.48, (u-r) 3.11, d 14.40 +6.78/-4.77, Re 14.57). m_r 18.53 +/- 0.20 BASS r. Li+2025: 10.99 +/- 0.12, 9.15 +/- 6.64, i = 65.18. RM 379.37 (Verdi).

**20231220A (Gemechu).** Connor Table 1/2 confirmed (z 0.3355 DBSP; DM 491.2; NE2001 49.9; m_r 20.49 +/- 0.04 PS1 r; E(B-V) 0.806 as printed, implausible). New: Li+2025 Bagpipes 10.35 +/- 0.10, 5.88 +/- 3.06, i 47.87. No offset/Re/class published. RM 37.78, tau 0.641 ms (Verdi). Foreground: NGC 2523 halo at ~70 kpc, UGC 4279 ~30 kpc (Connor Methods); Verdi: intervening galaxy at 78 kpc.

**20240119A (Nikhil).** z 0.3760 spectroscopic (LRIS 2024 Mar 31, ESM); Connor Table 1 "0.3700 (photometric)" is inconsistent with their own Table 2 and ESM - Hussaini copied 0.3700; Verdi and Li use 0.3760. m_r 21.20 +/- 0.02 BASS r. New: Li+2025 9.56 +/- 0.14, 2.51 +/- 1.54, i 31.73. No RM; tau < 3.132 ms.

**20240213A (Annie).** z 0.1185 (DBSP 2024 Apr 4); DM 357.4; NE2001 40.1; m_r 19.33 +/- 0.04 BASS r. New: Li+2025 10.13 +/- 0.11, 4.94 +/- 2.65, i 70.12 (edge-on-ish). RM -233.91 +/- 69.61; tau 6.533 ms.

**20240215A (Bubble).** z 0.2100 (DBSP 2024 Apr 4); DM 549.5; NE2001 48.0; m_r 20.36 +/- 0.02 BASS r. New: Li+2025 9.34 +/- 0.10, 2.60 +/- 1.39 (quality 1, no NIR), i 40.97. RM 374.32 +/- 65.21; tau 6.775 ms.

**20240229A (Casey).** z 0.2870 - spectroscopic (LRIS 2024 May 3, ESM) despite the dagger in Connor Table 1; seed source should be Connor+2025 (Hussaini is a compilation) and z_type spec. DM 491.15; NE2001 37.9 (Connor Table 1 DM_ex 453.20 is a rounding slip, 453.25). m_r 21.12 +/- 0.06 BASS r. New: Li+2025 10.08 +/- 0.20, 1.80 +/- 2.08, i 23.77. RM -6.15; tau < 0.196 ms.

**20190523A.** Ravi+2019 Table 1 confirmed (z 0.660 +/- 0.002, DM 760.8 +/- 0.6, NE2001 37, log M* 11.07 +/- 0.06, SFR < 1.3, tau 1.4 ms at 1 GHz). m_r 21.91 +/- 0.15 (Bhandari+2022 PS1 Kron) confirmed; alternatives 22.1 (Ravi) and 22.01 (Heintz PSF). Offset 27.2 +/- 22.6 kpc, Re 3.28 kpc (Bhandari). See Task 2 for SFR/mass status.

**20220912A.** Ravi+2023 Table 2 confirmed (z 0.0771, log M* 10.0 +/- 0.1, SFR > 0.1 H-alpha, Re 2.2 kpc, i 51, A_V 0.5, log Z -0.8, BPT composite); m_r 19.65 PS1 catalogue; offset 0.8 kpc from the EVN position (Hewitt+2024). DM_MW 125.24 (repo NE2001; no paper tabulates it). New SED SFR from Li+2025 (0.92 +/- 0.56) and a local radio SFR from Bhusare+2026 (1.3 +/- 0.1). Sheikh+2024 (ATA) has no host information.

## Anything odd

- Sharma Supp. Table 1 E(B-V) = 0.806 for 20230626A and Connor Table 2 E(B-V) = 0.806 for 20231220A are both implausible for $|b| > 30$ deg and identical; probably a placeholder that was never replaced. The tabulated extinction-corrected magnitudes look consistent with E(B-V) ~ 0.03-0.04, so m_r is kept.
- Connor Table 1 photometric daggers on 20240119A and 20240229A contradict Table 2 and the ESM spectroscopy log (and 20240119A has two different redshifts, 0.3700 and 0.3760, in the same paper).
- Hussaini+2025 cite Sharma+2024 for z = 0.0700 (20231120A) and 0.3700 (20240119A); both are Connor-compilation values, not Sharma's.
- Li+2025 Table 1 for 20231120A has SFR = 0.00 +/- 0.00 and log sSFR = -18.98 +/- 9.42: unconstrained, ignore.
- Li+2025 DM_exc for 20220912A (99.03) implies NE2001 = 120.4 whereas the repo NE2001 is 125.24.
- Verdi et al. 2026 and Athukoralalage et al. 2026 are arXiv only (noted in every row that uses them).
