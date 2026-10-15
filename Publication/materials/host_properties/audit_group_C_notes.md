# Audit group C: DSA-110 hosts in Sharma et al. 2024 / Law et al. 2024

Hosts: 20220207C, 20220307B, 20220310F, 20220319D, 20220418A, 20220509G, 20220825A, 20220914A, 20220920A, 20221012A, 20221219A, 20230124A, 20230307A.

Primary sources opened: Sharma et al. 2024 (Nature 635, 61; LaTeX source `papers/sharma2024_src/manuscript.tex`, Extended Data Tables 1-2, Supplementary Tables 2-3), Law et al. 2024 (ApJ 967, 29; Tables 1-3, Sec 3.2), Connor et al. 2025 (Nat. Astron.; Supplementary Table), Sharma et al. 2023 (Table 2), Connor et al. 2023 (Table 1), Ravi et al. 2025 (AJ; Tables 1, 3, Sec 3), Faber et al. 2024 (arXiv:2405.14182; Tables 1-3), Sherman et al. 2023/2024 (Tables 1-2), Verdi et al. 2026 (arXiv:2609.30616; Table A2), Hussaini et al. 2025 (Table 1), Bhardwaj et al. 2024 (Table 1). All bibliography keys already exist in `host_properties.bib`; `audit_group_C.bib` therefore adds no entries.

Nickname map used to read the Sharma 2024 tables (`preamble.tex`): Zach = 20220207C, Alex = 20220307B, Whitney = 20220310F, Mark = 20220319D, Quincy = 20220418A, Jackie = 20220509G, Ansel = 20220825A, Elektra = 20220914A, Etienne = 20220920A, Juan = 20221012A, Nihari = 20221219A, Leonidas = 20230124A, Phineas = 20230307A.

## Global findings

### Redshifts

- All 13 redshifts are spectroscopic pPXF fits. Sharma 2024 Extended Data Table 1 quotes $z$ to 4 decimals with errors; Law 2024 Table 3 quotes 6 decimals without errors. Instruments (Sharma 2024 Supplementary Table 2): Keck-I/LRIS for 20220207C, 20220310F, 20220418A, 20220509G, 20220825A, 20220914A, 20221012A, 20221219A; Keck-II/DEIMOS for 20220307B, 20230124A, 20230307A; Palomar P200/DBSP for 20220319D, 20220920A.
- The seed picks cited the Connor 2025 compilation for every $z$. Nine are rounded or differ at the 4th decimal (0.043, 0.479, 0.0111, 0.622, 0.1138, 0.1585, 0.284, 0.094, 0.271); all are replaced by the Sharma 2024 ED Table 1 values. Only 20220509G and 20220825A were already correct to 4 decimals.
- **20220307B settled**: $z = 0.2481$. Law 2024 Table 3 gives 0.248123 (Hale/DBSP, $\ge 3$ emission lines) and Sharma 2024 gives $0.2481 \pm 0.0001$ from an independent Keck-II/DEIMOS spectrum (2023 Jul 13 UT). The 0.2507 in Connor 2025 (copied by Hussaini 2025 and Verdi 2026) is a compilation entry with no spectrum behind it; $cz$ differs by about 620 km/s, i.e. a transcription error. Recorded as `corrected` plus a `conflict` row for the Connor value.
- Minor $z$ conflicts beyond quoted errors (same spectrum, different reduction): 20220207C (Sharma 0.0433 vs Law 0.0430), 20220418A (Sharma 0.6214 vs Law 0.622000, which looks rounded). 20221219A: Sharma 0.5530 (LRIS) vs Faber 2024 0.554 (DEIMOS, 3 d.p., arXiv only). In all three the Sharma 2024 value is recommended (newest refereed, 4-decimal precision).

### DM

- `dm_obs`: Connor 2025 values equal the Sherman 2024 structure-optimized DMs; Law 2024 Table 1 lists S/N-maximizing DMs that differ by $\le 0.2$ pc cm$^{-3}$. One conflict: 20221012A, where Connor 2025 and Sherman 2023 list 442.2 but Sherman 2024 (441.2 $\pm$ 0.26) and Law 2024 (441.08 $\pm$ 0.7) agree on 441.2; recommend 441.2. Verdi 2026 lists 590.0 and 610.15 for 20230124A/20230307A (Connor: 590.6, 608.9); Connor is refereed and retained.
- `dm_mw`: Law 2024 Table 1 NE2001 ISM values kept for the ten Law hosts; Connor 2025 NE2001 for the three later hosts. Connor's DM_MW column is NE2001 ISM only (their DM_ex = DM_obs $-$ DM_MW exactly; the 30 pc cm$^{-3}$ halo is added in their fits, not in the table). Connor's numbers differ from Law's by up to 6% (e.g. 76.0 vs 79.3, 128.2 vs 135.7, 139.8 vs 133.3), within the 10% NE2001 error Connor adopts; not treated as a conflict. YMW16 values added as `new` rows where published (20220319D: 187.7, Ravi 2025; 20220509G: 46 and 20220914A: 47, Connor 2023).
- 20220319D: NE2001 (133.3) exceeds DM_obs (110.95); Ravi 2025 Sec 3 discusses alternatives (pulsar DMs 46.7-65.3, H$\alpha$ EM 89).

### Stellar mass and SFR: Sharma 2024 vs Law 2024

Both papers use Prospector with the same model family (non-parametric continuity SFH, 7 bins, Chabrier IMF, redshift free within 1%, emcee; Law 2024 Sec 3.2 explicitly follows Sharma 2023). The differences are in the data and the dust/AGN treatment, not the code or SFH prior:

1. Law 2024 fits catalog photometry (PS1 grizy, Legacy, 2MASS, ALLWISE, WIRC JH where available) with no UV; it states that "the recent star formation history is not well constrained in a few cases ... primarily due to the lack of UV photometry".
2. Sharma 2024 re-measures aperture photometry on all archival images (GALEX FUV/NUV, SDSS, PS1, BASS, MzLS, 2MASS, ALLWISE; Supplementary Table 3), adds new spectra (e.g. Keck-II/DEIMOS for 20220307B) and routinely includes a Nenkova et al. 2008 AGN dust-torus component when IR data warrant (AGN column in ED Table 2: Y for all group-C hosts except 20220307B).
3. Sharma 2024 tabulates the posterior median and 68% interval; it is the newer refereed paper and explicitly re-fits the Law 2024 hosts (its tables cite `law2023deep` for the data).

Recommendation: adopt Sharma 2024 for `logmstar` and `sfr` throughout. Sharma 2024's SFR timescale: ED Table 2 has two columns, `SFR_20Myr` (7th data column) and `SFR_100Myr` (8th data column), both in $M_\odot$ yr$^{-1}$; the table note defines them as the SFR averaged over the last 20 and 100 Myr. The Table-1 target is `SFR_100Myr`; `SFR_20Myr` is recorded as a separate row with `method=prospector_sed_20Myr`.

Mass conflicts (beyond quoted errors): 20220307B (10.14 vs 10.04), 20220509G (10.70 vs 10.79; Sharma 2023 11.13), 20220914A (9.24 vs 9.48; Sharma 2023 9.99), 20220319D (Ravi 2025 9.93 vs 10.10/10.16). The other six Law hosts agree within errors.

SFR$_{100\,\rm Myr}$ conflicts: 20220207C (0.71 vs 1.16), 20220310F (0.97 vs 4.25), 20220319D (0.42 vs 0.21; Ravi 2025 1.8), 20220418A (1.31 vs 29.41), 20220825A (7.91 vs 1.62), 20220920A (1.62 vs 4.10), plus the superseded Sharma 2023 values for 20220509G (0.08) and 20220914A (1.45). 20220307B, 20220509G, 20220914A (0.46 vs 0.72, borderline) and 20221012A agree within errors. The 20220418A case is the clearest example of the UV problem: without GALEX the Law fit put $29^{+18}_{-11}$ $M_\odot$ yr$^{-1}$ into the youngest bin; with GALEX NUV Sharma finds 1.3.

### $r$-band magnitudes

Sharma 2024 ED Table 1 lists the PATH-filter magnitude (PS1 $r$ or BASS $r$), stated to be MW-extinction corrected. Two conflicts with Law 2024 Table 2:

- 20220418A: Sharma BASS $r$ = 21.22 $\pm$ 0.19 and PS1 $r$ = 21.27 $\pm$ 0.18 vs Law (Legacy) 22.22 $\pm$ 0.20. Exactly 1.0 mag apart; recommend 21.22.
- 20220825A: Sharma PS1 $r$ = 20.74 $\pm$ 0.16 vs Law (PS1) 19.36 $\pm$ 0.16, with the identical error. Sharma's full PS1 SED ($g$ = 21.36, $r$ = 20.74, $i$ = 20.34, $z$ = 19.97, $y$ = 19.56) is monotonic, whereas 19.36 would be brighter than $i$; recommend 20.74. Caveat: the 1.38 mag gap is close to $A_r = 1.24$ for $E(B-V) = 0.543$, so one entry may simply be uncorrected for extinction.

Caveat on the extinction correction: for the two heavily reddened hosts 20220207C ($E(B-V)$ = 0.813) and 20220319D (0.806), the Sharma ED Table 1 magnitudes (16.15, 13.25) are identical to the Law PS1 catalog values, so the correction stated in the table note may not have been applied to those entries. For 20220319D the seed had the SDSS aperture value 13.19 (Supp. Table 3); corrected to the tabulated host magnitude 13.25 (PS1 $r$) for consistency with the other hosts.

### Inclination

Bhardwaj 2024 Table 1 parametric-profile inclinations confirmed for nine hosts; Ravi 2025 Table 3 gives $23 \pm 3$ deg for 20220319D (PS1 $i$ isophotes). No inclination or axis ratio is published anywhere for 20220418A, 20221219A, 20230124A, 20230307A (checked Sharma 2024 tables, Law 2024, Faber 2024, Verdi 2026, Hussaini 2025).

## Per-host notes

**20220207C (Zach)** - $z$ = 0.0433 (Keck-I/LRIS). Mass agrees between papers (9.95 vs 9.91); SFR conflict (0.71 vs 1.16). Edge-on disk (Law Sec 4.1); Bhardwaj inclination 80 deg. RM = 162.48 rad m$^{-2}$ (Sherman 2023). $m_r$ = 16.15 (PS1) in both papers despite $E(B-V)$ = 0.81 (see caveat).

**20220307B (Alex)** - $z$ = 0.2481 settled (DBSP + DEIMOS; Connor's 0.2507 is wrong). Only group-C host fitted without an AGN component by Sharma 2024. Mass conflict (10.14 vs 10.04); SFR consistent (3.17 vs 2.71). Mass-weighted age differs strongly (3.64 vs 0.65 Gyr), another sign of the different photometric input. RM = $-947$ rad m$^{-2}$.

**20220310F (Whitney)** - $z$ = 0.4780 (Keck-I/LRIS; Law 0.477958). Mass consistent (9.98 vs 9.97); SFR conflict (0.97 vs 4.25). Large offset uncertainty ($4.4^{+19.9}_{-3.1}$ kpc).

**20220319D (Mark)** - $z$ = 0.0112 (P200/DBSP; Ravi 2025 0.0111 $\pm$ 0.0004 consistent). Host IRAS 02044+7048, barred spiral SBa with LINER-like nucleus (Ravi 2025). Three Prospector fits: Sharma 2024 (10.10, SFR 0.42), Law 2024 (10.16, 0.21), Ravi 2025 (9.93, 1.8; photometry-only MAP) - Sharma recommended. NE2001 exceeds DM_obs. Inclination 23 $\pm$ 3 deg (Ravi 2025). $m_r$ corrected from 13.19 (SDSS) to 13.25 (PS1).

**20220418A (Quincy)** - $z$ = 0.6214 (Keck-I/LRIS; Law lists 0.622000). Mass consistent (10.26 vs 10.31); SFR conflict 1.31 vs 29.41 (UV-driven). $m_r$ conflict 21.22 vs 22.22. No inclination published. RM = 6.1 $\pm$ 7.5 rad m$^{-2}$.

**20220509G (Jackie)** - $z$ = 0.0894 (Keck-I/LRIS; absorption-line redshift; also Sharma 2023). Quiescent early-type in Abell 2311 (2MASX J18504127+7014359). Mass conflict 10.70 (Sharma 2024) vs 10.79 (Law) vs 11.13 (Sharma 2023, superseded); SFR 0.25 vs 0.23 consistent, Sharma 2023 0.08 superseded. YMW16 = 46 (Connor 2023). $\tau$ = 0.08 ms at 1.4 GHz. Note: the user prompt cites Faber 2024 for 20220509G; Faber 2024 is the 20221219A paper, and no Faber paper covers 20220509G.

**20220825A (Ansel)** - $z$ = 0.2414 (Keck-I/LRIS). Mass consistent (10.01 vs 9.95); SFR conflict (7.91 vs 1.62; dusty, $A_{V,\rm young}$ = 1.5). $m_r$ conflict 20.74 vs 19.36 (recommend 20.74). RM = 750 rad m$^{-2}$. Highest DM_host in the group.

**20220914A (Elektra)** - $z$ = 0.1139 (Keck-I/LRIS; seed 0.1138 corrected). Mass conflict 9.24 (Sharma 2024) vs 9.48 (Law) vs 9.99 (Sharma 2023, superseded); SFR 0.46 vs 0.72 borderline-consistent; Sharma 2023 1.45 superseded. $m_r$ = 20.05 (PS1, Supp. Table 3 and Law); BASS 20.20. No polarization detected (Connor 2023). YMW16 = 47.

**20220920A (Etienne)** - $z$ = 0.1582 (P200/DBSP; seed 0.1585 corrected). Mass consistent (9.87 vs 9.81); SFR conflict (1.62 vs 4.10). $m_r$ = 19.46 (PS1 / Legacy). RM = $-830$ rad m$^{-2}$.

**20221012A (Juan)** - $z$ = 0.2847 (Keck-I/LRIS; seed 0.284 corrected). Quiescent, strong AGN component ($f_{\rm AGN}$ = 0.79). Mass and SFR consistent between papers. DM_obs conflict 442.2 vs 441.2 (recommend 441.2). RM = 166 rad m$^{-2}$.

**20221219A (Nihari)** - $z$ conflict 0.5530 (Sharma, LRIS) vs 0.554 (Faber 2024, DEIMOS, arXiv only); recommend 0.5530. Mass 10.21 and SFR 1.78 are the same Prospector fit in both papers. $m_r$ = 22.61 (BASS; Faber 22.6). Heavily scattered ($\tau$ = 19.2 ms at 1.4 GHz) with two intervening galaxies (Faber 2024). No inclination published.

**20230124A (Leonidas)** - $z$ = 0.0939 (Keck-II/DEIMOS; seed 0.094 corrected). Only Sharma 2024 has SED values (9.46, SFR 0.75); Athukoralalage 2026 compiles the same mass. RM = 1084 rad m$^{-2}$, $\tau < 7.3$ ms at 1 GHz (Verdi 2026, arXiv). No catalog host name, no inclination.

**20230307A (Phineas)** - $z$ = 0.2706 (Keck-II/DEIMOS; seed 0.271 corrected). Transition-class host, 10.76, SFR 0.46 (Sharma 2024 only). RM = $-476.67$ rad m$^{-2}$, $\tau$ = 0.39 ms (Verdi 2026, arXiv). No catalog host name, no inclination.
