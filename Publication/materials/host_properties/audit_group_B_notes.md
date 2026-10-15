# Audit group B: CRAFT/ASKAP newer hosts (13)

Hosts: 20220501C, 20220725A, 20221106A, 20230526A, 20230613A, 20230902A, 20230907D, 20231020B, 20231226A, 20240208A, 20240210A, 20240304A, 20240318A.

Companion files: `audit_group_B.csv` (231 rows), `audit_group_B.bib` (three arXiv-only papers checked during the gap search and not in `host_properties.bib`).

## Papers opened (all in `papers/`)

| key | paper | used for |
|---|---|---|
| shannon2025 | Shannon et al. 2025, PASA 42, e036 (arXiv:2408.02083) | Table 1 ($z$, DM, NE2001), Table 4 (PATH), Table 5 (FORS2/HAWK-I photometry), Sec. 5.2.x host descriptions |
| gordon2025 | Gordon et al. 2025, ApJ 993, 119 (arXiv:2506.06453) | Table 1 ($z$), Table 3 (GALFIT $R_e$, $n$, K-band $M_\star$), Table 4 (offsets), Table C2 (20240318A LRIS spectrum) |
| li2025 | Li et al. 2025, arXiv:2507.01270 | Table 2 (Bagpipes $M_\star$, SFR, sSFR, quality flag, PetroFit inclination) |
| leung2025 | Leung et al. 2025, arXiv:2507.16816 | Table 1 (CIGALE / NED-LVS masses) |
| muller2026 | Muller et al. 2026, ApJ 1001, 118 (arXiv:2506.20774) | Table 3 (X-Shooter pPXF $z$, H$\alpha$ SFR, extinction-corrected $m_R$, metallicity, BPT/Balmer flags) |
| pastormarazuela2026 | Pastor-Marazuela et al. 2026, MNRAS 545 (arXiv:2507.05982) | Table 1, Table D1, Table D3 (CIGALE) for the three MeerTRAP hosts |
| glowacki2025 | Glowacki et al. 2025, PASA 42, e157 (arXiv:2506.23403) | Tables 1-3 ($z$, SF flag, GALFIT $b/a$, $i$, offsets) |
| scott2025 | Scott et al. 2025, PASA 42, e133 (arXiv:2505.17497) | Tables 1-3 (HTR DM, scattering, RM) |
| hussaini2025 | Hussaini et al. 2025 (arXiv:2506.04186) | Table 1 ($DM_{\rm exgal}$, used to back out NE2001) |
| roxburgh2025 | Roxburgh et al. 2025, PASA (arXiv:2507.06594) | Table 3 (20240210A radio SFR, HI) |
| athukoralalage2026 | Athukoralalage et al. 2026 (arXiv:2608.23868) | Table 1 (copies Leung masses) |
| baptista2024 | Baptista et al. 2024, ApJ 965, 57 | Table 1 (cross-check of DM, NE2001) |
| james2026, james2026b, dial2026 | James et al. 2026a/b, Dial et al. 2026 (arXiv only) | gap search only; no host properties |

Independent NE2001 values quoted below were computed in this audit with `mwprop` (NE2001p v2.0.5, integrated to 30 kpc) at the Shannon 2025 Table 1 positions; they agree with `master_frb_localization.csv` and with Hussaini 2025 ($DM_{\rm obs}-DM_{\rm exgal}$) to 0.1 pc cm$^{-3}$.

## Cross-cutting findings

1. **Shannon 2025 Table 1 NE2001 column has four tabulation errors.** For 20230526A (printed 50, true 31.9), 20231226A (145 vs 38.1), 20240208A (98 vs 39.0) and 20240304A (30 vs 73.8) the printed value cannot be reproduced at the tabulated coordinates with NE2001 (or YMW16). Hussaini 2025 Table 1 independently implies 38.07 and 73.76 for 20231226A and 20240304A. The other nine Group-B sightlines reproduce to within rounding. Rows carry `verdict=conflict` with the recomputed value recommended (`source_short=this_audit` or `hussaini2025`).
2. **Li 2025 Table 2 (continued pages are headed "Table 5" in the arXiv PDF) was under-read in round 1.** It contains Bagpipes masses and SFRs for 20230526A, 20230902A, 20231226A and 20240210A in addition to the two already picked (20220725A, 20221106A). These are added as `new` rows and fill two of the six "gap" hosts.
3. **Inclinations: two independent measurements exist for four hosts and disagree.** Glowacki 2025 Table 3 (GALFIT on deep VLT/FORS2 imaging, fits from Marnoch et al. in prep) vs Li 2025 Table 2 (PetroFit on DESI Legacy Survey cut-outs). Differences: 20221106A 50.90 vs 39.69; 20230902A 69.67 vs 34.67; 20231226A 44.64 vs 35.63; 20240210A 54.35 vs 51.95 deg. Both convert $b/a$ with the Hubble formula ($q_0=0.2$); I reproduce Glowacki's $i$ from their $b/a$ to <0.2 deg. Recommendation: Glowacki 2025 (deeper, higher-resolution imaging; refereed) for all four; Li 2025 only where it is the sole value (20230526A).
4. **$m_r$:** Shannon 2025 Table 5 FORS2 $R$ magnitudes are *not* Galactic-extinction corrected; Muller 2026 Table 3 re-tabulates the same photometry extinction-corrected. Both are recorded; the Muller value is recommended where it exists. Round-1 picks left $m_r$ empty for 20220725A and 20221106A although Shannon Table 5 has them.
5. **Glowacki 2025 Table 2 "log SFR" column** (-0.98, -0.61, -1.05 for 20230526A/20230902A/20231226A) corresponds to an earlier FURBY draft of the Muller H$\alpha$ values and should not be cited; cite Muller 2026.
6. **Redshift of 20231226A:** Glowacki 2025 prints 0.1539 in Tables 1-3; Shannon, Muller, Leung, Hussaini, Athukoralalage, Li and Gordon give 0.1569/0.1570. Treated as a Glowacki typo; recommended 0.1569 (Shannon X-shooter). Glowacki kpc quantities for this host are ~2 % low as a result.
7. **Offsets for 20240210A and 20240304A:** Gordon 2025 Table 4 (4.38 and 8.73 kpc) vs Glowacki 2025 Table 3 (2.05 and 0.86 kpc) disagree far beyond errors. Gordon's angular offsets (9.0 and 2.26 arcsec) agree with Shannon's PATH $\theta$ for these hosts; Glowacki's do not. Recommend Gordon.
8. **Gordon 2025 K-band masses** ($M/L_K=0.6$, no errors) are systematically 0.5-0.7 dex above the SED masses for star-forming hosts in this group (and 1.4 dex for 20230526A, where the HAWK-I $K_s=17.1$ looks inconsistent with $R=21.1$). Recorded as alternates only; used as the headline mass only where nothing else exists (20240304A, 20240318A), flagged.

## Per-host notes

### 20220501C ($z=0.381$)
- $z$ 0.381 (Shannon Table 1; X-shooter H$\alpha$, [OIII]). DM 449.5(2); NE2001 31 (independent 30.6). All confirmed.
- No stellar mass or SFR anywhere: not in Gordon 2025 (host not imaged), Glowacki 2025 Table 2 (row empty), Li 2025 Table 2 (row has $z$/DM only), Leung 2025, Muller 2026 (not a FURBY target), Athukoralalage 2026. Only photometry: FORS2 $g=21.49(1)$, $I=20.47(1)$ (Shannon Table 5; not extinction corrected). No $r$-band magnitude published. Host is 15 arcsec from a $V=11.6$ star (TYC 7514), not catalogued in WISEA/DES.
- Inclination 36.65(4.83) deg, $b/a=0.81(5)$, offset 2.64(1.84) kpc, $R_e=2.50(10)$ kpc: Glowacki 2025 Table 3, confirmed. SF flag Y (basis not given).
- RM 35.5(3) rad m$^{-2}$, $\tau=0.35(25)$ ms at 864.5 MHz: Scott 2025.

### 20220725A ($z=0.1926$)
- $z$ 0.1926 (Shannon; SOAR), Muller pPXF 0.1926. DM 290.4(3); NE2001 31 (30.8). Confirmed.
- $M_\star$: Li 2025 Bagpipes 10.70(10) confirmed (quality 0). Alternates: Leung NED-LVS 10.71, Gordon K-band 11.26.
- SFR: Li Bagpipes 10.73(4.13) confirmed; Muller H$\alpha$ 3.97(1) added (`halpha`; host is a BPT AGN, so H$\alpha$ may be contaminated; Muller excludes it from the metallicity sample).
- $m_r$ was empty in the seed: FORS2 $R=17.843(4)$ (Shannon T5) / 17.806(4) extinction-corrected (Muller T3). Recommend 17.806.
- $i=58.58(17)$ deg, $b/a=0.548(2)$ (Glowacki T3) confirmed; offset 2.03$^{+1.30}_{-0.97}$ kpc (Gordon T4, R band) confirmed; $R_e$ 4.76 kpc (Gordon). Host WISEA J233315.68-355925.0. RM -26.3(7) (polarisation calibration flagged unreliable), $\tau=2.29(5)$ ms.

### 20221106A ($z=0.2044$)
- $z$ 0.2044 (Shannon), Muller 0.2043, Scott Table 1 misprints 2.044. DM 343.8(8); NE2001 35 (34.7). Confirmed.
- $M_\star$ 10.78(10) Bagpipes (Li, quality 0) confirmed; Gordon K 11.27 alternate. SFR 2.37(2.20) Bagpipes confirmed; Muller H$\alpha$ 0.535(8) added (BPT AGN; strong Balmer absorption).
- $m_r$ empty in seed: FORS2 $R=18.351(9)$ / 18.322(9) corrected. Recommend 18.322.
- **Inclination conflict**: seed 39.69 (Li, PetroFit/Legacy) vs Glowacki 50.90(16) ($b/a=0.65$). Recommend Glowacki. Offset 4.52$^{+2.29}_{-2.10}$ kpc (Gordon R) confirmed. Host WISEA J034649.07-253411.7. RM 444(1), $\tau=0.182(6)$ ms.

### 20230526A ($z=0.1570$)
- $z$ 0.1570 (Shannon X-shooter; Muller 0.1570). DM 361.4(2); Scott Table 1 prints 316.2 (digit transposition).
- **NE2001 conflict**: Shannon prints 50; recomputed 31.9. Recommend 31.9.
- $M_\star$ 9.31(16) CIGALE (Leung T1, flag b) confirmed; Li Bagpipes 9.44(13) (quality 1, no NIR) added, consistent; Gordon K 10.73 rejected (see above).
- SFR: seed 0.359(3) H$\alpha$ (Muller) confirmed; Li Bagpipes 0.59(62) added. Recommend the SED value as table entry with H$\alpha$ as alternate (per protocol SED > H$\alpha$), noting the large Bagpipes error.
- $m_r$ 21.08(1) (Shannon T5) confirmed, though Sec. 5.2.6 says the $R$ image is saturated at the nucleus; Muller 21.03(1) corrected value recommended.
- Inclination 20.63 deg (Li T2) added; only value (Glowacki row empty). Offset 1.76 kpc (Gordon) confirmed. 12+log(O/H)=8.5(1) O3N2. RM 613(2), $\tau=1.16(1)$ ms.

### 20230613A ($z=0.3923$, MeerTRAP)
- All seed values confirmed in Pastor-Marazuela 2026: $z$ 0.3923 (DEIMOS), DM 483.51(1) (Table D1 prints 483.82), NE2001 30 (30.0), CIGALE $\log M_\star=10.13(18)$, SFR 4.5 (0.28 dex), $m_r=20.132(4)$ (DESI DR10), offset 3.4(3.2) kpc. Added: age 2.0(8) Gyr, $R_{1/2}=1.39$ arcsec (DESI; not converted to kpc in paper). No inclination, axis ratio, host name or classification published.

### 20230902A ($z=0.3619$)
- $z$ 0.3619 (Shannon X-shooter; Muller 0.3619). DM 440.1(1); NE2001 34 (34.1). Confirmed.
- **Gap filled**: Li 2025 Bagpipes $\log M_\star=9.58(13)$ (quality 1), SFR 1.77(1.16). Gordon K-band 10.25 alternate (0.67 dex higher). Muller H$\alpha$ SFR 0.72(1) confirmed as alternate.
- $m_r$ 21.522(6) (Shannon) confirmed; 21.491(6) corrected (Muller) recommended.
- **Inclination conflict**: seed 69.67(77) (Glowacki, $b/a=0.39$) vs Li 34.67 (Legacy, 0.46 arcsec host). Recommend Glowacki. Offset 4.44 kpc (Gordon R) confirmed. 12+log(O/H)=8.5(1). RM 164.8(2), $\tau=0.123(2)$ ms.

### 20230907D ($z=0.4638$, MeerTRAP)
- Confirmed: $z$ 0.4638 (LRIS), DM 1030.79(4) (Table D1: 1031.11), NE2001 29 (28.7), CIGALE 10.88(20), SFR 14.91 (0.33 dex), offset 5.7(4.2) kpc. Added age 2.3(1.1) Gyr.
- **$m_r$ conflict within the paper**: Table D3 19.83(9) (authors' pipeline) vs Sec. 4.1.8 DESI DR10 catalogue 19.50(1). Both recorded; Table D3 value kept as seed since it is the one used in the SED. No inclination.

### 20231020B ($z=0.4775$, MeerTRAP)
- Confirmed: $z$ 0.4775 (DEIMOS), DM 952.2(3) (Table D1: 953.84), NE2001 34 (34.4), CIGALE 10.75(30), SFR 13.58 (0.99 dex, poorly constrained), $m_r=21.81(2)$ (text 21.79), offset 1.2(3.3) kpc. Added age 2.2(1.4) Gyr. No inclination.

### 20231226A ($z=0.1569$)
- $z$: Shannon 0.1569 confirmed; **Glowacki 0.1539 conflict** (typo). DM 329.9(1).
- **NE2001 conflict**: Shannon prints 145; recomputed 38.1; Hussaini implies 38.07. Recommend 38.1.
- $M_\star$ 9.59(16) CIGALE (Leung) confirmed; Li Bagpipes 9.71(8) (quality 0) added; Gordon K 10.34 alternate.
- SFR: seed 0.359(5) H$\alpha$ (Muller; strong Balmer absorption) confirmed; Li Bagpipes 2.43(74) added and recommended as the SED table value.
- $m_r$ 18.990(6) confirmed; 18.942(6) corrected (Muller) recommended.
- **Inclination conflict**: seed 44.64(22) (Glowacki, $b/a=0.725$) vs Li 35.63. Recommend Glowacki. Offset 4.84 kpc (Gordon R) confirmed. Host WISEA J102127.29+060634.5, grand-design spiral (Gordon Sec. 4.1). 12+log(O/H)=8.5(1) N2. RM 428.4(3), $\tau=0.1(7)$ ms (unconstrained).

### 20240208A ($z=0.2385$)
- $z$ 0.2385 appears only in Glowacki 2025 Tables 1-2 with no spectrum described (Shannon Table 1 lists $z$ as pending; Sec. 5.2.12 gives SDSS photo-$z$ 0.4(1)). Recorded as confirmed-spec with "unverified provenance"; the SDSS photo-$z$ is recorded as a separate `new` row.
- DM 260.2(3) confirmed. **NE2001 conflict**: Shannon 98 vs recomputed 39.0. Recommend 39.0.
- No mass, SFR, $m_r$, inclination, offset or $R_e$ in any paper (see gap search). Only SDSS $g=22.8$ (Shannon Sec. 5.2.12). Host SDSS J103654.96-005712.2 (coincidence within ~4 arcsec; 14-antenna localisation). Round-1 seed wrongly attributed this host to Gordon 2025. RM -73.7(1.4), $\tau=1.35(25)$ ms.

### 20240210A ($z=0.023686$)
- $z$ 0.023686 (archival SSRS, da Costa et al. 1998, via Shannon). DM 283.73(5); NE2001 31 (28.6; within rounding). Confirmed.
- $M_\star$ 9.67(16) CIGALE (Leung) confirmed; Li Bagpipes 9.71(5) (quality 0) added; Gordon K 10.33 alternate.
- **Gap filled**: SED SFR 0.81(13) (Li Bagpipes, quality 0; log sSFR -9.79). Alternate: Roxburgh 2025 1.4 GHz radio SFR 1.19$^{+0.66}_{-0.43}$. Caveat for both: Seyfert 1 nucleus.
- $m_r$ 14.919(9) (Shannon FORS2 $R$, not extinction corrected) confirmed; no corrected value published (not a FURBY target).
- **Inclination conflict (minor)**: Glowacki 54.35(4) ($b/a=0.6051$) vs Li 51.95. Recommend Glowacki.
- **Offset conflict**: Gordon 4.38(26) kpc (9.0 arcsec; two-component fit) vs Glowacki 2.05(41). Recommend Gordon. Host WISEA J003506.47-281619.1; HI-detected (log $M_{\rm HI}=9.61$). RM -325(1) (edge of beam), $\tau=0.10(3)$ ms.

### 20240304A ($z=0.2423$)
- $z$ 0.2423 (Shannon Sec. 5.2.14, X-shooter; Table 1 lists p). DM 652.6(5).
- **NE2001 conflict**: Shannon 30 vs recomputed 73.8; Hussaini implies 73.76. Recommend 73.8.
- Mass: only Gordon K-band 9.99 (no error) -> `new`, flagged as non-SED. SFR unavailable (X-shooter spectrum exists per Shannon but no flux or SFR published; not in FURBY).
- $m_r$ 20.723(7) (Shannon FORS2 $R$) confirmed. Inclination 51.77(42) deg ($b/a=0.638$, Glowacki) confirmed.
- **Offset conflict**: Gordon 8.73$^{+2.91}_{-3.02}$ kpc (2.26 arcsec, matches PATH) vs Glowacki 0.86(1.57). Recommend Gordon. RM 489.7(8) (Scott Table 3 caption says no RM fit; row lists a value), $\tau=2.51(12)$ ms.

### 20240318A ($z=0.1120$)
- $z$ 0.1120 from Gordon 2025 Table C2 (Keck/LRIS 2025 Jan 2; Shannon Table 1 lists p). DM 256.4(3); NE2001 37 (36.9). Confirmed.
- Mass: only Gordon K-band 10.31 (MMT/MMIRS; no error) -> `new`, flagged. SFR unavailable (Gordon lists emission lines but publishes no fluxes). $m_r$ unavailable (Gordon took LRIS $R$ imaging but tabulates no magnitude; Shannon Table 5 has no row; SDSS $g=19.5$ only).
- No inclination (Glowacki row empty; Gordon fits $b/a$ for the GC-ellipsoid model but does not tabulate it). Offset 2.63$^{+1.19}_{-1.14}$ kpc (Gordon R) confirmed; $R_e$ 4.95 kpc. Host WISEA J100134.45+373659.7. RM -48.5(3), $\tau=0.163(10)$ ms.

## Task 2: gap-host literature search log

Target: any 2025-2026 refereed or arXiv paper publishing a stellar mass, SFR or $r$-band magnitude for 20220501C, 20230902A, 20240208A, 20240304A, 20240318A, 20240210A.

### Queries run

ADS (full-text, via the ADS web UI; the API needs a token that is not available in this environment), sorted by date, 2026-10-07:

| query | hits (bibcodes) | relevant new host data |
|---|---|---|
| `full:"FRB 20220501C" OR full:"FRB 20220501"` | 2025PASA...42...36S, 2024MNRAS.528.1583H (Hoffmann FREDDA), 2023Sci...382..294R, 2022TNSFR1403 | none |
| `full:"FRB 20230902A" OR full:"FRB 20230902"` | 2025ApJ...993..119G, 2025PASA...42...36S, 2025PASA...42....5W (CRACO), 2023TNSFR2149 | none beyond Gordon |
| `full:"FRB 20240208A" OR full:"FRB 20240208" OR full:"J103654.96-005712.2"` | 2025PASA...42...36S, 2024TNSFR2412 | none |
| `full:"FRB 20240304A" OR full:"FRB 20240304"` | 2025ApJ...993..119G, 2025PASA...42..133S, 2025arXiv250801648C (Caleb, FRB 2024030**4B**, different source), 2025PASA...42...36S, 2024TNSFR1587 | none |
| `full:"FRB 20240318A" OR full:"FRB 20240318" OR full:"J100134.45+373659.7"` | 2026PASA...43...74B (Balzan, polarisation only), 2025ApJ...993..119G, 2025PASA...42...36S, 2024TNSFR1587 | none |
| `full:"FRB 20240210A" OR full:"FRB 20240210" OR full:"J003506.47-281619.1"` | 2026arXiv260918086B (Balzan, FRB 20250607A), 2026arXiv260904765J (James 2026b), 2025ApJ...993..119G, 2025MNRAS.538.2785H (Hackstein, local DM), 2025PASA...42...36S, 2024TNSFR1586 | none |
| `abs:"fast radio burst" abs:(host OR hosts) abs:(ASKAP OR CRAFT) abs:("stellar mass" OR "star formation") year:2025-2026` | 2026arXiv260810452J, 2026AAS...24840303K, 2026PASA...43..109J (James 2026a), 2025ApJ...993..119G, 2025ApJ...992..127L | none |

Note: ADS full-text does not return Glowacki 2025, Li 2025, Muller 2026 or Scott 2025 for several of these names even though those papers tabulate the hosts, so ADS full-text coverage of 2025 arXiv/PASA papers is incomplete; those papers were read directly instead.

arXiv (API `search_query=all:` and the listing search `arxiv.org/search/?query=...&searchtype=all`, `curl.exe -A "Mozilla/5.0" -L`), each of `"FRB 20220501C"`, `"FRB 20220501"`, `FRB20220501C` and the same three forms for 20230902A, 20240208A, 20240304A, 20240318A, 20240210A, plus `"WISEA J003506.47-281619.1"`, `"WISEA J100134.45+373659.7"`, `"SDSS J103654.96-005712.2"`, `CRAFT FRB host galaxies stellar mass`, `ASKAP FRB host galaxy catalogue 2025`, `Marnoch FRB host galaxy inclination`, `incoherent sum ASKAP FRB hosts stellar masses`: 0 hits for every FRB-name query (arXiv indexes titles/abstracts only); the topical queries returned only Shannon 2025, Gordon 2025, Glowacki 2025, Li 2025, James 2025/2026a/2026b.

Semantic Scholar Graph API (`/graph/v1/paper/search?query=...&fields=title,year,externalIds,venue&limit=30`), same query set with 8-40 s back-off: most calls returned HTTP 429 (no API key); the calls that succeeded (`"FRB 20220501C"`, `"FRB 20230902A"`, `"FRB 20230902"`) returned generic FRB papers with no match to the host names.

Candidate papers downloaded and grepped for the six names (`papers/`): James, Prochaska & Bera 2025 (arXiv:2502.16791, localisation precision vs host magnitude), Khrykin et al. 2025 (arXiv:2509.08896, magnetic fields), James et al. 2026a (arXiv:2606.10538, PATH update: lists 20220501C, 20230902A, 20240210A, 20240304A with PATH-input $I$/$R$ magnitudes 20.57, 21.52, 15.13, 21.08 only), James et al. 2026b (arXiv:2609.04765: notes 20240210A image omitted from Shannon 2025 and published in Gordon 2025), Dial et al. 2026 (arXiv:2609.21177, CRAFT HTR2 polarimetry: burst properties only, refers to an "updated CRAFT catalogue" that has no host table), Caleb et al. 2025 (arXiv:2508.01648: FRB 20240304**B**, MeerKAT, unrelated), Balzan et al. 2026 (arXiv:2609.18086: burst morphology only). None contains a stellar mass, SFR or $r$ magnitude for the six hosts. TNS discovery reports (2022TNSFR1403, 2023TNSFR2149, 2024TNSFR1586/1587/2412) contain burst parameters only. The Marnoch et al. host-morphology paper referenced by Glowacki 2025 is still "in prep" as of all 2026 papers checked.

### Outcome per host

| host | $\log M_\star$ | SFR | $m_r$ | status |
|---|---|---|---|---|
| 20220501C | none | none | none ($g=21.49$, $I=20.47$ FORS2) | **nothing published**; $i=36.65$ deg exists |
| 20230902A | 9.58(13) Bagpipes (Li 2025, Table 2, quality 1) | 1.77(1.16) Bagpipes; 0.72 H$\alpha$ (Muller) | 21.52 / 21.49 | **filled** (found in Li 2025 Table 2, missed in round 1) |
| 20240208A | none | none | none ($g=22.8$ SDSS) | **nothing published** beyond $z$ |
| 20240304A | 9.99 K-band only (Gordon) | none | 20.723 FORS2 $R$ | no SED mass/SFR exists |
| 20240318A | 10.31 K-band only (Gordon) | none | none ($g=19.5$ SDSS) | no SED mass/SFR/$m_r$ exists |
| 20240210A | 9.67 CIGALE (Leung) / 9.71 Bagpipes (Li) | 0.81(13) Bagpipes (Li); 1.19 radio (Roxburgh) | 14.919 FORS2 $R$ | **filled** |
