# Audit group E – notes

Hosts: 20200430A, 20220105A, 20220717A, 20230708A, 20240310A, 20240201A.
All rows in `audit_group_E.csv` carry `verdict=new` (none of these hosts was in the round-1 table). 314 rows, 18 distinct sources. New BibTeX entries (not in `host_properties.bib`) are in `audit_group_E.bib`: `rajwade2024`, `dial2025`, `bera2025`, `chang2015`, `skrutskie2006`.

Conventions: `DM` in pc cm$^{-3}$; magnitudes as published (AB unless stated; Galactic-extinction correction state recorded in `notes`); `logMstar` in $\log(M_*/M_\odot)$; inclinations from axis ratios use $\cos^2 i = (q^2-q_0^2)/(1-q_0^2)$ with $q_0=0.2$ where the source paper says so.

## Papers used (primary sources)

| key | paper | status | hosts in group E |
|---|---|---|---|
| heintz2020 | Heintz et al. 2020, ApJ 903, 152 | refereed | 20200430A |
| bhandari2022 | Bhandari et al. 2022, AJ 163, 69 | refereed | 20200430A |
| mannings2021 | Mannings et al. 2021, ApJ 917, 75 | refereed | 20200430A |
| gordon2023 | Gordon et al. 2023, ApJ 954, 80 | refereed | 20200430A, 20220105A |
| gordon2025 | Gordon et al. 2025, ApJ 993, 119 (arXiv:2506.06453) | refereed | all CRAFT hosts |
| shannon2025 | Shannon et al. 2025, PASA 42, e036 (arXiv:2408.02083) | refereed | all CRAFT hosts |
| scott2025 | Scott et al. 2025, PASA 42, e133 (arXiv:2505.17497) | refereed | all CRAFT hosts |
| glowacki2025 | Glowacki et al. 2025, PASA 42, e157 (arXiv:2506.23403) | refereed | all CRAFT hosts |
| muller2026 | Muller et al. 2026, ApJ 1001, 118 (arXiv:2506.20774; FURBY I) | refereed | 20220105A, 20230708A, 20240201A |
| leung2025 | Leung et al. 2025 (arXiv:2507.16816) | arXiv | 20200430A, 20240310A, 20240201A |
| li2025 | Li et al. 2025 (arXiv:2507.01270) | arXiv | all except 20240201A offsets (see below) |
| athukoralalage2026 | Athukoralalage et al. 2026 (arXiv:2608.23868) | arXiv | 20200430A, 20230708A, 20240201A |
| connor2025 | Connor et al. 2025 (arXiv:2409.16952) | refereed | DM compilation only |
| roxburgh2025 | Roxburgh et al. 2025 (arXiv:2507.06594) | refereed/accepted | 20240201A (HI) |
| rajwade2024 | Rajwade et al. 2024, MNRAS 532, 3881 | refereed | 20220717A |
| dial2025 | Dial et al. 2025, MNRAS 536, 3220 | refereed | 20230708A |
| bera2025 | Bera et al. 2025, ApJ 982, 119 (arXiv:2411.14784) | refereed | 20230708A (polarisation only; no host numbers, not in CSV) |
| chang2015 | Chang et al. 2015, ApJS 219, 8 (SDSS+WISE catalogue) | refereed | 20240201A |
| skrutskie2006 | Skrutskie et al. 2006, AJ 131, 1163 (2MASS XSC) | refereed | 20240201A |

Prompt correction: FRB 20220717A is **not** a CRAFT burst. It was detected and localised by MeerTRAP (MeerKAT transient buffer, UHF band) and published by Rajwade et al. (2024). None of the CRAFT papers (Shannon, Gordon 2025, Scott, Glowacki, Muller) contain it. Its only follow-up appearance is the Li et al. (2025) Bagpipes compilation.

## Independent search log

Queries were run on 2026-10-07. Semantic Scholar returned HTTP 429 after the first few hits, so the sweep was completed with Crossref, the arXiv API (`export.arxiv.org/api/query`) and the arXiv web search, plus ADS abstract pages for specific bibcodes. Each host was searched in both TNS-style (`FRB 20240201A`) and short-form (`FRB 240201`) naming, and for 20200430A also as `FRB 200430` / `HG200430`.

| host | queries (all/abs) | new hits beyond the primary list |
|---|---|---|
| 20200430A | `"FRB 20200430A"`, `"FRB 200430"`, `"HG200430"`, `"FRB200430"` | Mannings et al. 2021 (offsets, HST not available for this host – ground-based entry), Connor et al. 2025 (DM compilation), Eftekhari et al. 2023 (Hα flux quoted by Glowacki, not downloaded – value recorded in notes column only) |
| 20220105A | `"FRB 20220105A"`, `"FRB 220105"` | none beyond primary list; Connor et al. 2025 compilation |
| 20220717A | `"FRB 20220717A"`, `"FRB 220717"`, `"MeerTRAP transient buffer FRB"` | Rajwade et al. 2024 (discovery / host paper); Li et al. 2025 only compilation |
| 20230708A | `"FRB 20230708A"`, `"FRB 230708"`, `"quasi-periodic FRB CRAFT"` | Dial et al. 2025 (burst paper, host z), Bera et al. 2025 (polarisation; no host properties) |
| 20240310A | `"FRB 20240310A"`, `"FRB 240310"` | none; not in FURBY I (Muller 2026); Leung et al. 2025 excludes it as edge-on |
| 20240201A | `"FRB 20240201A"`, `"FRB 240201"`, `"WISEA J095937.44+140519.4"`, `"SDSS J095937.45+140519.5"`, `"FRB host HI MeerKAT"`, `"FRB host WALLABY"` | Roxburgh et al. 2025 (MeerKAT HI, environment), Leung et al. 2025 (inclination remark). Catalogue look-ups: VizieR `J/ApJS/219/8` (Chang et al. 2015, NYU-VAGC 2491950) and 2MASS XSC via NED/IRSA (Skrutskie et al. 2006). No WALLABY detection (field is at Dec $+14^\circ$, outside the WALLABY footprint), no kinematic / IFU paper, no independent photometric decomposition paper found. |

No ADS token was available, so ADS was used only via public abstract pages for verification of bibcodes/volumes/pages (Rajwade 2024, Dial 2025, Bera 2025, Chang 2015).

## Per-host summaries

### FRB 20200430A (HG200430, PSO J151849.52+122235.82)

- **Discovery/host**: Heintz et al. 2020 (Sec. 2.2.5), $z=0.160$ (NOT/ALFOSC); Bhandari et al. 2022 and Gordon et al. 2023 give $z=0.1608$ (Keck/DEIMOS). Shannon et al. 2025 give the final CRAFT position and $\mathrm{DM}=380.1\pm0.2$, $\mathrm{DM_{MW}}=27$ (NE2001). Scott et al. 2025: $\tau_{\rm scatt}=6.5\pm0.15$ ms at 863.5 MHz, $\mathrm{RM}=195.3\pm0.7$ rad m$^{-2}$, $L/I=0.43$.
- **Photometry**: three generations of Pan-STARRS magnitudes exist. Heintz et al. 2020 used PSF magnitudes ($r=21.51$); Bhandari et al. 2022 explicitly replaced them by Kron magnitudes ($r=21.25$); Gordon et al. 2023 used DECaLS DECam $r=21.05\pm0.02$ (extinction-corrected) which is the value Glowacki et al. 2025 adopt. All three are recorded; the DECaLS value is recommended as `m_r` (deepest, extinction-corrected). VLT/FORS2 $g$ and $I$ (Shannon 2025 Table 5) are not extinction-corrected.
- **Mass/SFR**: CIGALE (Heintz 2020) $\log M_*=9.11$, SFR $<0.2$; CIGALE re-fit with Kron photometry (Bhandari 2022) $\log M_*=9.32$, Hα SFR $0.26\pm0.08$; Prospector (Gordon 2023) $\log M_*=9.30^{+0.07}_{-0.10}$, SFR$_{0-100\,\rm Myr}=0.11^{+0.06}_{-0.04}$, $t_m=5.99$ Gyr, $A_V=1.08$. Prospector is the target column. Leung et al. 2025 adopt $\log M_*=9.51$, which is the Gordon 2023 *total mass formed*, not the surviving stellar mass – noted in the CSV. Li et al. 2025 Bagpipes fit has quality flag 4 (bad) and is recorded only for completeness.
- **Size/offset**: $R_e=0.57''$ / 1.6 kpc (Heintz 2020, Keck/LRIS I), $R_e=1.59$ kpc, $n=1.22$ (Gordon 2025, FORS2 I); Glowacki 2025 Table 3 gives $R_e=1.99\pm0.01$ kpc from Marnoch et al. (in prep) GALFIT, 25 % larger. Offsets range from 1.7 kpc (Bhandari 2022, old FRB position) to $4.17^{+1.17}_{-1.28}$ kpc (Gordon 2025, Shannon 2025 position). Recommend Gordon 2025 (uses the final CRAFT position and its full error ellipse).
- **Inclination**: Glowacki et al. 2025 Table 3: $b/a=0.316\pm0.006$, $i=75.56\pm0.45^\circ$. Li et al. 2025 Table 2 quote $i=43.88^\circ$ from a Sérsic fit to survey imaging – **conflict**. Glowacki's value comes from deep VLT imaging; Li's from shallow DESI/Pan-STARRS images of a 21-mag galaxy and their SED flag for this host is "bad". Recommend Glowacki.
- **Type**: star-forming (Gordon 2023, Tacchella criterion; Glowacki "SF"), although Heintz 2020 placed it in the green valley on the basis of $u-r=1.78$ and the CIGALE SFR limit.
- Mannings et al. 2021 Table 3 include this host with ground-based (Pan-STARRS) offset $1.3\pm0.79''$, $2.97\pm2.36$ kpc; no HST data.

### FRB 20220105A

- Localised by CRAFT (RACS high-band epoch, Shannon 2025 Sec. 5.1; first published in Gordon 2023). $z=0.2785$ (Keck/DEIMOS, Gordon 2023; X-shooter 0.2784, Muller 2026). $\mathrm{DM}=583\pm2$, $\mathrm{DM_{MW}}=22$. Scott 2025: $\tau=0.43$ ms at 1649.8 MHz, $\mathrm{RM}=-1312\pm8$ rad m$^{-2}$ (one of the largest |RM| in the CRAFT sample), $L/I=0.30$.
- No catalogue name is given in any paper; recorded as "(uncatalogued) host at 13:55:12.92 +22:27:59.4".
- **Photometry**: Pan-STARRS $r=21.19\pm0.08$ (Gordon 2023, extinction-corrected; adopted by Glowacki). VLT/FORS2 $R=21.317$ (Shannon 2025, uncorrected) and $R=21.27$ (Muller 2026, corrected). **Anomaly**: VLT/HAWK-I $K_s$ differs by 0.7 mag between Gordon 2023 ($19.76\pm0.02$, extinction-corrected) and Shannon 2025 Table 5 ($19.06\pm0.08$); same instrument, difference unexplained (aperture or zero-point). Both recorded in notes; the Gordon value is in the CSV row.
- **Mass/SFR**: Prospector (Gordon 2023) $\log M_*=10.01^{+0.05}_{-0.07}$, SFR $=0.42^{+0.31}_{-0.19}$ (0–100 Myr), $t_m=5.67$ Gyr, $A_V=1.15$. Hα SFR (Muller 2026, X-shooter) $0.425\pm0.008$, $12+\log(\mathrm{O/H})=8.5$ (N2). **Compilation error**: Glowacki et al. 2025 Table 2 list $\log \mathrm{SFR}=-0.93\pm0.01$ for this host; $\log(0.42)=-0.38$. We treat the Glowacki entry as a transcription error and cite Gordon 2023.
- **Size/offset**: $R_e=0.79''$ / 3.37 kpc, $n=0.75$ ($R$); 3.16 kpc, $n=1.10$ ($K$) (Gordon 2025). Glowacki: $3.01\pm0.03$ kpc. Offset $2.29''$ / $9.8^{+4.0}_{-3.7}$ kpc = 2.9 $R_e$ (Gordon 2025); large offset drives the globular-cluster-model preference in Gordon 2025 Sec. 5.2. Glowacki 7.91 kpc, Li 8.45 kpc.
- **Inclination**: Glowacki 2025 Table 3 gives $b/a=0.259\pm0.001$ but leaves the inclination cell blank (verified in the PDF). Li et al. 2025: $i=56.26^\circ$. With $q_0=0.2$, $b/a=0.259$ implies $i\approx80^\circ$, so the two are inconsistent; only Li's number is a published inclination.
- PATH: Shannon 2025 Table 3 $P(O|x)=1.0$ with prior $P(O)=0.57$; Muller 2026 0.9998.

### FRB 20220717A (PSO J293.3038−19.2876) – MeerTRAP, not CRAFT

- Rajwade et al. 2024 is the only primary paper. MeerKAT UHF (816 MHz) transient-buffer localisation, position 19:33:13.0 −19:17:15.8 ($\pm0.9''$). $z=0.36295\pm0.00018$ (Keck/DEIMOS; Fig. 6 caption says 0.3633). $\mathrm{DM}=637.34\pm3.52$ (S/N-maximising; scattering-corrected 634.69), $\mathrm{DM_{MW}}=118$ (NE2001) / 83 (YMW16), $\mathrm{DM_{halo}}=86$ (also 52 quoted), $\mathrm{DM_{host}}\sim100$ (from Hα EM, order of magnitude). $\tau_{1\,\rm GHz}=8.2\pm0.3$ ms with index $-3.7$; $\mathrm{RM}=385.7\pm0.4$; $L/I=0.30$; fluence $6.83$ Jy ms.
- Host: Pan-STARRS DR1 galaxy, $m=21.8$ (PATH table, filter not stated – assumed $r$), PATH $P(O|x)=0.969$, separation $0.6''$ (≈3.0 kpc at $z=0.363$, derived here). Hα flux $1.708\pm0.37\times10^{-16}$ (on a skyline, underestimated) gives SFR $\geq0.65\pm0.14$ (Kennicutt 1994) – authors call it a lower limit. No stellar mass, size, inclination or metallicity published by Rajwade et al.
- The only stellar mass is Li et al. 2025 Bagpipes: $\log M_*=10.84\pm0.60$, SFR $3.23\pm11.08$ (quality flag 1; essentially unconstrained). Recorded with that caveat. Li et al. leave inclination and offset blank for this host.

### FRB 20230708A (J201227.73−552122.72) – dwarf host

- Dial et al. 2025 (burst paper; $z=0.105$ from host spectrum), Muller et al. 2026 (FURBY host paper), Shannon 2025, Scott 2025, Gordon 2025, Glowacki 2025. $\mathrm{DM}=411.51\pm0.05$, $\mathrm{DM_{MW}}=50$. $\mathrm{DM_{host}}\approx200$ (Muller 2026 text) / $228^{+38}_{-46}$ rest-frame (Athukoralalage 2026).
- **Fluence conflict**: Shannon 2025 Table 1 $111\pm4$ Jy ms (ICS) vs Dial 2025 Table 2 $89.2\pm0.9$ Jy ms (HTR, 10 % systematic). Both CRAFT papers; recorded both. **Scattering**: Dial $0.17\pm0.02$ ms, Scott $0.24\pm0.02$ ms at ≈920 MHz, Muller $0.14$; differences stem from the multi-component fits of a five-sub-burst quasi-periodic profile. RM: $-6.9\pm0.04$ (Dial) / $-7.5\pm0.4$ (Scott). Bera et al. 2025 analyse the intra-burst polarisation state changes (no host information).
- **Photometry**: VLT/FORS2 $R=22.65$ (Shannon, uncorrected) / $22.53$ (Muller, corrected). Shannon Table 5 HAWK-I $K_s=18.8\pm0.1$ is implausibly bright/red for a $\log M_*\approx8$ dwarf ($R-K_s\approx3.8$); flagged in notes as likely contaminated or erroneous.
- **Mass/SFR**: Prospector (Muller 2026, Kroupa IMF, continuity SFH, joint photometry + X-shooter spectrum) $\log M_*=7.97^{+0.09}_{-0.08}$, SFR$_{0-100\,\rm Myr}=0.04^{+0.02}_{-0.01}$, $\log \mathrm{sSFR}=-9.62$, $t_m=5.82$ Gyr; Hα SFR $0.0262\pm0.0009$ (slit). Metallicity limits only: $8.0\lesssim12+\log(\mathrm{O/H})\lesssim8.3$ ([N II] undetected). Li et al. Bagpipes $\log M_*=7.78\pm0.14$, consistent. Lowest-luminosity non-repeater host ($\approx0.01\,L^*$); no PRS.
- **Size/offset**: $R_e=0.31''$ / 0.59 kpc with $n$ fixed to 1 ($R$, Gordon 2025); $K$-band free fit $n=2.69$, 0.78 kpc. Glowacki $0.58\pm0.08$ kpc. Offset $0.57''$ / $1.1^{+0.72}_{-0.54}$ kpc = 1.9 $R_e$ (Gordon); Glowacki $0.33\pm0.63$ kpc. Shannon PATH prior $P(O)=0.087$ but posterior 1.0; Dial 2025 false-association probability $8\times10^{-7}$.
- **Inclination**: Glowacki 2025 Table 3 $b/a=0.44\pm0.06$, $i=66.19\pm4.5^\circ$. Li et al. leave inclination blank. For a $0.3''$ dwarf this is a marginally resolved measurement.

### FRB 20240310A

- Shannon 2025 (Sec. 5.2.15; detected during an EMU observation, arrival time known to ~15 s), $z=0.127$ (X-shooter emission lines), $\mathrm{DM}=601.8\pm0.2$, $\mathrm{DM_{MW}}=36$. Scott 2025: $\tau=2.23$ ms at 846 MHz, $\mathrm{RM}=-1709.2\pm1.1$ rad m$^{-2}$ (largest |RM| in the CRAFT sample), $L/I=0.72$ (Glowacki Table 1 quotes 0.75). Scott Sec. 4.2 remark that scattering is small for the large implied $\mathrm{DM_{host}}$ ($\approx 600-36-\mathrm{DM_{cosmic}}(0.127)\sim 400$).
- No catalogue name; host "(uncatalogued) at 01:10:29.25 −44:26:21.5". Second PATH candidate ($3.65''$, mag 21.8) has $P(O|x)=0.012$; primary $P(O|x)=0.9884$.
- **Photometry**: VLT/FORS2 $R=20.143\pm0.005$, HAWK-I $K_s=18.51\pm0.07$ (Shannon, uncorrected). No $r$-band survey magnitude published.
- **Mass/SFR conflict**: Gordon 2025 give $\log M_*=9.97$ from a $K$-band mass-to-light ratio (0.6 $M_\odot/L_\odot$) with no uncertainty, not an SED fit; Li et al. 2025 Bagpipes give $9.14\pm0.12$ (quality flag 0) – 0.8 dex apart. No Prospector or Hα SFR exists (not in FURBY I); the only SFR is Li's Bagpipes $0.38\pm0.18$. Recommend quoting Li with the caveat, and flagging the mass as uncertain pending FURBY II.
- **Size/offset**: $R_e=0.92''$ / 2.11 kpc, $n=0.93$ (Gordon 2025 $R$); Glowacki $3.21\pm0.01$ kpc (50 % larger). Offset $0.95''$ / $2.18^{+1.24}_{-1.07}$ kpc = 1.0 $R_e$ (Gordon); Glowacki $4.64\pm0.26$ kpc (twice Gordon). Gordon note that the optical and NIR model preference differ for this host.
- **Inclination**: Glowacki 2025 Table 3 $b/a=0.22\pm0.002$, $i=84.58\pm0.27^\circ$; Leung et al. 2025 call it $i\sim85^\circ$ and exclude it as edge-on. Li et al. 2025 $i=56.27^\circ$ – **conflict**; recommend Glowacki (deep VLT imaging).
- Type: Glowacki "SF (Y)" from the emission-line spectrum; no SED-based classification.

### FRB 20240201A (WISEA J095937.44+140519.4; SDSS/NYU-VAGC 2491950) – bright edge-on disk

- Shannon 2025 Sec. 5.2.11: filler observation at $|b|=50^\circ$, FRB landed in a null of the primary beam, detected in an edge beam; fluence (47 Jy ms) and the polarisation fractions are flagged uncertain (Scott 2025 flag c, excluded from Glowacki's polarisation analysis). $\mathrm{DM}=374.5$, $\mathrm{DM_{MW}}=38$; Scott: $\tau=0.78$ ms at 915 MHz, $\mathrm{RM}=1275.0$ rad m$^{-2}$; Athukoralalage: $\mathrm{RM_{host}}=1380$, $\mathrm{DM_{host}}=270$ (rest frame).
- **Redshift typo**: Shannon 2025 Sec. 5.2.11 text gives $z=0.047279$; Table 1, Glowacki, Roxburgh and the SDSS/NYU-VAGC catalogue (Chang 2015) give $z=0.042729$. The text value is a digit transposition. CSV cites 0.042729 (Shannon Table 1) and the catalogue value 0.0427292 (Chang 2015).
- **Photometry**: VLT/FORS2 $R=16.97$ (Shannon, uncorrected), $16.91$ (Muller 2026, corrected; adopted by Glowacki). HAWK-I $K_s=15.41$. No SDSS model magnitude row was added because the CRAFT papers do not quote one and the galaxy is saturated/shredded in SDSS pipeline photometry; the Shannon PATH table "mag 26.17" for the automated centroid is a segmentation artefact (noted).
- **Mass/SFR – compilation error found**: Glowacki et al. 2025 Table 2 quote $\log M_*=10.21\pm0.02$, $\log \mathrm{SFR}=0.14\pm0.02$, $\log \mathrm{sSFR}=-10.07$ "from Chang et al. 2015". Querying VizieR `J/ApJS/219/8` for NYU-VAGC 2491950 shows these are the **97.5th-percentile columns** (`lmass97.5`, `lsfr97.5`). The medians are $\log M_*=10.04^{+0.08}_{-0.10}$ and $\log \mathrm{SFR}=-0.208$ (SFR $=0.62\,M_\odot$ yr$^{-1}$, 16–84 %: 0.42–0.93), $\log \mathrm{sSFR}=-10.22$. The CSV cites Chang 2015 with the median values; Athukoralalage 2026 inherit the 10.21 value. Gordon 2025 give $\log M_*=10.22$ from the $K_s$ mass-to-light method (no error) – coincidentally close to the Glowacki number but independent. Li et al. Bagpipes: $9.99\pm0.12$, SFR $1.18$. Muller 2026 Hα (X-shooter slit) SFR $0.106$ is aperture-limited (slit covers a small fraction of a $4.4''$-$R_e$ edge-on disk) and is a lower limit; $12+\log(\mathrm{O/H})=8.6$ (N2). Roxburgh 2025 radio-continuum SFR $1.04$ is for the interacting pair combined.
- **Size/offset – inconsistent literature**: Gordon 2025 $R_e=4.38''$ / 3.73 kpc, $n=0.68$ ($R$); 2.67 kpc, $n=1.23$ ($K$). Glowacki 4.13 kpc. Offsets: Gordon $3.3''$ / $2.81^{+0.44}_{-0.42}$ kpc (0.75 $R_e$); Glowacki $9.63$ kpc (3.4× Gordon, 2.33 $R_e$); Li $4.38$ kpc (identical to their 20240210A row – likely copy error); Shannon PATH table $6.25''$ (=5.3 kpc) to the automated (shredded) centroid. Recommend Gordon 2025. For an edge-on disk the projected offset along the minor axis is the physically relevant quantity; none of the papers decomposes it.
- **Inclination / axis ratio**: Glowacki 2025 Table 3 gives $b/a=0.166$ with the inclination cell blank, because $b/a<q_0=0.2$ makes the formula undefined (edge-on). Leung et al. 2025 state $i\approx80^\circ$ and exclude it ($>70^\circ$). Li et al. 2025 $i=82.45^\circ$. 2MASS XSC (Skrutskie 2006) gives a $K_s$ isophotal $b/a=0.44$, PA $=60^\circ$, semi-major $7.3''$ – NIR isophotes are much rounder than the optical 0.166, as expected for a bulge-dominated NIR light distribution plus dust lane in the optical. No independently published photometric or kinematic inclination exists beyond these.
- **HI / environment (Roxburgh 2025, MeerKAT)**: tidally interacting pair; neighbour at $z=0.04314$ ($\Delta v\approx300$ km s$^{-1}$); shared HI envelope, classified HI-disturbed. Pair-combined $S_{\rm HI}=1.01$ Jy km s$^{-1}$, $W_{50}=373\pm5$ km s$^{-1}$, $\log M_{\rm HI}=9.96$, 1.28 GHz continuum 0.42 mJy. Per-galaxy HI mass and HI kinematic inclination are not given. Not in the WALLABY footprint.
- Type: star-forming spiral (Chang 2015 Sp flag; Glowacki "SF (Y)"); Muller 2026 report no strong Balmer absorption.

## Cross-cutting anomalies

1. Glowacki et al. 2025 Table 2/3 contain three transcription issues affecting this group: 20240201A mass/SFR = Chang 2015 97.5th percentiles; 20220105A $\log \mathrm{SFR}=-0.93$ vs Gordon 2023 $-0.38$; $R_e$ and offsets from Marnoch et al. (in prep) differ by 25–240 % from Gordon et al. 2025 for the same hosts and same VLT imaging. CSV cites originals and records the Glowacki values as separate rows for traceability.
2. Li et al. 2025 inclinations (43.9°, 56.3°, 56.3°, 82.5°) disagree with the Glowacki/VLT axis ratios for 20200430A and 20240310A by 20–30°. Li's values come from shallow survey imaging; their SED quality flags are 4 (bad) for 20200430A. Recommend Glowacki for all CRAFT hosts.
3. Leung et al. 2025 use the Gordon 2023 *total mass formed* (9.51) rather than surviving mass (9.30) for 20200430A.
4. Shannon et al. 2025 Table 5 $K_s$ magnitudes are discrepant with other papers for 20220105A (0.7 mag) and implausible for 20230708A.
5. Shannon et al. 2025 Sec. 5.2.11 redshift typo (0.047279 → 0.042729) for 20240201A.
6. Athukoralalage et al. 2026 (arXiv only) host-DM/RM values depend on the adopted masses above; included as `DM_host`/`RM_host` rows with the caveat in notes.
