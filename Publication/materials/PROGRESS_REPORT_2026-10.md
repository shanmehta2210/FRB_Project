# FRB host inclinations: progress report (October 2026)

## Summary

- **Sample.** 56 confirmed, literature-localized, non-repeating FRB hosts with $m_r \le 22$ and $q > 0.2$. All 56 have a published DM. This replaces the 51-host sample in the committed paper results.
- **Main test.** The host inclinations are consistent with HSC-PDR2 field discs (KS $p = 0.45$, matched-$N$ $p = 0.44$; AD matched-$N$ $p = 0.55$). Hosts are, if anything, slightly *more face-on* than the field. The edge-on deficit claimed by Bhardwaj et al. (2024) is not reproduced: our edge-on count has $P(K \ge k) = 0.002$ under their 4/23 edge-on fraction.
- **DM_host.** No correlation with inclination: $r_S = -0.065$, $p = 0.64$, $N = 55$, 95% CI $[-0.32, +0.20]$. This excludes the $r_S \approx -0.44$ of Marnoch et al. (2026), which we would detect with about 91% probability.
- **Marnoch re-check.** Their numbers reproduce from their Table 1. Removing their four $q < 0.2$ hosts (set to $i = 90^\circ$) weakens the result to $p \approx 0.05$, and it fails in every subsample. On the 20 hosts both samples share, their own data give $r_S = -0.19$ ($p = 0.41$).
- **Other properties.** 41 host and burst properties were tested against inclination. None is significant after a false-discovery correction. The one physically expected trend, old-star (diffuse) dust attenuation rising towards edge-on ($r_S = -0.46$, $p = 0.010$), does appear.

## 1. Sample curation

### 1.1 Policies (set in this round)

| Policy | Reason | Effect |
|---|---|---|
| No repeaters | Repeaters may be a distinct population with different environments | 20220912A excluded (only repeater in the confirmed set) |
| Literature-claimed hosts only | A published association uses PATH plus redshift and follow-up; our own PATH alone is not an association | 9 Dial+2026 CRAFT HTR2 hosts (PATH-only, host "in prep") and 20250520A excluded |
| Published DM required | DM_host tests need it for every host | 20240615B excluded (Gordon+2025 host; radio properties not yet published) |

### 1.2 Sample history

| Step | N |
|---|---|
| Committed paper sample (included repeater 20220912A) | 51 |
| Earlier confirmations (20240203D, 20240224A, 20240312D, 20240525A, 20240615B) | 56 |
| Repeater 20220912A removed | 55 |
| October 2026 literature batch, 14 candidates, 11 pass triage | 66 |
| 9 PATH-only Dial+2026 hosts removed | 57 |
| 20240615B removed (no published DM) | **56** |

Literature batch: 14 candidates from Dial+2026, Pastor-Marazuela+2026, CHIME/FRB 2025, and DSA-110. Each got a production GALFIT fit and the triage protocol. Kept: **20231210F** (MeerTRAP; Pastor-Marazuela+2026 Sec. 4.1.11 name G1 the putative host, PATH 0.81) and **20250316A** (CHIME Outriggers; NGC 4141; CHIME/FRB 2025, ApJL 989, L48). Rejected: 20231010A (fails cuts), 20220208A (host undetected in PS1), 20250520A (no LS/DES imaging).

### 1.3 Final sample (N = 56)

| Property | Value |
|---|---|
| Surveys | ASKAP/CRAFT 25, DSA-110 25, MeerKAT 5, CHIME 1 |
| Redshift | 0.0067 to 0.66, median 0.238 |
| Fit leg used | production 42, $n=1$ + fixed sky 7, $n=1$ 5, fixed sky 1, PSF 1 |
| Cuts | confirmed host, $m_r \le 22$, $q > 0.2$, non-repeater, literature host, published DM |

Inclination: $\cos^2 i = (q^2 - q_0^2)/(1 - q_0^2)$ with $q_0 = 0.2$. Uncertainty (Protocol B): $\sigma_q$ is the GALFIT error combined in quadrature with the sky-variation term, plus a $5^\circ$ floor on $i$.

## 2. Dispersion measures

### 2.1 Coverage

- 56/56 hosts have a published DM.
- 17 rows of `master_frb_localization.csv` were missing DM, $z$ or DM$_{\rm MW}$ and were filled from the cited literature table.
- New-host DMs:

| FRB | DM (pc cm$^{-3}$) | $z$ | Source |
|---|---|---|---|
| 20231210F | 720.6 | photometric only | Pastor-Marazuela+2026 Table 2 |
| 20240203D | 272.60 | 0.0740 | Verdi+2026 Table A1 |
| 20240224A | 881.20 | 0.3726 | Verdi+2026 Table A1 |
| 20240312D | 334 (approx.) | 0.04986 | Jahns-Schindler+2026 |
| 20240525A | 491.23 | 0.3266 | Dial+2026 Table B1; Gordon+2025 |
| 20250316A | 161.82 | 0.0067 | Athukoralalage+2026; CHIME/FRB 2025 |

20240312D: the fitted DM row did not survive text extraction. The paper's DM budget sums to about 337, so the 334 dedispersion DM is adequate for rank tests.

### 2.2 DM_host method

$$
\mathrm{DM_{host,obs}} = \mathrm{DM_{obs}} - \mathrm{DM_{MW,NE2001}} - \widetilde{\mathrm{DM}}_{\rm cosmic}(z), \qquad
\mathrm{DM_{host,rest}} = \mathrm{DM_{host,obs}}\,(1+z)
$$

- **DM$_{\rm obs}$:** from the localization papers, chosen with a fixed source-priority list.
- **DM$_{\rm MW}$:** NE2001 (Cordes & Lazio 2002), recomputed for every sightline with `mwprop`. It reproduces the tabulated values to within 0.01 pc cm$^{-3}$ and caught misprints in Shannon+2025.
- **$\widetilde{\mathrm{DM}}_{\rm cosmic}$:** the *median* of the Macquart et al. (2020) distribution ($F = 0.2$, Planck 2018). The distribution is skewed, so the mean overestimates a typical sightline.
- **Halo:** no Milky Way halo term (following Bhardwaj+2024). It adds a near-constant 30–50 pc cm$^{-3}$ and cannot change rank correlations.
- **Redshift:** spectroscopic only, so 20231210F has no DM_host. **N = 55.**
- **Cross-check:** published DM_host values for 43 hosts, which mix methods.

Caveat: the scatter in DM$_{\rm cosmic}$ at fixed $z$ is of order 100 pc cm$^{-3}$ at $z \sim 0.3$. This noise affects any DM_host analysis, Marnoch's included. Rank tests are less sensitive to it.

## 3. Paper statistical tests (51 to 56 hosts)

Control: HSC-PDR2 (Kawinwanichakij+2021), $0.4 < n < 1.5$, `goodfits_flag = 1`, $r \le 22$, $q > 0.2$; $N = 24{,}450$. Matched-$N$ $p$ is the fraction of 10,000 HSC draws of the same size whose statistic against the pool is at least the FRB value.

| Test | 51 hosts | 56 hosts |
|---|---|---|
| Median $\cos i$, FRB / HSC | 0.610 / 0.575 | 0.626 / 0.575 |
| KS $D$; $p$ vs pool (matched-$N$) | 0.115; 0.476 (0.467) | 0.112; 0.446 (0.438) |
| Anderson–Darling matched-$N$ $p$ | 0.606 | 0.545 |
| Cramér–von Mises matched-$N$ $p$ | 0.476 | 0.419 |
| Kuiper matched-$N$ $p$ | 0.443 | 0.385 |
| Wasserstein-1 matched-$N$ $p$ | 0.522 | 0.473 |
| Mann–Whitney $p$ | 0.543 | 0.453 |
| One-sided KS / MWU $p$ (FRB more face-on) | 0.241 / 0.271 | 0.226 / 0.226 |
| $\chi^2$, HSC quartiles: counts edge-on to face-on; matched-$N$ $p$ | 11, 10, 19, 11; 0.254 | 13, 10, 19, 14; 0.404 |
| Edge-on $i > 60^\circ$: $k/N$, binomial $p$ | 15/51, 0.083 | 17/56, 0.094 |

Protocol B error Monte Carlo (2,000 draws):

| Quantity | 51 hosts | 56 hosts |
|---|---|---|
| Median $\cos i$ (16/50/84%) | 0.594 / 0.618 / 0.642 | 0.608 / 0.631 / 0.651 |
| AD $p$, median | 0.25 (SciPy cap) | 0.25 (cap) |
| Fraction of draws with KS $p < 0.05$ | 0.039 | 0.068 |
| Edge-on count $k$, median | 17 | 19 |
| $P(K \ge k)$ under the 2024 fraction 4/23, median | 0.004 | 0.002 |

Appendix controls, 56 hosts:

| Control | KS $p$ (2-sided) | KS $p$ (one-sided) | KS $p$ (matched-$N$) |
|---|---|---|---|
| $q_0 = 0.10$ ($N = 57$) | 0.437 | 0.221 | 0.430 |
| $q_0 = 0.15$ ($N = 57$) | 0.461 | 0.234 | 0.455 |
| $q_0 = 0.20$ ($N = 56$) | 0.446 | 0.226 | 0.437 |
| $q_0 = 0.25$ ($N = 52$) | 0.208 | 0.105 | 0.204 |
| HSC $n$ window 0.5–1.2 / 0.8–1.2 / 0.4–2.0 | matched-$N$ $p$ 0.33 / 0.43 / 0.46 | | |

Bug fixed: the error Monte Carlo hard-coded $N = 50$ in its binomial tests (the sample was 51). It now uses the actual sample size.

## 4. Host-property correlations with inclination

### 4.1 Literature survey

| Paper | Relation tested | Their result |
|---|---|---|
| Marnoch+2026 (CRAFT, N = 34) | DM_host vs $q$, $i$; $m < 21$, $n < 2$ subsamples; DM_host vs offset, colour | $r_S = -0.44$, $p = 0.01$ |
| Bhardwaj+2024 | DM_host vs $i$ | none |
| Glowacki+2025 | $\tau$ vs inclination, $\Sigma_\star$, age, metallicity | $\tau$–$\Sigma_\star$, $\tau$–age; none with $i$ |
| Leung+2025 | DM_host vs $M_\star$ | anticorrelation |
| Li+2025 | DM$_{\rm exc}$ vs sSFR, SFR, $M_\star$, $i$ | tight with sSFR |
| Verdi+2026 | RM$_{\rm host}$, $B_\parallel$, DM_host vs $z$ | DM–RM correlated; no $z$ evolution |
| Day+2020; Cordes+2022 | DM vs $\tau$ | — |
| Kovacs+2024; Orr+2024 (sims) | DM_host rises with $i$ | predicted |

### 4.2 Method

- Spearman, Kendall and Pearson against $\cos i$; a negative $r_S$ means larger values towards edge-on.
- Protocol B Monte Carlo (2,000 draws) on $r_S$.
- Benjamini–Hochberg false-discovery $q$ across all 41 properties.
- Edge-on ($i > 60^\circ$) vs face-on Mann–Whitney split.
- Literature values chosen by the host-table priority rules; masses fall back to K-band / compiled values where no SED mass exists.
- Polarization put on one scale (Sherman+2024 quote percent).

### 4.3 Results (Spearman $r_S$, $p$)

| Property | N | $r_S$ ($p$) | Note |
|---|---|---|---|
| DM_host, rest frame | 55 | −0.065 (0.64) | main test |
| DM_host, observer frame | 55 | −0.059 (0.67) | |
| DM_host, published | 43 | −0.131 (0.40) | |
| DM $-$ DM$_{\rm MW}$ | 56 | −0.051 (0.71) | |
| DM$_{\rm obs}$ | 56 | −0.060 (0.66) | |
| log $\tau$ (scattering) | 29 | −0.167 (0.39) | Glowacki: none |
| log $\tau$, rest 1 GHz (published) | 15 | +0.223 (0.42) | |
| log $\lvert$RM$\rvert$ | 41 | −0.018 (0.91) | |
| $\lvert$RM$_{\rm host}\rvert$, rest (published) | 21 | +0.058 (0.80) | |
| $B_\parallel$, host (published) | 17 | −0.321 (0.21) | signed: direction, not strength |
| $B_\parallel$ proxy $1.23\,\lvert{\rm RM}\rvert(1+z)^2/{\rm DM_{host}}$ | 39 | −0.079 (0.63) | RM$_{\rm MW}$ not removed |
| Linear polarization fraction | 24 | −0.331 (0.11) | direction opposite to depolarization |
| Circular polarization fraction | 17 | −0.189 (0.47) | |
| Projected offset (kpc) | 46 | −0.075 (0.62) | |
| Offset / $R_e$ | 44 | +0.163 (0.29) | |
| log $M_\star$ | 52 | −0.113 (0.43) | |
| SFR | 46 | +0.123 (0.42) | |
| log sSFR | 46 | +0.138 (0.36) | |
| log $\Sigma_\star$ | 40 | +0.091 (0.58) | |
| Mass-weighted age | 30 | −0.183 (0.33) | |
| Stellar metallicity | 31 | +0.138 (0.46) | |
| Gas metallicity | 11 | −0.055 (0.87) | |
| H$\alpha$ EW | 6 | −0.029 (0.96) | |
| **$A_V$, old stars (diffuse dust)** | 30 | **−0.460 (0.010)** | expected edge-on reddening; BH $q = 0.43$ |
| $A_V$, young stars (birth clouds) | 31 | −0.039 (0.84) | no geometric dependence expected |
| $g - r$ | 18 | −0.211 (0.40) | |
| $u - r$ | 25 | −0.326 (0.11) | |
| PAH fraction | 16 | +0.356 (0.18) | |
| $R_e$ (kpc) | 41 | −0.333 (0.033) | BH $q = 0.46$; projection effect |
| Sérsic $n$ (ours) | 56 | +0.302 (0.024) | BH $q = 0.46$; edge-on discs fit lower $n$ |
| Redshift | 55 | +0.025 (0.86) | |
| $m_r$ | 56 | −0.170 (0.21) | |
| PATH posterior | 49 | +0.057 (0.70) | |
| S/N, fluence, width | 33–35 | +0.13 to +0.19 (> 0.26) | |
| Null checks: DM$_{\rm MW}$, $\lvert b\rvert$, RA, Dec, scintillation bandwidth | 7–56 | all $p > 0.25$ | |
| By survey (Kruskal–Wallis, ASKAP / DSA-110 / MeerKAT) | 55 | $H = 0.16$, $p = 0.92$ | |

Old-star $A_V$ medians: 0.48 edge-on vs 0.13 face-on. 89% of Monte Carlo draws keep $p < 0.05$. Pearson is weaker ($r = -0.24$, $p = 0.20$), so treat this as consistent with dust geometry, not a detection.

### 4.4 DM_host: Marnoch-style subsamples (rest frame, vs $q$)

| Subsample | N | Spearman ($p$) | Partial given $z$ ($p$) |
|---|---|---|---|
| All | 55 | −0.065 (0.64) | −0.067 (0.63) |
| $m_r < 21$ | 45 | −0.093 (0.54) | −0.096 (0.54) |
| $n < 2$ | 46 | −0.037 (0.81) | −0.039 (0.80) |
| $m_r < 21$ and $n < 2$ | 36 | −0.071 (0.68) | −0.073 (0.68) |
| $z \le 0.3$ | 37 | −0.065 (0.70) | −0.030 (0.86) |

Interpreting the null result (rest frame, $N = 55$):

- **Shuffle test.** Randomly re-pairing DM_host with inclination 20,000 times gives $|r_S| \ge 0.065$ in 64% of shuffles. The 95% chance range is $\pm 0.27$.
- **Confidence interval.** 95% CI $[-0.32, +0.20]$ (Fisher and bootstrap).
- **Power.** A true $r_S = -0.44$ would give $p < 0.05$ in about 90% of samples of 55.
- Weak correlations ($|r_S| \lesssim 0.2$) are not excluded.

### 4.5 Relations not involving inclination

| Relation | Claim | N | $r_S$ ($p$) | Verdict |
|---|---|---|---|---|
| DM_host vs log $\tau$ | Cordes+22, Day+20 | 29 | +0.450 (0.014) | reproduced |
| DM_host vs log $\lvert$RM$\rvert$ | Verdi+26 | 41 | +0.490 (0.001) | reproduced |
| DM_host vs log $M_\star$ | Leung+25 anticorrelation | 51 | −0.007 (0.96) | not reproduced |
| DM$_{\rm exc}$ vs log sSFR | Li+25 | 46 | +0.257 (0.084) | weak |
| DM$_{\rm exc}$ vs SFR | Li+25 | 46 | +0.280 (0.060) | weak |
| log $\tau$ vs log $\Sigma_\star$ | Glowacki+25 | 20 | +0.261 (0.27) | not reproduced (small N) |
| log $\tau$ vs age / metallicity | Glowacki+25 | 13–14 | $\le 0.09$ (> 0.78) | not reproduced |
| DM_host vs offset / $R_e$ | Marnoch+26 | 44 | +0.234 (0.13) | not reproduced |
| DM_host vs $z$ | Verdi+26 (none) | 55 | +0.082 (0.55) | agrees |

## 5. Re-check of Marnoch et al. (2026)

Inputs: per-host $z$, DM_host (their observer-frame Macquart proxy), $q$, $i$, $n$, $m$ from their Table 1 (N = 34). Rest frame = observer × $(1+z)$.

| Subsample (rest frame vs $q$) | As published | Without 4 hosts $q < 0.2$ | Also without repeater 20190711A |
|---|---|---|---|
| Full | −0.43 (0.011), N = 34 | −0.36 (0.052), N = 30 | −0.40 (0.032), N = 29 |
| $m \le 21$ | −0.55 (0.019), N = 18 | −0.47 (0.068), N = 16 | −0.47 (0.068), N = 16 |
| $n < 2$ | −0.39 (0.038), N = 28 | −0.29 (0.170), N = 24 | −0.33 (0.119), N = 23 |
| No bulge fits | −0.42 (0.018), N = 31 | −0.36 (0.065), N = 27 | −0.40 (0.045), N = 26 |

- **Reproduction.** Our values match their published $r_P$ and $r_S$ to within about 0.01.
- **Low-$q$ hosts.** They set 20220105A ($q = 0.026$), 20231230D (0.101), 20240201A (0.167) and 20240203C (0.166) to $i = 90^\circ$. All four have DM_host,rest of 166–358 pc cm$^{-3}$, so they are high-leverage points.
- **Shared hosts (N = 20).** $q$ agrees ($r_S = 0.86$) and DM_host agrees ($r_S = 0.99$; ours about 60 higher, no halo term). On these hosts there is no correlation in any combination:
  - their $q$, their DM: −0.19 (0.41)
  - their $q$, our DM: −0.25 (0.29)
  - our $q$, their DM: −0.09 (0.71)
  - our $q$, our DM: −0.16 (0.50)
- **Conclusion.** The disagreement comes from sample composition, not methodology. Their signal is carried by the 14 hosts outside our sample:
  - the four with $q < 0.2$ and the repeater;
  - 20200430A and 20240310A, which we also measure at $q < 0.2$;
  - very faint hosts such as 20190611B ($m = 24$, DM_host = −110);
  - 20240615B (no published DM).

## 6. Outstanding work

- Manuscript text still quotes 50/51 hosts and the old sample counts (64/60/53/50).
- Paper figures (CDF envelopes, $q_0$, $n$-window panels) not regenerated for N = 56.
- `build_host_property_table.py` still seeded with the old 51 hosts; host table and `host_table_57.tex` need refreshing.
- 20240615B can return once its DM is published.

## 7. Files

| Content | Path |
|---|---|
| Sample (N = 56) | `plots/plots_null/v2/frb_vs_hsc_confirmed50/confirmed50_q.csv` |
| Paper tests | `plots/plots_null/v2/frb_vs_hsc_confirmed50/TEST_RESULTS.md` (`scripts/run_frb_hsc_confirmed50_tests.py`) |
| Error Monte Carlo | `scripts/run_frb_hsc_error_mc.py` |
| Appendix controls | `Publication/materials/appendix/CONTROLS_Q0_NRANGE.md` |
| Correlations | `plots/plots_null/v2/host_correlations/CORRELATIONS.md` (`scripts/run_host_property_correlations.py`) |
| Marnoch re-check | `Publication/materials/marnoch_recheck/MARNOCH_RECHECK.md` (`scripts/recheck_marnoch2026.py`) |
| New-host literature values | `Publication/materials/host_properties/group_6.csv` |
| Host decisions | `pipeline_scripts/verification/host_confirmation.csv`, `HOST_TRIAGE_CASES.md` |
| Sample changelog | `pipeline_scripts/docs/WEAK_ASSOCIATIONS_PRODUCTION67.md` |
