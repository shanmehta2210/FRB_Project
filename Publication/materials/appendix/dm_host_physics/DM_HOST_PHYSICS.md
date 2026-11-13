# What sets DM_host?

56 hosts; 55 with DM_host. DM_host = (DM - DM_MW - median DM_cosmic)(1 + z). Sigma = quantity / (2 pi R_e^2), with the literature R_e where published and our GALFIT R_e otherwise; offset / R_e uses the literature offset; Delta MS is log SFR minus the Speagle+14 main sequence at the host redshift. BH q over the 34 hypothesis tests (sanity checks excluded), separately for Spearman and AD. AD: two-sample Anderson-Darling p that y has the same distribution below and above the median of x (permutation p; inside the DM_cosmic MC, the fraction of draws with asymptotic AD p < 0.05).

## Redshift selection

| x | y | sample | claim | N | r_S (p) | AD p (BH q) | 95% CI | tau (p) | JK1 | JK3: p<0.05 (max p) | most influential (r_S without) | partial tau given z (p) | DM_cosmic MC: r_S, p<0.05 (AD) | foreground-corrected | BH q |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| z | log M* | All | flux limit: more massive hosts at higher z | 51 | +0.32 (0.024) | 0.022 (0.19) | [+0.05, +0.55] | +0.22 (0.023) | +0.27 to +0.35 | 0.83 (0.22) | 20190523A (+0.27) |  |  |  | 0.12 |
| z | log SFR | All | flux limit: more star-forming hosts at higher z | 46 | +0.34 (0.020) | 0.090 (0.23) | [+0.05, +0.58] | +0.22 (0.030) | +0.30 to +0.41 | 0.87 (0.15) | 20211127I (+0.41) |  |  |  | 0.12 |
| z | log sSFR | All | cosmic sSFR rises with z | 46 | +0.14 (0.35) | 0.96 (1.00) | [-0.14, +0.40] | +0.08 (0.45) | +0.10 to +0.20 | 0.00 (0.95) | 20211127I (+0.20) |  |  |  | 0.54 |
| z | R_e | All | larger hosts at higher z | 55 | +0.05 (0.73) | 0.99 (1.00) | [-0.27, +0.35] | +0.03 (0.75) | +0.01 to +0.09 | 0.00 (1.00) | 20231020B (+0.09) |  |  |  | 0.80 |
| z | offset (kpc) | All | localization in kpc degrades with z | 45 | +0.30 (0.042) | 0.036 (0.19) | [-0.02, +0.59] | +0.22 (0.037) | +0.26 to +0.37 | 0.52 (0.34) | 20231020B (+0.37) |  |  |  | 0.16 |
| z | DM_host (observer frame) | All | observer-frame DM_host diluted by 1/(1+z) | 55 | -0.10 (0.47) | 0.37 (0.59) | [-0.38, +0.19] | -0.07 (0.46) | -0.16 to -0.05 | 0.00 (1.00) | 20220319D (-0.16) |  |  |  | 0.67 |
| z | DM_host | All | Verdi+26: no evolution | 55 | +0.08 (0.55) | 0.92 (1.00) | [-0.21, +0.36] | +0.05 (0.56) | +0.03 to +0.14 | 0.00 (1.00) | 20220418A (+0.14) |  | +0.03 [-0.05, +0.10], 0.00 (0.01) | -0.14 (0.33) | 0.73 |

## Stars and star formation

| x | y | sample | claim | N | r_S (p) | AD p (BH q) | 95% CI | tau (p) | JK1 | JK3: p<0.05 (max p) | most influential (r_S without) | partial tau given z (p) | DM_cosmic MC: r_S, p<0.05 (AD) | foreground-corrected | BH q |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| log M* | DM_host | All | Leung+25: anticorrelation | 51 | -0.01 (0.96) | 0.91 (1.00) | [-0.29, +0.28] | -0.00 (0.98) | -0.06 to +0.04 | 0.00 (1.00) | 20230907D (-0.06) | -0.01 (0.91) | -0.03 [-0.12, +0.04], 0.00 (0.00) | -0.02 (0.89) | 0.97 |
| log SFR | DM_host | All | Li+25: positive | 46 | +0.24 (0.10) | 0.021 (0.19) | [-0.08, +0.53] | +0.18 (0.077) | +0.20 to +0.31 | 0.09 (0.51) | 20210410D (+0.31) | +0.19 (0.067) | +0.19 [+0.10, +0.26], 0.07 (0.23) | +0.12 (0.46) | 0.24 |
| log SFR | DM_host | z < 0.2 |  | 20 | +0.01 (0.97) | 0.50 (0.75) | [-0.45, +0.48] | +0.00 (1.00) | -0.11 to +0.15 | 0.00 (1.00) | 20210410D (+0.15) | -0.01 (0.95) | -0.06 [-0.15, +0.03], 0.00 (0.00) | -0.17 (0.49) | 0.97 |
| log SFR | DM_host | z >= 0.2 | earlier lead: r_S = +0.55 | 26 | +0.55 (0.004) | 0.050 (0.19) | [+0.16, +0.80] | +0.40 (0.003) | +0.50 to +0.67 | 0.99 (0.079) | 20230307A (+0.67) | +0.40 (0.005) | +0.46 [+0.33, +0.57], 0.71 (0.34) | +0.39 (0.063) | 0.07 |
| log SFR | DM_host | Prospector SFRs | one SFR method | 30 | +0.11 (0.56) | 0.13 (0.26) | [-0.25, +0.48] | +0.10 (0.44) | +0.06 to +0.22 | 0.00 (1.00) | 20210410D (+0.22) | +0.11 (0.35) | +0.07 [-0.03, +0.16], 0.00 (0.02) | -0.00 (0.98) | 0.73 |
| log SFR | DM_host | Prospector, z >= 0.2 | one SFR method | 17 | +0.48 (0.051) | 0.032 (0.19) | [-0.02, +0.83] | +0.35 (0.052) | +0.39 to +0.68 | 0.21 (0.42) | 20230307A (+0.68) | +0.32 (0.082) | +0.41 [+0.24, +0.55], 0.29 (0.33) | +0.31 (0.27) | 0.17 |
| log sSFR | DM_host | All | Li+25: tight positive | 46 | +0.28 (0.064) | 0.063 (0.21) | [-0.02, +0.53] | +0.18 (0.074) | +0.23 to +0.33 | 0.24 (0.29) | 20210410D (+0.33) | +0.18 (0.048) | +0.23 [+0.15, +0.30], 0.21 (0.23) | +0.18 (0.25) | 0.18 |
| log sSFR | DM_host | z >= 0.2 |  | 26 | +0.35 (0.084) | 0.049 (0.19) | [-0.07, +0.67] | +0.23 (0.10) | +0.28 to +0.47 | 0.15 (0.47) | 20230307A (+0.47) | +0.23 (0.052) | +0.30 [+0.18, +0.41], 0.22 (0.35) | +0.26 (0.21) | 0.22 |
| log Sigma_SFR | DM_host | All | H-alpha emission measure: DM ~ Sigma_SFR^0.5 | 46 | +0.34 (0.022) | 0.085 (0.23) | [+0.02, +0.61] | +0.25 (0.014) | +0.29 to +0.40 | 0.81 (0.21) | 20210410D (+0.40) | +0.26 (0.023) | +0.28 [+0.20, +0.36], 0.47 (0.19) | +0.35 (0.023) | 0.12 |
| log Sigma_SFR | DM_host | z >= 0.2 |  | 26 | +0.50 (0.009) | 0.076 (0.23) | [+0.09, +0.79] | +0.37 (0.008) | +0.44 to +0.59 | 0.89 (0.17) | 20230307A (+0.59) | +0.37 (0.012) | +0.45 [+0.31, +0.55], 0.67 (0.26) | +0.54 (0.006) | 0.11 |
| log Sigma_SFR | DM_host | Prospector SFRs | one SFR method | 30 | +0.28 (0.13) | 0.33 (0.59) | [-0.10, +0.61] | +0.21 (0.11) | +0.22 to +0.39 | 0.07 (0.65) | 20210410D (+0.39) | +0.21 (0.10) | +0.23 [+0.13, +0.32], 0.07 (0.03) | +0.35 (0.081) | 0.27 |
| log Sigma_SFR | DM_host | literature R_e | no R_e from our fits | 36 | +0.28 (0.092) | 0.68 (0.84) | [-0.07, +0.59] | +0.20 (0.081) | +0.22 to +0.35 | 0.13 (0.59) | 20230907D (+0.22) | +0.22 (0.078) | +0.26 [+0.16, +0.35], 0.23 (0.00) | +0.38 (0.032) | 0.22 |
| Delta MS | DM_host | All | gas-rich hosts above the main sequence | 46 | +0.31 (0.036) | 0.043 (0.19) | [+0.02, +0.56] | +0.22 (0.032) | +0.28 to +0.38 | 0.54 (0.18) | 20210410D (+0.38) | +0.22 (0.022) | +0.26 [+0.18, +0.33], 0.32 (0.30) | +0.26 (0.100) | 0.15 |
| log Sigma* | DM_host | All | denser discs | 51 | +0.17 (0.24) | 0.64 (0.84) | [-0.13, +0.44] | +0.12 (0.23) | +0.12 to +0.21 | 0.00 (0.94) | 20230907D (+0.12) | +0.11 (0.24) | +0.14 [+0.05, +0.22], 0.04 (0.01) | +0.24 (0.10) | 0.41 |
| R_e | DM_host | All | longer path through larger discs | 55 | -0.10 (0.45) | 0.52 (0.75) | [-0.38, +0.18] | -0.08 (0.39) | -0.15 to -0.06 | 0.00 (1.00) | 20231020B (-0.06) | -0.08 (0.39) | -0.13 [-0.20, -0.05], 0.03 (0.01) | -0.29 (0.045) | 0.67 |
| offset / R_e | DM_host | All | Mannings+21, Marnoch+26: less gas at large offset | 45 | +0.17 (0.28) | 0.69 (0.84) | [-0.15, +0.44] | +0.11 (0.28) | +0.11 to +0.21 | 0.00 (0.98) | 20230907D (+0.11) | +0.12 (0.21) | +0.12 [+0.02, +0.21], 0.03 (0.01) | +0.10 (0.51) | 0.45 |
| A_V (old stars) | DM_host | All | dust column traces gas column | 30 | +0.35 (0.055) | 0.10 (0.23) | [-0.04, +0.68] | +0.25 (0.054) | +0.30 to +0.48 | 0.29 (0.34) | 20230124A (+0.48) | +0.25 (0.039) | +0.26 [+0.14, +0.36], 0.15 (0.11) | +0.27 (0.18) | 0.17 |
| A_V (young stars) | DM_host | All | dust column traces gas column | 31 | +0.10 (0.60) | 0.53 (0.75) | [-0.30, +0.48] | +0.08 (0.55) | +0.02 to +0.20 | 0.00 (1.00) | 20230124A (+0.20) | +0.07 (0.56) | +0.05 [-0.05, +0.15], 0.00 (0.00) | +0.11 (0.59) | 0.76 |

## Host halo

| x | y | sample | claim | N | r_S (p) | AD p (BH q) | 95% CI | tau (p) | JK1 | JK3: p<0.05 (max p) | most influential (r_S without) | partial tau given z (p) | DM_cosmic MC: r_S, p<0.05 (AD) | foreground-corrected | BH q |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| predicted host-halo DM | DM_host | All | halo gas sets the floor of DM_host | 44 | -0.06 (0.72) | 0.80 (0.94) | [-0.35, +0.24] | -0.04 (0.69) | -0.11 to -0.00 | 0.00 (1.00) | 20230907D (-0.11) | -0.04 (0.64) | -0.08 [-0.17, +0.01], 0.01 (0.00) | -0.10 (0.55) | 0.80 |

## Propagation

| x | y | sample | claim | N | r_S (p) | AD p (BH q) | 95% CI | tau (p) | JK1 | JK3: p<0.05 (max p) | most influential (r_S without) | partial tau given z (p) | DM_cosmic MC: r_S, p<0.05 (AD) | foreground-corrected | BH q |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| DM_host | log tau | All | Cordes+22: scattering grows with DM | 29 | +0.45 (0.014) | 0.002 (0.07) | [+0.10, +0.71] | +0.30 (0.022) | +0.40 to +0.54 | 0.88 (0.13) | 20211212A (+0.54) | +0.30 (0.016) | +0.38 [+0.27, +0.47], 0.55 (0.47) | +0.29 (0.13) | 0.12 |
| DM_host | log |RM| | All | Verdi+26: RM grows with DM | 41 | +0.49 (0.001) | 0.036 (0.19) | [+0.20, +0.72] | +0.36 (<0.001) | +0.46 to +0.56 | 1.00 (0.017) | 20240229A (+0.56) | +0.37 (0.034) | +0.38 [+0.28, +0.48], 0.77 (0.66) | +0.54 (<0.001) | 0.04 |
| log Sigma_SFR | log tau | All | turbulent star-forming gas scatters | 23 | -0.04 (0.86) | 0.097 (0.23) | [-0.48, +0.46] | -0.00 (0.98) | -0.15 to +0.08 | 0.00 (1.00) | 20210410D (+0.08) | +0.03 (0.83) |  |  | 0.92 |
| log Sigma* | log tau | All | Glowacki+25: positive | 26 | +0.08 (0.68) | 0.20 (0.37) | [-0.36, +0.51] | +0.06 (0.64) | -0.03 to +0.18 | 0.00 (1.00) | 20211127I (-0.03) | +0.08 (0.52) |  |  | 0.80 |
| offset / R_e | log tau | All | less scattering at large offset | 24 | -0.08 (0.71) | 0.61 (0.83) | [-0.47, +0.34] | -0.04 (0.77) | -0.23 to -0.00 | 0.00 (1.00) | 20211127I (-0.23) | -0.04 (0.82) |  |  | 0.80 |
| log Sigma_SFR | log |RM| | All | stronger fields in star-forming gas | 33 | +0.21 (0.24) | 0.35 (0.59) | [-0.13, +0.53] | +0.15 (0.22) | +0.15 to +0.27 | 0.00 (0.82) | 20220418A (+0.15) | +0.16 (0.16) |  |  | 0.41 |
| offset / R_e | log |RM| | All | less RM at large offset | 33 | +0.23 (0.20) | 1.00 (1.00) | [-0.16, +0.55] | +0.14 (0.25) | +0.17 to +0.29 | 0.01 (0.80) | 20230124A (+0.29) | +0.19 (0.095) |  |  | 0.38 |
| log |RM| | linear polarization | All | depolarization by Faraday rotation | 23 | +0.32 (0.13) | 0.13 (0.26) | [-0.09, +0.65] | +0.21 (0.15) | +0.28 to +0.45 | 0.08 (0.50) | 20230526A (+0.45) | +0.23 (0.16) |  |  | 0.27 |

## Sanity checks

| x | y | sample | claim | N | r_S (p) | AD p (BH q) | 95% CI | tau (p) | JK1 | JK3: p<0.05 (max p) | most influential (r_S without) | partial tau given z (p) | DM_cosmic MC: r_S, p<0.05 (AD) | foreground-corrected | BH q |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| log M* | log SFR | All | star-forming main sequence | 46 | +0.25 (0.097) | 0.41 | [-0.09, +0.55] | +0.18 (0.080) | +0.20 to +0.31 | 0.12 (0.57) | 20221012A (+0.31) | +0.15 (0.11) |  |  |  |
| log M* | stellar metallicity | All | mass-metallicity relation | 31 | +0.27 (0.14) | 0.52 | [-0.09, +0.55] | +0.18 (0.15) | +0.20 to +0.31 | 0.01 (0.69) | 20231123B (+0.20) | +0.18 (0.12) |  |  |  |
| log M* | R_e | All | mass-size relation | 51 | +0.33 (0.018) | 0.028 | [+0.03, +0.60] | +0.24 (0.012) | +0.29 to +0.39 | 0.90 (0.13) | 20230907D (+0.39) | +0.25 (0.010) |  |  |  |
| log M* | mass-weighted age | All | downsizing: massive hosts older | 30 | +0.38 (0.036) | 0.040 | [+0.02, +0.66] | +0.25 (0.054) | +0.33 to +0.46 | 0.50 (0.30) | 20210410D (+0.46) | +0.32 (0.007) |  |  |  |
| log M* | A_V (young stars) | All | dustier massive hosts | 31 | -0.10 (0.59) | 0.89 | [-0.44, +0.24] | -0.08 (0.55) | -0.21 to -0.05 | 0.00 (1.00) | 20230124A (-0.21) | -0.07 (0.51) |  |  |  |

## Host-halo budget

Predicted rest-frame DM from the FRB to the edge of its own halo (Moster+13 halo mass, mNFW gas, half chord at the projected offset), N = 44: median 14 pc cm^-3 (16-84%: 12-19); median DM_host 194; median halo share of DM_host (DM_host > 0) 0.07; hosts with DM_host below the halo prediction: 2 / 44.

## Sigma_SFR slope

Theil-Sen slope of log DM_host vs log Sigma_SFR (DM_host > 0, N = 44): 0.17 [0.00, 0.27]; the H-alpha emission-measure argument predicts 0.5.

## Main-sequence offsets

Delta MS (N = 46): median +0.01 dex; within +-0.3 dex: 0.50; more than 0.6 dex below: 8.
