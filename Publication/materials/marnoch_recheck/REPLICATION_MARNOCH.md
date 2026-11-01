# Marnoch et al. (2026) analysis repeated on our hosts

Source: `scripts/marnoch_replication.py`. DM definitions as their eqs. 3–4 (NE2001, 40 pc cm$^{-3}$ Milky Way halo, mean Macquart $\langle$DM$_{\rm cosmic}\rangle$); LTSFIT, 2.6σ clip; $i$ from $q$ with $q_0 = 0.2$. Hosts with $z$ and DM: 55. Their bulge-fit hosts removed in the † subsample: 20190608B, 20220725A, 20231226A.

## DM vs axis ratio, inclination and cos i (their Tables 3–6)

| subsample | y | x | N | $r_P$ ($p$) | $r_S$ ($p$) | $r_S$ 95% CI | LTS slope | LTS intercept (pivot) | clipped |
|---|---|---|---|---|---|---|---|---|---|
| all | dm_mac_obs | q | 55 | 0.00 (0.990) | -0.04 (0.797) | [-0.30, 0.23] | -91.3 ± 76.7 | 85.9 (0.7) | 5 |
| all | dm_mac_rest | q | 55 | 0.02 (0.890) | -0.04 (0.765) | [-0.30, 0.23] | -1.0 ± 64.9 | 89.7 (0.7) | 13 |
| all | dm_exgal | q | 55 | 0.06 (0.650) | -0.01 (0.947) | [-0.27, 0.26] | 17.5 ± 107.8 | 344.1 (0.7) | 4 |
| all | dm_frb | q | 55 | 0.07 (0.625) | -0.02 (0.893) | [-0.28, 0.25] | -17.0 ± 110.3 | 421.8 (0.7) | 3 |
| all | dm_frb_minus_cosmic | q | 55 | 0.01 (0.960) | -0.02 (0.897) | [-0.28, 0.25] | -8.9 ± 65.6 | 160.8 (0.7) | 9 |
| all | dm_mac_rest | i | 55 | -0.03 (0.821) | 0.04 (0.765) | [-0.23, 0.30] | -0.1 ± 0.8 | 88.3 (50.1) | 13 |
| all | dm_mac_rest | cos_i | 55 | 0.02 (0.906) | -0.04 (0.765) | [-0.30, 0.23] | -4.4 ± 57.9 | 89.7 (0.6) | 13 |
| m_r < 21 | dm_mac_obs | q | 45 | 0.01 (0.961) | -0.03 (0.857) | [-0.32, 0.27] | -90.6 ± 79.3 | 98.0 (0.7) | 4 |
| m_r < 21 | dm_mac_rest | q | 45 | 0.03 (0.857) | -0.04 (0.806) | [-0.33, 0.26] | -105.6 ± 72.3 | 109.1 (0.7) | 9 |
| m_r < 21 | dm_exgal | q | 45 | 0.05 (0.748) | -0.05 (0.738) | [-0.34, 0.25] | -21.5 ± 110.6 | 321.6 (0.7) | 3 |
| m_r < 21 | dm_frb | q | 45 | 0.06 (0.712) | -0.06 (0.698) | [-0.35, 0.24] | -51.2 ± 113.6 | 400.7 (0.7) | 2 |
| m_r < 21 | dm_frb_minus_cosmic | q | 45 | 0.02 (0.914) | -0.03 (0.869) | [-0.32, 0.27] | -10.8 ± 60.0 | 187.9 (0.7) | 9 |
| m_r < 21 | dm_mac_rest | i | 45 | -0.04 (0.780) | 0.04 (0.806) | [-0.26, 0.33] | 1.2 ± 1.2 | 118.3 (48.4) | 4 |
| m_r < 21 | dm_mac_rest | cos_i | 45 | 0.02 (0.876) | -0.04 (0.806) | [-0.33, 0.26] | -96.0 ± 64.0 | 107.0 (0.7) | 9 |
| n < 2 | dm_mac_obs | q | 46 | 0.03 (0.848) | -0.01 (0.946) | [-0.30, 0.28] | -4.8 ± 61.3 | 68.7 (0.6) | 7 |
| n < 2 | dm_mac_rest | q | 46 | 0.04 (0.803) | -0.02 (0.902) | [-0.31, 0.27] | 16.9 ± 62.7 | 63.2 (0.6) | 19 |
| n < 2 | dm_exgal | q | 46 | 0.10 (0.524) | 0.02 (0.903) | [-0.27, 0.31] | 24.2 ± 113.1 | 332.3 (0.6) | 2 |
| n < 2 | dm_frb | q | 46 | 0.09 (0.540) | 0.01 (0.952) | [-0.28, 0.30] | 21.5 ± 113.2 | 419.4 (0.6) | 2 |
| n < 2 | dm_frb_minus_cosmic | q | 46 | 0.02 (0.885) | 0.01 (0.964) | [-0.28, 0.30] | -6.8 ± 62.4 | 154.0 (0.6) | 7 |
| n < 2 | dm_mac_rest | i | 46 | -0.05 (0.718) | 0.02 (0.902) | [-0.27, 0.31] | -0.3 ± 0.8 | 92.3 (51.2) | 9 |
| n < 2 | dm_mac_rest | cos_i | 46 | 0.03 (0.833) | -0.02 (0.902) | [-0.31, 0.27] | 12.2 ± 56.4 | 63.6 (0.6) | 19 |
| bulge fits removed | dm_mac_obs | q | 52 | 0.00 (0.993) | -0.04 (0.762) | [-0.31, 0.23] | -93.9 ± 78.8 | 95.0 (0.6) | 5 |
| bulge fits removed | dm_mac_rest | q | 52 | 0.02 (0.889) | -0.05 (0.739) | [-0.32, 0.23] | -109.3 ± 94.7 | 111.6 (0.6) | 5 |
| bulge fits removed | dm_exgal | q | 52 | 0.07 (0.604) | 0.01 (0.916) | [-0.26, 0.29] | 28.1 ± 110.4 | 348.4 (0.6) | 4 |
| bulge fits removed | dm_frb | q | 52 | 0.08 (0.573) | 0.00 (0.975) | [-0.27, 0.28] | -7.3 ± 113.0 | 429.7 (0.6) | 3 |
| bulge fits removed | dm_frb_minus_cosmic | q | 52 | 0.01 (0.956) | -0.03 (0.858) | [-0.30, 0.25] | -14.5 ± 67.6 | 161.3 (0.6) | 9 |
| bulge fits removed | dm_mac_rest | i | 52 | -0.03 (0.826) | 0.05 (0.739) | [-0.23, 0.32] | 1.2 ± 1.1 | 108.1 (51.2) | 5 |
| bulge fits removed | dm_mac_rest | cos_i | 52 | 0.02 (0.902) | -0.05 (0.739) | [-0.32, 0.23] | -100.0 ± 84.2 | 108.6 (0.6) | 5 |
| shared with Marnoch | dm_mac_obs | q | 20 | -0.01 (0.951) | -0.10 (0.677) | [-0.52, 0.36] | 196.4 ± 84.0 | 30.2 (0.7) | 8 |
| shared with Marnoch | dm_mac_rest | q | 20 | -0.04 (0.880) | -0.10 (0.663) | [-0.52, 0.36] | 218.6 ± 95.7 | 34.5 (0.7) | 8 |
| shared with Marnoch | dm_exgal | q | 20 | -0.25 (0.291) | -0.28 (0.235) | [-0.64, 0.19] | -521.6 ± 48.6 | 285.5 (0.7) | 8 |
| shared with Marnoch | dm_frb | q | 20 | -0.22 (0.348) | -0.32 (0.175) | [-0.67, 0.15] | -499.9 ± 45.4 | 358.8 (0.7) | 9 |
| shared with Marnoch | dm_frb_minus_cosmic | q | 20 | 0.03 (0.912) | 0.03 (0.905) | [-0.42, 0.47] | 197.7 ± 60.9 | 122.3 (0.7) | 9 |
| shared with Marnoch | dm_mac_rest | i | 20 | 0.01 (0.963) | 0.10 (0.663) | [-0.36, 0.52] | -2.6 ± 0.9 | 22.7 (49.7) | 8 |
| shared with Marnoch | dm_mac_rest | cos_i | 20 | -0.04 (0.875) | -0.10 (0.663) | [-0.52, 0.36] | 202.9 ± 90.9 | 17.7 (0.6) | 8 |
| all | dm_mac_rest | offset_kpc | 45 | -0.23 (0.123) | -0.13 (0.380) | [-0.41, 0.17] | — | — |  |
| m_r < 21 | dm_mac_rest | offset_kpc | 37 | -0.18 (0.273) | -0.04 (0.799) | [-0.36, 0.28] | — | — |  |
| all | dm_mac_rest | offset_normalized | 44 | 0.14 (0.365) | 0.10 (0.510) | [-0.20, 0.39] | — | — |  |
| m_r < 21 | dm_mac_rest | offset_normalized | 37 | 0.40 (0.013) | 0.11 (0.521) | [-0.22, 0.42] | — | — |  |
| all | dm_mac_rest | u_r_color | 25 | 0.11 (0.596) | 0.11 (0.614) | [-0.30, 0.48] | — | — |  |
| m_r < 21 | dm_mac_rest | u_r_color | 21 | 0.07 (0.763) | 0.01 (0.951) | [-0.42, 0.44] | — | — |  |
| all | dm_mac_rest | g_r_color | 18 | -0.22 (0.382) | -0.36 (0.142) | [-0.71, 0.13] | — | — |  |
| m_r < 21 | dm_mac_rest | g_r_color | 17 | -0.21 (0.414) | -0.34 (0.178) | [-0.71, 0.16] | — | — |  |

Their full-sample relation: slope $-309.3 \pm 69.5$, intercept $155.9 \pm 17.4$ at $q = 0.5$, $r_P = -0.40$, $r_S = -0.44$ (N = 34).

## Macquart-relation scatter (their Figure 13)

RMSE of DM$_{\rm exgal}$ about $\langle$DM$_{\rm cosmic}\rangle(z)$, i.e. of observer-frame DM$_{\rm host}^{\rm Macquart}$, before and after subtracting a host term converted to the observer frame.

| sample | correction | N | RMSE (pc cm$^{-3}$) |
|---|---|---|---|
| Marnoch Table 1 (34) | none | 34 | 181 |
| Marnoch Table 1 (34) | Marnoch model DM(q) | 34 | 107 |
| Marnoch Table 1 (34) | constant 155.9 (their intercept, no q) | 34 | 116 |
| Marnoch Table 1 (34) | constant, best-fitting mean (no q) | 34 | 115 |
| this work (55) | none | 55 | 186 |
| this work (55) | Marnoch model DM(q) | 55 | 149 |
| this work (55) | constant 155.9 (their intercept, no q) | 55 | 138 |
| this work (55) | constant, best-fitting mean (no q) | 55 | 138 |
| this work (55) | own LTS fit DM(q) | 55 | 147 |

## Which properties can exclude a correlation as strong as Marnoch's?

Fisher 95% CI on $r_S$ with $\cos i$; power = probability of $p < 0.05$ if the true $|r_S| = 0.44$.

| property | N | $r_S$ ($p$) | 95% CI | excludes $-0.44$ | excludes $\pm 0.44$ | power at $|r_S| = 0.44$ |
|---|---|---|---|---|---|---|
| DM_obs | 56 | -0.06 (0.658) | [-0.32, 0.21] | yes | yes | 0.93 |
| DM - DM_MW | 56 | -0.05 (0.707) | [-0.31, 0.21] | yes | yes | 0.93 |
| |b| Galactic | 56 | 0.05 (0.719) | [-0.22, 0.31] | yes | yes | 0.93 |
| m_r (ours) | 56 | -0.17 (0.211) | [-0.41, 0.10] | yes | yes | 0.93 |
| Sersic n (ours) | 56 | 0.30 (0.024) | [0.04, 0.52] | yes | no | 0.93 |
| DM_MW (NE2001) | 56 | -0.08 (0.557) | [-0.34, 0.19] | yes | yes | 0.93 |
| RA | 56 | 0.02 (0.910) | [-0.25, 0.28] | yes | yes | 0.93 |
| Dec | 56 | -0.05 (0.734) | [-0.31, 0.22] | yes | yes | 0.93 |
| DM_host rest (Macquart median, NE2001) | 55 | -0.06 (0.638) | [-0.32, 0.20] | yes | yes | 0.93 |
| redshift | 55 | 0.02 (0.858) | [-0.24, 0.29] | yes | yes | 0.93 |
| DM_host observer frame | 55 | -0.06 (0.670) | [-0.32, 0.21] | yes | yes | 0.93 |
| log M* | 52 | -0.11 (0.427) | [-0.37, 0.17] | yes | yes | 0.91 |
| PATH P(O|x) | 49 | 0.06 (0.699) | [-0.23, 0.33] | yes | yes | 0.89 |
| projected offset [kpc] | 46 | -0.08 (0.620) | [-0.36, 0.22] | yes | yes | 0.87 |
| SFR | 46 | 0.12 (0.417) | [-0.17, 0.40] | yes | yes | 0.87 |
| log sSFR | 46 | 0.14 (0.362) | [-0.16, 0.41] | yes | yes | 0.87 |
| offset / R_e | 44 | 0.16 (0.289) | [-0.14, 0.44] | yes | yes | 0.86 |
| DM_host as published | 43 | -0.13 (0.402) | [-0.42, 0.18] | yes | yes | 0.85 |
| log |RM| | 41 | -0.02 (0.912) | [-0.32, 0.29] | yes | yes | 0.83 |
| R_e [kpc] | 41 | -0.33 (0.033) | [-0.58, -0.03] | no | no | 0.83 |
| log Sigma* (M*/2 pi R_e^2) | 40 | 0.09 (0.576) | [-0.23, 0.39] | yes | yes | 0.82 |
| 1.23 |RM| (1+z)^2 / DM_host,rest | 39 | -0.08 (0.634) | [-0.38, 0.24] | yes | yes | 0.81 |
| S/N | 35 | 0.19 (0.269) | [-0.15, 0.49] | yes | no | 0.76 |
| width | 35 | 0.13 (0.455) | [-0.21, 0.44] | yes | no | 0.76 |
| fluence | 33 | 0.15 (0.407) | [-0.20, 0.47] | yes | no | 0.73 |
| stellar metallicity | 31 | 0.14 (0.460) | [-0.23, 0.47] | yes | no | 0.70 |
| A_V young stars (birth clouds + diffuse; Prospector) | 31 | -0.04 (0.837) | [-0.39, 0.32] | yes | yes | 0.70 |
| mass-weighted age | 30 | -0.18 (0.334) | [-0.51, 0.19] | no | no | 0.69 |
| A_V old stars (diffuse ISM; Prospector) | 30 | -0.46 (0.010) | [-0.70, -0.12] | no | no | 0.69 |
| log tau (scattering) | 29 | -0.17 (0.385) | [-0.50, 0.21] | no | no | 0.67 |
| u - r colour | 25 | -0.33 (0.112) | [-0.64, 0.08] | no | no | 0.60 |
| linear polarisation fraction | 24 | -0.33 (0.114) | [-0.65, 0.08] | no | no | 0.58 |
| |RM_host| rest (published) | 21 | 0.06 (0.803) | [-0.38, 0.48] | yes | no | 0.52 |
| g - r colour | 18 | -0.21 (0.401) | [-0.62, 0.28] | no | no | 0.45 |
| B_par host (published) | 17 | -0.32 (0.209) | [-0.69, 0.19] | no | no | 0.42 |
| circular polarisation fraction | 17 | -0.19 (0.468) | [-0.61, 0.32] | no | no | 0.42 |
| PAH fraction q_PAH | 16 | 0.36 (0.175) | [-0.17, 0.72] | yes | no | 0.40 |
| log tau rest 1 GHz (published) | 15 | 0.22 (0.423) | [-0.33, 0.66] | yes | no | 0.37 |
| gas metallicity 12+log O/H | 11 | -0.06 (0.872) | [-0.63, 0.56] | no | no | 0.27 |
| scintillation bandwidth | 7 | -0.50 (0.253) | [-0.91, 0.41] | no | no | 0.15 |
| H-alpha EW | 6 | -0.03 (0.957) | [-0.82, 0.80] | no | no | 0.13 |
