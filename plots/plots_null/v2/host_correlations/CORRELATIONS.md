# Host and burst properties vs inclination (confirmed sample, N = 57)

Source: `scripts/run_host_property_correlations.py`. Inclination is the winning-leg $q$ through $\cos^2 i = (q^2 - q_0^2)/(1 - q_0^2)$, $q_0 = 0.2$. Correlations are against $\cos i$ (Spearman and Kendall are identical in magnitude, opposite in sign, against $i$). MC: 2000 Protocol B draws ($q \sim N(q, \sigma_{q,B})$, $i$ += $N(0, 5^\circ)$); `frac p<0.05` uses the first 500 draws. BH $q$ is the Benjamini–Hochberg false-discovery value across the rows of Table 1. Edge-on = $i > 60^\circ$.

DM$_{\rm host}$ = (DM − DM$_{\rm MW,NE2001}$ − median Macquart DM$_{\rm cosmic}(z)$) [× (1+z) for rest frame], no halo term, spectroscopic $z$ only. Literature values are picked with the priority rules of `build_host_property_table.py`; masses fall back to K-band / compiled values where no SED mass exists.

## 1. Each property vs $\cos i$

| property | literature | N | Spearman $r_S$ (p) | Kendall $\tau$ (p) | Pearson $r$ (p) | MC $r_S$ 16/50/84 | MC frac p<0.05 | BH q | median edge / face (MWU p) |
|---|---|---|---|---|---|---|---|---|---|
| DM_host rest (Macquart median, NE2001) | Marnoch+26; Bhardwaj+24; Li+25 | 55 | -0.065 (0.638) | -0.045 (0.632) | 0.032 (0.816) | -0.10 / -0.05 / -6.49e-04 | 0.00 | 0.791 | 244.847 / 184.420 (0.410; 16/39) |
| DM_host observer frame | Marnoch+26 | 55 | -0.059 (0.670) | -0.041 (0.663) | 0.005 (0.970) | -0.10 / -0.04 / 0.01 | 0.00 | 0.791 | 171.806 / 146.773 (0.420; 16/39) |
| DM_host as published | Marnoch+26; Verdi+26 | 43 | -0.131 (0.402) | -0.101 (0.341) | -0.054 (0.733) | -0.17 / -0.12 / -0.06 | 0.00 | 0.791 | 208.510 / 155.790 (0.310; 12/31) |
| DM - DM_MW | Li+25 | 56 | -0.051 (0.707) | -0.038 (0.677) | 0.008 (0.954) | -0.08 / -0.04 / 0.01 | 0.00 | 0.791 | 441.320 / 386.840 (0.722; 17/39) |
| DM_obs | Bhardwaj+24 | 56 | -0.060 (0.658) | -0.040 (0.666) | 0.014 (0.916) | -0.09 / -0.04 / 8.64e-04 | 0.00 | 0.791 | 491.200 / 441.200 (0.762; 17/39) |
| log tau (scattering) | Glowacki+25; Ocker+22 | 29 | -0.167 (0.385) | -0.109 (0.409) | -0.158 (0.414) | -0.22 / -0.16 / -0.09 | 0.00 | 0.791 | 0.815 / -0.222 (0.125; 9/20) |
| log |RM| | Verdi+26; Mannings+23 | 41 | -0.018 (0.912) | -0.001 (0.991) | -0.061 (0.706) | -0.08 / -0.01 / 0.05 | 0.00 | 0.918 | 2.211 / 2.138 (0.900; 9/32) |
| projected offset [kpc] | Marnoch+26; Mannings+21 | 47 | -0.122 (0.413) | -0.080 (0.430) | -0.136 (0.363) | -0.16 / -0.11 / -0.06 | 0.00 | 0.791 | 5.035 / 4.260 (0.733; 12/35) |
| offset / R_e | Marnoch+26; Gordon+25 | 45 | 0.093 (0.543) | 0.073 (0.481) | 0.162 (0.287) | 0.05 / 0.10 / 0.14 | 0.00 | 0.791 | 0.634 / 1.180 (0.037; 10/35) |
| log M* | Leung+25; Li+25 | 53 | -0.054 (0.699) | -0.039 (0.684) | 8.53e-04 (0.995) | -0.09 / -0.05 / -0.00 | 0.00 | 0.791 | 10.170 / 10.140 (0.389; 16/37) |
| SFR | Li+25 | 46 | 0.123 (0.417) | 0.070 (0.495) | 0.233 (0.119) | 0.08 / 0.13 / 0.17 | 0.00 | 0.791 | 1.780 / 1.770 (0.348; 13/33) |
| log sSFR | Li+25 | 46 | 0.138 (0.362) | 0.091 (0.373) | 0.213 (0.155) | 0.09 / 0.14 / 0.18 | 0.00 | 0.791 | -9.944 / -9.713 (0.884; 13/33) |
| log Sigma* (M*/2 pi R_e^2) | Glowacki+25 | 41 | 0.141 (0.379) | 0.103 (0.345) | 0.186 (0.244) | 0.08 / 0.14 / 0.20 | 0.00 | 0.791 | 7.977 / 8.292 (0.092; 10/31) |
| mass-weighted age | Glowacki+25 | 30 | -0.183 (0.334) | -0.120 (0.353) | -0.223 (0.236) | -0.25 / -0.19 / -0.13 | 0.00 | 0.791 | 4.860 / 5.150 (0.839; 9/21) |
| stellar metallicity | Glowacki+25 | 31 | 0.138 (0.460) | 0.093 (0.464) | 0.205 (0.269) | 0.07 / 0.13 / 0.19 | 0.00 | 0.791 | -0.440 / -0.520 (0.663; 9/22) |
| u - r colour | Marnoch+26 (colour) | 25 | -0.326 (0.112) | -0.217 (0.129) | -0.364 (0.073) | -0.38 / -0.32 / -0.25 | 0.09 | 0.791 | 1.805 / 1.660 (0.475; 8/17) |
| R_e [kpc] | Li+25 (projected area) | 42 | -0.286 (0.067) | -0.211 (0.050) | -0.308 (0.047) | -0.32 / -0.27 / -0.22 | 0.27 | 0.791 | 5.215 / 3.775 (0.059; 10/32) |
| redshift | selection | 56 | -0.014 (0.918) | -0.011 (0.904) | 0.050 (0.713) | -0.05 / -0.01 / 0.04 | 0.00 | 0.918 | 0.223 / 0.242 (0.452; 16/40) |
| m_r (ours) | selection; Marnoch m<21 cut | 57 | -0.203 (0.129) | -0.148 (0.106) | -0.191 (0.154) | -0.24 / -0.19 / -0.14 | 0.03 | 0.791 | 20.060 / 19.485 (0.299; 17/40) |
| Sersic n (ours) | selection; Marnoch n<2 cut | 57 | 0.314 (0.017) | 0.220 (0.019) | 0.218 (0.103) | 0.26 / 0.30 / 0.35 | 0.84 | 0.482 | 0.700 / 1.000 (0.050; 17/40) |
| PATH P(O|x) | association | 49 | 0.057 (0.699) | 0.036 (0.732) | 0.225 (0.120) | 0.01 / 0.06 / 0.11 | 0.00 | 0.791 | 0.990 / 0.990 (0.224; 14/35) |
| S/N | selection | 35 | 0.192 (0.269) | 0.132 (0.267) | -0.010 (0.953) | 0.10 / 0.16 / 0.22 | 0.00 | 0.791 | 13.600 / 16.100 (0.298; 10/25) |
| fluence | selection | 33 | 0.149 (0.407) | 0.112 (0.360) | 0.159 (0.377) | 0.06 / 0.13 / 0.18 | 0.00 | 0.791 | 16.100 / 28.000 (0.462; 8/25) |
| width | selection; Glowacki+25 | 35 | 0.130 (0.455) | 0.086 (0.476) | 0.215 (0.215) | 0.05 / 0.12 / 0.19 | 0.00 | 0.791 | 2.600 / 5.300 (0.140; 9/26) |
| DM_MW (NE2001) | null check | 57 | -0.118 (0.381) | -0.072 (0.428) | 0.032 (0.811) | -0.16 / -0.11 / -0.06 | 0.00 | 0.791 | 40.250 / 38.685 (0.407; 17/40) |
| |b| Galactic | null check | 57 | 0.086 (0.525) | 0.038 (0.674) | 0.079 (0.557) | 0.03 / 0.08 / 0.13 | 0.00 | 0.791 | 38.324 / 44.606 (0.559; 17/40) |
| Dec | null check (survey footprint) | 57 | -0.061 (0.653) | -0.036 (0.695) | -0.041 (0.764) | -0.12 / -0.07 / -0.02 | 0.00 | 0.791 | 37.616 / 3.735 (0.427; 17/40) |
| RA | null check | 57 | -0.028 (0.834) | -0.009 (0.918) | -0.018 (0.893) | -0.07 / -0.03 / 0.02 | 0.00 | 0.898 | 183.980 / 204.512 (0.773; 17/40) |

## 2. Marnoch et al. (2026) DM$_{\rm host}$ test

Marnoch+2026 (CRAFT, N = 34): $r_S \approx -0.44$, $p \approx 0.01$ between rest-frame DM$_{\rm host}$ and $q$ (more DM in edge-on hosts). Partial $r_S$ controls for $z$.

| DM_host | vs | subsample | N | Spearman (p) | Pearson (p) | partial Spearman given z (p) |
|---|---|---|---|---|---|---|
| dm_host_rest | q | all | 55 | -0.065 (0.638) | 0.033 (0.813) | -0.067 (0.630) |
| dm_host_rest | q | m_r < 21 | 45 | -0.093 (0.543) | 0.037 (0.810) | -0.096 (0.536) |
| dm_host_rest | q | n < 2 | 46 | -0.037 (0.806) | 0.056 (0.713) | -0.039 (0.800) |
| dm_host_rest | q | m_r < 21 & n < 2 | 36 | -0.071 (0.682) | 0.074 (0.670) | -0.073 (0.675) |
| dm_host_rest | q | spectroscopic z <= 0.3 | 37 | -0.065 (0.700) | -0.045 (0.793) | -0.030 (0.864) |
| dm_host_rest | i | all | 55 | 0.065 (0.638) | -0.039 (0.775) | 0.067 (0.630) |
| dm_host_rest | i | m_r < 21 | 45 | 0.093 (0.543) | -0.047 (0.757) | 0.096 (0.536) |
| dm_host_rest | i | n < 2 | 46 | 0.037 (0.806) | -0.066 (0.661) | 0.039 (0.800) |
| dm_host_rest | i | m_r < 21 & n < 2 | 36 | 0.071 (0.682) | -0.090 (0.600) | 0.073 (0.675) |
| dm_host_rest | i | spectroscopic z <= 0.3 | 37 | 0.065 (0.700) | 0.042 (0.804) | 0.030 (0.864) |
| dm_host_obs | q | all | 55 | -0.059 (0.670) | 0.008 (0.954) | -0.057 (0.685) |
| dm_host_obs | q | m_r < 21 | 45 | -0.062 (0.687) | 0.012 (0.939) | -0.062 (0.689) |
| dm_host_obs | q | n < 2 | 46 | -0.028 (0.855) | 0.039 (0.798) | -0.022 (0.887) |
| dm_host_obs | q | m_r < 21 & n < 2 | 36 | -0.025 (0.886) | 0.059 (0.731) | -0.024 (0.892) |
| dm_host_obs | q | spectroscopic z <= 0.3 | 37 | -0.047 (0.784) | -0.043 (0.801) | -0.033 (0.847) |
| dm_host_obs | i | all | 55 | 0.059 (0.670) | -0.020 (0.885) | 0.057 (0.685) |
| dm_host_obs | i | m_r < 21 | 45 | 0.062 (0.687) | -0.027 (0.860) | 0.062 (0.689) |
| dm_host_obs | i | n < 2 | 46 | 0.028 (0.855) | -0.056 (0.714) | 0.022 (0.887) |
| dm_host_obs | i | m_r < 21 & n < 2 | 36 | 0.025 (0.886) | -0.081 (0.639) | 0.024 (0.892) |
| dm_host_obs | i | spectroscopic z <= 0.3 | 37 | 0.047 (0.784) | 0.035 (0.839) | 0.033 (0.847) |
| dm_host | q | all | 43 | -0.131 (0.402) | -0.053 (0.735) | -0.125 (0.432) |
| dm_host | q | m_r < 21 | 38 | -0.092 (0.581) | -0.038 (0.823) | -0.073 (0.669) |
| dm_host | q | n < 2 | 35 | -0.144 (0.409) | -0.099 (0.572) | -0.143 (0.419) |
| dm_host | q | m_r < 21 & n < 2 | 30 | -0.099 (0.604) | -0.074 (0.696) | -0.092 (0.635) |
| dm_host | q | spectroscopic z <= 0.3 | 32 | -0.103 (0.575) | -0.071 (0.701) | -0.051 (0.785) |
| dm_host | i | all | 43 | 0.131 (0.402) | 0.042 (0.790) | 0.125 (0.432) |
| dm_host | i | m_r < 21 | 38 | 0.092 (0.581) | 0.023 (0.891) | 0.073 (0.669) |
| dm_host | i | n < 2 | 35 | 0.144 (0.409) | 0.074 (0.671) | 0.143 (0.419) |
| dm_host | i | m_r < 21 & n < 2 | 30 | 0.099 (0.604) | 0.045 (0.812) | 0.092 (0.635) |
| dm_host | i | spectroscopic z <= 0.3 | 32 | 0.103 (0.575) | 0.065 (0.724) | 0.051 (0.785) |

## 3. Literature host-property relations on this sample (no inclination)

| x | y | claim | N | Spearman (p) | Pearson (p) |
|---|---|---|---|---|---|
| logmstar | dm_host_rest | Leung+25: DM_host anticorrelates with M* | 51 | -0.007 (0.963) | 0.050 (0.727) |
| log_ssfr | dm_exgal | Li+25: DM_exc tightly correlates with sSFR | 46 | 0.257 (0.084) | 0.177 (0.240) |
| sfr | dm_exgal | Li+25: DM_exc vs SFR | 46 | 0.280 (0.060) | 0.207 (0.167) |
| log_sigma_star | log_tau | Glowacki+25: tau correlates with Sigma* | 20 | 0.261 (0.266) | 0.098 (0.681) |
| tm_gyr | log_tau | Glowacki+25: tau correlates with mass-weighted age | 13 | 0.085 (0.782) | 0.161 (0.599) |
| metallicity_logz | log_tau | Glowacki+25: weak tau-metallicity | 14 | 0.008 (0.979) | -0.086 (0.769) |
| dm_host_rest | log_tau | Cordes+22 / Day+20: tau vs DM | 29 | 0.450 (0.014) | 0.450 (0.014) |
| dm_host_rest | log_abs_rm | Verdi+26: DM_host vs RM_host | 41 | 0.490 (0.001) | 0.423 (0.006) |
| offset_normalized | dm_host_rest | Marnoch+26 / Mannings+21: DM_host vs offset | 44 | 0.234 (0.126) | 0.237 (0.121) |
| z | dm_host_rest | Verdi+26: no z evolution of DM_host | 55 | 0.082 (0.551) | 0.126 (0.358) |

## 4. $\cos i$ by survey

| survey | N | median $\cos i$ |
|---|---|---|
| ASKAP | 26 | 0.647 |
| CHIME | 1 | 0.769 |
| DSA-110 | 25 | 0.610 |
| MeerKAT | 5 | 0.567 |

Kruskal–Wallis over surveys with N ≥ 5 (ASKAP, DSA-110, MeerKAT): H = 0.38, p = 0.826.

Hosts without DM$_{\rm host}$: 20231210F, 20240615B.
