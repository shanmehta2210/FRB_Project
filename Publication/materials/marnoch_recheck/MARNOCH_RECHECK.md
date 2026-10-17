# Re-check of Marnoch et al. (2026) DM_host vs axis ratio

Source: `scripts/recheck_marnoch2026.py`, using the per-host values of their Table 1 (N = 34). DM_host is their observer-frame Macquart proxy; rest frame = observer × (1+z). Hosts with $q < 0.2$ that they set to $i = 90^\circ$: 20220105A, 20231230D, 20240201A, 20240203C. Our sample cut is $q > 0.2$, so these are dropped in the second variant; the third also drops the repeater 20190711A. $i$ is their tabulated value (as published) or Hubble with $q_0 = 0.2$. 95% CI on $r_S$ from the Fisher transform.

## 1. Reproduction of the published numbers

| subsample | DM_host | vs | N | $r_P$ (p) | $r_S$ (p) | paper |
|---|---|---|---|---|---|---|
| full | dm_host_obs | q | 34 | -0.402 (0.019) | -0.427 (0.012) | r_P = -0.399 (p = 0.02), r_S = -0.425 (p = 0.012) |
| full | dm_host_obs | i_table | 34 | 0.383 (0.025) | 0.433 (0.011) |  |
| full | dm_host_rest | q | 34 | -0.406 (0.017) | -0.431 (0.011) | r_P = -0.402 (p = 0.018), r_S = -0.437 (p = 0.01) |
| full | dm_host_rest | i_table | 34 | 0.385 (0.024) | 0.438 (0.010) | p < 0.05 |
| m_galfit <= 21 | dm_host_obs | q | 18 | -0.539 (0.021) | -0.560 (0.016) |  |
| m_galfit <= 21 | dm_host_obs | i_table | 18 | 0.514 (0.029) | 0.560 (0.016) |  |
| m_galfit <= 21 | dm_host_rest | q | 18 | -0.542 (0.020) | -0.548 (0.019) | r_P = -0.538, r_S = -0.544 (p ~ 0.02) |
| m_galfit <= 21 | dm_host_rest | i_table | 18 | 0.510 (0.031) | 0.545 (0.019) |  |
| n < 2 | dm_host_obs | q | 28 | -0.427 (0.024) | -0.385 (0.043) |  |
| n < 2 | dm_host_obs | i_table | 28 | 0.405 (0.032) | 0.396 (0.037) |  |
| n < 2 | dm_host_rest | q | 28 | -0.432 (0.022) | -0.394 (0.038) | r_P = -0.43 (p = 0.022), r_S = -0.40 (p = 0.036) |
| n < 2 | dm_host_rest | i_table | 28 | 0.409 (0.031) | 0.404 (0.033) |  |
| no bulge fits | dm_host_obs | q | 31 | -0.404 (0.024) | -0.430 (0.016) |  |
| no bulge fits | dm_host_obs | i_table | 31 | 0.386 (0.032) | 0.437 (0.014) |  |
| no bulge fits | dm_host_rest | q | 31 | -0.407 (0.023) | -0.424 (0.018) | significant, moderately strong |
| no bulge fits | dm_host_rest | i_table | 31 | 0.387 (0.032) | 0.430 (0.016) |  |

## 2. drop q < 0.2

| subsample | DM_host | vs | N | $r_P$ (p) | $r_S$ (p) | $r_S$ 95% CI |
|---|---|---|---|---|---|---|
| full | dm_host_obs | q | 30 | -0.359 (0.051) | -0.347 (0.060) | [-0.64, 0.04] |
| full | dm_host_obs | i_hubble | 30 | 0.339 (0.067) | 0.347 (0.060) | [-0.04, 0.64] |
| full | dm_host_rest | q | 30 | -0.364 (0.048) | -0.358 (0.052) | [-0.65, 0.03] |
| full | dm_host_rest | i_hubble | 30 | 0.345 (0.062) | 0.358 (0.052) | [-0.03, 0.65] |
| m_galfit <= 21 | dm_host_obs | q | 16 | -0.405 (0.119) | -0.432 (0.094) | [-0.78, 0.11] |
| m_galfit <= 21 | dm_host_obs | i_hubble | 16 | 0.387 (0.139) | 0.432 (0.094) | [-0.11, 0.78] |
| m_galfit <= 21 | dm_host_rest | q | 16 | -0.430 (0.096) | -0.468 (0.068) | [-0.79, 0.07] |
| m_galfit <= 21 | dm_host_rest | i_hubble | 16 | 0.410 (0.115) | 0.468 (0.068) | [-0.07, 0.79] |
| n < 2 | dm_host_obs | q | 24 | -0.339 (0.105) | -0.267 (0.207) | [-0.62, 0.18] |
| n < 2 | dm_host_obs | i_hubble | 24 | 0.317 (0.132) | 0.267 (0.207) | [-0.18, 0.62] |
| n < 2 | dm_host_rest | q | 24 | -0.355 (0.089) | -0.290 (0.170) | [-0.64, 0.15] |
| n < 2 | dm_host_rest | i_hubble | 24 | 0.334 (0.110) | 0.290 (0.170) | [-0.15, 0.64] |
| no bulge fits | dm_host_obs | q | 27 | -0.366 (0.060) | -0.348 (0.075) | [-0.66, 0.06] |
| no bulge fits | dm_host_obs | i_hubble | 27 | 0.347 (0.076) | 0.348 (0.075) | [-0.06, 0.66] |
| no bulge fits | dm_host_rest | q | 27 | -0.372 (0.056) | -0.360 (0.065) | [-0.66, 0.05] |
| no bulge fits | dm_host_rest | i_hubble | 27 | 0.353 (0.071) | 0.360 (0.065) | [-0.05, 0.66] |

## 3. drop q < 0.2 and repeater

| subsample | DM_host | vs | N | $r_P$ (p) | $r_S$ (p) | $r_S$ 95% CI |
|---|---|---|---|---|---|---|
| full | dm_host_obs | q | 29 | -0.380 (0.042) | -0.386 (0.039) | [-0.67, 0.00] |
| full | dm_host_obs | i_hubble | 29 | 0.361 (0.054) | 0.386 (0.039) | [-0.00, 0.67] |
| full | dm_host_rest | q | 29 | -0.385 (0.039) | -0.398 (0.032) | [-0.68, -0.01] |
| full | dm_host_rest | i_hubble | 29 | 0.366 (0.051) | 0.398 (0.032) | [0.01, 0.68] |
| m_galfit <= 21 | dm_host_obs | q | 16 | -0.405 (0.119) | -0.432 (0.094) | [-0.78, 0.11] |
| m_galfit <= 21 | dm_host_obs | i_hubble | 16 | 0.387 (0.139) | 0.432 (0.094) | [-0.11, 0.78] |
| m_galfit <= 21 | dm_host_rest | q | 16 | -0.430 (0.096) | -0.468 (0.068) | [-0.79, 0.07] |
| m_galfit <= 21 | dm_host_rest | i_hubble | 16 | 0.410 (0.115) | 0.468 (0.068) | [-0.07, 0.79] |
| n < 2 | dm_host_obs | q | 23 | -0.364 (0.088) | -0.309 (0.152) | [-0.65, 0.14] |
| n < 2 | dm_host_obs | i_hubble | 23 | 0.342 (0.110) | 0.309 (0.152) | [-0.14, 0.65] |
| n < 2 | dm_host_rest | q | 23 | -0.377 (0.076) | -0.334 (0.119) | [-0.67, 0.12] |
| n < 2 | dm_host_rest | i_hubble | 23 | 0.358 (0.094) | 0.334 (0.119) | [-0.12, 0.67] |
| no bulge fits | dm_host_obs | q | 26 | -0.389 (0.050) | -0.383 (0.053) | [-0.68, 0.03] |
| no bulge fits | dm_host_obs | i_hubble | 26 | 0.370 (0.063) | 0.383 (0.053) | [-0.03, 0.68] |
| no bulge fits | dm_host_rest | q | 26 | -0.393 (0.047) | -0.397 (0.045) | [-0.69, 0.01] |
| no bulge fits | dm_host_rest | i_hubble | 26 | 0.375 (0.059) | 0.397 (0.045) | [-0.01, 0.69] |

## Dropped low-$q$ hosts

| FRB | q | z | DM_host obs | DM_host rest | m | n |
|---|---|---|---|---|---|---|
| 20220105A | 0.0259 | 0.2785 | 280 | 358 | 20.96 | 1.05 |
| 20231230D | 0.1010 | 0.5060 | 110 | 166 | 21.01 | 1.90 |
| 20240201A | 0.1666 | 0.0427 | 260 | 271 | 17.01 | 0.82 |
| 20240203C | 0.1660 | 0.2426 | 170 | 211 | 22.12 | 1.40 |

## Hosts shared with our sample (N = 20)

Subscript M = Marnoch Table 1, O = ours (DM_host without halo term, so offset by a near-constant).

| q | DM_host rest | N | $r_S$ (p) |
|---|---|---|---|
| q_M | dm_host_rest_M | 20 | -0.194 (0.413) |
| q_M | dm_host_rest_O | 20 | -0.251 (0.286) |
| q_O | dm_host_rest_M | 20 | -0.089 (0.710) |
| q_O | dm_host_rest_O | 20 | -0.160 (0.500) |

Agreement: $q$ $r_S$ = 0.857; DM_host $r_S$ = 0.988.
Their hosts not in our sample: 20190611B, 20190711A, 20191228A, 20200430A, 20210117A, 20220105A, 20220918A, 20230708A, 20231230D, 20240117B, 20240201A, 20240203C, 20240310A, 20240615B.
