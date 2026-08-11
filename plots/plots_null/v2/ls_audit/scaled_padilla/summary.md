# scaled_padilla

## Fit

| Parameter | Value |
|-----------|------:|
| N (fit pool) | 2,000,000 |
| chi2 | 1435307.24 |
| mu_g (mean thickness C/A) | 0.2645 |
| sig_g | 0.0250 |
| mu (mean ln eps) | -0.5214 |
| sig | 0.3036 |
| median eps = exp(mu) | 0.5937 |
| E0 (Padilla dust) | 1.5971 |

Literature seeds: Ryden (mu_g~0.22, ln eps~-1.85); Padilla (ln e~-2.33, E0~0.45).

## Mag-cut model cos(i) CDFs (strict b/a > q0=0.2)

| mag_limit | N | median cos(i) model | median cos(i) Hubble |
|----------:|--:|--------------------:|---------------------:|
| 20 | 25,274 | 0.6576 | 0.4964 |
| 21 | 138,605 | 0.6675 | 0.4554 |
| 22 | 524,666 | 0.6707 | 0.4253 |

Joint Padilla shape+E0 fit. LS has no redshifts; psi(theta) uses 10^{-0.4 E(theta)} only (not full 1/Vmax+LF). CDFs use observed mag cuts, strict b/a>q0. Per-galaxy cos(i) is DRAWN from the model posterior P(cos i | b/a) (one sample each), not the per-bin median.

Large E0 and large face-on ellipticity partly compensate for the EXP/REX selection (deficit of round systems) in this photometric adaptation.

