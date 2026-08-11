# scaled_ryden

## Fit

| Parameter | Value |
|-----------|------:|
| N (fit pool) | 2,000,000 |
| chi2 | 1968122.26 |
| mu_g (mean thickness C/A) | 0.2519 |
| sig_g | 0.0574 |
| mu (mean ln eps) | -0.5404 |
| sig | 0.3095 |
| median eps = exp(mu) | 0.5825 |
| E0 (Padilla dust) | 0.0000 |

Literature seeds: Ryden (mu_g~0.22, ln eps~-1.85); Padilla (ln e~-2.33, E0~0.45).

## Mag-cut model cos(i) CDFs (strict b/a > q0=0.2)

| mag_limit | N | median cos(i) model | median cos(i) Hubble |
|----------:|--:|--------------------:|---------------------:|
| 20 | 35,130 | 0.5053 | 0.424 |
| 21 | 185,477 | 0.5122 | 0.4047 |
| 22 | 641,099 | 0.5207 | 0.4008 |

Shape fit is Ryden photometry-only. Dust handled separately via Unterborn face-on mag re-cut for the CDFs in `cdfs/`. Per-galaxy cos(i) is DRAWN from the model posterior P(cos i | b/a) (one sample each), not the per-bin median: a given b/a maps to a distribution of inclinations, so using the median collapses every galaxy onto a narrow band and produces a degenerate near-vertical CDF.

Note: best-fit face-on ellipticity is much larger than Ryden's SDSS value (median eps~0.58 vs ~0.16). That is expected for Tractor EXP-only: REX removes near-round disks, so the apparent q distribution requires strong intrinsic ellipticity (and/or selection) to suppress q~1.

