# Confirmed-50 vs HSC MC cos(i) CDF

Paper sample: `in_53` intersect `confirmed` = **50**. Winning-leg q (13 re-fit, 37 production). Hubble q0 = 0.2.

HSC pool: Kawinwanichakij EXP analogue (`goodfits=1`, 0.4 < n < 1.5), photometric r <= 22, b/a > 0.2. **N = 24,450**. Each of 10,000 draws samples 50 galaxies without replacement.

| | |
|---|---|
| median q (FRB) | 0.630 |
| median cos(i) (point) | 0.610 |
| median sigma_q,A (10 x q_err) | 0.100 |
| median sigma_q,B (sky xor q_err) | 0.011 |
| hosts missing sky +/- | 0 (those get sigma_q,sky = 0) |

## Protocol A — `frb_vs_hsc_inflate10.png`

sigma_q = 10 * q_err, then i = Hubble(q) + N(0, 5 deg).

## Protocol B — `frb_vs_hsc_sky_quad.png`

sigma_q,sky = |q+ - q-| / 2
sigma_q = sqrt(q_err^2 + sigma_q,sky^2)
then i = Hubble(q) + N(0, 5 deg).

The 5 deg term is an inclination floor (Barnes & Sellwood 2003),
not a GALFIT q error. q_err is not multiplied by 10.

Table: `confirmed50_q.csv`. Tests: [`TEST_RESULTS.md`](TEST_RESULTS.md). Tracking note (50, Protocol B errors, and the 2024 hosts): [`STATISTICAL_FINDINGS.md`](STATISTICAL_FINDINGS.md).

Literature on GALFIT **axis-ratio** errors (not \(n\)/\(R_e\)/mag): [`GALFIT_Q_ERROR_LITERATURE.md`](GALFIT_Q_ERROR_LITERATURE.md).
