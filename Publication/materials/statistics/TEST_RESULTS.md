# Confirmed-50 vs HSC — statistical tests

Question: could the 51 confirmed FRB-host Hubble \(\cos(i)\) values be a random draw from the HSC disk pool (\(r\le 22\), \(b/a>0.2\), \(0.4<n<1.5\))?

Tests use **point estimates** (winning-leg \(q \to\) Hubble \(\cos(i)\), \(q_0=0.2\)). They do **not** fold in Protocol A/B error smearing. Same convention as [`run_sdss_frb_inclination_tests.py`](../../../scripts/run_sdss_frb_inclination_tests.py).

| | FRB | HSC |
|---|---:|---:|
| N | 51 | 24,450 |
| median \(\cos(i)\) | 0.610 | 0.575 |
| mean \(\cos(i)\) | 0.585 | 0.567 |
| FRB − HSC (median) | +0.035 | |
| FRB − HSC (mean) | +0.018 | |

## How to read a p-value

p is the chance of a discrepancy **at least this large** if FRB and HSC really share one continuous distribution. It is **not** a percent-similar score.

| p | Read |
|---|---|
| > 0.10 | Consistent with the same distribution (for that test). |
| 0.05–0.10 | Weak / suggestive only. |
| 0.01–0.05 | Significant at the usual 5% bar. |
| < 0.01 | Stronger evidence of a difference. |

Two p-values are listed:

- **p (vs pool)** — textbook two-sample test of the 51 vs all 24,450 HSC galaxies. Power is high: a *small* CDF offset can make p tiny.
- **p (matched N)** — among 10,000 random HSC draws of 51, the fraction whose statistic vs the pool is at least as large as the FRB’s. This matches the CDF envelope (same N). Prefer this when asking “is this 51 unusual among 51-galaxy HSC samples?”

SciPy AD p-values are **capped at 0.25** when samples look alike (true p is at least 0.25) and floored near 0.001 when they differ strongly.

**Takeaway.** Matched-N Monte Carlo p-values sit near 0.3-0.6 (typical HSC-50 draw). Location tests vs the pool agree. The CvM vs-pool p is a SciPy N2>>N1 artifact; use matched-N only.

## Summary

| Test | statistic | p (vs pool) | p (matched N) | notes |
|---|---:|---:|---:|---|
| Kolmogorov–Smirnov | 0.115 | 0.476 | 0.467 | D = max \|F−G\| |
| Anderson–Darling | -0.461 | >= 0.25 (SciPy cap) | 0.606 | tails weighted |
| Cramér–von Mises | 2.109 | unreliable | 0.476 | vs-pool p broken at N2 >> N1 |
| Kuiper | 0.169 | — | 0.443 | V = D+ + D−; p from MC only |
| Mann–Whitney U | 654199.0 | 0.543 | — | rank-biserial = -0.049 |
| Brunner–Munzel | -0.642 | 0.524 | — | location; unequal variance OK |
| Epps–Singleton | 3.43 | 0.489 | — | characteristic function |
| Welch t (means) | 0.601 | 0.551 | — | parametric; cos(i) is not Gaussian |
| Wasserstein-1 | 0.0325 | — | 0.522 | mean \|quantile gap\|; p from MC |

## Kolmogorov–Smirnov (two-sample)

**How it works.** Sort both samples and form empirical CDFs \(F_{50}\) and \(G_{\mathrm{HSC}}\). The statistic is \(D=\sup_x |F_{50}(x)-G_{\mathrm{HSC}}(x)|\), the largest vertical gap. SciPy’s two-sided p uses the asymptotic (or exact, when cheap) null distribution of \(D\) under a shared continuous law.

**How to read it.** \(D\) is a distance between CDFs on [0, 1]: here \(D=0.115\) means the curves never differ by more than that. p answers “would a gap this big show up often if both were draws from one law?” KS is most sensitive near the middle of the distribution and weaker in the tails than AD.

**Result.** \(D=0.115\), p (vs pool) = 0.476, p (matched N) = 0.467.

## Anderson–Darling (k-sample)

**How it works.** SciPy `anderson_ksamp` compares the two samples to the pooled ECDF with a weight \(1/[H(1-H)]\) that **boosts the tails**. The statistic is a standardized integral of the weighted gap. Null: one continuous distribution.

**How to read it.** Larger statistic → more evidence against a shared law, especially if the mismatch is at high/low \(\cos(i)\). Do **not** read the sign of the statistic (recent SciPy can return negative values when samples agree). Use the p-value. A printed p of 0.25 is a **cap** (true p ≥ 0.25).

**Result.** AD = -0.461, p (vs pool) = >= 0.25 (SciPy cap), p (matched N) = 0.606.

## Cramér–von Mises (two-sample)

**How it works.** Integrates \((F-G)^2\) over the pooled sample instead of taking only the single largest gap (KS). SciPy `cramervonmises_2samp` reports that integral and an approximate p.

**How to read it.** Middle ground between KS (one spike) and AD (tail emphasis): a modest gap over a wide stretch of \(\cos(i)\) can beat a tall narrow spike. Same p-scale as above.

SciPy vs-pool p is uncalibrated when N_HSC >> N_FRB. Use matched-N p.

**Result.** T=2.109, p (vs pool) = ignore (6.74e-06, uncalibrated), p (matched N) = 0.476.

## Kuiper (two-sample)

**How it works.** \(V=D^{+}+D^{-}\): the sum of the largest upward and downward ECDF gaps. Equal to KS \(D\) when the discrepancy is one-sided; larger when the curves cross. No SciPy two-sample p, so only the matched-N Monte Carlo p is quoted.

**How to read it.** Sensitive to **shape** differences that flip sign (excess edge-on *and* excess face-on). A large \(V\) with a small KS \(D\) means the CDFs cross. p (matched N) is the fraction of HSC-50 draws with \(V\) at least this large vs the pool.

**Result.** \(V=0.169\), p (matched N) = 0.443.

## Mann–Whitney U

**How it works.** Rank all \(N_1+N_2\) values. \(U\) counts how often an FRB \(\cos(i)\) beats an HSC one. Two-sided p tests \(P(\mathrm{FRB}>\mathrm{HSC})=1/2\) (no stochastic dominance / location shift). The rank-biserial correlation \(r=1-2U/(N_1 N_2)\) is the effect size on \([-1,1]\).

**How to read it.** p speaks only to a **shift**, not to variance or tail shape (AD/KS can fire when MWU does not). \(|r|\lesssim 0.1\) is a small location shift even if p is small against 24,450 HSC galaxies. Sign: \(r>0\) means FRB ranks **lower** \(\cos(i)\) (more edge-on) than HSC.

**Result.** \(U=654199.0\), p (vs pool) = 0.543, rank-biserial \(r=-0.049\) (median \(\Delta\cos(i)=+0.035\)).

## Brunner–Munzel

**How it works.** Same stochastic-dominance question as Mann–Whitney, but it does **not** assume equal variance (MWU’s rank null is touchy when the two spreads differ). SciPy `brunnermunzel`, two-sided.

**How to read it.** Treat it as a robustness check on MWU. If BM and MWU p-values agree, the location conclusion is not an equal-spread artifact. Same p-scale.

**Result.** \(W=-0.642\), p (vs pool) = 0.524.

## Epps–Singleton

**How it works.** Compares the empirical **characteristic functions** \(\langle e^{itX}\rangle\) of the two samples at a few \(t\) points (SciPy `epps_singleton_2samp`). Detects differences in location, scale, or shape that moment- or rank-only tests can miss. Needs \(N\gtrsim 5\) per sample; can fail on heavy ties.

**How to read it.** The statistic is a \(\chi^2\)-like distance between characteristic functions. Large stat / small p → the full laws differ, without saying *where*. Use the CDF plot for direction.

**Result.** \(W=3.43\), p (vs pool) = 0.489.

## Welch \(t\) (means)

**How it works.** Two-sample \(t\) on the means with unequal-variance standard errors. Assumes each sample mean is approximately normal (CLT is comfortable at \(N=51\), but \(\cos(i)\) is bounded).

**How to read it.** This is a **mean shift** test only. Do not prefer it over MWU/BM for the scientific claim; it is here because \(N=50\) is large enough that the mean is a stable summary. Mean \(\Delta\cos(i)=+0.018\).

**Result.** \(t=0.601\), p (vs pool) = 0.551.

## Wasserstein-1

**How it works.** Earth-mover distance between the two 1-D laws: \(\int |F^{-1}-G^{-1}|\,du\), the average absolute gap between quantile functions. Not a significance test by itself; p comes from the matched-N Monte Carlo.

**How to read it.** \(W_1\) is in \(\cos(i)\) units. \(W_1=0.0325\) is the typical horizontal shift you would need to morph one CDF into the other. Compare it to the plot’s 68% band width, not to a 0–1 “agreement” scale.

**Result.** \(W_1=0.0325\), p (matched N) = 0.522.

## Directional tests (fewer edge-on hosts)

Pre-specified direction: FRB hosts have fewer edge-on systems than the HSC disks, i.e. larger \(\cos(i)\). SciPy has no one-sided Anderson–Darling; AD stays the two-sided primary test above. One-sided KS and Mann–Whitney are the directional complements. Cramér–von Mises remains the two-sided optional test in the summary.

**One-sided KS** (`alternative='less'`): tests whether the FRB CDF lies below the HSC CDF. \(D^-=0.115\), p (vs pool) = 0.241, p (matched N) = 0.237.

**One-sided Mann–Whitney** (`alternative='greater'`): tests whether FRB \(\cos(i)\) is stochastically larger. \(U=654199.0\), p (vs pool) = 0.271, p (matched N) = 0.274.

### \(\chi^2\) in HSC quantile bins

Bins are equal-probability intervals of the HSC pool, so each bin expects \(N/n_{\mathrm{bins}}\) FRB hosts if the samples share one law. With \(N=51\), every split from 2 to 5 bins has expected counts above 5. Quartiles match the 2024 paper and are the adopted split (expected 12.8 per bin).

| bins | expected / bin | FRB counts (low \(\cos i\) → high) | \(\chi^2\) | p |
|---:|---:|---|---:|---:|
| 2 | 25.5 | 21,30 | 1.59 | 0.208 |
| 3 | 17.0 | 13,21,17 | 1.88 | 0.390 |
| 4 adopted | 12.8 | 11,10,19,11 | 4.14 | 0.247 |
| 5 | 10.2 | 8,8,12,15,8 | 4.00 | 0.406 |

Quartile matched-N p = 0.254. Counts (edge-on to face-on): 11,10,19,11.

### Binomial edge-on fraction

Edge-on means \(i\) above a fixed cut, equivalently \(\cos i < \cos i_{\mathrm{cut}}\). The null probability is the HSC fraction above that cut. The alternative is a smaller FRB fraction. \(i>60^\circ\) is the quoted cut; \(45^\circ\) and \(75^\circ\) are checks.

| cut | FRB \(k/N\) | FRB fraction | HSC fraction | one-sided p |
|---:|---:|---:|---:|---:|
| \(i>45^\circ\) | 34/51 | 0.667 | 0.681 | 0.465 |
| \(i>60^\circ\) | 15/51 | 0.294 | 0.398 | 0.083 |
| \(i>75^\circ\) | 4/51 | 0.078 | 0.102 | 0.396 |

For \(i>60^\circ\), matched-N p = 0.081 (fraction of HSC draws of 51 with an edge-on count at least as small).

## Caveats

1. Point estimates only — Protocol A/B error bands are not in these p-values.
2. HSC $q$ is from the Kawinwanichakij et al. Lenstronomy single-Sérsic fits; FRB $q$ is this pipeline’s GALFIT winning leg.
3. Pool tests have high power at \(N_{\mathrm{HSC}}=24{,}450\); lean on matched-N p and the CDF when judging *practical* difference.
4. Several tests share one sample; do not treat six small p-values as six independent discoveries.
