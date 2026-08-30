# GALFIT axis-ratio errors: literature census

Filter: a paper is included only if it treats uncertainties on the projected axis ratio \(q=b/a\), ellipticity \(\varepsilon=1-q\), or an equivalent shape parameter from GALFIT (or a GALFIT-class Levenberg–Marquardt Sérsic fit). Papers that discuss GALFIT errors only for \(n\), \(R_e\), or magnitude are listed in §5 and are **not** used to defend a \(q\) policy.

Context: confirmed-50 CDFs currently have two error protocols ([README.md](README.md)).

| | Protocol A | Protocol B |
|---|---|---|
| \(\sigma_q\) | \(10\times q_{\rm err}\) | \(\sqrt{q_{\rm err}^2+\sigma_{q,{\rm sky}}^2+\sigma_{q,5^\circ}^2}\) |
| extra term | \(i\leftarrow i+\mathcal{N}(0,5^\circ)\) after Hubble | \(5^\circ\) mapped into \(q\) *before* Hubble; no second floor |
| median \(\sigma_q\) (N=50) | 0.100 | 0.067 |

\(q_{\rm err}\) is GALFIT’s formal \(1\sigma\) from the covariance matrix. \(\sigma_{q,{\rm sky}}=|q_+-q_-|/2\) from \(\pm1\sigma\) sky re-fits. Verification §4.4 uses \(\Delta q_{\rm sky}=\max(|q_+-q_0|,|q_--q_0|)\), which is a different estimator of the same experiment.

---

## 1. What GALFIT itself claims about \(q\) errors

### Peng et al. (2002), AJ 124, 266 — GALFIT v1

Formal errors come from the Hessian / covariance of \(\chi^2\) at the minimum. The \(68\%\) interval on each free parameter (including \(q\)) is the \(\Delta\chi^2=1\) shell of that ellipsoid. They verified the recipe on artificial images with Poisson noise for one- and two-component models. That test applies to **all** free parameters, so it includes axis ratio — but only in the Poisson, perfect-model limit.

### Peng et al. (2010), AJ 139, 2097 — GALFIT v3

The covariance matrix is a **lower bound**. Quote (their §VI.5): when the model does not fit the data to Poisson residuals, “uncertainties inferred from covariance matrices [are] underestimated.” Realistic errors “are necessarily obtained by other processes, such as comparing fit results based on different assumptions.” They name profile mismatch, neighbours, and flat-field / sky as the dominant residuals. Axis ratio \(q\) is degenerate with a low-amplitude Fourier \(m=2\) mode; using both together is discouraged.

This is the strongest author-level statement: **do not quote `fit.log` \(\sigma_q\) as the science error**. Compare fits under different assumptions.

### GALFIT Technical FAQ (Peng)

Formal `fit.log` uncertainties are Poisson-only. Start-guess tolerance on axis ratio is \(\sim 0.5\) (that is a convergence window, not an error bar). Sky held fixed to a wrong value is listed as a leading cause of systematic parameter error.

---

## 2. Papers that actually treat \(\sigma_q\) / \(\sigma_\varepsilon\) / \(\sigma_{b/a}\)

### 2.1 Simulation scatter replacing the covariance (the GEMS method)

**Häussler et al. (2007), ApJS 172, 615.** Single-Sérsic GALFIT vs GIM2D on \(\sim 80{,}000\) simulated ACS galaxies. Input \(b/a\) is drawn from a Hubble disk with \(q_0=0.18\) (disks) or \(0.45\) (spheroids).

They **discard GALFIT/GIM2D formal errors for the catalog**. Replacement: at the object’s \(\mu\) and \(n\), take the simulation residual scatter, linearly interpolate between the \(n=1\) and \(n=4\) grids, and apply a floor \(\sigma_{b/a}=0.001\). That recipe is applied to magnitude, \(R_e\), \(n\), **\(b/a\)**, and PA.

Fig. 14 (formal \(\sigma/\Delta\) vs truth) is shown for magnitude, \(n\), and \(R_e\) only. Caption: magnitude and \(n\) “dramatically underestimated”; \(R_e\) “significantly better estimated.” They do **not** plot \(b/a\) on that figure. The catalog method nevertheless treats \(b/a\) the same way as the other four parameters.

Sky-estimator test on **\(n=4\)** simulations, recovered/input ellipticity (“\(e\) ratio”):

| sample | sky | \(n\) | \(r_{50}\) ratio | \(\Delta\)mag | **\(e\) ratio** |
|---|---|---|---|---|---|
| bright (\(\mu_{\rm in}<22.5\), \(m<22.5\)) | isophotal | \(3.99\pm0.27\) | \(1.00\pm0.05\) | \(0.00\pm0.03\) | **\(0.99\pm0.04\)** |
| | SExtractor | \(3.79\pm0.29\) | \(0.96\pm0.07\) | \(0.03\pm0.04\) | **\(1.00\pm0.04\)** |
| | GALFIT free | \(3.94\pm0.24\) | \(0.99\pm0.06\) | \(0.01\pm0.03\) | **\(1.00\pm0.04\)** |
| faint (\(23.5<\mu_{\rm in}<26\)) | isophotal | \(3.95\pm1.13\) | \(0.97\pm0.44\) | \(0.02\pm0.42\) | **\(1.04\pm0.18\)** |
| | SExtractor | \(3.05\pm0.98\) | \(0.61\pm0.27\) | \(0.41\pm0.39\) | **\(1.05\pm0.18\)** |
| | GALFIT free | \(3.78\pm1.16\) | \(0.94\pm0.47\) | \(0.09\pm0.44\) | **\(1.04\pm0.19\)** |

Wrong sky wrecks \(n\), \(R_e\), and mag. **It does not shift \(q\)**. Even at faint \(\mu\), the \(e\)-ratio mean stays \(\approx1.04\)–\(1.05\) across all three sky methods; only the scatter grows (\(\sigma\approx0.18\)). For bright galaxies the \(q\) scatter is \(0.04\).

Häussler’s prose conclusion is that exponential disks are “well fit and have small measurement errors.” That is the relevant regime for our \(n\sim1\), \(m\le22\) hosts.

**How they treat \(q\):** simulation residual scatter vs \(\mu\), interpolated in \(n\); formal GALFIT \(\sigma_q\) unused; floor \(0.001\).

### 2.2 Repeated real-data measurements (CANDELS)

**van der Wel et al. (2012), ApJS 203, 24.** GALFIT + GALAPAGOS on CANDELS. Random errors from **deep vs shallow** reductions of the same 6,492 GOODS-S objects, not from `fit.log`. Uncertainties are the half-width of the \(16\)–\(84\) percentile range of those differences, nearest-neighbour interpolated in \((m,n,R_e)\) and scaled as \(1/(S/N)\). **\(\delta q\) is one of the five reported parameters** (their Eq. 2). They explicitly state they see **no correlation between \(q\) (or PA) and the uncertainties in the other parameters** after \(S/N\) normalisation.

Bottom line (their Table 3 / §V.1): \(m\), \(R_e\), and **\(q\)** reach random accuracy of **\(20\%\) or better** for \(H_{\rm F160W}\lesssim24.5\); \(n\) only to \(H\lesssim23.5\). Typical faint high-\(z\) galaxies (small, low \(n\)) reach **\(\sim10\%\)** on size and shape to \(H\sim24.5\). Systematic offsets from simulations are smaller than random errors over most of the science range.

Sky contribution: they re-ran GALFIT on 1,000 objects with scrambled GALAPAGOS backgrounds. Background is **at most \(25\)–\(30\%\) of the total error budget**, except for \(H>25.5\) and \(R_e>0.4''\), where it dominates. An object at \(H\sim22\) has \(S/N\sim100\)–\(200\) in CANDELS-wide — our \(m\le22\) hosts sit in that high-\(S/N\) regime.

**How they treat \(q\):** empirical \(\delta q(S/N,m,n,R_e)\) from deep–shallow repeats; formal GALFIT unused; sky is a minority term.

**Davari et al. (2016), ApJ 787, 69.** Compact high-\(z\) GALFIT simulations. For \(23<H_{160}<24\): uncertainties on \(R_e\), \(m_H\), \(n\), and **ellipticity \(e\)** are \(<20\%\), \(0.2\) mag, \(0.2\), and **\(10\%\)**, respectively. Brighter galaxies are better. They cite van der Wel (2012) Table 3 as the parent result.

**How they treat \(q\):** simulation residual on \(e\); quote \(10\%\) at \(H\sim23\)–\(24\).

### 2.3 Simulations that bin in ellipticity but do not publish \(\sigma_\varepsilon\)

**Hoyos et al. (2011), MNRAS 411, 2439.** Coma ACS, GALFIT and GIM2D. \(2\times10^5\) GALFIT simulations, **binned in ellipticity** (2 bins, \(0\le\varepsilon\le0.8\)) as well as mag, \(\log R_e\), and \(n\). They state GALFIT “tends to underestimate errors” (citing Häussler). The **published error formula** (their Table 3 / Eq. 5) is only for magnitude, \(R_e\), and \(n\). The catalog column \(\sigma_{\rm Ellip}\) (GF) is the **formal GALFIT ellipticity error**.

**How they treat \(q\):** ellipticity is a fitted parameter and a simulation bin; the science error recipe they ship is **not** for ellipticity. Formal \(\sigma_\varepsilon\) is what appears in the table.

### 2.4 Sky \(\pm1\sigma\) re-fits (the Huang / Gao method) — \(\varepsilon\) included for bulges

**Huang et al. (2013), ApJ 766, 47** (CGS III), as used by **Gao & Ho (2017), ApJ 845, 114**. Formal GALFIT errors “do not properly capture” sky systematics. Empirical fix: **re-fit with sky held at best-fit \(\pm1\sigma\)** and take the variation of the **bulge parameters**. Gao & Ho list those bulge parameters as \(m\), \(\mu_e\), \(n\), \(r_e\), **and apparent ellipticity \(\varepsilon\)**. Additional model-mismatch error from the range of still-plausible multi-component models. For 1-D fits they add sky and “excluded-range” terms **in quadrature**.

This is the literature source of Protocol B’s sky term. Caveat: Gao & Ho apply it to **bulge** \(\varepsilon\), not disk \(q\). Bulge \(n\) and \(R_e\) couple to sky through the high-\(n\) wings; disk \(q\) does not (Häussler Tables 4–5). Using the same experiment on host \(q\) is still valid as a systematic, but the expected \(\sigma_{q,{\rm sky}}\) is small.

**How they treat \(q\):** sky \(\pm1\sigma\) re-fit \(\rightarrow\Delta\varepsilon\) (bulge); plus model-to-model range; 1-D: quadrature of sky and analysis choices. Formal covariance not used as the science error.

### 2.5 Formal errors published with an explicit “underestimated” flag

**Almaini et al. (2017) UDS DR11 GALA catalog** (GALAPAGOS-2 + GALFIT3). Columns `Q_GALFIT_BAND` and `QERR_GALFIT_BAND`. README: “errors on the GALFIT structural parameters have been provided, these are known to be underestimated (Häussler et al. 2007), so should be used with caution.” No \(q\)-specific inflation factor.

### 2.6 Bayesian / MCMC posteriors on \(q\) (not GALFIT, but the HSC comparison catalog)

**Kawinwanichakij et al.** HSC PDR2 structural catalog (the parent of our EXP analogue). `fitted_q` from lenstronomy; `fitted_x_err` from **MCMC**. Bias corrections (`corrected_x`) from image simulations. This is a posterior width on \(q\), not a GALFIT Hessian. Our HSC Monte Carlo currently uses point \(q\) with no error smear — that is **more conservative against finding a difference** than smearing HSC too.

### 2.7 Method comparison as a systematic on \(i\), not a GALFIT \(\sigma_q\)

**Bhardwaj et al. (2024), Nature 634, 1065.** Not GALFIT. AutoProf Sérsic (or Superellipse/Gaussian/Spline for 3 hosts) and Photutils isophotes \(\rightarrow\) Hubble \(q_0=0.2\). Covariance on \(i\) is \(\approx1^\circ\) from both packages. **RMSD between the two methods is \(3.7^\circ\)**, which they treat as a systematic **not captured by the covariance**. They do **not** multiply covariance errors by 10. Science conclusions are unchanged if Photutils \(i\) is used instead of AutoProf.

This is the closest published number to Protocol A/B’s \(5^\circ\) floor for *this* science case.

---

## 3. The \(5^\circ\) term: inclination, not GALFIT

No GALFIT paper assigns a \(5^\circ\) floor to \(q\). That number lives in the Hubble-formula / Tully–Fisher literature:

- Hubble (1926) formula with finite thickness \(q_0\). We and Bhardwaj (2024) use \(q_0=0.2\) (Vallejo, Unterborn & Ryden 2008).
- Photometric \(i\) is biased by bars, arms, dust, and a non-universal \(q_0\) (Bhardwaj Methods). Comparing FRB hosts to a survey catalog measured the same way cancels a large part of that bias — which is why both Bhardwaj and our HSC test compare CDFs rather than to \(\mathrm{CDF}=\cos i\).
- Tully–Fisher practice commonly adopts \(\sigma_i\sim5^\circ\) as a typical photometric inclination uncertainty (Tully & Fouqué 1985; Giovanelli et al. 1997; later TF samples). That is a floor on **\(i\)**, not a multiplier on \(\sigma_q\).
- Bhardwaj’s empirical method-mismatch is \(3.7^\circ\), not \(5^\circ\).

Protocol A adds \(\mathcal{N}(0,5^\circ)\) **after** Hubble. Protocol B converts \(\pm5^\circ\) into \(\sigma_q\) **before** Hubble. Adding both is double-counting the same geometric systematic.

---

## 4. The \(10\times q_{\rm err}\) factor: no \(q\)-specific literature support

I did not find a paper that multiplies GALFIT \(\sigma_q\) (or \(\sigma_\varepsilon\), \(\sigma_{b/a}\)) by 10.

What exists:

- Häussler (2007) Fig. 14: formal errors on **magnitude and \(n\)** are low by “a large factor”; **\(R_e\) is better**. \(q\) is not in the figure.
- Peng (2010): covariance is a lower bound for **all** parameters when the model is imperfect. No factor given.
- Survey catalogs (UDS GALA, and the usual GALAPAGOS READMEs): “use GALFIT errors with caution (Häussler 2007).” No \(\times10\).
- van der Wel / Davari: replace the covariance with an empirical \(\delta q\sim10\)–\(20\%\) **relative** error at faint \(H\). For a typical host \(q=0.63\), \(20\%\) is \(\sigma_q=0.13\); at \(m\le22\) and \(n\sim1\), their own scaling says the random error is much smaller than that faint-end 20%.

Protocol A’s median \(\sigma_{q,A}=0.100\) is numerically close to the **faint, high-\(z\), mixed-\(n\)** 10–20% numbers, applied indiscriminately to **nearby, mag\(\le22\), exponential** hosts whose Häussler scatter on \(q\) is \(0.04\) (bright) and whose van der Wel \(S/N\) is in the \(100\)–\(200\) regime. The \(\times10\) is a transfer of the mag/\(n\) result onto \(q\).

Some hosts have \(10\times q_{\rm err}\sim1\), which is why Protocol A smears toward \(i=90^\circ\). No paper endorses drawing \(q\sim\mathcal{N}(q,1)\).

---

## 5. Papers checked and **not** used for \(\sigma_q\)

These discuss GALFIT (or equivalent) uncertainties but the published error analysis is \(n\), \(R_e\), mag, and/or \(B/T\) — not axis ratio. Citing them in defence of a \(q\) policy would be incorrect.

| Paper | What they actually error-budget |
|---|---|
| Hoyos et al. (2011) *error formula* | mag, \(R_e\), \(n\) only (ellipticity is a bin, not a delivered \(\sigma\)) |
| Meert, Vikram & Bernardi (2013) PyMorph simulations, as summarised in Meert et al. (2015) | mag, size, \(n\), \(B/T\), sky; total axis ratio is measured from the model image but the quoted 5% accuracy is for mag and \(R_e\) |
| Kelvin et al. (2012) GAMA SIGMA | GALFIT single-Sérsic catalog includes \(q\); the Häussler caveat is stated generally; I did not find a \(q\)-specific scatter recipe in the error discussion |
| Salo et al. (2015) S4G | multi-component GALFIT; science product includes disk \(q\); uncertainties are not a published \(\sigma_q\) calibration |
| Simard et al. (2011) GIM2D SDSS | Metropolis posteriors including ellipticity — **not GALFIT**; Häussler showed GIM2D still underestimates |
| Lackner & Gunn (2012) | SDSS bulge+disk; not used here as a GALFIT-\(q\) calibration |
| Häussler Fig. 14 *as plotted* | mag, \(n\), \(R_e\) only |
| Guo et al. (2009), Bruce et al. (2012) | background \(\rightarrow\) size / \(n\) / mag (van der Wel cites them for sky, then measures \(\delta q\) himself) |
| Odewahn et al. (1997) | **isophotal** \(b/a\) vs \(S/N\), not GALFIT; faint images measured too round — a CDF systematic, not a GALFIT \(\sigma_q\) |

---

## 6. Approaches, collapsed

Every paper that actually treats \(q\) uses one (or a combination) of these. None uses \(\times10\) on formal \(\sigma_q\).

| Approach | Who | What it captures | Misses |
|---|---|---|---|
| Hessian / `fit.log` | Peng 2002 (Poisson tests); Hoyos catalog \(\sigma_{\rm Ellip}\); UDS GALA with warning | photon noise, parameter covariance | sky, neighbours, model mismatch |
| Simulation residual scatter vs \(\mu,n\) | Häussler 2007 (includes \(b/a\)); Davari 2016 (\(e\)) | Poisson + crowding + sky structure in the mock | real non-Sérsic structure (Häussler and van der Wel both say this underestimates) |
| Repeat measurements (deep/shallow or independent reductions) | van der Wel 2012 (\(\delta q\)) | noise + sky + reduction; real galaxy structure | model choice (still single Sérsic) |
| Sky \(\pm1\sigma\) re-fit | Huang 2013; Gao & Ho 2017 (\(\varepsilon_{\rm bulge}\)); our §4.4 / Protocol B | sky systematic on that object | Poisson, model mismatch, Hubble \(q_0\) |
| Compare different models / codes | Peng 2010 (recommended); Gao & Ho model range; Bhardwaj AutoProf vs Photutils (\(3.7^\circ\)) | model mismatch | needs more than one fit |
| MCMC / bootstrap posterior | Kawinwanichakij HSC; IMFIT (Erwin 2015) and ProFit (Robotham et al. 2017) recommend this because LM covariance is low | full posterior under one model | still conditional on the model |
| Floor on \(i\) (\(\sim5^\circ\)) | TF tradition; Bhardwaj \(3.7^\circ\) RMSD | Hubble \(q_0\), bars, dust, method mismatch | not a GALFIT error; double-counts if stacked on a large \(\sigma_q\) |
| Quadrature of independent terms | Gao & Ho (1-D: sky \(\oplus\) excluded-range) | combines named systematics | only as good as the term list |

Peng (2010) is explicit that the right operation for galaxy fitting is **the last two rows**, not a scalar times the Hessian.

---

## 7. What this means for Protocols A and B

**Protocol A (\(10\times q_{\rm err}\) then \(+5^\circ\)).**

- Defensible only as a blunt reading of “Häussler: formal errors are too small” plus a TF \(i\) floor.
- Not defensible as a \(q\) policy: Häussler’s large-factor result is for mag and \(n\); their own \(e\)-ratio scatter at bright \(\mu\) is \(0.04\), and \(q\) does not move when sky is wrong. van der Wel’s \(20\%\) is a faint-\(H\) relative error, not a licence to draw \(q\sim\mathcal{N}(q,1)\). Bhardwaj (2024) used \(\sim1^\circ\) covariance plus \(3.7^\circ\) method RMSD — the opposite of \(\times10\).
- Adding \(5^\circ\) *after* a \(\sigma_q\) that is already \(\sim0.1\) double-smears edge-on.

**Protocol B (formal \(\oplus\) sky half-range \(\oplus\) \(5^\circ\) in \(q\)).**

- Matches Peng (2010): combine named assumptions (sky, Hubble geometry) rather than trust the Hessian.
- Sky term is exactly Huang / Gao & Ho’s experiment, applied to host \(q\). Häussler Tables 4–5 say that term should be **small** for disks, which is what we want the data to decide per host (`flag_sky_sensitive` is \(\Delta q_{\rm sky}>0.05\)).
- Formal \(q_{\rm err}\) is kept as the Poisson piece (Peng 2002, valid when \(S/N\) is high and the model is not terrible). It is **not** multiplied.
- \(5^\circ\) is mapped into \(q\) so it is not applied twice. \(3.7^\circ\) (Bhardwaj) or \(5^\circ\) (TF) are the published scales for that term.
- van der Wel: at our \(S/N\), sky should be \(\lesssim30\%\) of the budget. Quadrature with a non-inflated \(q_{\rm err}\) is the right shape; \(\times10\) is not.

**Neither protocol** does the two things the large-survey papers actually did for \(q\): (i) replace \(\sigma_q\) with simulation or deep–shallow scatter at the object’s \(\mu,n,S/N\), or (ii) MCMC. For N=50 nearby disks, (i) would give \(\sigma_q\sim0.04\) (Häussler bright) to \(\lesssim0.10\) (van der Wel 10–20% relative at much fainter \(H\)). Protocol B’s median \(0.067\) sits in that band. Protocol A’s median \(0.100\) is at the faint-end edge, with a heavy tail that those papers do not justify.

---

## 8. Ingredients a literature-faithful policy would use

Not a new implementation — the pieces the papers actually combine:

1. **Keep formal \(q_{\rm err}\)** as the Poisson/covariance term. Do not \(\times10\). (Peng 2002; van der Wel’s high-\(S/N\) limit; Häussler \(e\)-ratio at bright \(\mu\).)
2. **Keep sky \(\pm1\sigma\) re-fits**, object by object. Prefer a single estimator and stick to it: either Protocol B’s \(|q_+-q_-|/2\) or verification’s \(\max(|q_\pm-q_0|)\). (Huang 2013; Gao & Ho 2017; Peng 2010 “different assumptions.”)
3. **Keep a Hubble-geometry / method term of a few degrees**, applied **once**, in \(q\) or in \(i\), not both. Calibrate to \(3.7^\circ\) (Bhardwaj same-sample RMSD) or \(5^\circ\) (TF convention). (Bhardwaj 2024; Tully–Fisher.)
4. **Optional, if we want Häussler/van der Wel literally:** a floor \(\sigma_q\gtrsim0.02\)–\(0.04\) so that hosts with \(q_{\rm err}\sim0.002\) are not treated as infinitely precise. That floor is the bright-disk simulation scatter, not \(10\times\) the Hessian.
5. **Do not smear the HSC catalog** unless we also smear with a published HSC \(q\) error (Kawinwanichakij MCMC). Point-estimate HSC vs smeared FRB is conservative against claiming a difference.

Protocol B is already the closest published analogue. Protocol A is the mag/\(n\) Häussler result applied to the wrong parameter.
