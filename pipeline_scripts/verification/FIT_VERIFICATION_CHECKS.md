# GALFIT fit verification — method, implementation and results

**Created 2026-08-05.** A complete walkthrough of every diagnostic used to verify
the production Sérsic fits: why each test exists, what it is physically
measuring, how it is implemented, what can make it lie, and what it returned on
the first full run.

> **Equations will not render in Cursor's native Preview tab** (known Cursor
> bug — KaTeX is not wired there). Use either:
>
> - **`FIT_VERIFICATION_CHECKS.html`** — open in a browser (or regenerate with
>   `python render_checks_html.py`). This is the reliable viewer.
> - Classic preview: Command Palette → **Markdown: Open Preview to the Side**
>   (`Ctrl+Shift+V`). Requires `markdown.math.enabled: true` (now set in your
>   user settings).
>
> Math delimiters in this file are `$…$` / `$$…$$` only.

**Contents**

| § | |
|---|---|
| [1](#1-scope-and-ground-rules) | Scope, cohort, ground rules |
| [2](#2-running-the-suite) | Running the suite |
| [3](#3-the-shared-layer) | The shared layer: geometry, masks, model rebuild |
| [4](#4-the-checks) | The nine checks, in full |
| [5](#5-aggregation-flags-and-trust-tiers) | Aggregation, flags, trust tiers |
| [6](#6-data-walkthrough) | Data walkthrough — every file the suite writes |
| [7](#7-results--first-full-run) | Results of the first full run |
| [8](#8-deferred) | Deferred checks |
| [9](#9-open-questions) | Open questions |

---

## 1. Scope and ground rules

**What is being verified.** The pipeline reports an axis ratio $q = b/a$ per host
and converts it to an inclination. Everything downstream — the inclination
distribution that is the science result — rests on $q$ being right. These checks
exist to attack $q$ from as many independent directions as possible and see
whether it survives.

**Cohort.** All **64** hosts in
[`pipeline_galfit_results.csv`](../../pipeline_galfit_results.csv). Each row is
tagged `in_53` — magnitude $\le 22$ **and** $b/a > 0.2$ — so the 53-host science
cut is applied at *analysis* time, never at run time. The 11 hosts outside it are
still measured, because 9 of them are excluded (at least in part) *by the very
quantity under study*; see [§9](#9-open-questions).

**Read-only.** The suite never writes to `pipeline_scripts/Output/`. The one
check that re-runs GALFIT (§4.4) stages copies of the Phase 3a inputs into
`outputs/per_host/<FRB>/galfit_sky_*/`. Production `fit.log` and `out.fits` are
never touched.

**Baseline provenance.** Host geometry is read from the `out.fits` model-HDU
header, which is guaranteed consistent with the model and residual planes in the
same file. It is independently cross-checked against
[`scripts/galfit_fitlog_parse.py`](../../scripts/galfit_fitlog_parse.py) with
`sersic_component_index=0` — the production policy — and any disagreement above
0.02 in any parameter is recorded as a note. Every number therefore reconciles
with `pipeline_galfit_results.csv`.

**Design rule.** Each check is independent and idempotent, writes its own JSON,
and never raises: a failure in one check leaves the other eight intact. That is
what makes the suite re-runnable check-by-check as thresholds get calibrated.

---

## 2. Running the suite

```bash
cd pipeline_scripts/verification

python run_verification.py --checks all --jobs 4          # everything, 64 hosts
python run_verification.py --checks fourier isophote      # a subset
python run_verification.py --frb 20240210A --checks all --force
python run_verification.py --aggregate-only               # rebuild tables only
python -m pytest tests/test_fourier_recovery.py -q        # estimator unit tests
```

| Flag | Effect |
|---|---|
| `--checks` | `all` or any of `chi2 rff fourier psf mag isophote sky astrophot visual`. Order is fixed internally; `visual` always runs last because it consumes the others' products. |
| `--frb` | Restrict to named hosts. Names not in the cohort are warned about, not silently dropped. |
| `--cohort` | `all64` (default) or `53`. |
| `--jobs` | Process-pool width. Checks are CPU-bound and independent per host. |
| `--force` | Recompute even when a result JSON exists. Without it, existing results are reported as `cached`. |
| `--no-aggregate` / `--aggregate-only` | Skip or isolate the table-building step. |

**Runtime** (64 hosts, `--jobs 4`): the image-space checks are milliseconds;
`fourier` ~4 s/host (it rebuilds and convolves the model five times to pick the
sampling, then twice more for the response kernels); `astrophot` ~6 s/host;
`sky` ~4 s/host (two GALFIT runs); `isophote` is the slowest at ~28 s/host
because of the strategy search in §4.5. Full suite: about 25 minutes.

**Related runbooks** (not part of this check bible): index
[`VERIFICATION_README.md`](VERIFICATION_README.md); re-fits
[`REFIT_AND_REJECT_GRID.md`](REFIT_AND_REJECT_GRID.md); sky protocol
[`SKY_PROTOCOL.md`](SKY_PROTOCOL.md); sandbox [`SANDBOX.md`](SANDBOX.md);
confirmation [`HOST_CONFIRMATION_WORKFLOW.md`](HOST_CONFIRMATION_WORKFLOW.md).

---

## 3. The shared layer

Shared code lives in `vercommon.py`.

Everything a check needs about a host is assembled once by `load_host(frb)` and
passed around as a `HostData`.

### 3.1 What `load_host` assembles

| From | What |
|---|---|
| `out.fits` HDU 1/2/3 | data, **PSF-convolved model**, residual |
| `out.fits` HDU 2 header | final $x_c, y_c, q, \mathrm{PA}, R_e, n, m$ with GALFIT's errors, sky level, $\chi^2/\nu$, zero point |
| `host_sigma.fits` | the per-pixel $\sigma$ GALFIT weighted with |
| `host_mask.fits` | pixels GALFIT ignored |
| `proto_image.fits` | the PSF stamp GALFIT convolved with |
| `psfex.xml` | PSF FWHM, Stokes ellipticities, star count |
| `galfit.feedme` | plate scale, zero point, fit region |
| `cutout_meta.json`, `pipeline_summary.json`, `sky_fit_audit.json` | ROI policy, phase status, sky provenance |
| `fit.log` | independent parse for the baseline cross-check |

Two provenance assertions run on load. **Residual closure**: the suite verifies
$\lvert \mathrm{data} - \mathrm{model} - \mathrm{residual}\rvert \approx 0$
across the stamp and stores the maximum as `residual_closure`; a nonzero value
would mean the three planes are not from the same fit. **Coordinate convention**:
GALFIT reports 1-based FITS pixel coordinates, so centres are converted to
0-based NumPy indices on load. Getting this wrong shifts every elliptical
aperture by one pixel, which for a barely-resolved host is a large error.

### 3.2 Elliptical coordinates

Every localized metric needs "how far along the galaxy's own ellipse is this
pixel, and at what angle". For a pixel offset $(\Delta x, \Delta y)$ from the
fitted centre, rotate into the major-axis frame and define

$$
a = \sqrt{x'^2 + (y'/q)^2},
\qquad
x' = a\cos\theta,
\qquad
y' = a\,q\sin\theta
$$

so $a$ is the **semi-major axis (SMA)** of the ellipse passing through that
pixel and $\theta$ is the in-ellipse azimuth measured from the major axis.
Contours of constant $a$ are exactly the model's own isophotes. GALFIT measures
PA counter-clockwise from $+y$, so the rotation uses $\mathrm{PA} + 90^\circ$.

**$a$ vs $R_e$ (why plots use $a/R_e$).** $R_e$ is **not** a circular radius.
In GALFIT it is the SMA of the elliptical isophote that encloses half the total
*model* light — the half-light ellipse, stretched by the fitted $q$. The
coordinate $a$ is that same SMA measured for an arbitrary pixel (or fitted
isophote). So $a/R_e = 1$ is “on the half-light ellipse,” $a/R_e = 2$ is an
ellipse twice as large in every linear dimension (homothetic), not a circle of
radius $2R_e$.

```
elliptical_coords(shape, xc, yc, q, pa):
    dx, dy   = pixel grid - (xc, yc)
    phi      = radians(pa)
    x_maj    = -dx sin(phi) + dy cos(phi)      # along the major axis
    y_min    =  dx cos(phi) + dy sin(phi)      # along the minor axis
    a        = hypot(x_maj, y_min / max(q, 1e-3))
    theta    = atan2(y_min / q, x_maj)
    return a, theta
```

### 3.3 Valid pixels

A pixel enters any metric only if data, model and $\sigma$ are finite,
$\sigma > 0$, and `host_mask == 0`. When `cutout_meta.json` reports
`n_fit_components > 1`, pixels within $1 R_e$ of any *other* fitted component are
dropped as well. Without that, a neighbour whose own model is slightly wrong
would dump its residual into the host's aperture and every metric below would
inherit it.

### 3.4 Sigma calibration and a short $\chi^2/\nu$ primer

**What $\chi^2/\nu$ is.** GALFIT
([Peng et al. 2002](https://doi.org/10.1086/340952);
[Peng et al. 2010](https://doi.org/10.1088/0004-6256/139/6/2097))
and the suite score a fit by

$$
\chi^2 = \sum_j\left(\frac{I_j - I^{\rm model}_j}{\sigma_j}\right)^2,
\qquad
\frac{\chi^2}{\nu}
\quad\text{with}\quad
\nu = N_{\rm pix} - k_{\rm free}.
$$

Each term is the residual in units of the per-pixel uncertainty $\sigma_j$ from
`host_sigma.fits`. If the model is correct and $\sigma$ is honest Gaussian noise,
$E[\chi^2/\nu] = 1$. Three things break that reading: (i) a wrong model
(bright/resolved hosts *should* sit above 1 — a single Sérsic is never perfect;
see [Andrae et al. 2010](https://arxiv.org/abs/1012.3754) on misusing reduced
$\chi^2$), (ii) a wrong absolute $\sigma$ scale ($\chi^2$ moves as $1/\sigma^2$), and
(iii) dilution of the global value by empty sky after the Re-separation ROI
change. Localized $\chi^2/\nu$ inside $2R_e$ addresses (iii); this section
addresses (ii).

**Production rescale gate (Phase 3a, pre-fit).** The pipeline builds
$\sigma_{\rm invvar} = 1/\sqrt{\mathrm{invvar}}$ and compares it to the sky of
the *data*:

$$
k = \frac{1.4826\times\mathrm{MAD}(\text{sky pixels of the cutout})}
{\mathrm{median}(\sigma_{\rm invvar}\text{ on those same pixels})}.
$$

$\mathrm{MAD}$ is the median absolute deviation,
$\mathrm{median}_i\lvert x_i - \mathrm{median}(x)\rvert$
([Hampel 1974](https://doi.org/10.1080/01621459.1974.10482962);
[Rousseeuw & Croux 1993](https://doi.org/10.1080/01621459.1993.10476408)).
For Gaussian noise of width $\sigma$,
$\mathrm{MAD} = \sigma\cdot\Phi^{-1}(0.75)$ with
$\Phi^{-1}(0.75)\approx 0.6745$, so

$$
\sigma \approx \frac{\mathrm{MAD}}{0.6745} = 1.4826\times\mathrm{MAD}.
$$

That factor is **not** $\sqrt{2}$. ($\sqrt{2/\pi}\approx 0.80$ appears elsewhere,
in the RFF noise subtraction, because that uses the *mean* absolute deviation
rather than the median.) The denominator is simply the median of the invvar-derived
sigma map on the same sky pixels — the map's absolute scale, ignoring its spatial
structure.

Rescale fires **only** when $k\notin[0.5,\,2.0]$. Mild mismatches are left alone;
only vacuous unit-scale failures (typically Legacy Survey flux/invvar unit
mismatches with $k\sim10^3$–$10^4$;
[Dey et al. 2019](https://doi.org/10.3847/1538-3881/ab089d))
get a global multiply of the usable sigma map by $k$. On the current tree:
52/64 hosts `scale OK` at $k\sim0.9$; 12/64 rescaled.

**Verification ratio (post-fit audit).** After the fit, the suite measures how the
*residual* sky compares to the sigma map GALFIT actually used (already
gate-corrected if it was):

$$
r_\sigma = \frac{1.4826\times\mathrm{MAD}(\text{residual on sky})}
{\mathrm{median}(\mathtt{host\_sigma}\text{ on sky})}
$$

on unmasked pixels beyond $3R_e$ (fallback $1.5R_e$; `nan` below 20 pixels).
Cohort median $r_\sigma = 0.892$ [0.839, 0.959] — a mild, expected offset from
different sky definitions, model-wing leakage, and the fact that invvar is an
estimate. Nobody sits outside $(0.5,\,2)$. This is **not** evidence that everyday
fits had their $\chi^2$ forced to 1; it is a calibration thermometer.

**Corrected $\chi^2/\nu$.** Because $\chi^2\propto 1/\sigma^2$ and
$r_\sigma = \sigma_{\rm emp}/\sigma_{\rm map}$,

$$
\left.\frac{\chi^2}{\nu}\right|_{\rm corr}
= \frac{(\chi^2/\nu)}{r_\sigma^2}.
$$

With $r_\sigma = 0.89$ this raises the value by $\sim 1.26\times$, putting the
comparison to 1 back on an empirically calibrated noise scale. Applied to both
the global (whole-image) and localized ($2R_e$) values; see §4.1.

### 3.5 Rebuilding GALFIT's model analytically

The Fourier check (§4.3) needs $\partial M/\partial q$ of the **PSF-convolved**
model. The residual is defined against that plane, so the derivative has to be
too. The only honest route is to rebuild $M$, perturb $q$, and finite-difference —
which is only legal if the rebuild matches GALFIT's plane to much better than the
effect being measured.

**The continuous Sérsic**
([Sérsic 1963](https://ui.adsabs.harvard.edu/abs/1963BAAA....6...41S/abstract);
[Graham & Driver 2005](https://doi.org/10.1071/AS05001)):

$$
I(a) = I_e \exp\!\left[-b_n\!\left(\left(\frac{a}{R_e}\right)^{1/n} - 1\right)\right],
\qquad
F_{\rm tot} = 2\pi n R_e^2\, q\, e^{b_n} b_n^{-2n}\, \Gamma(2n)\, I_e
$$

with the Ciotti & Bertin $b_n$
([Ciotti & Bertin 1999](https://arxiv.org/abs/astro-ph/9911078)).
Magnitude → total flux via the zero point; the stamp is then scaled by
$F/F_{\rm analytic}(I_e=1)$, **not** by the stamp sum, so truncated wings cannot
change the total light.

**How GALFIT itself samples pixels**
([Peng et al. 2002](https://doi.org/10.1086/340952), §3;
[GALFIT Technical FAQ](https://users.obs.carnegiescience.edu/peng/work/galfit/TFAQ.html)).
GALFIT
does **not** uniformly oversample the whole stamp. It decides per pixel, by
distance from the component centroid:

- far from centre → evaluate the analytic function once, at the **pixel centre**;
- near the centre, where the second derivative is large → subdivide the pixel into
  a square $k\times k$ grid, evaluate on each sub-cell, and sum (i.e. integrate
  over the pixel area).

For the pathological Nuker cusp it switches to an elliptical polar grid inside
$r<3$ px; for Sérsic profiles the square sub-grid is what GALFIT uses. Convolution
with the PSF is then an FFT. The central few pixels carry a large fraction of the
light at high $n$, so under-sampling them before convolution permanently mis-states
the convolved core — which is why a naive point-sampled rebuild fails on compact
hosts.

**What “chasing GALFIT” means.** The mathematically exact pixel integral of a
Sérsic is *not* the target. GALFIT's adaptive, finite-$k$ scheme is. Pushing a
uniform high-order integrator harder can therefore move the rebuild *away* from
`out.fits` HDU 2 even as it gets closer to mathematical truth — observed directly
on `20230708A` ($R_e=1.5$ px, $n=6$), where some large $k$ values are worse than
smaller ones. The match-to-GALFIT error surface is **non-monotonic** in $k$.

**What the suite does.** Mirror GALFIT's distance dependence: point-sample the
whole stamp, and replace only an inner circular core with a $k\times k$ block
average. Neighbours get the same treatment about their own centres. Then
convolve with `proto_image.fits` and add sky. Because the error is non-monotonic
in both the subdivision and the core size, the suite searches a small product
grid — dense in $k$ (including even factors; block averaging does not require
odd $k$) and three core radii ($5$ px, $8$ px, and the whole stamp):

```
model_reconstruction_error(host):
    for core_radius in (5, 8, whole_stamp):
        for k in (1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 15, 21):
            recon = build_model(host, oversample=k, core_radius=...)
            err = max|recon - galfit_model| / peak(galfit_model - sky)
    keep the (k, core_radius) that minimises err
```

This is still an approximation of GALFIT (we do not re-implement the polar Nuker
integrator or GALFIT's exact internal $k$ schedule), but it is the right *shape*
of approximation, and the per-host comparison against HDU 2 is the ground truth
for whatever residual remains. The Fourier calibration refuses the rebuild when
the residual exceeds 3% of peak.

**Why step 6 needs this.** Finite-differencing $M(q\pm\varepsilon)$ only yields a
usable $\partial M/\partial q$ if rebuild error $\ll$ the finite-difference signal.
With a $0.1\%$ rebuild and $\varepsilon_q = 0.02$, that headroom is what lets
§4.3 convert residual $A_2$ into a number in units of $q$ itself.

---

## 4. The checks

Each check below is described as: **why it exists**, **what it measures**,
**how it is implemented**, **what can make it lie**, and **what it outputs**.

---

### 4.1 Chi-squared per degree of freedom — `chi2_local.py`

**Why.** The obvious first question, and the one most often misused
([Andrae et al. 2010](https://arxiv.org/abs/1012.3754)).

**The problem it solves.** The global $\chi^2/\nu$ in `fit.log`
([Peng et al. 2002](https://doi.org/10.1086/340952)) is computed over
the whole fitting region. After the Re-separation ROI change the host can occupy
a small fraction of the stamp, so sky pixels — which fit perfectly, being pure
noise — dilute the host residuals toward 1. A compact galaxy can be fit badly and
still report $\chi^2/\nu \approx 1$. Localizing the statistic to the host removes
that dilution.

**Definition.**

$$
\left.\frac{\chi^2}{\nu}\right|_{\rm loc} = \frac{1}{N_{\rm pix} - k}
\sum_{a \le 2R_e}\left(\frac{I_j - I^{\rm model}_j}{\sigma_j}\right)^2
$$

over valid pixels, with $k$ = 7 free parameters per Sérsic component plus 1 for
sky. Computed at both $1R_e$ and $2R_e$.

**Implementation.**

```
run(host):
    global:     last "Chi^2 = ..., ndof = ..." pair in fit.log   (regex)
    localized:  a, _ = elliptical_coords(...)
                sel  = valid_mask(host) & (a <= n_re * host.re)
                chi2 = sum(resid[sel]^2 / sigma[sel]^2)
                nu   = npix - (7 * n_components + 1)
    also:       sigma_calibration_ratio(host)            # §3.4
                chi2nu_global_corrected     = chi2nu_global     / r_sigma^2
                chi2nu_local_2re_corrected  = chi2nu_local_2re  / r_sigma^2
```

**What makes it lie.** Two things, both reported alongside. Sigma miscalibration
($r_\sigma \ne 1$) rescales it directly. And SNR: a bright, well-resolved galaxy
*should* show a large $\chi^2/\nu$, because a single Sérsic is wrong at the
percent level and there are enough photons to see it. The population plot of
$\chi^2/\nu$ against `snr_win` exists to confirm that trend — the interesting
host is one whose high $\chi^2/\nu$ is *not* explained by SNR. A second trend
against $R_e/\mathrm{FWHM}$ checks the same idea along resolution:
$\chi^2/\nu$ tracks whether model mismatch is *visible after PSF convolution*,
not whether the fit is trustworthy — unresolved hosts can look clean in
$\chi^2$ while $q$ is still weakly constrained. Neither trend is a pass/fail
gate; $\chi^2/\nu$ here is a visual / ranking diagnostic only.

**Outputs.** `chi2nu_global`, `chi2nu_local_1re`, `chi2nu_local_2re`,
`chi2nu_global_corrected`, `chi2nu_local_2re_corrected`, `chi2_local_npix_*`,
`chi2_local_nu_*`, `sigma_calibration_ratio`, `sigma_cal_npix`, `sky_mad_adu`,
`sigma_map_median_adu`, `residual_closure`, `re_over_fwhm`.

**Final document.** Report the corrected values
`chi2nu_global_corrected` (whole image) and `chi2nu_local_2re_corrected`
(within $2R_e$).

---

### 4.2 Localized residual flux fraction — `rff.py`

**Why.** $\chi^2/\nu$ asks whether residuals match the *noise model* and weights
bright pixels hard ($r^2/\sigma^2$). RFF asks a different question: what
*fraction of the galaxy's light* is left unmodelled, after subtracting the
noise you would see even for a perfect fit. That fraction is comparable across
bright and faint hosts in a way raw $\chi^2$ is not.

**Definition**
([Hoyos et al. 2012](https://doi.org/10.1111/j.1365-2966.2011.19918.x);
[arXiv:1109.6828](https://arxiv.org/abs/1109.6828)),
over an aperture $A$:

$$
\mathrm{RFF}
= \frac{\displaystyle\sum_{j\in A}\lvert I_j - M_j\rvert
- \sqrt{2/\pi}\,\sum_{j\in A}\sigma_j}
{\displaystyle\sum_{j\in A} M_j^{\rm gal}}.
$$

Three pieces, separately:

1. **Numerator, first term.** $\sum\lvert I-M\rvert$ is the total absolute
   residual flux in $A$. Absolute value so positive and negative lobes do not
   cancel (a wrong $q$ produces opposite-signed lobes that would vanish in a
   signed sum).
2. **Numerator, second term.** For zero-mean Gaussian noise,
   $E[\lvert x\rvert]=\sigma\sqrt{2/\pi}\simeq 0.798\,\sigma$. Even a perfect
   model therefore contributes $\sim 0.8\sum\sigma$ to the absolute residual.
   Subtracting that debiasing term makes pure noise give
   $\mathrm{RFF}\approx 0$ at any depth — the point of the statistic.
3. **Denominator.** $\sum M^{\rm gal}$ is the model *galaxy* flux in the same
   aperture (sky stripped from `out.fits` HDU 2). RFF is then a dimensionless
   fraction of host light. Including sky in the denominator would let a large
   empty aperture drive RFF to 0 by dilution.

**Uncertainty.** $\mathrm{Var}(\lvert x\rvert)=(1-2/\pi)\sigma^2$, so under
independent pixels

$$
\sigma_{\rm RFF}
= \frac{\sqrt{(1-2/\pi)\sum_{j\in A}\sigma_j^2}}{\sum_{j\in A} M_j^{\rm gal}}.
$$

At low SNR the debiasing subtraction is noisy and RFF can go **negative**. That
is expected arithmetic: $-0.05\pm 0.06$ means "consistent with a perfect model,"
not "better than perfect."

**Locality is an assumption.** The suite evaluates $A$ as elliptical apertures
tied to the *fitted* geometry: $a\le 1R_e$, $a\le 2R_e$, and the
$1$–$2R_e$ annulus, with the same valid-pixel / neighbour mask as §3.3.
That does **not** guarantee the whole host is inside $A$:

- $R_e$ is the half-light radius of the fitted Sérsic, so even a perfect match
  leaves $\sim$ half the light outside $1R_e$ and a non-negligible wing outside
  $2R_e$ (especially at high $n$).
- If the true galaxy is more extended than the fit, or the fit pinned $R_e$ low,
  flux (and residual structure) outside $2R_e$ is simply not scored.
- Conversely, if $R_e$ is over-large, $A$ includes more sky and the denominator
  can pick up model wing that is mostly prior, not data.

So $2R_e$ is a *convention* matched to where most of the constrained light lives,
not a proof that every host photon is counted. The three apertures exist partly
to expose that: core vs $2R_e$ vs annulus.

**Implementation.**

```
run(host):
    model_gal = host.model - sky_level        # galaxy only in the denominator
    a, _      = elliptical_coords(...)        # fitted q, PA, Re
    ok        = valid_mask(host)
    for aperture in (a <= 1Re, a <= 2Re, 1Re < a <= 2Re):
        sel   = ok & aperture
        rff   = (sum|resid[sel]| - sqrt(2/pi)*sum(sigma[sel])) / sum(model_gal[sel])
        err   = sqrt((1 - 2/pi) * sum(sigma[sel]^2)) / sum(model_gal[sel])
    rff_outer_minus_inner = rff_annulus_1_2re - rff_1re
```

**Why the annulus separately.** The $1$–$2R_e$ annulus is low surface brightness
and large area, so a small sky offset moves it while barely touching the core.
A large `rff_outer_minus_inner` points at background error rather than a wrong
model shape.

**$\chi^2/\nu$ vs RFF — why both, and which to keep if forced.** They are not
substitutes:

| | $\chi^2/\nu$ | RFF |
|---|---|---|
| Weights | $r^2/\sigma^2$ — brightest, noisiest-weighted pixels dominate | $\lvert r\rvert$ — every ADU of residual counts once |
| Scale | absolute, tied to the noise map | fraction of model galaxy light |
| Depth dependence | rises with SNR when model error becomes visible | debiased toward 0 for pure noise at any depth |
| Cancelling lobes | signed squares still accumulate | absolute value keeps opposite-signed geometry errors |
| Best at answering | "are residuals consistent with $\sigma$?" | "how much of the host light is unmodelled?" |

Keep both when possible: $\chi^2/\nu$ is the right thermometer for noise-model
and SNR trends (§4.1); RFF is the right single number for "how dirty is the
residual relative to the galaxy," and is the better triage scalar across a
magnitude-mixed cohort. **If only one were allowed for visual confirmation of
fit dirtiness, choose RFF** (prefer `rff_2re` with its error bar) — it is
directly in units of host light and does not require a mental SNR correction.
Keep $\chi^2/\nu$ when the question is specifically about the noise map or the
population SNR/resolution trends.

#### How to read $\chi^2/\nu$ and RFF (good / bad / expected)

Neither number has a universal “pass at $X$” cut independent of depth and
morphology. Use the rules of thumb below, then compare to **this cohort**.

**$\chi^2/\nu$ (prefer the sigma-corrected values on the panel).**

| value | reading |
|---|---|
| $\sim 1$ after correction | residuals consistent with the noise map (ideal Gaussian case) |
| $\gg 1$ on a bright, resolved host | usually **real structure** the single Sérsic cannot absorb (arms, clumps, neighbours) — expected, not automatically a failed $q$ |
| $\ll 1$ | almost always a **noise-map** problem (wrong units / overestimated $\sigma$), not a “better than perfect” fit |
| rises with SNR / falls for faint hosts | normal: model mismatch becomes visible only when $S/N$ is high |

Peng’s own FAQ is blunt: in galaxy fitting the absolute $\chi^2/\nu$ is
“mostly meaningless” because residuals are dominated by imperfect models, not
Gaussian noise
([GALFIT $\chi^2$ FAQ](https://users.obs.carnegiescience.edu/peng/work/galfit/CHI2.html)).
Use it **relatively** (host vs host, global vs $2R_e$, before/after a change)
and always prefer `chi2nu_*_corrected` so a wrong $\sigma$ scale is not mistaken
for a dirty residual. On the 53-cut, median
$\chi^2/\nu|_{2R_e,{\rm corr}}\approx 1.35$ overall, but
$\sim 3.2$ for $m<19$ vs $\sim 1.0$ for $m>20.5$ — brighter hosts look
“worse” in $\chi^2$ because structure is resolved, not because the fainter ones
are better models.

**RFF (prefer `rff_2re` $\pm$ error).**

| value | reading |
|---|---|
| $\approx 0$ (within $\sim\sigma_{\rm RFF}$) | residual consistent with noise after the Hoyos debias — the design point of RFF |
| slightly negative | common at low SNR; **not** “better than perfect,” just noisy debiasing |
| $\lvert\mathrm{RFF}\rvert \lesssim 0.05$ | clean single-Sérsic residual for survey work |
| $\lvert\mathrm{RFF}\rvert \gtrsim 0.10$ | Hoyos et al. (2012) regime where residuals often flag mergers / strong disturbances; our `flag_rff_high` cut |
| $\mathrm{RFF} \gtrsim 0.5$ | commonly used in modern JWST/HST pipelines as a **hard reject** for unusable single-Sérsic fits (e.g. Lee et al. 2024; Ormerod et al. 2024) |

Literature medians for “typical” galaxies after a single Sérsic are often
$\mathrm{RFF}\sim 0.1$–$0.12$ in deep HST/JWST samples
([Lee et al. 2024](https://arxiv.org/abs/2312.04899)), not zero — real disks
and irregularities leave light. There is **no** published universal
RFF$(m)$ or $\chi^2(m)$ curve to import: both scale with depth, PSF sampling,
and morphology. The practical calibration is **internal**: on this 53-cut,
median $\mathrm{RFF}_{2R_e}\approx -0.002$ with p84 $\approx +0.07$, and
median $\lvert\mathrm{RFF}\rvert$ is higher for bright hosts ($m<19$:
$\sim 0.06$) than for $19<m<20.5$ ($\sim 0.006$), the same SNR/structure
trend as $\chi^2$. So: treat $\lvert\mathrm{RFF}\rvert>0.10$ as “look at the
panel,” not as automatic science rejection; treat $\chi^2/\nu_{\rm corr}>1$ on
a mag$\lesssim 19$ host as the default, not a surprise.

**Outputs.** `rff_1re`, `rff_2re`, `rff_annulus_1_2re`, each with `_err`, `_sig`,
`_npix`, `_model_flux`; plus `rff_outer_minus_inner` and `sky_level_adu`.

---

### 4.3 Fourier decomposition of the residual — `fourier.py`

#### Why $m=2$ matters for axis ratios at all

An ellipse is a quadrupole distortion of a circle: stretching one axis and
squeezing the other is exactly an $m=2$ perturbation in polar angle. That is
why classical isophote work expands surface brightness as
$I(\theta)=\sum (A_m\cos m\theta+B_m\sin m\theta)$ and treats the $m=2$ terms
as the carriers of ellipticity and PA
([Jedrzejewski 1987](https://doi.org/10.1093/mnras/226.4.747);
[Bender et al. 1988](https://doi.org/10.1051/aas:1988002);
[Lauer 1985](https://doi.org/10.1086/163219)).
Higher even modes ($m=4$) describe boxy/disky deviations *from* an ellipse;
odd modes describe lopsidedness. So if the fitted $q$ (or PA) is wrong, the
residual after subtracting that elliptical model is not random — it must show
up primarily as $m=2$. That is the sole reason this check stares at $A_2$ and
$B_2$ rather than at $\chi^2$ or RFF.

#### What this check is even trying to do

GALFIT returns a best-fit axis ratio $q_{\rm fit}$. The science question is
whether that number is right. $\chi^2$ and RFF only say “the residual has some
size.” They do not say “the residual looks like a wrong $q$.”

$\delta q$ is the answer to a sharper question:

> If the true axis ratio were $q_{\rm fit}+\delta q$ instead of $q_{\rm fit}$,
> and that were the *only* thing wrong, how large would $\delta q$ have to be
> to produce the $m=2$ pattern we actually see in the residual?

So $\delta q$ is **not** GALFIT’s formal $\sigma_q$, and not a Monte Carlo error
bar. It is a *residual-inferred correction* (or inconsistency): a signed number
in the same units as $q$, read off the shape of $I-M$. If
$|\delta q|\ll \sigma_q$ and the pattern is reliable, the residual is not asking
for a different axis ratio. If $|\delta q|$ is large and the reliability gate
passes, the residual *is* asking for a different $q$ — either the fit is biased
or real $m=2$ structure (bar/spiral) is present (phase test below).

This is the same spirit as classical Fourier analysis of isophotes
([Jedrzejewski 1987](https://doi.org/10.1093/mnras/226.4.747);
[Bender et al. 1988](https://doi.org/10.1051/aas:1988002)), turned into a
direct estimator of a parameter error rather than a shape coefficient alone.

#### Why a pure continuous analytic picture is not enough

It is tempting to stay in the Ciotti–Bertin continuum: differentiate the Sérsic,
never touch a pixel. That closed form is derived below and still computed as
`fourier_dq_radial`. It is **not** what the calibrated estimator uses as truth,
because the residual that exists in the real world is

$$
R = I_{\rm data} - M_{\rm GALFIT},
$$

and both sides are **detector pixel arrays**. $I_{\rm data}$ was never a
continuous function in our pipeline; $M_{\rm GALFIT}$ is GALFIT’s
pixel-integrated, PSF-convolved model plane
([Peng et al. 2002](https://doi.org/10.1086/340952)). Masks, $\sigma_j$, and
annulus membership are pixel sets. The Fourier step is a weighted least-squares
sum over those pixels. There is no continuous $R(a,\theta)$ to expand until you
commit to the same sampling.

What §3.5’s oversampling is doing is **not** “4.3 needs a fine continuum for
Fourier mathematics.” Fourier mathematics here is discrete by construction.
Oversampling is only so that when we build

$$
K_q \approx \frac{M_{\rm rebuild}(q+\varepsilon)-M_{\rm rebuild}(q-\varepsilon)}{2\varepsilon}
$$

by finite difference, $M_{\rm rebuild}$ lives on the **same pixelization
convention as $M_{\rm GALFIT}$**. Without that, $K_q$ would be the derivative of
a different function than the one that defined $R$, and
$\delta q = A_2/c_q$ would be nonsense. The continuous analytic derivative is
exactly the approximation that under-recovered injected $\delta q$ by $\sim 18\%$
once a real PSF and pixel grid are present.

In short: pixelate because $R$ is pixelated; oversample when rebuilding so the
response kernel matches GALFIT’s pixels, not because Fourier requires a
continuum.

#### Why elliptical annuli (not one global Fourier)

A geometry error is a function of elliptical radius $a$: bars live in a limited
radial range, spirals wind with $a$, and a wrong global $q$ should be roughly
constant across $a$. A single stamp-wide $m=2$ coefficient would mix those
behaviours into one number and erase the phase test that separates them.
Elliptical annuli (geometry of §3.2) keep the radial coordinate the model
actually uses, so each ring’s $A_2(a),B_2(a)$ is local.

**How Fourier picks its rings (not photutils).** §4.5’s free isophotes and
§4.3’s Fourier rings are **different samplings**. Photutils steps geometrically
in SMA ($\times 1.1$) and fits a new ellipse each time. Fourier fixes GALFIT’s
$(q,\mathrm{PA})$, builds the elliptical-radius map $a$, and slices it into
**fixed-width** elliptical annuli via `azimuthal_annuli`:

$$
\mathrm{width} = \max\!\big(1\,\mathrm{px},\,\mathrm{FWHM}/2\big),
\quad
a \in [\!\max(0.2R_e,\,\tfrac12\mathrm{width}),\; 3R_e].
$$

On a barely-resolved host ($\mathrm{FWHM}\!\sim\!R_e$) that width is $\sim 2$ px,
so you get only a handful of fat rings (e.g. midpoints at
$0.63,\,1.26,\,1.89,\,2.52,\,3.15\,R_e$) — not the ~20 thin isophotes on the
middle panel. Rings with fewer than **20 valid pixels** get no $\delta q$ (the
azimuthal WLS has 9 coefficients; $<20$ pixels is underdetermined / pure noise).
That is why the Fourier panel can skip the first midpoint and start near
$\sim 1.3R_e$: the ring exists, but `npix < 20` → $\delta q = \mathrm{NaN}$.

On each usable annulus, decompose the residual azimuthally by inverse-variance
least squares on the valid pixels (not an FFT — masking makes $\theta$-sampling
non-uniform; [Jedrzejewski 1987](https://doi.org/10.1093/mnras/226.4.747)):

$$
R(a,\theta) = A_0(a) + \sum_{m\ge1}\Big[A_m(a)\cos m\theta + B_m(a)\sin m\theta\Big].
$$

#### Why the m = 2 mode is the diagnostic

Suppose the true axis ratio is $q + \delta q$. To first order

$$
R \approx \frac{\partial M}{\partial q}\,\delta q.
$$

With $M = f(a)$ and the elliptical parametrization of §3.2,

$$
\frac{\partial a}{\partial q} = -\frac{a}{q}\sin^2\theta
= -\frac{a}{2q}\big(1 - \cos 2\theta\big)
\quad\Longrightarrow\quad
\frac{\partial M}{\partial q} = -\frac{a f'(a)}{2q}\big(1 - \cos 2\theta\big).
$$

So an axis-ratio error puts power into $\cos 2\theta$ — four alternating lobes
**aligned with the axes**. For a PA error $\delta\phi$ (radians),

$$
\frac{\partial a}{\partial \phi} = \frac{a}{2}\Big(\frac{1}{q} - q\Big)\sin 2\theta
\quad\Longrightarrow\quad
\frac{\partial M}{\partial \phi} = \frac{a f'(a)}{2}\Big(\frac1q - q\Big)\sin 2\theta,
$$

the same quadrupole rotated $45^\circ$. Because $\cos 2\theta$ and $\sin 2\theta$
are orthogonal, $A_2$ and $B_2$ separate the two errors. Matching coefficients
and using $f'(a)<0$:

$$
\delta q = -\frac{2q\,A_2(a)}{a\,\lvert f'(a)\rvert},
\qquad
\delta\phi = -\frac{2\,B_2(a)}{a\,\lvert f'(a)\rvert\,(1/q - q)}.
$$

A **positive** $\delta q$ means the true galaxy is rounder than the model
($q_{\rm true}\approx q_{\rm fit}+\delta q$). The factor $(1/q-q)$ makes
$\delta\phi$ diverge as $q\to 1$: a round galaxy has no measurable PA.

#### Calibrating the conversion (PSF + pixels)

The closed form assumes the PSF does not change the $A_2\to\delta q$ map. It
does: convolution damps a quadrupole more than the radial gradient, so the
analytic conversion **under-recovers an injected $\delta q$ by $\sim 18\%$** at
this cohort’s resolution. The suite therefore builds $K_q$ and $K_{\rm PA}$ from
the §3.5 rebuild (same `proto_image.fits`, same pixel grid), projects them with
the **identical** design matrix and weights as $R$, and sets

$$
K_q = \frac{M_{\rm conv}(q+\epsilon) - M_{\rm conv}(q-\epsilon)}{2\epsilon},
\qquad
c_q(a) = \big[K_q\big]_{\cos 2\theta},
\qquad
\delta q(a) = \frac{A_2(a)}{c_q(a)}.
$$

Unbiased by construction for that projector. The analytic form remains as
`fourier_dq_radial` / fallback.

**What “injected $\delta q$” means (and what it does *not*).** Nothing is injected
into production science images. `tests/test_fourier_recovery.py` builds a
*synthetic* Sérsic, convolves it with a real `proto_image.fits`, differences
two copies that differ by a known $\delta q$ (or $\delta\mathrm{PA}$), and asks
the estimator to recover that known number. That is a unit-test ground truth,
not a science-sample perturbation. “Calibrated” = $A_2/c_q$ from the numerical
PSF-convolved kernels; “analytic” = the thin-shell formula from $A_2$ and
$d\langle M\rangle/da$. They differ because the analytic map ignores PSF damping
of the quadrupole relative to the radial gradient (~18% under-recovery here).

**Even / odd and crosstalk, in plain language.** Under reflection
$\theta\to-\theta$, $\cos 2\theta$ is even and $\sin 2\theta$ is odd. A pure
$\delta q$ response $K_q$ should therefore project almost entirely onto
$\cos 2\theta$; a pure $\delta\mathrm{PA}$ response $K_{\rm PA}$ onto
$\sin 2\theta$. Crosstalk is the fraction of $K_q$ that leaked into the PA
channel on that annulus:

$$
\texttt{fourier\_kernel\_crosstalk}(a)
= \Big\lvert\frac{[K_q]_{\sin 2\theta}}{[K_q]_{\cos 2\theta}}\Big\rvert.
$$

Ideal value $\approx 0$. A large value means the finite-difference kernels (or
the local sampling) are not cleanly separating the two channels — e.g. extreme
$q$, heavy masking, or a rebuild mismatch. It is **reported always** and is a
useful QA field. It does **not** gate or hide $\delta q$: if PA and $q$ are
mixed, the geometry is usually untrustworthy anyway (except near-degenerate
$180^\circ$ flips), so suppressing $\delta q$ would hide the hosts that most need
inspection. Flagging still uses $|\delta q|$ against its threshold; read
crosstalk alongside when interpreting.

#### Step D in detail — the ratio on one annulus

Fix an annulus. After the least-squares fit you have a scalar $A_2$ (ADU): how
much $\cos 2\theta$ is in the residual on those pixels. After projecting $K_q$
with the **same** design matrix and weights you have a scalar $c_q$ (ADU per
unit $q$): how much $\cos 2\theta$ appears in the model when $q$ moves by 1.

If the residual’s quadrupole were caused only by a wrong $q$, then
$A_2 = c_q\cdot\delta q$, hence

$$
\delta q(a)=\frac{A_2(a)}{c_q(a)}.
$$

That is ordinary linear response: observed amplitude divided by response per
unit parameter. Using the same projector on both sides is what makes PSF
damping, pixelation, masking, and annulus width cancel — they affect $A_2$ and
$c_q$ equally. The PA channel is identical with $B_2$ and $c_{\rm PA}$. Errors
propagate as $\sigma_{\delta q}=|A_2^{\rm err}/c_q|$ from the LSQ covariance on
$A_2$ (treating $c_q$ as known from the noise-free rebuild).

#### Step E in detail — “collapse” = one summary number from many rings

After Step D you still have a **profile** $\delta q(a_i)$ — one estimate per
annulus. “Collapse” here is not a physical collapse of the galaxy; it is just
reducing that curve to the single scalar `fourier_dq` used in tables and flags.
The reduction is an **inverse-variance weighted (IVW) mean** over annuli with
**$a\le 2R_e$ only** (points plotted beyond $2R_e$ are diagnostic, not in the
mean):

$$
\delta q
= \frac{\sum_i w_i\,\delta q(a_i)}{\sum_i w_i},
\qquad
w_i = 1/\sigma_{\delta q,i}^2,
\qquad
\sigma_{\delta q}
= \big(\sum_i w_i\big)^{-1/2}.
$$

**Why cut at $2R_e$ (physical).** For a Sérsic, most of the light that actually
constrains $(q,\mathrm{PA},R_e,n)$ lives inside $\sim 2R_e$; beyond that the
pixels are sky-, neighbour-, and mask-dominated, so a large outer $\delta q$ is
usually not a statement about the fitted geometry. The same $2R_e$ aperture is
used for localized $\chi^2/\nu$, RFF, and `iso_dq_2re`, so the headline numbers
answer the same question. Opposite-signed rings inside $2R_e$ **cancel** in the
IVW mean — e.g. $\delta q(1.26R_e)\approx -0.19$ and
$\delta q(1.89R_e)\approx +0.05$ can yield `fourier_dq` $\approx -0.05$ — which
is correct for a *global* geometry error (should be roughly constant) and a
warning to read the curve when the profile flips sign (structure / noise, not
one $\Delta q$).

**Why inverse-variance weights (and not a plain / flux / equal mean).**
Treat each annulus estimate as an independent measurement
$\hat\delta q_i = \delta q^\star + \varepsilon_i$ with
$\varepsilon_i\sim\mathcal{N}(0,\sigma_i^2)$ and $\sigma_i\equiv\sigma_{\delta q,i}$
from the WLS covariance on $A_2$ (Step D). The log-likelihood for a constant
$\delta q^\star$ is

$$
\ln\mathcal{L}(\delta q^\star)
= -\tfrac12\sum_i\frac{(\hat\delta q_i-\delta q^\star)^2}{\sigma_i^2}+\mathrm{const}.
$$

Maximising gives exactly the IVW mean above; the Fisher information is
$\sum_i\sigma_i^{-2}$, so
$\mathrm{Var}(\widehat{\delta q})=(\sum_i\sigma_i^{-2})^{-1}$.
Among all unbiased linear estimators
$\sum_i\alpha_i\hat\delta q_i$ with $\sum_i\alpha_i=1$, the Gauss–Markov /
Aitken theorem says this choice uniquely minimises the variance
([Aitken 1935](https://doi.org/10.1017/S0370164600014346)).
Equivalently: any other weights (equal, flux-weighted, …) are a
**strictly worse** estimator of a constant $\delta q$ under these assumptions —
they throw away information encoded in $\sigma_i$.

Important: $\sigma_i$ is the uncertainty on **$\delta q$**, not the pixel ADU
noise and not “brighter ⇒ larger $\sigma$.” Brighter, better-sampled annuli
almost always have *smaller* $\sigma_{\delta q}$ (more photons constrain $A_2$),
so IVW **up-weights** the luminous rings and **down-weights** faint noisy ones.
That is the opposite of unfairly ignoring the bright galaxy — it is optimal
use of the well-measured core. If $\delta q(a)$ is *not* constant (bar, spiral,
sign flip), the IVW scalar is still the ML constant fit, but the profile itself
is the science product; that is why the panel plots $\delta q(a)$.

That is `fourier_dq` / `fourier_dq_err`. The full $\delta q(a)$, phase, and
amplitudes stay in `fourier_profiles.npz`.

**Implementation (full path the code takes).**

```
run(host):
    recon = model_reconstruction_error(host)   # §3.5: pick (k, core_radius)
    K_q, K_PA, eps = response_kernels(...)     # finite-diff rebuild ±eps, no sky/neighbours
    model_gal = host.model - sky
    a_map, theta = elliptical_coords(...)
    annuli = 0.2 Re -> 3 Re, width = max(1 px, FWHM/2)

    for each annulus with >= 20 valid pixels:
        design = [1, cos mθ, sin mθ] for m=1..4     # 9 columns
        w = 1/sigma^2
        coef, coef_err = WLS(design, resid, w)     # -> A0,A1,B1,...,A4,B4
        kq_coef = WLS(design, K_q,  w)            # same design, same w
        kpa_coef = WLS(design, K_PA, w)
        dq(a)  = A2 / kq_coef[cos2]
        dPA(a) = B2 / kpa_coef[sin2]
        also analytic dq_radial from A2 and d<M>/da
        psi2 = 0.5 * atan2(B2, A2)
        crosstalk = |kq_coef[sin2] / kq_coef[cos2]|   # reported, not a gate

    fourier_dq = IVW mean of dq(a) for a <= 2 Re     # Step E ("collapse")
    reliability gate (recon, n_annuli, Re/FWHM)      # not crosstalk
    write fourier.json + fourier_profiles.npz
```

Least squares rather than an FFT, because masked pixels make the azimuthal
sampling non-uniform and an FFT would alias the gaps into the coefficients.
Annulus width never goes below the resolution element.

#### Distinguishing a geometry error from bars / spirals (we do not subtract them)

A nonzero $m=2$ does not by itself mean the geometry is wrong — bars and spiral
arms produce $m=2$ too
([Elmegreen & Elmegreen 1985](https://doi.org/10.1086/163147);
[Rix & Zaritsky 1995](https://doi.org/10.1086/175701)).
The suite does **not** remove bar flux from the image. It **classifies** the
$m=2$ pattern using the phase
$\psi_2(a) = \tfrac12\arctan(B_2/A_2)$, fitted linearly against $a$:

| phase behaviour | interpretation | what to trust |
|---|---|---|
| flat $\psi_2$, $|\delta q(a)|$ smooth across $a$ | global geometry error | headline `fourier_dq` |
| flat $\psi_2$, $|\delta q(a)|$ peaked in a limited $a$ range | bar-like | report phase + radial peak; do not treat collapsed $\delta q$ as a $q$ correction |
| $\psi_2$ winds with $a$ (nonzero slope) | spiral | report `fourier_m2_phase_slope_*`; $\delta q$ is structure, not a bias |

**Good reporting when $m=2$ is interesting.** Publish the triad
(`fourier_dq`, `fourier_dq_err`, `fourier_m2_phase_slope_deg_per_re` /
`_sig` / `_scatter`), plus a one-line morphological note from the table above,
and for flagged hosts the $\delta q(a)$ curve from the NPZ (or the visual
panel). That is how a bar is “caught”: not by excision, but by refusing to
reinterpret a localized or winding quadrupole as a single global $\Delta q$.

**Other modes (panel scalars).** Amplitudes on the panel are
$\sqrt{A_m^2+B_m^2}\,/\langle M\rangle$ medians over usable annuli
(`fourier_m{1,2,3,4}_amp_median`), plus the $m=2$ phase slope
$\psi_2'=$ `fourier_m2_phase_slope_deg_per_re`:

| scalar | meaning |
|---|---|
| $m=1$ | lopsidedness / centroid drift — large if the nucleus is offset or one side is brighter |
| $m=2$ | quadrupole: geometry error **or** bar/spiral (use $\psi_2'$ to tell them apart) |
| $m=3$ | triangular / odd asymmetry (interactions, uneven spiral arms) |
| $m=4$ | boxy/disky isophotes ([Bender et al. 1988](https://doi.org/10.1051/aas:1988002)) — physical, not a fit failure |
| $\psi_2'\approx 0$ | flat phase → global $\delta q$ or bar; **winding** $\psi_2'$ → spiral (do not treat $\delta q$ as a $q$ correction) |

There is no hard “bad if $m_k>X$” cut: bright structured disks routinely show
medians $\sim 0.05$–$0.2$. Use modes to *classify* residual structure next to
$\delta q$. $A_0(a)$ beyond $2R_e$ remains the sky cross-check (§4.4), reported
in units of the sky MAD as a median over outer annuli.

#### The reliability gate

$\delta q$ is only meaningful if the estimator had leverage. It is marked
unreliable, and excluded from **flagging** (the number is still written), when
any of:

| gate | threshold | why |
|---|---|---|
| `model_recon_max_frac` | $> 0.03$ | the rebuilt model is not GALFIT's, so the kernel is wrong |
| `fourier_n_annuli_inner` | $< 4$ | too few rings inside $2R_e$ to average |
| `re_over_fwhm` | $< 0.8$ | host smaller than the PSF; see below |

`fourier_kernel_crosstalk` is always stored but is **not** a reliability gate
(see even/odd paragraph above).

**Why unresolved is fatal for this check (not for $q$ in general).** The PSF
does give the image *some* structure — a roundish blob — but that structure is
almost entirely the PSF itself. After convolution, a wide range of intrinsic
$(q,n,R_e)$ collapse to nearly the same observed shape. Consequences for §4.3
specifically:

1. **Too few independent annuli.** Annulus width is $\ge\mathrm{FWHM}/2$, and
   the useful radial range is only a couple of resolution elements across, so
   you cannot assemble $\ge 4$ rings inside $2R_e$ with $\ge 20$ pixels each.
2. **$c_q$ has no leverage.** Moving intrinsic $q$ barely changes the
   convolved image, so $K_q$ (and $c_q$) → 0; then $\delta q=A_2/c_q$ blows up
   or becomes noise-dominated. That is the mathematical meaning of “nothing to
   decompose.”
3. **Whatever $m=2$ you measure is often PSF ellipticity / pixel noise**, not
   host geometry — which is exactly the systematic §4.6 is built to catch at
   population level.

So unresolved does **not** mean “the galaxy has no shape” or “GALFIT’s $q$ is
automatically wrong.” It means **this residual→$\delta q$ converter has no
information**, and reporting a large $|\delta q|$ would be a false condemnation.
The gate exists to refuse the number, not to discard the host from science.

Without this gate a barely-resolved host would be condemned by arithmetic that
means nothing. On the science cut, 24/53 hosts pass the remaining gates — a
statement about the *imaging*, not about the estimator (see §7).

#### Statmorph / Gini–M20 (Lotz classification)

Checks 1–9 do **not** re-run or re-measure non-parametric morphology. Production
Phase Statmorph (`galfit_fitting/run_statmorph_pipeline.py`) already writes
`Output/<FRB>_all/statmorph_results.json` from `host_cutout.fits` (+ optional
`host_sigma.fits`, `proto_image.fits`) via the [statmorph](https://statmorph.readthedocs.io/)
package (Rodriguez-Gomez et al. 2019). Verification only **reads** those JSON
fields and applies the Lotz et al. (2008) Gini–$M_{20}$ plane cuts.

**Provenance.** For each host:

1. Phase 3a builds the cutout / segmentation inputs.
2. Phase Statmorph calls `statmorph.source_morphology` on the background-subtracted
   cutout (central segment only) and writes `statmorph_results.json`.
3. `run_gini_m20.py` loads `gini` and `m20` from that JSON for the science-cut
   cohort and writes `outputs/tables/gini_m20_53.csv`.
4. Visual panels (§4.9) show ${\rm Lotz}={\rm late}|{\rm early}|{\rm merger}$ in
   the header from the same JSON (same cut equations as `run_gini_m20.py`).

**Fields queried from `statmorph_results.json`:**

| JSON key | Quantity | Role here |
|---|---|---|
| `gini` | Gini coefficient $G$ | Lotz plane $x$-axis partner |
| `m20` | $M_{20}$ | second-order moment of the brightest 20% of light |
| `concentration`, `asymmetry`, `smoothness` | CAS | reported in CSV; not used for Lotz class |
| `gini_m20_merger`, `gini_m20_bulge` | $F$, $S$ | Rodriguez-Gomez / Lotz $S$–$F$ statistics from statmorph |
| `rpetro_circ_px`, `r20_px`, `r80_px` | Petrosian / $r_{20}$ / $r_{80}$ | radii used internally by CAS/Gini–$M_{20}$ |
| `sn_per_pixel`, `flag`, `flag_sersic` | quality | `flag>0` ⇒ measurement unreliable |

**Metric definitions (compact).** Following Abraham et al. (2003), Lotz et al.
(2004), and Conselice (2003):

- **Gini $G$** — inequality of the pixel flux distribution within the segmentation
  map ($G=0$ equal pixels; $G\to 1$ light in few pixels). Sensitive to clumpiness
  and central concentration without assuming elliptical isophotes.
- **$M_{20}$** — normalized second-order moment of the pixels that contain the
  brightest 20% of the galaxy’s flux; more negative for centrally concentrated
  light, less negative when bright clumps sit off-centre (mergers / double nuclei).
- **CAS** — Concentration $C=5\log_{10}(r_{80}/r_{20})$; Asymmetry $A$ from
  $180^\circ$ residual; Smoothness $S$ from high-frequency residual after boxcar
  smoothing (Conselice 2003).

**Lotz et al. (2008) classification** used in this repo
([ApJ 672, 177](https://doi.org/10.1086/523659); calibrations also in
[Lotz et al. 2004, ApJ 613, 898](https://doi.org/10.1086/421849)):

$$
\begin{aligned}
{\rm merger\ (ULIRG/merger\ region):\quad}
  & G > -0.14\,M_{20} + 0.33 \\
{\rm early\ (E/S0/Sa):\quad}
  & G \le -0.14\,M_{20} + 0.33
    \;\land\; G > 0.14\,M_{20} + 0.80 \\
{\rm late\ (Sb{-}Irr):\quad}
  & G \le -0.14\,M_{20} + 0.33
    \;\land\; G \le 0.14\,M_{20} + 0.80
\end{aligned}
$$

CSV labels are `merger` / `early` / `late` (plus `missing` / `unknown` if JSON
or $G,M_{20}$ absent). Panel header uses the same three science labels.

**statmorph $S$ / $F$ (secondary).** Rodriguez-Gomez et al. (2019) expose
`gini_m20_bulge` ($S$) and `gini_m20_merger` ($F$) as signed distances to the
Lotz dividing lines. `run_gini_m20.py` reports `bulge_side = early` if $S>0$,
else `late`. This is a convenient summary, **not** a second independent
classifier — the primary label for triage / panels is the Lotz+2008 region above.

**How to use it here.** Lotz class is **context** for dirty residuals and $m=2$
structure (e.g. true mergers often fail a clean single-Sérsic), not a
replacement for $\delta q$, RFF, or confirmation. Neighbours in the stamp can
drive a spurious `merger` label even when the host is a quiet disk.

**Citations.**

- Lotz, J. M., Primack, J., & Madau, P. 2004, ApJ, 613, 898
  ([10.1086/421849](https://doi.org/10.1086/421849)) — Gini–$M_{20}$ definitions.
- Lotz, J. M., et al. 2008, ApJ, 672, 177
  ([10.1086/523659](https://doi.org/10.1086/523659)) — $G$–$M_{20}$ early / late / merger cuts used here.
- Conselice, C. J. 2003, ApJS, 147, 1
  ([10.1086/375001](https://doi.org/10.1086/375001)) — CAS.
- Abraham, R. G., van den Bergh, S., & Nair, P. 2003, ApJ, 588, 218
  ([10.1086/373919](https://doi.org/10.1086/373919)) — Gini for galaxy morphology.
- Rodriguez-Gomez, V., et al. 2019, MNRAS, 483, 4140
  ([10.1093/mnras/sty3345](https://doi.org/10.1093/mnras/sty3345)) — `statmorph` implementation ($S$, $F$, CAS, Gini–$M_{20}$).

#### Validation

`tests/test_fourier_recovery.py` builds an analytic Sérsic at $q$, rebuilds at
$q + \delta q$, convolves both with the host's real `proto_image.fits`,
differences them, and checks what comes back:

| injected $\delta q$ | calibrated | analytic |
|---|---|---|
| 0.020 | 0.0193 | 0.0169 |
| 0.050 | 0.0459 | 0.0402 |
| −0.040 | −0.0428 | −0.0376 |

The few percent that remains is genuine second-order nonlinearity
($\delta q = 0.05$ is an 8% change in $q$) and shrinks with the perturbation,
which is the regime real residuals live in. An injected $\delta\mathrm{PA}$ of 3.000° comes
back as 2.996°. The suite also asserts that an injected $\delta q$ produces no
spurious $\delta\mathrm{PA}$ and vice versa, that the null case returns zero,
and that the calibration beats the analytic form. Six tests, all passing.

**Outputs.** `fourier_dq`, `fourier_dq_err`, `fourier_dq_sig`,
`fourier_dpa_deg`, `fourier_dpa_err_deg`, `fourier_dq_radial`,
`fourier_m{1,2,3,4}_amp_{max,median}`, `fourier_m2_phase_slope_deg_per_re`,
`fourier_m2_phase_slope_sig`, `fourier_m2_phase_scatter_deg`,
`fourier_a0_outer_adu`, `fourier_a0_sky_offset_sigma`, `fourier_reliable`,
`fourier_unreliable_reasons`, `fourier_kernel_crosstalk`, `fourier_n_annuli`,
`q_fourier_corrected`, `model_recon_*`; full radial profiles in
`fourier_profiles.npz`.

---

### 4.4 Sky perturbation — `sky_perturb.py`

**Why.** Background error is the dominant systematic for faint and extended
galaxies. It moves outer-disk flux, and outer-disk flux is what sets $q$ — the
core is nearly round under any PSF, so the axis ratio is carried by the faintest
part of the image, which is exactly the part a sky error corrupts.

**Choosing the perturbation amplitude.** The relevant quantity is the
**large-scale** sky uncertainty, not the per-pixel noise. Taken as

$$
\sigma_{\rm sky} = 1.4826 \times \mathrm{MAD}\big(\texttt{BACKGROUND}\big)
$$

over `image.cat` sources with `FLAGS == 0`. SExtractor writes a *local* background
estimate for every detected source, so the scatter of that column across the
field is a direct empirical measure of large-scale structure — gradients,
scattered light, flat-field residuals. The per-pixel sky RMS is recorded for
reference but deliberately **not** used: it is much larger than the uncertainty
on the *mean* sky, and using it would turn the test into an implausible worst
case rather than a realistic systematic. Cohort median
`sky_sigma_over_pixel_rms` = 0.085, i.e. the large-scale term is about 12×
smaller than per-pixel noise, as it should be.

**Method.** Two re-fits per host, both minimal-diff copies of production.

```
for label, value in (plus, sky_best + sigma), (minus, sky_best - sigma):
    stage outdir/galfit_sky_<label>/ with copies of:
        host_cutout.fits  host_sigma.fits  host_mask.fits  proto_image.fits
        constraints.txt   galfit.feedme
    edit galfit.feedme, minimally:
        every Sersic parameter -> reseeded at its best-fit value
        sky component line  " 1) <sky_best> 1"  ->  " 1) <value> 0"     # fixed
    strip any sky constraint from constraints.txt
    run_galfit(wkdir)                       # production's own runner
    parse out.fits header for q, PA, Re, n, mag
```

**Why reseed at the best fit.** Starting the optimizer from the answer isolates
the sky effect. If it started from the original guesses, part of the observed
$\Delta q$ would be optimizer path dependence rather than sky sensitivity, and
the two would be inseparable.

**Why strip the sky constraint.** GALFIT will happily override a fixed value if
`constraints.txt` still bounds it, which would silently turn the perturbation
into a no-op. The sky QA retry loop is disabled for the same reason. Everything
else — data, sigma, PSF, mask, region, zero point — is byte-identical to
production, and the run uses `run_galfit()` from
[`run_galfit_fitting.py`](../galfit_fitting/run_galfit_fitting.py) so the two
paths share crash handling.

**What it produces.** Report the full triad for every host:

| symbol | field | meaning |
|---|---|---|
| $q_0$ | `q_sky_0` (= production `b_a`) | baseline best-fit axis ratio |
| $q_{+1}$ | `q_sky_plus` | re-fit with sky fixed at $\mathrm{sky}+\sigma_{\rm sky}$ |
| $q_{-1}$ | `q_sky_minus` | re-fit with sky fixed at $\mathrm{sky}-\sigma_{\rm sky}$ |

Then
$\Delta q_{\rm sky} = \max\big(\lvert q_{+1}-q_0\rvert,\,\lvert q_{-1}-q_0\rvert\big)$
(`dq_sky`) is the per-host *systematic*, compared to GALFIT's *statistical*
`b_a_err` via `dq_sky_over_q_err`. Always look at $q_{-1},\,q_0,\,q_{+1}$
together: a one-sided jump (only plus or only minus moves) is as informative as
a symmetric swing.

**Outputs.** `q_sky_0`, `q_sky_plus`, `q_sky_minus`, `dq_sky`, `dq_sky_signed`,
`dq_sky_over_q_err`, the same for PA / $R_e$ / $n$ / mag, `sky_sigma_adu`,
`sky_sigma_source`, `sky_pixel_rms_adu`, `sky_sigma_over_pixel_rms`,
`sky_cat_nsources`, `converged_sky_{plus,minus}`, `status_sky_{plus,minus}`.

---

### 4.5 Isophote comparison — `isophotes.py`

**Why.** An independent, **non-parametric** view of the same geometry. GALFIT
imposes a single $q$ on the whole galaxy by construction; free isophotes do not.
So this is the only check that can answer *where in radius* a single-Sérsic
description breaks.

**Method.** `photutils.isophote.Ellipse` free fits — $q$ and PA free at every
isophote — on **both** the data and the **PSF-convolved model** (`out.fits`
HDU 2), using the same mask and same starting geometry.

**Ellipse is a parametrization, not a claim that SB is elliptical.** Real
isophotes need not be ellipses (bars, spirals, boxy/disky deviations). The
algorithm still fits the *best* ellipse at each $a$ (Jedrzejewski harmonic
iteration: sample intensity on a trial ellipse, adjust centre/$q$/PA until the
$m=1,2$ harmonics vanish). Non-elliptical structure is absorbed into higher
harmonics / fit residuals — it is not modelled as free shape. So
$q_{\rm data}(a)$ and $q_{\rm model}(a)$ are “best-fitting ellipse” axis ratios,
on both images, by construction.

**Why the convolved model and not the data alone.** Isophotal $q$ measured on
data is biased round by the PSF. Comparing it directly against GALFIT's
*deconvolved* $q$ would charge the PSF's rounding to the fit as an error.
Running the identical isophote machinery on the convolved model applies the
identical bias to both sides, so the difference isolates the fit.

**Sampling in $a$ (semi-major axis).** Here $a$ **is** the SMA of each fitted
ellipse (photutils `sma`, pixels). Growth is geometric with default
`linear=False`, `step=0.1`: $a_{n+1}=a_n\times 1.1$, from
`minsma = max(2, 0.2 Re)` out to `maxsma = min(3 Re, 0.45×stamp)` (and inward
from the seed). Typical recovery on the 53-cut: ~20 isophotes per host
(`iso_n_data` median 20). That is a discrete sampling of the continuum of
brightness contours — denser steps mostly correlate neighbours.

**The strategy search.** `fit_image` grows from a seed `sma0` and aborts if that
first ellipse fails, so one guess is not a strategy (on `20240119A`, `sma0=8`
gives 21 isophotes; 5 and 12 give none):

```
_search(...):
    seeds = 12 geom-spaced sma0 in [1.2*minsma, 0.9*maxsma]
    for (step, fix_center) in ((0.10, True), (0.15, True), (0.10, False)):
        for sma0 in seeds:
            fit_image(minsma, maxsma, step, nclip=3, fix_pa/eps free)
            keep the run with the most isophotes; stop early if >= 8
    reuse that strategy on the model (iso_same_strategy)
```

**What “profile” and $q(a)$ mean.** GALFIT’s $q_{\rm GALFIT}$ is **one number**.
A profile is $q$, PA, or $I$ **vs** $a$: $q_{\rm data}(a)$ on the cutout,
$q_{\rm model}(a)$ from the same free-ellipse fit on the convolved model (not
$q_{\rm GALFIT}$ copied). PSF convolution alone gives a mild $q(a)$ even for
constant-$q$ Sérsics; differencing cancels that.

**$\Delta q(a)$ and $\sigma_{\Delta q}$.**
$\Delta q(a)=q_{\rm data}(a)-q_{\rm model}(a)$, with

$$
\sigma_{\Delta q}(a)
= \sqrt{\sigma_{q,{\rm data}}(a)^2 + \sigma_{q,{\rm model}}(a)^2}
$$

from photutils ellipticity errors — **not** GALFIT `b_a_err`.
`iso_break_radius_re` = first $a/R_e$ with $|\Delta q|>3\sigma_{\Delta q}$.

**Scalars, IVW, and the flag.** Collapse $\Delta q(a)$ for $a\le 2R_e$ by
inverse-variance mean (not a plain mean: outer/noisy isophotes have larger
$\sigma_{\Delta q}$ and would dominate an unweighted average):

$$
\texttt{iso\_dq\_2re}
= \frac{\sum_i w_i\Delta q(a_i)}{\sum_i w_i},
\quad w_i=1/\sigma_{\Delta q,i}^2,
\quad
\texttt{iso\_dq\_2re\_err}=(\sum_i w_i)^{-1/2},
\quad
\texttt{iso\_dq\_2re\_sig}
= \texttt{iso\_dq\_2re}/\texttt{iso\_dq\_2re\_err}.
$$

`flag_iso_dq` is $|{\tt iso\_dq\_2re}| > 0.10$ — an absolute cut on that
collapsed difference, not a cut on `iso_dq_2re_sig` and not using GALFIT
$\sigma_q$. Other scalars: `iso_q_at_{1,2}re_*` (interpolated profile values),
`iso_dq_max_abs_inner` / `iso_frac_discrepant_inner` (how often $|\Delta q|$
exceeds $3\sigma$ inside $2R_e$), `iso_dpa_2re_*`, `iso_n_*`, strategy fields.

Also the major-axis SB profile in mag arcsec$^{-2}$ (data vs model), kept where
`intens > 0`, `mu_err > 0`, and $\sigma_\mu\le 0.3$ mag.

**Guards.** Search fail with $R_e/\mathrm{FWHM}<0.8$ → `unresolved`; tiny stamp →
`stamp_too_small`. Both benign in the aggregator.

**Population cross-check.** $q_{\rm GALFIT}$ vs $q_{\rm isophote}$ coloured by
$\mathrm{FWHM}/R_e$: isophotal $q$ should look rounder when unresolved; a
smooth trend validates both methods, outliers are real discrepancies.

**Outputs.** scalars above; full curves in `isophote_profiles.npz`.

---

### 4.6 PSF leakage — `psf_leakage.py`

**Why.** The direct test of whether the measured shapes are partly tracking the
*instrument* rather than the galaxies — the classic weak-lensing systematic. If
$q$ correlated with PSF ellipticity, the inclination distribution would be
partly an instrument artefact and the science result would be void.

**Method.** PSFEx writes Stokes ellipticity components to `psfex.xml`:

$$
e_{\rm PSF} = \sqrt{e_1^2 + e_2^2},
\qquad
\mathrm{PA}_{\rm PSF} = \tfrac12\operatorname{atan2}(e_2, e_1)
$$

The production parser in `master_run.py` keeps only the scalar
`Ellipticity_Mean`, which has no direction, so the suite carries its own parser
rather than modifying production code. An independent estimate from the
flux-weighted second moments of `proto_image.fits` — the actual stamp GALFIT
convolved with — is computed as a cross-check and used as a fallback if the XML
lacks the Stokes columns:

$$
e_1 = \frac{Q_{xx} - Q_{yy}}{Q_{xx} + Q_{yy}},
\qquad
e_2 = \frac{2Q_{xy}}{Q_{xx} + Q_{yy}}
$$

Both are converted into GALFIT's PA convention (from $+y$) and wrapped to
$[-90^\circ, 90^\circ)$, since PA is defined modulo 180.

**Population tests** (in `aggregate.py`, since a leakage test is meaningless one
galaxy at a time):

1. Spearman correlation of $q_{\rm host}$ with $e_{\rm PSF}$.
2. Distribution of $\lvert\mathrm{PA}_{\rm host} - \mathrm{PA}_{\rm PSF}\rvert$
   wrapped to $[0, 90^\circ)$, tested against uniform with a KS test. **This is
   the sharper of the two**: an amplitude correlation dilutes easily across a
   heterogeneous cohort, but alignment is unambiguous — there is no astrophysical
   reason for galaxies to align with a telescope's optics.
3. $q_{\rm host}$ against $\mathrm{FWHM}/R_e$ as the resolution control.

**Outputs.** `psf_e1`, `psf_e2`, `psf_ellipticity`, `psf_ellipticity_xml`,
`psf_ellipticity_moments`, `psf_pa_deg`, `psf_pa_moments_deg`, `psf_pa_source`,
`psf_fwhm_xml_px`, `psf_fwhm_moments_px`, `psf_nstars`, `psfex_chi2`,
`re_over_fwhm`, `fwhm_over_re`, `dpa_host_psf_deg`, `q_host`, `pa_host_deg`.

---

### 4.7 Magnitude leakage — `mag_leakage.py`

**Why.** A model can be wrong in a way $\chi^2$ tolerates but total flux does
not. Sky over-subtraction eats the outer disk; a runaway $n$ inflates it. Both
show up in the integrated magnitude even when the residual map looks clean,
because they are large-area, low-amplitude errors.

**Definition.** $\Delta m = m_{\rm GALFIT} - m_{\rm ref}$, using the
`ref_survey` / `ref_mag` columns already in `pipeline_galfit_results.csv`.

**Tests.** $\Delta m$ regressed against $R_e$, $n$, the fitted sky offset, and
$\chi^2/\nu$, each with slope, standard error and significance. The
interpretation differs by regressor:

| trend | fingerprint of |
|---|---|
| $\Delta m$ vs $R_e$ | sky error (a sky offset integrates over area, so it scales with size) — **or** aperture-mismatch in the reference survey |
| $\Delta m$ vs $n$ | runaway Sérsic index inflating the extrapolated wing |
| $\Delta m$ vs sky offset | direct confirmation of the above |
| $\Delta m$ vs $\chi^2/\nu$ | fit quality; a *flat* trend argues the others are not fit failures |

**Constraint bounds.** Bounds are parsed for component 1 from `constraints.txt`
as GALFIT was given them, and a parameter is reported as pinned when it sits
within 5% **of the bound value**, not of the allowed range:

```
_at_bound(value, (lo, hi)):
    return value <= lo * 1.05  or  value >= hi * 0.95
```

The range-relative version originally used called $n = 0.62$ "at the bound" with
$n$ free over 0.5–6.0, which it plainly is not; that inflated the flag count
from 15 to 27 on the science cut. `n_at_floor` and `n_at_ceiling` are reported
separately, because they mean opposite things.

**Caveat.** Hosts whose Phase 2 zero point failed carry an untrustworthy raw
`mag`. The check records `zp_ok` and `mag_final_source` so those can be
separated rather than silently polluting the trend fits.

**Outputs.** `dmag_ref`, `dmag_flag`, `ref_mag`, `ref_mag_err`,
`ref_sep_arcsec`, `mag_galfit`, `mag_galfit_err`, `mag_final_source`,
`n_at_bound`, `n_at_floor`, `n_at_ceiling`, `re_at_bound`, `mag_at_bound`,
`sky_fitted_adu`, `sky_offset_adu`, `sky_offset_sigma`,
`flux_outside_stamp_frac`.

---

### 4.8 AstroPhot cross-fit — `astrophot_refit.py`

**Why.** Agreement between two independent codes — different optimizer, different
likelihood implementation, different pixel-integration scheme — rules out a whole
class of implementation-specific and local-minimum failures that no internal
diagnostic can see. If GALFIT had a systematic quirk in how it handles
barely-resolved sources, every check above would inherit it; this one would not.

**Method.** AstroPhot v0.16.13 fits a Sérsic to `host_cutout.fits` with
`host_sigma.fits` as the variance source and `host_mask.fits` applied. The PSF is
`proto_image.fits` — **the same stamp GALFIT convolved with** — which makes this
a like-for-like comparison rather than a comparison of two PSF treatments.
Sérsic index free, matching the GALFIT configuration.

**AstroPhot finds its own starting point.** Seeding it with GALFIT's answer
would make agreement much cheaper than it looks; the whole value of the check is
that it converges independently.

**Sky is the one quantity held common.** GALFIT's fitted level is subtracted from
the data and AstroPhot fits the galaxy alone. The reason is a parametrization
trap: AstroPhot's flat sky is
$\log_{10}(\mathrm{flux}\,/\,\mathrm{arcsec}^{2})$, which cannot represent a
zero-or-negative background — left free it runs to $-\infty$, and the resulting
disagreement would be about that parametrization rather than about our fits. The
sky systematic is measured properly and separately by §4.4.

**Convention reconciliation.** AstroPhot shares GALFIT's PA convention (from
$+y$) but reports radians; its $R_e$ is in arcsec, not pixels. Its LM optimizer's
`res_loss()` is already the **reduced** $\chi^2$, so it is directly comparable to
GALFIT's `Chi^2/nu` — dividing by `ndf` again (as an older repo script does) is
wrong by a factor of $\nu$.

**Outputs.** `ap_q`, `ap_q_err`, `ap_pa_deg`, `ap_re_px`, `ap_re_arcsec`,
`ap_n`, `ap_chi2nu`, `ap_message`, `ap_runtime_s`, `ap_sky_fixed_adu`, and the
comparisons `dq_astrophot`, `dq_astrophot_over_q_err`, `dpa_astrophot_deg`,
`dn_astrophot`, `dre_astrophot_frac`. The normalized version matters: two codes
differing by 0.01 in $q$ is unremarkable unless GALFIT claims a 0.001 error.

---

### 4.9 Visual panel — `visual.py`

**Why.** No metric replaces looking at the residual. Top-row **data** and
**model** use an asinh stretch on sky-subtracted flux over the data’s
$[1,99]$ percentile window (grayscale, soft factor 10) so cores and envelopes
are both visible. The **residual** stays in $\sigma$ units on diverging
`RdBu_r`, clipped at $\pm5$. (Ops detail / regen / pptx:
[`VISUAL_PANELS.md`](VISUAL_PANELS.md).)

**Per-host panel** (`panel.png`, 2×3). The two black ellipses on the images
are **GALFIT geometry at $1R_e$ and $2R_e$** (scale markers), not the free
isophote fits of §4.5 (those are the curves in the bottom-centre plot; often
$\sim 20$ of them).

| | role | why this, not something else |
|---|---|---|
| top: data / model / resid$/\sigma$ | see morphology + residual | data+model: asinh 1–99%; resid: $/\sigma$, $\pm5$. Overlays = GALFIT ellipses at $1R_e$ and $2R_e$ (homothetic; $R_e$ = **semi-major**) |
| bottom left: $\mu(a)$ | major- **and** minor-axis SB cuts vs elliptical $a/R_e$, data and model | flux/scale check on both axes; for a correct ellipse the major and minor curves coincide at fixed $a$ |
| bottom centre: $q(a)$ | free isophotes on data **and** convolved model; blue dotted = **intrinsic** GALFIT $q$ | only plot that shows *where* shape disagrees (PA twin dropped — triage does not need it) |
| bottom right: Fourier $\delta q(a)$ | residual→geometry converter of §4.3 | answers “does the residual look like wrong $q$?” — orthogonal to isophotes |

Those three bottom plots are the geometry triage set: **flux** ($\mu$ on both
axes), **observed shape vs radius** (isophotes), **residual-implied $\Delta q$**
(Fourier). Header columns: geometry; \(n\)+photometry+\(\delta m\) errs; corrected
\(\chi^2/\nu\)+RFF+Lotz+2008 Gini–M20 class (`late`/`early`/`merger`; see
§4.3 Statmorph / Gini–M20); then
\(\delta q_{\rm Fou}\), \(\psi_2'\), sky \(q_\pm\), \(\Delta q_{\rm AP-G}\).

**Reading trap — GALFIT $q$ vs isophotal $q$.** The blue line is the
*deconvolved* Sérsic $q$. Isophotal $q(a)$ is measured on **pixel images after
the PSF**. For $R_e\lesssim\mathrm{FWHM}$ the PSF rounds observed isophotes, so
both $q_{\rm data}(a)$ and $q_{\rm model}(a)$ sit well above GALFIT $q$ even
when the fit is right — that gap is expected, not a smoking gun. Compare data
vs model isophotes to each other; use GALFIT $q$ only as the intrinsic
reference. Fourier $\delta q(a)$ may start at $a/R_e\gtrsim 1$ when inner
annuli lack $\ge 20$ pixels (width $\sim\mathrm{FWHM}/2$).

**$\chi^2/\nu|_{2R_e}$ vs an outer residual blob.** Localized $\chi^2/\nu$ sums
only pixels with elliptical $a\le 2R_e$. A bright unmodeled spot *outside*
that ellipse does **not** enter it. Global $\chi^2/\nu$ is usually *lower*
because thousands of sky pixels (with $\sigma$ slightly overestimated) dilute
the mean toward $<1$; the local value is higher when the galaxy region itself
is imperfect — that ordering is expected, not a bug.

**Contact sheet** (`outputs/plots/contact_sheet.png`): $\sigma$-unit residuals
for all 64 hosts on a common scale.

---

## 5. Aggregation, flags and trust tiers

`aggregate.py` collects every per-host JSON into one row per host, runs the
population-level tests (which by nature cannot live in a per-host check), and
assigns flags.

**Thresholds.** Deliberately loose and explicitly provisional until calibrated
against the real distribution of each metric:

| Threshold | Value |
|---|---|
| `rff_2re` | 0.10 |
| `fourier_dq` | 0.05 |
| `iso_dq_2re` | 0.10 |
| `dq_sky` | 0.05 |
| `dq_astrophot` | 0.05 |
| `dmag_ref` | 0.5 mag |
| `sigma_calibration_ratio` | outside (0.5, 2.0) |

**Flags.**

| Flag | Condition |
|---|---|
| `flag_sigma_miscalibrated` | $r_\sigma$ outside (0.5, 2.0) |
| `flag_rff_high` | $\lvert\mathrm{RFF}_{2R_e}\rvert > 0.10$ |
| `flag_rff_sky` | $\lvert\mathrm{RFF}_{1{-}2R_e} - \mathrm{RFF}_{1R_e}\rvert > 0.10$ |
| `flag_fourier_dq` | reliable **and** $\lvert\delta q\rvert > 0.05$ |
| `flag_fourier_dq_significant` | reliable **and** $\lvert\delta q/\sigma_{\delta q}\rvert > 3$ |
| `flag_fourier_unusable` | failed the §4.3 reliability gate |
| `flag_iso_dq` | $\lvert\Delta q_{\rm iso}\rvert > 0.10$ |
| `flag_sky_sensitive` | $\Delta q_{\rm sky} > 0.05$ |
| `flag_astrophot_disagrees` | $\lvert q_{\rm AP} - q_{\rm GALFIT}\rvert > 0.05$ |
| `flag_dmag` | $\lvert\Delta m\rvert > 0.5$ |
| `flag_param_at_bound` | $n$ or $R_e$ pinned (§4.7) |
| `flag_q_near_floor` | $q \le 0.25$ |
| `flag_unresolved` | $R_e/\mathrm{FWHM} < 1$ |
| `flag_check_error` | any check status not in {`ok`, `unresolved`, `stamp_too_small`} |

**Trust tier.** Only the four flags that speak to whether $q$ *itself* is
trustworthy count — `flag_fourier_dq`, `flag_iso_dq`, `flag_sky_sensitive`,
`flag_astrophot_disagrees`. A high RFF on a bright spiral is astrophysics, not a
bad fit, and does not belong in the tier.

| Tier | Meaning |
|---|---|
| **A** | no geometry flags |
| **B** | exactly one |
| **C** | two or more |
| **?** | a check did not complete — the host *cannot* be graded |

The `?` tier exists so an un-run host cannot silently score as clean. Note the
gating on `fourier_reliable`: without it, every barely-resolved host would be
condemned by a $\delta q$ that is arithmetic rather than measurement.

---

## 6. Data walkthrough

Suite products live under `verification/outputs/`. Re-fits, sandboxes, and
sky-protocol JSON live under `verification/Re-fits/` (see
[`VERIFICATION_README.md`](VERIFICATION_README.md)).

```
outputs/
├── per_host/<FRB>/                 64 directories
│   ├── chi2.json  rff.json  fourier.json  psf.json  mag.json
│   ├── isophote.json  sky.json  astrophot.json  visual.json
│   ├── fourier_profiles.npz        radial profiles, check 3
│   ├── isophote_profiles.npz       radial profiles, check 5
│   ├── panel.png                   the 2x3 diagnostic panel
│   └── galfit_sky_{plus,minus}/    full staged GALFIT re-runs
├── panels/
│   ├── <FRB>.png                   published production panel
│   └── <FRB>_{n1,sky,n1_sky,psf}.png   confirmed alternate legs
├── tables/
│   ├── fit_verification_metrics.csv    64 x 212
│   ├── fit_verification_flags.csv      64 x 20
│   ├── gini_m20_53.csv                 Lotz+2008 class from production statmorph
│   └── population_summary.json         116 cohort-level statistics
├── plots/
│   ├── population_diagnostics.png  mag_leakage.png
│   ├── dq_comparison.png           contact_sheet.png
├── confirmed_fit_panels.pptx       one confirmed panel per slide
└── logs/

Re-fits/<FRB>/                      staged re-fits (outside outputs/)
├── sky_protocol.json
├── n1/  sky/  n1_sky/  sandbox/
└── panel_{production,n1,sky,n1_sky}.png
```

### 6.1 Per-host JSON

One file per check, sorted keys, `NaN`/`Inf` written as `null` so the files stay
valid JSON. Every file carries `_check`, `_frb`, `_runtime_s` and `status`.
Key counts and typical runtimes, from `20240210A`:

| File | Keys | Runtime | Status vocabulary |
|---|---|---|---|
| `chi2.json` | 29 | 0.02 s | `ok` |
| `rff.json` | 22 | 0.01 s | `ok` |
| `fourier.json` | 41 | 4.1 s | `ok` |
| `psf.json` | 27 | 0.00 s | `ok` |
| `mag.json` | 33 | 0.03 s | `ok` |
| `isophote.json` | 29 | 28.4 s | `ok`, `unresolved`, `stamp_too_small`, `isophote_fit_failed` |
| `sky.json` | 43 | 3.9 s | `ok` |
| `astrophot.json` | 25 | 6.4 s | `ok` |
| `visual.json` | 7 | 1.0 s | `ok` |

Any exception is caught and written as `status: "error"` with the exception text
and an 8-frame traceback, so a failure is inspectable without re-running.

### 6.2 `fourier_profiles.npz`

33 annuli for `20240210A`. Everything the $m=2$ analysis rests on, kept so the
summary scalars can be re-derived or re-weighted without re-running.

| Key | Shape | Contents |
|---|---|---|
| `a_mid`, `a_mid_re` | (N,) | annulus midpoint in px and in $R_e$ |
| `npix` | (N,) | valid pixels in the annulus |
| `coef`, `coef_err` | (N, 9) | $[A_0, A_1, B_1, A_2, B_2, A_3, B_3, A_4, B_4]$ |
| `amp` | (5, N) | $\sqrt{A_m^2 + B_m^2}$ normalized by the local model, $m = 0..4$ |
| `dq`, `dq_err` | (N,) | calibrated $\delta q(a)$ |
| `dq_radial`, `dq_radial_err` | (N,) | analytic $\delta q(a)$, for comparison |
| `dpa_deg`, `dpa_err_deg` | (N,) | calibrated $\delta\mathrm{PA}(a)$ |
| `psi2_rad`, `psi2_err_rad` | (N,) | $m=2$ phase |
| `kernel_q_coef`, `kernel_pa_coef` | (N, 9) | the response kernels, same basis |
| `crosstalk` | (N,) | leakage between the two kernels |
| `dmda`, `model_profile` | (N,) | azimuthally averaged model and its gradient |
| `calibrated` | scalar | whether kernels or the analytic fallback were used |

### 6.3 `isophote_profiles.npz`

29 isophotes for `20240210A`, interpolated onto a common $a/R_e$ grid so data and
model are directly differenceable.

| Key | Contents |
|---|---|
| `a_re`, `sma`, `sma_model` | radius grid, and the native grids of each fit |
| `q_data`, `q_data_err`, `q_model`, `q_model_err` | axis ratio profiles |
| `q_model_native` | model $q$ before interpolation |
| `pa_data_deg`, `pa_model_deg` (+ errors) | PA profiles, GALFIT convention |
| `pa_model_native_deg` | as above, pre-interpolation |
| `dq`, `dq_err`, `dpa_deg`, `dpa_err_deg` | data − model differences |
| `intens_data`, `intens_data_err`, `intens_model` | isophotal intensities |
| `mu_data`, `mu_model`, `mu_err` | surface brightness, mag arcsec$^{-2}$ |
| `mu_valid` | where $I > 0$ and $\sigma_\mu > 0$ — the only points safe to plot |

### 6.4 `galfit_sky_{plus,minus}/`

A complete, self-contained GALFIT run per leg — 128 across the cohort:
`host_cutout.fits`, `host_sigma.fits`, `host_mask.fits`, `proto_image.fits`,
`galfit.feedme`, `constraints.txt` (sky constraint stripped), plus GALFIT's own
`out.fits`, `fit.log`, `galfit.01`–`.03` restart files and `galfit_stdout.log`.
Because the inputs are staged copies, each directory can be re-run by hand for
forensics without touching production.

### 6.5 `fit_verification_metrics.csv` — 64 × 212

One row per host. Columns are grouped by originating check, in run order; each
block includes that check's `<check>_status` column.

| Block | Cols | Examples |
|---|---|---|
| production baseline | 13 | `frb`, `in_53`, `snr_win`, `mag`, `b_a`, `b_a_err`, `re`, `n`, `pa`, `chi2nu`, `zp_ok`, `ref_survey`, `n_sersic_components` |
| check 1 — `chi2_*` | 24 | `chi2nu_global`, `chi2nu_local_2re`, `sigma_calibration_ratio`, `residual_closure`, `re_over_fwhm` |
| check 2 — `rff_*` | 18 | `rff_1re`, `rff_2re`, `rff_annulus_1_2re`, `rff_outer_minus_inner` |
| check 3 — `fourier_*`, `model_recon_*` | 36 | `fourier_dq`, `fourier_reliable`, `fourier_unreliable_reasons`, `model_recon_max_frac` |
| check 6 — `psf_*` | 21 | `psf_ellipticity`, `psf_pa_deg`, `dpa_host_psf_deg`, `fwhm_over_re` |
| check 7 — `mag_*` | 20 | `dmag_ref`, `n_at_floor`, `n_at_ceiling`, `sky_offset_sigma` |
| check 5 — `iso_*` | 26 | `iso_dq_2re`, `iso_q_at_1re_data`, `iso_break_radius_re`, `iso_same_strategy` |
| check 4 — `*_sky*` | 37 | `dq_sky`, `dq_sky_over_q_err`, `sky_sigma_adu`, `q_sky_plus`, `q_sky_minus` |
| check 8 — `ap_*` | 17 | `ap_q`, `ap_chi2nu`, `dq_astrophot`, `dq_astrophot_over_q_err` |

Production columns take precedence on merge, so a check re-reporting `mag` or
`b_a` for convenience cannot shadow the authoritative value or produce `_x`/`_y`
suffixes.

### 6.6 `fit_verification_flags.csv` — 64 × 20

`frb`, `in_53`, the 14 `flag_*` booleans of §5, `fourier_reliable`,
`n_geometry_flags`, `n_flags_total`, `trust_tier`. Deliberately small: this is
the triage table, meant to be read directly.

### 6.7 `population_summary.json` — 116 keys

Every cohort-level statistic, computed twice — suffix `_all64` and `_in53` — so
the effect of the science cut is always visible. Families:

| Prefix | Contents |
|---|---|
| `psf_q_vs_epsf_*` | Spearman $r$, $p$, $n$ for the leakage test |
| `psf_dpa_ks_*` | KS statistic, $p$, median, $n$ for PA alignment |
| `psf_q_vs_fwhm_over_re_*` | the resolution control |
| `chi2nu_local_vs_{snr,re_over_fwhm,mag}_*` | is $\chi^2/\nu$ an SNR meter? |
| `dmag_vs_{re,n,sky,chi2nu}_*` | slope, error, significance, $r$, $p$, $n$ |
| `q_vs_{isophote,astrophot}_*` | median offset, MAD scatter, $n$ |
| `<metric>_{median,p16,p84}_in53` | distribution summaries for 8 headline metrics |

### 6.8 Plots

| File | Contents |
|---|---|
| `population_diagnostics.png` | 2×3: $\chi^2/\nu$ vs SNR; RFF vs $R_e$/FWHM; $q$ vs $e_{\rm PSF}$; PA-alignment histogram against uniform; $q_{\rm AP}$ vs $q_{\rm GALFIT}$; $q_{\rm iso}$ vs $q_{\rm GALFIT}$ |
| `mag_leakage.png` | 1×4: $\Delta m$ against $R_e$, $n$, sky offset, $\chi^2/\nu$ |
| `dq_comparison.png` | 1×3 histograms of the three independent $\delta q$ handles, with medians |
| `contact_sheet.png` | 64 $\sigma$-unit residual thumbnails |

### 6.9 Reading it back

```python
import pandas as pd, numpy as np, json

m = pd.read_csv("outputs/tables/fit_verification_metrics.csv", dtype={"frb": str})
f = pd.read_csv("outputs/tables/fit_verification_flags.csv",   dtype={"frb": str})
pop = json.load(open("outputs/tables/population_summary.json"))

sci = m[m.in_53]                                    # the 53-host science cut
sci[sci.fourier_reliable == True].fourier_dq        # only where the estimator had leverage
f[f.trust_tier.isin(["C", "?"])].frb                # hosts needing a look

z = np.load("outputs/per_host/20240210A/fourier_profiles.npz")
z["a_mid_re"], z["dq"], z["dq_err"]                 # the delta-q profile
```

---

## 7. Results — first full run

Date: 2026-08-05.

All nine checks over all 64 hosts. Numbers are the 53-host science cut, quoted as
median [p16, p84].

| Metric | Median [p16, p84] | Reading |
|---|---|---|
| `chi2nu_global` | 0.889 [0.743, 1.707] | $\sigma$ mildly over-estimated, not manufactured |
| `sigma_calibration_ratio` | 0.892 [0.839, 0.959] | blank-sky MAD $\approx 0.89\times$ `host_sigma` |
| `rff_2re` | −0.002 [−0.034, +0.071] | residual consistent with noise inside $2R_e$ |
| `iso_dq_2re` | −0.026 [−0.076, +0.008] | GALFIT $q$ slightly rounder than free isophotes |
| `dq_astrophot` | −0.006 [−0.020, −0.000] | scatter 0.008 against an independent fitter |
| `dq_sky` | 0.004 [0.001, 0.055] | $q$ insensitive to a $\pm1\sigma$ sky shift |
| `dmag_ref` | −0.038 [−0.355, +0.112] | photometry ties to the reference surveys |
| `model_recon_max_frac` | 0.0013 [0.0007, 0.0090] | the analytic rebuild reproduces GALFIT to 0.1% |
| `re_over_fwhm` | 1.03 [0.53, 4.18] | **half the cohort is at or below the PSF scale** |

### 7.1 The four falsification tests all pass

1. **No PSF leakage.** $q_{\rm host}$ versus $e_{\rm PSF}$: Spearman $r = +0.04$,
   $p = 0.77$. $q$ versus FWHM/$R_e$: $r = -0.045$, $p = 0.75$.
   $\mathrm{PA}_{\rm host} - \mathrm{PA}_{\rm PSF}$ is consistent with uniform
   (KS $p = 0.85$). The axis ratios are not tracking the instrument.
2. **An independent fitter agrees.** AstroPhot reproduces GALFIT's $q$ with a
   median offset of $-0.006$ and MAD scatter $0.008$ — far below the $\sim0.05$
   that would matter for inclination. Three hosts disagree (§7.3).
3. **Sky is not driving $q$.** Median $\lvert\Delta q\rvert$ under a
   $\pm1\sigma$ sky shift is 0.004; only 9/53 exceed 0.05.
4. **Residuals carry no $m=2$ signal** where the estimator is usable: median
   $\delta q = -0.035$, $\max\lvert\delta q\rvert = 0.121$ over the 20 hosts
   that pass the reliability gate.

### 7.2 What the run exposed

**Resolution, not fit quality, is the limiting factor.** 25/53 hosts have
$R_e < \mathrm{FWHM}$. This is why the Fourier estimator is usable on only 24/53
(`too_few_annuli` and `unresolved` dominate; crosstalk is recorded but no longer
gates — hosts can fail more than one gate): a source smaller than the PSF has
no annuli to decompose.
The gate is doing its job; reporting $\delta q$ for those hosts would be
reporting noise.

**$\chi^2/\nu$ is an SNR meter, as anticipated.** Localized $\chi^2/\nu$
correlates with `snr_win` ($r = +0.41$, $p = 0.003$) and with $R_e/\mathrm{FWHM}$
($r = +0.31$, $p = 0.023$). It ranks brightness, not correctness, and must not be
used as a pass/fail gate.

**A real $\Delta m$–$R_e$ trend**, slope $-0.150$ mag arcsec$^{-1}$ at
$6.3\sigma$: larger fits are systematically brighter than the reference aperture
magnitude. Expected — the reference surveys use fixed apertures that miss
extended wings — but it is also the classic sky-error fingerprint, so it is worth
separating before the photometry is quoted. $\Delta m$ versus $\chi^2/\nu$ is
flat ($0.3\sigma$), which argues against a fit-quality origin.

**15/53 hosts have Sérsic $n$ pinned at a constraint** (9 at the 0.5 floor, 6 at
the 6.0 ceiling), and 3 have $R_e$ at a bound. Given the $n$–$q$ degeneracy this
is the strongest argument for settling the fixed-versus-free $n$ question in §9.

### 7.3 Trust tiers

**53-host cut: A = 34, B = 13, C = 6, ? = 0.** The six tier-C hosts:

| FRB | $b/a$ | $R_e$/FWHM | Issue |
|---|---|---|---|
| `20190711A` | 0.56 | 0.26 | $\chi^2/\nu = 347$; AstroPhot gives $q = 0.055$; unresolved, $n$ at ceiling |
| `20220501C` | 0.52 | 0.47 | most sky-sensitive, $\Delta q = 0.171$; AstroPhot $\Delta q = -0.30$; too compact for isophotes |
| `20190523A` | 0.73 | 0.33 | $\Delta q_{\rm sky} = 0.155$; AstroPhot $\Delta q = -0.26$; unresolved |
| `20230930A` | 0.84 | 18.7 | Phase 2 ZP failure — `mag = 8.9` is not physical, $\Delta m = -8.6$ |
| `20220825A` | 0.62 | 1.42 | isophote $\Delta q = +0.161$ alongside a reliable Fourier $\delta q = -0.121$ |
| `20220912A` | 0.63 | 1.26 | isophote $\Delta q = -0.127$; everything else clean |

The first four are the same two root causes seen from different angles — three
are unresolved, one is a known zero-point failure — rather than six independent
problems. `20230930A` should not carry a magnitude at all.

---

## 8. Deferred

Specified here so they can be switched on without redesign.

### 8.1 Profile likelihood chi-squared of q

Fix $q$ on a grid (0.05 to 0.95, finer near the minimum) and refit **all** other
parameters at each node. The result is a curve $\chi^2_{\min}(q)$, not two
points. Three things follow that nothing else provides:

1. **The real error bar.** $\Delta\chi^2 \le 1$ gives the $1\sigma$ interval,
   $\Delta\chi^2 \le 4$ the $2\sigma$. This is correct under degeneracy because
   the other parameters re-adjust at each node — precisely what a
   covariance-matrix error cannot capture when the surface is non-quadratic.
   Comparing this width with GALFIT's reported $\sigma_q$ gives the inflation
   factor, which matters directly: `master_run._mc_inclination` propagates
   `b_a_err` through $10^4$ draws to build every inclination confidence interval
   the pipeline quotes.
2. **Whether $q$ is measured at all.** A flat-bottomed curve means the data do
   not constrain $q$; a double minimum means the fit could have landed elsewhere.
3. **It subsumes multi-start convergence** along the axis of interest: a smooth
   single-minimum curve *is* the proof that there is no local-minimum problem in
   $q$.

Cost is roughly 19 re-fits per host, about an hour for the cohort. The
$n$-at-bound result in §7.2 strengthens the case for running it.

### 8.2 Injection–recovery into real images

Inject synthetic galaxies of known parameters into blank regions of the *same*
image and run the full pipeline on them, measuring bias and scatter in $q$ at
each host's true PSF, sky and SNR. Rejected as too computationally intensive:
useful scatter needs many realizations per host.

### 8.3 Multi-start convergence

Re-fit from deliberately bad initial guesses to map local minima. Rejected;
largely subsumed by the $\chi^2(q)$ scan.

---

## 9. Open questions

Raised, not settled. Not implemented.

- **Sérsic $n$ fixed versus free.** $n$ and $q$ are correlated in practice, and
  §7.2 shows 15/53 hosts with $n$ pinned at a bound. Re-fitting at $n = 1$ and
  $n = 4$ and reporting $\Delta q$ would test whether the axis ratios survive the
  degeneracy.
- **Fitted-centre drift.** Fitted centre against the `host_components.csv`
  centroid. Distinct from association: this catches a fit wandering off-centre
  during optimization even when the literature host is correct, which no residual
  metric flags.
- **The $b/a > 0.2$ selection cut.** The 53-host sample is selected on the very
  quantity whose distribution is the science result, and the cut removes the
  edge-on tail. Of the eleven excluded hosts, seven fail on $b/a$ alone
  (`20200430A`, `20220105A`, `20220717A`, `20221116A`, `20230708A`, `20240201A`,
  `20240310A`), two on magnitude alone (`20190611B`, `20211203C`), and two on
  both (`20220918A`, `20230712A`). The right treatment is to establish, per host,
  whether the low $b/a$ is a genuine edge-on system or a fit pinned near the
  constraint floor, and then either model the censoring or justify the cut. The
  all-64 run exists so that evidence is available.
- **How to present the 25 unresolved hosts.** They pass every test that *can* be
  run on them, but several tests cannot be run at all. Silently reporting them
  alongside the resolved hosts overstates the evidence; excluding them
  re-introduces a selection on size, which correlates with $q$.
