# Capped-Ryden experiment: bake the REX 0.8 ceiling into the intrinsic shape law

**Question (user):** instead of the ad-hoc `scaled` renormalisation, build a *genuine*
Ryden distribution whose intrinsic face-on b/a is uniform below 0.8 and falls off
exponentially above it. Does that reproduce `scaled`?

**Answer: No.** It produces a **near-uniform** cos(i) CDF, not the `scaled` curve
(max CDF difference ~0.14 from scaled), and it fits the LS b/a histogram poorly. The
reason is exactly the mechanism established in `SCALED_IS_DEGENERATE_RYDEN.md` and
`../scaled_ryden/CIRCULARITY_CHECK.md`: any intrinsic shape distribution *with real
spread* leaves a broad `P(cos i | b/a)`, and marginalising an isotropic model over the
data drives the recovered cos(i) back to uniform. The cap bounds the b/a ceiling but
does not remove the spread, so it cannot recover the deterministic `scaled` curve.

## Model (implemented)

`CappedShapeParams` in `scripts/elliptical_disk_model.py`. Intrinsic face-on axis ratio
`f = 1 - eps`:

    p(f) ∝ 1                         for 0 <= f <= cap   (uniform "otherwise")
    p(f) ∝ exp(-(f - cap)/lam)       for cap < f <= 1    (controllable fall-off)

with `cap = 0.8`. Thickness `gamma ~ N(mu_g, sig_g)` (truncated) and isotropic
orientation, projected with the Binney/Ryden equations (`generate_q_capped`). Because a
disk's projected b/a never exceeds its face-on b/a (= 1 - eps; proven: face-on q = 1-eps),
capping `f <= ~0.8` enforces observed `b/a <= ~0.8` *physically*, with genuine cos(i)
spread. `intrinsic_faceon_law.png` shows the drawn law (flat to 0.8, sharp cutoff).

## Fit to LS EXP b/a (`ba_hist_data_vs_model.png`)

Fitted `(mu_g, sig_g, lam)` with `cap=0.8` fixed (differential_evolution + chi2):

| param | value |
|-------|------:|
| mu_g (thickness) | 0.292 |
| sig_g | 0.023 |
| lam (fall-off) | 0.005 (=> nearly hard cap at 0.8) |
| chi2 (b/a hist) | 3.49e6 |

The fit is **poor**: LS b/a peaks at ~0.42-0.47, but a uniform face-on law + isotropic
projection peaks too low (~0.32) and overproduces low-b/a galaxies. For reference the
free 4-param Ryden refit reached chi2 ~1.97e6, so the uniform-face-on assumption is a
markedly worse description of LS. `lam` collapses to ~0.005, i.e. the fit wants an almost
hard cap at 0.8.

## cos(i) result vs `scaled` (`cdf_compare.png`, `cdfs/mag*.png`)

| mag | capped-Ryden median | scaled median | max CDF diff (capped vs scaled) |
|----:|--------------------:|--------------:|--------------------------------:|
| 20 | 0.536 | 0.616 | 0.137 |
| 21 | 0.517 | 0.573 | 0.136 |
| 22 | 0.495 | 0.537 | 0.149 |

The capped-Ryden CDF hugs the **uniform diagonal** (medians ~0.5); `scaled` is the
curved S. They differ by KS ~0.14 - not "very similar".

## Why, and what *would* reproduce scaled

At a fixed observed b/a, the uniform face-on law admits many `(eps, cos theta)` combos
(a face-on elongated disk vs an inclined rounder disk), so `P(cos i | b/a)` is broad.
Averaged over the sample under isotropy, the recovered cos(i) marginalises back toward
uniform - the same degeneracy that makes any fitted Ryden model ~uniform.

`scaled` is recovered only when the intrinsic face-on b/a is a **delta at the cap**
(eps = 0.2 for every disk), which is the *degenerate* elliptical-disk limit with the
face-on edge at 0.8 (see `../scaled_ryden_fixed/SCALED_IS_DEGENERATE_RYDEN.md`). A delta
kills the spread and restores the deterministic 1:1 map; a *uniform* law keeps the spread
and gives ~uniform. So along the "width of the intrinsic shape law" axis:

    delta at 0.8  ->  scaled (deterministic curve)
    uniform < 0.8 ->  ~uniform (this experiment)

There is no intermediate genuine-distribution setting that looks like `scaled`; the curve
is a property of the zero-spread limit, not of the ceiling.

## Takeaway

The capped model is physically reasonable for enforcing the REX ceiling and gives a
sensible median (~0.5), but it does **not** emulate `scaled`; it emulates the isotropic
uniform null. Practically this means the two defensible endpoints remain:
`scaled` (= degenerate delta, max power) and a genuine shape fit (= ~uniform null).
Reproduced by `scripts/build_ryden_capped_null.py`; params in `fit_params.json`.
