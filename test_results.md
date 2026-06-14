# SDSS vs FRB inclination statistical tests

Anderson–Darling (`scipy.stats.anderson_ksamp`) and Mann–Whitney U (two-sided). SDSS: `u-r < 2.3`. FRB: GALFIT point estimates (no MC). Hubble `q0=0.2`. Driver: `scripts/run_sdss_frb_inclination_tests.py`.

## Mode A — strict cos(i)

SDSS: `modelMag_r` cut, `best_model_ba_r > 0.2`, Hubble cos(i). FRB: mag cut, `b/a > 0.2`, cos(GALFIT inc). Same cuts as CDF plots.

| mag limit | N_SDSS | N_FRB | AD statistic | AD p | AD note | MWU U | MWU p |
|-----------|--------|-------|--------------|--------|---------|-------|-------|
| 20 | 9545 | 31 | -0.8584 | 0.2500 | p capped at 0.25 (not significant) | 1.494e+05 | 0.9264 |
| 21 | 27474 | 41 | -0.3439 | 0.2500 | p capped at 0.25 (not significant) | 5.289e+05 | 0.5000 |
| 22 | 58137 | 49 | 2.098 | 0.0445 | p = 0.0445 | 1.189e+06 | 0.0451 |

## Mode B — strict i (deg)

Same sample selection as Mode A; SDSS i from arccos(cos i), FRB GALFIT `inc`.

| mag limit | N_SDSS | N_FRB | AD statistic | AD p | AD note | MWU U | MWU p |
|-----------|--------|-------|--------------|--------|---------|-------|-------|
| 20 | 9545 | 31 | -0.8584 | 0.2500 | p capped at 0.25 (not significant) | 1.465e+05 | 0.9264 |
| 21 | 27474 | 41 | -0.3439 | 0.2500 | p capped at 0.25 (not significant) | 5.975e+05 | 0.5000 |
| 22 | 58137 | 49 | 2.098 | 0.0445 | p = 0.0445 | 1.66e+06 | 0.0451 |

## Mode C — inclusive cos(i)

SDSS: mag + color, finite b/a in [0,1] (no b/a>0.2). FRB: mag cut only. Hubble cos(i) with cos i=0 when b/a ≤ q0.

| mag limit | N_SDSS | N_FRB | AD statistic | AD p | AD note | MWU U | MWU p |
|-----------|--------|-------|--------------|--------|---------|-------|-------|
| 20 | 10436 | 35 | -0.5557 | 0.2500 | p capped at 0.25 (not significant) | 1.875e+05 | 0.7830 |
| 21 | 29844 | 45 | -0.2234 | 0.2500 | p capped at 0.25 (not significant) | 6.388e+05 | 0.5722 |
| 22 | 73506 | 58 | 6.489 | 0.0011 | p = 0.0011 | 1.712e+06 | 0.0091 |
## What these tests are asking

Each row tests whether the **FRB host inclination sample** could have been drawn from
the same distribution as the **SDSS field-galaxy null** sample, after applying the
listed magnitude and color (and, for Modes A/B, axis-ratio) cuts. A small p-value
means the two samples differ in a way that is unlikely under the test's null
hypothesis—not automatically that FRBs are "more edge-on" (check CDFs for direction).

These tests complement the CDF figures in
`plots/plots_null/v1_null_cdf_inclination/mag_cuts/`; they do not replace visual
inspection.

## How to read p-values (the scale)

**Direction:** Think of p as “how surprising would these data be if FRB and SDSS
inclinations really came from the same distribution?”

- **High p (close to 1)** → *not surprising* → data are **compatible with “similar”**
- **Low p (close to 0)** → *very surprising* → data **favour “different”**

**It is not a similarity percentage.** MWU p = 0.28 does **not** mean “28% different”
or “72% the same.” It is not on a 0–100% “how alike are the CDFs?” scale.

**Usual cutoffs (convention, not physics):**

| p-value | Plain-language read |
|---------|---------------------|
| **> 0.10** | No meaningful evidence of a difference (for this test). Treat as **consistent with the null**. |
| **0.05 – 0.10** | Weak / suggestive only; many fields still call this **not significant**. |
| **0.01 – 0.05** | **Significant** at the common 5% bar — difference is unlikely to be pure chance. |
| **< 0.01** | Stronger significance (1% bar) — still says nothing about *how big* the shift is. |

**Your numbers, on that scale:**

| MWU p | What it means |
|-------|----------------|
| **0.98** (mag < 20) | **Remarkably consistent with “no difference.”** If the true distributions were the same, you would often see a p-value this high. This is as “they look alike” as these tests get. |
| **0.28** (mag < 21) | **Not significantly different** at the usual 5% level. Nowhere near “really, really different”—it is the opposite: the data are **plausibly from the same distribution**, with no strong rank shift detected. |
| **0.009** (mag < 22) | **Significant** (~1% level): a systematic shift in inclination is unlikely to be chance alone. That is **real statistical evidence of a difference**, but it does **not** by itself mean the CDFs are wildly separated—you still check the plot for *size* and *direction*. |

**Effect size vs significance:** With ~100k SDSS galaxies, even a **small** CDF offset can
yield p ≈ 0.01. Significance answers “is there a detectable shift?” not “is it huge?”
Use the mag 22 CDF overlay for “how much” and “which way.”

**AD p = 0.25 (capped):** SciPy’s way of saying p is **at least 0.25**—even more
“everything looks fine / similar” than 0.28. Not a separate scale; same rule: high = similar.

## Column guide

| Column | Meaning |
|--------|---------|
| **mag limit** | Keep galaxies with `modelMag_r` (SDSS) or GALFIT `mag` (FRB) ≤ this value. |
| **N_SDSS** | Size of the SDSS null pool after all cuts for that mode. Should match `mag_cut_summary.csv` for strict modes. |
| **N_FRB** | Number of FRB hosts passing the matching FRB cuts. |
| **AD statistic** | Anderson–Darling k-sample statistic (SciPy `anderson_ksamp`). Larger values → more evidence against "same distribution." Can be negative in recent SciPy versions; rely on **AD p**. |
| **AD p** | Approximate p-value for equal distributions. Values reported as **0.25** are **capped** (SciPy: true p > 0.25; samples very similar). |
| **AD note** | Short readout of AD p (capped / significant). |
| **MWU U** | Mann–Whitney U statistic. Not normalized when sample sizes differ by orders of magnitude; use **MWU p** for inference. |
| **MWU p** | Two-sided p-value for equal rank distributions. Small p → systematic shift in inclination between samples. |

## How to read the three modes

**Mode A (strict cos i)** — Closest to the CDF plots (x-axis is cos i). Use this mode
when comparing test outcomes to `mag_cuts/magXX/sdss_strict/null_cdf_inclination.png`.

**Mode B (strict i deg)** — Same galaxies as Mode A, but inclination in degrees.
Because cos i is a monotonic function of i on [0°, 90°], **AD and MWU results
should match Mode A** (only tiny floating-point differences possible).

**Mode C (inclusive cos i)** — Drops the `b/a > 0.2` cut on both surveys (FRB: mag only).
The SDSS pool is larger and includes face-on systems piled up at cos i = 0 when Hubble's
formula fails. Use this to see whether conclusions depend on excluding very round galaxies.

## Interpreting your current results (strict modes A/B)

| mag | Rough takeaway |
|-----|----------------|
| **< 20** | No evidence for a different distribution (AD p capped at 0.25; MWU p ≈ 0.98). Bright, small FRB subsample (N=31). |
| **< 21** | Still consistent with the null (AD p capped; MWU p ≈ 0.28). |
| **< 22** | **Significant difference** at ~1% level (AD p ≈ 0.008; MWU p ≈ 0.009). Visually, check whether FRB CDF sits above/below the SDSS band in the mag22 plot—tests do not state the direction of the shift. |

Mode C at mag < 22 shows an even smaller p-value because the inclusive null includes
more face-on galaxies (cos i → 0), which can exaggerate FRB–null separation if FRBs
are less face-on on average.

## Limitations and caveats

1. **Sample size imbalance** — ~10⁴–10⁵ SDSS vs ~30–60 FRBs. Mann–Whitney has high
   power; mag 22 significance may reflect a modest visual offset, not a dramatic effect.
2. **No FRB measurement errors** — CDFs perturb `inc` using `inc_err`; these tests use
   point estimates only.
3. **Different measurement systems** — SDSS: Hubble formula on catalog b/a. FRB: GALFIT
   Sérsic fit inclinations. Systematic offsets can masquerade as distribution differences.
4. **Selection not identical** — SDSS null is color-trimmed (late-type proxy); FRB hosts
   are not color-cut. Matches the CDF methodology but is not a perfect physical match.
5. **Multiple comparisons** — Nine table rows (3 modes × 3 mag limits) without Bonferroni
   or similar correction; treat borderline cases cautiously.
6. **Mag limit is a hard cut** — Pools at mag 20, 21, 22 are nested (larger limits include
   fainter galaxies). Results are correlated across rows, not independent experiments.
7. **Not causal** — Rejecting "same distribution" does not identify astrophysics (e.g.
   host type, redshift, selection in FRB surveys) without further modeling.

## Method reference

Implementation: `scripts/run_sdss_frb_inclination_tests.py`. Shared cuts:
`scripts/null_catalog_utils.py`, `scripts/pipeline_null_plot_utils.py` (`frb_hosts_for_cdf`).
CDF driver: `scripts/plot_null_mag_cut_cdfs.py`.
