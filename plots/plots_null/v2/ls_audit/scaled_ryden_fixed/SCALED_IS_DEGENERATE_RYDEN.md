# `scaled` is the elliptical-disk (Ryden) model in its degenerate limit

Proof + the exact physical approximations, and why full Ryden can never "look" like
`scaled` unless you make it degenerate. Figure: `ryden_equals_scaled.png`.
Script: `scripts/ryden_reduces_to_scaled.py`.

## 1. Analytic reduction

Ryden (2004) projected axis ratio (Binney 1985, eqs 12-15) for a triaxial disk with
thickness gamma=C/A, face-on ellipticity eps=1-B/A, viewed at (theta, phi):

    e2 = eps(2-eps)
    A  = (1 - e2 sin^2 phi) cos^2 theta + gamma^2 sin^2 theta
    B  = 4 e2^2 cos^2 theta sin^2 phi cos^2 phi
    C  = 1 - e2 cos^2 phi
    q^2 = (A + C - sqrt((A-C)^2 + B)) / (A + C + sqrt((A-C)^2 + B))

Set **eps = 0** (a perfectly circular disk). Then e2 = 0, so B = 0 and the phi terms
drop out:

    A = cos^2 theta + gamma^2 sin^2 theta,   C = 1
    q^2 = A / C = cos^2 theta (1 - gamma^2) + gamma^2

which inverts to a single, deterministic value (no azimuthal scatter):

    cos(theta) = sqrt( (q^2 - gamma^2) / (1 - gamma^2) ).

That is **exactly the Hubble formula with q0 = gamma**. The elliptical-disk model
therefore collapses onto the Hubble deprojection iff:

| # | approximation | parameter | meaning |
|---|---------------|-----------|---------|
| 1 | eps = 0            | mu -> -inf (we use mu=-12) | intrinsic disks are **perfectly circular** |
| 2 | sigma_gamma = 0    | sig_g = 0 (we use 1e-4)    | a **single universal thickness** |
| 3 | gamma = q0 = 0.2   | mu_g = 0.2                 | thickness = assumed edge-on axis ratio |

`scaled` adds one non-shape step - the **REX renormalisation** cos i -> cos i / H(0.8),
i.e. relocate the face-on edge to the empirically-observed b/a=0.8 cap (see
`REX_AND_ELLIPTICAL_DISK.md`; the EXP b/a distribution truncates at the 99.9th pct
= 0.803). H(0.8) = sqrt((0.8^2 - 0.2^2)/(1 - 0.2^2)) = 0.7906.

## 2. Numerical confirmation

- Degenerate-Ryden sampled cos(theta) vs Hubble(q): median |diff| = 1.8e-5 (exact).
- Degenerate-Ryden + /0.8 renorm vs `scaled`, LS EXP: max CDF difference
  **0.005 / 0.004 / 0.003** at mag 20 / 21 / 22; medians agree to 0.001.

So: **scaled == elliptical-disk model with {eps=0, sigma_gamma=0, gamma=q0=0.2} + REX /0.8.**

## 3. What `scaled` physically assumes (and why it "looks good")

`scaled` treats every galaxy as an **identical, perfectly circular, single-thickness
oblate disk**, and asserts the roundest observable object (b/a=0.8, set by REX) is
exactly face-on. It has **zero free shape distribution**, so nothing can absorb the
inclination signal - the b/a -> cos i map is a fixed, monotonic 1:1 curve. That is
precisely why it is stable and preserves contrast.

Full Ryden with real eps > 0 does the physically correct thing: a given b/a maps to a
*distribution* of cos i (azimuthal + ellipticity scatter). Consequences, all observed:
- **Fitting** eps to LS lets the shape distribution soak up the REX-truncated b/a
  structure -> cos i collapses to ~uniform (`scaled_ryden/CIRCULARITY_CHECK.md`).
- **Fixed literature** eps reads the REX round-deficit as real inclination -> too
  edge-on (~0.41).
- **Degenerate** eps=0 is the *only* Ryden setting that reproduces scaled - because it
  throws away the very feature (intrinsic ellipticity) that makes it "Ryden".

## 4. Literature on the REX / linear-stretch question

- No paper publishes a "renormalise-to-max-b/a" de-REX correction. DESI-LS morphology
  work (e.g. the LS4 star-galaxy paper) states shape params "are only calculated for
  sources best fit to an extended profile" and **excludes** them; REX shapes are treated
  as missing, not corrected.
- The inclination/Tully-Fisher community handles unreliable photometric inclinations by
  **discarding low-inclination (face-on) galaxies** (i < 40-45 deg), using **kinematic
  inclinations**, or **statistical (Lucy) rectification** that bypasses per-galaxy
  b/a -> i (typical q0 = 0.13-0.2).

The linear stretch is therefore an ad-hoc device; its justification here is (i) the
empirical b/a=0.8 truncation, and (ii) the exact reduction above - scaled is not
"un-physical", it is the maximally-constrained (circular, single-thickness) member of
the same elliptical-disk family, with the face-on edge fixed by REX.

## 5. Two defensible pipelines for the FRB comparison

A. **Simple scaled** (= degenerate Ryden + /0.8). Apply the identical map to FRB host
   b/a and to the LS control; two-sample KS/AD. Max power, strongest shape assumption.
B. **Full shape fit** (fit eps, gamma to the control, sample P(cos i|b/a)); apply the
   SAME fitted model to FRB hosts and compare. Because a good fit drives the control to
   ~uniform, this is effectively "FRB vs uniform" but with self-consistent error
   propagation. Weakest assumptions, least power.

Recommend carrying **A as primary** (physically the circular-disk limit, REX-anchored)
and **B as a cross-check**.
