# scaled_ryden_fixed — Ryden (2004) frozen shapes + Unterborn A1 dust

## The REX trilemma

REX (round-object excision) removes near-round disks from LS Tractor EXP, so LS
has a deficit of face-on/round galaxies. This forces a three-way choice:

1. **Refit shapes to LS** -> the shape distribution absorbs all the b/a
   structure and cos(i) collapses to ~uniform (circular; isotropy in = out).
   See `scaled_ryden/CIRCULARITY_CHECK.md`.
2. **Use SDSS Ryden (2004) shapes directly** (`pure_lit`) -> the REX round-
   deficit is misread as real inclination, pushing the median too edge-on
   (~0.41 after dust).
3. **Use SDSS Ryden shapes + assume b/a=0.8 is face-on** (`ceiling08`;
   de-REX via b/a -> b/a/0.8) -> curved, median ~0.5, matching ad-hoc scaled.

Assuming the 0.8 ceiling is the price of directly using literature Ryden shapes
on REX-truncated data. Per-galaxy cos(i) is drawn from P(cos i | b/a) of the
FROZEN model (never refit); dust via the Unterborn face-on mag re-cut.

## Median cos(i) (Unterborn A1 dust)

| mag | pure-lit + A1 | ceiling-0.8 + A1 | scaled + A1 (ref) |
|----:|--------------:|-----------------:|------------------:|
| 20 | 0.438 | 0.5747 | 0.5271 |
| 21 | 0.4122 | 0.5524 | 0.5094 |
| 22 | 0.4058 | 0.5478 | 0.5062 |

Apply the identical frozen sampler + same b/a treatment (fit_params.json) to
FRB host b/a and compare FRB-vs-null (two-sample KS/AD), not FRB-vs-uniform.
