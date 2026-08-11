# Independent sky protocol (B / E / F)

Implemented in `sky_protocol.py`. Used to seed **fixed-sky** GALFIT re-fits —
not to replace production Phase sky QA.

## Purpose

GALFIT’s free sky (constant ± optional tip/tilt) is degenerate with outer disk
light and \(q\). The protocol estimates sky **independently** on a large Legacy
Survey cutout, then re-fits with that sky held fixed so \(\Delta q\) can be
judged under a controlled background.

Production `sky_final_adu` is recorded for comparison only; it **does not**
enter the consensus.

## Inputs

| input | path |
|---|---|
| Flux stamp | `large_cutouts/{FRB}_flux.fits` (multi-band; **r** plane) |
| Invvar (optional) | `large_cutouts/{FRB}_invvar.fits` |
| Host geometry | production `Output/{FRB}_all/` (`cutout_meta.json`, segmap, …) |

Analysis is restricted to a box around the host (default half-width **90″**),
not the full ~10′ stamp — avoids distant gradients dominating the estimate.

Band: **r** (`BANDS=griz` → plane index 1). Science fits are r-band.

## Methods

### B — annular growth plateau (`B_growth`)

GALAPAGOS-style: annular median vs radius; sky = median on the outer plateau
where the growth curve has flattened.

### E — mask growth (`E_mask_growth`)

Dilate the segmented source mask; take a robust mean of remaining pixels;
report the stable dilation range (Ji-style practical asymptote).

### F — empty patches (`F_empty_patches`, optional `F_plane`)

Medians of many small patches free of segmented sources inside the analysis
box. Optional local plane evaluated at the host.

## Consensus

- Collect finite estimates from B, E, F (and F_plane when used).
- **`sky_adu`** = median of those values.
- **`sky_sigma_adu`** = MAD-based scatter of the method set.
- **`agree`** ≈ method span is small vs that scatter (implementation:
  span ≲ a few × σ; see JSON `consensus.agree` / `span_adu`).
- Box stability: also evaluate at **60″ / 90″ / 120″** half-boxes
  (`box_stability` in the JSON).

`agree=False` with methods still the same order of magnitude can be acceptable
with a caveat; wild disagreement or unphysical (e.g. large negative) skies are
a reject signal for inclination work.

## Outputs

```
Re-fits/<FRB>/
  sky_protocol.json
  sky_protocol.png          # diagnostic figure (if written)
```

Key JSON fields: `estimates_adu`, `consensus.{sky_adu,sky_sigma_adu,agree,…}`,
`production` comparison / `delta_vs_production_adu`, `box_stability`,
`host_xy_full`, `halfbox_arcsec`, `band`.

## CLI

```bash
python sky_protocol.py 20181112A
# options: --halfbox-arcsec, --band, …
```

## Consumers

| consumer | how |
|---|---|
| `run_refit.py --sky-from-protocol` | reads `consensus.sky_adu`, freezes sky in feedme |
| `run_reject_grid.py` | runs protocol then legs `sky` / `n1_sky` |
| Triage notes | cite B/E/F span and whether `agree` |

## Relation to check 4 (`sky_perturb`)

Check 4 (`checks/sky_perturb.py`) is different: it holds sky at
\(S\pm\sigma_{\rm sky}\) where \(\sigma_{\rm sky}\) is the **catalog BACKGROUND
MAD**, and refits. Protocol sky is an independent absolute level for re-fits.
Both probe sky systematics on \(q\); they are complementary.

See also [`FIT_VERIFICATION_CHECKS.md` §4.4](FIT_VERIFICATION_CHECKS.md).
