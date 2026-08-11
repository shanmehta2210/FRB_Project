# COSMOS Zurich morphology catalog — notes

**Catalog:** `data/cosmos_morph_zurich_1.0.tbl` (COSMOS Zurich morphology v1.0)  
**Used for:** HST arm of the COSMOS b/a vs magnitude audit (PSF-safe reference vs SDSS)

---

## Column map (production)

| Role | Zurich column | Derived / output |
|------|---------------|------------------|
| ID | `SequentialID` | preserved |
| Position | `RA`, `DEC` | `ra`, `dec` (degrees J2000) |
| Magnitude | `ACS_MAG_AUTO` | ACS I-band AUTO magnitude (AB) |
| Ellipticity | `ELL_GIM2D` | `b_a = 1 - ELL_GIM2D` |
| Size | `R_0P5_GIM2D` | `Re_arcsec` (GIM2D half-light radius, arcsec) |
| Galaxy class | `ACS_MU_CLASS` | require `== 1` |
| Star/galaxy | `STELLARITY` | require `== 0` (galaxy) |
| Quality | `JUNKFLAG` | require `== 0` |
| GIM2D flux | `FLUX_GIM2D` | require `> 0` (successful fit) |
| Morphology type | `TYPE` | preserved (informational) |

**Sentinel:** `-999999` marks missing/invalid values; treated as NaN in [`scripts/build_cosmos_hst_zurich_catalog.py`](../scripts/build_cosmos_hst_zurich_catalog.py).

---

## Depth limit

GIM2D structural parameters are reliable to **I_AB ≲ 22.5** (Sargent et al. 2007; see also CATALOG_DECISIONS §2).  
Pass-1 overlay plots include a dedicated **mag ≤ 22.5** panel alongside the full SDSS curve.

---

## Strict cut funnel (ACS footprint)

Applied in order by `build_cosmos_hst_zurich_catalog.py`:

1. ACS/WFC mosaic box (RA 149.43–150.80°, Dec 1.57–2.89°)
2. `ACS_MU_CLASS == 1` (galaxy)
3. `STELLARITY == 0`
4. `JUNKFLAG == 0`
5. Finite `ACS_MAG_AUTO`, `FLUX_GIM2D > 0`, `ELL_GIM2D ∈ [0, 1)`
6. `b_a > 0.2`, finite `Re_arcsec > 0`

**Production pool (this run):** N ≈ 34,232 galaxies  
Funnel log: `cosmos_hst_cut_funnel.csv`

---

## Comparison policy (this audit)

- **No colour cut** on HST or SDSS (ablation vs v2 full-sky strict null with u−r < 2.3)
- Compare **trend shapes** of median b/a vs mag, not absolute magnitudes (ACS I vs SDSS r)
- Re slices are **survey-native only** — do not compare absolute Re across HST and SDSS

---

## References

- Sargent et al. 2007 — GIM2D COSMOS morphology
- Scoville et al. 2007 — ACS mosaic footprint
