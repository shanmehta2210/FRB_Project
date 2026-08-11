# Production-62 GTC science review

Tiered review of the **62** hosts in `pipeline_galfit_results.csv` for cases where
deeper GTC imaging may be scientifically justified (association ambiguity or
depth/SNR-limited fits — not intrinsic morphology or low-SNR alone).

**Does not modify** `pipeline_scripts/` or re-run the imaging pipeline.

## Refresh

```bash
python "GTC data/production_62/build_production_gtc_review.py"
```

Re-run after updating production fits, `pipeline_unphysical_fits_review.csv`, or
GTC visibility summaries under `GTC data/visibility/summaries/`.

## Files

| File | Purpose |
|------|---------|
| `literature_confident_hosts.csv` | Parsed from `Archive/notes/new_confident_hosts.txt` (Sharma/Verdi table) |
| `gtc_science_candidates.csv` | All Tier A/B flagged FRBs with metrics and visibility |
| `gtc_science_candidates.md` | Human-readable summary tables |
| `gtc_visible_intersection.csv` | Tier A∪B with `n_pass_nights >= 1` (scheduling list) |
| `build_production_gtc_review.py` | Build script |

## Tier rules (summary)

| Tier | Trigger |
|------|---------|
| **A** | `P(O) < 0.85` and/or host offset > 2″ from published localization. **Downgraded** when the FRB appears in `new_confident_hosts.txt` with `P_host ≥ 0.85` and archival GALFIT is acceptable — literature association (Sharma/Verdi) overrides AstroPath mismatch alone. |
| **B** | Degenerate/suspect fit from depth/SNR/calibration; excludes face-on pins, n/Re ceilings, secure multi-Sérsic blends |
| **C** (notes only) | Face-on geometry, secure multi-Sérsic — reported in `tier_c_note`, not promoted to GTC |

Visibility join: `GTC data/visibility/summaries/gtc_availability_by_frb_2026-06-24_to_2026-07-24.csv`.
