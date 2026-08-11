# Production-62 GTC science review

Last built: **2026-06-24 11:36:07 UTC**
(`python "GTC data/production_62/build_production_gtc_review.py"`)

Machine-readable: `GTC data/production_62/gtc_science_candidates.csv`

## Summary

| Metric | Count |
|--------|------:|
| Production hosts (`pipeline_galfit_results.csv`) | 62 |
| GTC-visible (Jun 24–Jul 24 2026) | 27 |
| Tier A or B flagged | 7 |
| — Tier A (association) | 2 |
| — Tier B (depth/SNR degenerate) | 5 |
| Flagged and GTC-visible | 3 |
| Clean (not flagged) and GTC-visible | 24 |
| Literature confident hosts (`new_confident_hosts.txt`) | 109 |
| — overlap 46-host cohort | 44 |
| Tier A downgraded (lit host + good GALFIT) | 0 |

## Tier definitions

- **Tier A**: `P(O) < 0.85` and/or AstroPath host > 2″ from published localization.
  Downgraded when FRB is in `Archive/notes/new_confident_hosts.txt` with
  `P_host ≥ 0.85` and archival GALFIT is acceptable (AP mismatch alone
  is not a GTC driver when Sharma/Verdi association is already secure).
- **Tier B**: Degenerate/suspect archival fit driven by depth/SNR/calibration,
  excluding intrinsic morphology (face-on pins, n/Re ceilings, secure multi-Sérsic blends).
- **Tier C** (notes only): face-on pins, secure multi-Sérsic — not GTC drivers.

## Tier A + B candidates

| FRB | Tiers | Reasons | P(O) | Lit P_host | Offset″ | mag | χ²/ν | GTC nights | Visible |
|-----|-------|---------|------|------------|---------|-----|------|------------|---------|
| 20220307B | B | degenerate_depth_snr | 0.999 | 0.983 | 0.75 | 19.86 | 0.97 | 31 | yes |
| 20220912A | B | degenerate_depth_snr | 0.999 | 0.95 | 0.00 | 19.52 | 0.94 | 31 | yes |
| 20220319D | A | host_offset_gt2as | 0.991 | — | 9.21 | 15.88 | 1.55 | 2 | yes |
| 20240304A | A | low_P_O | 0.261 | — | 0.04 | 21.42 | 0.71 | 0 | no |
| 20211127I | B | degenerate_depth_snr | 0.996 | 1 | 0.95 | 14.95 | 0.49 | 0 | no |
| 20211212A | B | degenerate_depth_snr | 0.996 | 1 | 0.65 | 16.51 | 6.34 | 0 | no |
| 20220725A | B | degenerate_depth_snr | 1 | 1 | 0.10 | 18.12 | 3.72 | 0 | no |

## GTC-visible intersection (scheduling list)

See `gtc_visible_intersection.csv`.

- **20220307B** [B] degenerate_depth_snr (31 nights)
- **20220912A** [B] degenerate_depth_snr (31 nights)
- **20220319D** [A] host_offset_gt2as (2 nights)
