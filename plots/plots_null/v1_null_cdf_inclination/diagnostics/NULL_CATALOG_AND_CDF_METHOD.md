# Null catalog and inclination CDF methodology (v1) — stub

**Canonical document (full provenance, cut order, funnel, columns):**  
[`NULL_CATALOG_DATA_AND_INCLINATION_AUDIT.md`](../../../../NULL_CATALOG_DATA_AND_INCLINATION_AUDIT.md)  
(see **Appendix A** for the former contents of this file).

**Local artifacts (stay here):**

| Path | Role |
|------|------|
| [`cut_funnel.csv`](cut_funnel.csv) | Machine-readable funnel counts |
| [`sdss_cut_evolution/AUDIT_VERIFICATION.md`](sdss_cut_evolution/AUDIT_VERIFICATION.md) | Pass/fail check of root-audit claims vs v1 plots |
| [`legacy_morphology/`](legacy_morphology/) | Legacy morphology diagnostics |
| [`sdss_profile_winner/`](sdss_profile_winner/) | lnL exp-winner validation |

Regenerate funnel: `python scripts/audit_and_plot_null_v1_diagnostics.py --full`.

Pre-morphology-cut plots: `Archive/plots_null_pre_morphology_cut/`.
