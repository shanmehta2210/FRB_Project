# Materials

Copies of the inclination-test results. The working copies stay in `plots/plots_null/v2/frb_vs_hsc_confirmed50/`. Scripts stay in `scripts/`.

| Path | Contents |
|---|---|
| `statistics/STATISTICAL_FINDINGS.md` | Confirmed-50 tests, Protocol B error draws, and the 2024 hosts remeasured (22 of 23) |
| `statistics/TEST_RESULTS.md` | Point-estimate test dump, including the appendix tests |
| `statistics/confirmed50_q.csv` | Winning-leg $q$, $\sigma_q$, and $\cos i$ for the 50 |
| `statistics/nature23_tests.json` | Machine-readable output of the 2024-host run |
| `figures/frb_vs_hsc_protocolB.png` | CDF with sky $\oplus$ $q$ error, then $5^\circ$ on $i$ |
| `notes/GALFIT_Q_ERROR_LITERATURE.md` | Why the $5^\circ$ floor sits on $i$, and the citations |
| `notes/protocol_B_figure.md` | Sample definition and the two error protocols |
| `appendix/STATMORPH_MORPHOLOGY.md` | statmorph on the 51: GALFIT/statmorph/AstroPhot inclination RMS (basis of the $5^\circ$ floor), $F(G,M_{20})$ disk classification, per-host table; `statmorph_51.csv`, `statmorph_51_table.tex` |
| `appendix/CONTROLS_Q0_NRANGE.md` | $q_0 \in \{0.10, 0.15, 0.20, 0.25\}$ and HSC Sérsic-$n$ window controls with KS/MWU $p$-values; `controls_q0.csv`, `controls_nrange.csv`; figures `figures/figA_q0.*`, `figures/figA_nrange.*` |
| `host_properties/` | Literature sweep of the host properties (51 + 6 pre-cut hosts), Table-1 clone and extended tables; see `SWEEP_PROTOCOL.md`, `AUDIT_PROTOCOL.md` |

`statmorph_51_full.py` + `appendix_statmorph_section.py` build the statmorph appendix; `appendix_controls_tables.py` builds the $q_0$ / $n$-window tables. `run_frb_hsc_confirmed50_tests.py` writes the point-estimate file. `run_frb_hsc_error_mc.py` is the 2,000-draw Protocol B rerun of the 50. `run_frb_hsc_nature23_tests.py` reruns the same tests on the 2024 hosts.
