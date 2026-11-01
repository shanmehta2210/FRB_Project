# Appendix figures

Source: `scripts/appendix_figures.py`. Style inherited from Figure 1 (`plot_paper_inclination_figures.py` rcParams: serif/STIX, 9 pt axis labels, 8 pt ticks, inward ticks; control `#1F77B4`, FRB hosts `#D62728`, reference lines black). The two side-by-side panels are saved at exactly 3.4 × 2.8 in with identical axes margins (`subplots_adjust`, no bbox cropping) so they line up in one row; the multi-panel figures use `bbox_inches='tight', pad_inches=0.02`. PNGs at 200 dpi. $p_{\rm m}$ annotations are read from `robustness_results.json` (run of 10,000 draws).

## figA_control_hist (3.4 × 2.8 in, left panel)

Histogram of $\cos i$ for the 24,450 control galaxies (adopted pool: `goodfits_flag = 1`, $r \le 22$, $0.4 < n < 1.5$, $q > 0.2$, $q_0 = 0.2$), 20 bins of width 0.05, `density=True` (uniform = 1). Filled blue ($\alpha = 0.5$) with a darker blue step outline; dashed black line at 1 labelled "Random orientation". Top axis: inclination at 90, 75, 60, 45, 30, 0 deg via `top_inclination`. y from 0 to 1.15 × tallest bar (1.42). No median line, quartile edges or FRB hosts. Median $i = 54.9^\circ$, $f(i > 60^\circ) = 39.8\%$ go in the caption.

## figB_gini_m20 (3.4 × 2.8 in, right panel)

$G$ and $M_{20}$ from `gini_m20_56.csv`. Reliable = `flag = 0`, $\langle S/N\rangle_{\rm pix} \ge 2.5$, $R_e/{\rm FWHM} \ge 1.0$: 33 filled red circles; 23 unreliable as open red circles (s = 18). $x$ reversed, -0.5 to -3.0; $y$ 0.30 to 0.70. Lotz et al. (2008) boundaries (eqs. 4–6 of that paper, checked): merger line $G = -0.14\,M_{20} + 0.33$ across the full range; E/S0/Sa–Sb/Irr line $G = 0.14\,M_{20} + 0.80$ drawn only for $M_{20} \le -1.68$ (the lines meet at $G = 0.565$). Grey italic region labels. Only the four hosts outside Sb–Irr are labelled: 20190608B, 20220319D, 20220509G, 20231120A. Legend lower right (fewest points).

## Control-variant CDFs: figA_ctrl_sersic, figA_ctrl_mag, figA_ctrl_color

These replace `figA_control_variants`. One control variant per panel, in the style of `figC_q0`: uniform diagonal (dashed black, lw 0.6); the 68% band of 10,000 samples drawn from *that panel's* control pool, each of that panel's FRB sample size, on the 500-point $\cos i$ grid (light blue); that variant's control CDF (solid blue, lw 1.0); the FRB step CDF (red, lw 1.5, `where='post'`, legend gives the panel's $N$); "AD $p_{\rm m}$" lower right from `robustness_results.json` (tables A2, A3, A5; the adopted configuration reproduces the main-text run, seed 42, AD $p_{\rm m} = 0.54$); short title above. `sharex`/`sharey`; $y$ label on the left column, $x$ label on the bottom row, top inclination axis on the top row only.

`figA_ctrl_sersic` (7.1 × 5.6 in, 2 × 2): $0.4 < n < 2.0$ | $0.4 < n < 2.5$ / $n < 2.5$ | no morphology cut.

`figA_ctrl_mag` (7.1 × 5.6 in, 2 × 2): left column the limit on the control only (all 56 hosts), right column on both samples (the band then uses samples of 45 and 32); top row $r \le 21$, bottom row $r \le 20$.

`figA_ctrl_color` (7.1 × 2.6 in, 1 × 3): star-forming (`quiescent_flag = 0`) | $g - i < 1.2$ | $g - i \ge 1.2$.

| figure | panel | $N_{\rm FRB}$ | $N_{\rm pool}$ | AD $p_{\rm m}$ |
|---|---|---|---|---|
| figA_ctrl_sersic | $0.4 < n < 2.0$ | 56 | 31,212 | 0.64 |
| figA_ctrl_sersic | $0.4 < n < 2.5$ | 56 | 36,700 | 0.68 |
| figA_ctrl_sersic | $n < 2.5$ | 56 | 37,670 | 0.69 |
| figA_ctrl_sersic | No morphology cut | 56 | 52,682 | 0.74 |
| figA_ctrl_mag | $r \leq 21$ (control only) | 56 | 7,423 | 0.74 |
| figA_ctrl_mag | $r \leq 21$ (both samples) | 45 | 7,423 | 0.44 |
| figA_ctrl_mag | $r \leq 20$ (control only) | 56 | 1,426 | 0.84 |
| figA_ctrl_mag | $r \leq 20$ (both samples) | 32 | 1,426 | 0.75 |
| figA_ctrl_color | Star-forming | 56 | 23,331 | 0.57 |
| figA_ctrl_color | $g - i < 1.2$ | 56 | 11,720 | 0.63 |
| figA_ctrl_color | $g - i \geq 1.2$ | 56 | 12,730 | 0.45 |

## figC_q0 (7.1 × 6.2 in, 2 × 2)

Revision of the old `figA_q0`: per panel the $q > q_0$ cut is applied to both samples and $\cos i$ recomputed with that $q_0$. Added (1) the 68% band of 10,000 draws from that panel's control pool, each of that panel's FRB sample size; (2) the AD $p_{\rm m}$ from the C1 table, lower right; (3) the top inclination axis; legend entry renamed "Control sample". Panel title replaced by the in-panel $q_0$ label.

| $q_0$ | $N_{\rm FRB}$ | $N_{\rm pool}$ | AD $p_{\rm m}$ |
|---|---|---|---|
| 0.10 | 57 | 24,976 | 0.57 |
| 0.15 | 57 | 24,870 | 0.59 |
| 0.20 | 56 | 24,450 | 0.55 |
| 0.25 | 52 | 23,699 | 0.21 |

LaTeX placement (two panels in one row):

```latex
\begin{figure*}
\centering
\includegraphics[width=3.4in]{figA_control_hist.pdf}\hfill
\includegraphics[width=3.4in]{figB_gini_m20.pdf}
\caption{...}
\end{figure*}
```

Control-variant figures (full width, one per float):

```latex
\begin{figure*}\centering\includegraphics[width=7.1in]{figA_ctrl_sersic.pdf}\caption{...}\end{figure*}
\begin{figure*}\centering\includegraphics[width=7.1in]{figA_ctrl_mag.pdf}\caption{...}\end{figure*}
\begin{figure*}\centering\includegraphics[width=7.1in]{figA_ctrl_color.pdf}\caption{...}\end{figure*}
```
