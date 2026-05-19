# Preliminary Photutils vs GALFIT Comparison

Rows compared: 23
Valid inclination pairs: 23
Valid b/a pairs: 23

## Inclination summary
- Mean Photutils inclination: 43.895 deg
- Mean GALFIT inclination: 56.246 deg
- Mean delta (Photutils - GALFIT): -12.351 deg
- Median delta: -9.702 deg
- Std delta: 16.022 deg
- RMSE delta: 19.952 deg
- Correlation (inc): 0.502

## Axis-ratio summary
- Mean Photutils b/a: 0.7179
- Mean GALFIT b/a: 0.5678
- Mean delta (Photutils - GALFIT): 0.1501
- Median delta: 0.1347
- Std delta: 0.1917
- RMSE delta: 0.2402
- Correlation (b/a): 0.508

## Photutils Metrics For Fit Quality
- GALFIT provides catalog chi2nu directly as chi2nu_psf.
- This Photutils pipeline does not currently output a reduced-chi2 (chi2nu) value.
- Available Photutils diagnostics in this run include n_isophotes, fit_strategy_index, b_a_err_photutils, inc_err_photutils, and bkg_std.
- If needed, we can add a photutils chi2nu-like metric by rebuilding an ellipse model and computing residual/sigma over valid pixels.

## FRBs Where Photutils Inclination > GALFIT
- Count: 5
- Saved table: photutils_higher_inclination_frbs.csv

## Output files
- photutils_master_summary.csv
- photutils_vs_galfit_comparison.csv
- photutils_vs_galfit_prelim_stats.csv
- photutils_higher_inclination_frbs.csv
- plots/plots_photutils/photutils_vs_galfit_inc_scatter.png
- plots/plots_photutils/photutils_vs_galfit_ba_scatter.png
- plots/plots_photutils/photutils_minus_galfit_inc_hist.png
- plots/plots_photutils/photutils_minus_galfit_ba_hist.png
- plots/plots_photutils/delta_inc_galfit_minus_photutils_vs_galfit_reff_arcsec.png
- plots/plots_photutils/delta_inc_galfit_minus_photutils_vs_photutils_sma_proxy_arcsec.png