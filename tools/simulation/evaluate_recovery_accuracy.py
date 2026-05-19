import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def hubble_inc_from_q(q, q0=0.2):
    if not np.isfinite(q) or q < 0: return np.nan
    q = np.clip(q, 0.0, 1.0)
    if q <= q0: return 90.0
    val = (q**2 - q0**2) / (1.0 - q0**2)
    val = np.clip(val, 0.0, 1.0)
    return float(np.degrees(np.arccos(np.sqrt(val))))

def evaluate():
    base = "tools/simulation"
    true_cat = pd.read_csv(f"{base}/mock_catalog.csv")
    galfit = pd.read_csv(f"{base}/mock_galfit_results.csv")
    phot = pd.read_csv(f"{base}/mock_photutils_results.csv")

    galfit['inc_galfit'] = galfit['b_a'].apply(hubble_inc_from_q)
    phot['inc_photutils'] = phot['b_over_a'].apply(hubble_inc_from_q)

    # Merge
    merged = pd.merge(true_cat, galfit[['FRB', 'inc_galfit', 're', 'status']], on='FRB', how='left')
    merged = pd.merge(merged, phot[['FRB', 'inc_photutils', 'sma_rep_pix', 'status']], on='FRB', suffixes=('_galfit', '_photutils'))
    
    # Exclude failed fits
    merged = merged[(merged['status_galfit'] == 'ok') & (merged['status_photutils'] == 'ok')]

    merged['err_galfit'] = merged['inc_galfit'] - merged['true_inc']
    merged['err_photutils'] = merged['inc_photutils'] - merged['true_inc']

    merged.to_csv(f"{base}/mock_recovery_summary.csv", index=False)

    # PLOT: Error in inclination vs True Re, segregated by Photutils and GALFIT
    res_list = sorted(merged['true_re_pix'].unique())
    n_res = len(res_list)
    fig, axes = plt.subplots(1, n_res, figsize=(n_res * 4.5, 4), sharey=True)
    if n_res == 1: axes = [axes]
    
    for ax, r in zip(axes, res_list):
        sub = merged[merged['true_re_pix'] == r]
        ax.scatter(sub['true_inc'], sub['err_galfit'], c='r', alpha=0.6, label='GALFIT Error')
        ax.scatter(sub['true_inc'], sub['err_photutils'], c='b', alpha=0.6, marker='x', label='Photutils Error')
        ax.axhline(0, color='k', linestyle='--')
        ax.set_title(f"True Re = {r} pix")
        ax.set_xlabel("True Inclination (deg)")
        if ax == axes[0]:
            ax.set_ylabel(r"$\Delta i$ (Recovered - True)")
            ax.legend()
        ax.grid(True, alpha=0.3)

    fig_path = f"{base}/inclination_recovery_vs_re.png"
    plt.tight_layout()
    plt.savefig(fig_path)
    
    # Summary texts
    print(f"Evaluated {len(merged)} successful mock recoveries.")
    for r in res_list:
        sub = merged[merged['true_re_pix'] == r]
        rmse_g = np.sqrt(np.mean(sub['err_galfit']**2))
        rmse_p = np.sqrt(np.mean(sub['err_photutils']**2))
        print(f"Re={r} | GALFIT RMSE: {rmse_g:.2f} deg | Photutils RMSE: {rmse_p:.2f} deg")

if __name__ == "__main__":
    evaluate()
