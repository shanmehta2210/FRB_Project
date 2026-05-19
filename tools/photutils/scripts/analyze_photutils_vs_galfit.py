import argparse
import math
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits


def hubble_inclination_from_q(q, q0=0.2):
    if not np.isfinite(q) or not np.isfinite(q0):
        return np.nan
    if q <= q0:
        return 90.0
    if q >= 1.0:
        return 0.0
    val = (q * q - q0 * q0) / (1.0 - q0 * q0)
    val = np.clip(val, 0.0, 1.0)
    return float(np.degrees(np.arccos(np.sqrt(val))))


def hubble_inclination_err_from_qerr(q, q_err, q0=0.2):
    if not np.isfinite(q) or not np.isfinite(q_err) or not np.isfinite(q0):
        return np.nan
    if q_err <= 0:
        return 0.0
    q_lo = max(0.0, q - q_err)
    q_hi = min(1.0, q + q_err)
    i_lo = hubble_inclination_from_q(q_lo, q0=q0)
    i_hi = hubble_inclination_from_q(q_hi, q0=q0)
    if not np.isfinite(i_lo) or not np.isfinite(i_hi):
        return np.nan
    return abs(i_hi - i_lo) / 2.0


def corr_safe(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(m) < 2:
        return np.nan
    return float(np.corrcoef(x[m], y[m])[0, 1])


def pixel_scale_arcsec_from_fits(fits_path):
    try:
        with fits.open(fits_path) as hdul:
            h = hdul[0].header
    except Exception:
        return np.nan

    # Prefer explicit CDELT keywords when present.
    cdelt_vals = [abs(float(h[k])) for k in ("CDELT1", "CDELT2") if k in h and np.isfinite(h[k])]
    if len(cdelt_vals) > 0:
        return float(np.nanmean(cdelt_vals) * 3600.0)

    # Fallback to CD matrix scale, robust to simple rotations.
    cd_keys = ["CD1_1", "CD1_2", "CD2_1", "CD2_2"]
    if all(k in h for k in cd_keys):
        cd11 = float(h["CD1_1"])
        cd12 = float(h["CD1_2"])
        cd21 = float(h["CD2_1"])
        cd22 = float(h["CD2_2"])
        sx = np.hypot(cd11, cd21) * 3600.0
        sy = np.hypot(cd12, cd22) * 3600.0
        vals = [v for v in (sx, sy) if np.isfinite(v) and v > 0]
        if len(vals) > 0:
            return float(np.nanmean(vals))

    return np.nan


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--photutils-csv", default="tools/photutils/results/photutils_ellipse_inclination_angles_nopsf.csv")
    parser.add_argument("--galfit-csv", default="galfit_sigma_metrics_summary.csv")
    parser.add_argument("--master-csv", default="master_frb_summary.csv")
    parser.add_argument("--out-root", default=".")
    parser.add_argument("--plots-dir", default="plots/plots_photutils")
    args = parser.parse_args()

    os.makedirs(args.plots_dir, exist_ok=True)

    phot = pd.read_csv(args.photutils_csv)
    gal = pd.read_csv(args.galfit_csv)
    master = pd.read_csv(args.master_csv)

    # Build root-level photutils master table for easy access.
    phot_ok = phot.copy()
    phot_ok["b_over_a"] = pd.to_numeric(phot_ok["b_over_a"], errors="coerce")
    phot_ok["b_over_a_err"] = pd.to_numeric(phot_ok["b_over_a_err"], errors="coerce")
    phot_ok["inclination_deg"] = pd.to_numeric(phot_ok["inclination_deg"], errors="coerce")
    phot_ok["inclination_err_deg"] = pd.to_numeric(phot_ok["inclination_err_deg"], errors="coerce")

    phot_master = phot_ok.rename(
        columns={
            "b_over_a": "b_a_photutils",
            "b_over_a_err": "b_a_err_photutils",
            "inclination_deg": "inc_photutils",
            "inclination_err_deg": "inc_err_photutils",
            "pa_deg": "pa_photutils",
            "pa_err_deg": "pa_err_photutils",
        }
    )

    phot_master_cols = [
        "FRB",
        "status",
        "q0",
        "b_a_photutils",
        "b_a_err_photutils",
        "inc_photutils",
        "inc_err_photutils",
        "pa_photutils",
        "pa_err_photutils",
        "n_isophotes",
        "sma_rep_pix",
        "half_light_flux_proxy",
        "fit_strategy_index",
        "bkg_median",
        "bkg_std",
        "source_flux",
        "source_sigma",
    ]
    phot_master = phot_master[phot_master_cols].sort_values("FRB").reset_index(drop=True)
    out_phot_master = os.path.join(args.out_root, "photutils_master_summary.csv")
    phot_master.to_csv(out_phot_master, index=False)

    # Prepare GALFIT inclination with per-FRB q0 from master when needed.
    q0_map = dict(zip(master["FRB"].astype(str), pd.to_numeric(master["q0"], errors="coerce")))

    gal["b_a_psf"] = pd.to_numeric(gal["b_a_psf"], errors="coerce")
    gal["b_a_err_psf"] = pd.to_numeric(gal["b_a_err_psf"], errors="coerce")
    gal["re_psf"] = pd.to_numeric(gal["re_psf"], errors="coerce")

    gal_inc = []
    gal_inc_err = []
    for _, row in gal.iterrows():
        frb = str(row["FRB"])
        q0 = q0_map.get(frb, np.nan)
        if not np.isfinite(q0):
            q0 = 0.2
        q = row["b_a_psf"]
        q_err = row["b_a_err_psf"]
        gal_inc.append(hubble_inclination_from_q(q, q0=q0))
        gal_inc_err.append(hubble_inclination_err_from_qerr(q, q_err, q0=q0))

    gal["inc_galfit"] = gal_inc
    gal["inc_err_galfit"] = gal_inc_err

    # Merge for comparison.
    cmp_df = pd.merge(
        phot_master,
        gal[["FRB", "b_a_psf", "b_a_err_psf", "inc_galfit", "inc_err_galfit", "chi2nu_psf", "re_psf"]],
        on="FRB",
        how="inner",
    )

    cmp_df["delta_inc_phot_minus_galfit"] = cmp_df["inc_photutils"] - cmp_df["inc_galfit"]
    cmp_df["delta_inc_galfit_minus_photutils"] = cmp_df["inc_galfit"] - cmp_df["inc_photutils"]
    cmp_df["delta_ba_phot_minus_galfit"] = cmp_df["b_a_photutils"] - cmp_df["b_a_psf"]

    cmp_df["pix_scale_arcsec"] = cmp_df["source_flux"].apply(
        lambda p: pixel_scale_arcsec_from_fits(os.path.normpath(p)) if isinstance(p, str) else np.nan
    )
    cmp_df["re_psf_arcsec"] = cmp_df["re_psf"] * cmp_df["pix_scale_arcsec"]
    cmp_df["sma_rep_arcsec"] = cmp_df["sma_rep_pix"] * cmp_df["pix_scale_arcsec"]

    out_cmp = os.path.join(args.out_root, "photutils_vs_galfit_comparison.csv")
    cmp_df.sort_values("FRB").to_csv(out_cmp, index=False)

    higher_mask = np.isfinite(cmp_df["delta_inc_phot_minus_galfit"]) & (cmp_df["delta_inc_phot_minus_galfit"] > 0)
    higher_df = cmp_df.loc[higher_mask, [
        "FRB",
        "inc_photutils",
        "inc_galfit",
        "delta_inc_phot_minus_galfit",
        "chi2nu_psf",
        "n_isophotes",
        "fit_strategy_index",
        "sma_rep_pix",
    ]].sort_values("delta_inc_phot_minus_galfit", ascending=False)
    out_higher = os.path.join(args.out_root, "photutils_higher_inclination_frbs.csv")
    higher_df.to_csv(out_higher, index=False)

    # Preliminary summary stats.
    valid_inc = np.isfinite(cmp_df["inc_photutils"]) & np.isfinite(cmp_df["inc_galfit"])
    valid_ba = np.isfinite(cmp_df["b_a_photutils"]) & np.isfinite(cmp_df["b_a_psf"])

    d_inc = cmp_df.loc[valid_inc, "delta_inc_phot_minus_galfit"].to_numpy(dtype=float)
    d_ba = cmp_df.loc[valid_ba, "delta_ba_phot_minus_galfit"].to_numpy(dtype=float)

    stats = {
        "n_total": int(len(cmp_df)),
        "n_valid_inc": int(np.count_nonzero(valid_inc)),
        "n_valid_ba": int(np.count_nonzero(valid_ba)),
        "mean_inc_photutils": float(np.nanmean(cmp_df["inc_photutils"])),
        "mean_inc_galfit": float(np.nanmean(cmp_df["inc_galfit"])),
        "mean_delta_inc_phot_minus_galfit": float(np.nanmean(d_inc)) if len(d_inc) else np.nan,
        "median_delta_inc_phot_minus_galfit": float(np.nanmedian(d_inc)) if len(d_inc) else np.nan,
        "std_delta_inc_phot_minus_galfit": float(np.nanstd(d_inc, ddof=1)) if len(d_inc) > 1 else np.nan,
        "rmse_delta_inc_phot_minus_galfit": float(np.sqrt(np.nanmean(d_inc**2))) if len(d_inc) else np.nan,
        "corr_inc_phot_vs_galfit": corr_safe(cmp_df["inc_photutils"], cmp_df["inc_galfit"]),
        "mean_ba_photutils": float(np.nanmean(cmp_df["b_a_photutils"])),
        "mean_ba_galfit": float(np.nanmean(cmp_df["b_a_psf"])),
        "mean_delta_ba_phot_minus_galfit": float(np.nanmean(d_ba)) if len(d_ba) else np.nan,
        "median_delta_ba_phot_minus_galfit": float(np.nanmedian(d_ba)) if len(d_ba) else np.nan,
        "std_delta_ba_phot_minus_galfit": float(np.nanstd(d_ba, ddof=1)) if len(d_ba) > 1 else np.nan,
        "rmse_delta_ba_phot_minus_galfit": float(np.sqrt(np.nanmean(d_ba**2))) if len(d_ba) else np.nan,
        "corr_ba_phot_vs_galfit": corr_safe(cmp_df["b_a_photutils"], cmp_df["b_a_psf"]),
    }

    stats_df = pd.DataFrame({"metric": list(stats.keys()), "value": list(stats.values())})
    out_stats = os.path.join(args.out_root, "photutils_vs_galfit_prelim_stats.csv")
    stats_df.to_csv(out_stats, index=False)

    # Quick diagnostic plots.
    # Inclination scatter
    plt.figure(figsize=(6, 6))
    plt.scatter(cmp_df["inc_galfit"], cmp_df["inc_photutils"], s=35)
    plt.plot([0, 90], [0, 90], "k--", linewidth=1)
    plt.xlim(0, 90)
    plt.ylim(0, 90)
    plt.xlabel("GALFIT inclination (deg)")
    plt.ylabel("Photutils inclination (deg)")
    plt.title("Photutils vs GALFIT inclinations")
    plt.tight_layout()
    out_sc_inc = os.path.join(args.plots_dir, "photutils_vs_galfit_inc_scatter.png")
    plt.savefig(out_sc_inc, dpi=220)
    plt.close()

    # b/a scatter
    plt.figure(figsize=(6, 6))
    plt.scatter(cmp_df["b_a_psf"], cmp_df["b_a_photutils"], s=35)
    plt.plot([0, 1], [0, 1], "k--", linewidth=1)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.xlabel("GALFIT b/a")
    plt.ylabel("Photutils b/a")
    plt.title("Photutils vs GALFIT axis ratio")
    plt.tight_layout()
    out_sc_ba = os.path.join(args.plots_dir, "photutils_vs_galfit_ba_scatter.png")
    plt.savefig(out_sc_ba, dpi=220)
    plt.close()

    # Delta histograms
    plt.figure(figsize=(7, 4.5))
    plt.hist(d_inc, bins=10, alpha=0.7)
    plt.axvline(np.nanmean(d_inc), color="r", linestyle="--", linewidth=1.2)
    plt.xlabel("Photutils - GALFIT inclination (deg)")
    plt.ylabel("Count")
    plt.title("Inclination difference distribution")
    plt.tight_layout()
    out_hist_inc = os.path.join(args.plots_dir, "photutils_minus_galfit_inc_hist.png")
    plt.savefig(out_hist_inc, dpi=220)
    plt.close()

    plt.figure(figsize=(7, 4.5))
    plt.hist(d_ba, bins=10, alpha=0.7)
    plt.axvline(np.nanmean(d_ba), color="r", linestyle="--", linewidth=1.2)
    plt.xlabel("Photutils - GALFIT b/a")
    plt.ylabel("Count")
    plt.title("Axis-ratio difference distribution")
    plt.tight_layout()
    out_hist_ba = os.path.join(args.plots_dir, "photutils_minus_galfit_ba_hist.png")
    plt.savefig(out_hist_ba, dpi=220)
    plt.close()

    # Requested: delta inclination (GALFIT - Photutils) vs GALFIT Reff in arcsec.
    m_re_galfit = np.isfinite(cmp_df["re_psf_arcsec"]) & np.isfinite(cmp_df["delta_inc_galfit_minus_photutils"])
    plt.figure(figsize=(7, 5))
    plt.scatter(cmp_df.loc[m_re_galfit, "re_psf_arcsec"], cmp_df.loc[m_re_galfit, "delta_inc_galfit_minus_photutils"], s=38)
    plt.axhline(0.0, color="k", linestyle="--", linewidth=1)
    plt.xlabel("GALFIT Reff (arcsec)")
    plt.ylabel("GALFIT - Photutils inclination (deg)")
    plt.title("Delta inclination vs GALFIT Reff (arcsec)")
    plt.tight_layout()
    out_sc_dinc_re_galfit = os.path.join(args.plots_dir, "delta_inc_galfit_minus_photutils_vs_galfit_reff_arcsec.png")
    plt.savefig(out_sc_dinc_re_galfit, dpi=220)
    plt.close()

    # Optional companion: use Photutils half-light isophote SMA proxy in arcsec.
    m_re_phot = np.isfinite(cmp_df["sma_rep_arcsec"]) & np.isfinite(cmp_df["delta_inc_galfit_minus_photutils"])
    plt.figure(figsize=(7, 5))
    plt.scatter(cmp_df.loc[m_re_phot, "sma_rep_arcsec"], cmp_df.loc[m_re_phot, "delta_inc_galfit_minus_photutils"], s=38)
    plt.axhline(0.0, color="k", linestyle="--", linewidth=1)
    plt.xlabel("Photutils half-light SMA proxy (arcsec)")
    plt.ylabel("GALFIT - Photutils inclination (deg)")
    plt.title("Delta inclination vs Photutils size proxy (arcsec)")
    plt.tight_layout()
    out_sc_dinc_re_phot = os.path.join(args.plots_dir, "delta_inc_galfit_minus_photutils_vs_photutils_sma_proxy_arcsec.png")
    plt.savefig(out_sc_dinc_re_phot, dpi=220)
    plt.close()

    # Human-readable markdown report.
    report_lines = [
        "# Preliminary Photutils vs GALFIT Comparison",
        "",
        f"Rows compared: {stats['n_total']}",
        f"Valid inclination pairs: {stats['n_valid_inc']}",
        f"Valid b/a pairs: {stats['n_valid_ba']}",
        "",
        "## Inclination summary",
        f"- Mean Photutils inclination: {stats['mean_inc_photutils']:.3f} deg",
        f"- Mean GALFIT inclination: {stats['mean_inc_galfit']:.3f} deg",
        f"- Mean delta (Photutils - GALFIT): {stats['mean_delta_inc_phot_minus_galfit']:.3f} deg",
        f"- Median delta: {stats['median_delta_inc_phot_minus_galfit']:.3f} deg",
        f"- Std delta: {stats['std_delta_inc_phot_minus_galfit']:.3f} deg",
        f"- RMSE delta: {stats['rmse_delta_inc_phot_minus_galfit']:.3f} deg",
        f"- Correlation (inc): {stats['corr_inc_phot_vs_galfit']:.3f}",
        "",
        "## Axis-ratio summary",
        f"- Mean Photutils b/a: {stats['mean_ba_photutils']:.4f}",
        f"- Mean GALFIT b/a: {stats['mean_ba_galfit']:.4f}",
        f"- Mean delta (Photutils - GALFIT): {stats['mean_delta_ba_phot_minus_galfit']:.4f}",
        f"- Median delta: {stats['median_delta_ba_phot_minus_galfit']:.4f}",
        f"- Std delta: {stats['std_delta_ba_phot_minus_galfit']:.4f}",
        f"- RMSE delta: {stats['rmse_delta_ba_phot_minus_galfit']:.4f}",
        f"- Correlation (b/a): {stats['corr_ba_phot_vs_galfit']:.3f}",
        "",
        "## Photutils Metrics For Fit Quality",
        "- GALFIT provides catalog chi2nu directly as chi2nu_psf.",
        "- This Photutils pipeline does not currently output a reduced-chi2 (chi2nu) value.",
        "- Available Photutils diagnostics in this run include n_isophotes, fit_strategy_index, b_a_err_photutils, inc_err_photutils, and bkg_std.",
        "- If needed, we can add a photutils chi2nu-like metric by rebuilding an ellipse model and computing residual/sigma over valid pixels.",
        "",
        "## FRBs Where Photutils Inclination > GALFIT",
        f"- Count: {int(len(higher_df))}",
        "- Saved table: photutils_higher_inclination_frbs.csv",
        "",
        "## Output files",
        "- photutils_master_summary.csv",
        "- photutils_vs_galfit_comparison.csv",
        "- photutils_vs_galfit_prelim_stats.csv",
        "- photutils_higher_inclination_frbs.csv",
        "- plots/plots_photutils/photutils_vs_galfit_inc_scatter.png",
        "- plots/plots_photutils/photutils_vs_galfit_ba_scatter.png",
        "- plots/plots_photutils/photutils_minus_galfit_inc_hist.png",
        "- plots/plots_photutils/photutils_minus_galfit_ba_hist.png",
        "- plots/plots_photutils/delta_inc_galfit_minus_photutils_vs_galfit_reff_arcsec.png",
        "- plots/plots_photutils/delta_inc_galfit_minus_photutils_vs_photutils_sma_proxy_arcsec.png",
    ]
    out_report = os.path.join(args.out_root, "Reports", "photutils_vs_galfit_prelim.md")
    os.makedirs(os.path.dirname(out_report), exist_ok=True)
    with open(out_report, "w", encoding="utf-8") as f:
        f.write("\n".join(report_lines))

    print(f"Saved: {out_phot_master}")
    print(f"Saved: {out_cmp}")
    print(f"Saved: {out_stats}")
    print(f"Saved: {out_higher}")
    print(f"Saved: {out_sc_inc}")
    print(f"Saved: {out_sc_ba}")
    print(f"Saved: {out_hist_inc}")
    print(f"Saved: {out_hist_ba}")
    print(f"Saved: {out_sc_dinc_re_galfit}")
    print(f"Saved: {out_sc_dinc_re_phot}")
    print(f"Saved: {out_report}")


if __name__ == "__main__":
    main()
