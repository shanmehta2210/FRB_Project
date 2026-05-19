"""
Generate CDF bias plots for sigma-weighted GALFIT results.
Produces plots in plots_nopsf_sigma/ and plots_psf_sigma/,
matching the format of the existing plots_psf/ and plots_nopsf/ directories.
"""

import os
import csv
import math
import random
import re
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from scipy import interpolate
from astropy.io import ascii


def parse_float(value):
    if value is None:
        return float('nan')
    if isinstance(value, (int, float, np.floating, np.integer)):
        return float(value)
    text = str(value).strip()
    if text == '':
        return float('nan')
    m = re.search(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', text)
    return float(m.group(0)) if m else float('nan')


def build_mc_cosi_draws(frb_list, ba_map, ba_err_map, n_draws=500):
    draws = []
    for _ in range(n_draws):
        vals = []
        for frb in frb_list:
            mu = parse_float(ba_map.get(frb))
            sig = parse_float(ba_err_map.get(frb))
            if not np.isfinite(mu):
                continue
            if not np.isfinite(sig) or sig < 0:
                sig = 0.0
            sampled = np.random.normal(mu, sig) if sig > 0 else mu
            sampled = min(1.0, max(0.0, float(sampled)))
            vals.append(hubble_cosi(sampled))
        draws.append(sorted(vals))
    return draws


def hubble_inclination_angle(b_over_a, q0=0.2):
    """Return inclination angle in degrees from axis ratio."""
    if b_over_a <= q0:
        return 90.0
    val = (b_over_a**2 - q0**2) / (1 - q0**2)
    if val < 0:
        val = 0
    if val > 1:
        val = 1
    return math.degrees(math.acos(math.sqrt(val)))


def hubble_cosi(b_over_a, q0=0.2):
    """Return cos(i) from axis ratio (for CDF)."""
    val = (b_over_a**2 - q0**2) / (1 - q0**2)
    if val < 0:
        return 0.0
    if val > 1:
        return 1.0
    return math.sqrt(val)


def generate_plots(
    frb_inc_angles,
    out_dir,
    label_suffix,
    b_over_a_sdss,
    frb_sample,
    frb_ba_map=None,
    frb_ba_err_map=None,
    mc_draws=500,
    mc_alpha=0.03,
):
    """Generate the full set of bias plots for a given set of FRB inclination angles."""
    os.makedirs(out_dir, exist_ok=True)

    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    font_prop = font_manager.FontProperties(family='Arial', style='normal', size=7)

    # Compute cos(i) for FRB sample
    frb_list = frb_sample['FRB'].tolist()
    incl_ang = np.array([frb_inc_angles[f] for f in frb_list])
    cosis = np.cos(np.radians(incl_ang))
    cosis_sorted = sorted([c for c in cosis if not np.isnan(c)])

    # --- SDSS bin sampling ---
    bins_sdss = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    bins2 = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    bins_collected = [[] for _ in range(10)]
    means = []
    stdevs = [[], []]

    print(f"  Generating SDSS samples for {label_suffix}...")
    for _ in range(10000):
        new_bins = np.histogram(random.sample(b_over_a_sdss, len(cosis_sorted)), bins=bins_sdss)
        for j in range(10):
            bins_collected[j].append(new_bins[0][j] / len(cosis_sorted))

    for i in range(10):
        means.append(np.mean(bins_collected[i]))
        bins_collected[i].sort()
        stdevs[0].append(abs(means[i] - bins_collected[i][1600]))
        stdevs[1].append(abs(bins_collected[i][8400] - means[i]))

    # --- CDF sampling ---
    sdss_cdf_list = sorted(b_over_a_sdss + [1])
    indices_norm = [0] + [i / len(b_over_a_sdss) for i in range(len(b_over_a_sdss))]
    sdss_cdf_list_full = [0] + sdss_cdf_list

    total_samples = []
    print(f"  Generating CDFs for {label_suffix}...")
    for _ in range(10000):
        new_sample = sorted(random.sample(b_over_a_sdss, len(cosis_sorted)))
        idx_norm = [0] + [j / len(cosis_sorted) for j in range(1, len(cosis_sorted) + 1)] + [1]
        new_sample = [0] + new_sample + [1]
        total_samples.append(interpolate.interp1d(new_sample, idx_norm))

    means_cdf = []
    sigma_down = []
    sigma_up = []
    for value in np.linspace(0, 1, 100):
        index_sample = sorted([s(value) for s in total_samples])
        means_cdf.append(np.mean(index_sample))
        sigma_down.append(index_sample[1600])
        sigma_up.append(index_sample[8400])

    # --- Plot 1: SDSS Distribution ---
    print(f"  Plotting sdss_dist_frb for {label_suffix}...")
    plt.figure(figsize=(10, 6))
    plt.bar(bins2, means, yerr=stdevs, align='edge', width=0.1, zorder=3, alpha=0.5,
            edgecolor='black', capsize=10, label="SDSS Distribution ($m_r < 21$, q0 = 0.2)")
    plt.plot((0, 1), (0.1, 0.1), color='black', label="Uniform Distribution")
    ax = plt.gca()
    ax.axvspan(min(cosis_sorted), max(cosis_sorted), color='gray', alpha=0.3,
               label="FRB host galaxies sample range", zorder=1)
    plt.ylabel('Probability', fontproperties=font_prop, fontsize=8)
    plt.xlabel('cos(i)', fontproperties=font_prop, fontsize=8)
    plt.ylim(0, 0.4)
    plt.xlim(0, 1)
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(np.cos(np.radians([90, 78, 66, 53, 37, 0])))
    ax2.set_xticklabels(['90°', '78°', '66°', '53°', '37°', '0°'], fontproperties=font_prop)
    ax2.set_xlabel('Inclination angle i (degrees)', fontproperties=font_prop, fontsize=8)
    ax.tick_params(axis='y', labelsize=8)
    ax.tick_params(axis='x', labelsize=8)
    ax2.tick_params(axis='x', labelsize=8)
    for spine in ax.spines.values(): spine.set_linewidth(1)
    for spine in ax2.spines.values(): spine.set_linewidth(1)
    ax.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5, zorder=0)
    plt.legend(fontsize=7)
    plt.savefig(os.path.join(out_dir, "sdss_dist_frb.pdf"), bbox_inches="tight", dpi=300)
    plt.close()

    # --- Plot 2: CDF Bias ---
    print(f"  Plotting CDF_bias for {label_suffix}...")
    plt.figure()
    y = [0] + [i / len(cosis_sorted) for i in range(1, len(cosis_sorted) + 1)] + [1]
    if frb_ba_map is not None and frb_ba_err_map is not None:
        mc_samples = build_mc_cosi_draws(frb_list, frb_ba_map, frb_ba_err_map, n_draws=mc_draws)
        funcs = []
        for draw in mc_samples:
            x_draw = [0.0] + draw + [1.0]
            y_draw = [0.0] + [j / len(draw) for j in range(1, len(draw) + 1)] + [1.0]
            plt.step(x_draw, y_draw, where='mid', color='red', linewidth=0.9, alpha=mc_alpha)
            idx_norm = [0.0] + [j / len(draw) for j in range(1, len(draw) + 1)] + [1.0]
            funcs.append(interpolate.interp1d(x_draw, idx_norm))
        x_grid = np.linspace(0, 1, 100)
        y_mean = [np.mean([float(f(v)) for f in funcs]) for v in x_grid]
        plt.plot(x_grid, y_mean, color='red', linewidth=2.0, label="FRB host galaxies (MC, errors)")
    else:
        x = [0] + cosis_sorted + [1]
        plt.step(x, y, where='mid', color='red', label="FRB host galaxies")
    plt.plot(np.linspace(0, 1, 100), means_cdf, color='black', label='SDSS distribution (mr < 21)')
    plt.fill_between(np.linspace(0, 1, 100), sigma_down, sigma_up, color='gray', alpha=0.4,
                     label='68% Confidence Interval')
    plt.plot((0, 1), (0, 1), color='black', linestyle='--', label="Uniform Distribution")
    plt.ylabel("Cumulative distribution", fontproperties=font_prop)
    plt.xlabel("cos(i)", fontproperties=font_prop)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    ax = plt.gca()
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(np.cos(np.radians([90, 78, 66, 53, 37, 0])))
    ax2.set_xticklabels(['90°', '78°', '66°', '53°', '37°', '0°'], fontproperties=font_prop)
    ax2.set_xlabel('Inclination angle i (degrees)', fontproperties=font_prop)
    plt.legend(fontsize=7)
    plt.savefig(os.path.join(out_dir, "CDF_bias.pdf"), bbox_inches="tight", dpi=300)
    plt.close()

    # --- Plot 3: CDF number of hosts ---
    print(f"  Plotting cdf_number_of_host for {label_suffix}...")
    plt.figure(figsize=(8, 8))
    x = [0] + cosis_sorted + [1]
    plt.step(x, y, where='mid', color='red', label="FRB host sample")
    plt.ylabel("Cumulative distribution")
    plt.xlabel("cos(i)")
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    ax = plt.gca()
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(np.cos(np.radians([90, 78, 66, 53, 37, 0])))
    ax2.set_xticklabels(['90°', '78°', '66°', '53°', '37°', '0°'])
    ax2.set_xlabel('Inclination angle i (degrees)')
    plt.legend(fontsize=12)
    plt.savefig(os.path.join(out_dir, "cdf_number_of_host.pdf"), bbox_inches="tight", dpi=300)
    plt.close()

    # --- Plot 4: Figure 2 (bin barchart + boxplots) ---
    print(f"  Plotting Figure_2 for {label_suffix}...")
    b_val = 59.5; c_val = 44.7; d_val = 28.3; a_val = 90; e_val = 0
    bin_edges = [[e_val, d_val], [d_val, c_val], [c_val, b_val], [b_val, a_val]]

    def values_for_bins(data):
        bins_data = [[], [], [], []]
        for i, value in enumerate(data):
            for idx, (lo, hi) in enumerate(bin_edges):
                if lo <= incl_ang[i] < hi:
                    bins_data[idx].append(value)
                elif hi == 90 and incl_ang[i] >= 90:
                    bins_data[idx].append(value)
        return bins_data

    DMh = frb_sample['DM_Host'].values
    z = frb_sample['z'].values
    SFR = frb_sample['SFR'].values
    sSFR = frb_sample['log(sSFR)'].values

    sfr_bins = values_for_bins(SFR)
    ssfr_bins = values_for_bins(sSFR)
    z_bins = values_for_bins(z)
    dmh_bins = values_for_bins(DMh)
    for i in range(4):
        if not sfr_bins[i]: sfr_bins[i] = [np.nan]
        if not ssfr_bins[i]: ssfr_bins[i] = [np.nan]
        if not z_bins[i]: z_bins[i] = [np.nan]
        if not dmh_bins[i]: dmh_bins[i] = [np.nan]

    bin_labels = ['Bin 1', 'Bin 2', 'Bin 3', 'Bin 4']
    bin_counts = [len([v for v in b_ if not np.isnan(v)]) for b_ in values_for_bins(incl_ang)]

    plt.figure(figsize=(5, 6))

    ax_a = plt.subplot2grid((3, 4), (0, 1), colspan=2)
    ax_a.bar(bin_labels, bin_counts, zorder=3, color='orange')
    ax_a.set_ylabel('Number of FRB host galaxies', fontproperties=font_prop)
    ax_a.text(0.95, 0.95, "(a)", fontsize=8, weight='bold', ha='right', va='top', transform=ax_a.transAxes)
    ax_a.tick_params(axis='both', which='major', labelsize=7)
    for spine in ax_a.spines.values(): spine.set_linewidth(1)
    ax2_a = ax_a.twiny()
    ax2_a.set_xlim(ax_a.get_xlim())
    ax2_a.set_xticks([i for i in range(4)])
    ax2_a.set_xticklabels([f'{a}°' for a in [21, 51, 68, 83]], fontproperties=font_prop)
    ax2_a.set_xlabel('Inclination angle (degrees)', fontproperties=font_prop, fontsize=7)
    ax2_a.tick_params(axis='x', labelsize=7)

    def make_boxplot(ax, data, ylabel, panel_label):
        valid = [[x for x in b if not np.isnan(x)] for b in data]
        if any(valid):
            ax.boxplot(valid, boxprops=dict(linewidth=0.5), whiskerprops=dict(linewidth=0.5),
                       capprops=dict(linewidth=0.5), medianprops=dict(linewidth=0.5),
                       flierprops=dict(marker='o', markersize=2, linestyle='none', linewidth=0.25))
        ax.set_xticks(ticks=[1, 2, 3, 4])
        ax.set_xticklabels(bin_labels, fontproperties=font_prop)
        ax.set_ylabel(ylabel, fontproperties=font_prop)
        ax.text(0.95, 0.95, panel_label, fontsize=8, weight='bold', ha='right', va='top', transform=ax.transAxes)
        ax.tick_params(axis='both', labelsize=7)

    make_boxplot(plt.subplot2grid((3, 4), (1, 0), colspan=2), dmh_bins,
                 r'DM$_{\mathdefault{host}}$ (pc cm$^{-3}$)', "(b)")
    make_boxplot(plt.subplot2grid((3, 4), (1, 2), colspan=2), sfr_bins,
                 'SFR (M$\\odot$/yr)', "(c)")
    make_boxplot(plt.subplot2grid((3, 4), (2, 0), colspan=2), ssfr_bins,
                 'Log(sSFR)', "(d)")
    make_boxplot(plt.subplot2grid((3, 4), (2, 2), colspan=2), z_bins,
                 'Redshift', "(e)")

    plt.subplots_adjust(wspace=0.5, hspace=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "Figure_2.pdf"), format='pdf', bbox_inches="tight", dpi=300)
    plt.close()

    # --- Plot 5: Bin barchart standalone ---
    print(f"  Plotting bin_barchart for {label_suffix}...")
    plt.figure(figsize=(6, 4))
    plt.bar(bin_labels, bin_counts, color='orange', edgecolor='black', zorder=3)
    plt.ylabel('Number of FRB host galaxies', fontproperties=font_prop, fontsize=9)
    plt.xlabel('Inclination Bin', fontproperties=font_prop, fontsize=9)
    for i, v in enumerate(bin_counts):
        plt.text(i, v + 0.2, str(v), ha='center', fontsize=8)
    plt.grid(axis='y', alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "bin_barchart.pdf"), bbox_inches="tight", dpi=300)
    plt.close()

    print(f"  Done! All plots saved to {out_dir}/")


def main():
    random.seed(42)

    # Load sigma GALFIT results
    sigma_data = {}
    with open('galfit_sigma_metrics_summary.csv') as f:
        for row in csv.DictReader(f):
            sigma_data[row['FRB']] = row

    # Compute inclination angles for nopsf and psf sigma cases
    inc_nopsf = {}
    inc_psf = {}
    ba_nopsf_map = {}
    ba_psf_map = {}
    ba_err_nopsf_map = {}
    ba_err_psf_map = {}
    for frb, row in sigma_data.items():
        ba_nopsf = parse_float(str(row.get('b_a_nopsf', '')).replace('*', ''))
        ba_psf = parse_float(str(row.get('b_a_psf', '')).replace('*', ''))
        ba_err_nopsf = parse_float(row.get('b_a_err_nopsf'))
        ba_err_psf = parse_float(row.get('b_a_err_psf'))
        ba_nopsf_map[frb] = ba_nopsf
        ba_psf_map[frb] = ba_psf
        ba_err_nopsf_map[frb] = ba_err_nopsf
        ba_err_psf_map[frb] = ba_err_psf
        try:
            inc_nopsf[frb] = hubble_inclination_angle(float(ba_nopsf))
        except (ValueError, TypeError):
            inc_nopsf[frb] = float('nan')
        try:
            inc_psf[frb] = hubble_inclination_angle(float(ba_psf))
        except (ValueError, TypeError):
            inc_psf[frb] = float('nan')

    # Load SDSS background
    print("Loading SDSS catalog...")
    sdss_data = ascii.read("SDSS_catalog.csv")
    b_over_a_sdss = []
    for i in range(len(sdss_data)):
        try:
            if sdss_data['petroMag_r'][i] <= 21 and sdss_data['expAB_r'][i] > 0.2:
                b_over_a_sdss.append(hubble_cosi(sdss_data['expAB_r'][i]))
        except:
            pass
    print(f"  SDSS sample: {len(b_over_a_sdss)} galaxies")

    # Load FRB sample
    frb_sample = pd.read_csv('frb_sample.txt', sep='\t')

    # Generate plots
    print("\n=== No-PSF + Sigma ===")
    generate_plots(
        inc_nopsf,
        "plots_nopsf_sigma",
        "nopsf_sigma",
        b_over_a_sdss,
        frb_sample,
        frb_ba_map=ba_nopsf_map,
        frb_ba_err_map=ba_err_nopsf_map,
    )

    print("\n=== With-PSF + Sigma ===")
    generate_plots(
        inc_psf,
        "plots_psf_sigma",
        "psf_sigma",
        b_over_a_sdss,
        frb_sample,
        frb_ba_map=ba_psf_map,
        frb_ba_err_map=ba_err_psf_map,
    )

    print("\nAll done!")


if __name__ == '__main__':
    main()
