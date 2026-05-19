import os
import argparse
import pandas as pd
import numpy as np
import random
import re
import matplotlib.pyplot as plt
from scipy import interpolate
import matplotlib.font_manager as font_manager
from astropy.io import ascii
import matplotlib
import sys


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


def build_mc_cosi_draws(frb_list, ba_map, ba_err_map, hubble_func, n_draws=500):
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
            vals.append(hubble_func(sampled))
        draws.append(sorted(vals))
    return draws

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--type', choices=['nopsf', 'psf'], required=True)
    parser.add_argument('--error-source', default='galfit_sigma_metrics_summary.csv')
    parser.add_argument('--mc-draws', type=int, default=500)
    parser.add_argument('--mc-alpha', type=float, default=0.03)
    args = parser.parse_args()

    out_dir = f"plots_{args.type}"
    os.makedirs(out_dir, exist_ok=True)
    # Removed custom style usage

    matplotlib.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    font_prop = font_manager.FontProperties(family='Arial', style='normal', size=7)

    def hubble_inclination(b_over_a):
        try:
            val = (b_over_a**2 - 0.2**2) / (1 - 0.2**2)
            if val < 0: val = 0
            return np.sqrt(val)
        except:
            return np.nan

    # Load SDSS
    print("Loading SDSS...")
    b_over_a_sdss = []
    data = ascii.read("SDSS_catalog.csv")
    for i in range(0, len(data)):
        try:
            if data['petroMag_r'][i] <= 21 and data['expAB_r'][i] > 0.2:
                b_over_a_sdss.append(hubble_inclination(data['expAB_r'][i]))
        except:
            pass

    # Read FRB sample
    print("Loading FRB sample...")
    df_frb = pd.read_csv('frb_sample.txt', sep='\t')
    frb_list = df_frb['FRB'].tolist()

    df_inc = pd.read_csv('master_frb_summary.csv')
    inc_map = dict(zip(df_inc['frb_name'], df_inc[f'inc_{args.type}']))

    ba_map = {}
    ba_err_map = {}
    if os.path.exists(args.error_source):
        err_df = pd.read_csv(args.error_source)
        ba_col = f'b_a_{args.type}'
        err_col = f'b_a_err_{args.type}'
        if ba_col in err_df.columns and err_col in err_df.columns:
            for _, row in err_df.iterrows():
                frb = str(row.get('FRB', ''))
                ba_map[frb] = parse_float(str(row.get(ba_col, '')).replace('*', ''))
                ba_err_map[frb] = parse_float(row.get(err_col))

    incl_ang = np.array([inc_map[f] for f in frb_list])
    cosis = np.cos(np.radians(incl_ang))
    cosis_sorted = pd.Series(cosis).dropna().sort_values().tolist()

    bins_sdss = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
    bins2 = [0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
    bins_collected = [[] for _ in range(10)]
    means = []
    stdevs = [[], []]

    print("Generating SDSS samples...")
    for i in range(10000):
        new_bins = np.histogram(random.sample(b_over_a_sdss, len(cosis_sorted)), bins=bins_sdss)
        for j in range(10):
            bins_collected[j].append(new_bins[0][j] / len(cosis_sorted))

    for i in range(10):
        means.append(np.mean(bins_collected[i]))
        bins_collected[i].sort()
        stdevs[0].append(abs(means[i] - bins_collected[i][1600]))
        stdevs[1].append(abs(bins_collected[i][8400] - means[i]))

    b_over_a_sdss.append(1)
    b_over_a_sdss.sort()
    
    indices_normalized = [0]
    for i in range(0, len(b_over_a_sdss)-1):
        indices_normalized.append(i/len(b_over_a_sdss))
    
    b_over_a_sdss.insert(0,0)
    
    total_samples = []
    sdss_raw = []
    
    print("Generating CDFs...")
    for i in range(10000):
        new_sample = random.sample(b_over_a_sdss, len(cosis_sorted))
        indices_normalized_3 = []
        for j in range(1, len(cosis_sorted)+1):
            indices_normalized_3.append(j/len(cosis_sorted))
        new_sample.insert(0,0)
        new_sample.append(1)
        indices_normalized_3.append(1)
        indices_normalized_3.insert(0,0)
        new_sample.sort()
        sdss_raw.append(new_sample)
        total_samples.append(interpolate.interp1d(new_sample, indices_normalized_3))
        
    means_cdf = []
    sigma_down = []
    sigma_up = []
    for value in np.linspace(0, 1, 100):
        index_sample = []
        for sample in total_samples:
            index_sample.append(sample(value))
        means_cdf.append(np.mean(index_sample))
        index_sample.sort()
        sigma_down.append(index_sample[1600])
        sigma_up.append(index_sample[8400])

    print("Plotting Figure 1 (SDSS Distribution)...")
    plt.figure(figsize=(10,6))
    plt.bar(bins2, means, yerr=stdevs, align='edge', width=0.1, zorder=3, alpha=0.5, edgecolor='black', capsize=10, label="SDSS Distribution ($m_r < 21$, q0 = 0.2)")
    plt.plot((0,1), (0.1,0.1), color='black', label="Uniform Distribution")
    ax = plt.gca()
    p1 = ax.axvspan(min(cosis_sorted), max(cosis_sorted), color='gray', alpha=0.3, label="FRB host galaxies sample range", zorder=1)
    plt.ylabel('Probability', fontproperties=font_prop, fontsize=8)
    plt.xlabel('cos(i)', fontproperties=font_prop, fontsize=8)
    plt.ylim(0, 0.4)
    plt.xlim(0, 1)
    
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(np.cos(np.radians([90,78,66,53,37,0])))
    ax2.set_xticklabels(['90°','78°','66°','53°','37°','0°'], fontproperties=font_prop)
    ax2.set_xlabel('Inclination angle i (degrees)', fontproperties=font_prop, fontsize=8)
    
    ax.tick_params(axis='y', labelsize=8)
    ax.tick_params(axis='x', labelsize=8)
    ax2.tick_params(axis='x', labelsize=8)
    for spine in ax.spines.values(): spine.set_linewidth(1)
    for spine in ax2.spines.values(): spine.set_linewidth(1)
    ax.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5, zorder=0)
    plt.legend(prop=font_prop)
    plt.savefig(os.path.join(out_dir, "sdss_dist_frb.pdf"), bbox_inches="tight", dpi=300)
    plt.close()

    print("Plotting CDF Bias...")
    plt.figure()
    y = [0] + [i/len(cosis_sorted) for i in range(1, len(cosis_sorted)+1)] + [1]
    if ba_map and ba_err_map:
        mc_samples = build_mc_cosi_draws(
            frb_list,
            ba_map,
            ba_err_map,
            hubble_inclination,
            n_draws=args.mc_draws,
        )
        funcs = []
        for draw in mc_samples:
            x_draw = [0.0] + draw + [1.0]
            y_draw = [0.0] + [j / len(draw) for j in range(1, len(draw) + 1)] + [1.0]
            plt.step(x_draw, y_draw, where='mid', color='red', linewidth=0.9, alpha=args.mc_alpha)
            idx_norm = [0.0] + [j / len(draw) for j in range(1, len(draw) + 1)] + [1.0]
            funcs.append(interpolate.interp1d(x_draw, idx_norm))
        x_grid = np.linspace(0, 1, 100)
        y_mean = [np.mean([float(f(v)) for f in funcs]) for v in x_grid]
        plt.plot(x_grid, y_mean, color='red', linewidth=2.0, label="FRB host galaxies (MC, errors)")
    else:
        x = [0] + cosis_sorted + [1]
        plt.step(x, y, where='mid', color='red', label="FRB host galaxies")
    plt.plot(np.linspace(0,1,100), means_cdf, color='black', label='SDSS distribution (mr < 21)')
    plt.fill_between(np.linspace(0,1,100), sigma_down, sigma_up, color='gray', alpha=0.4, label='68% Confidence Interval')
    plt.plot((0,1), (0,1), color='black', linestyle='--', label="Uniform Distribution")
    plt.ylabel("Cumulative distribution", fontproperties=font_prop)
    plt.xlabel("cos(i)", fontproperties=font_prop)
    plt.xlim(0,1)
    plt.ylim(0,1)
    ax = plt.gca()
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(np.cos(np.radians([90,78,66,53,37,0])))
    ax2.set_xticklabels(['90°','78°','66°','53°','37°','0°'], fontproperties=font_prop)
    ax2.set_xlabel('Inclination angle i (degrees)', fontproperties=font_prop)
    plt.legend(prop=font_prop)
    plt.savefig(os.path.join(out_dir, "CDF_bias.pdf"), bbox_inches="tight", dpi=300)
    plt.close()
    
    print("Plotting CDF number of host...")
    plt.figure(figsize=(8,8))
    x = [0] + cosis_sorted + [1]
    plt.step(x, y, where='mid', color='red', label="FRB host sample")
    plt.ylabel("Cumulative distribution")
    plt.xlabel("cos(i)")
    plt.xlim(0,1)
    plt.ylim(0,1)
    ax = plt.gca()
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks(np.cos(np.radians([90,78,66,53,37,0])))
    ax2.set_xticklabels(['90°','78°','66°','53°','37°','0°'])
    ax2.set_xlabel('Inclination angle i (degrees)')
    plt.legend(prop={'size': 16})
    plt.savefig(os.path.join(out_dir, "cdf_number_of_host.pdf"), bbox_inches="tight", dpi=300)
    plt.close()

    print("Generating Figure 2 components...")
    b=59.5; c=44.7; d=28.3; a=90; e=0
    bin_ = [[e,d], [d,c], [c,b], [b,a]]
    
    def values_for_bins(data):
        bins_data = [[], [], [], []]
        for i, value in enumerate(data):
            for idx, (min_val, max_val) in enumerate(bin_):
                if min_val <= incl_ang[i] < max_val:
                    bins_data[idx].append(value)
                elif a == 90 and max_val == 90 and incl_ang[i] >= 90:
                    bins_data[idx].append(value)
        return bins_data

    DMh = df_frb['DM_Host'].values
    z = df_frb['z'].values
    SFR = df_frb['SFR'].values
    sSFR = df_frb['log(sSFR)'].values

    sfr_bins = values_for_bins(SFR)
    sffr_bins = values_for_bins(sSFR)
    z_bins = values_for_bins(z)
    dmh_bins = values_for_bins(DMh)
    
    for i in range(4):
        if not sfr_bins[i]: sfr_bins[i] = [np.nan]
        if not sffr_bins[i]: sffr_bins[i] = [np.nan]
        if not z_bins[i]: z_bins[i] = [np.nan]
        if not dmh_bins[i]: dmh_bins[i] = [np.nan]

    print("Plotting Figure 2...")
    plt.figure(figsize=(5, 6))
    ax_a = plt.subplot2grid((3, 4), (0, 1), colspan=2)
    bin_labels = ['Bin 1', 'Bin 2', 'Bin 3', 'Bin 4']
    bin_counts = [len([v for v in b_ if not np.isnan(v)]) for b_ in values_for_bins(incl_ang)]
    bar_chart = ax_a.bar(bin_labels, bin_counts, zorder=3, color='orange')
    ax_a.set_ylabel('Number of FRB host galaxies', fontproperties=font_prop)
    ax_a.text(0.95, 0.95, "(a)", fontsize=8, weight='bold', horizontalalignment='right', verticalalignment='top', transform=ax_a.transAxes)
    ax_a.tick_params(axis='both', which='major', labelsize=7)
    for spine in ax_a.spines.values(): spine.set_linewidth(1)
    
    ax2_a = ax_a.twiny()
    ax2_a.set_xlim(ax_a.get_xlim())
    mean_bin = [21, 51, 68, 83]
    ax2_a.set_xticks([i for i in range(len(mean_bin))])
    ax2_a.set_xticklabels([f'{angle}°' for angle in mean_bin], fontproperties=font_prop)
    ax2_a.set_xlabel('Inclination angle (degrees)', fontproperties=font_prop, fontsize=7)
    ax2_a.tick_params(axis='x', labelsize=7)

    ax_b = plt.subplot2grid((3, 4), (1, 0), colspan=2)
    valid_dmh = [[x for x in b if not np.isnan(x)] for b in dmh_bins]
    if any(valid_dmh):
        ax_b.boxplot(valid_dmh, boxprops=dict(linewidth=0.5), whiskerprops=dict(linewidth=0.5), capprops=dict(linewidth=0.5), medianprops=dict(linewidth=0.5), flierprops=dict(marker='o', markersize=2, linestyle='none', linewidth=0.25))
    ax_b.set_xticks(ticks=[1, 2, 3, 4])
    ax_b.set_xticklabels(bin_labels, fontproperties=font_prop)
    ax_b.set_ylabel(r'DM$_{\mathdefault{host}}$ (pc cm$^{-3}$)', fontproperties=font_prop)
    ax_b.text(0.95, 0.95, "(b)", fontsize=8, weight='bold', horizontalalignment='right', verticalalignment='top', transform=ax_b.transAxes)
    ax_b.tick_params(axis='both', labelsize=7)
    
    ax_c = plt.subplot2grid((3, 4), (1, 2), colspan=2)
    valid_sfr = [[x for x in b if not np.isnan(x)] for b in sfr_bins]
    if any(valid_sfr):
        ax_c.boxplot(valid_sfr, boxprops=dict(linewidth=0.5), whiskerprops=dict(linewidth=0.5), capprops=dict(linewidth=0.5), medianprops=dict(linewidth=0.5), flierprops=dict(marker='o', markersize=2, linestyle='none', linewidth=0.25))
    ax_c.set_xticks(ticks=[1, 2, 3, 4])
    ax_c.set_xticklabels(bin_labels, fontproperties=font_prop)
    ax_c.set_ylabel('SFR (M$\\odot$/yr)', fontproperties=font_prop)
    ax_c.text(0.95, 0.95, "(c)", fontsize=8, weight='bold', horizontalalignment='right', verticalalignment='top', transform=ax_c.transAxes)
    ax_c.tick_params(axis='both', labelsize=7)

    ax_d = plt.subplot2grid((3, 4), (2, 0), colspan=2)
    valid_ssfr = [[x for x in b if not np.isnan(x)] for b in sffr_bins]
    if any(valid_ssfr):
        ax_d.boxplot(valid_ssfr, boxprops=dict(linewidth=0.5), whiskerprops=dict(linewidth=0.5), capprops=dict(linewidth=0.5), medianprops=dict(linewidth=0.5), flierprops=dict(marker='o', markersize=2, linestyle='none', linewidth=0.25))
    ax_d.set_xticks(ticks=[1, 2, 3, 4])
    ax_d.set_xticklabels(bin_labels, fontproperties=font_prop)
    ax_d.set_ylabel('Log(sSFR)', fontproperties=font_prop)
    ax_d.text(0.95, 0.95, "(d)", fontsize=8, weight='bold', horizontalalignment='right', verticalalignment='top', transform=ax_d.transAxes)
    ax_d.tick_params(axis='both', labelsize=7)

    ax_e = plt.subplot2grid((3, 4), (2, 2), colspan=2)
    valid_z = [[x for x in b if not np.isnan(x)] for b in z_bins]
    if any(valid_z):
        ax_e.boxplot(valid_z, boxprops=dict(linewidth=0.5), whiskerprops=dict(linewidth=0.5), capprops=dict(linewidth=0.5), medianprops=dict(linewidth=0.5), flierprops=dict(marker='o', markersize=2, linestyle='none', linewidth=0.25))
    ax_e.set_xticks(ticks=[1, 2, 3, 4])
    ax_e.set_xticklabels(bin_labels, fontproperties=font_prop)
    ax_e.set_ylabel('Redshift', fontproperties=font_prop)
    ax_e.text(0.95, 0.95, "(e)", fontsize=8, weight='bold', horizontalalignment='right', verticalalignment='top', transform=ax_e.transAxes)
    ax_e.tick_params(axis='both', labelsize=7)

    plt.subplots_adjust(wspace=0.5, hspace=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "Figure_2.pdf"), format='pdf', bbox_inches="tight", dpi=300)
    plt.close()

    print("Success! Plots generated in", out_dir)

if __name__ == '__main__':
    main()
