import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import matplotlib

matplotlib.use('Agg')

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--type', choices=['nopsf', 'psf'], required=True)
    args = parser.parse_args()

    out_dir = f"plots_{args.type}"
    os.makedirs(out_dir, exist_ok=True)

    font_prop = font_manager.FontProperties(family='Arial', style='normal', size=7)

    df_frb = pd.read_csv('frb_sample.txt', sep='\t')
    df_inc = pd.read_csv('master_frb_summary.csv')
    inc_map = dict(zip(df_inc['frb_name'], df_inc[f'inc_{args.type}']))
    frb_list = df_frb['FRB'].tolist()
    incl_ang = np.array([inc_map.get(f, np.nan) for f in frb_list])

    b = 59.5; c = 44.7; d = 28.3; a = 90; e = 0
    bin_ = [[e, d], [d, c], [c, b], [b, a]]

    def values_for_bins(data):
        bins_data = [[], [], [], []]
        for i, value in enumerate(data):
            for idx, (min_val, max_val) in enumerate(bin_):
                if min_val <= incl_ang[i] < max_val:
                    bins_data[idx].append(value)
                elif max_val == 90 and incl_ang[i] >= 90:
                    bins_data[idx].append(value)
        return bins_data

    bin_labels = [f'0°–28°', f'28°–45°', f'45°–60°', f'60°–90°']
    bin_counts = [len([v for v in b_ if not np.isnan(v)]) for b_ in values_for_bins(incl_ang)]

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(bin_labels, bin_counts, color='steelblue', edgecolor='black', linewidth=0.6)
    ax.set_ylabel('Number of FRB host galaxies', fontsize=9)
    ax.set_xlabel('Inclination angle bin', fontsize=9)
    ax.set_title('FRB Host Galaxies per Inclination Bin', fontsize=10)
    
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    ax.grid(axis='y', linestyle='--', linewidth=0.4, alpha=0.6)
    ax.tick_params(axis='both', labelsize=8)
    
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "bin_barchart.pdf"), bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Bar chart saved to {out_dir}/bin_barchart.pdf")

if __name__ == '__main__':
    main()
