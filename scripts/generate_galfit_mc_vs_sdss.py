import argparse
import math
import os
import sys
from pathlib import Path

import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import cosi_array_from_df, prepare_null_sample

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SDSS = REPO_ROOT / "catalog/SDSS_catalog_v1_allsky_modelmr.csv"
DEFAULT_OUT_SUBDIR = "plots/plots_null/v1_null_plots"

np.random.seed(42)

def cdf_envelope_null(reference_vals: np.ndarray, n_sample: int, n_draws: int = 10000):
    rng = np.random.default_rng(42)
    ref = np.asarray(reference_vals, dtype=float)
    ref = ref[np.isfinite(ref)]

    if len(ref) < n_sample:
        raise RuntimeError(f"Reference sample too small: {len(ref)} < {n_sample}")

    total_samples = []
    for _ in range(n_draws):
        draw = np.sort(rng.choice(ref, size=n_sample, replace=False))
        idx_norm = [0.0] + [j / n_sample for j in range(1, n_sample + 1)] + [1.0]
        draw_ext = [0.0] + draw.tolist() + [1.0]
        total_samples.append(interpolate.interp1d(draw_ext, idx_norm))

    x = np.linspace(0, 1, 500)
    means, down, up = [], [], []
    lo = int(0.16 * n_draws)
    hi = int(0.84 * n_draws)
    for value in x:
        idx_sample = sorted(float(s(value)) for s in total_samples)
        means.append(float(np.mean(idx_sample)))
        down.append(float(idx_sample[lo]))
        up.append(float(idx_sample[hi]))

    return x, np.array(means), np.array(down), np.array(up)

def cdf_envelope_mc_measured(inc_deg: np.ndarray, inc_err_deg: np.ndarray, n_draws: int = 10000):
    rng = np.random.default_rng(42)
    # Cos(i) mapping logic
    x = np.linspace(0, 1, 500)
    
    n_sample = len(inc_deg)
    total_samples = []
    
    for _ in range(n_draws):
        # Sample inclination from gaussian using original errors
        mc_inc = rng.normal(inc_deg, inc_err_deg)
        mc_inc = np.clip(mc_inc, 0, 90) # physical boundaries
        
        # Convert to cosi
        mc_cosi = np.cos(np.radians(mc_inc))
        
        # Calculate empirical CDF for this specific draw
        draw = np.sort(mc_cosi) # sort cosi
        idx_norm = [0.0] + [j / n_sample for j in range(1, n_sample + 1)] + [1.0]
        draw_ext = [0.0] + draw.tolist() + [1.0]
        
        # Create an interpolator. To avoid duplicated points which break interpolation:
        # We ensure monotonicity. Because of clipping, we might have multiple exactly 0s or 1s.
        # So we just use np.searchsorted later to evaluate the CDF directly rather than interp1d.
        total_samples.append(draw)

    # Evaluate the exact empirical CDF on the fixed grid x for each draw
    # For a sorted array `draw`, searchsorted gives the count of elements <= value
    cdf_evals = np.zeros((n_draws, len(x)))
    for i, draw in enumerate(total_samples):
        # fraction of samples <= x_val
        y_vals = np.searchsorted(draw, x, side='right') / n_sample
        cdf_evals[i, :] = y_vals
        
    means = np.mean(cdf_evals, axis=0)
    # For 68% interval (roughly 1 sigma), get 16th and 84th percentiles
    lo_bound = np.percentile(cdf_evals, 16, axis=0)
    hi_bound = np.percentile(cdf_evals, 84, axis=0)
    
    return x, means, lo_bound, hi_bound

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdss-csv", default=str(DEFAULT_SDSS))
    parser.add_argument("--inc-csv", default="legacy_vs_galfit_two_inclinations.csv")
    parser.add_argument("--mag-limit", type=float, default=21.0)
    parser.add_argument("--mag-column", default="rmag")
    parser.add_argument(
        "--sample-mode",
        choices=("strict",),
        default="strict",
    )
    parser.add_argument("--q0", type=float, default=0.2)
    parser.add_argument("--mc-draws", type=int, default=10000)
    parser.add_argument("--out-dir", default=DEFAULT_OUT_SUBDIR)
    parser.add_argument("--tag", default="v1_allsky_modelmr")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    sdss = pd.read_csv(args.sdss_csv)
    inc_df = pd.read_csv(args.inc_csv)

    galfit_inc = pd.to_numeric(inc_df['galfit_inc_deg'], errors='coerce').values
    galfit_err = pd.to_numeric(inc_df['galfit_inc_err_deg'], errors='coerce').values
    
    valid_galfit = ~np.isnan(galfit_inc) & ~np.isnan(galfit_err)
    galfit_inc = galfit_inc[valid_galfit]
    galfit_err = galfit_err[valid_galfit]
    n_frb = len(galfit_inc)

    print(f"Using {n_frb} valid GALFIT measurements for Monte Carlo...")

    sdss_cut = prepare_null_sample(
        sdss,
        sample_mode=args.sample_mode,
        mag_column=args.mag_column,
        mag_limit=args.mag_limit,
        q0=args.q0,
        is_legacy=False,
    )
    sdss_cosi = cosi_array_from_df(sdss_cut, q0=args.q0)

    # 1. SDSS Null Envelope
    x_s, mean_s, lo_s, hi_s = cdf_envelope_null(sdss_cosi, n_sample=n_frb, n_draws=int(args.mc_draws/2)) # Null uses resampling
    
    # 2. GALFIT Error Monte Carlo Envelope
    x_mc, mean_mc, lo_mc, hi_mc = cdf_envelope_mc_measured(galfit_inc, galfit_err, n_draws=args.mc_draws)

    # 3. Original exactly measured Galfit CDF
    orig_galfit_cosi = np.cos(np.radians(galfit_inc))
    sorted_orig_g = np.sort(orig_galfit_cosi)
    y_orig_g = np.arange(1, n_frb + 1) / n_frb

    font_prop = font_manager.FontProperties(family="Arial", style="normal", size=8)
    fig, ax = plt.subplots(figsize=(8, 8))

    # SDSS Plot
    ax.plot(
        x_s,
        mean_s,
        color="#377eb8",
        linewidth=2.0,
        label=f"SDSS model $m_r$ null ({args.sample_mode})",
    )
    ax.fill_between(x_s, lo_s, hi_s, color="#377eb8", alpha=0.20, label="SDSS 68% CI")

    # GALFIT MC Error envelope
    ax.plot(x_mc, mean_mc, color="#e41a1c", linewidth=2.0, linestyle='--', label=f"GALFIT Monte Carlo Mean")
    ax.fill_between(x_mc, lo_mc, hi_mc, color="#e41a1c", alpha=0.25, label="GALFIT 68% Uncertainty Band")
    
    # GALFIT original line
    # Prepend 0,0 to make step look right from origin
    sg = np.concatenate([[0.0], sorted_orig_g, [1.0]])
    yg = np.concatenate([[0.0], y_orig_g, [1.0]])
    ax.step(sg, yg, where='post', color='darkred', lw=2.0, label='GALFIT Point Estimates')

    ax.plot((0, 1), (0, 1), color="black", linestyle=":", linewidth=1.0, label="Uniform distribution")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("cos(i)", fontproperties=font_prop, fontsize=11)
    ax.set_ylabel("Cumulative distribution", fontproperties=font_prop, fontsize=11)
    ax.set_title(
        f"GALFIT vs SDSS model $m_r$ null ({args.sample_mode}, joint footprint)",
        fontsize=12,
    )
    ax.legend(fontsize=9, loc="upper left")

    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    tick_vals = np.cos(np.radians([90, 78, 66, 53, 37, 0]))
    ax_top.set_xticks(tick_vals)
    ax_top.set_xticklabels(["90", "78", "66", "53", "37", "0"], fontproperties=font_prop)
    ax_top.set_xlabel("Inclination angle i (degrees)", fontproperties=font_prop, fontsize=10)

    plt.tight_layout()
    out_png = os.path.join(args.out_dir, f"{args.tag}.png")
    out_pdf = os.path.join(args.out_dir, f"{args.tag}.pdf")
    
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, bbox_inches="tight")
    plt.close()
    print(f"Generated GALFIT MC vs SDSS plots at:\n{out_png}\n{out_pdf}")


if __name__ == '__main__':
    main()