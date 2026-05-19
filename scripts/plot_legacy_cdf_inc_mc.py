import argparse
import math
import os
import random
import re
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

from null_catalog_utils import hubble_cosi_from_ba, prepare_null_sample

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_LEGACY = REPO_ROOT / "LS_catalog_v1_allsky_modelmr.csv"
DEFAULT_OUT_SUBDIR = "plots/plots_legacy_cdf/v1_null_plots"

random.seed(42)
np.random.seed(42)

def legacy_null_cdf_envelope_from_shape_mc(legacy_df, n_sample, n_draws=10000, q0=0.2):
    e1 = pd.to_numeric(legacy_df["shape_e1"], errors="coerce").to_numpy(dtype=float)
    e2 = pd.to_numeric(legacy_df["shape_e2"], errors="coerce").to_numpy(dtype=float)
    e1_ivar = pd.to_numeric(legacy_df["shape_e1_ivar"], errors="coerce").to_numpy(dtype=float)
    e2_ivar = pd.to_numeric(legacy_df["shape_e2_ivar"], errors="coerce").to_numpy(dtype=float)

    e1_sigma = np.zeros_like(e1, dtype=float)
    e2_sigma = np.zeros_like(e2, dtype=float)
    good_e1 = np.isfinite(e1_ivar) & (e1_ivar > 0)
    good_e2 = np.isfinite(e2_ivar) & (e2_ivar > 0)
    e1_sigma[good_e1] = 1.0 / np.sqrt(e1_ivar[good_e1])
    e2_sigma[good_e2] = 1.0 / np.sqrt(e2_ivar[good_e2])

    valid = np.isfinite(e1) & np.isfinite(e2)
    if int(np.count_nonzero(valid)) < n_sample:
        raise RuntimeError("Not enough valid Legacy shape entries for null sampling.")

    e1, e2 = e1[valid], e2[valid]
    e1_sigma, e2_sigma = e1_sigma[valid], e2_sigma[valid]

    rng = np.random.default_rng(42)
    total_samples = []

    for _ in range(n_draws):
        idx = rng.choice(len(e1), size=n_sample, replace=False)
        e1_draw = rng.normal(e1[idx], e1_sigma[idx])
        e2_draw = rng.normal(e2[idx], e2_sigma[idx])
        eabs = np.hypot(e1_draw, e2_draw)
        eabs = np.clip(eabs, 0.0, 0.999999)
        q = (1.0 - eabs) / (1.0 + eabs)

        cosi = [hubble_cosi_from_ba(float(qv), q0=q0) for qv in q]
        cosi_sorted = sorted(cosi)

        idx_norm = [0.0] + [j / n_sample for j in range(1, n_sample + 1)] + [1.0]
        draw_ext = [0.0] + cosi_sorted + [1.0]
        total_samples.append(interpolate.interp1d(draw_ext, idx_norm))

    x = np.linspace(0, 1, 100)
    means, down, up = [], [], []
    lo, hi = int(0.16 * n_draws), int(0.84 * n_draws)

    for value in x:
        idx_sample = sorted(float(s(value)) for s in total_samples)
        means.append(float(np.mean(idx_sample)))
        down.append(float(idx_sample[lo]))
        up.append(float(idx_sample[hi]))

    return x, np.array(means), np.array(down), np.array(up)

def mc_mean_cdf(cosi_draws):
    n_sample = len(cosi_draws[0])
    funcs = []
    for draw in cosi_draws:
        idx_norm = [0.0] + [j / n_sample for j in range(1, n_sample + 1)] + [1.0]
        draw_ext = [0.0] + draw + [1.0]
        funcs.append(interpolate.interp1d(draw_ext, idx_norm))

    x = np.linspace(0, 1, 100)
    mean = []
    for value in x:
        vals = [float(f(value)) for f in funcs]
        mean.append(float(np.mean(vals)))
    return x, np.array(mean)

def build_frb_mc_draws_inc(galfit_df, n_draws=500):
    draws = []
    inc_vals = []
    err_vals = []
    
    for _, row in galfit_df.iterrows():
        try:
            mu = float(row.get("inc_psf", np.nan))
            sigma = float(row.get("inc_err_psf", 0.0))
        except:
            continue
            
        if not np.isfinite(mu) or pd.isna(mu):
            continue
        if not np.isfinite(sigma) or pd.isna(sigma) or sigma < 0.0:
            sigma = 0.0
            
        mu = min(90.0, max(0.0, mu))
        inc_vals.append((mu, sigma))

    if len(inc_vals) == 0:
        raise RuntimeError("No valid FRBs found for MC sampling.")

    for _ in range(n_draws):
        cosi = []
        for mu, sigma in inc_vals:
            sampled_inc = np.random.normal(mu, sigma) if sigma > 0 else mu
            sampled_inc = min(90.0, max(0.0, float(sampled_inc)))
            # angle is in degrees
            cosi.append(math.cos(math.radians(sampled_inc)))
        draws.append(sorted(cosi))

    return draws

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--legacy-csv", default=str(DEFAULT_LEGACY))
    parser.add_argument("--galfit-csv", default="new_16_frbs_galfit_results.csv")
    parser.add_argument("--mag-limit", type=float, default=21.0)
    parser.add_argument("--mag-column", default="rmag")
    parser.add_argument(
        "--sample-mode",
        choices=("strict", "inclusive"),
        default="strict",
    )
    parser.add_argument("--q0", type=float, default=0.2)
    parser.add_argument("--mc-draws-frb", type=int, default=500)
    parser.add_argument("--mc-draws-null", type=int, default=10000)
    parser.add_argument("--mc-alpha", type=float, default=0.03)
    parser.add_argument("--out-dir", default=DEFAULT_OUT_SUBDIR)
    parser.add_argument("--tag", default="v1_allsky_modelmr")
    parser.add_argument("--exclude-types", default="REX", help="Comma-separated types to exclude")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    legacy = pd.read_csv(args.legacy_csv)
    galfit = pd.read_csv(args.galfit_csv)

    # 1. Filter out missing or bad FRBs
    galfit = galfit[galfit['inc_psf'].notna()].copy()
    
    # Exclude the 4 currently problematic FRBs
    exclude_frbs = ["20190611B", "20190711A", "20230526A", "20240310A"]
    galfit = galfit[~galfit['FRB'].isin(exclude_frbs)].copy()

    frb_cosi_draws = build_frb_mc_draws_inc(galfit, n_draws=args.mc_draws_frb)
    n_frb = len(frb_cosi_draws[0])

    legacy_cut = prepare_null_sample(
        legacy,
        sample_mode=args.sample_mode,
        mag_column=args.mag_column,
        mag_limit=args.mag_limit,
        q0=args.q0,
        exclude_legacy_types=args.exclude_types,
        is_legacy=True,
    )

    # 3. Simulate Null Cdf
    x_null, mean_null, lo_null, hi_null = legacy_null_cdf_envelope_from_shape_mc(
        legacy_cut, n_sample=n_frb, n_draws=args.mc_draws_null, q0=args.q0
    )

    y_steps = [0.0] + [i / n_frb for i in range(1, n_frb + 1)] + [1.0]
    x_mc, y_mc = mc_mean_cdf(frb_cosi_draws)

    # 4. Plot
    font_prop = font_manager.FontProperties(family="Arial", style="normal", size=8)
    fig, ax = plt.subplots(figsize=(8, 8))
    
    for draw in frb_cosi_draws:
        x_draw = [0.0] + draw + [1.0]
        ax.step(x_draw, y_steps, where="mid", color="red", linewidth=0.9, alpha=args.mc_alpha)

    ax.plot(x_mc, y_mc, color="red", linewidth=2.0, label="FRB host sample (MC, errors)")
    ax.plot(
        x_null,
        mean_null,
        color="#377eb8",
        linewidth=2.0,
        label=f"Legacy Tractor model $m_r$ null ({args.sample_mode})",
    )
    ax.fill_between(x_null, lo_null, hi_null, color="gray", alpha=0.35, label="68% confidence interval")
    ax.plot((0, 1), (0, 1), color="black", linestyle="--", linewidth=1.0, label="Uniform distribution")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("cos(i)", fontproperties=font_prop, fontsize=11)
    ax.set_ylabel("Cumulative distribution", fontproperties=font_prop, fontsize=11)
    ax.set_title(
        f"CDF Null: Legacy Tractor model $m_r$ ({args.sample_mode}, joint footprint)",
        fontsize=12,
    )
    ax.legend(fontsize=8, loc="upper left")

    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    tick_vals = np.cos(np.radians([90, 78, 66, 53, 37, 0]))
    ax_top.set_xticks(tick_vals)
    ax_top.set_xticklabels(["90", "78", "66", "53", "37", "0"], fontproperties=font_prop)
    ax_top.set_xlabel("Inclination angle i (degrees)", fontproperties=font_prop, fontsize=10)

    plt.tight_layout()
    out_pdf = os.path.join(args.out_dir, f"CDF_bias_legacy_{args.tag}.pdf")
    plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Plot saved to: {out_pdf} (N={n_frb})")

if __name__ == "__main__":
    main()
