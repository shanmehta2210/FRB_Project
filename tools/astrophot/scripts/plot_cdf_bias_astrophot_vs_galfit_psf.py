import math
import os
import random
import argparse
import re
import sys
from pathlib import Path

import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate

REPO_ROOT = Path(__file__).resolve().parents[3]
_SCRIPTS = REPO_ROOT / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import cosi_array_from_df, hubble_cosi_from_ba, prepare_null_sample

DEFAULT_SDSS = REPO_ROOT / "catalog/SDSS_catalog_v1_allsky_modelmr.csv"
DEFAULT_OUT_SUBDIR = "plots/plots_astrophot_psf_sigma/v1_allsky_modelmr"

random.seed(42)
np.random.seed(42)


def build_sdss_reference(sdss_csv, mag_column="rmag", mag_limit=21.0, sample_mode="strict", q0=0.2):
    sdss = pd.read_csv(sdss_csv)
    cut = prepare_null_sample(
        sdss,
        sample_mode=sample_mode,
        mag_column=mag_column,
        mag_limit=mag_limit,
        q0=q0,
        is_legacy=False,
    )
    return cosi_array_from_df(cut, q0=q0).tolist()


def cdf_envelope(reference_vals, n_sample, n_draws=10000):
    total_samples = []
    for _ in range(n_draws):
        draw = sorted(random.sample(reference_vals, n_sample))
        idx_norm = [0.0] + [j / n_sample for j in range(1, n_sample + 1)] + [1.0]
        draw_ext = [0.0] + draw + [1.0]
        total_samples.append(interpolate.interp1d(draw_ext, idx_norm))

    x = np.linspace(0, 1, 100)
    means, down, up = [], [], []
    lo = int(0.16 * n_draws)
    hi = int(0.84 * n_draws)
    for value in x:
        idx_sample = sorted(float(s(value)) for s in total_samples)
        means.append(float(np.mean(idx_sample)))
        down.append(float(idx_sample[lo]))
        up.append(float(idx_sample[hi]))

    return x, np.array(means), np.array(down), np.array(up)


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


def build_mc_cosi_draws(shared, ba_map, ba_err_map, n_draws=500):
    draws = []
    for _ in range(n_draws):
        vals = []
        for frb in shared:
            mu = parse_float(ba_map.get(frb))
            sig = parse_float(ba_err_map.get(frb))
            if not np.isfinite(mu):
                continue
            if not np.isfinite(sig) or sig < 0:
                sig = 0.0
            sampled = np.random.normal(mu, sig) if sig > 0 else mu
            sampled = min(1.0, max(0.0, float(sampled)))
            vals.append(hubble_cosi_from_ba(sampled))
        draws.append(sorted(vals))
    return draws


def mean_cdf_from_draws(draws):
    n_sample = len(draws[0])
    funcs = []
    for draw in draws:
        x_draw = [0.0] + draw + [1.0]
        idx_norm = [0.0] + [j / n_sample for j in range(1, n_sample + 1)] + [1.0]
        funcs.append(interpolate.interp1d(x_draw, idx_norm))

    x_grid = np.linspace(0, 1, 100)
    y_mean = [np.mean([float(f(v)) for f in funcs]) for v in x_grid]
    return x_grid, y_mean


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--astrophot-csv",
        default="AstroPhot_Analysis/results/astrophot_psf_sigma_inclination_angles.csv",
    )
    parser.add_argument("--galfit-csv", default="galfit_sigma_metrics_summary.csv")
    parser.add_argument("--frb-sample", default="frb_sample.txt")
    parser.add_argument("--sdss-csv", default=str(DEFAULT_SDSS))
    parser.add_argument("--mag-limit", type=float, default=21.0)
    parser.add_argument("--mag-column", default="rmag")
    parser.add_argument(
        "--sample-mode",
        choices=("strict",),
        default="strict",
    )
    parser.add_argument("--q0", type=float, default=0.2)
    parser.add_argument("--out-dir", default=DEFAULT_OUT_SUBDIR)
    parser.add_argument("--tag", default="v1_allsky_modelmr")
    parser.add_argument("--astrophot-label", default="AstroPhot (PSF + sigma)")
    parser.add_argument("--title", default="CDF Bias Comparison: AstroPhot vs GALFIT (PSF + sigma)")
    parser.add_argument("--mc-draws", type=int, default=500)
    parser.add_argument("--mc-alpha", type=float, default=0.03)
    args = parser.parse_args()

    astrophot_csv = args.astrophot_csv
    galfit_csv = args.galfit_csv
    frb_sample_path = args.frb_sample
    sdss_csv = args.sdss_csv
    out_dir = args.out_dir
    os.makedirs(out_dir, exist_ok=True)

    ap = pd.read_csv(astrophot_csv)
    gf = pd.read_csv(galfit_csv)
    frb_sample = pd.read_csv(frb_sample_path, sep="\t")

    frb_list = frb_sample["FRB"].astype(str).tolist()

    ap_inc = dict(zip(ap["frb_name"], ap["inclination_angle"]))
    ap_q = {
        str(row["frb_name"]): parse_float(row.get("q"))
        for _, row in ap.iterrows()
    }
    ap_q_err = {
        str(row["frb_name"]): parse_float(row.get("error_q"))
        for _, row in ap.iterrows()
    }
    gf_ba = {
        str(row["FRB"]): float(str(row["b_a_psf"]).replace("*", ""))
        for _, row in gf.iterrows()
    }
    gf_ba_err = {
        str(row["FRB"]): parse_float(row.get("b_a_err_psf"))
        for _, row in gf.iterrows()
    }

    shared = [f for f in frb_list if f in ap_inc and f in gf_ba]
    if len(shared) == 0:
        raise RuntimeError("No shared FRBs found between AstroPhot and GALFIT tables.")

    ap_cosi = sorted([math.cos(math.radians(float(ap_inc[f]))) for f in shared])
    gf_cosi = sorted([hubble_cosi_from_ba(gf_ba[f]) for f in shared])

    n_frb = len(shared)
    sdss_vals = build_sdss_reference(
        sdss_csv,
        mag_column=args.mag_column,
        mag_limit=args.mag_limit,
        sample_mode=args.sample_mode,
        q0=args.q0,
    )
    x_ref, mean_cdf, lo_cdf, hi_cdf = cdf_envelope(sdss_vals, n_sample=n_frb)

    y_steps = [0.0] + [i / n_frb for i in range(1, n_frb + 1)] + [1.0]
    x_ap = [0.0] + ap_cosi + [1.0]
    x_gf = [0.0] + gf_cosi + [1.0]

    font_prop = font_manager.FontProperties(family="Arial", style="normal", size=8)

    fig, ax = plt.subplots(figsize=(7, 5))
    ap_mc_samples = build_mc_cosi_draws(shared, ap_q, ap_q_err, n_draws=args.mc_draws)
    for draw in ap_mc_samples:
        x_draw = [0.0] + draw + [1.0]
        ax.step(x_draw, y_steps, where="mid", color="red", linewidth=0.9, alpha=args.mc_alpha)
    x_ap_mc, y_ap_mc = mean_cdf_from_draws(ap_mc_samples)
    ax.plot(x_ap_mc, y_ap_mc, color="red", linewidth=2.0, label=f"{args.astrophot_label} (MC errors)")

    gf_mc_samples = build_mc_cosi_draws(shared, gf_ba, gf_ba_err, n_draws=args.mc_draws)
    for draw in gf_mc_samples:
        x_draw = [0.0] + draw + [1.0]
        ax.step(x_draw, y_steps, where="mid", color="blue", linewidth=0.9, alpha=args.mc_alpha)
    x_gf_mc, y_gf_mc = mean_cdf_from_draws(gf_mc_samples)
    ax.plot(x_gf_mc, y_gf_mc, color="blue", linewidth=2.0, label="GALFIT (PSF + sigma, MC errors)")
    ax.plot(x_ref, mean_cdf, color="black", linewidth=1.2, label="SDSS distribution ($m_r < 21$)")
    ax.fill_between(x_ref, lo_cdf, hi_cdf, color="gray", alpha=0.35, label="68% confidence interval")
    ax.plot((0, 1), (0, 1), color="black", linestyle="--", linewidth=1.0, label="Uniform distribution")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("cos(i)", fontproperties=font_prop, fontsize=10)
    ax.set_ylabel("Cumulative distribution", fontproperties=font_prop, fontsize=10)
    ax.set_title(args.title, fontsize=11)
    ax.legend(fontsize=8)

    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    tick_vals = np.cos(np.radians([90, 78, 66, 53, 37, 0]))
    ax_top.set_xticks(tick_vals)
    ax_top.set_xticklabels(["90", "78", "66", "53", "37", "0"], fontproperties=font_prop)
    ax_top.set_xlabel("Inclination angle i (degrees)", fontproperties=font_prop, fontsize=10)

    plt.tight_layout()
    out_pdf = os.path.join(out_dir, f"CDF_bias_comparison_{args.tag}.pdf")
    out_png = os.path.join(out_dir, f"CDF_bias_comparison_{args.tag}.png")
    plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Shared FRBs: {n_frb}")
    print(f"Saved {out_pdf}")
    print(f"Saved {out_png}")


if __name__ == "__main__":
    main()
