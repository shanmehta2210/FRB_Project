import argparse
import math
import os
import random
import re

import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate


random.seed(42)
np.random.seed(42)


def hubble_cosi_from_ba(b_over_a, q0=0.2):
    val = (b_over_a**2 - q0**2) / (1 - q0**2)
    if val < 0:
        return 0.0
    if val > 1:
        return 1.0
    return math.sqrt(val)


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


def _parse_float(value):
    if value is None:
        return np.nan
    if isinstance(value, (float, int, np.floating, np.integer)):
        return float(value)
    text = str(value).strip()
    if text == "":
        return np.nan
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", text)
    if not m:
        return np.nan
    return float(m.group(0))


def build_frb_mc_draws(galfit_df, q0=0.2, n_draws=500):
    draws = []
    mu_vals = []
    for _, row in galfit_df.iterrows():
        mu = _parse_float(row.get("b_a_psf"))
        sigma = _parse_float(row.get("b_a_err_psf"))
        if not np.isfinite(mu):
            continue
        if not np.isfinite(sigma) or sigma < 0:
            sigma = 0.0
        mu = min(1.0, max(0.0, mu))
        mu_vals.append(mu)

    if len(mu_vals) == 0:
        raise RuntimeError("No valid GALFIT b/a values found for Monte Carlo sampling.")

    for _ in range(n_draws):
        cosi = []
        for _, row in galfit_df.iterrows():
            mu = _parse_float(row.get("b_a_psf"))
            sigma = _parse_float(row.get("b_a_err_psf"))
            if not np.isfinite(mu):
                continue
            if not np.isfinite(sigma) or sigma < 0:
                sigma = 0.0
            sampled = np.random.normal(mu, sigma) if sigma > 0 else mu
            sampled = min(1.0, max(0.0, float(sampled)))
            cosi.append(hubble_cosi_from_ba(sampled, q0=q0))
        draws.append(sorted(cosi))

    return draws


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


def sdss_cosi_by_band(sdss, band_col, mag_limit=21.0, use_b_proxy=False):
    vals = []
    for _, row in sdss.iterrows():
        try:
            expab = float(row["expAB_r"])
            if use_b_proxy:
                gmag = float(row["gmag"])
                rmag = float(row["rmag"])
                mag = gmag + 0.3130 * (gmag - rmag) + 0.2271
            else:
                mag = float(row[band_col])
        except Exception:
            continue

        if mag <= mag_limit and expab > 0.2:
            vals.append(hubble_cosi_from_ba(expab))

    return vals


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdss-csv", default="SDSS_catalog.csv")
    parser.add_argument("--master-summary", default="master_frb_summary.csv")
    parser.add_argument("--galfit-csv", default="galfit_sigma_metrics_summary.csv")
    parser.add_argument("--null-bands", default="rgb", help="Combination of null lines to plot: r, g, b (e.g. r, rg, rgb)")
    parser.add_argument("--show-ci", choices=["auto", "on", "off"], default="auto")
    parser.add_argument("--mag-limit", type=float, default=21.0)
    parser.add_argument("--mc-draws", type=int, default=500)
    parser.add_argument("--mc-alpha", type=float, default=0.03)
    parser.add_argument("--out-dir", default="plots/plots_multiband_cdf")
    parser.add_argument("--tag", default="rgb_bands")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    sdss = pd.read_csv(args.sdss_csv)
    master = pd.read_csv(args.master_summary)
    galfit = pd.read_csv(args.galfit_csv)

    master_frbs = set(master["FRB"].astype(str))
    galfit = galfit[galfit["FRB"].astype(str).isin(master_frbs)].copy()

    frb_cosi_draws = build_frb_mc_draws(galfit, n_draws=args.mc_draws)
    n_frb = len(frb_cosi_draws[0])
    if n_frb == 0:
        raise RuntimeError("No FRB galaxies available for Monte Carlo CDF.")

    selected = [b for b in ["r", "g", "b"] if b in args.null_bands.lower()]
    if not selected:
        raise RuntimeError("No null bands selected. Use --null-bands with any of r, g, b.")

    band_data = {}
    if "r" in selected:
        band_data["r"] = sdss_cosi_by_band(sdss, "rmag", mag_limit=args.mag_limit)
    if "g" in selected:
        band_data["g"] = sdss_cosi_by_band(sdss, "gmag", mag_limit=args.mag_limit)
    if "b" in selected:
        band_data["b"] = sdss_cosi_by_band(sdss, "", mag_limit=args.mag_limit, use_b_proxy=True)

    if min(len(v) for v in band_data.values()) < n_frb:
        raise RuntimeError("At least one selected band has fewer galaxies than FRB sample size.")

    show_ci = args.show_ci == "on" or (args.show_ci == "auto" and len(selected) == 1)

    cdf_lines = {}
    for band in selected:
        x_band, mean_band, lo_band, hi_band = cdf_envelope(band_data[band], n_sample=n_frb)
        cdf_lines[band] = (x_band, mean_band, lo_band, hi_band)

    y_steps = [0.0] + [i / n_frb for i in range(1, n_frb + 1)] + [1.0]
    x_mc, y_mc = mc_mean_cdf(frb_cosi_draws)

    font_prop = font_manager.FontProperties(family="Arial", style="normal", size=8)

    fig, ax = plt.subplots(figsize=(8, 8))
    for draw in frb_cosi_draws:
        x_draw = [0.0] + draw + [1.0]
        ax.step(x_draw, y_steps, where="mid", color="red", linewidth=0.9, alpha=args.mc_alpha)
    ax.plot(x_mc, y_mc, color="red", linewidth=2.0, label="FRB host sample (MC, errors)")
    if "r" in cdf_lines:
        x_r, mean_r, lo_r, hi_r = cdf_lines["r"]
        ax.plot(x_r, mean_r, color="#377eb8", linewidth=1.8, label=f"r-band null ($m_r < {args.mag_limit:g}$)")
    if "g" in cdf_lines:
        x_g, mean_g, lo_g, hi_g = cdf_lines["g"]
        ax.plot(x_g, mean_g, color="#4daf4a", linewidth=1.8, label=f"g-band null ($m_g < {args.mag_limit:g}$)")
    if "b" in cdf_lines:
        x_b, mean_b, lo_b, hi_b = cdf_lines["b"]
        ax.plot(x_b, mean_b, color="#984ea3", linewidth=1.8, label=f"b-band null ($B < {args.mag_limit:g}$)")
    if show_ci:
        ci_band = selected[0]
        x_ci, _, lo_ci, hi_ci = cdf_lines[ci_band]
        ax.fill_between(x_ci, lo_ci, hi_ci, color="gray", alpha=0.3, label="68% confidence interval")
    ax.plot((0, 1), (0, 1), color="black", linestyle="--", linewidth=1.0, label="Uniform distribution")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("cos(i)", fontproperties=font_prop, fontsize=11)
    ax.set_ylabel("Cumulative distribution", fontproperties=font_prop, fontsize=11)
    ax.set_title("CDF Null Comparison: r/g/b Bands", fontsize=12)
    ax.legend(fontsize=8, loc="upper left")

    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    tick_vals = np.cos(np.radians([90, 78, 66, 53, 37, 0]))
    ax_top.set_xticks(tick_vals)
    ax_top.set_xticklabels(["90", "78", "66", "53", "37", "0"], fontproperties=font_prop)
    ax_top.set_xlabel("Inclination angle i (degrees)", fontproperties=font_prop, fontsize=10)

    plt.tight_layout()
    out_pdf = os.path.join(args.out_dir, f"CDF_bias_multiband_{args.tag}.pdf")
    out_png = os.path.join(args.out_dir, f"CDF_bias_multiband_{args.tag}.png")
    plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"FRB sample size: {n_frb}")
    print(f"FRB MC draws: {args.mc_draws}")
    print("Selected null bands:", ",".join(selected))
    print("Show CI:", show_ci)
    print("SDSS sample sizes:", {k: len(v) for k, v in band_data.items()})
    print(f"Saved {out_pdf}")
    print(f"Saved {out_png}")


if __name__ == "__main__":
    main()
