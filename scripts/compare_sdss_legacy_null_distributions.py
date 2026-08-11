import argparse
import os
import sys
from pathlib import Path

import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate
from scipy.stats import ks_2samp, wasserstein_distance

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (
    cosi_array_from_df,
    LEGACY_GR_MAX_CDF,
    SDSS_UR_MAX_CDF,
    prepare_null_sample,
)

np.random.seed(42)

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SDSS = REPO_ROOT / "catalog/SDSS_catalog_v1_allsky_modelmr.csv"
DEFAULT_LEGACY = REPO_ROOT / "catalog/LS_catalog_v1_allsky_modelmr.csv"
DEFAULT_OUT_SUBDIR = "plots/plots_null/v1_null_plots"


def cdf_envelope(reference_vals: np.ndarray, n_sample: int, n_draws: int = 10000):
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


def summary_stats(values: np.ndarray, label: str) -> dict:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    return {
        "sample": label,
        "n": int(len(arr)),
        "mean": float(np.mean(arr)),
        "median": float(np.median(arr)),
        "std": float(np.std(arr, ddof=1)),
        "min": float(np.min(arr)),
        "q10": float(np.quantile(arr, 0.10)),
        "q25": float(np.quantile(arr, 0.25)),
        "q75": float(np.quantile(arr, 0.75)),
        "q90": float(np.quantile(arr, 0.90)),
        "max": float(np.max(arr)),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sdss-csv", default=str(DEFAULT_SDSS))
    parser.add_argument("--legacy-csv", default=str(DEFAULT_LEGACY))
    parser.add_argument("--master-summary", default="master_frb_summary.csv")
    parser.add_argument("--galfit-csv", default="Archive/csv/galfit/galfit_metrics_summary.csv")
    parser.add_argument("--mag-limit", type=float, default=21.0)
    parser.add_argument("--mag-column", default="rmag", help="Model r: rmag (SDSS) / tractor_mag_r (Legacy).")
    parser.add_argument(
        "--sample-mode",
        choices=("strict",),
        default="strict",
        help="strict: b/a > q0 before CDF (inclusive mode removed from plots).",
    )
    parser.add_argument(
        "--sdss-mag-column",
        default="modelMag_r",
        help="SDSS brightness column for mag cut (best single-profile r).",
    )
    parser.add_argument(
        "--sdss-q-column",
        default="expAB_r",
        help="SDSS axis ratio for Hubble i (default expAB_r after lnL exp-winner cut).",
    )
    parser.add_argument("--q0", type=float, default=0.2)
    parser.add_argument("--mc-draws", type=int, default=10000)
    parser.add_argument(
        "--exclude-types",
        default="REX",
        help="Comma-separated Tractor types to drop from Legacy null (empty string = keep REX).",
    )
    parser.add_argument(
        "--include-rex",
        action="store_true",
        help="Include Legacy REX (sets --exclude-types to empty).",
    )
    parser.add_argument("--out-dir", default=DEFAULT_OUT_SUBDIR)
    parser.add_argument("--tag", default="v1_allsky_modelmr")
    args = parser.parse_args()
    if args.include_rex:
        args.exclude_types = ""

    os.makedirs(args.out_dir, exist_ok=True)

    sdss = pd.read_csv(args.sdss_csv)
    legacy = pd.read_csv(args.legacy_csv)
    master = pd.read_csv(args.master_summary)
    galfit = pd.read_csv(args.galfit_csv)

    master_frbs = set(master["FRB"].astype(str))
    galfit = galfit[galfit["FRB"].astype(str).isin(master_frbs)].copy()
    n_frb = int(galfit["b_a_psf"].notna().sum())
    if n_frb <= 0:
        raise RuntimeError("No FRB sample size found from galfit_sigma_metrics_summary.csv")

    sdss_cut = prepare_null_sample(
        sdss,
        sample_mode=args.sample_mode,
        mag_column=args.sdss_mag_column,
        mag_limit=args.mag_limit,
        q0=args.q0,
        q_column=args.sdss_q_column,
        is_legacy=False,
        sdss_ur_max=SDSS_UR_MAX_CDF,
    )
    legacy_cut = prepare_null_sample(
        legacy,
        sample_mode=args.sample_mode,
        mag_column=args.mag_column,
        mag_limit=args.mag_limit,
        q0=args.q0,
        exclude_legacy_types=args.exclude_types,
        is_legacy=True,
        legacy_gr_max=LEGACY_GR_MAX_CDF,
    )

    sdss_cosi = cosi_array_from_df(sdss_cut, q_col=args.sdss_q_column, q0=args.q0)
    legacy_cosi = cosi_array_from_df(legacy_cut, q0=args.q0)

    if len(sdss_cosi) < n_frb or len(legacy_cosi) < n_frb:
        raise RuntimeError("Null pool is smaller than FRB sample size.")

    x_s, mean_s, lo_s, hi_s = cdf_envelope(sdss_cosi, n_sample=n_frb, n_draws=args.mc_draws)
    x_l, mean_l, lo_l, hi_l = cdf_envelope(legacy_cosi, n_sample=n_frb, n_draws=args.mc_draws)

    font_prop = font_manager.FontProperties(family="Arial", style="normal", size=8)
    fig, ax = plt.subplots(figsize=(8, 8))

    mode_lbl = "q>q0" if args.sample_mode == "strict" else "inclusive"
    sdss_q_lbl = (
        "exp disk $b/a$"
        if args.sdss_q_column == "expAB_r"
        else "best deV/exp $b/a$"
    )
    ax.plot(
        x_s,
        mean_s,
        color="#377eb8",
        linewidth=2.0,
        label=f"SDSS {sdss_q_lbl} ({mode_lbl}, $m<{args.mag_limit:g}$)",
    )
    ax.fill_between(x_s, lo_s, hi_s, color="#377eb8", alpha=0.20, label="SDSS 68% CI")

    ax.plot(
        x_l,
        mean_l,
        color="#e41a1c",
        linewidth=2.0,
        label=(
            f"Legacy Tractor $m_r$ null ({mode_lbl}, incl. REX)"
            if args.include_rex or not str(args.exclude_types).strip()
            else f"Legacy Tractor $m_r$ null ({mode_lbl}, no REX)"
        ),
    )
    ax.fill_between(x_l, lo_l, hi_l, color="#e41a1c", alpha=0.18, label="Legacy 68% CI")

    ax.plot((0, 1), (0, 1), color="black", linestyle="--", linewidth=1.0, label="Uniform distribution")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("cos(i)", fontproperties=font_prop, fontsize=11)
    ax.set_ylabel("Cumulative distribution", fontproperties=font_prop, fontsize=11)
    ax.set_title(
        f"Null CDF: SDSS vs Legacy model $m_r$ ({args.sample_mode}, joint footprint)",
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
    out_png = os.path.join(args.out_dir, f"null_overlay_{args.tag}.png")
    out_pdf = os.path.join(args.out_dir, f"null_overlay_{args.tag}.pdf")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)

    s_sdss = summary_stats(sdss_cosi, "SDSS")
    s_legacy = summary_stats(legacy_cosi, "Legacy_no_REX")

    stats_df = pd.DataFrame([s_sdss, s_legacy])
    stats_out = os.path.join(args.out_dir, f"null_distribution_stats_{args.tag}.csv")
    stats_df.to_csv(stats_out, index=False)

    ks = ks_2samp(sdss_cosi, legacy_cosi)
    wass = wasserstein_distance(sdss_cosi, legacy_cosi)

    test_df = pd.DataFrame(
        [
            {
                "test": "ks_2samp",
                "statistic": float(ks.statistic),
                "pvalue": float(ks.pvalue),
            },
            {
                "test": "wasserstein",
                "statistic": float(wass),
                "pvalue": np.nan,
            },
            {
                "test": "mean_diff_legacy_minus_sdss",
                "statistic": float(s_legacy["mean"] - s_sdss["mean"]),
                "pvalue": np.nan,
            },
            {
                "test": "median_diff_legacy_minus_sdss",
                "statistic": float(s_legacy["median"] - s_sdss["median"]),
                "pvalue": np.nan,
            },
        ]
    )
    tests_out = os.path.join(args.out_dir, f"null_distribution_tests_{args.tag}.csv")
    test_df.to_csv(tests_out, index=False)

    print(f"FRB sample size used for envelopes: {n_frb}")
    print(f"SDSS rows after cut: {len(sdss_cut)}")
    print(f"Legacy rows after cut/filter: {len(legacy_cut)}")
    print(f"Saved {out_png}")
    print(f"Saved {out_pdf}")
    print(f"Saved {stats_out}")
    print(f"Saved {tests_out}")


if __name__ == "__main__":
    main()
