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

REPO_ROOT = Path(__file__).resolve().parents[3]
_SCRIPTS = REPO_ROOT / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import cosi_array_from_df, prepare_null_sample

DEFAULT_SDSS = REPO_ROOT / "catalog/SDSS_catalog_v1_allsky_modelmr.csv"
DEFAULT_OUT_SUBDIR = "plots/plots_photutils/v1_allsky_modelmr"


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


def hubble_cosi_from_ba(b_over_a, q0=0.2):
    val = (b_over_a**2 - q0**2) / (1 - q0**2)
    if val < 0:
        return 0.0
    if val > 1:
        return 1.0
    return math.sqrt(val)


def cdf_envelope_null(reference_vals, n_sample, n_draws=10000, seed=42):
    rng = np.random.default_rng(seed)
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

    x = np.linspace(0.0, 1.0, 500)
    mean = []
    lo = []
    hi = []
    for value in x:
        vals = sorted(float(f(value)) for f in total_samples)
        mean.append(float(np.mean(vals)))
        lo.append(float(np.percentile(vals, 16)))
        hi.append(float(np.percentile(vals, 84)))

    mean = np.array(mean)
    lo = np.array(lo)
    hi = np.array(hi)
    return x, mean, lo, hi


def cdf_envelope_mc_measured(inc_deg, inc_err_deg, n_draws=10000, seed=42):
    rng = np.random.default_rng(seed)

    inc_deg = np.asarray(inc_deg, dtype=float)
    inc_err_deg = np.asarray(inc_err_deg, dtype=float)
    valid = np.isfinite(inc_deg) & np.isfinite(inc_err_deg)
    inc_deg = inc_deg[valid]
    inc_err_deg = inc_err_deg[valid]

    if len(inc_deg) == 0:
        raise RuntimeError("No valid measured inclinations for MC envelope.")

    n_sample = len(inc_deg)
    x = np.linspace(0.0, 1.0, 500)
    cdf_evals = np.zeros((n_draws, len(x)))

    for i in range(n_draws):
        mc_inc = rng.normal(inc_deg, inc_err_deg)
        mc_inc = np.clip(mc_inc, 0.0, 90.0)
        mc_cosi = np.cos(np.radians(mc_inc))
        mc_cosi = np.sort(mc_cosi)
        cdf_evals[i, :] = np.searchsorted(mc_cosi, x, side="right") / n_sample

    mean = np.mean(cdf_evals, axis=0)
    lo = np.percentile(cdf_evals, 16, axis=0)
    hi = np.percentile(cdf_evals, 84, axis=0)
    return x, mean, lo, hi


def _parse_float(value):
    if value is None:
        return np.nan
    if isinstance(value, (float, int, np.floating, np.integer)):
        return float(value)
    text = str(value).strip()
    if text == "":
        return np.nan
    return pd.to_numeric(text, errors="coerce")


def build_method_mc_draws(df, mu_col, sigma_col, q0=0.2, n_draws=500):
    draws = []
    has_valid = False

    for _ in range(n_draws):
        cosi = []
        for _, row in df.iterrows():
            mu = _parse_float(row.get(mu_col))
            sigma = _parse_float(row.get(sigma_col))
            if not np.isfinite(mu):
                continue
            if not np.isfinite(sigma) or sigma < 0:
                sigma = 0.0
            sampled = np.random.normal(mu, sigma) if sigma > 0 else mu
            sampled = min(1.0, max(0.0, float(sampled)))
            cosi.append(hubble_cosi_from_ba(sampled, q0=q0))

        if len(cosi) > 0:
            has_valid = True
        draws.append(sorted(cosi))

    if not has_valid:
        raise RuntimeError(f"No valid values found in {mu_col} for Monte Carlo sampling.")
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
    median = []
    for value in x:
        vals = [float(f(value)) for f in funcs]
        mean.append(float(np.mean(vals)))
        median.append(float(np.median(vals)))
    return x, np.array(mean), np.array(median)


def sample_cosi_draws_measured(inc_deg, inc_err_deg, n_draws=1000, seed=42):
    rng = np.random.default_rng(seed)

    inc_deg = np.asarray(inc_deg, dtype=float)
    inc_err_deg = np.asarray(inc_err_deg, dtype=float)
    valid = np.isfinite(inc_deg) & np.isfinite(inc_err_deg)
    inc_deg = inc_deg[valid]
    inc_err_deg = inc_err_deg[valid]

    if len(inc_deg) == 0:
        raise RuntimeError("No valid measured inclinations for MC sampling.")

    draws = []
    for _ in range(n_draws):
        mc_inc = rng.normal(inc_deg, inc_err_deg)
        mc_inc = np.clip(mc_inc, 0.0, 90.0)
        mc_cosi = np.sort(np.cos(np.radians(mc_inc)))
        draws.append(mc_cosi)
    return draws


def clean_err(err, default_floor=0.0):
    if not np.isfinite(err):
        return default_floor
    return max(default_floor, float(err))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--photutils-csv", default="tools/photutils/results/photutils_ellipse_inclination_angles_nopsf.csv")
    parser.add_argument("--galfit-csv", default="galfit_sigma_metrics_summary.csv")
    parser.add_argument("--sdss-csv", default=str(DEFAULT_SDSS))
    parser.add_argument("--master-csv", default="master_frb_summary.csv")
    parser.add_argument("--mag-limit", type=float, default=21.0)
    parser.add_argument("--mag-column", default="rmag")
    parser.add_argument(
        "--sample-mode",
        choices=("strict",),
        default="strict",
    )
    parser.add_argument("--q0", type=float, default=0.2)
    parser.add_argument("--mc-draws", type=int, default=10000)
    parser.add_argument("--mc-draws-spaghetti", type=int, default=500)
    parser.add_argument("--spaghetti-alpha", type=float, default=0.04)
    parser.add_argument("--out-dir", default=DEFAULT_OUT_SUBDIR)
    parser.add_argument("--tag", default="v1_allsky_modelmr")
    args = parser.parse_args()

    np.random.seed(42)

    os.makedirs(args.out_dir, exist_ok=True)

    phot = pd.read_csv(args.photutils_csv)
    galfit = pd.read_csv(args.galfit_csv)
    master = pd.read_csv(args.master_csv)
    sdss = pd.read_csv(args.sdss_csv)

    # Photutils sample
    phot_ok = phot[phot["status"].astype(str) == "ok"].copy()
    phot_ok["inclination_deg"] = pd.to_numeric(phot_ok["inclination_deg"], errors="coerce")
    phot_ok["inclination_err_deg"] = pd.to_numeric(phot_ok["inclination_err_deg"], errors="coerce")
    phot_ok = phot_ok[np.isfinite(phot_ok["inclination_deg"])].copy()
    phot_ok["inclination_err_deg"] = phot_ok["inclination_err_deg"].apply(clean_err)

    # Build GALFIT inclinations (prefer provided inc_psf, else derive from b_a_psf)
    gal = galfit.copy()
    gal["b_a_psf"] = pd.to_numeric(gal.get("b_a_psf"), errors="coerce")
    gal["b_a_err_psf"] = pd.to_numeric(gal.get("b_a_err_psf"), errors="coerce")

    if "inc_psf" in gal.columns:
        gal["inclination_deg"] = pd.to_numeric(gal["inc_psf"], errors="coerce")
    else:
        gal["inclination_deg"] = gal["b_a_psf"].apply(lambda q: hubble_inclination_from_q(q, q0=args.q0))

    if "inc_err_psf" in gal.columns:
        gal["inclination_err_deg"] = pd.to_numeric(gal["inc_err_psf"], errors="coerce")
    else:
        gal["inclination_err_deg"] = gal.apply(
            lambda r: hubble_inclination_err_from_qerr(r.get("b_a_psf"), r.get("b_a_err_psf"), q0=args.q0),
            axis=1,
        )

    gal["inclination_err_deg"] = gal["inclination_err_deg"].apply(clean_err)

    # Keep FRBs from master and shared with photutils
    master_frb = set(master["FRB"].astype(str))
    phot_ok = phot_ok[phot_ok["FRB"].astype(str).isin(master_frb)].copy()
    gal = gal[gal["FRB"].astype(str).isin(master_frb)].copy()

    shared = sorted(set(phot_ok["FRB"].astype(str)).intersection(set(gal["FRB"].astype(str))))
    if len(shared) == 0:
        raise RuntimeError("No shared FRBs between Photutils and GALFIT after filtering.")

    phot_shared = phot_ok[phot_ok["FRB"].astype(str).isin(shared)].copy()
    gal_shared = gal[gal["FRB"].astype(str).isin(shared)].copy()

    # Align by FRB for deterministic sampling
    phot_shared = phot_shared.sort_values("FRB").reset_index(drop=True)
    gal_shared = gal_shared.sort_values("FRB").reset_index(drop=True)

    # Ensure both methods have valid b/a values for exactly the same FRBs (multiband MC policy).
    phot_shared["b_over_a"] = pd.to_numeric(phot_shared["b_over_a"], errors="coerce")
    phot_shared["b_over_a_err"] = pd.to_numeric(phot_shared["b_over_a_err"], errors="coerce")
    gal_shared["b_a_psf"] = pd.to_numeric(gal_shared["b_a_psf"], errors="coerce")
    gal_shared["b_a_err_psf"] = pd.to_numeric(gal_shared["b_a_err_psf"], errors="coerce")

    valid_phot = set(phot_shared.loc[np.isfinite(phot_shared["b_over_a"]), "FRB"].astype(str))
    valid_gal = set(gal_shared.loc[np.isfinite(gal_shared["b_a_psf"]), "FRB"].astype(str))
    shared_ba = sorted(valid_phot.intersection(valid_gal))
    if len(shared_ba) == 0:
        raise RuntimeError("No shared FRBs with valid b/a in both Photutils and GALFIT.")

    phot_shared = phot_shared[phot_shared["FRB"].astype(str).isin(shared_ba)].sort_values("FRB").reset_index(drop=True)
    gal_shared = gal_shared[gal_shared["FRB"].astype(str).isin(shared_ba)].sort_values("FRB").reset_index(drop=True)

    sdss_cut = prepare_null_sample(
        sdss,
        sample_mode=args.sample_mode,
        mag_column=args.mag_column,
        mag_limit=args.mag_limit,
        q0=args.q0,
        is_legacy=False,
    )
    sdss_cosi = cosi_array_from_df(sdss_cut, q0=args.q0)

    n_frb = len(shared_ba)
    x_null, mean_null, lo_null, hi_null = cdf_envelope_null(sdss_cosi, n_sample=n_frb, n_draws=args.mc_draws)

    # Photutils/GALFIT MC from b/a + b/a_err (same policy as multiband null script).
    phot_draws = build_method_mc_draws(
        phot_shared,
        mu_col="b_over_a",
        sigma_col="b_over_a_err",
        q0=args.q0,
        n_draws=args.mc_draws_spaghetti,
    )
    gal_draws = build_method_mc_draws(
        gal_shared,
        mu_col="b_a_psf",
        sigma_col="b_a_err_psf",
        q0=args.q0,
        n_draws=args.mc_draws_spaghetti,
    )
    x_p, mean_p, med_p = mc_mean_cdf(phot_draws)
    x_g, mean_g, med_g = mc_mean_cdf(gal_draws)

    y_step = np.arange(1, n_frb + 1) / n_frb

    font_prop = font_manager.FontProperties(family="Arial", style="normal", size=8)

    fig, ax = plt.subplots(figsize=(8, 8))

    # SDSS null
    ax.plot(x_null, mean_null, color="#377eb8", linewidth=2.0, label=f"SDSS null ($m_r < {args.mag_limit:g}$)")
    ax.fill_between(x_null, lo_null, hi_null, color="#377eb8", alpha=0.20, label="SDSS 68% CI")

    # GALFIT (multiband-style): draw traces + mean.
    for draw in gal_draws:
        ax.step(
            np.concatenate([[0.0], draw, [1.0]]),
            np.concatenate([[0.0], y_step, [1.0]]),
            where="mid",
            color="#1b4f72",
            linewidth=0.9,
            alpha=args.spaghetti_alpha,
        )
    ax.plot(x_g, med_g, color="#1b4f72", linewidth=2.0, linestyle="--", label="GALFIT MC mean")

    # Photutils (same mechanism as GALFIT): draw traces + mean.
    for draw in phot_draws:
        ax.step(
            np.concatenate([[0.0], draw, [1.0]]),
            np.concatenate([[0.0], y_step, [1.0]]),
            where="mid",
            color="#e41a1c",
            linewidth=0.9,
            alpha=args.spaghetti_alpha,
        )
    ax.plot(x_p, med_p, color="#e41a1c", linewidth=2.0, linestyle="-", label="Photutils MC mean")

    ax.plot((0, 1), (0, 1), color="black", linestyle=":", linewidth=1.0, label="Uniform distribution")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("cos(i)", fontproperties=font_prop, fontsize=11)
    ax.set_ylabel("Cumulative distribution", fontproperties=font_prop, fontsize=11)
    ax.set_title("CDF Overlay: Photutils vs GALFIT vs SDSS null", fontsize=12)
    ax.legend(fontsize=8, loc="upper left")

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
    plt.close(fig)

    # Secondary plot (same policy, reduced clutter title for side-by-side comparison).
    fig2, ax2 = plt.subplots(figsize=(8, 8))

    ax2.plot(x_null, mean_null, color="#377eb8", linewidth=2.0, label=f"SDSS null ($m_r < {args.mag_limit:g}$)")
    ax2.fill_between(x_null, lo_null, hi_null, color="#377eb8", alpha=0.20, label="SDSS 68% CI")

    for draw in gal_draws:
        ax2.step(
            np.concatenate([[0.0], draw, [1.0]]),
            np.concatenate([[0.0], y_step, [1.0]]),
            where="post",
            color="#1b4f72",
            linewidth=0.9,
            alpha=args.spaghetti_alpha,
        )
    ax2.plot(x_g, med_g, color="#1b4f72", linewidth=2.0, linestyle="--", label="GALFIT MC mean")

    for draw in phot_draws:
        ax2.step(
            np.concatenate([[0.0], draw, [1.0]]),
            np.concatenate([[0.0], y_step, [1.0]]),
            where="post",
            color="#e41a1c",
            linewidth=0.9,
            alpha=args.spaghetti_alpha,
        )
    ax2.plot(x_p, med_p, color="#e41a1c", linewidth=2.0, linestyle="-", label="Photutils MC mean")

    ax2.plot((0, 1), (0, 1), color="black", linestyle=":", linewidth=1.0, label="Uniform distribution")

    ax2.set_xlim(0, 1)
    ax2.set_ylim(0, 1)
    ax2.set_xlabel("cos(i)", fontproperties=font_prop, fontsize=11)
    ax2.set_ylabel("Cumulative distribution", fontproperties=font_prop, fontsize=11)
    ax2.set_title("CDF Overlay (MC draw traces): Photutils vs GALFIT vs SDSS null", fontsize=12)
    ax2.legend(fontsize=8, loc="upper left")

    ax2_top = ax2.twiny()
    ax2_top.set_xlim(ax2.get_xlim())
    ax2_top.set_xticks(tick_vals)
    ax2_top.set_xticklabels(["90", "78", "66", "53", "37", "0"], fontproperties=font_prop)
    ax2_top.set_xlabel("Inclination angle i (degrees)", fontproperties=font_prop, fontsize=10)

    plt.tight_layout()
    out_spaghetti_png = os.path.join(args.out_dir, f"{args.tag}_spaghetti.png")
    out_spaghetti_pdf = os.path.join(args.out_dir, f"{args.tag}_spaghetti.pdf")
    fig2.savefig(out_spaghetti_png, dpi=300, bbox_inches="tight")
    fig2.savefig(out_spaghetti_pdf, bbox_inches="tight")
    plt.close(fig2)

    print(f"Shared FRBs used: {n_frb}")
    print(f"Photutils rows (ok): {len(phot_ok)}")
    print(f"GALFIT rows (usable): {len(gal_shared)}")
    print(f"SDSS rows after cut: {len(sdss_cut)}")
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")
    print(f"Saved: {out_spaghetti_png}")
    print(f"Saved: {out_spaghetti_pdf}")


if __name__ == "__main__":
    main()
