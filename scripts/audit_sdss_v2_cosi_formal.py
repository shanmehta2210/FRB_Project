#!/usr/bin/env python3
"""Formal cos(i) hypothesis tests and SDSS magnitude benchmark (v2 audit)."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

_SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (  # noqa: E402
    Q0,
    SDSS_UR_MAX_CDF,
    apply_strict_q_cut,
    cosi_array_from_df,
    filter_sdss_drop_dev_winners,
    filter_sdss_ur,
    prepare_null_strict_color_base,
    read_sdss_null_catalog,
)
from pipeline_null_plot_utils import DEFAULT_SDSS_V2, REPO_ROOT  # noqa: E402

OUT_DIR = REPO_ROOT / "plots" / "plots_null" / "v2_sdss_audit" / "formal"
REF_DIR = OUT_DIR / "reference"
MODELMAG_REF_CSV = REF_DIR / "sdss_modelmag_r_counts.csv"
YASUDA_CSV = REF_DIR / "yasuda2001_r_counts.csv"
MAG_LIMITS = (19.0, 20.0, 21.0, 22.0)
# Typical r* − modelMag_r offset for SDSS pipeline galaxies (late-type mix).
PETRO_TO_MODELMAG_DELTA = 0.055


def bootstrap_stat(
    values: np.ndarray,
    stat_fn,
    n_boot: int = 10_000,
    seed: int = 42,
    batch: int = 64,
) -> tuple[float, float, float]:
    rng = np.random.default_rng(seed)
    n = len(values)
    draws = np.empty(n_boot, dtype=float)
    for start in range(0, n_boot, batch):
        b = min(batch, n_boot - start)
        idx = rng.integers(0, n, size=(b, n))
        chunk = values[idx]
        if stat_fn is np.median:
            draws[start : start + b] = np.median(chunk, axis=1)
        elif stat_fn is np.mean:
            draws[start : start + b] = np.mean(chunk, axis=1)
        else:
            draws[start : start + b] = [stat_fn(row) for row in chunk]
    point = float(stat_fn(values))
    return point, float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


def permutation_mag_cut_pvalue(
    mag: np.ndarray,
    cosi: np.ndarray,
    mag_limit: float,
    observed_median: float,
    n_perm: int = 5_000,
    seed: int = 42,
) -> float:
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(n_perm):
        perm_mag = rng.permutation(mag)
        mask = perm_mag <= mag_limit
        if mask.sum() == 0:
            continue
        if np.median(cosi[mask]) >= observed_median:
            count += 1
    return (count + 1) / (n_perm + 1)


def spearman_with_ci(x: np.ndarray, y: np.ndarray, n_boot: int = 500, seed: int = 42) -> dict:
    rho, p = stats.spearmanr(x, y)
    rng = np.random.default_rng(seed)
    n = len(x)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, size=n)
        boots[i], _ = stats.spearmanr(x[idx], y[idx])
    return {
        "rho": float(rho),
        "p_value": float(p),
        "ci_lo": float(np.percentile(boots, 2.5)),
        "ci_hi": float(np.percentile(boots, 97.5)),
    }


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> dict:
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    rz = stats.rankdata(z)
    ry_res = ry - np.polyval(np.polyfit(rz, ry, 1), rz)
    rx_res = rx - np.polyval(np.polyfit(rz, rx, 1), rz)
    rho, p = stats.spearmanr(rx_res, ry_res)
    return {"rho_partial": float(rho), "p_value": float(p)}


def weighted_median(values: np.ndarray, weights: np.ndarray) -> float:
    order = np.argsort(values)
    v = values[order]
    w = weights[order]
    cw = np.cumsum(w) / np.sum(w)
    idx = int(np.searchsorted(cw, 0.5))
    return float(v[min(idx, len(v) - 1)])


def mag_bin_table(mag: np.ndarray, cosi: np.ndarray, step: float = 0.5) -> pd.DataFrame:
    rows = []
    lo = 15.0
    while lo < 28.0:
        hi = lo + step
        if lo == 15.0:
            mask = mag <= hi
        else:
            mask = (mag > lo) & (mag <= hi)
        n = int(mask.sum())
        if n < 30:
            lo = hi
            continue
        rows.append(
            {
                "mag_lo": lo,
                "mag_hi": hi,
                "n": n,
                "median_cosi": float(np.median(cosi[mask])),
                "mean_cosi": float(np.mean(cosi[mask])),
                "frac_pool": n / len(mag),
            }
        )
        lo = hi
    return pd.DataFrame(rows)


def load_pools() -> tuple[pd.DataFrame, pd.DataFrame]:
    df = read_sdss_null_catalog(DEFAULT_SDSS_V2)
    base = prepare_null_strict_color_base(
        df,
        mag_column="modelMag_r",
        q0=Q0,
        q_column="expAB_r",
        is_legacy=False,
        sdss_ur_max=SDSS_UR_MAX_CDF,
        sdss_exp_winner_only=True,
    )
    return df, base


def run_hypothesis_tests(base: pd.DataFrame, n_boot: int, n_perm: int) -> pd.DataFrame:
    mag = pd.to_numeric(base["modelMag_r"], errors="coerce").to_numpy()
    cosi = cosi_array_from_df(base, q_col="expAB_r", q0=Q0)
    rows = []

    ks_u, p_u = stats.kstest(cosi, "uniform")
    one_minus = 1.0 - cosi
    hist, _ = np.histogram(one_minus, bins=15, range=(0, 1))
    _, p_chi2 = stats.chisquare(hist)

    med, med_lo, med_hi = bootstrap_stat(cosi, np.median, n_boot=n_boot)
    rows.append(
        {
            "test": "pre_mag_isotropy",
            "stage": "production_no_mag_cut",
            "n": len(cosi),
            "statistic": float(med),
            "stat_lo": med_lo,
            "stat_hi": med_hi,
            "p_value": float(p_u),
            "note": f"KS vs U(0,1) p={p_u:.4g}; chi2 flat(1-cosi) p={p_chi2:.4g}",
        }
    )

    full_median = float(np.median(cosi))
    for mlim in MAG_LIMITS:
        mask = mag <= mlim
        sub = cosi[mask]
        med_o = float(np.median(sub))
        mean_o = float(np.mean(sub))
        med_b, lo_b, hi_b = bootstrap_stat(sub, np.median, n_boot=n_boot)
        mean_b, mean_lo, mean_hi = bootstrap_stat(sub, np.mean, n_boot=n_boot)
        delta = med_o - full_median
        p_perm = permutation_mag_cut_pvalue(mag, cosi, mlim, med_o, n_perm=n_perm)

        rows.append(
            {
                "test": "mag_cut_median",
                "stage": f"mag_le_{mlim:g}",
                "n": int(mask.sum()),
                "statistic": med_o,
                "stat_lo": lo_b,
                "stat_hi": hi_b,
                "p_value": p_perm,
                "note": f"perm p={p_perm:.4g}; delta_median_vs_full={delta:+.4f}",
            }
        )
        rows.append(
            {
                "test": "mag_cut_mean",
                "stage": f"mag_le_{mlim:g}",
                "n": int(mask.sum()),
                "statistic": mean_o,
                "stat_lo": mean_lo,
                "stat_hi": mean_hi,
                "p_value": p_perm,
                "note": "bootstrap mean; same permutation as median",
            }
        )
    return pd.DataFrame(rows)


def run_correlation_tests(base: pd.DataFrame) -> pd.DataFrame:
    mag = pd.to_numeric(base["modelMag_r"], errors="coerce").to_numpy()
    ba = pd.to_numeric(base["expAB_r"], errors="coerce").to_numpy()
    cosi = cosi_array_from_df(base, q_col="expAB_r", q0=Q0)
    ok = np.isfinite(mag) & np.isfinite(ba) & np.isfinite(cosi)

    s_mag = spearman_with_ci(mag[ok], cosi[ok])
    s_ba = spearman_with_ci(ba[ok], cosi[ok])
    s_part = partial_spearman(mag[ok], cosi[ok], ba[ok])
    slope_mag, intercept, r_lin, p_lin, _ = stats.linregress(mag[ok], cosi[ok])
    slope_ba, intercept_ba, r_ba, p_ba, _ = stats.linregress(ba[ok], cosi[ok])

    X = np.column_stack([np.ones(ok.sum()), mag[ok], ba[ok]])
    coef, _, _, _ = np.linalg.lstsq(X, cosi[ok], rcond=None)

    rows = [
        {"relation": "spearman_mag_cosi", **s_mag},
        {"relation": "spearman_expAB_cosi", **s_ba},
        {"relation": "partial_spearman_mag_cosi_given_expAB", **s_part},
        {
            "relation": "ols_slope_cosi_vs_mag",
            "rho": float(slope_mag),
            "p_value": float(p_lin),
            "ci_lo": float(intercept),
            "ci_hi": float(r_lin),
        },
        {
            "relation": "ols_slope_cosi_vs_expAB",
            "rho": float(slope_ba),
            "p_value": float(p_ba),
            "ci_lo": float(intercept_ba),
            "ci_hi": float(r_ba),
        },
        {
            "relation": "ols_cosi_intercept_mag_coef_expAB_coef",
            "rho": float(coef[1]),
            "p_value": np.nan,
            "ci_lo": float(coef[2]),
            "ci_hi": float(coef[0]),
        },
    ]
    return pd.DataFrame(rows)


def run_mixture_and_sim(base: pd.DataFrame, seed: int = 42) -> tuple[pd.DataFrame, pd.DataFrame]:
    mag = pd.to_numeric(base["modelMag_r"], errors="coerce").to_numpy()
    cosi = cosi_array_from_df(base, q_col="expAB_r", q0=Q0)
    bins = mag_bin_table(mag, cosi)
    rng = np.random.default_rng(seed)

    w_uniform = np.ones_like(mag, dtype=float)
    for _, brow in bins.iterrows():
        if brow["mag_lo"] == 15.0:
            mask = mag <= brow["mag_hi"]
        else:
            mask = (mag > brow["mag_lo"]) & (mag <= brow["mag_hi"])
        if mask.sum():
            w_uniform[mask] = 1.0 / mask.sum()

    obs21 = float(np.median(cosi[mag <= 21]))
    mix_rows = [
        {"scenario": "observed_mag_le_21", "median_cosi": obs21, "note": "direct"},
        {
            "scenario": "uniform_mag_bin_reweight_mag_le_21",
            "median_cosi": weighted_median(cosi[mag <= 21], w_uniform[mag <= 21]),
            "note": "equal weight per 0.5 mag bin then median",
        },
        {
            "scenario": "bin_median_composition_mag_le_21",
            "median_cosi": weighted_median(
                bins.loc[bins["mag_hi"] <= 21, "median_cosi"].to_numpy(),
                bins.loc[bins["mag_hi"] <= 21, "n"].to_numpy(),
            ),
            "note": "weighted mean of per-bin medians",
        },
    ]

    shuffled = rng.permutation(cosi)
    sim_cosi = rng.uniform(0, 1, size=50_000)
    sim_rows = [
        {
            "null": "B_isotropic_uniform_cosi",
            "median_full": float(np.median(sim_cosi)),
            "median_mag21": float(np.median(sim_cosi)),
            "observed_mag21": obs21,
        },
        {
            "null": "A_shuffled_cosi_vs_mag",
            "median_full": float(np.median(shuffled)),
            "median_mag21": float(np.median(shuffled[mag <= 21])),
            "observed_mag21": obs21,
        },
        {
            "null": "observed",
            "median_full": float(np.median(cosi)),
            "median_mag21": obs21,
            "observed_mag21": obs21,
        },
    ]
    return pd.DataFrame(mix_rows), pd.DataFrame(sim_rows)


def cut_survival_by_mag(df: pd.DataFrame) -> pd.DataFrame:
    """Per-bin cut survival on the 1.9M catalog (boolean masks computed once)."""
    mag = pd.to_numeric(df["modelMag_r"], errors="coerce").to_numpy()
    after_ur = filter_sdss_ur(df, SDSS_UR_MAX_CDF)
    after_lnl = filter_sdss_drop_dev_winners(after_ur)
    after_strict = apply_strict_q_cut(after_lnl, q_col="expAB_r", q0=Q0)

    ur_mask = np.zeros(len(df), dtype=bool)
    lnl_mask = np.zeros(len(df), dtype=bool)
    strict_mask = np.zeros(len(df), dtype=bool)
    pos = {idx: i for i, idx in enumerate(df.index)}
    for idx in after_ur.index:
        ur_mask[pos[idx]] = True
    for idx in after_lnl.index:
        lnl_mask[pos[idx]] = True
    for idx in after_strict.index:
        strict_mask[pos[idx]] = True

    rows = []
    for lo in np.arange(15, 24.5, 0.5):
        hi = lo + 0.5
        if lo == 15.0:
            mask = mag <= hi
        else:
            mask = (mag > lo) & (mag <= hi)
        n_raw = int(mask.sum())
        if n_raw < 50:
            continue
        rows.append(
            {
                "mag_lo": lo,
                "mag_hi": hi,
                "n_catalog": n_raw,
                "frac_pass_ur": float(ur_mask[mask].mean()),
                "frac_pass_lnl": float(lnl_mask[mask].mean()),
                "frac_pass_strict": float(strict_mask[mask].mean()),
            }
        )
    return pd.DataFrame(rows)


def _mag_bin_counts(mag: pd.Series, ref: pd.DataFrame) -> pd.Series:
    """Count objects in each half-mag bin defined by ref[['mag_lo','mag_hi']]."""
    m = pd.to_numeric(mag, errors="coerce")
    counts = []
    for _, row in ref.iterrows():
        lo, hi = row["mag_lo"], row["mag_hi"]
        counts.append(int(((m > lo) & (m <= hi)).sum()))
    return pd.Series(counts, index=ref.index, dtype=int)


def compare_mag_sky_validation(df: pd.DataFrame, base: pd.DataFrame) -> pd.DataFrame:
    """Compare strict-pool N(modelMag_r) to full catalog and literature sky counts."""
    ref = pd.read_csv(MODELMAG_REF_CSV)
    mag_pool = pd.to_numeric(base["modelMag_r"], errors="coerce")
    mag_raw = pd.to_numeric(df["modelMag_r"], errors="coerce")

    comp = ref.copy()
    comp["n_pool"] = _mag_bin_counts(mag_pool, ref).to_numpy()
    comp["n_raw_catalog"] = _mag_bin_counts(mag_raw, ref).to_numpy()
    comp["frac_pool"] = comp["n_pool"] / comp["n_pool"].sum()
    comp["frac_raw_catalog"] = comp["n_raw_catalog"] / comp["n_raw_catalog"].sum()
    comp["frac_sky_ref"] = (
        comp["n_per_deg2_per_half_mag"] / comp["n_per_deg2_per_half_mag"].sum()
    )
    comp["ratio_pool_over_raw"] = comp["frac_pool"] / comp["frac_raw_catalog"]
    comp["ratio_pool_over_sky_ref"] = comp["frac_pool"] / comp["frac_sky_ref"]

    if YASUDA_CSV.is_file():
        yas = pd.read_csv(YASUDA_CSV)
        yas_map = yas.set_index(
            yas["mag_lo"].astype(str) + "_" + yas["mag_hi"].astype(str)
        )["n_per_deg2_per_half_mag"]
        key = comp["mag_lo"].astype(str) + "_" + comp["mag_hi"].astype(str)
        comp["yasuda_petrosian_n_per_deg2"] = key.map(yas_map)
        comp["frac_yasuda_petrosian"] = (
            comp["yasuda_petrosian_n_per_deg2"]
            / comp["yasuda_petrosian_n_per_deg2"].sum()
        )

    obs = comp["n_pool"].to_numpy(dtype=float)
    exp_raw = obs.sum() * (
        comp["frac_raw_catalog"].to_numpy() / comp["frac_raw_catalog"].sum()
    )
    chi2_raw, p_raw = stats.chisquare(obs, exp_raw)
    exp_sky = obs.sum() * (
        comp["frac_sky_ref"].to_numpy() / comp["frac_sky_ref"].sum()
    )
    chi2_sky, p_sky = stats.chisquare(obs, exp_sky)
    comp.attrs["chi2_vs_raw"] = float(chi2_raw)
    comp.attrs["chi2_p_vs_raw"] = float(p_raw)
    comp.attrs["chi2_vs_sky_ref"] = float(chi2_sky)
    comp.attrs["chi2_p_vs_sky_ref"] = float(p_sky)
    comp.attrs["chi2_dof"] = len(obs) - 1
    return comp


def plot_mag_counts(comp: pd.DataFrame, out_png: Path) -> None:
    x = 0.5 * (comp["mag_lo"] + comp["mag_hi"])
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(x, comp["frac_pool"], "o-", label="Strict null pool (normalized)")
    ax.plot(x, comp["frac_raw_catalog"], "^-", alpha=0.85, label="Full v2 catalog (normalized)")
    ax.plot(
        x,
        comp["frac_sky_ref"],
        "s--",
        label="SDSS sky ref: modelMag_r (Yasuda+Stoughton offset)",
    )
    ax.set_xlabel("modelMag_r (bin center)")
    ax.set_ylabel("Fraction per 0.5 mag bin")
    ax.set_title("Sky magnitude distribution validation (modelMag_r)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_cut_survival(surv: pd.DataFrame, out_png: Path) -> None:
    x = 0.5 * (surv["mag_lo"] + surv["mag_hi"])
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x, surv["frac_pass_ur"], label="pass u-r")
    ax.plot(x, surv["frac_pass_lnl"], label="pass lnL exp-wins")
    ax.plot(x, surv["frac_pass_strict"], label="pass expAB_r > 0.2")
    ax.set_xlabel("modelMag_r")
    ax.set_ylabel("Survival fraction")
    ax.set_title("Cut efficiency vs magnitude (from 1.9M catalog)")
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(0, 1.05)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_joint_panel(bins: pd.DataFrame, out_png: Path) -> None:
    x = 0.5 * (bins["mag_lo"] + bins["mag_hi"])
    fig, ax1 = plt.subplots(figsize=(9, 5))
    ax1.plot(x, bins["median_cosi"], "o-", color="C0", label="median cos(i)")
    ax1.axhline(0.5, color="0.5", ls=":", label="isotropic 0.5")
    ax1.set_xlabel("modelMag_r")
    ax1.set_ylabel("median cos(i)", color="C0")
    ax1.tick_params(axis="y", labelcolor="C0")
    ax2 = ax1.twinx()
    ax2.plot(x, bins["frac_pool"], "s--", color="C1", alpha=0.8, label="N(m) fraction")
    ax2.set_ylabel("pool fraction per bin", color="C1")
    ax2.tick_params(axis="y", labelcolor="C1")
    ax1.set_title("cos(i) vs magnitude composition (production pool, no mag cut)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def write_report_stub(out_md: Path, hyp: pd.DataFrame) -> None:
    """Do not overwrite the narrative FORMAL_COSI_AUDIT.md if it already exists."""
    if out_md.exists():
        return
    pre = hyp[hyp["test"] == "pre_mag_isotropy"].iloc[0]
    mag21 = hyp[(hyp["test"] == "mag_cut_median") & (hyp["stage"] == "mag_le_21")].iloc[0]
    out_md.write_text(
        f"# Formal cos(i) audit\n\nPre-mag median={pre['statistic']:.4f}; "
        f"mag≤21={mag21['statistic']:.4f}. Expand this file with full analysis.\n",
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-bootstrap", type=int, default=2_000)
    parser.add_argument("--n-perm", type=int, default=1_000)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading pools...")
    df, base = load_pools()
    mag = pd.to_numeric(base["modelMag_r"], errors="coerce").to_numpy()
    cosi = cosi_array_from_df(base, q_col="expAB_r", q0=Q0)
    bins = mag_bin_table(mag, cosi)

    print("A1-A2: hypothesis tests...")
    hyp = run_hypothesis_tests(base, n_boot=args.n_bootstrap, n_perm=args.n_perm)
    hyp.to_csv(args.out_dir / "cosi_hypothesis_tests.csv", index=False)

    print("A3: correlations...")
    corr = run_correlation_tests(base)
    corr.to_csv(args.out_dir / "cosi_mag_correlation.csv", index=False)

    print("A4-A5: mixture + simulation...")
    mix, sim = run_mixture_and_sim(base)
    mix.to_csv(args.out_dir / "cosi_mixture_decomposition.csv", index=False)
    sim.to_csv(args.out_dir / "cosi_simulation_null.csv", index=False)

    print("B: sky N(m) validation + cut survival...")
    comp = compare_mag_sky_validation(df, base)
    comp.to_csv(args.out_dir / "mag_counts_sky_validation.csv", index=False)
    surv = cut_survival_by_mag(df)
    surv.to_csv(args.out_dir / "cut_survival_vs_mag.csv", index=False)

    print("Plots...")
    plot_mag_counts(comp, args.out_dir / "mag_counts_comparison.png")
    plot_cut_survival(surv, args.out_dir / "cut_survival_vs_mag.png")
    plot_joint_panel(bins, args.out_dir / "cosi_mag_joint_panel.png")
    bins.to_csv(args.out_dir / "cosi_per_mag_bin.csv", index=False)

    write_report_stub(args.out_dir / "FORMAL_COSI_AUDIT.md", hyp)

    print(hyp.to_string(index=False))
    print(
        f"\nchi2 vs raw catalog: {comp.attrs['chi2_vs_raw']:.2f}, "
        f"p={comp.attrs['chi2_p_vs_raw']:.4g}"
    )
    print(
        f"chi2 vs modelMag sky ref: {comp.attrs['chi2_vs_sky_ref']:.2f}, "
        f"p={comp.attrs['chi2_p_vs_sky_ref']:.4g}"
    )
    print(f"Wrote outputs to {args.out_dir}")


if __name__ == "__main__":
    main()
