#!/usr/bin/env python3
"""Bootstrap median cos(i) vs subsample size for SDSS v2 strict mag<=21 pool."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (
    Q0,
    SDSS_UR_MAX_CDF,
    cosi_array_from_df,
    prepare_null_strict_color_base,
    read_sdss_null_catalog,
    slice_null_base_by_mag,
)
from pipeline_null_plot_utils import DEFAULT_SDSS_V2, REPO_ROOT

OUT_DIR = REPO_ROOT / "plots" / "plots_null" / "v2_sdss_audit" / "formal"
OUT_CSV = OUT_DIR / "cosi_subsample_stability.csv"
OUT_PNG = OUT_DIR / "cosi_subsample_stability.png"

MAG_LIMIT = 21.0
MIN_SUBSAMPLE = 2000
DEFAULT_N_BOOTSTRAP = 500
DEFAULT_SEED = 42


def load_mag21_pool() -> tuple[np.ndarray, int]:
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
    pool = slice_null_base_by_mag(base, mag_column="modelMag_r", mag_limit=MAG_LIMIT)
    cosi = cosi_array_from_df(pool, q_col="expAB_r", q0=Q0)
    return cosi, len(pool)


def default_subsample_sizes(n_pool: int, min_n: int = MIN_SUBSAMPLE) -> np.ndarray:
    lo = min(min_n, n_pool)
    if n_pool <= lo:
        return np.array([n_pool], dtype=int)
    sizes = np.unique(np.round(np.logspace(np.log10(lo), np.log10(n_pool), 22)).astype(int))
    if sizes[-1] != n_pool:
        sizes = np.append(sizes, n_pool)
    return sizes


def parse_sizes_arg(sizes_str: str | None, n_pool: int) -> np.ndarray:
    if not sizes_str:
        return default_subsample_sizes(n_pool)
    vals = [int(x.strip()) for x in sizes_str.split(",") if x.strip()]
    sizes = np.unique(np.clip(vals, 1, n_pool))
    if sizes[-1] != n_pool:
        sizes = np.append(sizes, n_pool)
    return sizes.astype(int)


def bootstrap_median_cosi(
    cosi: np.ndarray,
    sizes: np.ndarray,
    *,
    n_bootstrap: int,
    seed: int,
) -> pd.DataFrame:
    n_pool = len(cosi)
    full_median = float(np.median(cosi))
    rng = np.random.default_rng(seed)
    rows: list[dict] = []

    for n in sizes:
        n = int(n)
        medians = np.empty(n_bootstrap, dtype=float)
        for b in range(n_bootstrap):
            idx = rng.choice(n_pool, size=n, replace=False)
            medians[b] = np.median(cosi[idx])
        rows.append(
            {
                "n_subsample": n,
                "n_bootstrap": n_bootstrap,
                "mean_median_cosi": float(np.mean(medians)),
                "std_median_cosi": float(np.std(medians, ddof=1)),
                "p16_median_cosi": float(np.percentile(medians, 16)),
                "p84_median_cosi": float(np.percentile(medians, 84)),
                "full_pool_median_cosi": full_median,
                "n_pool": n_pool,
            }
        )
    return pd.DataFrame(rows)


def plot_subsample_stability(df: pd.DataFrame, out_png: Path) -> None:
    x = df["n_subsample"].to_numpy(dtype=float)
    y = df["mean_median_cosi"].to_numpy()
    p16 = df["p16_median_cosi"].to_numpy()
    p84 = df["p84_median_cosi"].to_numpy()
    full_med = float(df["full_pool_median_cosi"].iloc[0])

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.fill_between(x, p16, p84, alpha=0.25, color="C0", label="16–84% of subsample medians")
    ax.plot(x, y, "o-", color="C0", lw=1.5, ms=4, label="Mean of subsample medians")
    ax.axhline(full_med, color="C1", ls="--", lw=1.2, label=f"Full pool median ({full_med:.3f})")
    ax.axhline(0.5, color="0.5", ls=":", lw=1.0, label="Isotropic cos(i)=0.5")

    ax.set_xscale("log")
    ax.set_xlabel("Subsample size n")
    ax.set_ylabel("cos(i)")
    ax.set_title(
        f"SDSS v2 strict pool (modelMag_r ≤ {MAG_LIMIT:g}): "
        f"median cos(i) vs subsample size\n"
        f"B={int(df['n_bootstrap'].iloc[0])} bootstrap draws per n; "
        f"N_pool={int(df['n_pool'].iloc[0]):,}"
    )
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.3, which="both")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-bootstrap", type=int, default=DEFAULT_N_BOOTSTRAP)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument(
        "--sizes",
        type=str,
        default=None,
        help="Comma-separated subsample sizes (default: log-spaced 2000..N)",
    )
    parser.add_argument("--out-csv", type=Path, default=OUT_CSV)
    parser.add_argument("--out-png", type=Path, default=OUT_PNG)
    args = parser.parse_args()

    cosi, n_pool = load_mag21_pool()
    sizes = parse_sizes_arg(args.sizes, n_pool)
    print(f"Pool N={n_pool:,}; full-pool median cos(i)={np.median(cosi):.4f}")
    print(f"Bootstrap B={args.n_bootstrap}; {len(sizes)} sample sizes: {sizes.min()}..{sizes.max()}")

    result = bootstrap_median_cosi(
        cosi,
        sizes,
        n_bootstrap=args.n_bootstrap,
        seed=args.seed,
    )
    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.out_csv, index=False)
    plot_subsample_stability(result, args.out_png)

    last = result.iloc[-1]
    print(result.to_string(index=False))
    print(f"\nWrote {args.out_csv}")
    print(f"Wrote {args.out_png}")
    print(
        f"Sanity @ n={int(last['n_subsample']):,}: "
        f"mean_median={last['mean_median_cosi']:.4f} "
        f"(full pool {last['full_pool_median_cosi']:.4f})"
    )


if __name__ == "__main__":
    main()
