"""
FRB inflated-error cos(i) CDF overlaid on SDSS null (mag21, full strict cuts).

SDSS: 10_000 draws of N_FRB galaxies from the strict late-type pool
(modelMag_r <= 21, u-r < 2.3, lnL exp-wins, expAB_r > q0) → median + 16–84%.

FRB: same inflated b/a MC as ``plot_frb_inflated_error_cdf.py``
(q ~ N(b/a, k*b_a_err), Hubble, +5° floor) → median + 16–84%.

Writes separate panels for k=3 and k=10.

Run from repo root::

    python scripts/plot_frb_vs_sdss_inflated_cdf.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (  # noqa: E402
    Q0,
    SDSS_Q_COLUMN_CDF,
    SDSS_UR_MAX_CDF,
    cosi_array_from_df,
    prepare_null_strict_color_base,
    read_sdss_null_catalog,
    slice_null_base_by_mag,
)
from pipeline_null_plot_utils import (  # noqa: E402
    DEFAULT_PIPELINE,
    DEFAULT_SDSS,
    PLOTS_NULL,
    add_inclination_top_axis,
    default_font,
    save_figure,
)
from plot_frb_inflated_error_cdf import (  # noqa: E402
    BA_ERR_INFLATE_DEFAULTS,
    INC_FLOOR_DEG,
    N_DRAWS,
    SEED,
    X,
    mc_inflated_cdf_envelope,
    select_hosts,
)

OUT_DIR = PLOTS_NULL / "v2" / "frb_inflated_error_cdf"


def mc_subsample_envelope(
    cosi: np.ndarray,
    *,
    n_sample: int,
    n_draws: int,
    x: np.ndarray,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Median + 16/84 ECDF envelope from resampling n_sample from a null pool."""
    ref = np.asarray(cosi, dtype=float)
    ref = ref[np.isfinite(ref)]
    if len(ref) < n_sample:
        raise RuntimeError(f"SDSS pool N={len(ref)} < n_sample={n_sample}")

    rng = np.random.default_rng(seed)
    y_corners = np.concatenate([[0.0], np.arange(1, n_sample + 1) / n_sample, [1.0]])
    ecdfs = np.empty((n_draws, len(x)), dtype=float)

    for i in range(n_draws):
        draw = np.sort(rng.choice(ref, size=n_sample, replace=False))
        x_corners = np.concatenate([[0.0], draw, [1.0]])
        for j in range(1, len(x_corners)):
            if x_corners[j] <= x_corners[j - 1]:
                x_corners[j] = x_corners[j - 1] + 1e-12
        f = interpolate.interp1d(
            x_corners,
            y_corners,
            kind="linear",
            bounds_error=False,
            fill_value=(0.0, 1.0),
        )
        ecdfs[i] = f(x)

    return (
        np.median(ecdfs, axis=0),
        np.percentile(ecdfs, 16.0, axis=0),
        np.percentile(ecdfs, 84.0, axis=0),
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sdss-csv", type=Path, default=DEFAULT_SDSS)
    parser.add_argument("--mag-limit", type=float, default=21.0)
    parser.add_argument("--n-draws", type=int, default=N_DRAWS)
    parser.add_argument("--q0", type=float, default=Q0)
    parser.add_argument(
        "--inflate",
        type=float,
        nargs="+",
        default=list(BA_ERR_INFLATE_DEFAULTS),
    )
    parser.add_argument("--floor-deg", type=float, default=INC_FLOOR_DEG)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()

    hosts = select_hosts(DEFAULT_PIPELINE, mag_limit=args.mag_limit, q0=args.q0)
    ba = hosts["b_a"].to_numpy(dtype=float)
    ba_err = hosts["b_a_err"].to_numpy(dtype=float)
    n_frb = len(hosts)

    print(f"[*] Loading SDSS null from {args.sdss_csv}", flush=True)
    sdss = read_sdss_null_catalog(args.sdss_csv)
    base = prepare_null_strict_color_base(
        sdss,
        mag_column="modelMag_r",
        q0=args.q0,
        q_column=SDSS_Q_COLUMN_CDF,
        is_legacy=False,
        sdss_ur_max=SDSS_UR_MAX_CDF,
        sdss_exp_winner_only=True,
    )
    sdss_cut = slice_null_base_by_mag(
        base, mag_column="modelMag_r", mag_limit=args.mag_limit
    )
    sdss_cosi = cosi_array_from_df(sdss_cut, q_col=SDSS_Q_COLUMN_CDF, q0=args.q0)
    print(
        f"[*] SDSS pool N={len(sdss_cosi):,}  "
        f"(modelMag_r<={args.mag_limit:g}, u-r<{SDSS_UR_MAX_CDF}, "
        f"lnL exp-wins, b/a>{args.q0:g})",
        flush=True,
    )
    print(f"[*] FRB N={n_frb}; MC draws={args.n_draws}", flush=True)

    print("[*] SDSS subsample envelope...", flush=True)
    sdss_med, sdss_lo, sdss_hi = mc_subsample_envelope(
        sdss_cosi,
        n_sample=n_frb,
        n_draws=args.n_draws,
        x=X,
        seed=SEED,
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    font_prop = default_font()
    sdss_col = "#4daf4a"
    frb_col = "#d62728"

    for inflate in args.inflate:
        tag = f"x{int(inflate) if float(inflate).is_integer() else inflate:g}"
        print(f"[*] FRB inflate={inflate:g}x overlay...", flush=True)
        frb_med, frb_lo, frb_hi = mc_inflated_cdf_envelope(
            ba,
            ba_err,
            n_draws=args.n_draws,
            inflate=inflate,
            floor_deg=args.floor_deg,
            q0=args.q0,
            x=X,
            seed=SEED,
        )

        fig, ax = plt.subplots(figsize=(8, 8))
        # Uniform — no legend
        ax.plot((0, 1), (0, 1), "k--", lw=1.2, zorder=1)

        ax.fill_between(
            X, sdss_lo, sdss_hi, color=sdss_col, alpha=0.22, linewidth=0, zorder=2
        )
        ax.plot(
            X,
            sdss_med,
            color=sdss_col,
            lw=2.0,
            zorder=3,
            label=rf"SDSS null  ($N={n_frb}\times{args.n_draws // 1000}$k, 68% CI)",
        )

        ax.fill_between(
            X, frb_lo, frb_hi, color=frb_col, alpha=0.22, linewidth=0, zorder=4
        )
        ax.plot(
            X,
            frb_med,
            color=frb_col,
            lw=2.0,
            zorder=5,
            label=(
                rf"FRB hosts  ($N={n_frb}$, "
                rf"${inflate:g}\times\sigma_{{b/a}}+{args.floor_deg:g}^\circ$)"
            ),
        )

        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel(
            rf"$\cos(i)$ (Hubble, $q_0={args.q0:g}$)",
            fontproperties=font_prop,
            fontsize=11,
        )
        ax.set_ylabel("Cumulative distribution", fontproperties=font_prop, fontsize=11)
        ax.set_title(
            rf"FRB vs SDSS $\cos(i)$ CDF  "
            rf"($m_r\leq{args.mag_limit:g}$, strict+color, "
            rf"${inflate:g}\times\sigma_{{b/a}}+{args.floor_deg:g}^\circ$)",
            fontsize=11,
        )
        ax.legend(loc="upper left", fontsize=8, frameon=True)
        ax.grid(True, alpha=0.3)
        add_inclination_top_axis(ax, font_prop)
        fig.tight_layout()

        stem = args.out_dir / f"frb_vs_sdss_inflated_cdf_{tag}"
        save_figure(fig, stem)
        plt.close(fig)


if __name__ == "__main__":
    main()
