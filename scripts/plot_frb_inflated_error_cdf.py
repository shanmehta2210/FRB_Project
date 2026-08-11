"""
Error-aware FRB host cos(i) CDF with inflated GALFIT uncertainties.

For each of 10_000 draws, every host is resampled as:
  q ~ N(b/a, k * b_a_err)          # formal GALFIT 1σ inflated by k (default: 3 and 10)
  i = Hubble(q; q0=0.2)
  i = i + N(0, 5°)                  # photometric inclination floor
  cos(i) = cos(i)

Each draw builds a continuous ECDF (linear interp through order statistics),
then the plot shows the median CDF and 16–84% (68% CI) band across draws —
same style as ``plot_jimin_cdf_overlay_mc.py``.

Sample: pipeline hosts with ``mag_final <= 21``, ``b/a > q0``, excluding weak
associations (``WEAK_ASSOCIATIONS_PRODUCTION67.md`` §A/§B).

Run from repo root::

    python scripts/plot_frb_inflated_error_cdf.py
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import Q0  # noqa: E402
from pipeline_null_plot_utils import (  # noqa: E402
    DEFAULT_PIPELINE,
    PLOTS_NULL,
    add_inclination_top_axis,
    default_font,
    frb_hosts_for_cdf,
    load_pipeline_hosts,
    save_figure,
)

# Weak / Verdi / no-mag hosts — not for secure inclination CDFs.
WEAK_FRBS = frozenset(
    {
        "20210214G",
        "20221116A",
        "20230913",
        "20240104A",
        "20240203",
        "20250518",
    }
)

N_DRAWS = 10_000
SEED = 42
BA_ERR_INFLATE_DEFAULTS = (3.0, 10.0)
INC_FLOOR_DEG = 5.0
X = np.linspace(0.0, 1.0, 200)


def _hubble_inc_deg(q: np.ndarray, q0: float = Q0) -> np.ndarray:
    """Vectorized Hubble b/a -> i (degrees)."""
    val = (q**2 - q0**2) / (1.0 - q0**2)
    val = np.clip(val, 0.0, 1.0)
    return np.degrees(np.arccos(np.sqrt(val)))


def mc_inflated_cdf_envelope(
    ba: np.ndarray,
    ba_err: np.ndarray,
    *,
    n_draws: int = N_DRAWS,
    inflate: float = 10.0,
    floor_deg: float = INC_FLOOR_DEG,
    q0: float = Q0,
    x: np.ndarray = X,
    seed: int = SEED,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Median + 16/84 CDF envelope from inflated b/a sampling + i floor.

    Returns (median, lo, hi) evaluated on ``x``.
    """
    ba = np.asarray(ba, dtype=float)
    ba_err = np.asarray(ba_err, dtype=float)
    n = len(ba)
    if n == 0:
        raise RuntimeError("No hosts for MC CDF.")

    sigma_q = np.maximum(inflate * np.nan_to_num(ba_err, nan=0.0), 0.0)
    rng = np.random.default_rng(seed)

    y_corners = np.concatenate([[0.0], np.arange(1, n + 1) / n, [1.0]])
    ecdfs = np.empty((n_draws, len(x)), dtype=float)

    for i in range(n_draws):
        q = rng.normal(ba, sigma_q)
        q = np.clip(q, q0 + 1e-6, 1.0)
        inc = _hubble_inc_deg(q, q0=q0)
        inc = inc + rng.normal(0.0, floor_deg, size=n)
        inc = np.clip(inc, 0.0, 90.0)
        cosi = np.sort(np.cos(np.radians(inc)))

        x_corners = np.concatenate([[0.0], cosi, [1.0]])
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

    med = np.median(ecdfs, axis=0)
    lo = np.percentile(ecdfs, 16.0, axis=0)
    hi = np.percentile(ecdfs, 84.0, axis=0)
    return med, lo, hi


def select_hosts(
    pipeline_csv: Path,
    *,
    q0: float = Q0,
    mag_limit: float = 21.0,
    mag_column: str = "mag_final",
) -> pd.DataFrame:
    hosts = load_pipeline_hosts(pipeline_csv)
    hosts = frb_hosts_for_cdf(
        hosts,
        sample_mode="strict",
        q0=q0,
        ba_col="b_a",
        mag_limit=mag_limit,
        frb_mag_column=mag_column,
    )
    hosts = hosts.loc[~hosts["frb"].astype(str).isin(WEAK_FRBS)].copy()
    if hosts.empty:
        raise RuntimeError(
            "No FRB hosts left after mag / strict b/a / weak-association cuts."
        )
    return hosts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pipeline-csv", type=Path, default=DEFAULT_PIPELINE)
    parser.add_argument("--n-draws", type=int, default=N_DRAWS)
    parser.add_argument("--q0", type=float, default=Q0)
    parser.add_argument(
        "--inflate",
        type=float,
        nargs="+",
        default=list(BA_ERR_INFLATE_DEFAULTS),
        help="One or more b/a error inflation factors (default: 3 and 10).",
    )
    parser.add_argument("--floor-deg", type=float, default=INC_FLOOR_DEG)
    parser.add_argument("--mag-limit", type=float, default=21.0)
    parser.add_argument(
        "--mag-column",
        default="mag_final",
        help="Magnitude column for the mr cut (default: mag_final).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=PLOTS_NULL / "v2" / "frb_inflated_error_cdf",
    )
    args = parser.parse_args()

    hosts = select_hosts(
        args.pipeline_csv,
        q0=args.q0,
        mag_limit=args.mag_limit,
        mag_column=args.mag_column,
    )
    ba = pd.to_numeric(hosts["b_a"], errors="coerce").to_numpy(dtype=float)
    ba_err = pd.to_numeric(hosts["b_a_err"], errors="coerce").to_numpy(dtype=float)
    n = len(hosts)

    # Point-estimate median cos(i) for legend (Hubble on measured b/a).
    cosi0 = np.array(
        [
            max(0.0, min(1.0, float(np.sqrt(max(0.0, (q**2 - args.q0**2) / (1.0 - args.q0**2))))))
            for q in ba
        ]
    )
    med_cosi = float(np.median(cosi0))

    n_zero_err = int(np.sum(~np.isfinite(ba_err) | (ba_err <= 0)))
    print(
        f"[*] N={n} hosts ({args.mag_column}<={args.mag_limit:g}, "
        f"strict b/a>{args.q0:g}, weak excluded)  "
        f"med_cosi(point)={med_cosi:.3f}  MC={args.n_draws}  "
        f"floor={args.floor_deg:g} deg",
        flush=True,
    )
    if n_zero_err:
        print(
            f"[*] {n_zero_err}/{n} hosts have b_a_err<=0; only {args.floor_deg:g} deg floor applies",
            flush=True,
        )

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    mag_cols = [c for c in ("mag_final", "mag") if c in hosts.columns]
    summary = hosts[["frb", "b_a", "b_a_err", "inc", "inc_err", *mag_cols]].copy()
    summary["cosi_point"] = cosi0
    summary_path = out_dir / "frb_inflated_error_cdf_hosts.csv"
    summary.to_csv(summary_path, index=False)
    print(f"[*] Wrote {summary_path}", flush=True)

    font_prop = default_font()
    col = "#d62728"

    for inflate in args.inflate:
        tag = f"x{int(inflate) if float(inflate).is_integer() else inflate:g}"
        print(f"[*] Running inflate={inflate:g}x -> *_{tag}.*", flush=True)

        med, lo, hi = mc_inflated_cdf_envelope(
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
        ax.plot((0, 1), (0, 1), "k--", lw=1.2, zorder=1)
        ax.fill_between(X, lo, hi, color=col, alpha=0.22, linewidth=0, zorder=2, label="68% CI")
        ax.plot(
            X,
            med,
            color=col,
            lw=2.0,
            zorder=3,
            label=rf"FRB hosts  ($N={n}$, med$\approx${med_cosi:.3f})",
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
            rf"FRB GALFIT $\cos(i)$ CDF  "
            rf"($m_r\leq{args.mag_limit:g}$, $N={n}\times{args.n_draws // 1000}$k, "
            rf"${inflate:g}\times\sigma_{{b/a}}+{args.floor_deg:g}^\circ$ floor, 68% CI)",
            fontsize=11,
        )
        ax.legend(loc="upper left", fontsize=9, frameon=True)
        ax.grid(True, alpha=0.3)
        add_inclination_top_axis(ax, font_prop)
        fig.tight_layout()

        stem = out_dir / f"frb_inflated_error_cdf_{tag}"
        save_figure(fig, stem)
        plt.close(fig)


if __name__ == "__main__":
    main()