"""
Regenerate Jimin cdf_overlay with MC subsample envelopes (advisor-style).

Bug fix vs first attempt: evaluating the *step* ECDF makes CDF values only
{0, 1/23, ..., 1}, so the median across draws is also quantized → staircase.
Advisor-style (and our compare_sdss_legacy_null_distributions.py) linearly
interpolates each draw's ECDF onto a continuous grid first.

Default: N=23 galaxies × 10_000 draws; median + 16–84% (68% CI) band.

Run from repo root::

    python scripts/plot_jimin_cdf_overlay_mc.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import Q0, cosi_array_from_df  # noqa: E402
from pipeline_null_plot_utils import REPO_ROOT  # noqa: E402

JIMIN = REPO_ROOT / "plots" / "plots_null" / "v2" / "sdss_audit" / "Jimin"
N_SAMPLE = 23
N_DRAWS = 10_000
SEED = 42
X = np.linspace(0.0, 1.0, 200)

VARIANTS = [
    {
        "csv": JIMIN / "catalog" / "v1_fracDev0_and_lnL.csv",
        "short": "fracDev0+lnL",
        "color": "#e41a1c",
    },
    {
        "csv": JIMIN / "catalog" / "v2_lnL_exp.csv",
        "short": "lnL only",
        "color": "#377eb8",
    },
]


def strict_cosi(df: pd.DataFrame) -> np.ndarray:
    ba = pd.to_numeric(df["expAB_r"], errors="coerce")
    ok = ba.notna() & (ba > Q0) & (ba <= 1.0)
    return cosi_array_from_df(df.loc[ok], q_col="expAB_r", q0=Q0)


def mc_envelope(
    cosi: np.ndarray,
    *,
    n_sample: int,
    n_draws: int,
    x: np.ndarray,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Advisor / repo style: for each draw, build a continuous ECDF via linear
    interpolation through the order statistics (plus endpoints at 0 and 1),
    evaluate on ``x``, then take median + 16/84 percentiles across draws.
    """
    ref = np.asarray(cosi, dtype=float)
    ref = ref[np.isfinite(ref)]
    if len(ref) < n_sample:
        raise RuntimeError(f"pool N={len(ref)} < n_sample={n_sample}")

    rng = np.random.default_rng(seed)
    # y-grid for one draw's ECDF corners: F=0 at cos=0, then k/n at order stats, F=1 at cos=1
    y_corners = np.concatenate([[0.0], np.arange(1, n_sample + 1) / n_sample, [1.0]])
    ecdfs = np.empty((n_draws, len(x)), dtype=float)

    for i in range(n_draws):
        draw = np.sort(rng.choice(ref, size=n_sample, replace=False))
        x_corners = np.concatenate([[0.0], draw, [1.0]])
        # ensure strictly increasing x for interp (rare ties)
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


def main() -> None:
    fig, ax = plt.subplots(figsize=(7.2, 7.2))
    ax.plot((0, 1), (0, 1), "k--", lw=1.2, label="Uniform")

    for i, meta in enumerate(VARIANTS):
        df = pd.read_csv(meta["csv"])
        if "expAB_r" not in df.columns and "expab_r" in df.columns:
            df = df.rename(columns={"expab_r": "expAB_r"})
        cosi = strict_cosi(df)
        pool_med = float(np.median(cosi))
        print(
            f"[*] {meta['short']}: pool N={len(cosi):,}  med_cosi={pool_med:.3f}  "
            f"MC n={N_SAMPLE} x {N_DRAWS}",
            flush=True,
        )
        med, lo, hi = mc_envelope(
            cosi,
            n_sample=N_SAMPLE,
            n_draws=N_DRAWS,
            x=X,
            seed=SEED + i,
        )
        col = meta["color"]
        ax.fill_between(X, lo, hi, color=col, alpha=0.22, linewidth=0, zorder=2)
        ax.plot(
            X,
            med,
            color=col,
            lw=2.0,
            zorder=3,
            label=f"{meta['short']}  (med={pool_med:.3f})",
        )

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(rf"$\cos(i)$ (Hubble, $q_0={Q0:g}$)")
    ax.set_ylabel("Cumulative distribution")
    ax.set_title(rf"Jimin MC CDF  ($N={N_SAMPLE}\times{N_DRAWS // 1000}$k, 68% CI)", fontsize=11)
    ax.legend(loc="upper left", fontsize=8, frameon=True)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = JIMIN / "plots" / "cdf_overlay.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(JIMIN / "plots" / "cdf_overlay_mc.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[*] Wrote {out}", flush=True)


if __name__ == "__main__":
    main()
