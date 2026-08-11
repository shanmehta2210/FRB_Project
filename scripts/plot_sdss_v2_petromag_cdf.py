"""
SDSS v2 strict null cos(i) CDF with petroMag_r <= 21.

Production cuts (same as modelMag mag-cut CDFs):
  u-r < 2.3, lnLExp_r > lnLDeV_r, expAB_r > q0 (=0.2), then mag cut.
cos(i) via Hubble formula with q0=0.2.

Also overlays the usual modelMag_r <= 21 pool for comparison.

Outputs under plots/plots_null/v2/sdss_audit/formal/:
  - cdf_petromag21_strict.png
  - cdf_petromag21_vs_modelmag21_overlay.png
  - cdf_petromag21_strict_summary.csv

Run from repo root::

    python scripts/plot_sdss_v2_petromag_cdf.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (  # noqa: E402
    Q0,
    SDSS_UR_MAX_CDF,
    apply_strict_q_cut,
    cosi_array_from_df,
    ensure_sdss_colors,
    filter_sdss_drop_dev_winners,
    filter_sdss_ur,
    slice_null_base_by_mag,
)
from pipeline_null_plot_utils import DEFAULT_SDSS_V2, REPO_ROOT  # noqa: E402

OUT_DIR = REPO_ROOT / "plots" / "plots_null" / "v2" / "sdss_audit" / "formal"
MAG_LIMIT = 21.0
USECOLS = (
    "modelMag_r",
    "petroMag_r",
    "modelMag_u",
    "modelMag_g",
    "expAB_r",
    "lnLDeV_r",
    "lnLExp_r",
    "model_winner_is_exp",
    "u_r",
    "g_r",
)


def empirical_cdf(vals: np.ndarray, x: np.ndarray) -> np.ndarray:
    s = np.sort(vals[np.isfinite(vals)])
    return np.searchsorted(s, x, side="right") / max(1, len(s))


def production_base(df: pd.DataFrame) -> pd.DataFrame:
    """u-r + lnL exp-wins + expAB_r > q0; keep both mag columns."""
    out = ensure_sdss_colors(df)
    out = filter_sdss_ur(out, SDSS_UR_MAX_CDF)
    out = filter_sdss_drop_dev_winners(out)
    out = apply_strict_q_cut(out, q_col="expAB_r", q0=Q0)
    keep = [c for c in ("modelMag_r", "petroMag_r", "expAB_r") if c in out.columns]
    return out.loc[:, keep].reset_index(drop=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sdss-csv", type=Path, default=DEFAULT_SDSS_V2)
    p.add_argument("--out-dir", type=Path, default=OUT_DIR)
    p.add_argument("--mag-limit", type=float, default=MAG_LIMIT)
    args = p.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[*] Loading {args.sdss_csv.name} ...", flush=True)
    df = pd.read_csv(args.sdss_csv, usecols=lambda c: c in USECOLS)
    base = production_base(df)
    print(
        f"    production base (ur+lnL+strict ba): N={len(base):,}",
        flush=True,
    )

    petro = slice_null_base_by_mag(
        base, mag_column="petroMag_r", mag_limit=args.mag_limit
    )
    model = slice_null_base_by_mag(
        base, mag_column="modelMag_r", mag_limit=args.mag_limit
    )
    cosi_p = cosi_array_from_df(petro, q_col="expAB_r", q0=Q0)
    cosi_m = cosi_array_from_df(model, q_col="expAB_r", q0=Q0)
    med_p = float(np.median(cosi_p)) if len(cosi_p) else float("nan")
    med_m = float(np.median(cosi_m)) if len(cosi_m) else float("nan")
    print(
        f"    petroMag_r <= {args.mag_limit:g}: N={len(cosi_p):,}  med cos(i)={med_p:.3f}",
        flush=True,
    )
    print(
        f"    modelMag_r <= {args.mag_limit:g}: N={len(cosi_m):,}  med cos(i)={med_m:.3f}",
        flush=True,
    )

    x = np.linspace(0, 1, 401)

    # Standalone petro CDF
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot((0, 1), (0, 1), "k--", lw=1.2, label="Uniform")
    ax.plot(
        x,
        empirical_cdf(cosi_p, x),
        color="#4daf4a",
        lw=2.0,
        label=f"SDSS strict (petroMag_r ≤ {args.mag_limit:g})",
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"$\cos(i)$ (Hubble, $q_0=0.2$)")
    ax.set_ylabel("Cumulative distribution")
    ax.set_title(
        f"SDSS v2 strict null — petroMag$_r$ ≤ {args.mag_limit:g}\n"
        f"u−r < {SDSS_UR_MAX_CDF:g}, lnL exp-wins, expAB$_r$ > {Q0:g}\n"
        f"N = {len(cosi_p):,}  |  median cos(i) = {med_p:.3f}"
    )
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out_petro = args.out_dir / "cdf_petromag21_strict.png"
    fig.savefig(out_petro, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Overlay petro vs model
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot((0, 1), (0, 1), "k--", lw=1.2, label="Uniform")
    ax.plot(
        x,
        empirical_cdf(cosi_p, x),
        color="#e41a1c",
        lw=2.2,
        label=f"petroMag_r ≤ {args.mag_limit:g}  (N={len(cosi_p):,}, med={med_p:.3f})",
    )
    ax.plot(
        x,
        empirical_cdf(cosi_m, x),
        color="#377eb8",
        lw=2.0,
        ls="--",
        label=f"modelMag_r ≤ {args.mag_limit:g}  (N={len(cosi_m):,}, med={med_m:.3f})",
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"$\cos(i)$ (Hubble, $q_0=0.2$)")
    ax.set_ylabel("Cumulative distribution")
    ax.set_title(
        f"SDSS v2 strict null CDFs — petroMag vs modelMag cut at {args.mag_limit:g}\n"
        f"u−r < {SDSS_UR_MAX_CDF:g}, lnL exp-wins, expAB$_r$ > {Q0:g}"
    )
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out_ov = args.out_dir / "cdf_petromag21_vs_modelmag21_overlay.png"
    fig.savefig(out_ov, dpi=300, bbox_inches="tight")
    plt.close(fig)

    pd.DataFrame(
        [
            {
                "mag_column": "petroMag_r",
                "mag_limit": args.mag_limit,
                "n": len(cosi_p),
                "median_cosi": round(med_p, 4),
                "ur_max": SDSS_UR_MAX_CDF,
                "q0": Q0,
                "lnL_exp_wins": True,
            },
            {
                "mag_column": "modelMag_r",
                "mag_limit": args.mag_limit,
                "n": len(cosi_m),
                "median_cosi": round(med_m, 4),
                "ur_max": SDSS_UR_MAX_CDF,
                "q0": Q0,
                "lnL_exp_wins": True,
            },
        ]
    ).to_csv(args.out_dir / "cdf_petromag21_strict_summary.csv", index=False)

    print(f"[*] Wrote {out_petro.name}", flush=True)
    print(f"[*] Wrote {out_ov.name}", flush=True)


if __name__ == "__main__":
    main()
