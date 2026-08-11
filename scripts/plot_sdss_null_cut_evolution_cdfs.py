#!/usr/bin/env python3
"""
SDSS null cos(i) diagnostics: cut evolution, final mag slices, MC null sensitivity.

Strict b/a > q0 is applied first (Hubble-valid pool only).

Outputs under plots/plots_null/v1_null_cdf_inclination/diagnostics/sdss_cut_evolution/
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
    apply_mag_cut,
    apply_strict_q_cut,
    cosi_array_from_df,
    filter_sdss_drop_dev_winners,
    filter_frb_hosts_mag,
    filter_frb_hosts_strict_ba,
    filter_sdss_ur,
    prepare_null_strict_color_base,
    read_sdss_null_catalog,
    slice_null_base_by_mag,
)
from pipeline_null_plot_utils import (  # noqa: E402
    DEFAULT_PIPELINE,
    cdf_envelope,
    load_pipeline_hosts,
)

REPO = _SCRIPTS.parent
DEFAULT_OUT = (
    REPO
    / "plots"
    / "plots_null"
    / "v1_null_cdf_inclination"
    / "diagnostics"
    / "sdss_cut_evolution"
)

# Production overlay defaults (plot_null_mag_cut_cdfs.py)
MC_DRAWS_DEFAULT = 10_000
MC_SAMPLE_FRB_DEFAULT = 41


def empirical_cdf(values: np.ndarray, x: np.ndarray) -> np.ndarray:
    v = np.sort(values[np.isfinite(values)])
    if len(v) == 0:
        return np.zeros_like(x)
    return np.searchsorted(v, x, side="right") / len(v)


def x_at_cdf(y: np.ndarray, x: np.ndarray, target: float = 0.5) -> float:
    idx = np.searchsorted(y, target)
    return float(x[min(idx, len(x) - 1)])


def build_pre_mag_stages(
    df: pd.DataFrame,
    *,
    ur_max: float,
    q0: float,
) -> list[tuple[str, str, pd.DataFrame]]:
    q_tag = f"{q0:g}".replace(".", "p")
    s = apply_strict_q_cut(df, q_col="expAB_r", q0=q0)
    stages: list[tuple[str, str, pd.DataFrame]] = [
        (
            f"00_strict_ba_gt_{q_tag}",
            f"Strict pool: expAB_r > {q0:g} (Hubble-valid)",
            s,
        ),
    ]
    s = filter_sdss_ur(s, ur_max)
    stages.append((f"01_ur_lt_{ur_max:g}".replace(".", "p"), f"+ u-r < {ur_max:g}", s))
    s = filter_sdss_drop_dev_winners(s)
    stages.append(("02_lnl_exp_wins", "+ lnLExp_r > lnLDeV_r", s))
    return stages


def production_pool_at_mag(
    df: pd.DataFrame,
    *,
    ur_max: float,
    mag_limit: float,
) -> pd.DataFrame:
    """Same pool as plot_null_mag_cut_cdfs (ur → lnL → strict → mag)."""
    base = prepare_null_strict_color_base(
        df,
        mag_column="modelMag_r",
        q_column="expAB_r",
        is_legacy=False,
        sdss_ur_max=ur_max,
        sdss_exp_winner_only=True,
    )
    return slice_null_base_by_mag(base, mag_column="modelMag_r", mag_limit=mag_limit)


def mag_tag(mag_limit: float) -> str:
    return str(int(mag_limit)) if mag_limit == int(mag_limit) else f"{mag_limit:g}".replace(".", "p")


def plot_single_cdf(
    x: np.ndarray,
    y: np.ndarray,
    *,
    title: str,
    n: int,
    median_cosi: float,
    out_path: Path,
    ylabel: str = "Cumulative distribution",
    show_uniform: bool = False,
) -> None:
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot(x, y, color="#4daf4a", linewidth=2.0, label="SDSS null")
    if show_uniform:
        ax.plot((0, 1), (0, 1), "k--", linewidth=1.2, label="Uniform")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("cos(i)")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title}\nN = {n:,}  |  median cos(i) = {median_cosi:.3f}")
    if show_uniform:
        ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.25)
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_final_mag_panel(
    panels: list[tuple[float, np.ndarray, np.ndarray, int, float]],
    *,
    out_path: Path,
    catalog_label: str = "SDSS v2",
) -> None:
    """Side-by-side empirical null CDFs with uniform isotropic reference (dashed)."""
    n_panels = len(panels)
    fig, axes = plt.subplots(1, n_panels, figsize=(5.5 * n_panels, 5.5), squeeze=False)
    for ax, (mag_limit, x, y, n_pool, med_cosi) in zip(axes.ravel(), panels):
        ax.plot(x, y, color="#4daf4a", linewidth=2.0, label="SDSS null")
        ax.plot((0, 1), (0, 1), "k--", linewidth=1.2, label="Uniform")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel("cos(i)")
        ax.set_ylabel("Cumulative distribution")
        ax.set_title(
            f"modelMag$_r$ $\\leq$ {mag_limit:g}\n"
            f"N = {n_pool:,}  |  median cos(i) = {med_cosi:.3f}"
        )
        ax.legend(fontsize=7, loc="upper left")
        ax.grid(True, alpha=0.25)
    fig.suptitle(
        f"{catalog_label} strict null ECDF (expAB$_r$, lnL exp-wins, u-r cut)",
        fontsize=12,
        y=1.02,
    )
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_path}")


def plot_mc_sensitivity(
    *,
    cosi_pool: np.ndarray,
    mag_limit: float,
    n_frb: int,
    mc_draws: int,
    mc_sample_sizes: list[int],
    out_path: Path,
) -> pd.DataFrame:
    """Full-pool ECDF vs cdf_envelope at several n_sample (each: n_sample galaxies, mc_draws times)."""
    x = np.linspace(0, 1, 201)
    y_full = empirical_cdf(cosi_pool, x)
    n_pool = len(cosi_pool)

    fig, ax = plt.subplots(figsize=(9, 8))
    ax.plot(
        x,
        y_full,
        color="#4daf4a",
        linewidth=2.4,
        label=f"Full pool ECDF (N={n_pool:,})",
        zorder=5,
    )

    rows = []
    colors = plt.cm.viridis(np.linspace(0.15, 0.85, len(mc_sample_sizes)))
    for color, n_s in zip(colors, mc_sample_sizes):
        ns = min(int(n_s), n_pool)
        x_mc, y_mc, lo, hi = cdf_envelope(cosi_pool, n_sample=ns, n_draws=mc_draws)
        j05 = int(np.argmin(np.abs(x_mc - 0.5)))
        i50 = int(np.argmin(np.abs(y_mc - 0.5)))
        rows.append(
            {
                "mag_limit": mag_limit,
                "n_pool": n_pool,
                "n_sample": ns,
                "n_draws": mc_draws,
                "cdf_at_cosi_0.5": float(y_mc[j05]),
                "x_at_cdf_0.5": float(x_mc[i50]),
                "max_abs_diff_vs_full_ecdf": float(
                    np.max(np.abs(y_full - np.interp(x, x_mc, y_mc)))
                ),
            }
        )
        lbl = (
            f"MC mean: {ns} gal × {mc_draws:,} draws"
            if ns < n_pool
            else f"MC mean: full N (={ns}) × {mc_draws:,} draws"
        )
        ax.plot(x_mc, y_mc, color=color, linewidth=1.8, linestyle="--", label=lbl)
        if ns == n_frb:
            ax.fill_between(x_mc, lo, hi, color=color, alpha=0.15)

    ax.plot((0, 1), (0, 1), "k:", linewidth=1.0, alpha=0.5, label="Uniform")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("cos(i)")
    ax.set_ylabel("Cumulative distribution")
    ax.set_title(
        f"SDSS strict pool @ modelMag_r <= {mag_limit:g}\n"
        f"Read carefully: y @ x=0.5 ≠ x where CDF=0.5  |  production uses n={n_frb}, "
        f"{mc_draws:,} draws"
    )
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(True, alpha=0.25)
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_path}")
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sdss-csv", default=str(REPO / "catalog/SDSS_catalog_v1_allsky_modelmr.csv"))
    parser.add_argument(
        "--final-mag-limits",
        type=float,
        nargs="+",
        default=[20.0, 21.0, 22.0],
        help="Final-stage empirical CDFs (production pools at each mag cut).",
    )
    parser.add_argument(
        "--evolution-mag-limit",
        type=float,
        default=21.0,
        help="Mag limit for 00–02 + 03 stage-evolution series.",
    )
    parser.add_argument("--ur-max", type=float, default=SDSS_UR_MAX_CDF)
    parser.add_argument("--q0", type=float, default=Q0)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--pipeline-csv", default=str(DEFAULT_PIPELINE))
    parser.add_argument("--mc-draws", type=int, default=MC_DRAWS_DEFAULT)
    parser.add_argument(
        "--mc-sample-sizes",
        type=int,
        nargs="+",
        default=[41, 100, 500, 2000, 5000],
        help="Subsample sizes for MC sensitivity (41 = production N_FRB).",
    )
    parser.add_argument("--no-panel", action="store_true")
    parser.add_argument(
        "--no-mc",
        action="store_true",
        help="Skip MC sensitivity plots (faster; still writes final CDF panel).",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    final_dir = out_dir / "final_by_mag"
    mc_dir = out_dir / "mc_sensitivity"
    final_dir.mkdir(parents=True, exist_ok=True)
    mc_dir.mkdir(parents=True, exist_ok=True)

    df = read_sdss_null_catalog(args.sdss_csv)
    hosts = load_pipeline_hosts(args.pipeline_csv)

    # --- Stage evolution (fixed mag for step 03) ---
    pre_mag = build_pre_mag_stages(df, ur_max=args.ur_max, q0=args.q0)
    evo_stages = list(pre_mag)
    s_mag = apply_mag_cut(pre_mag[-1][2], mag_column="modelMag_r", limit=args.evolution_mag_limit)
    evo_stages.append(
        (
            f"03_mag_lt_{mag_tag(args.evolution_mag_limit)}",
            f"+ modelMag_r <= {args.evolution_mag_limit:g}",
            s_mag,
        )
    )

    x = np.linspace(0, 1, 201)
    evo_rows = []
    cdf_out = pd.DataFrame({"cos_i": x})
    for tag, label, sub in evo_stages:
        cosi = cosi_array_from_df(sub, q_col="expAB_r", q0=args.q0)
        y = empirical_cdf(cosi, x)
        cdf_out[tag] = y
        med = float(np.median(cosi))
        evo_rows.append(
            {
                "stage_tag": tag,
                "stage_label": label,
                "n": len(sub),
                "median_cosi": med,
                "mean_cosi": float(np.mean(cosi)),
                "cdf_at_cosi_0.5_y": float(empirical_cdf(cosi, np.array([0.5]))[0]),
                "x_at_cdf_0.5": x_at_cdf(y, x, 0.5),
                "median_expAB_r": float(pd.to_numeric(sub["expAB_r"], errors="coerce").median()),
            }
        )
        plot_single_cdf(
            x,
            y,
            title=label,
            n=len(sub),
            median_cosi=med,
            out_path=out_dir / f"cdf_{tag}.png",
        )

    pd.DataFrame(evo_rows).to_csv(out_dir / "stage_summary.csv", index=False)
    cdf_out.to_csv(out_dir / "cdf_curves.csv", index=False)

    if not args.no_panel:
        fig, axes = plt.subplots(2, 2, figsize=(10, 10))
        for ax, (tag, label, sub) in zip(axes.ravel(), evo_stages):
            cosi = cosi_array_from_df(sub, q_col="expAB_r", q0=args.q0)
            ax.plot(x, empirical_cdf(cosi, x), color="#4daf4a", lw=1.8)
            ax.set_title(f"{label}\nN={len(sub):,}", fontsize=9)
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.grid(True, alpha=0.2)
        fig.suptitle(f"Cut evolution (final step m<{args.evolution_mag_limit:g})", fontsize=11)
        plt.tight_layout()
        p = out_dir / "cdf_all_stages_panel.png"
        fig.savefig(p, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved {p}")

    # --- Final empirical CDFs at mag 20 / 21 / 22 ---
    final_rows = []
    mc_audit_rows = []
    panel_data: list[tuple[float, np.ndarray, np.ndarray, int, float]] = []
    for mag_limit in sorted(args.final_mag_limits):
        pool = production_pool_at_mag(df, ur_max=args.ur_max, mag_limit=mag_limit)
        cosi = cosi_array_from_df(pool, q_col="expAB_r", q0=args.q0)
        y = empirical_cdf(cosi, x)
        med = float(np.median(cosi))
        mt = mag_tag(mag_limit)
        final_rows.append(
            {
                "mag_limit": mag_limit,
                "n": len(pool),
                "median_cosi": med,
                "cdf_at_cosi_0.5_y": float(empirical_cdf(cosi, np.array([0.5]))[0]),
                "x_at_cdf_0.5": x_at_cdf(y, x, 0.5),
            }
        )
        plot_single_cdf(
            x,
            y,
            title=f"Production pool: modelMag_r <= {mag_limit:g} (strict+color)",
            n=len(pool),
            median_cosi=med,
            out_path=final_dir / f"cdf_final_empirical_mag{mt}.png",
            show_uniform=True,
        )
        panel_data.append((mag_limit, x, y, len(pool), med))

        if not args.no_mc:
            frb = filter_frb_hosts_strict_ba(
                filter_frb_hosts_mag(hosts, mag_limit=mag_limit, mag_column="mag")
            )
            n_frb = len(frb)
            sizes = sorted(set(args.mc_sample_sizes + [n_frb]))

            mc_df = plot_mc_sensitivity(
                cosi_pool=cosi,
                mag_limit=mag_limit,
                n_frb=n_frb,
                mc_draws=args.mc_draws,
                mc_sample_sizes=sizes,
                out_path=mc_dir / f"mc_sensitivity_mag{mt}.png",
            )
            mc_audit_rows.append(mc_df)

            # Production-exact MC audit row
            x_mc, y_mc, _, _ = cdf_envelope(cosi, n_sample=n_frb, n_draws=args.mc_draws)
            j05 = int(np.argmin(np.abs(x_mc - 0.5)))
            i50 = int(np.argmin(np.abs(y_mc - 0.5)))
            mc_audit_rows.append(
                pd.DataFrame(
                    [
                        {
                            "mag_limit": mag_limit,
                            "n_pool": len(pool),
                            "n_sample": n_frb,
                            "n_draws": args.mc_draws,
                            "cdf_at_cosi_0.5": float(y_mc[j05]),
                            "x_at_cdf_0.5": float(x_mc[i50]),
                            "note": "production_overlay_exact",
                        }
                    ]
                )
            )

    pd.DataFrame(final_rows).to_csv(final_dir / "final_by_mag_summary.csv", index=False)
    if mc_audit_rows:
        pd.concat(mc_audit_rows, ignore_index=True).to_csv(mc_dir / "mc_audit.csv", index=False)

    if len(panel_data) >= 2:
        mags = "_".join(mag_tag(m) for m, *_ in panel_data)
        plot_final_mag_panel(
            panel_data,
            out_path=final_dir / f"cdf_final_empirical_panel_mag{mags}.png",
            catalog_label="SDSS null",
        )

    # README snippet as MC_AUDIT.txt
    audit_txt = out_dir / "MC_AUDIT.txt"
    audit_txt.write_text(
        f"Production null curve: cdf_envelope(pool, n_sample=N_FRB, n_draws={args.mc_draws})\n"
        f"  Each draw: WITHOUT replacement, sample n_sample cos(i) values from pool, "
        f"build step-function CDF, evaluate on 100-point grid; mean over {args.mc_draws} draws.\n"
        f"  Default N_FRB = {MC_SAMPLE_FRB_DEFAULT} (actual per mag in mc_audit.csv).\n\n"
        "Two ways to read the plot (do not confuse):\n"
        "  cdf_at_cosi_0.5  = y-value at x=0.5  (~0.38 for mag21: fraction with cos i < 0.5)\n"
        "  x_at_cdf_0.5     = x-value where CDF crosses 0.5 (~0.58: median cos i)\n\n"
        "Subsampling 41 from 27k does NOT move cdf_at_cosi_0.5 to ~0.65; "
        "that would be mis-reading x_at_cdf_0.5 as the first quantity.\n",
        encoding="utf-8",
    )
    print(f"Wrote {audit_txt}")
    print(pd.DataFrame(final_rows).to_string(index=False))
    print(f"\nWrote {out_dir}")


if __name__ == "__main__":
    main()
