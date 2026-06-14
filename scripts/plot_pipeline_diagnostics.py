"""
Pipeline diagnostic plots: inclination histogram, Re/n CDFs, random-host test, master deltas.

Default null sample: **strict** at $m_r<21$ with SDSS $u-r<2.3$ and Legacy $g-r<0.75$.
Sky maps and Re/n histograms are not generated.

Run from repo root:
    python scripts/plot_pipeline_diagnostics.py --section all
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (  # noqa: E402
    LEGACY_GR_MAX_CDF,
    SDSS_UR_MAX_CDF,
    cosi_array_from_df,
    n_from_legacy_df,
    n_from_sdss_df,
    re_arcsec_from_legacy_df,
    re_arcsec_from_sdss_df,
)
from pipeline_null_plot_utils import (  # noqa: E402
    DEFAULT_LEGACY,
    DEFAULT_PIPELINE,
    DEFAULT_SDSS,
    PLOTS_NULL,
    build_frb_mc_draws_inc,
    frb_hosts_for_cdf,
    frb_n,
    frb_re_arcsec,
    load_pipeline_hosts,
    load_null_cuts,
    mc_mean_cdf_from_draws,
    plot_delta_cleveland,
    plot_inclination_bin_fractions,
    plot_inclination_cdf_overlay,
    plot_inclination_cdf_dual_null_overlay,
    plot_random_host_inclination,
    plot_scalar_cdf_overlay,
    cdf_envelope,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DIFF = REPO_ROOT / "pipeline_vs_master_galfit_diff.csv"
DEFAULT_TAG = "v1_allsky_modelmr_strict_color"


def section_hist(
    hosts: pd.DataFrame,
    sdss_cut: pd.DataFrame,
    *,
    tag: str,
    sample_mode: str,
    sdss_q_column: str = "expAB_r",
    q0: float = 0.2,
    mag_limit: float = 21.0,
) -> None:
    frb_plot = frb_hosts_for_cdf(
        hosts,
        sample_mode=sample_mode,
        q0=q0,
        mag_limit=mag_limit,
    )
    inc_frb = pd.to_numeric(frb_plot["inc"], errors="coerce").to_numpy()
    sdss_inc = np.degrees(
        np.arccos(np.clip(cosi_array_from_df(sdss_cut, q_col=sdss_q_column, q0=q0), 0, 1))
    )
    plot_inclination_bin_fractions(
        frb_inc=inc_frb,
        sdss_inc=sdss_inc,
        n_frb=len(frb_plot),
        n_sdss=len(sdss_cut),
        out_stem=PLOTS_NULL / "v1_hist_inclination" / f"hist_inclination_frb_vs_sdss_{tag}",
        sample_mode=sample_mode,
    )


def section_cdf(
    hosts: pd.DataFrame,
    legacy_cut: pd.DataFrame,
    sdss_cut: pd.DataFrame,
    *,
    tag: str,
    mag_limit: float,
    mc_draws_null: int,
    mc_alpha: float,
    sample_mode: str,
    sdss_q_column: str = "expAB_r",
    q0: float = 0.2,
    sdss_ur_max: float = SDSS_UR_MAX_CDF,
    legacy_gr_max: float = LEGACY_GR_MAX_CDF,
) -> None:
    frb_plot = frb_hosts_for_cdf(
        hosts,
        sample_mode=sample_mode,
        q0=q0,
        mag_limit=mag_limit,
    )
    n_frb = len(frb_plot)
    frb_draws = build_frb_mc_draws_inc(frb_plot)
    x_frb, y_frb = mc_mean_cdf_from_draws(frb_draws)
    legacy_color_lbl = f"strict, $g-r<{legacy_gr_max:g}$"
    sdss_color_lbl = f"strict, $u-r<{sdss_ur_max:g}$"
    mag_tag = int(mag_limit)
    out_base = PLOTS_NULL / "v1_null_cdf_inclination" / "mag_cuts" / f"mag{mag_tag}"
    frb_title = f"FRB hosts (N={n_frb}, $m_r<{mag_limit:g}$"
    if sample_mode == "strict":
        frb_title += f", $b/a>{q0:g}$"
    frb_title += ")"

    for label, cut, color, stem, q_col, shape_lbl, color_lbl in (
        (
            "Legacy",
            legacy_cut,
            "#377eb8",
            "legacy",
            "expAB_r",
            "Tractor $b/a$",
            legacy_color_lbl,
        ),
        (
            "SDSS",
            sdss_cut,
            "#4daf4a",
            "sdss",
            sdss_q_column,
            "best-model $b/a$",
            sdss_color_lbl,
        ),
    ):
        cosi = cosi_array_from_df(cut, q_col=q_col, q0=q0)
        x_n, mn, lo, hi = cdf_envelope(cosi, n_sample=n_frb, n_draws=mc_draws_null)
        out_dir = out_base / f"{stem}_{sample_mode}"
        out_dir.mkdir(parents=True, exist_ok=True)
        plot_inclination_cdf_overlay(
            null_label=f"{label} null ({color_lbl}, N={len(cut)})",
            null_color=color,
            x_null=x_n,
            mean_null=mn,
            lo_null=lo,
            hi_null=hi,
            frb_draws=frb_draws,
            x_frb=x_frb,
            y_frb=y_frb,
            n_frb=n_frb,
            title=(
                f"{frb_title} vs {label} null "
                f"($m_r<{mag_limit:g}$, {shape_lbl}, {color_lbl})"
            ),
            out_stem=out_dir / "null_cdf_inclination",
            mc_alpha=mc_alpha,
        )

    combined_dir = out_base / f"legacy_sdss_{sample_mode}_combined"
    combined_dir.mkdir(parents=True, exist_ok=True)
    null_layers = []
    for label, cut, color, _stem, q_col, color_lbl in (
        ("Legacy null", legacy_cut, "#377eb8", "legacy", "expAB_r", legacy_color_lbl),
        ("SDSS null", sdss_cut, "#4daf4a", "sdss", sdss_q_column, sdss_color_lbl),
    ):
        cosi = cosi_array_from_df(cut, q_col=q_col, q0=q0)
        x_n, mn, lo, hi = cdf_envelope(cosi, n_sample=n_frb, n_draws=mc_draws_null)
        null_layers.append(
            (f"{label} ({color_lbl})", color, x_n, mn, lo, hi)
        )
    plot_inclination_cdf_dual_null_overlay(
        nulls=tuple(null_layers),
        frb_draws=frb_draws,
        x_frb=x_frb,
        y_frb=y_frb,
        n_frb=n_frb,
        title=(
            f"{frb_title} vs Legacy + SDSS nulls "
            f"($m_r<{mag_limit:g}$, Legacy {legacy_color_lbl}, SDSS {sdss_color_lbl})"
        ),
        out_stem=combined_dir / "null_cdf_inclination",
        mc_alpha=mc_alpha,
    )

    frb_re = frb_re_arcsec(frb_plot)
    frb_nvals = frb_n(frb_plot)
    out_re = PLOTS_NULL / "v1_null_cdf_re"
    out_n = PLOTS_NULL / "v1_null_cdf_n"

    plot_scalar_cdf_overlay(
        null_label=f"Legacy null $R_e$ ({legacy_color_lbl})",
        null_color="#377eb8",
        null_vals=re_arcsec_from_legacy_df(legacy_cut),
        frb_vals=frb_re,
        n_frb=n_frb,
        xlabel=r"$R_e$ (arcsec)",
        title=f"FRB vs Legacy null $R_e$ ({legacy_color_lbl}; per-FRB plate scale)",
        out_stem=out_re / f"null_cdf_pipeline_legacy_re_{tag}",
        n_draws_null=mc_draws_null,
    )
    if "best_model_re_r" in sdss_cut.columns and sdss_cut["best_model_re_r"].notna().any():
        plot_scalar_cdf_overlay(
            null_label=f"SDSS null $R_e$ ({sdss_color_lbl})",
            null_color="#4daf4a",
            null_vals=re_arcsec_from_sdss_df(sdss_cut),
            frb_vals=frb_re,
            n_frb=n_frb,
            xlabel=r"$R_e$ (arcsec)",
            title=f"FRB vs SDSS null $R_e$ ({sdss_color_lbl})",
            out_stem=out_re / f"null_cdf_pipeline_sdss_re_{tag}",
            n_draws_null=mc_draws_null,
        )

    plot_scalar_cdf_overlay(
        null_label=f"Legacy null Sérsic $n$ ({legacy_color_lbl})",
        null_color="#377eb8",
        null_vals=n_from_legacy_df(legacy_cut),
        frb_vals=frb_nvals,
        n_frb=n_frb,
        xlabel="Sérsic n",
        title=f"FRB vs Legacy null $n$ ({legacy_color_lbl})",
        out_stem=out_n / f"null_cdf_pipeline_legacy_n_{tag}",
        n_draws_null=mc_draws_null,
    )
    plot_scalar_cdf_overlay(
        null_label=f"SDSS $n_{{eff}}$ proxy ({sdss_color_lbl})",
        null_color="#4daf4a",
        null_vals=n_from_sdss_df(sdss_cut),
        frb_vals=frb_nvals,
        n_frb=n_frb,
        xlabel=r"$n_{\rm eff}$ (from fracDeV)",
        title=f"FRB vs SDSS null $n$ ({sdss_color_lbl})",
        out_stem=out_n / f"null_cdf_pipeline_sdss_n_{tag}",
        n_draws_null=mc_draws_null,
    )


def section_random(
    hosts: pd.DataFrame,
    legacy_cut: pd.DataFrame,
    sdss_cut: pd.DataFrame,
    *,
    tag: str,
    mag_limit: float,
    n_random: int,
    mc_draws_null: int,
    sample_mode: str,
    sdss_q_column: str = "expAB_r",
    q0: float = 0.2,
) -> None:
    frb_plot = frb_hosts_for_cdf(
        hosts,
        sample_mode=sample_mode,
        q0=q0,
        mag_limit=mag_limit,
    )
    n_frb = len(frb_plot)
    frb_draws = build_frb_mc_draws_inc(frb_plot)
    x_frb, y_frb = mc_mean_cdf_from_draws(frb_draws)
    out = PLOTS_NULL / "v1_random_host_inclination"
    rows = []
    for label, cut, stem, q_col in (
        ("Legacy", legacy_cut, "legacy", "expAB_r"),
        ("SDSS", sdss_cut, "sdss", sdss_q_column),
    ):
        cosi = cosi_array_from_df(cut, q_col=q_col, q0=q0)
        p = plot_random_host_inclination(
            null_cosi=cosi,
            frb_draws=frb_draws,
            x_frb=x_frb,
            y_frb=y_frb,
            n_frb=n_frb,
            null_label=f"{label} null",
            title=f"Random-host inclination test vs {label} ($m<{mag_limit:g}$, {sample_mode})",
            out_stem=out / f"random_host_inclination_{stem}_{tag}",
            n_random=n_random,
            n_draws_null=mc_draws_null,
        )
        rows.append({"survey": stem, "tag": tag, "p_mean_within_0.02": p})
    pd.DataFrame(rows).to_csv(out / f"random_host_stats_{tag}.csv", index=False)


def section_master(diff_csv: Path) -> None:
    diff = pd.read_csv(diff_csv)
    out = PLOTS_NULL / "v1_pipeline_vs_master"
    for col, label in (
        ("inc_delta", "inclination (deg)"),
        ("re_delta", r"$R_e$ (px)"),
        ("n_delta", "Sérsic n"),
    ):
        plot_delta_cleveland(diff, col, label, out / f"delta_{col.replace('_delta', '')}_pipeline_vs_master")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--section",
        choices=("all", "hist", "cdf", "random", "master"),
        default="all",
    )
    parser.add_argument("--pipeline-csv", type=Path, default=DEFAULT_PIPELINE)
    parser.add_argument("--diff-csv", type=Path, default=DEFAULT_DIFF)
    parser.add_argument("--legacy-csv", type=Path, default=DEFAULT_LEGACY)
    parser.add_argument("--sdss-csv", type=Path, default=DEFAULT_SDSS)
    parser.add_argument(
        "--sample-mode",
        choices=("strict",),
        default="strict",
        help="Null catalog cut (strict b/a > q0 only).",
    )
    parser.add_argument("--sdss-ur-max", type=float, default=SDSS_UR_MAX_CDF)
    parser.add_argument("--legacy-gr-max", type=float, default=LEGACY_GR_MAX_CDF)
    parser.add_argument("--tag", default=DEFAULT_TAG)
    parser.add_argument("--mag-limit", type=float, default=21.0)
    parser.add_argument("--mag-column", default="rmag")
    parser.add_argument("--q0", type=float, default=0.2)
    parser.add_argument("--exclude-types", default="REX,DEV")
    parser.add_argument("--sdss-mag-column", default="modelMag_r")
    parser.add_argument("--sdss-q-column", default="expAB_r")
    parser.add_argument("--mc-draws-null", type=int, default=10000)
    parser.add_argument("--mc-alpha", type=float, default=0.03)
    parser.add_argument("--n-random", type=int, default=200)
    args = parser.parse_args()

    hosts = load_pipeline_hosts(args.pipeline_csv)
    print("[*] Loading null catalogs once (extended usecols for Re/n)...")
    legacy_cut, sdss_cut = load_null_cuts(
        args.legacy_csv,
        args.sdss_csv,
        sample_mode=args.sample_mode,
        mag_limit=args.mag_limit,
        mag_column=args.mag_column,
        q0=args.q0,
        exclude_types=args.exclude_types,
        sdss_q_column=args.sdss_q_column,
        sdss_mag_column=args.sdss_mag_column,
        sdss_ur_max=args.sdss_ur_max,
        legacy_gr_max=args.legacy_gr_max,
        extended_columns=True,
    )

    sections = {args.section} if args.section != "all" else {"hist", "cdf", "random", "master"}
    if "hist" in sections:
        section_hist(
            hosts,
            sdss_cut,
            tag=args.tag,
            sample_mode=args.sample_mode,
            sdss_q_column=args.sdss_q_column,
            q0=args.q0,
            mag_limit=args.mag_limit,
        )
    if "cdf" in sections:
        section_cdf(
            hosts,
            legacy_cut,
            sdss_cut,
            tag=args.tag,
            mag_limit=args.mag_limit,
            mc_draws_null=args.mc_draws_null,
            mc_alpha=args.mc_alpha,
            sample_mode=args.sample_mode,
            sdss_q_column=args.sdss_q_column,
            q0=args.q0,
            sdss_ur_max=args.sdss_ur_max,
            legacy_gr_max=args.legacy_gr_max,
        )
    if "random" in sections:
        section_random(
            hosts,
            legacy_cut,
            sdss_cut,
            tag=args.tag,
            mag_limit=args.mag_limit,
            n_random=args.n_random,
            mc_draws_null=args.mc_draws_null,
            sample_mode=args.sample_mode,
            sdss_q_column=args.sdss_q_column,
            q0=args.q0,
        )
    if "master" in sections:
        section_master(args.diff_csv)


if __name__ == "__main__":
    main()
