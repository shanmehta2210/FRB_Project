#!/usr/bin/env python3
"""
Generate strict null inclination CDFs (Legacy + SDSS) with late-type color cuts.

SDSS null: u-r < 2.3, lnL exp-wins, expAB_r > 0.2, modelMag_r <= mag limit.
Legacy null: g-r < 0.75, EXP or n in [0.75,2], no REX/DEV, expAB_r > 0.2.

Layout:
    mag_cuts/mag21/legacy_strict/null_cdf_inclination.png
    mag_cuts/mag21/sdss_strict/null_cdf_inclination.png

Run from repo root:
    python scripts/plot_null_mag_cut_cdfs.py
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

import pandas as pd

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (  # noqa: E402
    LEGACY_CDF_TYPE_EXCLUDE,
    LEGACY_GR_MAX_CDF,
    SDSS_Q_COLUMN_CDF,
    SDSS_UR_MAX_CDF,
    cosi_array_from_df,
    prepare_null_strict_color_base,
    read_legacy_null_catalog,
    read_sdss_null_catalog,
    slice_null_base_by_mag,
)
from pipeline_null_plot_utils import (  # noqa: E402
    DEFAULT_LEGACY,
    DEFAULT_PIPELINE,
    DEFAULT_SDSS,
    PLOTS_NULL,
    build_frb_mc_draws_inc,
    cdf_envelope,
    frb_hosts_for_cdf,
    load_pipeline_hosts,
    mc_mean_cdf_from_draws,
    plot_inclination_cdf_overlay,
)

MAG_CUTS_DEFAULT = [24, 23, 22, 21, 20, 19, 18, 17, 16, 15]
OUT_ROOT = PLOTS_NULL / "v1_null_cdf_inclination" / "mag_cuts"
SAMPLE_MODE = "strict"

SURVEYS = {
    "legacy": {"display": "Legacy", "color": "#377eb8"},
    "sdss": {"display": "SDSS", "color": "#4daf4a"},
}


def plot_one(
    *,
    hosts: pd.DataFrame,
    null_cut: pd.DataFrame,
    survey_key: str,
    mag_limit: float,
    mag_label: str,
    shape_label: str,
    q_column: str,
    color_note: str,
    mc_draws_null: int,
    mc_alpha: float,
    out_dir: Path,
    q0: float = 0.2,
    frb_mag_limit: float | None = None,
) -> dict:
    meta = SURVEYS[survey_key]
    frb_plot = frb_hosts_for_cdf(
        hosts,
        sample_mode=SAMPLE_MODE,
        q0=q0,
        mag_limit=frb_mag_limit,
    )
    n_frb = len(frb_plot)
    frb_draws = build_frb_mc_draws_inc(frb_plot)
    x_frb, y_frb = mc_mean_cdf_from_draws(frb_draws)

    cosi = cosi_array_from_df(null_cut, q_col=q_column, q0=q0)
    x_n, mn, lo, hi = cdf_envelope(cosi, n_sample=n_frb, n_draws=mc_draws_null)

    out_dir.mkdir(parents=True, exist_ok=True)
    out_stem = out_dir / "null_cdf_inclination"

    frb_lbl = f"FRB hosts (N={n_frb}"
    if frb_mag_limit is not None:
        frb_lbl += f", GALFIT $m_r < {frb_mag_limit:g}$"
    frb_lbl += f", $b/a > {q0:g}$)"

    plot_inclination_cdf_overlay(
        null_label=(
            f"{meta['display']} null (strict, {color_note}, N={len(null_cut)})"
        ),
        null_color=meta["color"],
        x_null=x_n,
        mean_null=mn,
        lo_null=lo,
        hi_null=hi,
        frb_draws=frb_draws,
        x_frb=x_frb,
        y_frb=y_frb,
        n_frb=n_frb,
        title=(
            f"{frb_lbl} vs {meta['display']} null "
            f"({mag_label} $<$ {mag_limit:g}, {shape_label}, strict)"
        ),
        out_stem=out_stem,
        mc_alpha=mc_alpha,
    )
    return {
        "survey": survey_key,
        "sample_mode": SAMPLE_MODE,
        "mag_limit": mag_limit,
        "n_null_pool": len(null_cut),
        "n_frb": n_frb,
        "output_dir": str(out_dir),
    }


def clear_out_root(out_root: Path) -> None:
    if out_root.is_dir():
        for child in out_root.iterdir():
            if child.is_dir():
                shutil.rmtree(child)
            else:
                child.unlink()
    else:
        out_root.mkdir(parents=True, exist_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pipeline-csv", type=Path, default=DEFAULT_PIPELINE)
    parser.add_argument("--legacy-csv", type=Path, default=DEFAULT_LEGACY)
    parser.add_argument("--sdss-csv", type=Path, default=DEFAULT_SDSS)
    parser.add_argument(
        "--mag-limits",
        type=float,
        nargs="+",
        default=MAG_CUTS_DEFAULT,
    )
    parser.add_argument("--mag-column", default="rmag")
    parser.add_argument("--q0", type=float, default=0.2)
    parser.add_argument("--exclude-types", default=LEGACY_CDF_TYPE_EXCLUDE)
    parser.add_argument("--sdss-mag-column", default="modelMag_r")
    parser.add_argument("--sdss-q-column", default=SDSS_Q_COLUMN_CDF)
    parser.add_argument("--sdss-ur-max", type=float, default=SDSS_UR_MAX_CDF)
    parser.add_argument("--legacy-gr-max", type=float, default=LEGACY_GR_MAX_CDF)
    parser.add_argument("--mc-draws-null", type=int, default=10000)
    parser.add_argument("--mc-alpha", type=float, default=0.03)
    parser.add_argument("--out-root", type=Path, default=OUT_ROOT)
    parser.add_argument(
        "--no-clear",
        action="store_true",
        help="Do not delete existing contents of out-root before running.",
    )
    args = parser.parse_args()

    if not args.no_clear:
        clear_out_root(args.out_root)
    else:
        args.out_root.mkdir(parents=True, exist_ok=True)

    print("[*] Loading null catalogs once (minimal usecols)...")
    legacy_raw = read_legacy_null_catalog(args.legacy_csv)
    sdss_raw = read_sdss_null_catalog(args.sdss_csv)

    legacy_base = prepare_null_strict_color_base(
        legacy_raw,
        mag_column=args.mag_column,
        q0=args.q0,
        q_column="expAB_r",
        exclude_legacy_types=args.exclude_types,
        is_legacy=True,
        legacy_gr_max=args.legacy_gr_max,
    )
    sdss_base = prepare_null_strict_color_base(
        sdss_raw,
        mag_column=args.sdss_mag_column,
        q0=args.q0,
        q_column=args.sdss_q_column,
        is_legacy=False,
        sdss_ur_max=args.sdss_ur_max,
        sdss_exp_winner_only=True,
    )
    del legacy_raw, sdss_raw

    print(
        f"[*] Strict+color base pools: Legacy N={len(legacy_base)}, "
        f"SDSS N={len(sdss_base)}"
    )

    hosts = load_pipeline_hosts(args.pipeline_csv)
    rows = []

    sdss_color = f"$u-r<{args.sdss_ur_max:g}$"
    legacy_color = f"$g-r<{args.legacy_gr_max:g}$"
    survey_meta = {
        "legacy": (
            "tractor $m_r$",
            "Tractor $b/a$",
            "expAB_r",
            f"{legacy_color}, EXP or $n\\in[0.75,2]$, no REX/DEV",
        ),
        "sdss": (
            "modelMag$_r$",
            "exp $b/a$",
            args.sdss_q_column,
            f"{sdss_color}, lnL exp-wins",
        ),
    }

    for mag_limit in sorted(args.mag_limits, reverse=True):
        mag_tag = int(mag_limit)
        frb_mag = mag_limit
        n_frb_need = len(
            frb_hosts_for_cdf(
                hosts,
                sample_mode=SAMPLE_MODE,
                q0=args.q0,
                mag_limit=frb_mag,
            )
        )

        cuts = {
            "legacy": slice_null_base_by_mag(
                legacy_base,
                mag_column=args.mag_column,
                mag_limit=mag_limit,
            ),
            "sdss": slice_null_base_by_mag(
                sdss_base,
                mag_column=args.sdss_mag_column,
                mag_limit=mag_limit,
            ),
        }

        for survey_key, cut in cuts.items():
            if len(cut) < n_frb_need:
                print(
                    f"[!] {survey_key} n={len(cut)} < n_frb={n_frb_need} "
                    f"(mag<{mag_limit}, strict+color)"
                )
            out_dir = args.out_root / f"mag{mag_tag}" / f"{survey_key}_{SAMPLE_MODE}"
            try:
                mag_lbl, shape_lbl, q_col, color_note = survey_meta[survey_key]
                row = plot_one(
                    hosts=hosts,
                    null_cut=cut,
                    survey_key=survey_key,
                    mag_limit=mag_limit,
                    mag_label=mag_lbl,
                    shape_label=shape_lbl,
                    q_column=q_col,
                    color_note=color_note,
                    q0=args.q0,
                    mc_draws_null=args.mc_draws_null,
                    mc_alpha=args.mc_alpha,
                    out_dir=out_dir,
                    frb_mag_limit=frb_mag,
                )
                rows.append(row)
                print(f"[*] {out_dir}")
            except RuntimeError as exc:
                print(f"[FAIL] mag{mag_tag}/{survey_key}_{SAMPLE_MODE}: {exc}")

    summary_path = args.out_root / "mag_cut_summary.csv"
    new_df = pd.DataFrame(rows)
    if args.no_clear and summary_path.is_file():
        old = pd.read_csv(summary_path)
        new_df = pd.concat([old, new_df], ignore_index=True)
        new_df = new_df.drop_duplicates(
            subset=["survey", "sample_mode", "mag_limit"],
            keep="last",
        )
    new_df.to_csv(summary_path, index=False)
    print(f"[*] Done: {len(rows)} plots under {args.out_root}")


if __name__ == "__main__":
    main()
