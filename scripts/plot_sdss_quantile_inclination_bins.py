#!/usr/bin/env python3
"""
Equal-count SDSS quantile bins on cos(i), with FRB hosts mapped into those bins.

Uses the same strict + color + mag selection as mag-cut CDF plots at m_r < 21:
  SDSS: u-r < 2.3, lnL exp-wins, expAB_r > 0.2, modelMag_r <= mag limit
  FRB:  GALFIT mag <= 21, b/a > 0.2

Default: 8 bins (~12.5% of SDSS per bin). Edges are quantiles of SDSS cos(i).

Outputs under plots/plots_null/v1_hist_inclination/quantile8_mag21_sdss_strict/:
  - quantile_inclination_bins.png
  - bin_edges_summary.csv
  - frb_bin_assignments.csv

Run from repo root:
    python scripts/plot_sdss_quantile_inclination_bins.py
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
    SDSS_UR_MAX_CDF,
    cosi_array_from_df,
    equal_count_quantile_edges,
    prepare_null_strict_color_base,
    read_sdss_null_catalog,
    slice_null_base_by_mag,
)
from pipeline_null_plot_utils import (  # noqa: E402
    DEFAULT_PIPELINE,
    DEFAULT_SDSS,
    PLOTS_NULL,
    frb_hosts_for_cdf,
    load_pipeline_hosts,
    plot_frb_in_sdss_quantile_cosi_bins,
)

MAG_LIMIT_DEFAULT = 21.0
N_BINS_DEFAULT = 8
def quantile_out_dir(mag_limit: float) -> Path:
    tag = int(mag_limit) if mag_limit == int(mag_limit) else mag_limit
    return PLOTS_NULL / "v1_hist_inclination" / f"quantile8_mag{tag}_sdss_strict"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pipeline-csv", type=Path, default=DEFAULT_PIPELINE)
    parser.add_argument("--sdss-csv", type=Path, default=DEFAULT_SDSS)
    parser.add_argument("--mag-limit", type=float, default=MAG_LIMIT_DEFAULT)
    parser.add_argument("--n-bins", type=int, default=N_BINS_DEFAULT)
    parser.add_argument("--q0", type=float, default=0.2)
    parser.add_argument("--sdss-mag-column", default="modelMag_r")
    parser.add_argument("--sdss-q-column", default="expAB_r")
    parser.add_argument("--sdss-ur-max", type=float, default=SDSS_UR_MAX_CDF)
    parser.add_argument("--out-dir", type=Path, default=None)
    args = parser.parse_args()
    if args.out_dir is None:
        args.out_dir = quantile_out_dir(args.mag_limit)

    print("[*] Loading SDSS null (strict + color)...")
    sdss_raw = read_sdss_null_catalog(args.sdss_csv)
    sdss_base = prepare_null_strict_color_base(
        sdss_raw,
        mag_column=args.sdss_mag_column,
        q0=args.q0,
        q_column=args.sdss_q_column,
        is_legacy=False,
        sdss_ur_max=args.sdss_ur_max,
    )
    sdss_cut = slice_null_base_by_mag(
        sdss_base,
        mag_column=args.sdss_mag_column,
        mag_limit=args.mag_limit,
    )
    sdss_cosi = cosi_array_from_df(sdss_cut, q_col=args.sdss_q_column, q0=args.q0)
    print(f"[*] SDSS pool N={len(sdss_cut)} (mag<{args.mag_limit:g})")

    edges = equal_count_quantile_edges(sdss_cosi, n_bins=args.n_bins)
    n_bins_eff = len(edges) - 1
    print(f"[*] {n_bins_eff} equal-count bins on cos(i); edges={np.round(edges, 4).tolist()}")

    hosts = load_pipeline_hosts(args.pipeline_csv)
    frb_plot = frb_hosts_for_cdf(
        hosts,
        sample_mode="strict",
        q0=args.q0,
        mag_limit=args.mag_limit,
    )
    inc = np.clip(
        np.asarray(pd.to_numeric(frb_plot["inc"], errors="coerce"), dtype=float),
        0.0,
        90.0,
    )
    frb_cosi = np.cos(np.radians(inc))
    print(f"[*] FRB hosts N={len(frb_plot)} (strict, mag<{args.mag_limit:g})")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_stem = args.out_dir / "quantile_inclination_bins"

    summary_df, frb_df = plot_frb_in_sdss_quantile_cosi_bins(
        sdss_cosi=sdss_cosi,
        frb_hosts=frb_plot,
        frb_cosi=frb_cosi,
        edges=edges,
        n_bins=args.n_bins,
        mag_limit=args.mag_limit,
        out_stem=out_stem,
    )

    summary_path = args.out_dir / "bin_edges_summary.csv"
    frb_path = args.out_dir / "frb_bin_assignments.csv"
    summary_df.to_csv(summary_path, index=False)
    frb_df.to_csv(frb_path, index=False)
    print(f"[*] Wrote {summary_path}")
    print(f"[*] Wrote {frb_path}")
    print(summary_df.to_string(index=False, float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()
