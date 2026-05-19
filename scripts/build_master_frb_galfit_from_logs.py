"""
Build master_frb_galfit_from_logs.csv by parsing tools/galfit/runs/<FRB>/*/fit.log.

Parsers live in galfit_fitlog_parse.py (robust dash splitting + second-to-last
single-sersic block policy, aligned with append_old_frbs_galfit_results.py).

Run from repo root:
  python scripts/build_master_frb_galfit_from_logs.py
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import pandas as pd

# Allow `python scripts/build_master_frb_galfit_from_logs.py` from repo root
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from galfit_fitlog_parse import (
    inclination_err_from_b_a_err,
    inclination_from_b_a,
    parse_fitlog_file,
)

_METRIC_KEYS = [
    "chi2nu",
    "mag",
    "mag_err",
    "re",
    "re_err",
    "n",
    "n_err",
    "b_a",
    "b_a_err",
    "pa",
    "pa_err",
    "x",
    "x_err",
    "y",
    "y_err",
]

# Column order aligned with new_16_frbs_galfit_results.csv: all no_psf_sigma, then all with_psf_sigma.
_ORDERED_COLS = (
    ["frb"]
    + [f"{k}_nopsf" for k in _METRIC_KEYS]
    + [f"{k}_psf" for k in _METRIC_KEYS]
    + [
        "parse_strategy_nopsf",
        "parse_strategy_psf",
        "fit_log_path_nopsf",
        "fit_log_path_psf",
        "inc_nopsf",
        "inc_err_nopsf",
        "inc_psf",
        "inc_err_psf",
    ]
)


def _cell(val):
    if val is None:
        return ""
    if isinstance(val, float) and (math.isnan(val) or not math.isfinite(val)):
        return ""
    return val


def extract_metrics(log_path: Path | None, repo_root: Path) -> tuple[dict, str, str]:
    if log_path is None or not log_path.is_file():
        return {}, "", ""
    data, strategy = parse_fitlog_file(str(log_path))
    try:
        rel_s = str(log_path.resolve().relative_to(repo_root)).replace("\\", "/")
    except ValueError:
        rel_s = str(log_path).replace("\\", "/")
    return data, strategy, rel_s


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--frb-table",
        type=Path,
        default=Path("master_frb_localization.csv"),
        help="CSV whose `frb` column defines row order",
    )
    parser.add_argument(
        "--runs-dir",
        type=Path,
        default=Path("tools/galfit/runs"),
        help="Directory containing per-FRB run folders",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("master_frb_galfit_from_logs.csv"),
        help="Output CSV path",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    frb_csv = args.frb_table if args.frb_table.is_absolute() else root / args.frb_table
    runs_dir = args.runs_dir if args.runs_dir.is_absolute() else root / args.runs_dir
    out_path = args.output if args.output.is_absolute() else root / args.output

    base = pd.read_csv(frb_csv)
    rows = []

    for _, r in base.iterrows():
        frb = str(r["frb"]).strip()
        nopsf_p = runs_dir / frb / "no_psf_sigma" / "fit.log"
        psf_p = runs_dir / frb / "with_psf_sigma" / "fit.log"

        d_n, st_n, path_n = extract_metrics(nopsf_p if nopsf_p.is_file() else None, root)
        d_p, st_p, path_p = extract_metrics(psf_p if psf_p.is_file() else None, root)

        row = {"frb": frb}
        for k in _METRIC_KEYS:
            row[f"{k}_nopsf"] = _cell(d_n.get(k))
            row[f"{k}_psf"] = _cell(d_p.get(k))

        row["parse_strategy_nopsf"] = st_n
        row["parse_strategy_psf"] = st_p
        row["fit_log_path_nopsf"] = path_n
        row["fit_log_path_psf"] = path_p

        ba_n = d_n.get("b_a")
        be_n = d_n.get("b_a_err")
        ba_p = d_p.get("b_a")
        be_p = d_p.get("b_a_err")

        row["inc_nopsf"] = (
            _cell(inclination_from_b_a(ba_n)) if ba_n is not None else ""
        )
        row["inc_err_nopsf"] = (
            _cell(inclination_err_from_b_a_err(ba_n, be_n))
            if ba_n is not None
            else ""
        )
        row["inc_psf"] = _cell(inclination_from_b_a(ba_p)) if ba_p is not None else ""
        row["inc_err_psf"] = (
            _cell(inclination_err_from_b_a_err(ba_p, be_p))
            if ba_p is not None
            else ""
        )

        rows.append(row)

    out_df = pd.DataFrame(rows)
    out_df = out_df[[c for c in _ORDERED_COLS if c in out_df.columns]]
    out_df.to_csv(out_path, index=False)
    print(f"Wrote {out_path} ({len(out_df)} rows)")


if __name__ == "__main__":
    main()
