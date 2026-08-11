"""
One-shot audit: master_frb_galfit_from_logs.csv vs fresh fit.log parse and new_16 table.

Round-trip: every row with fit_log_path_* is re-parsed; all metrics must match.

Cross-check new_16_frbs_galfit_results.csv only for the 16 expansion FRBs (hard-coded
list matching scripts/run_galfit_16_expansion.py). The full new_16 file also contains
legacy-appended rows that may still carry stale chi^2 blow-ups or NaNs — do not expect
those rows to match without rebuilding new_16 from logs.

Run from repo root: python scripts/audit_master_frb_galfit_csv.py
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pandas as pd

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(_ROOT / "scripts"))

from galfit_fitlog_parse import (  # noqa: E402
    inclination_from_b_a,
    parse_fitlog_file,
)

# Expansion-set FRBs — must match scripts/run_galfit_16_expansion.py `targets`
EXPANSION_16_FRBS = (
    "20190611B",
    "20190711A",
    "20200430A",
    "20220105A",
    "20220725A",
    "20221106A",
    "20230526A",
    "20230708A",
    "20230902A",
    "20231226A",
    "20240201A",
    "20240208A",
    "20240210A",
    "20240304A",
    "20240310A",
    "20240318A",
)

METRICS = [
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


def as_float(x) -> float | None:
    if x is None:
        return None
    if isinstance(x, float) and (math.isnan(x) or not math.isfinite(x)):
        return None
    if isinstance(x, str) and not str(x).strip():
        return None
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def close(a, b, tol: float = 5e-4) -> bool:
    fa, fb = as_float(a), as_float(b)
    if fa is None and fb is None:
        return True
    if fa is None or fb is None:
        return False
    if abs(fa - fb) < 1e-6:
        return True
    return abs(fa - fb) <= tol * max(1.0, abs(fa))


def main() -> None:
    master_path = _ROOT / "master_frb_galfit_from_logs.csv"
    new16_path = _ROOT / "Archive" / "csv" / "galfit" / "new_16_frbs_galfit_results.csv"

    df = pd.read_csv(master_path)
    issues: list[tuple] = []
    roundtrip_fail: list[tuple] = []

    for _, row in df.iterrows():
        frb = row["frb"]
        for suf in ("nopsf", "psf"):
            rel = row.get(f"fit_log_path_{suf}")
            if not isinstance(rel, str) or not str(rel).strip():
                continue
            path = _ROOT / str(rel).replace("/", "\\")
            if not path.is_file():
                issues.append((frb, suf, "missing_file", rel))
                continue
            data, st = parse_fitlog_file(str(path))
            exp_st = row.get(f"parse_strategy_{suf}")
            if str(exp_st) != str(st):
                issues.append((frb, suf, f"strategy csv={exp_st!r} fresh={st!r}", ""))
            for k in METRICS:
                v_csv = row.get(f"{k}_{suf}")
                v_raw = data.get(k)
                if not close(v_csv, v_raw if v_raw is not None else None):
                    roundtrip_fail.append((frb, suf, k, v_csv, v_raw))
            ba = data.get("b_a")
            inc_e = inclination_from_b_a(ba) if ba is not None else None
            inc_csv = as_float(row.get(f"inc_{suf}"))
            if inc_e is not None and inc_csv is not None:
                if not close(inc_e, inc_csv, tol=2e-3):
                    roundtrip_fail.append((frb, suf, "inc_recalc", inc_csv, inc_e))

    print("=== Round-trip / strategy audit ===")
    print("missing_log_files:", len([x for x in issues if x[2] == "missing_file"]))
    print("strategy_mismatch:", len([x for x in issues if "strategy" in str(x[2])]))
    print("value_mismatches:", len(roundtrip_fail))
    for x in issues:
        print("  ISSUE", x)
    for x in roundtrip_fail[:40]:
        print("  MISMATCH", x)
    if len(roundtrip_fail) > 40:
        print("  ...", len(roundtrip_fail) - 40, "more mismatches")

    if new16_path.is_file():
        new16 = pd.read_csv(new16_path)
        m = df.set_index("frb")
        ncheck: list[tuple] = []
        for _, nr in new16.iterrows():
            frb = nr["FRB"]
            if frb not in EXPANSION_16_FRBS:
                continue
            if frb not in m.index:
                continue
            mr = m.loc[frb]
            for col in new16.columns:
                if col == "FRB":
                    continue
                if col not in mr.index:
                    continue
                a, b = nr[col], mr[col]
                if pd.isna(a) and pd.isna(b):
                    continue
                if pd.isna(a) or pd.isna(b):
                    ncheck.append((frb, col, a, b))
                    continue
                fa, fb = float(a), float(b)
                if abs(fa - fb) > max(0.002, 1e-3 * max(1.0, abs(fa))):
                    ncheck.append((frb, col, fa, fb))

        print("\n=== new_16 (expansion 16 only) vs master ===")
        print("mismatch_count:", len(ncheck))
        for x in ncheck:
            print("  ", x)


if __name__ == "__main__":
    main()
