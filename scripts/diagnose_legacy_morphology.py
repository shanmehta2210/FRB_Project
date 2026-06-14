#!/usr/bin/env python3
"""
Legacy Tractor morphology diagnostic for null-catalog cuts.

Reports tractor_type and Sérsic n (rdVrad) counts for full catalog and
strict + g-r + mag < 21 reference pool. Flags re-query if union pool is small.

Run from repo root:
    python scripts/diagnose_legacy_morphology.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (  # noqa: E402
    LEGACY_GR_MAX_CDF,
    LEGACY_MORPH_N_MAX,
    LEGACY_MORPH_N_MIN,
    Q0,
    apply_legacy_type_cut,
    apply_strict_q_cut,
    filter_legacy_gr,
    filter_legacy_spiral_morph,
    legacy_spiral_morph_mask,
    read_legacy_null_catalog,
    resolve_mag_column,
)
from pipeline_null_plot_utils import DEFAULT_LEGACY, PLOTS_NULL  # noqa: E402

OUT_DIR = (
    PLOTS_NULL
    / "v1_null_cdf_inclination"
    / "diagnostics"
    / "legacy_morphology"
)
MAG_REF = 21.0
MIN_UNION_POOL = 5_000
MIN_FRAC_OF_STRICT = 0.35


def _type_counts(df: pd.DataFrame) -> pd.Series:
    if "tractor_type" not in df.columns:
        return pd.Series(dtype=int)
    return df["tractor_type"].astype(str).str.upper().value_counts()


def _summary_block(df: pd.DataFrame, label: str) -> list[dict]:
    n = len(df)
    types = _type_counts(df)
    n_mask = legacy_spiral_morph_mask(df) if n else pd.Series(dtype=bool)
    n_in_range = (
        pd.to_numeric(df.get("rdVrad"), errors="coerce").between(
            LEGACY_MORPH_N_MIN, LEGACY_MORPH_N_MAX
        )
        if n
        else pd.Series(dtype=bool)
    )
    is_exp = df["tractor_type"].astype(str).str.upper() == "EXP" if n else pd.Series(dtype=bool)
    rows = [
        {"pool": label, "metric": "n_total", "value": n},
        {"pool": label, "metric": "n_exp_type", "value": int(is_exp.sum()) if n else 0},
        {
            "pool": label,
            "metric": "n_n_in_range",
            "value": int(n_in_range.sum()) if n else 0,
        },
        {
            "pool": label,
            "metric": "n_exp_or_n_range",
            "value": int(n_mask.sum()) if n else 0,
        },
    ]
    for t, c in types.items():
        rows.append({"pool": label, "metric": f"type_{t}", "value": int(c)})
    return rows


def reference_pool(df: pd.DataFrame, mag_col: str) -> pd.DataFrame:
    mag_c = resolve_mag_column(df, mag_col)
    mag = pd.to_numeric(df[mag_c], errors="coerce")
    out = df.loc[mag < MAG_REF].copy()
    out = apply_strict_q_cut(out, q_col="expAB_r", q0=Q0)
    out = filter_legacy_gr(out, LEGACY_GR_MAX_CDF)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-csv", type=Path, default=DEFAULT_LEGACY)
    parser.add_argument("--mag-column", default="rmag")
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    df = read_legacy_null_catalog(args.legacy_csv, extended=True)
    ref = reference_pool(df, args.mag_column)
    ref_no_rex = apply_legacy_type_cut(ref, exclude="REX")
    ref_cdf = apply_legacy_type_cut(ref, exclude="REX,DEV")
    ref_spiral = filter_legacy_spiral_morph(ref_cdf)

    rows: list[dict] = []
    rows.extend(_summary_block(df, "full_catalog"))
    rows.extend(_summary_block(ref, f"strict_gr_mag_lt_{MAG_REF:g}"))
    rows.extend(_summary_block(ref_spiral, "cdf_eligible_exp_or_n"))

    summary = pd.DataFrame(rows)
    summary.to_csv(args.out_dir / "morphology_counts.csv", index=False)

    n_strict = len(ref)
    n_union = len(ref_spiral)
    frac = n_union / n_strict if n_strict else 0.0
    requery = n_union < MIN_UNION_POOL or frac < MIN_FRAC_OF_STRICT

    md = [
        "# Legacy morphology diagnostic",
        "",
        f"Reference pool: strict b/a > {Q0}, g-r < {LEGACY_GR_MAX_CDF}, mag < {MAG_REF:g}.",
        "",
        f"- N strict+color @ m<{MAG_REF:g}: **{n_strict}**",
        f"- N after REX+DEV drop + EXP or n in [{LEGACY_MORPH_N_MIN}, {LEGACY_MORPH_N_MAX}]: **{n_union}** ({frac:.1%} of strict)",
        "",
        f"Re-query gate (&lt; {MIN_UNION_POOL} or &lt; {MIN_FRAC_OF_STRICT:.0%} of strict): "
        f"**{'YES — consider build_legacy_catalog_csv --top' if requery else 'NO — proceed'}**",
    ]
    (args.out_dir / "legacy_morphology_summary.md").write_text(
        "\n".join(md), encoding="utf-8"
    )
    print("\n".join(md))
    print(f"Wrote {args.out_dir}")


if __name__ == "__main__":
    main()
