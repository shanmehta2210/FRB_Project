"""
Build master_frb_localization.csv from master_frb_summary.csv.

Coordinate semantics:
  - First N_ORIGINAL_HOST_FRBS rows (canonical original sample): RA_deg/DEC_deg are
    host positions.
  - Later rows: if notes indicate host association is done, RA_deg/DEC_deg are host;
    otherwise they are treated as burst/signal localization centers.

survey: inferred from short `notes` tags — LS, PS1, LS+PS1, unsure, or empty when
unknown. Original "perfect match." sample is labeled LS (Legacy imaging cutouts).

Run from repo root: python scripts/build_master_frb_localization.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

N_ORIGINAL_HOST_FRBS = 23

HOST_ASSOC_MARKERS = (
    "host association done",
    "host association complete",
)


def _combined_notes(row: pd.Series) -> str:
    a = row.get("notes", "")
    b = row.get("Notes", "")
    if pd.isna(a):
        a = ""
    if pd.isna(b):
        b = ""
    return f"{a} {b}"


def host_association_done(row: pd.Series) -> bool:
    text = _combined_notes(row).lower()
    return any(m in text for m in HOST_ASSOC_MARKERS)


def coord_semantics(row: pd.Series, index: int) -> str:
    if index < N_ORIGINAL_HOST_FRBS:
        return "host"
    return "host" if host_association_done(row) else "signal"


def infer_survey(row: pd.Series) -> str:
    """
    Primary imaging / fetch source tag for science cutouts.
    Uses bracket tags from master summary when present; otherwise tags the original
    vetted sample as LS (Legacy Survey–based cutouts).
    """
    text = _combined_notes(row)
    low = text.lower()

    if "[legacy surveys only]" in low:
        return "LS"
    if "[pan-starrs only]" in low:
        return "PS1"
    if "[both ls and ps1]" in low:
        return "LS+PS1"
    if "[none (outside footprint or too faint)]" in low:
        return "unsure"

    # Original 23 + workflow hosts: Legacy-oriented pipeline cutouts
    if "perfect match" in low:
        return "LS"

    return ""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("master_frb_summary.csv"),
        help="Source master summary CSV",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("master_frb_localization.csv"),
        help="Output localization truth table",
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    src = args.input if args.input.is_absolute() else root / args.input
    out = args.output if args.output.is_absolute() else root / args.output

    df = pd.read_csv(src)

    semantics = [coord_semantics(df.iloc[i], i) for i in range(len(df))]
    surveys = [infer_survey(df.iloc[i]) for i in range(len(df))]

    out_df = pd.DataFrame(
        {
            "frb": df["FRB"],
            "ra_deg": df["RA_deg"],
            "dec_deg": df["DEC_deg"],
            "ra_hms": df["RA_hms"],
            "dec_dms": df["DEC_dms"],
            "coord_semantics": semantics,
            "survey": surveys,
            "ra_err_as": df["ra_err_as"],
            "dec_err_as": df["dec_err_as"],
            "major_sigma_as": df["major_sigma_as"],
            "minor_sigma_as": df["minor_sigma_as"],
            "pa_deg": df["pa_deg"],
            "z": df["z"],
            "DM": df["DM"],
            "status": df["status"],
        }
    )

    out_df.to_csv(out, index=False)
    print(f"Wrote {out} ({len(out_df)} rows)")


if __name__ == "__main__":
    main()
