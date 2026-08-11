"""
Merge confident-host positions from new_confident_hosts.txt (LaTeX tables)
into master_frb_localization.csv.

Rules:
  - Skip Survey == CHIME (case-insensitive).
  - Match FRB ids with optional trailing-letter normalization.
  - On match: apply TeX host coords, z, DM, survey, DM_MW, DM_exgal;
    coord_semantics=host.
  - Append new non-CHIME FRBs not already in CSV.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord

ERR_COLS = (
    "ra_err_as",
    "dec_err_as",
    "major_sigma_as",
    "minor_sigma_as",
    "pa_deg",
)

_ROW_RE = re.compile(
    r"^\s*(\d{8}[A-Z]?)\s*&\s*"
    r"([-\d.]+)\s*&\s*"
    r"([-\d.]+)\s*&\s*"
    r"([\d.]+)\s*&\s*"
    r"([\d.]+)\s*&\s*"
    r"([\d.]+)\s*&\s*"
    r"([\d.]+)\s*&\s*"
    r"([\d.]+)\s*&\s*"
    r"([^&\\]+)",
    re.IGNORECASE,
)

OUTPUT_COLS = [
    "frb",
    "ra_deg",
    "dec_deg",
    "ra_hms",
    "dec_dms",
    "coord_semantics",
    "survey",
    *ERR_COLS,
    "z",
    "DM",
    "DM_MW",
    "DM_exgal",
    "status",
]


@dataclass
class TexRow:
    frb: str
    ra_deg: float
    dec_deg: float
    dm_obs: float
    dm_mw: float
    z: float
    dm_exgal: float
    survey: str


def normalize_frb(name: str) -> str:
    return str(name).strip().upper()


def frb_date_prefix(name: str) -> str:
    n = normalize_frb(name)
    return n[:8] if len(n) >= 8 else n


def has_trailing_letter(name: str) -> bool:
    n = normalize_frb(name)
    return len(n) == 9 and n[8].isalpha()


def parse_tex_file(path: Path) -> list[TexRow]:
    rows: list[TexRow] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        m = _ROW_RE.match(line)
        if not m:
            continue
        rows.append(
            TexRow(
                frb=normalize_frb(m.group(1)),
                ra_deg=float(m.group(2)),
                dec_deg=float(m.group(3)),
                dm_obs=float(m.group(5)),
                dm_mw=float(m.group(6)),
                z=float(m.group(7)),
                dm_exgal=float(m.group(8)),
                survey=m.group(9).strip(),
            )
        )
    return rows


def deg_to_hms_dms(ra_deg: float, dec_deg: float) -> tuple[str, str]:
    c = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
    ra_hms = c.ra.to_string(unit=u.hourangle, sep=":", precision=2, pad=True)
    dec_dms = c.dec.to_string(sep=":", precision=1, alwayssign=True, pad=True)
    return ra_hms, dec_dms


def build_csv_index(frb_series: pd.Series) -> dict[str, list[int]]:
    idx: dict[str, list[int]] = {}
    for i, name in enumerate(frb_series):
        idx.setdefault(normalize_frb(name), []).append(i)
    return idx


def resolve_csv_row(
    tex_frb: str, csv_frb: pd.Series, csv_index: dict[str, list[int]]
) -> tuple[int | None, str, str]:
    """Return (row_index, matched_csv_name, note)."""
    key = normalize_frb(tex_frb)
    if key in csv_index:
        if len(csv_index[key]) == 1:
            i = csv_index[key][0]
            return i, normalize_frb(csv_frb.iloc[i]), ""
        return None, "", f"ambiguous exact match for {key}"

    if has_trailing_letter(key):
        return None, "", "no match"

    prefix = frb_date_prefix(key)
    candidates: list[int] = []
    for csv_key, indices in csv_index.items():
        if frb_date_prefix(csv_key) == prefix and has_trailing_letter(csv_key):
            candidates.extend(indices)
    if len(candidates) == 1:
        i = candidates[0]
        return i, normalize_frb(csv_frb.iloc[i]), f"date alias {key}"
    if len(candidates) > 1:
        names = [normalize_frb(csv_frb.iloc[i]) for i in candidates]
        return None, "", f"ambiguous alias {key} -> {names}"
    return None, "", "no match"


def apply_tex(row: pd.Series, tex: TexRow) -> None:
    row["ra_deg"] = tex.ra_deg
    row["dec_deg"] = tex.dec_deg
    row["ra_hms"], row["dec_dms"] = deg_to_hms_dms(tex.ra_deg, tex.dec_deg)
    row["coord_semantics"] = "host"
    row["z"] = tex.z
    row["DM"] = tex.dm_obs
    row["DM_MW"] = tex.dm_mw
    row["DM_exgal"] = tex.dm_exgal
    row["survey"] = tex.survey


def new_row_from_tex(tex: TexRow) -> dict:
    ra_hms, dec_dms = deg_to_hms_dms(tex.ra_deg, tex.dec_deg)
    return {
        "frb": tex.frb,
        "ra_deg": tex.ra_deg,
        "dec_deg": tex.dec_deg,
        "ra_hms": ra_hms,
        "dec_dms": dec_dms,
        "coord_semantics": "host",
        "survey": tex.survey,
        "ra_err_as": "",
        "dec_err_as": "",
        "major_sigma_as": "",
        "minor_sigma_as": "",
        "pa_deg": "",
        "z": tex.z,
        "DM": tex.dm_obs,
        "DM_MW": tex.dm_mw,
        "DM_exgal": tex.dm_exgal,
        "status": "pending",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tex", type=Path, default=Path("Archive/notes/new_confident_hosts.txt"))
    parser.add_argument("--csv", type=Path, default=Path("master_frb_localization.csv"))
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    tex_path = args.tex if args.tex.is_absolute() else root / args.tex
    csv_path = args.csv if args.csv.is_absolute() else root / args.csv

    tex_rows = parse_tex_file(tex_path)
    usable = [r for r in tex_rows if r.survey.upper() != "CHIME"]

    df = pd.read_csv(csv_path)
    for col in ERR_COLS:
        if col in df.columns:
            df[col] = df[col].astype(object)
    for col in ("DM_MW", "DM_exgal"):
        if col not in df.columns:
            df[col] = ""

    csv_index = build_csv_index(df["frb"])
    new_rows: list[dict] = []
    ambiguous: list[tuple[str, str]] = []

    stats = {"parsed": len(tex_rows), "skipped_chime": len(tex_rows) - len(usable), "updated": 0, "added": 0}

    for tex in usable:
        row_i, csv_name, note = resolve_csv_row(tex.frb, df["frb"], csv_index)
        if row_i is None:
            if note.startswith("ambiguous"):
                ambiguous.append((tex.frb, note))
            else:
                new_rows.append(new_row_from_tex(tex))
                stats["added"] += 1
                print(f"[+] New: {tex.frb} ({tex.survey})")
            continue

        stats["updated"] += 1
        row = df.iloc[row_i].copy()
        if note:
            print(f"[*] {tex.frb} -> {csv_name} ({note})")
        elif csv_name != tex.frb:
            print(f"[*] {tex.frb} -> {csv_name}")
        apply_tex(row, tex)
        df.iloc[row_i] = row

    if new_rows:
        df = pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)

    for c in OUTPUT_COLS:
        if c not in df.columns:
            df[c] = ""
    df = df[OUTPUT_COLS]

    print("\n=== Merge summary ===")
    print(f"TeX rows parsed:     {stats['parsed']}")
    print(f"Skipped (CHIME):     {stats['skipped_chime']}")
    print(f"Updated existing:    {stats['updated']}")
    print(f"New rows added:      {stats['added']}")
    print(f"Output rows:         {len(df)}")

    if ambiguous:
        print("\n=== Ambiguous (skipped) ===")
        for frb, msg in ambiguous:
            print(f"  {frb}: {msg}")

    if not args.dry_run:
        df.to_csv(csv_path, index=False)
        print(f"\nWrote {csv_path}")
    else:
        print("\n(dry-run: CSV not written)")


if __name__ == "__main__":
    main()
