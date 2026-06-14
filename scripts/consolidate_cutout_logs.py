#!/usr/bin/env python3
"""One-off: merge cutout logs into large_cutouts/cutout_registry.csv."""

from __future__ import annotations

import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
CUTOUT = REPO / "large_cutouts"
LOC = REPO / "master_frb_localization.csv"
MASTER_CSV = REPO / "pipeline_scripts" / "new_hosts_master.csv"
REGISTRY = CUTOUT / "cutout_registry.csv"

NO_CUTOUT = frozenset(
    {"20230930A", "20230125D", "20230718A", "20201123A", "20230731A"}
)
NO_CUTOUT_NOTE = (
    "No Legacy/PS1/DES coverage at host position; manual cutout required."
)

COLS = [
    "frb",
    "cohort_46",
    "ra_deg",
    "dec_deg",
    "status",
    "source",
    "layer",
    "resampled",
    "notes",
    "updated_utc",
]


def _ts() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")


def paired(frb: str) -> bool:
    return (CUTOUT / f"{frb}_flux.fits").is_file() and (
        CUTOUT / f"{frb}_invvar.fits"
    ).is_file()


def _load_old() -> dict[str, dict]:
    meta: dict[str, dict] = {}
    for name in ("downloaded_cutouts.csv", "refetch_results.csv"):
        p = CUTOUT / name
        if not p.is_file():
            continue
        df = pd.read_csv(p)
        for _, r in df.iterrows():
            frb = str(r.get("frb", "")).strip()
            if not frb:
                continue
            meta[frb] = {
                "source": str(r.get("source", "") or "").strip(),
                "layer": str(r.get("layer", "") or "").strip(),
                "resampled": r.get("resampled", ""),
                "notes": str(r.get("notes", "") or "").strip(),
                "updated_utc": str(
                    r.get("downloaded_utc") or r.get("utc") or ""
                ).strip(),
            }
    return meta


def main() -> None:
    if not MASTER_CSV.is_file():
        raise SystemExit(f"Missing {MASTER_CSV}; run consolidate_new_hosts_logs.py first.")
    names_46 = pd.read_csv(MASTER_CSV)["frb"].astype(str).tolist()
    loc = pd.read_csv(LOC)
    old = _load_old()
    now = _ts()
    rows: list[dict] = []

    for frb in names_46:
        m = loc.loc[loc["frb"] == frb]
        ra = float(m.iloc[0]["ra_deg"]) if len(m) else ""
        dec = float(m.iloc[0]["dec_deg"]) if len(m) else ""
        prev = old.get(frb, {})
        if frb in NO_CUTOUT:
            rows.append(
                {
                    "frb": frb,
                    "cohort_46": True,
                    "ra_deg": ra,
                    "dec_deg": dec,
                    "status": "no_coverage",
                    "source": "",
                    "layer": "",
                    "resampled": "",
                    "notes": NO_CUTOUT_NOTE,
                    "updated_utc": prev.get("updated_utc") or now,
                }
            )
        elif paired(frb):
            rows.append(
                {
                    "frb": frb,
                    "cohort_46": True,
                    "ra_deg": ra,
                    "dec_deg": dec,
                    "status": "ok",
                    "source": prev.get("source", ""),
                    "layer": prev.get("layer", ""),
                    "resampled": prev.get("resampled", ""),
                    "notes": prev.get("notes", ""),
                    "updated_utc": prev.get("updated_utc") or now,
                }
            )
        else:
            rows.append(
                {
                    "frb": frb,
                    "cohort_46": True,
                    "ra_deg": ra,
                    "dec_deg": dec,
                    "status": "missing",
                    "source": "",
                    "layer": "",
                    "resampled": "",
                    "notes": "cutout pair not on disk",
                    "updated_utc": now,
                }
            )

    pat = re.compile(r"^(\d{4}[0-9A-Za-z]+)_(flux|invvar)\.fits$")
    cohort = set(names_46)
    for p in CUTOUT.glob("*_flux.fits"):
        m = pat.match(p.name)
        if not m:
            continue
        frb = m.group(1)
        if frb in cohort or not paired(frb):
            continue
        lm = loc.loc[loc["frb"] == frb]
        ra = float(lm.iloc[0]["ra_deg"]) if len(lm) else ""
        dec = float(lm.iloc[0]["dec_deg"]) if len(lm) else ""
        prev = old.get(frb, {})
        rows.append(
            {
                "frb": frb,
                "cohort_46": False,
                "ra_deg": ra,
                "dec_deg": dec,
                "status": "on_disk",
                "source": prev.get("source", ""),
                "layer": prev.get("layer", ""),
                "resampled": prev.get("resampled", ""),
                "notes": "",
                "updated_utc": prev.get("updated_utc") or "",
            }
        )

    df = pd.DataFrame(rows, columns=COLS).sort_values(["cohort_46", "frb"], ascending=[False, True])
    df.to_csv(REGISTRY, index=False)
    print(f"Wrote {REGISTRY} ({len(df)} rows, {sum(df.cohort_46)} in 46-host cohort)")


if __name__ == "__main__":
    main()
