#!/usr/bin/env python3
"""Download cutouts for FRBs in new_hosts_master.csv that are not yet on disk."""

from __future__ import annotations

import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from cutout_download import (
    AUDIT_CSV,
    CUTOUT_DIR,
    LOC_CSV,
    fetch_frb,
    load_registry,
    paired_on_disk,
    save_registry,
    upsert_row,
)
from cutout_fetch_common import PS1_DEC_MIN

REPO = _SCRIPTS.parent
MASTER_CSV = REPO / "pipeline_scripts" / "new_hosts_master.csv"
def _ts() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")


def _likely_fast_fail(frb: str, dec: float) -> bool:
    """Southern / audit-no-PS1 fields finish in ~30s; high-north PS1 ~10+ min each."""
    if dec <= PS1_DEC_MIN:
        return True
    if not AUDIT_CSV.is_file():
        return False
    m = pd.read_csv(AUDIT_CSV)
    m = m.loc[m["frb"] == frb]
    if m.empty:
        return False
    row = m.iloc[0]
    ps1_ok = bool(row.get("ps1_ok", False)) and dec > PS1_DEC_MIN
    leg_ok = bool(row.get("legacy_global_ok", False)) or bool(
        row.get("legacy_hemisphere_ok", False)
    )
    return not ps1_ok and not leg_ok


def missing_from_list() -> list[tuple[str, float, float]]:
    if MASTER_CSV.is_file():
        names = pd.read_csv(MASTER_CSV)["frb"].astype(str).tolist()
    else:
        raise SystemExit(f"Missing {MASTER_CSV}")
    loc = pd.read_csv(LOC_CSV)
    out = []
    for frb in names:
        if frb in paired_on_disk():
            continue
        m = loc.loc[loc["frb"] == frb]
        if m.empty:
            continue
        out.append((frb, float(m.iloc[0]["ra_deg"]), float(m.iloc[0]["dec_deg"])))
    # Quick failures first so the batch log shows progress before long PS1 jobs.
    out.sort(key=lambda t: (not _likely_fast_fail(t[0], t[2]), -t[2]))
    return out


def main():
    targets = missing_from_list()
    print(f"[batch] {len(targets)} FRB(s) missing cutouts")
    if not targets:
        return

    reg = load_registry()
    ok = fail = 0
    for i, (frb, ra, dec) in enumerate(targets, 1):
        print(f"\n{'='*60}\n[{i}/{len(targets)}] {frb}\n{'='*60}")
        try:
            success, source, layer, resampled = fetch_frb(frb, ra, dec, force=True)
            on_disk = frb in paired_on_disk()
            if success and on_disk:
                reg = upsert_row(
                    reg, frb, ra, dec,
                    source=source or "unknown",
                    layer=layer,
                    resampled=resampled,
                    status="ok",
                )
                ok += 1
                print(f">>> SUCCESS {frb} ({source})")
            else:
                reg = upsert_row(
                    reg, frb, ra, dec, source="none", layer="",
                    resampled=False, status="failed",
                )
                fail += 1
                print(f">>> FAIL {frb}")
        except Exception as exc:
            fail += 1
            print(f">>> ERROR {frb}: {exc}")
        save_registry(reg)

    print(f"\n[batch] done: {ok} ok, {fail} fail | registry {CUTOUT_DIR / 'cutout_registry.csv'}")


if __name__ == "__main__":
    main()
