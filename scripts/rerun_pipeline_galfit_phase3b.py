"""
Re-run Phase 3b (GALFIT fit + sky QA) on existing pipeline Output/*_all folders.

Uses each folder's zero_points.json zp_aper and writes galfit_config.yaml
before invoking run_galfit_fitting.py.

    python scripts/rerun_pipeline_galfit_phase3b.py
    python scripts/rerun_pipeline_galfit_phase3b.py --frb 20200906A
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUT_ROOT = REPO_ROOT / "pipeline_scripts" / "Output"
PHASE3B = REPO_ROOT / "pipeline_scripts" / "galfit_fitting" / "run_galfit_fitting.py"
DEFAULT_GALFIT_CFG = REPO_ROOT / "pipeline_scripts" / "galfit_fitting" / "galfit_config.yaml"


def write_galfit_config(odir: Path) -> float | None:
    cfg = {}
    if DEFAULT_GALFIT_CFG.is_file():
        with open(DEFAULT_GALFIT_CFG, encoding="utf-8") as f:
            loaded = yaml.safe_load(f) or {}
        if isinstance(loaded, dict):
            cfg = loaded
    zp_path = odir / "zero_points.json"
    if zp_path.is_file():
        zp_data = json.loads(zp_path.read_text(encoding="utf-8"))
        zp = zp_data.get("zp_aper") or zp_data.get("zp_aper_40px")
        if zp is not None and math.isfinite(float(zp)):
            cfg["mag_zeropoint"] = float(zp)
    with open(odir / "galfit_config.yaml", "w", encoding="utf-8") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)
    return cfg.get("mag_zeropoint")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--frb", action="append", help="Only re-run these FRB names")
    args = parser.parse_args()
    only = set(args.frb) if args.frb else None

    ok, fail, skip = 0, 0, 0
    for odir in sorted(OUTPUT_ROOT.glob("*_all")):
        frb = odir.name.replace("_all", "")
        if only and frb not in only:
            continue
        if not (odir / "host_cutout.fits").is_file():
            print(f"[skip] {frb}: no host_cutout.fits")
            skip += 1
            continue
        zp = write_galfit_config(odir)
        zp_str = f"{zp:.4f}" if zp is not None else "?"
        print(f"\n=== {frb} (J)={zp_str}) ===")
        rc = subprocess.run(
            [sys.executable, str(PHASE3B), "--dir", str(odir)],
            cwd=str(odir),
        ).returncode
        if rc == 0:
            ok += 1
        else:
            fail += 1
            print(f"[!] {frb}: exit {rc}")

    print(f"\nDone: ok={ok}, failed={fail}, skipped={skip}")
    return 1 if fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
