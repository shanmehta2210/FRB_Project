#!/usr/bin/env python3
"""
Archive pre-morphology-cut null plots before regenerating.

Moves plots/plots_null/ to Archive/plots_null_pre_morphology_cut/<timestamp>/.
Recreates an empty plots/plots_null/ directory.

Run from repo root:
    python scripts/archive_null_plots_pre_morphology.py
"""

from __future__ import annotations

import argparse
import shutil
import sys
from datetime import datetime
from pathlib import Path

_REPO = Path(__file__).resolve().parents[1]
PLOTS_NULL = _REPO / "plots" / "plots_null"
ARCHIVE_ROOT = _REPO / "Archive" / "plots_null_pre_morphology_cut"
TEST_RESULTS = _REPO / "test_results.md"

# Subdirs that must not remain after a successful archive+recreate.
LEGACY_SUBDIRS = (
    "mag_cuts",
    "v1_null_cdf_inclination",
    "v1_hist_inclination",
    "v1_null_plots",
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned moves without writing.",
    )
    args = parser.parse_args()

    if not PLOTS_NULL.is_dir():
        print(f"[!] {PLOTS_NULL} does not exist; creating empty tree.")
        if not args.dry_run:
            PLOTS_NULL.mkdir(parents=True, exist_ok=True)
        return

    has_content = any(PLOTS_NULL.iterdir())
    if not has_content:
        print(f"{PLOTS_NULL} already empty.")
        return

    stamp = datetime.now().strftime("%Y%m%d_%H%M")
    dest_root = ARCHIVE_ROOT / stamp
    dest_plots = dest_root / "plots_null"

    print(f"Archive destination: {dest_plots}")
    if args.dry_run:
        return

    dest_root.mkdir(parents=True, exist_ok=False)
    shutil.move(str(PLOTS_NULL), str(dest_plots))
    PLOTS_NULL.mkdir(parents=True, exist_ok=True)

    if TEST_RESULTS.is_file():
        shutil.move(str(TEST_RESULTS), str(dest_root / "test_results_pre_morphology.md"))

    readme = ARCHIVE_ROOT / "README.md"
    readme.write_text(
        "\n".join(
            [
                "# Archived null plots (pre morphology cut)",
                "",
                f"Latest archive: `{stamp}/`",
                "",
                "Contents: SDSS `best_model_ba_r` / no lnL exp-winner cut; "
                "Legacy REX-only type exclusion (no DEV drop, no EXP∪n filter).",
                "",
                "Regenerated outputs live under `plots/plots_null/` after morphology cuts.",
            ]
        ),
        encoding="utf-8",
    )

    remaining = [d for d in LEGACY_SUBDIRS if (PLOTS_NULL / d).exists()]
    if remaining:
        raise SystemExit(
            f"Safety check failed: plots/plots_null still contains {remaining}"
        )

    print(f"Archived to {dest_root}")
    print(f"Fresh {PLOTS_NULL} ready for regeneration.")


if __name__ == "__main__":
    main()
