#!/usr/bin/env python3
"""Re-scan pipeline outputs and refresh pipeline_scripts/new_hosts_master.{csv,md}."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
CONSOLIDATE = REPO / "scripts" / "consolidate_new_hosts_logs.py"


def main() -> None:
    subprocess.run([sys.executable, str(CONSOLIDATE)], check=True)


if __name__ == "__main__":
    main()
