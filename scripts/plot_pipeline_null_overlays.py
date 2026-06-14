"""
Deprecated: use ``plot_pipeline_diagnostics.py`` (inclination CDFs as ``null_cdf_pipeline_*``).

This script previously wrote duplicate ``null_overlay_pipeline_*`` figures identical to
``null_cdf_pipeline_*`` in ``plots/plots_null/v1_null_cdf_inclination/``.
"""

from __future__ import annotations

import sys


def main() -> None:
    print(
        "plot_pipeline_null_overlays.py is deprecated.\n"
        "Run: python scripts/plot_pipeline_diagnostics.py --section cdf\n"
        "Outputs: plots/plots_null/v1_null_cdf_inclination/null_cdf_pipeline_{legacy|sdss}_*.png"
    )
    sys.exit(0)


if __name__ == "__main__":
    main()
