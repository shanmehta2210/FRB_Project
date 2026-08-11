"""Shared Sérsic initial-guess helpers for Phase 3a / 3b.

Keep Phase 3a neighbor-Re tests and Phase 3b GALFIT ``re`` seeds on the same
``FLUX_RADIUS`` recipe so the two cannot drift.
"""

from __future__ import annotations

import math
from typing import Any, Mapping


def effective_re_px(row: Mapping[str, Any] | Any) -> float:
    """Ordinary GALFIT ``re`` seed (pixels) from a SExtractor catalog row.

    Matches ``run_galfit_fitting.py`` for non-extended-host components:
    ``FLUX_RADIUS``, floored to ``1.0`` when missing / non-positive.

    The extended-host ``max(re, sqrt(AWIN*BWIN))`` boost stays in Phase 3b only
    (host-specific); neighbor ROI tests use this ordinary seed.
    """
    try:
        re = float(row["FLUX_RADIUS"])
    except (KeyError, TypeError, ValueError):
        return 1.0
    if not math.isfinite(re) or re <= 0.0:
        return 1.0
    return re
