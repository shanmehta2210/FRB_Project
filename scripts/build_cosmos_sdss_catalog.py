#!/usr/bin/env python3
"""Build SDSS DR17 pools for COSMOS ACS footprint (b/a audit).

Writes entire (all footprint galaxies with valid shape) and strict (b/a > 0.2).

Run from repo root::

    python scripts/build_cosmos_sdss_catalog.py --no-color-cut
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astroquery.sdss import SDSS

_SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS))

from cosmos_ba_audit_utils import (  # noqa: E402
    COSMOS_BA_MIN,
    COSMOS_DEC_MAX,
    COSMOS_DEC_MIN,
    COSMOS_RA_MAX,
    COSMOS_RA_MIN,
    COSMOS_UR_MAX,
    SDSS_CSV,
    SDSS_DISK_CSV,
    SDSS_DISK_ENTIRE_CSV,
    SDSS_ENTIRE_CSV,
    apply_ba_strict,
)
from null_catalog_utils import (  # noqa: E402
    ensure_sdss_colors,
    filter_sdss_drop_dev_winners,
    filter_sdss_ur,
)

SDSS_DR = 17

SQL = f"""
SELECT
    p.objID,
    p.ra,
    p.dec,
    p.modelMag_r,
    p.modelMag_u,
    p.modelMag_g,
    p.expAB_r,
    p.expRad_r,
    p.lnLExp_r,
    p.lnLDeV_r,
    p.type,
    p.clean,
    p.mode
FROM PhotoObj AS p
WHERE p.ra BETWEEN {COSMOS_RA_MIN} AND {COSMOS_RA_MAX}
  AND p.dec BETWEEN {COSMOS_DEC_MIN} AND {COSMOS_DEC_MAX}
  AND p.type = 3
  AND p.clean = 1
  AND p.mode = 1
  AND p.cmodelMag_r > 0 AND p.cmodelMag_r < 90
  AND p.deVAB_r > 0 AND p.deVAB_r <= 1
  AND p.expAB_r > 0 AND p.expAB_r <= 1
"""


def apply_cosmos_sdss_entire(
    df: pd.DataFrame, *, color: bool, disk_only: bool
) -> tuple[pd.DataFrame, list[dict]]:
    """Footprint galaxies with finite mag and expAB_r; optional exp-winner disk cut."""
    funnel: list[dict] = []
    n0 = len(df)
    funnel.append({"stage": "footprint_sql", "n_remaining": n0})

    out = ensure_sdss_colors(df.copy())

    mag = pd.to_numeric(out["modelMag_r"], errors="coerce")
    ba = pd.to_numeric(out["expAB_r"], errors="coerce")
    ok = mag.notna() & ba.notna() & (ba > 0) & (ba <= 1)
    out = out.loc[ok].copy()
    funnel.append({"stage": "finite_mag_expAB", "n_remaining": len(out)})

    if color:
        out = filter_sdss_ur(out, COSMOS_UR_MAX)
        funnel.append({"stage": f"u_r_lt_{COSMOS_UR_MAX:g}", "n_remaining": len(out)})

    if disk_only:
        out = filter_sdss_drop_dev_winners(out)
        funnel.append({"stage": "lnL_exp_wins_disk", "n_remaining": len(out)})

    return out.reset_index(drop=True), funnel


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-entire", type=Path, default=SDSS_ENTIRE_CSV)
    parser.add_argument("--out-strict", type=Path, default=SDSS_CSV)
    parser.add_argument(
        "--no-color-cut",
        action="store_true",
        help="Skip u-r cut (default for Zurich audit).",
    )
    parser.add_argument("--timeout", type=int, default=900)
    args = parser.parse_args()

    print(f"[*] Querying SDSS DR{SDSS_DR} PhotoObj in COSMOS ACS box ...")
    tbl = SDSS.query_sql(SQL, data_release=SDSS_DR, timeout=args.timeout)
    raw = tbl.to_pandas()
    print(f"[*] SQL returned N={len(raw):,}")

    entire, funnel = apply_cosmos_sdss_entire(raw, color=not args.no_color_cut, disk_only=False)
    strict = apply_ba_strict(entire, "expAB_r", COSMOS_BA_MIN)
    funnel.append({"stage": f"expAB_r_gt_{COSMOS_BA_MIN:g}", "n_remaining": len(strict)})

    disk_entire, funnel_disk = apply_cosmos_sdss_entire(
        raw, color=not args.no_color_cut, disk_only=True
    )
    disk_strict = apply_ba_strict(disk_entire, "expAB_r", COSMOS_BA_MIN)
    funnel_disk.append({"stage": f"expAB_r_gt_{COSMOS_BA_MIN:g}", "n_remaining": len(disk_strict)})

    print(f"[*] Entire pool N={len(entire):,}; strict pool N={len(strict):,}")
    print(f"[*] Disk entire N={len(disk_entire):,}; disk strict N={len(disk_strict):,}")

    args.out_entire.parent.mkdir(parents=True, exist_ok=True)
    entire.to_csv(args.out_entire, index=False)
    strict.to_csv(args.out_strict, index=False)
    disk_entire.to_csv(SDSS_DISK_ENTIRE_CSV, index=False)
    disk_strict.to_csv(SDSS_DISK_CSV, index=False)
    print(f"[*] Wrote {args.out_entire}")
    print(f"[*] Wrote {args.out_strict}")
    print(f"[*] Wrote {SDSS_DISK_ENTIRE_CSV}")
    print(f"[*] Wrote {SDSS_DISK_CSV}")

    pd.DataFrame(funnel).to_csv(
        args.out_strict.with_name("cosmos_sdss_dr17_cut_funnel.csv"), index=False
    )
    pd.DataFrame(funnel_disk).to_csv(
        SDSS_DISK_CSV.with_name("cosmos_sdss_dr17_disk_cut_funnel.csv"), index=False
    )

    mag = pd.to_numeric(entire["modelMag_r"], errors="coerce")
    ba = pd.to_numeric(entire["expAB_r"], errors="coerce")
    print(
        f"[*] median modelMag_r={np.median(mag):.3f}, "
        f"median expAB_r={np.median(ba):.3f}"
    )


if __name__ == "__main__":
    main()
