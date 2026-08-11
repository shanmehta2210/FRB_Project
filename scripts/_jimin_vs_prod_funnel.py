"""
Contrast production SDSS v2 cuts vs Jimin cuts in the Jimin 4x4 deg box.

Uses live DR16 COUNTs (not the HTM-random production CSV, which only sparsely
samples any given patch).
"""
from __future__ import annotations

import csv
import time
from pathlib import Path

from astroquery.sdss import SDSS

DR = 16
BOX = "p.ra BETWEEN 148.0 AND 152.0 AND p.dec BETWEEN 0.0 AND 4.0"
AREA = 16.0
OUT = (
    Path(__file__).resolve().parents[1]
    / "plots"
    / "plots_null"
    / "v2"
    / "sdss_audit"
    / "Jimin"
    / "prod_vs_jimin_funnel.csv"
)


def count(sql: str, label: str) -> int:
    print(f"\n=== {label} ===", flush=True)
    for attempt in range(1, 4):
        try:
            tbl = SDSS.query_sql(sql, data_release=DR, timeout=900)
            n = int(tbl[0][0]) if tbl is not None else -1
            print(f"N={n:,}  dens={n/AREA:,.1f}/deg2", flush=True)
            return n
        except Exception as exc:
            print(f"fail {attempt}: {exc}", flush=True)
            time.sleep(2 * attempt)
    return -1


# Shared minimal galaxy base matching production fetch + mag window
# Production fetch: type=3, clean=1, mode=1 (no Photoz, no lnLStar, no score)
# Mag: model r in [12,21] to match Jimin window for fair contrast
BASE = f"""
FROM PhotoObj AS p
WHERE {BOX}
AND p.mode=1 AND p.clean=1 AND p.type=3
AND p.r BETWEEN 12 AND 21
"""

STEPS = [
    (
        "P0_prod_style_base_type3_mag12_21",
        f"SELECT COUNT(*) {BASE}",
        "Production-like: type=3, mode, clean, model r [12,21]",
    ),
    (
        "P0b_type_r3_instead",
        f"""SELECT COUNT(*) FROM PhotoObj AS p
        WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
        AND p.r BETWEEN 12 AND 21""",
        "Swap type->type_r=3 (Jimin)",
    ),
    (
        "J1_plus_lnLStar",
        f"""SELECT COUNT(*) FROM PhotoObj AS p
        WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
        AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10""",
        "Jimin add: lnLStar_r < -10",
    ),
    (
        "J2_plus_Photoz_nnAvgZ",
        f"""SELECT COUNT(*) FROM PhotoObj AS p
        JOIN Photoz AS pz ON pz.objid=p.objid
        WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
        AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10 AND pz.nnAvgZ > 0""",
        "Jimin add: Photoz nnAvgZ > 0",
    ),
    (
        "J3_plus_score",
        f"""SELECT COUNT(*) FROM PhotoObj AS p
        JOIN Photoz AS pz ON pz.objid=p.objid
        WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
        AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10 AND pz.nnAvgZ > 0
        AND p.score > 0.8""",
        "Jimin add: score > 0.8",
    ),
    (
        "J4_V2_plus_lnL_exp",
        f"""SELECT COUNT(*) FROM PhotoObj AS p
        JOIN Photoz AS pz ON pz.objid=p.objid
        WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
        AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10 AND pz.nnAvgZ > 0
        AND p.score > 0.8 AND p.lnLDeV_r < p.lnLExp_r""",
        "Jimin V2: + lnLExp wins",
    ),
    (
        "J5_V1_plus_fracDeV0",
        f"""SELECT COUNT(*) FROM PhotoObj AS p
        JOIN Photoz AS pz ON pz.objid=p.objid
        WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
        AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10 AND pz.nnAvgZ > 0
        AND p.score > 0.8 AND p.fracDeV_r=0 AND p.lnLDeV_r < p.lnLExp_r""",
        "Jimin V1: + fracDeV=0",
    ),
    # Production post-hoc style from P0 (no Jimin extras)
    (
        "P1_prod_plus_lnL_exp",
        f"""SELECT COUNT(*) {BASE} AND p.lnLDeV_r < p.lnLExp_r""",
        "Production CDF morph: lnLExp wins (no Jimin extras)",
    ),
    (
        "P2_prod_plus_lnL_and_ur_lt_2p3",
        f"""SELECT COUNT(*) {BASE}
        AND p.lnLDeV_r < p.lnLExp_r
        AND (p.u - p.r) < 2.3""",
        "Production CDF: lnL + u-r < 2.3",
    ),
    (
        "P3_prod_lnL_ur_expAB_gt_0p2",
        f"""SELECT COUNT(*) {BASE}
        AND p.lnLDeV_r < p.lnLExp_r
        AND (p.u - p.r) < 2.3
        AND p.expAB_r > 0.2""",
        "Production strict null @ mag<=21 in box",
    ),
    (
        "J4b_V2_plus_expAB",
        f"""SELECT COUNT(*) FROM PhotoObj AS p
        JOIN Photoz AS pz ON pz.objid=p.objid
        WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
        AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10 AND pz.nnAvgZ > 0
        AND p.score > 0.8 AND p.lnLDeV_r < p.lnLExp_r AND p.expAB_r > 0.2""",
        "Jimin V2 strict ba",
    ),
    (
        "J5b_V1_plus_expAB",
        f"""SELECT COUNT(*) FROM PhotoObj AS p
        JOIN Photoz AS pz ON pz.objid=p.objid
        WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
        AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10 AND pz.nnAvgZ > 0
        AND p.score > 0.8 AND p.fracDeV_r=0 AND p.lnLDeV_r < p.lnLExp_r
        AND p.expAB_r > 0.2""",
        "Jimin V1 strict ba",
    ),
    # Value tests: apply Jimin extras one-at-a-time on prod base+lnL
    (
        "V_lnLStar_alone_on_P1",
        f"""SELECT COUNT(*) {BASE}
        AND p.lnLDeV_r < p.lnLExp_r AND p.lnLStar_r < -10""",
        "Value: lnLStar on prod+lnL",
    ),
    (
        "V_score_Photoz_on_P1",
        f"""SELECT COUNT(*) FROM PhotoObj AS p
        JOIN Photoz AS pz ON pz.objid=p.objid
        WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type=3
        AND p.r BETWEEN 12 AND 21 AND p.lnLDeV_r < p.lnLExp_r
        AND pz.nnAvgZ > 0 AND p.score > 0.8""",
        "Value: Photoz+score on prod+lnL (no lnLStar)",
    ),
    (
        "V_fracDeV0_on_P1",
        f"""SELECT COUNT(*) {BASE}
        AND p.lnLDeV_r < p.lnLExp_r AND p.fracDeV_r = 0""",
        "Value: fracDeV=0 on prod+lnL",
    ),
    (
        "V_ur_on_JiminV2",
        f"""SELECT COUNT(*) FROM PhotoObj AS p
        JOIN Photoz AS pz ON pz.objid=p.objid
        WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
        AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10 AND pz.nnAvgZ > 0
        AND p.score > 0.8 AND p.lnLDeV_r < p.lnLExp_r
        AND (p.u - p.r) < 2.3""",
        "What we add that Jimin lacks: u-r<2.3 on Jimin V2",
    ),
]


def main() -> None:
    rows = []
    for key, sql, note in STEPS:
        n = count(sql, f"{key}: {note}")
        rows.append(
            {
                "stage": key,
                "note": note,
                "n": n,
                "per_deg2": round(n / AREA, 2) if n >= 0 else "",
            }
        )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["stage", "note", "n", "per_deg2"])
        w.writeheader()
        w.writerows(rows)
    print(f"\n[*] Wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
