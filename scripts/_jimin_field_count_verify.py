"""
Careful DR16 COUNT funnel for the Jimin 4x4 deg box.

Area = 4 deg x 4 deg = 16 deg^2.
Mag: both inclusive BETWEEN 12 AND 21 (historical) and r < 21 / r <= 21 variants.
"""
from __future__ import annotations

import csv
import time
from pathlib import Path

from astroquery.sdss import SDSS

DR = 16
BOX = "p.ra BETWEEN 148.0 AND 152.0 AND p.dec BETWEEN 0.0 AND 4.0"
AREA_DEG2 = 16.0
OUT = (
    Path(__file__).resolve().parents[1]
    / "plots"
    / "plots_null"
    / "v2"
    / "sdss_audit"
    / "Jimin"
    / "field_count_verify.csv"
)


def count(sql: str, label: str) -> int:
    print(f"\n=== {label} ===", flush=True)
    for attempt in range(1, 4):
        try:
            tbl = SDSS.query_sql(sql, data_release=DR, timeout=900)
            n = int(tbl[0][0]) if tbl is not None else -1
            dens = n / AREA_DEG2 if n >= 0 else float("nan")
            print(f"N = {n:,}   density = {dens:,.1f} / deg^2", flush=True)
            return n
        except Exception as exc:
            print(f"  fail {attempt}: {exc}", flush=True)
            time.sleep(2.0 * attempt)
    return -1


# Ordered funnel + ablations. Labels must stay unique.
STEPS: list[tuple[str, str]] = [
    (
        "01_box_all_PhotoObj",
        f"SELECT COUNT(*) FROM PhotoObj AS p WHERE {BOX}",
    ),
    (
        "02_box_mode1",
        f"SELECT COUNT(*) FROM PhotoObj AS p WHERE {BOX} AND p.mode=1",
    ),
    (
        "03_box_mode1_clean1",
        f"SELECT COUNT(*) FROM PhotoObj AS p WHERE {BOX} AND p.mode=1 AND p.clean=1",
    ),
    (
        "04_mode_clean_type3",
        f"SELECT COUNT(*) FROM PhotoObj AS p WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type=3",
    ),
    (
        "05_mode_clean_type_r3",
        f"SELECT COUNT(*) FROM PhotoObj AS p WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3",
    ),
    # Mag variants on type_r=3 + mode + clean
    (
        "06_type_r3_r_BETWEEN_12_21",
        f"""SELECT COUNT(*) FROM PhotoObj AS p WHERE {BOX}
        AND p.mode=1 AND p.clean=1 AND p.type_r=3 AND p.r BETWEEN 12 AND 21""",
    ),
    (
        "07_type_r3_r_lt_21",
        f"""SELECT COUNT(*) FROM PhotoObj AS p WHERE {BOX}
        AND p.mode=1 AND p.clean=1 AND p.type_r=3 AND p.r < 21""",
    ),
    (
        "08_type_r3_r_le_21",
        f"""SELECT COUNT(*) FROM PhotoObj AS p WHERE {BOX}
        AND p.mode=1 AND p.clean=1 AND p.type_r=3 AND p.r <= 21""",
    ),
    (
        "09_type_r3_r_BETWEEN_12_21_no_bright_floor",
        f"""SELECT COUNT(*) FROM PhotoObj AS p WHERE {BOX}
        AND p.mode=1 AND p.clean=1 AND p.type_r=3 AND p.r <= 21 AND p.r >= 0""",
    ),
    (
        "10_type3_r_BETWEEN_12_21",
        f"""SELECT COUNT(*) FROM PhotoObj AS p WHERE {BOX}
        AND p.mode=1 AND p.clean=1 AND p.type=3 AND p.r BETWEEN 12 AND 21""",
    ),
    (
        "11_type_r3_petro_BETWEEN_12_21",
        f"""SELECT COUNT(*) FROM PhotoObj AS p WHERE {BOX}
        AND p.mode=1 AND p.clean=1 AND p.type_r=3 AND p.petroMag_r BETWEEN 12 AND 21""",
    ),
    # Continue historical funnel from 06
    (
        "12_plus_lnLStar_lt_m10",
        f"""SELECT COUNT(*) FROM PhotoObj AS p WHERE {BOX}
        AND p.mode=1 AND p.clean=1 AND p.type_r=3 AND p.r BETWEEN 12 AND 21
        AND p.lnLStar_r < -10""",
    ),
    (
        "13_plus_Photoz_nnAvgZ_gt0",
        f"""SELECT COUNT(*) FROM PhotoObj AS p
        JOIN Photoz AS pz ON pz.objid=p.objid
        WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
        AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10 AND pz.nnAvgZ > 0""",
    ),
    (
        "14_plus_score_gt_0p8_BASE",
        f"""SELECT COUNT(*) FROM PhotoObj AS p
        JOIN Photoz AS pz ON pz.objid=p.objid
        WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
        AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10 AND pz.nnAvgZ > 0
        AND p.score > 0.8""",
    ),
    (
        "15_BASE_plus_lnL_V2",
        f"""SELECT COUNT(*) FROM PhotoObj AS p
        JOIN Photoz AS pz ON pz.objid=p.objid
        WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
        AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10 AND pz.nnAvgZ > 0
        AND p.score > 0.8 AND p.lnLDeV_r < p.lnLExp_r""",
    ),
    (
        "16_BASE_plus_fracDeV0_and_lnL_V1",
        f"""SELECT COUNT(*) FROM PhotoObj AS p
        JOIN Photoz AS pz ON pz.objid=p.objid
        WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
        AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10 AND pz.nnAvgZ > 0
        AND p.score > 0.8 AND p.fracDeV_r = 0 AND p.lnLDeV_r < p.lnLExp_r""",
    ),
    (
        "17_BASE_plus_fracDeV0_only",
        f"""SELECT COUNT(*) FROM PhotoObj AS p
        JOIN Photoz AS pz ON pz.objid=p.objid
        WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
        AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10 AND pz.nnAvgZ > 0
        AND p.score > 0.8 AND p.fracDeV_r = 0""",
    ),
    (
        "18_V2_plus_expAB_gt_0p2",
        f"""SELECT COUNT(*) FROM PhotoObj AS p
        JOIN Photoz AS pz ON pz.objid=p.objid
        WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
        AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10 AND pz.nnAvgZ > 0
        AND p.score > 0.8 AND p.lnLDeV_r < p.lnLExp_r AND p.expAB_r > 0.2""",
    ),
    (
        "19_V1_plus_expAB_gt_0p2",
        f"""SELECT COUNT(*) FROM PhotoObj AS p
        JOIN Photoz AS pz ON pz.objid=p.objid
        WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
        AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10 AND pz.nnAvgZ > 0
        AND p.score > 0.8 AND p.fracDeV_r = 0 AND p.lnLDeV_r < p.lnLExp_r
        AND p.expAB_r > 0.2""",
    ),
    # Sanity: stars + galaxies mix without type cut, mag limited
    (
        "20_mode_clean_r_BETWEEN_12_21_any_type",
        f"""SELECT COUNT(*) FROM PhotoObj AS p WHERE {BOX}
        AND p.mode=1 AND p.clean=1 AND p.r BETWEEN 12 AND 21""",
    ),
    (
        "21_mode_clean_type3_r_lt_21",
        f"""SELECT COUNT(*) FROM PhotoObj AS p WHERE {BOX}
        AND p.mode=1 AND p.clean=1 AND p.type=3 AND p.r < 21""",
    ),
]


def main() -> None:
    rows: list[dict] = []
    prev: int | None = None
    print(f"Field: RA[148,152] x Dec[0,4] = {AREA_DEG2:g} deg^2  DR{DR}", flush=True)
    for label, sql in STEPS:
        n = count(sql, label)
        dropped = "" if prev is None else str(prev - n)
        rows.append(
            {
                "stage": label,
                "n": n,
                "per_deg2": round(n / AREA_DEG2, 2) if n >= 0 else "",
                "delta_from_prev": dropped,
            }
        )
        # Only chain deltas along the main funnel stages that share parentage loosely
        if label.startswith(("01", "02", "03", "05", "06", "12", "13", "14", "15", "16")):
            prev = n
        elif label.startswith(("18", "19")):
            prev = n

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["stage", "n", "per_deg2", "delta_from_prev"])
        w.writeheader()
        w.writerows(rows)
    print(f"\n[*] Wrote {OUT}", flush=True)


if __name__ == "__main__":
    main()
