"""Diagnostic counts for Jimin-box SDSS DR16 selection (no TOP)."""
from __future__ import annotations

import time
from astroquery.sdss import SDSS

DR = 16
BOX = """
p.ra BETWEEN 148.0 AND 152.0
AND p.dec BETWEEN 0.0 AND 4.0
"""


def run(sql: str, label: str) -> int:
    print(f"\n=== {label} ===", flush=True)
    print(sql.strip()[:200], "...", flush=True)
    for attempt in range(3):
        try:
            tbl = SDSS.query_sql(sql, data_release=DR, timeout=900)
            n = int(tbl[0][0]) if tbl is not None else -1
            print(f"N = {n:,}", flush=True)
            return n
        except Exception as e:
            print(f"fail {attempt+1}: {e}", flush=True)
            time.sleep(2 * (attempt + 1))
    return -1


queries = [
    ("A: box only, mode=1 clean=1 type=3", f"""
SELECT COUNT(*) FROM PhotoObj AS p
WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type=3
"""),
    ("B: box + type_r=3 (Jimin)", f"""
SELECT COUNT(*) FROM PhotoObj AS p
WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
"""),
    ("C: B + r BETWEEN 12 AND 21", f"""
SELECT COUNT(*) FROM PhotoObj AS p
WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
AND p.r BETWEEN 12 AND 21
"""),
    ("D: C + lnLStar < -10", f"""
SELECT COUNT(*) FROM PhotoObj AS p
WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10
"""),
    ("E: D + INNER JOIN Photoz nnAvgZ>0", f"""
SELECT COUNT(*) FROM PhotoObj AS p
JOIN Photoz AS pz ON pz.objid=p.objid
WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10 AND pz.nnAvgZ > 0
"""),
    ("F: E + score>0.8 (Jimin base, no morph)", f"""
SELECT COUNT(*) FROM PhotoObj AS p
JOIN Photoz AS pz ON pz.objid=p.objid
WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10 AND pz.nnAvgZ > 0
AND p.score > 0.8
"""),
    ("G: F + lnLExp wins (V2 full, no TOP)", f"""
SELECT COUNT(*) FROM PhotoObj AS p
JOIN Photoz AS pz ON pz.objid=p.objid
WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10 AND pz.nnAvgZ > 0
AND p.score > 0.8 AND p.lnLDeV_r < p.lnLExp_r
"""),
    ("H: F + fracDeV=0 AND lnL (V1 full, no TOP)", f"""
SELECT COUNT(*) FROM PhotoObj AS p
JOIN Photoz AS pz ON pz.objid=p.objid
WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10 AND pz.nnAvgZ > 0
AND p.score > 0.8 AND p.fracDeV_r = 0 AND p.lnLDeV_r < p.lnLExp_r
"""),
    ("I: F + fracDeV < 0.01 (soft) AND lnL", f"""
SELECT COUNT(*) FROM PhotoObj AS p
JOIN Photoz AS pz ON pz.objid=p.objid
WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10 AND pz.nnAvgZ > 0
AND p.score > 0.8 AND p.fracDeV_r < 0.01 AND p.lnLDeV_r < p.lnLExp_r
"""),
    ("J: V2 morph WITHOUT Photoz join", f"""
SELECT COUNT(*) FROM PhotoObj AS p
WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10
AND p.score > 0.8 AND p.lnLDeV_r < p.lnLExp_r
"""),
    ("K: V2 morph WITHOUT score", f"""
SELECT COUNT(*) FROM PhotoObj AS p
JOIN Photoz AS pz ON pz.objid=p.objid
WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10 AND pz.nnAvgZ > 0
AND p.lnLDeV_r < p.lnLExp_r
"""),
    ("L: V1 WITHOUT Photoz (fracDeV=0 + lnL)", f"""
SELECT COUNT(*) FROM PhotoObj AS p
WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type_r=3
AND p.r BETWEEN 12 AND 21 AND p.lnLStar_r < -10
AND p.score > 0.8 AND p.fracDeV_r = 0 AND p.lnLDeV_r < p.lnLExp_r
"""),
    ("M: box + type=3 + r 12-21 + lnL only (minimal)", f"""
SELECT COUNT(*) FROM PhotoObj AS p
WHERE {BOX} AND p.mode=1 AND p.clean=1 AND p.type=3
AND p.r BETWEEN 12 AND 21 AND p.lnLDeV_r < p.lnLExp_r
"""),
]

if __name__ == "__main__":
    for label, sql in queries:
        run(sql, label)
