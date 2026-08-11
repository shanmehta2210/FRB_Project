"""Four-panel re-fit grid for the 15 in-cut rejected hosts.

For each FRB under ``Re-fits/<FRB>/``:

* ``panel_production.png`` — **byte copy** of ``outputs/panels/<FRB>.png``
  (production: n free, sky free; never regenerated)
* ``panel_n1.png`` — GALFIT with ``n=1`` fixed, sky free + full verification
* ``panel_sky.png`` — GALFIT with protocol sky fixed, n free + full verification
* ``panel_n1_sky.png`` — GALFIT with ``n=1`` and protocol sky fixed + full verification

Full check products live under ``Re-fits/<FRB>/{n1,sky,n1_sky}/``.
Protocol sky is written to ``Re-fits/<FRB>/sky_protocol.json``.

Example
-------
    python run_reject_grid.py
    python run_reject_grid.py --frb 20181112A
"""

from __future__ import annotations

import argparse
import os
import sys
import time
import traceback

import pandas as pd

VER_DIR = os.path.dirname(os.path.abspath(__file__))
if VER_DIR not in sys.path:
    sys.path.insert(0, VER_DIR)

import vercommon as vc  # noqa: E402
import sky_protocol  # noqa: E402
import run_refit as rr  # noqa: E402

# Science-cut rejects from triage (mag ≤ 22, b/a > 0.2, confirmed=False).
IN53_REJECTS = [
    "20181112A",
    "20190102C",
    "20190523A",
    "20190711A",
    "20210410D",
    "20210807D",
    "20211127I",
    "20220501C",
    "20220509G",
    "20220725A",
    "20221219A",
    "20230526A",
    "20230626A",
    "20230930A",
    "20231220A",
]


def in53_rejects_from_csv() -> list[str]:
    """Recompute the 15 from host_confirmation + results (fallback to constant)."""
    hc_path = os.path.join(VER_DIR, "host_confirmation.csv")
    if not os.path.isfile(hc_path) or not os.path.isfile(vc.RESULTS_CSV):
        return list(IN53_REJECTS)
    hc = pd.read_csv(hc_path, dtype={"frb": str})
    res = vc.cohort("all64")
    m = res.merge(hc, on="frb", how="left")
    rej = m[m["in_53"] & (m["confirmed"] == False)]  # noqa: E712
    frbs = sorted(rej["frb"].astype(str).tolist())
    return frbs if frbs else list(IN53_REJECTS)


def _publish_panel(frb_root: str, label: str, wk_panel: str) -> str | None:
    if not os.path.isfile(wk_panel):
        return None
    dst = os.path.join(frb_root, f"panel_{label}.png")
    import shutil
    shutil.copy2(wk_panel, dst)
    return dst


def run_one_frb(frb: str, *, force_sky: bool = False) -> dict:
    frb_root = os.path.join(rr.REFITS_ROOT, frb)
    os.makedirs(frb_root, exist_ok=True)
    out: dict = {"frb": frb, "legs": {}}

    # 1) production panel — byte copy only
    prod_dst = os.path.join(frb_root, "panel_production.png")
    copied = rr.copy_production_panel(frb, prod_dst)
    out["panel_production"] = copied
    if copied is None:
        out["error"] = "missing outputs/panels/<FRB>.png"
        return out

    # 2) protocol sky
    sky_json = os.path.join(frb_root, "sky_protocol.json")
    if force_sky or not os.path.isfile(sky_json):
        sky_protocol.run_one(frb)
    sky = rr._sky_from_protocol(frb)
    out["sky_protocol_adu"] = sky

    # 3) three re-fits
    legs = (
        ("n1", {"fix_n": 1.0, "sky_adu": None}),
        ("sky", {"fix_n": None, "sky_adu": sky}),
        ("n1_sky", {"fix_n": 1.0, "sky_adu": sky}),
    )
    for label, kw in legs:
        t0 = time.time()
        try:
            summary = rr.run_refit(frb, label=label, **kw)
            pub = _publish_panel(frb_root, label, summary.get("panel_png") or "")
            out["legs"][label] = {
                "status": (summary.get("fit") or {}).get("status"),
                "fit": summary.get("fit"),
                "verification": (summary.get("verification") or {}).get("status"),
                "panel": pub,
                "runtime_s": round(time.time() - t0, 1),
            }
        except Exception as exc:
            out["legs"][label] = {
                "status": "error",
                "error": f"{type(exc).__name__}: {exc}",
                "traceback": traceback.format_exc(limit=6),
            }
    vc.write_json(os.path.join(frb_root, "reject_grid_summary.json"), out)
    return out


def main(argv: list[str] | None = None) -> int:
    # Windows consoles are often cp1252; keep progress prints ASCII-safe and live.
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")  # type: ignore[attr-defined]
    except Exception:
        pass
    try:
        sys.stdout.reconfigure(line_buffering=True)  # type: ignore[attr-defined]
    except Exception:
        pass

    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frb", nargs="+", default=None,
                    help="subset of FRBs (default: all 15 in-cut rejects)")
    ap.add_argument("--force-sky", action="store_true",
                    help="recompute sky_protocol.json even if present")
    args = ap.parse_args(argv)

    frbs = args.frb or in53_rejects_from_csv()
    print(f"Reject grid: {len(frbs)} hosts × (production copy + n1 + sky + n1_sky)")
    t0 = time.time()
    for i, frb in enumerate(frbs, 1):
        print(f"\n[{i}/{len(frbs)}] {frb}")
        res = run_one_frb(frb, force_sky=args.force_sky)
        if res.get("error"):
            print(f"  ERROR: {res['error']}")
            continue
        sky = res.get("sky_protocol_adu")
        print(f"  sky_protocol = {sky:.6g} ADU" if sky is not None else "  sky_protocol = ?")
        print(f"  panel_production = {res.get('panel_production')}")
        for label, leg in (res.get("legs") or {}).items():
            st = leg.get("status")
            fit = leg.get("fit") or {}
            print(f"  {label:8s}  fit={st}  q={fit.get('q')}  n={fit.get('n')}  "
                  f"Re={fit.get('re')}  ({leg.get('runtime_s')}s)")
            if leg.get("error"):
                print(f"           {leg['error']}")
    print(f"\nDone in {(time.time() - t0) / 60:.1f} min")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
