"""Orchestrator for the GALFIT fit verification suite.

Every check is idempotent and writes its own JSON under
``outputs/per_host/<FRB>/``, so a failure in one check never blocks the rest and
re-runs are cheap. See ``FIT_VERIFICATION_CHECKS.md`` for what each check does.

Examples
--------
    python run_verification.py --checks chi2 rff fourier
    python run_verification.py --checks all --jobs 4
    python run_verification.py --frb 20240210A --checks fourier --force
    python run_verification.py --aggregate-only
"""

from __future__ import annotations

import argparse
import importlib
import os
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed

VER_DIR = os.path.dirname(os.path.abspath(__file__))
if VER_DIR not in sys.path:
    sys.path.insert(0, VER_DIR)

import vercommon as vc  # noqa: E402

CHECK_MODULES = {
    "chi2": "checks.chi2_local",
    "rff": "checks.rff",
    "fourier": "checks.fourier",
    "psf": "checks.psf_leakage",
    "mag": "checks.mag_leakage",
    "isophote": "checks.isophotes",
    "sky": "checks.sky_perturb",
    "astrophot": "checks.astrophot_refit",
    "visual": "checks.visual",
}

# visual consumes the products of the others, so it must come last.
CHECK_ORDER = ["chi2", "rff", "fourier", "psf", "mag", "isophote", "sky",
               "astrophot", "visual"]


def _load_check(name: str):
    return importlib.import_module(CHECK_MODULES[name])


def run_host_dir(
    frb: str,
    host_dir: str,
    outdir: str,
    checks: list[str],
    force: bool,
) -> dict:
    """Run checks for a host loaded from ``host_dir``, writing into ``outdir``.

    Used for production (``host_dir = Output/<FRB>_all``) and for Re-fits
    (``host_dir = outdir = Re-fits/<FRB>/<label>``). Never raises.
    """
    status: dict[str, str] = {}
    host = None
    for name in checks:
        json_path = os.path.join(outdir, f"{name}.json")
        if os.path.isfile(json_path) and not force:
            status[name] = "cached"
            continue
        try:
            if host is None:
                host = vc.load_host_from_dir(host_dir, frb=frb)
            module = _load_check(name)
            t0 = time.time()
            payload = module.run(host, outdir)
            payload["_check"] = name
            payload["_frb"] = frb
            payload["_runtime_s"] = round(time.time() - t0, 3)
            payload.setdefault("status", "ok")
            vc.write_json(json_path, payload)
            status[name] = payload["status"]
        except Exception as exc:  # a broken check must not kill the host
            vc.write_json(
                json_path,
                {
                    "_check": name,
                    "_frb": frb,
                    "status": "error",
                    "error": f"{type(exc).__name__}: {exc}",
                    "traceback": traceback.format_exc(limit=8),
                },
            )
            status[name] = "error"
    return {"frb": frb, "status": status,
            "notes": list(host.notes) if host is not None else [],
            "host_dir": host_dir, "outdir": outdir}


def run_host(frb: str, checks: list[str], force: bool) -> dict:
    """Run the requested checks for one production host. Never raises."""
    return run_host_dir(frb, vc.host_dir(frb), vc.per_host_dir(frb), checks, force)


def _worker(args: tuple[str, list[str], bool]) -> dict:
    return run_host(*args)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--checks", nargs="+", default=["all"],
                    choices=["all", *CHECK_ORDER])
    ap.add_argument("--frb", nargs="+", default=None, help="restrict to these FRBs")
    ap.add_argument("--cohort", choices=["all64", "53"], default="all64")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--force", action="store_true",
                    help="recompute checks even when a JSON already exists")
    ap.add_argument("--no-aggregate", action="store_true")
    ap.add_argument("--aggregate-only", action="store_true")
    args = ap.parse_args(argv)

    checks = CHECK_ORDER if "all" in args.checks else [
        c for c in CHECK_ORDER if c in args.checks
    ]

    df = vc.cohort(args.cohort)
    frbs = list(df["frb"].astype(str))
    if args.frb:
        wanted = set(args.frb)
        frbs = [f for f in frbs if f in wanted]
        missing = wanted - set(frbs)
        if missing:
            print(f"[warn] not in cohort: {', '.join(sorted(missing))}")

    os.makedirs(vc.TABLES_ROOT, exist_ok=True)
    os.makedirs(vc.LOGS_ROOT, exist_ok=True)

    if not args.aggregate_only:
        print(f"Verification suite: {len(frbs)} hosts x {len(checks)} checks "
              f"({', '.join(checks)}), jobs={args.jobs}")
        t0 = time.time()
        results: list[dict] = []
        if args.jobs > 1:
            payload = [(f, checks, args.force) for f in frbs]
            with ProcessPoolExecutor(max_workers=args.jobs) as pool:
                futures = {pool.submit(_worker, p): p[0] for p in payload}
                for i, fut in enumerate(as_completed(futures), 1):
                    res = fut.result()
                    results.append(res)
                    _report(res, i, len(frbs))
        else:
            for i, frb in enumerate(frbs, 1):
                res = run_host(frb, checks, args.force)
                results.append(res)
                _report(res, i, len(frbs))
        print(f"Done in {time.time() - t0:.1f}s")

        errors = [
            (r["frb"], k) for r in results for k, v in r["status"].items() if v == "error"
        ]
        if errors:
            print(f"\n{len(errors)} check errors:")
            for frb, check in sorted(errors):
                print(f"  {frb:12s} {check}")

    if not args.no_aggregate:
        import aggregate

        aggregate.main(cohort=args.cohort)
    return 0


def _report(res: dict, i: int, total: int) -> None:
    bad = [k for k, v in res["status"].items() if v == "error"]
    skipped = [k for k, v in res["status"].items() if v == "skipped"]
    tag = "ERROR " + ",".join(bad) if bad else ("skip " + ",".join(skipped) if skipped else "ok")
    print(f"[{i:3d}/{total}] {res['frb']:12s} {tag}")
    for note in res.get("notes", []):
        print(f"           note: {note}")


if __name__ == "__main__":
    raise SystemExit(main())
