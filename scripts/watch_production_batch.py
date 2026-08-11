"""Live per-FRB validator for a ``run_all_frbs.py`` production batch.

Polls ``pipeline_scripts/Output/<FRB>_all/`` for the cohort under test and
validates each host as soon as its ``pipeline_summary.json`` is refreshed, so a
broken run is reported immediately instead of after the whole batch.

Emitted sentinels (stable, greppable):

    [watch] PASS   <frb>  <detail>
    [watch] WARN   <frb>  <reason>            # tolerated: Phase 2 ZP failed, fit is fine
    [watch] BROKEN <frb>  <reason>; <reason>
    [watch] DONE   pass=<n> warn=<n> broken=<n> pending=<n>

A Phase 2 failure alone is a **warning**, not a break: with a sparse calibration
field (`Too few calibration matches`) the pipeline still fits in
localization-host mode and `reference_photometry.py` falls back to the survey
magnitude (`mag_final_source = reference_*`). Only the raw GALFIT `mag` is
untrusted. Anything touching Phase 1 / 3a / 3b is a genuine break.

Usage:
    python scripts/watch_production_batch.py --list-file production_confirmed_lit_hosts.csv \
        --since-epoch <unix_ts> [--timeout-min 240]
"""
import argparse
import json
import time
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
OUTPUT_DIR = REPO / "pipeline_scripts" / "Output"

EXPECTED_POLICY = "re_separation"
REQUIRED_FILES = (
    "pipeline_summary.json",
    "cutout_meta.json",
    "host_components.csv",
    "host_cutout.fits",
    "host_sigma.fits",
    "fit.log",
    "galfit.feedme",
)


TOLERATED_PHASES = {"phase2_photometry_astropath"}


def validate(frb: str) -> tuple[list[str], list[str]]:
    """Return (problems, warnings); both empty means the host is healthy."""
    out = OUTPUT_DIR / f"{frb}_all"
    problems: list[str] = []
    warnings: list[str] = []
    if not out.is_dir():
        return [f"missing {out.name}/"], warnings

    for name in REQUIRED_FILES:
        if not (out / name).is_file():
            problems.append(f"missing {name}")

    summary_path = out / "pipeline_summary.json"
    if summary_path.is_file():
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            problems.append(f"unreadable pipeline_summary.json ({exc})")
            summary = {}
        for phase, rc in (summary.get("phase_exit_codes") or {}).items():
            if rc in (0, -1):
                continue
            if phase in TOLERATED_PHASES:
                ref = summary.get("reference_photometry") or {}
                if ref.get("status") == "ok" and ref.get("mag") is not None:
                    fallback = f"survey mag {ref['survey']} r={ref['mag']}"
                else:
                    fallback = "NO survey fallback — host has no usable magnitude"
                warnings.append(f"{phase} rc={rc}; {fallback}")
            else:
                problems.append(f"{phase} rc={rc}")
        host = summary.get("galfit_host") or {}
        chi2nu = host.get("chi2nu")
        if chi2nu is None:
            problems.append("no galfit_host.chi2nu")

    meta = {}
    meta_path = out / "cutout_meta.json"
    if meta_path.is_file():
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as exc:
            problems.append(f"unreadable cutout_meta.json ({exc})")
        policy = meta.get("neighbor_policy")
        if policy != EXPECTED_POLICY:
            problems.append(f"neighbor_policy={policy!r} (expected {EXPECTED_POLICY!r})")
        if meta.get("host_pad") is None or meta.get("re_sep_factor") is None:
            problems.append("cutout_meta missing host_pad / re_sep_factor")

    comps_path = out / "host_components.csv"
    if comps_path.is_file():
        try:
            comps = pd.read_csv(comps_path)
        except Exception as exc:  # noqa: BLE001 - report any parse failure
            problems.append(f"unreadable host_components.csv ({exc})")
            comps = None
        if comps is not None:
            if comps.empty:
                problems.append("host_components.csv empty")
            elif meta.get("host_number") is not None:
                first = int(comps.iloc[0]["NUMBER"])
                if first != int(meta["host_number"]):
                    problems.append(
                        f"host is not component 1 (first NUMBER={first}, "
                        f"host_number={meta['host_number']})"
                    )
    return problems, warnings


def detail(frb: str) -> str:
    out = OUTPUT_DIR / f"{frb}_all"
    try:
        summary = json.loads((out / "pipeline_summary.json").read_text(encoding="utf-8"))
        meta = json.loads((out / "cutout_meta.json").read_text(encoding="utf-8"))
    except Exception:  # noqa: BLE001 - detail line is best-effort
        return ""
    host = summary.get("galfit_host") or {}
    return (
        f"chi2nu={host.get('chi2nu')} b/a={host.get('b_a')} "
        f"fit={meta.get('n_fit_components')} mask={meta.get('n_mask_objects')} "
        f"roi={meta.get('cutout_bounds')}"
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--list-file", required=True, help="CSV with an 'frb' column")
    ap.add_argument("--since-epoch", type=float, required=True,
                    help="Only validate hosts whose pipeline_summary.json is newer than this")
    ap.add_argument("--timeout-min", type=float, default=300.0)
    ap.add_argument("--poll-sec", type=float, default=20.0)
    args = ap.parse_args()

    frbs = [str(f) for f in pd.read_csv(args.list_file, dtype={"frb": str})["frb"]]
    print(f"[watch] tracking {len(frbs)} FRB(s) from {Path(args.list_file).name}", flush=True)

    seen: dict[str, str] = {}
    deadline = time.time() + args.timeout_min * 60.0

    while time.time() < deadline and len(seen) < len(frbs):
        for frb in frbs:
            if frb in seen:
                continue
            summary = OUTPUT_DIR / f"{frb}_all" / "pipeline_summary.json"
            if not summary.is_file() or summary.stat().st_mtime < args.since_epoch:
                continue
            time.sleep(1.0)  # let the writer finish flushing
            problems, warnings = validate(frb)
            if problems:
                seen[frb] = "BROKEN"
                print(f"[watch] BROKEN {frb}  " + "; ".join(problems + warnings), flush=True)
            elif warnings:
                seen[frb] = "WARN"
                print(f"[watch] WARN   {frb}  " + "; ".join(warnings)
                      + f"  |  {detail(frb)}", flush=True)
            else:
                seen[frb] = "PASS"
                print(f"[watch] PASS   {frb}  {detail(frb)}", flush=True)
        time.sleep(args.poll_sec)

    counts = {k: sum(1 for v in seen.values() if v == k) for k in ("PASS", "WARN", "BROKEN")}
    pending = [f for f in frbs if f not in seen]
    if pending:
        print(f"[watch] pending (never refreshed): {', '.join(pending)}", flush=True)
    print(
        f"[watch] DONE pass={counts['PASS']} warn={counts['WARN']} "
        f"broken={counts['BROKEN']} pending={len(pending)}",
        flush=True,
    )
    return 1 if (counts["BROKEN"] or pending) else 0


if __name__ == "__main__":
    raise SystemExit(main())
