"""Stage a GALFIT re-fit and re-run the full verification suite on it.

A re-fit is not just a new ``out.fits`` + panel. The verification checks
(RFF, isophotes, Fourier, sky perturbation, AstroPhot, visual, …) are re-run
against the new model and written under ``Re-fits/<FRB>/<label>/``, with the
production panel copied alongside for comparison.

Examples
--------
    python run_refit.py 20181112A --sky-adu 6.317e-5 --fix-n 1 --label n1_sky
    python run_refit.py 20181112A --sky-from-protocol --fix-n 1
"""

from __future__ import annotations

import argparse
import contextlib
import math
import os
import re
import shutil
import sys
import time
from typing import Any

VER_DIR = os.path.dirname(os.path.abspath(__file__))
if VER_DIR not in sys.path:
    sys.path.insert(0, VER_DIR)

import vercommon as vc  # noqa: E402
from checks import sky_perturb as sp  # noqa: E402
import run_verification as rv  # noqa: E402

REFITS_ROOT = os.path.join(VER_DIR, "Re-fits")

_COMP_HEADER = re.compile(r"^\s*#\s*Component number:\s*(\d+)", re.IGNORECASE)
_PARAM = re.compile(r"^(\s*)(\d+)\)\s*(.*?)(\s*#.*)?$")

# Sidecars needed so verification checks see the same supporting files as
# production (catalog sky MAD, PSFEx XML, constraints, …).
_SIDECARS = (
    "host_cutout.fits",
    "host_sigma.fits",
    "host_mask.fits",
    "proto_image.fits",
    "image.cat",
    "psfex.xml",
    "constraints.txt",
    "cutout_meta.json",
    "sky_fit_audit.json",
    "pipeline_summary.json",
    "reference_photometry.json",
    "zero_points.json",
)


def _edit_feedme_fixed(
    lines: list[str],
    comps: list[dict],
    *,
    sky_value: float | None,
    sky_comp: int,
    fix_n: float | None,
    reseed: bool = True,
    sky_seed_free: float | None = None,
) -> list[str]:
    """Reseed Sérsics at best-fit; optionally fix sky and/or host n.

    ``sky_value`` → sky held fixed (flag 0).
    ``sky_seed_free`` → sky free (flag 1) but reseeding at that ADU (used when
    only ``n`` is fixed so sky can still float from the production level).
    """
    out = list(lines)
    comp_id = 0
    comp_type = ""
    sersic_seen = 0
    params: dict | None = None

    for i, line in enumerate(lines):
        head = _COMP_HEADER.match(line)
        if head:
            comp_id = int(head.group(1))
            comp_type = ""
            params = None
            continue
        m = _PARAM.match(line)
        if not m:
            continue
        indent, num, body, comment = m.group(1), m.group(2), m.group(3), m.group(4) or ""
        if num == "0":
            comp_type = body.split()[0].lower() if body.split() else ""
            if comp_type == "sersic":
                params = comps[sersic_seen] if sersic_seen < len(comps) else None
                sersic_seen += 1
            continue

        if comp_type == "sky" and comp_id == sky_comp and num == "1":
            if sky_value is not None:
                out[i] = f"{indent}1) {sky_value:.8f} 0{comment}"
            elif sky_seed_free is not None:
                out[i] = f"{indent}1) {sky_seed_free:.8f} 1{comment}"
            continue

        if comp_type != "sersic" or params is None:
            continue

        # Host = first Sérsic. Fix n only on the host component.
        is_host = sersic_seen == 1
        if num == "5" and is_host and fix_n is not None:
            out[i] = f"{indent}5) {float(fix_n):.4f} 0{comment}"
            continue

        if not reseed:
            continue
        new = {
            "1": f"{params['xc']:.4f} {params['yc']:.4f} 1 1",
            "3": f"{params['mag']:.4f} 1",
            "4": f"{params['re']:.4f} 1",
            "5": f"{params['n']:.4f} 1",
            "9": f"{params['q']:.4f} 1",
            "10": f"{params['pa']:.4f} 1",
        }.get(num)
        # Don't overwrite a just-fixed n with the free reseed.
        if num == "5" and is_host and fix_n is not None:
            continue
        if new is not None:
            out[i] = f"{indent}{num}) {new}{comment}"
    return out


def _strip_param_constraints(src: str, dst: str, drop: set[str]) -> None:
    """Drop constraint rows for named params (e.g. {'n'} when n is fixed)."""
    try:
        with open(src, encoding="utf-8", errors="replace") as f:
            lines = f.read().splitlines()
    except OSError:
        lines = []
    kept = []
    for ln in lines:
        toks = ln.split()
        if len(toks) >= 2 and toks[1].lower() in drop:
            continue
        kept.append(ln)
    with open(dst, "w", newline="\n", encoding="utf-8") as f:
        f.write("\n".join(kept) + ("\n" if kept else ""))


def stage_refit(
    frb: str,
    label: str,
    *,
    sky_adu: float | None = None,
    fix_n: float | None = None,
) -> str:
    """Build ``Re-fits/<FRB>/<label>/`` ready for GALFIT + verification."""
    prod = vc.load_host(frb)
    wkdir = os.path.join(REFITS_ROOT, frb, label)
    if os.path.isdir(wkdir):
        shutil.rmtree(wkdir)
    os.makedirs(wkdir, exist_ok=True)

    for name in _SIDECARS:
        src = os.path.join(prod.dir, name)
        if os.path.isfile(src):
            shutil.copy2(src, os.path.join(wkdir, name))

    hdr = vc.parse_out_header(os.path.join(prod.dir, "out.fits"))
    comps = hdr["components"]
    sky_comp = int(hdr["sky"].get("comp", len(comps) + 1))

    feed = vc.parse_feedme(os.path.join(prod.dir, "galfit.feedme"))
    # If sky is left free, reseed it at the production best-fit (still flag 1).
    prod_sky = float(hdr["sky"].get("level", float("nan")))
    sky_seed = prod_sky if (sky_adu is None and math.isfinite(prod_sky)) else None
    edited = _edit_feedme_fixed(
        feed["lines"], comps,
        sky_value=sky_adu, sky_comp=sky_comp, fix_n=fix_n, reseed=True,
        sky_seed_free=sky_seed,
    )
    with open(os.path.join(wkdir, "galfit.feedme"), "w", newline="\n",
              encoding="utf-8") as f:
        f.write("\n".join(edited) + "\n")

    # Sky fixed → drop sky constraints; n fixed → drop n bounds row.
    drop: set[str] = set()
    if sky_adu is not None:
        # sky_perturb helper drops by component id; here rewrite constraints.
        csrc = os.path.join(prod.dir, "constraints.txt")
        try:
            with open(csrc, encoding="utf-8", errors="replace") as f:
                clines = f.read().splitlines()
        except OSError:
            clines = []
        kept = []
        for ln in clines:
            toks = ln.split()
            if not toks:
                kept.append(ln)
                continue
            if sky_adu is not None and toks[0] == str(sky_comp):
                continue
            if fix_n is not None and len(toks) >= 2 and toks[1].lower() == "n":
                continue
            kept.append(ln)
        with open(os.path.join(wkdir, "constraints.txt"), "w", newline="\n",
                  encoding="utf-8") as f:
            f.write("\n".join(kept) + ("\n" if kept else ""))
    elif fix_n is not None:
        _strip_param_constraints(
            os.path.join(prod.dir, "constraints.txt"),
            os.path.join(wkdir, "constraints.txt"),
            {"n"},
        )

    meta = {
        "frb": frb,
        "label": label,
        "sky_fixed_adu": sky_adu,
        "n_fixed": fix_n,
        "production_sky_adu": prod.sky_level,
        "production": {
            "q": prod.q, "n": prod.n, "re": prod.re,
            "mag": prod.mag, "pa": prod.pa,
        },
    }
    vc.write_json(os.path.join(wkdir, "refit_meta.json"), meta)
    return wkdir


def run_galfit_refit(wkdir: str) -> dict[str, Any]:
    run_galfit = sp._run_galfit()
    for stale in ("fit.log", "out.fits"):
        path = os.path.join(wkdir, stale)
        if os.path.isfile(path):
            os.remove(path)
    ok = False
    try:
        with open(os.path.join(wkdir, "galfit_stdout.log"), "w",
                  encoding="utf-8") as fh, contextlib.redirect_stdout(fh):
            ok = bool(run_galfit(wkdir))
    except Exception as exc:
        return {"status": f"launch_failed: {type(exc).__name__}: {exc}",
                "converged": False}
    fit = sp._read_result(wkdir)
    fit["converged"] = ok
    return fit


def production_panel_src(frb: str) -> str | None:
    """Canonical production panel — ``outputs/panels/<FRB>.png`` only."""
    path = os.path.join(vc.OUT_ROOT, "panels", f"{frb}.png")
    return path if os.path.isfile(path) else None


def copy_production_panel(frb: str, dest_path: str) -> str | None:
    """Byte-copy the production panel. Never regenerate it."""
    src = production_panel_src(frb)
    if src is None:
        return None
    os.makedirs(os.path.dirname(os.path.abspath(dest_path)) or ".", exist_ok=True)
    shutil.copy2(src, dest_path)
    return dest_path


def run_refit(
    frb: str,
    *,
    label: str,
    sky_adu: float | None = None,
    fix_n: float | None = None,
    checks: list[str] | None = None,
    force: bool = True,
) -> dict[str, Any]:
    """Stage → GALFIT → full verification suite for one re-fit label."""
    checks = checks or list(rv.CHECK_ORDER)
    wkdir = stage_refit(frb, label, sky_adu=sky_adu, fix_n=fix_n)
    fit = run_galfit_refit(wkdir)
    if fit.get("status") != "ok":
        return {"frb": frb, "label": label, "wkdir": wkdir, "fit": fit,
                "verification": None}

    t0 = time.time()
    ver = rv.run_host_dir(frb, wkdir, wkdir, checks, force=force)

    summary = {
        "frb": frb,
        "label": label,
        "wkdir": wkdir,
        "sky_fixed_adu": sky_adu,
        "n_fixed": fix_n,
        "fit": fit,
        "verification": ver,
        "verification_runtime_s": round(time.time() - t0, 3),
        "panel_png": os.path.join(wkdir, "panel.png"),
    }
    vc.write_json(os.path.join(wkdir, "refit_summary.json"), summary)
    return summary


def _sky_from_protocol(frb: str) -> float:
    path = os.path.join(REFITS_ROOT, frb, "sky_protocol.json")
    data = vc.read_json(path)
    sky = (data.get("consensus") or {}).get("sky_adu")
    if sky is None or not math.isfinite(float(sky)):
        raise SystemExit(f"no consensus sky in {path}; run sky_protocol.py first")
    return float(sky)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("frb")
    ap.add_argument("--label", default=None,
                    help="subdir under Re-fits/<FRB>/ (default derived from flags)")
    ap.add_argument("--sky-adu", type=float, default=None)
    ap.add_argument("--sky-from-protocol", action="store_true",
                    help="use Re-fits/<FRB>/sky_protocol.json consensus sky")
    ap.add_argument("--fix-n", type=float, default=None,
                    help="hold host Sérsic index fixed at this value")
    ap.add_argument("--checks", nargs="+", default=["all"],
                    choices=["all", *rv.CHECK_ORDER])
    ap.add_argument("--force", action="store_true", default=True)
    args = ap.parse_args(argv)

    sky = args.sky_adu
    if args.sky_from_protocol:
        sky = _sky_from_protocol(args.frb)

    parts = []
    if sky is not None:
        parts.append("sky")
    if args.fix_n is not None:
        parts.append(f"n{args.fix_n:g}")
    label = args.label or ("_".join(parts) if parts else "refit")

    checks = rv.CHECK_ORDER if "all" in args.checks else [
        c for c in rv.CHECK_ORDER if c in args.checks
    ]
    print(f"Re-fit {args.frb}  label={label}  sky={sky}  fix_n={args.fix_n}")
    print(f"  checks: {', '.join(checks)}")
    summary = run_refit(
        args.frb, label=label, sky_adu=sky, fix_n=args.fix_n,
        checks=checks, force=True,
    )
    fit = summary.get("fit") or {}
    print(f"  GALFIT: status={fit.get('status')}  q={fit.get('q')}  "
          f"n={fit.get('n')}  Re={fit.get('re')}  m={fit.get('mag')}")
    ver = summary.get("verification") or {}
    print(f"  verification: {ver.get('status')}")
    print(f"  panel: {summary.get('panel_png')}")
    return 0 if fit.get("status") == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
