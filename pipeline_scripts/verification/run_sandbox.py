"""Manual GALFIT sandbox for weird hosts.

Edit ``Re-fits/<FRB>/sandbox/galfit.feedme`` by hand, then::

    python run_sandbox.py 20220509G
    python run_sandbox.py 20230930A

Runs GALFIT in that sandbox and writes ``panel.png`` (data | model | resid).

Keeps every ``galfit.NN`` restart file and any renamed archives (e.g.
``cent_fix``) — only refreshes ``fit.log`` / ``out.fits`` / ``panel.png``.

Setup (once)::

    python run_sandbox.py --init

Panel stretch (see ``make_panel``):
  * data/model — grayscale asinh of the [1, 99] percentile window
  * resid — linear ``RdBu_r``, clipped to ±98th percentile of |resid|
"""

from __future__ import annotations

import argparse
import contextlib
import os
import shutil
import sys

import numpy as np
from astropy.io import fits
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

VER_DIR = os.path.dirname(os.path.abspath(__file__))
if VER_DIR not in sys.path:
    sys.path.insert(0, VER_DIR)

import vercommon as vc  # noqa: E402
from checks import sky_perturb as sp  # noqa: E402

REFITS_ROOT = os.path.join(VER_DIR, "Re-fits")
WEIRD_DEFAULT = ("20220509G", "20230930A")

_INPUTS = (
    "host_cutout.fits",
    "host_sigma.fits",
    "host_mask.fits",
    "proto_image.fits",
    "constraints.txt",
)


def sandbox_dir(frb: str) -> str:
    return os.path.join(REFITS_ROOT, frb, "sandbox")


def init_sandbox(frb: str, *, force_feedme: bool = False) -> str:
    """Copy production inputs + seed feedme (centroid free: ``1 1``)."""
    prod = vc.host_dir(frb)
    if not os.path.isdir(prod):
        raise FileNotFoundError(f"production dir missing: {prod}")
    sb = sandbox_dir(frb)
    os.makedirs(sb, exist_ok=True)

    for name in _INPUTS:
        src = os.path.join(prod, name)
        dst = os.path.join(sb, name)
        if os.path.isfile(src):
            shutil.copy2(src, dst)
        elif name == "constraints.txt":
            open(dst, "w", encoding="utf-8").close()
        else:
            raise FileNotFoundError(f"missing {src}")

    feedme_dst = os.path.join(sb, "galfit.feedme")
    feedme_src = os.path.join(prod, "galfit.feedme")
    if force_feedme or not os.path.isfile(feedme_dst):
        shutil.copy2(feedme_src, feedme_dst)
        print(f"[{frb}] seeded galfit.feedme from production")
    else:
        print(f"[{frb}] keeping existing galfit.feedme")

    print(f"[{frb}] sandbox ready: {sb}")
    return sb


def _asinh_stretch(img: np.ndarray, lo: float, hi: float) -> np.ndarray:
    x = np.clip((img - lo) / (hi - lo + 1e-30), 0.0, 1.0)
    return np.arcsinh(10.0 * x) / np.arcsinh(10.0)


def make_panel(wkdir: str, out_png: str | None = None) -> str:
    """3-panel data / model / resid from GALFIT ``out.fits``."""
    out_fits = os.path.join(wkdir, "out.fits")
    if not os.path.isfile(out_fits):
        raise FileNotFoundError(f"no out.fits in {wkdir}")
    out_png = out_png or os.path.join(wkdir, "panel.png")

    with fits.open(out_fits) as hdul:
        if len(hdul) < 4:
            raise RuntimeError(f"out.fits has {len(hdul)} HDUs; need ≥4")
        data = np.asarray(hdul[1].data, float)
        model = np.asarray(hdul[2].data, float)
        resid = np.asarray(hdul[3].data, float)

    finite = np.isfinite(data)
    lo, hi = (np.nanpercentile(data[finite], [1, 99]) if finite.any()
              else (0.0, 1.0))
    r_ok = np.isfinite(resid)
    rlim = (np.nanpercentile(np.abs(resid[r_ok]), 98) if r_ok.any() else 1.0)
    rlim = max(float(rlim), 1e-8)

    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.6), constrained_layout=True)
    for ax, arr, title in (
        (axes[0], data, "data"),
        (axes[1], model, "model"),
        (axes[2], resid, "resid"),
    ):
        if title != "resid":
            ax.imshow(_asinh_stretch(arr, lo, hi), origin="lower", cmap="gray",
                      vmin=0, vmax=1, interpolation="nearest")
        else:
            ax.imshow(arr, origin="lower", cmap="RdBu_r",
                      vmin=-rlim, vmax=rlim, interpolation="nearest")
        ax.set_title(title, fontsize=11)
        ax.set_xticks([])
        ax.set_yticks([])

    frb = os.path.basename(os.path.dirname(wkdir.rstrip("\\/")))
    fig.suptitle(f"{frb}  sandbox", fontsize=12)
    fig.savefig(out_png, dpi=140)
    plt.close(fig)
    return out_png


def run_sandbox(frb: str) -> dict:
    sb = sandbox_dir(frb)
    feedme = os.path.join(sb, "galfit.feedme")
    if not os.path.isfile(feedme):
        raise FileNotFoundError(f"no feedme — run with --init first: {feedme}")

    run_galfit = sp._run_galfit()
    # Do NOT wipe galfit.NN or named archives (cent_fix, …) — iterate keeps them.
    for stale in ("fit.log", "out.fits"):
        path = os.path.join(sb, stale)
        if os.path.isfile(path):
            os.remove(path)

    log_path = os.path.join(sb, "galfit_stdout.log")
    with open(log_path, "w", encoding="utf-8") as fh, contextlib.redirect_stdout(fh):
        ok = bool(run_galfit(sb))

    fit = sp._read_result(sb)
    fit["converged"] = ok
    panel = None
    if os.path.isfile(os.path.join(sb, "out.fits")):
        panel = make_panel(sb)
    print(f"[{frb}] GALFIT ok={ok} status={fit.get('status')} panel={panel}")
    if fit.get("status") == "ok":
        print(
            f"  mag={fit.get('mag')} Re={fit.get('re')} n={fit.get('n')} "
            f"q={fit.get('q')} PA={fit.get('pa')} sky={fit.get('sky')}"
        )
    return {"frb": frb, "wkdir": sb, "fit": fit, "panel": panel}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("frb", nargs="?", default=None,
                    help="FRB to run (required unless --init)")
    ap.add_argument("--init", action="store_true",
                    help="create/refresh sandbox inputs (+ seed feedme if missing)")
    ap.add_argument("--reseed-feedme", action="store_true",
                    help="with --init, overwrite galfit.feedme from production")
    ap.add_argument("--panel-only", action="store_true",
                    help="rebuild panel.png from existing out.fits (no GALFIT)")
    ns = ap.parse_args()

    if ns.init and ns.frb is None:
        frbs = list(WEIRD_DEFAULT)
    elif ns.frb:
        frbs = [ns.frb]
    else:
        ap.error("pass an FRB, or --init")

    rc = 0
    for frb in frbs:
        try:
            if ns.init:
                init_sandbox(frb, force_feedme=ns.reseed_feedme)
                if ns.frb is None:
                    # --init alone: setup both, do not run
                    continue
            if ns.panel_only:
                print(f"[{frb}] wrote {make_panel(sandbox_dir(frb))}")
                continue
            if ns.init and ns.frb is None:
                continue
            # --init FRB without wanting a run: only init
            if ns.init and "--init" in sys.argv and ns.frb:
                continue
            run_sandbox(frb)
        except Exception as exc:
            print(f"[{frb}] ERROR: {exc}", file=sys.stderr)
            rc = 1
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
