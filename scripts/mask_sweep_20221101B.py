"""Star-mask sensitivity sweep for FRB 20221101B.

The production mask is the SExtractor segmentation footprint of the bright
Milky Way star 5.5 arcsec NE of the host (392 px, r_eq ~ 11 px). This script
grows / shrinks that footprint by a fixed number of pixels (disk structuring
element), re-runs GALFIT with the production feedme, and tabulates the host
parameters. A second family replaces the footprint with a circular mask of
radius R around the star peak.

Output: pipeline_scripts/verification/Re-fits/20221101B/mask_sweep/
    <variant>/  (GALFIT work dir)
    mask_sweep_results.csv
"""
from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy import ndimage

ROOT = Path(__file__).resolve().parents[1]
PROD = ROOT / "pipeline_scripts/Output/20221101B_all"
OUT = ROOT / "pipeline_scripts/verification/Re-fits/20221101B/mask_sweep"
Q0 = 0.2
STAR_XY = (33, 50)  # 1-based peak of the star in host_cutout.fits

DILATIONS = [-2, -1, 0, 1, 2, 4, 6, 8]
DISK_RADII = [8, 11, 14, 17, 20]


def disk(r: int) -> np.ndarray:
    y, x = np.ogrid[-r:r + 1, -r:r + 1]
    return (x * x + y * y) <= r * r


def star_blob(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    lab, n = ndimage.label(mask > 0)
    sizes = ndimage.sum(mask > 0, lab, range(1, n + 1))
    star = lab == (1 + int(np.argmax(sizes)))
    other = (mask > 0) & ~star
    return star, other


def make_mask(kind: str, val: int) -> np.ndarray:
    m0 = fits.getdata(PROD / "host_mask.fits")
    star, other = star_blob(m0)
    if kind == "dil":
        if val > 0:
            star = ndimage.binary_dilation(star, structure=disk(val))
        elif val < 0:
            star = ndimage.binary_erosion(star, structure=disk(-val))
    elif kind == "disk":
        yy, xx = np.mgrid[0:m0.shape[0], 0:m0.shape[1]]
        cx, cy = STAR_XY[0] - 1, STAR_XY[1] - 1
        star = (xx - cx) ** 2 + (yy - cy) ** 2 <= val * val
    return ((star | other).astype(np.int16))


def run_galfit(wk: Path) -> str:
    proc = subprocess.run(["wsl", "galfit", "galfit.feedme"], cwd=wk, capture_output=True,
                          text=True, encoding="utf-8", errors="replace", stdin=subprocess.DEVNULL, timeout=600)
    (wk / "galfit_stdout.log").write_text((proc.stdout or "") + (proc.stderr or ""), encoding="utf-8")
    return (wk / "fit.log").read_text(encoding="utf-8", errors="replace") if (wk / "fit.log").exists() else ""


def parse(fitlog: str) -> dict:
    """Parse the last sersic/sky block in fit.log."""
    blocks = fitlog.strip().split("-----------------------------------------------------------------------------")
    txt = [b for b in blocks if "sersic" in b][-1]
    lines = txt.splitlines()
    out: dict = {}
    for k, ln in enumerate(lines):
        if ln.strip().startswith("sersic"):
            nums = [float(v) for v in re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", ln.split(":", 1)[1].replace("*", ""))]
            errs = [float(v) for v in re.findall(r"[-+]?\d*\.\d+|[-+]?\d+", lines[k + 1].replace("*", ""))]
            out.update(dict(x=nums[0], y=nums[1], mag=nums[2], re=nums[3], n=nums[4], q=nums[5], pa=nums[6],
                            mag_err=errs[2], re_err=errs[3], n_err=errs[4], q_err=errs[5], pa_err=errs[6]))
            out["flagged"] = "*" in ln
        if ln.strip().startswith("sky"):
            nums = [float(v) for v in re.findall(r"[-+]?\d*\.\d+(?:e[-+]?\d+)?|[-+]?\d+", ln.split(":", 1)[1].replace("*", ""))]
            out["sky"] = nums[2]
        m = re.search(r"Chi\^2/nu = ([\d.]+)", ln)
        if m:
            out["chi2nu"] = float(m.group(1))
    return out


def hubble_i(q: float) -> float:
    if not (Q0 < q <= 1):
        return float("nan")
    return float(np.degrees(np.arccos(np.sqrt((q * q - Q0 * Q0) / (1 - Q0 * Q0)))))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    feedme = (PROD / "galfit.feedme").read_text(encoding="utf-8")
    rows = []
    variants = [("dil", d) for d in DILATIONS] + [("disk", r) for r in DISK_RADII]
    for kind, val in variants:
        name = f"{kind}_{val:+d}" if kind == "dil" else f"{kind}_r{val}"
        wk = OUT / name
        wk.mkdir(exist_ok=True)
        for f in ("host_cutout.fits", "host_sigma.fits", "proto_image.fits", "constraints.txt"):
            shutil.copy(PROD / f, wk / f)
        mask = make_mask(kind, val)
        fits.writeto(wk / "host_mask.fits", mask, overwrite=True)
        (wk / "galfit.feedme").write_text(feedme, encoding="utf-8")
        for stale in ("fit.log", "out.fits", "galfit.01"):
            (wk / stale).unlink(missing_ok=True)
        log = run_galfit(wk)
        r = {"variant": name, "kind": kind, "param": val, "masked_px": int((mask > 0).sum())}
        if log:
            r.update(parse(log))
            r["i_deg"] = hubble_i(r["q"])
        rows.append(r)
        print(name, {k: r.get(k) for k in ("masked_px", "q", "q_err", "pa", "re", "n", "mag", "chi2nu", "i_deg", "flagged")})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "mask_sweep_results.csv", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
