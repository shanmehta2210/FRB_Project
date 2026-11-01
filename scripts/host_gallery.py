"""Morphology gallery for the hosts whose Gini-M20 class is unreliable or not Sb-Irr.

Per host: colour image (LS DR10 / DR9 grz, else Pan-STARRS gri), and the r-band data,
GALFIT model and residual of the adopted fit. Residuals show spiral arms, bars and
disks that the single-Sersic model does not capture.

Outputs (Overleaf-ready folder): Publication/manuscript/galaxy/
    {FRB}.png              per-host 4-panel strip
    gallery_{k}.pdf/.png   pages of 8 hosts used in the paper
"""
from __future__ import annotations

import io
import json
import sys
import time
import urllib.request
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.visualization import AsinhStretch, ImageNormalize, PercentileInterval
from astropy.wcs import WCS
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Publication" / "manuscript" / "galaxy"
CACHE = ROOT / "Publication" / "materials" / "galaxy_cache"
GINI = ROOT / "Publication" / "materials" / "appendix" / "robustness" / "gini_m20_56.csv"
WIN = ROOT / "plots" / "plots_null" / "v2" / "frb_vs_hsc_confirmed50" / "confirmed50_q.csv"
PER_PAGE = 8
SEED = pd.read_csv(ROOT / "Publication" / "materials" / "host_properties" / "seed_56.csv", dtype={"frb": str}).set_index("frb")


def flagged_hosts() -> list[str]:
    g = pd.read_csv(GINI, dtype={"frb": str})
    return g.loc[(~g.reliable) | (g.lotz2008 != "late"), "frb"].tolist()


def _comp1(h) -> dict:
    def v(k):
        s = str(h.get(k, "nan")).strip("[]*").split("+/-")[0].strip("[]* ")
        try:
            return float(s)
        except ValueError:
            return float("nan")
    return {k: v(f"1_{k}") for k in ("XC", "YC", "RE", "N", "AR", "PA")}


def host_radec(frb: str) -> tuple[float, float, float]:
    prod = ROOT / "pipeline_scripts" / "Output" / f"{frb}_all"
    meta = json.loads((prod / "cutout_meta.json").read_text())
    c = _comp1(fits.getheader(prod / "out.fits", 2))
    hdr = fits.getheader(ROOT / "large_cutouts" / f"{frb}_flux.fits")
    w = WCS(hdr).celestial
    x0, _, y0, _ = meta["cutout_bounds"]
    ra, dec = w.pixel_to_world_values(x0 + c["XC"] - 1, y0 + c["YC"] - 1)
    ps = abs(float(hdr.get("CD1_1") or hdr.get("CDELT1"))) * 3600
    return float(ra), float(dec), c["RE"] * ps


def _get(url: str) -> bytes | None:
    for _ in range(3):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read()
        except Exception:  # noqa: BLE001
            time.sleep(2)
    return None


def colour_image(frb: str, ra: float, dec: float, size_arcsec: float) -> tuple[np.ndarray | None, str]:
    CACHE.mkdir(parents=True, exist_ok=True)
    npx = int(np.clip(size_arcsec / 0.262, 48, 400))
    for layer, label in (("ls-dr10", "LS DR10"), ("ls-dr9", "LS DR9")):
        p = CACHE / f"{frb}_{layer}.jpg"
        if not p.exists():
            b = _get(f"https://www.legacysurvey.org/viewer/cutout.jpg?ra={ra:.6f}&dec={dec:.6f}"
                     f"&layer={layer}&pixscale=0.262&size={npx}")
            if b:
                p.write_bytes(b)
        if p.exists():
            im = np.asarray(Image.open(p).convert("RGB"))
            if im.std() > 3:
                return im, label
    p = CACHE / f"{frb}_ps1.jpg"
    if not p.exists():
        lst = _get(f"https://ps1images.stsci.edu/cgi-bin/ps1filenames.py?ra={ra}&dec={dec}&filters=gri")
        if lst:
            rows = [r.split() for r in lst.decode().strip().splitlines()[1:]]
            fn = {r[4]: r[7] for r in rows if len(r) > 7}
            if all(b in fn for b in "gri"):
                n = int(np.clip(size_arcsec / 0.25, 48, 400))
                url = (f"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?ra={ra}&dec={dec}&size={n}&format=jpg"
                       f"&red={fn['i']}&green={fn['r']}&blue={fn['g']}")
                b = _get(url)
                if b:
                    p.write_bytes(b)
    if p.exists():
        return np.asarray(Image.open(p).convert("RGB")), "Pan-STARRS"
    return None, ""


def stamps(frb: str, wkdir: str) -> tuple[list[np.ndarray], np.ndarray, dict]:
    with fits.open(ROOT / wkdir / "out.fits") as h:
        data, model, resid = (np.asarray(h[k].data, float) for k in (1, 2, 3))
        c = _comp1(h[2].header)
    mask = np.zeros(data.shape, bool)
    for d in (ROOT / wkdir, ROOT / "pipeline_scripts" / "Output" / f"{frb}_all"):
        p = d / "host_mask.fits"
        if p.exists():
            m = np.asarray(fits.getdata(p)) > 0
            if m.shape == data.shape:
                mask = m
                break
    return [data, model, resid], mask, c


def draw_row(axes, frb: str, wkdir: str, note: str) -> None:
    ra, dec, _ = host_radec(frb)
    (data, model, resid), mask, c = stamps(frb, wkdir)
    hdr = fits.getheader(ROOT / "large_cutouts" / f"{frb}_flux.fits")
    ps = abs(float(hdr.get("CD1_1") or hdr.get("CDELT1"))) * 3600
    col, lab = colour_image(frb, ra, dec, max(data.shape) * ps)
    ax = axes[0]
    if col is not None:
        ax.imshow(col, origin="upper")
        ax.text(0.03, 0.80, lab, transform=ax.transAxes, color="w", fontsize=5, va="top")
    ax.text(0.03, 0.96, f"FRB {frb}", transform=ax.transAxes, color="w", fontsize=6, va="top", weight="bold")
    yy, xx = np.indices(data.shape)
    near = np.hypot(xx + 1 - c["XC"], yy + 1 - c["YC"]) < max(3 * c["RE"], 4)
    good = near & ~mask
    norm = ImageNormalize(data[~mask], interval=PercentileInterval(99.5), stretch=AsinhStretch(0.1))
    for a, img in zip(axes[1:3], (data, model)):
        a.imshow(img, origin="lower", cmap="gray_r", norm=norm)
    r = np.nanpercentile(np.abs(resid[good]), 98) if good.any() else np.nanpercentile(np.abs(resid), 98)
    shown = np.ma.masked_array(resid, mask)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("0.75")
    axes[3].imshow(shown, origin="lower", cmap=cmap, vmin=-r, vmax=r)
    n_free = SEED.loc[frb, "n_galfit"]
    axes[2].text(0.03, 0.96, f"$n_{{\\rm free}}={n_free:.1f}$, $q={c['AR']:.2f}$", transform=axes[2].transAxes,
                 fontsize=5.5, va="top")
    axes[0].text(0.97, 0.04, note, transform=axes[0].transAxes, fontsize=5, va="bottom", ha="right",
                 color="yellow", weight="bold")
    for a in axes:
        a.set_xticks([]); a.set_yticks([])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    win = pd.read_csv(WIN, dtype={"frb": str}).set_index("frb")
    ev = pd.read_csv(ROOT / "Publication" / "materials" / "appendix" / "morphology" / "disk_evidence.csv",
                     dtype={"frb": str}).set_index("frb") if len(sys.argv) < 2 else None
    hosts = flagged_hosts()
    titles = ("Colour", "Data ($r$)", "Model", "Residual")
    for frb in hosts:
        note = ev.loc[frb, "class"] if ev is not None and frb in ev.index else ""
        fig, axes = plt.subplots(1, 4, figsize=(4.4, 1.15))
        draw_row(axes, frb, win.loc[frb, "wkdir"], note)
        fig.subplots_adjust(0.005, 0.01, 0.995, 0.99, 0.03, 0)
        fig.savefig(OUT / f"{frb}.png", dpi=250)
        plt.close(fig)
        print(frb, "ok", flush=True)
    for k in range(0, len(hosts), PER_PAGE):
        page = hosts[k:k + PER_PAGE]
        ncol = 2
        nrow = int(np.ceil(len(page) / ncol))
        fig, axes = plt.subplots(nrow, 4 * ncol, figsize=(7.1, 0.92 * nrow + 0.2))
        axes = np.atleast_2d(axes)
        for j, frb in enumerate(page):
            row, cc = divmod(j, ncol)
            note = ev.loc[frb, "class"] if ev is not None and frb in ev.index else ""
            draw_row(axes[row, 4 * cc:4 * cc + 4], frb, win.loc[frb, "wkdir"], note)
        for j in range(len(page), nrow * ncol):
            row, cc = divmod(j, ncol)
            for a in axes[row, 4 * cc:4 * cc + 4]:
                a.axis("off")
        for cc in range(ncol):
            for t, a in zip(titles, axes[0, 4 * cc:4 * cc + 4]):
                a.set_title(t, fontsize=6, pad=2)
        fig.subplots_adjust(0.005, 0.005, 0.995, 1 - 0.2 / (0.92 * nrow + 0.2), 0.03, 0.03)
        for ext in ("pdf", "png"):
            fig.savefig(OUT / f"gallery_{k // PER_PAGE + 1}.{ext}", dpi=250)
        plt.close(fig)
        print("page", k // PER_PAGE + 1, flush=True)


if __name__ == "__main__":
    main()
