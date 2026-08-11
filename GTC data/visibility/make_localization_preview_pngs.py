#!/usr/bin/env python3
"""1 arcmin PNG previews centered on FRB localization from large_cutouts flux FITS."""

from __future__ import annotations

import argparse
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.visualization import ZScaleInterval
from astropy.wcs import WCS

REPO = Path(__file__).resolve().parents[2]
CUTOUT_DIR = REPO / "large_cutouts"
LOC_CSV = REPO / "master_frb_localization.csv"
DEFAULT_OUT = Path(__file__).resolve().parent / "preview_cutouts_1arcmin"

FRBS_GTC = [
    "20210117A",
    "20210214G",
    "20210809C",
    "20220204A",
    "20220506D",
    "20221116A",
    "20230501A",
    "20230521A",
    "20230521B",
    "20230814B",
    "20230913",
    "20230930A",
    "20240203",
]


def _flux_2d(data: np.ndarray, header) -> np.ndarray:
    data = np.squeeze(data)
    if data.ndim == 2:
        return data.astype(float)
    if data.ndim == 3:
        bands = (header.get("BANDS") or "").lower()
        if "r" in bands:
            return data[bands.index("r")].astype(float)
        for i in range(4):
            if (header.get(f"BAND{i}") or "").lower() == "r":
                return data[i].astype(float)
        return data[0].astype(float)
    raise ValueError(f"unexpected flux shape {data.shape}")


def make_png(frb: str, ra: float, dec: float, flux_path: Path, out_path: Path) -> None:
    with fits.open(flux_path) as hdul:
        header = hdul[0].header
        data = _flux_2d(hdul[0].data, header)
        wcs = WCS(header).celestial

    center = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
    cut = Cutout2D(data, center, (1.0, 1.0) * u.arcmin, wcs=wcs)
    vmin, vmax = ZScaleInterval().get_limits(cut.data)

    fig = plt.figure(figsize=(5, 5), dpi=120)
    ax = fig.add_subplot(1, 1, 1, projection=cut.wcs)
    ax.imshow(cut.data, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)
    cx, cy = cut.wcs.world_to_pixel(center)
    ax.plot(cx, cy, "+", color="red", markersize=14, markeredgewidth=1.5)
    ax.set_title(f"{frb}  ({ra:.4f}, {dec:.4f})  1′ zscale", fontsize=10)
    ax.coords[0].set_axislabel("RA")
    ax.coords[1].set_axislabel("Dec")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frb", nargs="*", help="FRB names (default: GTC month-visible 13)")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    frbs = args.frb or FRBS_GTC
    loc = pd.read_csv(LOC_CSV).set_index("frb")
    ok = 0
    for frb in frbs:
        flux = CUTOUT_DIR / f"{frb}_flux.fits"
        if not flux.is_file():
            print(f"[skip] {frb}: missing {flux.name}")
            continue
        row = loc.loc[frb]
        out = args.out_dir / f"{frb}_1arcmin.png"
        make_png(frb, float(row["ra_deg"]), float(row["dec_deg"]), flux, out)
        print(f"[ok] {out}")
        ok += 1
    print(f"Wrote {ok}/{len(frbs)} PNGs -> {args.out_dir}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
