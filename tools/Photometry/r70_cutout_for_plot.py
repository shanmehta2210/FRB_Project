"""
Resolve which FITS image to use for r70 AstroPath overlay plots.

Preference:
1. ``coadded_astrometrically_corrected_rband_r70_cutout.fits`` if its on-sky
   field of view is approximately 1 arcminute (see tolerance).
2. Otherwise, if the full coadd exists, extract a 1' x 1' cutout centered on
   the FRB position and write
   ``coadded_astrometrically_corrected_rband_r70_1arcmin_frbcenter.fits``.

FOV is derived from the image WCS (not assumed to match Legacy Survey pixel scale).
"""
from __future__ import annotations

import os
from typing import Tuple

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

CUTOUT_CANDIDATE = "coadded_astrometrically_corrected_rband_r70_cutout.fits"
FULL_COADD = "coadded_astrometrically_corrected_rband_r70.fits"
OUT_ONE_ARCMIN = "coadded_astrometrically_corrected_rband_r70_1arcmin_frbcenter.fits"

# Accept ~1' cutouts (WCS-based; allows small non-square residuals)
FOV_TOL_ARCMIN = 0.12


def fov_arcmin_from_hdu(hdu) -> Tuple[float, float]:
    """Return approximate (width, height) of the image in arcminutes from WCS."""
    wcs = WCS(hdu.header)
    data = np.asarray(hdu.data)
    if data.ndim != 2:
        data = data[0]
    ny, nx = data.shape
    cx, cy = (nx - 1) / 2.0, (ny - 1) / 2.0
    c0 = wcs.pixel_to_world(0, cy)
    c1 = wcs.pixel_to_world(nx - 1, cy)
    c2 = wcs.pixel_to_world(cx, 0)
    c3 = wcs.pixel_to_world(cx, ny - 1)
    width_as = c0.separation(c1).to(u.arcmin).value
    height_as = c2.separation(c3).to(u.arcmin).value
    return float(width_as), float(height_as)


def _is_about_one_arcmin(w: float, h: float) -> bool:
    return (
        abs(w - 1.0) <= FOV_TOL_ARCMIN
        and abs(h - 1.0) <= FOV_TOL_ARCMIN
    )


def extract_one_arcmin_cutout(
    full_fits_path: str,
    center: SkyCoord,
    out_path: str,
) -> str:
    with fits.open(full_fits_path) as hdul:
        hdu = hdul[0]
        wcs = WCS(hdu.header)
        data = np.asarray(hdu.data)
        if data.ndim != 2:
            data = np.squeeze(data)
            if data.ndim != 2:
                data = data[0]
        size = u.Quantity((1.0, 1.0), u.arcmin)
        cut = Cutout2D(data, center, size, wcs=wcs, mode="partial")
        phdu = fits.PrimaryHDU(data=cut.data, header=cut.wcs.to_header())
        phdu.writeto(out_path, overwrite=True)
    return out_path


def resolve_r70_plot_image(
    work_dir: str,
    center: SkyCoord,
) -> Tuple[str, str]:
    """
    Return (path_to_fits, human-readable log line about FOV / which file was chosen).
    """
    work_dir = os.path.abspath(work_dir)
    cut_cand = os.path.join(work_dir, CUTOUT_CANDIDATE)
    full_path = os.path.join(work_dir, FULL_COADD)
    out_1am = os.path.join(work_dir, OUT_ONE_ARCMIN)

    if os.path.isfile(cut_cand):
        with fits.open(cut_cand) as hdul:
            w, h = fov_arcmin_from_hdu(hdul[0])
        msg = f"Found {CUTOUT_CANDIDATE}: FOV ~ {w:.3f}' x {h:.3f}' (from WCS)."
        if _is_about_one_arcmin(w, h):
            return cut_cand, msg + " Using as 1' display image."
        msg += f" Not within ~1' (+/-{FOV_TOL_ARCMIN}')."
        if os.path.isfile(full_path):
            extract_one_arcmin_cutout(full_path, center, out_1am)
            return (
                out_1am,
                msg
                + f" Extracted 1' cutout from {FULL_COADD} -> {OUT_ONE_ARCMIN}.",
            )
        return cut_cand, msg + f" {FULL_COADD} missing; using cutout anyway."

    if os.path.isfile(full_path):
        extract_one_arcmin_cutout(full_path, center, out_1am)
        return (
            out_1am,
            f"No {CUTOUT_CANDIDATE}; extracted 1' cutout from {FULL_COADD} -> {OUT_ONE_ARCMIN}.",
        )

    raise FileNotFoundError(
        f"Need either {CUTOUT_CANDIDATE} or {FULL_COADD} under {work_dir}."
    )
