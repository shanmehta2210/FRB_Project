"""Field depth helpers for Phase 2 (5σ aperture limiting magnitude)."""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS

DEFAULT_APERTURE_DIAMS_PX = [
    4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 15.0, 20.0, 25.0, 30.0, 40.0,
]
PRODUCTION_APERTURE_INDEX = 14


def parse_phot_apertures_from_sex(sex_path: Path | None) -> list[float]:
    if sex_path is None or not Path(sex_path).is_file():
        return list(DEFAULT_APERTURE_DIAMS_PX)
    text = Path(sex_path).read_text(encoding="utf-8")
    match = re.search(r"^PHOT_APERTURES\s+(.+)$", text, re.MULTILINE)
    if not match:
        return list(DEFAULT_APERTURE_DIAMS_PX)
    return [float(x.strip()) for x in match.group(1).split(",")]


def pixel_scale_from_wcs(image_path: Path) -> float:
    with fits.open(image_path) as hdul:
        w = WCS(hdul[0].header).celestial
    scales = w.proj_plane_pixel_scales()
    return float(np.median([s.to(u.arcsec).value for s in scales]))


def measure_sky_sigma(image_path: Path, seg_path: Path | None) -> float:
    with fits.open(image_path) as hdul:
        img = np.asarray(hdul[0].data, dtype=float)
    if seg_path is not None and Path(seg_path).is_file():
        with fits.open(seg_path) as hdul:
            seg = np.asarray(hdul[0].data)
        sky_pixels = img[seg == 0]
    else:
        sky_pixels = img.ravel()
    med = float(np.nanmedian(sky_pixels))
    mad = float(np.nanmedian(np.abs(sky_pixels - med)))
    sigma_sky = 1.4826 * mad
    if sigma_sky < 0.1:
        sigma_sky, _, _ = sigma_clipped_stats(sky_pixels, sigma=3.0)
        sigma_sky = float(sigma_sky)
    return max(float(sigma_sky), 0.1)


def m_lim_5sigma_aperture(zp: float, sigma_sky: float, diam_px: float) -> tuple[float, float, float]:
    area = math.pi * (diam_px / 2.0) ** 2
    flux_5sigma = 5.0 * sigma_sky * math.sqrt(area)
    if flux_5sigma <= 0:
        raise ValueError(f"Non-positive 5σ flux ({flux_5sigma})")
    m_lim = zp - 2.5 * math.log10(flux_5sigma)
    return float(m_lim), float(flux_5sigma), float(area)


def measure_seeing_fwhm_arcsec(proto_path: Path, px_scale: float) -> float:
    if not Path(proto_path).is_file():
        return float("nan")
    with fits.open(proto_path) as hdul:
        data = hdul[0].data
    prof = np.nanmedian(data, axis=0)
    prof = prof - np.nanmin(prof)
    peak = np.nanmax(prof)
    half = 0.5 * peak
    above = np.where(prof >= half)[0]
    if len(above) < 2:
        return float("nan")
    fwhm_px = float(above[-1] - above[0])
    return fwhm_px * px_scale


def production_aperture_diameter_px(
    config_diam: float | None,
    sex_path: Path | None,
) -> float:
    if config_diam is not None and config_diam > 0:
        return float(config_diam)
    apertures = parse_phot_apertures_from_sex(sex_path)
    if len(apertures) > PRODUCTION_APERTURE_INDEX:
        return float(apertures[PRODUCTION_APERTURE_INDEX])
    return 40.0


def measure_field_depth(
    image_path: Path,
    zp_aper_40px: float,
    *,
    seg_path: Path | None = None,
    proto_path: Path | None = None,
    sex_path: Path | None = None,
    aperture_diameter_px: float | None = None,
) -> dict:
    """
    5σ limiting AB magnitude at the production circular aperture.

    Uses pre-computed ZP (no extra Vizier query).
    """
    image_path = Path(image_path)
    diam = production_aperture_diameter_px(aperture_diameter_px, sex_path)
    px_scale = pixel_scale_from_wcs(image_path)
    sigma_sky = measure_sky_sigma(image_path, seg_path)
    m_lim, flux5, area = m_lim_5sigma_aperture(zp_aper_40px, sigma_sky, diam)
    seeing = measure_seeing_fwhm_arcsec(proto_path, px_scale) if proto_path else float("nan")

    return {
        "m_lim_5sigma_ab": round(m_lim, 3),
        "zp_aper_40px": round(float(zp_aper_40px), 6),
        "sigma_sky_adu_per_pix": round(sigma_sky, 4),
        "aperture_diameter_px": diam,
        "aperture_area_pix2": round(area, 1),
        "flux_5sigma_adu": round(flux5, 4),
        "seeing_fwhm_arcsec": round(seeing, 3) if np.isfinite(seeing) else None,
        "pixel_scale_arcsec": round(px_scale, 6),
        "formula_note": "5σ in pi*(d/2)^2 px^2 at production aperture; not seeing disk",
    }
