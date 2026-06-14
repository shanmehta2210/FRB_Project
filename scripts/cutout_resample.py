"""Resample flux / invvar cutouts to the pipeline standard grid (0.262 arcsec/px, 2290 px)."""

from __future__ import annotations

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import zoom

TARGET_PIXSCALE = 0.262  # arcsec/pixel
TARGET_SIZE = 2290
FOV_ARCMIN = 10.0


def _native_pixscale(header) -> float:
    w = WCS(header)
    if w.is_celestial and w.pixel_scale_matrix is not None:
        scales = w.pixel_scale_matrix.diagonal()[:2]
        return float(np.mean(np.abs(scales)))
    cdelt = header.get("CDELT1") or header.get("CD1_1")
    if cdelt is not None:
        return abs(float(cdelt)) * 3600.0
    return TARGET_PIXSCALE


def _squeeze2d(data: np.ndarray) -> np.ndarray:
    data = np.squeeze(np.asarray(data, dtype=np.float64))
    while data.ndim > 2:
        data = data[0]
    return data


def _crop_zoom_centered(
    data: np.ndarray,
    header,
    ra: float,
    dec: float,
    native_scale: float,
    *,
    variance: bool,
    center_on_array: bool = False,
) -> np.ndarray:
    """Crop/zoom so the target lands on the output stamp centre.

    PS1 (and other SIAP-style) cutouts may have CRVAL offset from the stamp
    centre; default ``center_on_array=False`` uses WCS(ra,dec) for the crop centre.
    """
    data = _squeeze2d(data)
    if center_on_array:
        xc, yc = (data.shape[1] + 1) / 2.0, (data.shape[0] + 1) / 2.0
    else:
        w = WCS(header)
        try:
            xc, yc = w.wcs_world2pix(ra, dec, 0)
        except Exception:
            xc, yc = (data.shape[1] + 1) / 2.0, (data.shape[0] + 1) / 2.0

    factor = native_scale / TARGET_PIXSCALE
    if abs(factor - 1.0) > 1e-6:
        data = zoom(data, factor, order=1)
        if variance:
            data = data * (factor ** 2)
        xc *= factor
        yc *= factor

    h, w = data.shape
    x0 = int(round(xc - TARGET_SIZE / 2.0))
    y0 = int(round(yc - TARGET_SIZE / 2.0))
    out = np.zeros((TARGET_SIZE, TARGET_SIZE), dtype=np.float64)
    # paste with clipping
    sy0 = max(0, -y0)
    sx0 = max(0, -x0)
    sy1 = min(h, TARGET_SIZE - y0)
    sx1 = min(w, TARGET_SIZE - x0)
    dy0 = max(0, y0)
    dx0 = max(0, x0)
    if sy1 > sy0 and sx1 > sx0:
        out[dy0 : dy0 + (sy1 - sy0), dx0 : dx0 + (sx1 - sx0)] = data[sy0:sy1, sx0:sx1]
    return out


def _output_header(ra: float, dec: float, template=None) -> fits.Header:
    hdr = template.copy() if template is not None else fits.Header()
    hdr["NAXIS"] = 2
    hdr["NAXIS1"] = TARGET_SIZE
    hdr["NAXIS2"] = TARGET_SIZE
    hdr["CDELT1"] = -TARGET_PIXSCALE / 3600.0
    hdr["CDELT2"] = TARGET_PIXSCALE / 3600.0
    hdr["CRVAL1"] = ra
    hdr["CRVAL2"] = dec
    hdr["CRPIX1"] = (TARGET_SIZE + 1) / 2.0
    hdr["CRPIX2"] = (TARGET_SIZE + 1) / 2.0
    hdr["CTYPE1"] = "RA---TAN"
    hdr["CTYPE2"] = "DEC--TAN"
    hdr["BAND"] = "r"
    hdr["PIXSCALE"] = TARGET_PIXSCALE
    return hdr


def adu_to_nanomaggies(flux_adu: np.ndarray, zeropoint_ab: float) -> np.ndarray:
    """Convert SkyMapper image ADU to Legacy-like nanomaggies.

    SkyMapper: mag = ZP - 2.5*log10(counts).  Legacy: mag = 22.5 - 2.5*log10(nanomaggies).
    """
    scale = 10.0 ** ((22.5 - float(zeropoint_ab)) / 2.5)
    return np.asarray(flux_adu, dtype=np.float64) * scale


def invvar_from_sky(
    flux_nmgy: np.ndarray,
    good: np.ndarray,
    *,
    floor_sigma: float = 1e-6,
) -> np.ndarray:
    """Flat invvar map from robust sky MAD on good pixels (Legacy-style scale)."""
    out = np.zeros_like(flux_nmgy, dtype=np.float64)
    if not np.any(good):
        return out
    sky = flux_nmgy[good]
    med = float(np.median(sky))
    mad = float(np.median(np.abs(sky - med)))
    sigma = max(1.4826 * mad, floor_sigma)
    iv = 1.0 / (sigma * sigma)
    out[good] = iv
    return out


def standardize_pair(
    flux: np.ndarray,
    flux_header,
    invvar: np.ndarray,
    invvar_header,
    ra: float,
    dec: float,
    native_pixscale: float | None = None,
    *,
    center_on_array: bool = False,
):
    """Return (flux, flux_hdr, invvar, invvar_hdr) on the standard grid."""
    scale = native_pixscale or _native_pixscale(flux_header)
    flux_out = _crop_zoom_centered(
        flux, flux_header, ra, dec, scale, variance=False, center_on_array=center_on_array
    )
    inv_out = _crop_zoom_centered(
        invvar, invvar_header, ra, dec, scale, variance=True, center_on_array=center_on_array
    )
    hdr = _output_header(ra, dec, flux_header)
    inv_hdr = hdr.copy()
    return flux_out, hdr, inv_out, inv_hdr


def write_standardized(
    flux_path: str,
    invvar_path: str,
    flux: np.ndarray,
    flux_header,
    invvar: np.ndarray,
    invvar_header,
    ra: float,
    dec: float,
    native_pixscale: float | None = None,
    *,
    center_on_array: bool = False,
) -> None:
    f, fh, iv, ivh = standardize_pair(
        flux,
        flux_header,
        invvar,
        invvar_header,
        ra,
        dec,
        native_pixscale,
        center_on_array=center_on_array,
    )
    fits.writeto(flux_path, f, fh, overwrite=True)
    fits.writeto(invvar_path, iv, ivh, overwrite=True)
