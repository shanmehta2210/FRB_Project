"""Shared r-band cutout fetch helpers (Legacy Survey, PS1, DES via Legacy viewer)."""

from __future__ import annotations

import time
from io import BytesIO
from typing import Optional

import numpy as np
import requests
from astropy.io import fits

_HTTP_RETRIES = 5
_HTTP_BACKOFF_S = 20
PROBE_TIMEOUT = 25
PROBE_RETRIES = 2
VIEWER_TIMEOUT = 360
PS1_TIMEOUT = 600

LS_BASE = "https://www.legacysurvey.org/viewer/cutout.fits"
LS_LAYERS_GLOBAL = ("ls-dr10", "ls-dr9")
LS_LAYERS_SOUTH = ("ls-dr10-south", "ls-dr9-south")
LS_LAYERS_NORTH = ("ls-dr10-north", "ls-dr9-north")
# DES DR1 coadds served through the Legacy Survey viewer (southern footprint).
DES_LAYERS = ("des-dr1",)
DEC_SPLIT = 32.375
PS1_DEC_MIN = -30.0
PS1_PIXSCALE = 0.25

# Standard pipeline cutout: 10 arcmin square at 0.262 arcsec/pixel.
FOV_ARCMIN = 10.0
LS_PIXSCALE = 0.262
PROBE_SIZE = 32
CUTOUT_SIZE = int(round(FOV_ARCMIN * 60.0 / LS_PIXSCALE))  # 2290 px


def _viewer_cutout_url(ra: float, dec: float, layer: str, size: int = CUTOUT_SIZE) -> str:
    return (
        f"{LS_BASE}?ra={ra}&dec={dec}&size={size}&layer={layer}"
        f"&pixscale={LS_PIXSCALE}&bands=r&invvar=True"
    )


def _legacy_url(ra: float, dec: float, layer: str, size: int = CUTOUT_SIZE) -> str:
    return _viewer_cutout_url(ra, dec, layer, size=size)


def _flux_ok(data) -> bool:
    if data is None:
        return False
    arr = np.asarray(data)
    if arr.size == 0:
        return False
    return bool(np.any(np.isfinite(arr)) and np.nanmax(np.abs(arr)) > 0)


def _http_get(url: str, timeout: int = 180, *, max_retries: int | None = None) -> Optional[requests.Response]:
    n_try = _HTTP_RETRIES if max_retries is None else max_retries
    for attempt in range(1, n_try + 1):
        try:
            resp = requests.get(url, timeout=timeout)
            if resp.status_code == 200:
                return resp
        except (requests.RequestException, OSError):
            pass
        if attempt < n_try:
            time.sleep(_HTTP_BACKOFF_S * attempt)
    return None


def fetch_viewer_pair(
    ra: float,
    dec: float,
    layers: tuple[str, ...],
    size: int = CUTOUT_SIZE,
    timeout: int = VIEWER_TIMEOUT,
) -> tuple[Optional[np.ndarray], Optional[object], Optional[np.ndarray], Optional[object], str]:
    """Download flux + invvar from legacysurvey.org/viewer (Legacy or DES layers)."""
    for layer in layers:
        url = _viewer_cutout_url(ra, dec, layer, size=size)
        resp = _http_get(url, timeout=timeout)
        if resp is None:
            continue
        try:
            with fits.open(BytesIO(resp.content)) as hdul:
                if len(hdul) < 2:
                    continue
                flux = hdul[0].data
                invv = hdul[1].data
                if not _flux_ok(flux):
                    continue
                return (
                    np.squeeze(flux),
                    hdul[0].header.copy(),
                    np.squeeze(invv),
                    hdul[1].header.copy(),
                    layer,
                )
        except Exception:
            continue
    return None, None, None, None, ""


def fetch_legacy_pair(
    ra: float,
    dec: float,
    layers: tuple[str, ...],
    size: int = CUTOUT_SIZE,
    timeout: int = 180,
) -> tuple[Optional[np.ndarray], Optional[object], Optional[np.ndarray], Optional[object], str]:
    return fetch_viewer_pair(ra, dec, layers, size=size, timeout=timeout)


def fetch_des_pair(
    ra: float,
    dec: float,
    size: int = CUTOUT_SIZE,
    timeout: int = 180,
) -> tuple[Optional[np.ndarray], Optional[object], Optional[np.ndarray], Optional[object], str]:
    return fetch_viewer_pair(ra, dec, DES_LAYERS, size=size, timeout=timeout)


def legacy_layers_for_dec(dec: float) -> tuple[str, ...]:
    if dec < DEC_SPLIT:
        return LS_LAYERS_SOUTH
    return LS_LAYERS_NORTH


def probe_viewer_layer(
    ra: float, dec: float, layer: str, timeout: int = PROBE_TIMEOUT
) -> bool:
    for layer_name in (layer,) if isinstance(layer, str) else layer:
        url = _viewer_cutout_url(ra, dec, layer_name, size=PROBE_SIZE)
        resp = _http_get(url, timeout=timeout, max_retries=PROBE_RETRIES)
        if resp is None:
            continue
        try:
            with fits.open(BytesIO(resp.content)) as hdul:
                if len(hdul) < 1:
                    continue
                if _flux_ok(hdul[0].data):
                    return True
        except Exception:
            continue
    return False


def preflight_coverage(
    ra: float,
    dec: float,
    *,
    allow_ps1: bool = True,
) -> dict[str, bool]:
    """Fast 32px probes — avoid multi-minute full cutouts when a tier has no data."""
    hem = legacy_layers_for_dec(dec)
    leg_g = probe_viewer_layer(ra, dec, LS_LAYERS_GLOBAL[0])
    leg_h = False if leg_g else probe_viewer_layer(ra, dec, hem[0])
    ps1 = bool(
        allow_ps1
        and dec > PS1_DEC_MIN
        and get_ps1_r_filename(ra, dec, timeout=15) is not None
    )
    des = probe_viewer_layer(ra, dec, DES_LAYERS[0])
    any_ok = leg_g or leg_h or ps1 or des
    return {
        "legacy_global": leg_g,
        "legacy_hemisphere": leg_h,
        "ps1": ps1,
        "des": des,
        "any": any_ok,
        "hemisphere_layer": hem[0],
    }


def probe_legacy(ra: float, dec: float, layer: str, timeout: int = 60) -> bool:
    return probe_viewer_layer(ra, dec, layer, timeout=timeout)


def probe_des(ra: float, dec: float, timeout: int = PROBE_TIMEOUT) -> bool:
    return probe_viewer_layer(ra, dec, DES_LAYERS[0], timeout=timeout)


def get_ps1_r_filename(ra: float, dec: float, timeout: int = 30) -> Optional[str]:
    if dec <= PS1_DEC_MIN:
        return None
    url = (
        "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py?"
        f"ra={ra}&dec={dec}&filters=r&sep=comma"
    )
    try:
        resp = requests.get(url, timeout=timeout)
        if resp.status_code != 200:
            return None
        lines = resp.text.strip().split("\n")
        if len(lines) < 2:
            return None
        keys = lines[0].split(",")
        vals = lines[1].split(",")
        row = dict(zip(keys, vals))
        return row.get("filename") or None
    except Exception:
        return None


def fetch_ps1_pair(
    ra: float,
    dec: float,
    size: int = CUTOUT_SIZE,
    timeout: int = PS1_TIMEOUT,
) -> tuple[Optional[np.ndarray], Optional[object], Optional[np.ndarray], Optional[object]]:
    if dec <= PS1_DEC_MIN:
        return None, None, None, None
    filename = get_ps1_r_filename(ra, dec, timeout=timeout)
    if not filename:
        return None, None, None, None
    ps1_size = int(round(size * LS_PIXSCALE / PS1_PIXSCALE))
    base = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"

    def _get(red: str):
        url = f"{base}?ra={ra}&dec={dec}&size={ps1_size}&format=fits&red={red}"
        resp = _http_get(url, timeout=timeout)
        if resp is None or b"<html" in resp.content[:200].lower():
            return None, None
        with fits.open(BytesIO(resp.content)) as hdul:
            return hdul[0].data, hdul[0].header.copy()

    flux, fh = _get(filename)
    if flux is None:
        return None, None, None, None
    wt_name = filename.replace(".fits", ".wt.fits")
    inv, ih = _get(wt_name)
    if inv is None:
        return None, None, None, None
    return np.squeeze(flux), fh, np.squeeze(inv), ih
