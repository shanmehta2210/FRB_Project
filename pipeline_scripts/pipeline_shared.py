"""Shared helpers for the FRB host pipeline.

Single source of truth for:
  * the SExtractor convolution kernel + neural-net star/galaxy weights
    (identical between Phase 1 and Phase 2),
  * photometric-aperture resolution (YAML-configurable array + production
    aperture, defaulting to the largest aperture), and
  * a consistently formatted logger.

Stdlib-only so it stays importable from every Windows-side phase script.
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
_LOG_FORMAT = "[%(asctime)s] %(levelname)s %(name)s | %(message)s"
_DATE_FORMAT = "%H:%M:%S"
_ROOT_NAME = "frb_pipeline"
_CONFIGURED = False


def get_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """Return a child logger under a single shared stdout handler."""
    global _CONFIGURED
    root = logging.getLogger(_ROOT_NAME)
    if not _CONFIGURED:
        handler = logging.StreamHandler(sys.stdout)
        handler.setFormatter(logging.Formatter(_LOG_FORMAT, _DATE_FORMAT))
        root.addHandler(handler)
        root.setLevel(level)
        root.propagate = False
        _CONFIGURED = True
    return root.getChild(name)


# ---------------------------------------------------------------------------
# Shared SExtractor static templates
# ---------------------------------------------------------------------------
TEMPLATE_CONV = """CONV NORM
# 3x3 ``all-ground'' convolution mask with FWHM = 2 pixels.
1 2 1
2 4 2
1 2 1
"""

TEMPLATE_NNW = """NNW
# Neural Network Weights for the SExtractor star/galaxy classifier (V1.3)
# inputs:	9 for profile parameters + 1 for seeing.
# outputs:	``Stellarity index'' (0.0 to 1.0)
# Seeing FWHM range: from 0.025 to 5.5'' (images must have 1.5 < FWHM < 5 pixels)
# Optimized for Moffat profiles with 2<= beta <= 4.

 3 10 10  1

-1.56604e+00 -2.48265e+00 -1.44564e+00 -1.24675e+00 -9.44913e-01 -5.22453e-01  4.61342e-02  8.31957e-01  2.15505e+00  2.64769e-01
 3.03477e+00  2.69561e+00  3.16188e+00  3.34497e+00  3.51885e+00  3.65570e+00  3.74856e+00  3.84541e+00  4.22811e+00  3.27734e+00

-3.22480e-01 -2.12804e+00  6.50750e-01 -1.11242e+00 -1.40683e+00 -1.55944e+00 -1.84558e+00 -1.18946e-01  5.52395e-01 -4.36564e-01 -5.30052e+00
 4.62594e-01 -3.29127e+00  1.10950e+00 -6.01857e-01  1.29492e-01  1.42290e+00  2.90741e+00  2.44058e+00 -9.19118e-01  8.42851e-01 -4.69824e+00
-2.57424e+00  8.96469e-01  8.34775e-01  2.18845e+00  2.46526e+00  8.60878e-02 -6.88080e-01 -1.33623e-02  9.30403e-02  1.64942e+00 -1.01231e+00
 4.81041e+00  1.53747e+00 -1.12216e+00 -3.16008e+00 -1.67404e+00 -1.75767e+00 -1.29310e+00  5.59549e-01  8.08468e-01 -1.01592e-02 -7.54052e+00
 1.01933e+01 -2.09484e+01 -1.07426e+00  9.87912e-01  6.05210e-01 -6.04535e-02 -5.87826e-01 -7.94117e-01 -4.89190e-01 -8.12710e-02 -2.07067e+01
-5.31793e+00  7.94240e+00 -4.64165e+00 -4.37436e+00 -1.55417e+00  7.54368e-01  1.09608e+00  1.45967e+00  1.62946e+00 -1.01301e+00  1.13514e-01
 2.20336e-01  1.70056e+00 -5.20105e-01 -4.28330e-01  1.57258e-03 -3.36502e-01 -8.18568e-02 -7.16163e+00  8.23195e+00 -1.71561e-02 -1.13749e+01
 3.75075e+00  7.25399e+00 -1.75325e+00 -2.68814e+00 -3.71128e+00 -4.62933e+00 -2.13747e+00 -1.89186e-01  1.29122e+00 -7.49380e-01  6.71712e-01
-8.41923e-01  4.64997e+00  5.65808e-01 -3.08277e-01 -1.01687e+00  1.73127e-01 -8.92130e-01  1.89044e+00 -2.75543e-01 -7.72828e-01  5.36745e-01
-3.65598e+00  7.56997e+00 -3.76373e+00 -1.74542e+00 -1.37540e-01 -5.55400e-01 -1.59195e-01  1.27910e-01  1.91906e+00  1.42119e+00 -4.35502e+00

-1.70059e+00 -3.65695e+00  1.22367e+00 -5.74367e-01 -3.29571e+00  2.46316e+00  5.22353e+00  2.42038e+00  1.22919e+00 -9.22250e-01 -2.32028e+00


 0.00000e+00 
 1.00000e+00 
"""


# ---------------------------------------------------------------------------
# Photometric apertures
# ---------------------------------------------------------------------------
# Standard 15-aperture diameter ladder (px). The largest entry is the default
# production / calibration aperture.
DEFAULT_APERTURE_DIAMS_PX: list[float] = [
    4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 15.0, 20.0, 25.0, 30.0, 40.0,
]
# Legacy base (everything except the production aperture) — used to interpret an
# old scalar ``phot_apertures_px`` without changing historical behaviour.
_LEGACY_BASE_APERTURES: list[float] = DEFAULT_APERTURE_DIAMS_PX[:-1]


def resolve_apertures(sex_cfg: dict) -> tuple[list[float], int, float]:
    """Resolve the aperture configuration.

    Reads from a SExtractor config block (``sextractor`` / ``sextractor_psf``):

      phot_apertures_px      : list of aperture *diameters* in px. A bare scalar
                               is still accepted (legacy) and is appended to the
                               standard base ladder.
      production_aperture_px : aperture used for the zero point AND the
                               production magnitude. If omitted / null, the
                               LARGEST aperture in the array is used.

    Returns ``(aperture_list, production_index, production_diameter_px)``.
    """
    raw = sex_cfg.get("phot_apertures_px", None)
    prod = sex_cfg.get("production_aperture_px", None)

    if isinstance(raw, (list, tuple)) and len(raw) > 0:
        apertures = [float(x) for x in raw]
    elif raw is None:
        apertures = list(DEFAULT_APERTURE_DIAMS_PX)
    else:  # legacy scalar: base ladder + this production aperture
        apertures = list(_LEGACY_BASE_APERTURES) + [float(raw)]
        if prod is None:
            prod = float(raw)

    if not apertures:
        apertures = list(DEFAULT_APERTURE_DIAMS_PX)

    if prod is None:
        prod_index = max(range(len(apertures)), key=lambda i: apertures[i])
    else:
        target = float(prod)
        prod_index = min(range(len(apertures)), key=lambda i: abs(apertures[i] - target))

    return apertures, prod_index, apertures[prod_index]


def format_phot_apertures(apertures) -> str:
    """Render an aperture list as a SExtractor ``PHOT_APERTURES`` value."""
    return ", ".join(f"{float(a):g}" for a in apertures)


def render_param_template(base_param: str, n_aper: int) -> str:
    """Substitute the aperture multiplicity ``{NAPER}`` into a .param template."""
    return base_param.replace("{NAPER}", str(int(n_aper)))


def header_mag_zeropoint_from_fits(path: str | Path) -> float | None:
    """Survey photometric ZP from a flux-image FITS header (PS1, etc.).

    Used when Phase 2 calibration is unavailable so GALFIT ``J)`` matches
    SExtractor instrumental mags (``MAG_ZEROPOINT=0`` in the pipeline).
    """
    import math
    from pathlib import Path as _Path

    p = _Path(path)
    if not p.is_file():
        return None
    try:
        from astropy.io import fits

        with fits.open(p) as hdul:
            header = hdul[0].header
        for key in ("HIERARCH FPA.ZP", "MAGZPT", "ZPT", "PHOTZP"):
            raw = header.get(key)
            if raw is None:
                continue
            val = float(raw)
            if math.isfinite(val) and 15.0 <= val <= 30.0:
                return val
    except (OSError, ValueError, TypeError, ImportError):
        return None
    return None
