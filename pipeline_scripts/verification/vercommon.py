"""Shared IO, geometry, masking and cohort handling for the fit verification suite.

Read-only with respect to ``pipeline_scripts/Output/``. Everything the checks
need about a host is assembled once by :func:`load_host` and passed around as a
:class:`HostData`.
"""

from __future__ import annotations

import json
import math
import os
import re
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd
from astropy.io import fits

VER_DIR = os.path.dirname(os.path.abspath(__file__))
PIPELINE_DIR = os.path.dirname(VER_DIR)
REPO = os.path.dirname(PIPELINE_DIR)
OUTPUT_ROOT = os.path.join(PIPELINE_DIR, "Output")
RESULTS_CSV = os.path.join(REPO, "pipeline_galfit_results.csv")

OUT_ROOT = os.path.join(VER_DIR, "outputs")
PER_HOST_ROOT = os.path.join(OUT_ROOT, "per_host")
TABLES_ROOT = os.path.join(OUT_ROOT, "tables")
LOGS_ROOT = os.path.join(OUT_ROOT, "logs")

if REPO not in sys.path:
    sys.path.insert(0, REPO)

from scripts.galfit_fitlog_parse import parse_fitlog_file  # noqa: E402

MAG_CUT = 22.0
BA_CUT = 0.2


# --------------------------------------------------------------------------
# cohort
# --------------------------------------------------------------------------

def cohort(which: str = "all64") -> pd.DataFrame:
    """Production cohort with the ``in_53`` science-cut tag attached."""
    df = pd.read_csv(RESULTS_CSV, dtype={"frb": str})
    mag = pd.to_numeric(df["mag"], errors="coerce")
    ba = pd.to_numeric(df["b_a"], errors="coerce")
    df["in_53"] = (mag <= MAG_CUT) & (ba > BA_CUT)
    if which == "53":
        df = df[df["in_53"]].reset_index(drop=True)
    return df


def host_dir(frb: str) -> str:
    return os.path.join(OUTPUT_ROOT, f"{frb}_all")


_FEEDME_COMP = re.compile(r"^\s*#\s*Component number:\s*(\d+)", re.IGNORECASE)
_FEEDME_PARAM = re.compile(r"^(\s*)(\d+)\)\s*(.*?)(\s*#.*)?$")


def host_n_held_fixed(feedme_path: str) -> float | None:
    """If the first Sérsic in ``galfit.feedme`` has n frozen (flag 0), return it."""
    if not os.path.isfile(feedme_path):
        return None
    try:
        with open(feedme_path, encoding="utf-8", errors="replace") as f:
            lines = f.read().splitlines()
    except OSError:
        return None
    comp_type = ""
    sersic_seen = 0
    for line in lines:
        if _FEEDME_COMP.match(line):
            comp_type = ""
            continue
        m = _FEEDME_PARAM.match(line)
        if not m:
            continue
        num, body = m.group(2), m.group(3)
        if num == "0":
            comp_type = body.split()[0].lower() if body.split() else ""
            if comp_type == "sersic":
                sersic_seen += 1
            continue
        if comp_type == "sersic" and sersic_seen == 1 and num == "5":
            toks = body.split()
            if len(toks) >= 2 and toks[1] == "0":
                try:
                    return float(toks[0])
                except ValueError:
                    return None
            return None
    return None


def per_host_dir(frb: str, create: bool = True) -> str:
    path = os.path.join(PER_HOST_ROOT, frb)
    if create:
        os.makedirs(path, exist_ok=True)
    return path


# --------------------------------------------------------------------------
# JSON helpers
# --------------------------------------------------------------------------

def _jsonable(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, (np.floating, float)):
        val = float(obj)
        return val if math.isfinite(val) else None
    if isinstance(obj, (np.integer, int)) and not isinstance(obj, bool):
        return int(obj)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return _jsonable(obj.tolist())
    return obj


def write_json(path: str, payload: dict) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(_jsonable(payload), f, indent=2, sort_keys=True)
        f.write("\n")


def read_json(path: str) -> dict:
    try:
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
    except (OSError, json.JSONDecodeError):
        return {}
    return data if isinstance(data, dict) else {}


def _safe_json(path: str) -> dict:
    return read_json(path) if os.path.isfile(path) else {}


# --------------------------------------------------------------------------
# GALFIT header / feedme parsing
# --------------------------------------------------------------------------

def _clean_num(token: str) -> float:
    s = str(token)
    for ch in "*[]{}":
        s = s.replace(ch, " ")
    s = s.strip()
    try:
        return float(s)
    except ValueError:
        return float("nan")


def _hdr_value(header, key: str) -> tuple[float, float]:
    """``('166.02 +/- 0.02',)`` -> ``(166.02, 0.02)``; missing -> ``(nan, nan)``."""
    if key not in header:
        return float("nan"), float("nan")
    raw = str(header[key])
    if "+/-" in raw:
        lhs, rhs = raw.split("+/-", 1)
        return _clean_num(lhs), _clean_num(rhs)
    return _clean_num(raw), float("nan")


def parse_out_header(path: str) -> dict:
    """Final GALFIT parameters straight from the ``out.fits`` model HDU header.

    These are guaranteed consistent with the model and residual planes in the
    same file, which ``fit.log`` block selection is not.
    """
    with fits.open(path, memmap=False) as hdul:
        header = hdul[2].header
    comps: list[dict] = []
    sky: dict = {}
    idx = 1
    while f"COMP_{idx}" in header:
        ctype = str(header[f"COMP_{idx}"]).strip().lower()
        if ctype == "sky":
            val, err = _hdr_value(header, f"{idx}_SKY")
            sky = {"level": val, "level_err": err, "comp": idx}
        else:
            xc, xc_err = _hdr_value(header, f"{idx}_XC")
            yc, yc_err = _hdr_value(header, f"{idx}_YC")
            mag, mag_err = _hdr_value(header, f"{idx}_MAG")
            re_px, re_err = _hdr_value(header, f"{idx}_RE")
            n, n_err = _hdr_value(header, f"{idx}_N")
            q, q_err = _hdr_value(header, f"{idx}_AR")
            pa, pa_err = _hdr_value(header, f"{idx}_PA")
            fwhm, fwhm_err = _hdr_value(header, f"{idx}_FWHM")
            c_moff, c_err = _hdr_value(header, f"{idx}_C")
            comps.append(
                {
                    "comp": idx,
                    "type": ctype,
                    "xc": xc, "xc_err": xc_err,
                    "yc": yc, "yc_err": yc_err,
                    "mag": mag, "mag_err": mag_err,
                    "re": re_px, "re_err": re_err,
                    "n": n, "n_err": n_err,
                    "q": q, "q_err": q_err,
                    "pa": pa, "pa_err": pa_err,
                    "fwhm": fwhm, "fwhm_err": fwhm_err,
                    "c": c_moff, "c_err": c_err,
                }
            )
        idx += 1
    chi2nu, _ = _hdr_value(header, "CHI2NU")
    return {"components": comps, "sky": sky, "chi2nu": chi2nu,
            "magzpt": _hdr_value(header, "MAGZPT")[0]}


_FEEDME_KEY = re.compile(r"^\s*([A-Z0-9]+)\)\s*(.*?)\s*(?:#.*)?$")


def parse_feedme(path: str) -> dict:
    """Control-block values plus the byte offsets of the sky component lines."""
    out: dict[str, Any] = {"lines": []}
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            lines = f.read().splitlines()
    except OSError:
        return out
    out["lines"] = lines
    for line in lines:
        m = _FEEDME_KEY.match(line)
        if not m:
            continue
        key, val = m.group(1), m.group(2)
        if key == "J":
            out["magzpt"] = _clean_num(val.split()[0]) if val.split() else float("nan")
        elif key == "K":
            toks = val.split()
            if len(toks) >= 2:
                out["plate_scale"] = (_clean_num(toks[0]), _clean_num(toks[1]))
        elif key == "H":
            out["region"] = [int(_clean_num(t)) for t in val.split()[:4]]
    return out


def plate_scale_arcsec(frb: str) -> float:
    feed = parse_feedme(os.path.join(host_dir(frb), "galfit.feedme"))
    ps = feed.get("plate_scale")
    if ps and math.isfinite(ps[0]) and ps[0] > 0:
        return float(ps[0])
    return 0.262


# --------------------------------------------------------------------------
# PSF metrics
# --------------------------------------------------------------------------

_PSFEX_KEYS = (
    "FWHM_FromFluxRadius_Mean",
    "FWHM_FromFluxRadius_StDev",
    "Ellipticity_Mean",
    "Ellipticity_StDev",
    "Ellipticity1_Mean",
    "Ellipticity2_Mean",
    "Ellipticity1_PixelFree_Mean",
    "Ellipticity2_PixelFree_Mean",
    "NStars_Accepted_Total",
    "Chi2_Mean",
)


def parse_psfex_xml(path: str) -> dict:
    """PSF metrics including the Stokes components production drops.

    ``master_run._parse_psfex_xml`` keeps only the scalar ``Ellipticity_Mean``;
    the leakage test needs ``e1``/``e2`` to get a position angle, so the suite
    carries its own parser rather than changing production behaviour.
    """
    if not os.path.isfile(path):
        return {}
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError:
        return {}
    fields = root.findall(".//{*}FIELD") or root.findall(".//FIELD")
    tds = root.findall(".//{*}TR/{*}TD") or root.findall(".//TR/TD")
    if not fields or not tds:
        return {}
    field_map: dict[str, str] = {}
    for i, fld in enumerate(fields):
        name = fld.attrib.get("name", "")
        if name and i < len(tds) and name not in field_map:
            field_map[name] = tds[i].text or ""
    out: dict[str, float] = {}
    for key in _PSFEX_KEYS:
        raw = field_map.get(key)
        if raw is None:
            continue
        val = _clean_num(raw)
        if math.isfinite(val):
            out[key] = val
    return out


def psf_second_moments(psf: np.ndarray) -> dict:
    """Flux-weighted second moments of the PSF stamp GALFIT actually convolved with.

    Independent of PSFEx's own model fit, so it doubles as a cross-check and a
    fallback when the XML lacks the Stokes columns.
    """
    img = np.asarray(psf, dtype=float)
    img = np.where(np.isfinite(img), img, 0.0)
    img = np.clip(img, 0.0, None)
    total = img.sum()
    if total <= 0:
        return {}
    ny, nx = img.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    xbar = float((img * xx).sum() / total)
    ybar = float((img * yy).sum() / total)
    dx, dy = xx - xbar, yy - ybar
    qxx = float((img * dx * dx).sum() / total)
    qyy = float((img * dy * dy).sum() / total)
    qxy = float((img * dx * dy).sum() / total)
    denom = qxx + qyy
    if denom <= 0:
        return {}
    e1 = (qxx - qyy) / denom
    e2 = 2.0 * qxy / denom
    ellip = math.hypot(e1, e2)
    # Stokes -> position angle measured from +y (GALFIT convention).
    pa_x = 0.5 * math.degrees(math.atan2(e2, e1))
    sigma = (max(qxx * qyy - qxy * qxy, 0.0)) ** 0.25
    return {
        "e1": e1,
        "e2": e2,
        "ellipticity": ellip,
        "pa_deg": _wrap_pa(pa_x - 90.0),
        "fwhm_px": 2.0 * math.sqrt(2.0 * math.log(2.0)) * sigma,
    }


def _wrap_pa(pa_deg: float) -> float:
    """Wrap a position angle into ``[-90, 90)``; PA is defined modulo 180."""
    if not math.isfinite(pa_deg):
        return float("nan")
    return (float(pa_deg) + 90.0) % 180.0 - 90.0


def wrap_pa(pa_deg: float) -> float:
    return _wrap_pa(pa_deg)


# --------------------------------------------------------------------------
# host bundle
# --------------------------------------------------------------------------

@dataclass
class HostData:
    frb: str
    dir: str
    data: np.ndarray
    model: np.ndarray
    resid: np.ndarray
    sigma: np.ndarray
    mask: np.ndarray                 # True where GALFIT ignored the pixel
    xc: float                        # 0-based numpy column of the host centre
    yc: float
    q: float
    pa: float                        # degrees, GALFIT convention (from +y)
    re: float                        # pixels
    n: float
    mag: float
    q_err: float
    pa_err: float
    re_err: float
    mag_err: float
    chi2nu_log: float
    sky_level: float
    neighbours: list[dict]
    meta: dict
    summary: dict
    sky_audit: dict
    plate_scale: float
    magzpt: float
    psf: np.ndarray
    psf_fwhm: float
    residual_closure: float          # max |data - model - resid|
    log_vs_header: dict              # baseline reconciliation
    notes: list[str] = field(default_factory=list)
    # Sub-pixel sampling that best reproduces GALFIT's model plane for this
    # host; set by model_reconstruction_error(), consumed by build_model().
    oversample: int | None = None
    core_radius: float | None = None
    # GALFIT mask ∪ production neighbour mask. Residual metrics (RFF, χ², …)
    # use this so an unmasked companion does not leak into host scores.
    metric_mask: np.ndarray | None = None
    # Every extra GALFIT component (Sérsic / PSF / Moffat), 0-based xy.
    other_components: list[dict] = field(default_factory=list)

    @property
    def shape(self) -> tuple[int, int]:
        return self.data.shape

    @property
    def re_arcsec(self) -> float:
        return self.re * self.plate_scale


def _read_image(path: str) -> np.ndarray | None:
    if not os.path.isfile(path):
        return None
    with fits.open(path, memmap=False) as hdul:
        for hdu in hdul:
            if getattr(hdu, "data", None) is not None and hdu.data.ndim >= 2:
                return np.squeeze(np.asarray(hdu.data, dtype=float))
    return None


def load_host_from_dir(hdir: str, frb: str | None = None) -> HostData:
    """Load a HostData from any GALFIT workdir (production or Re-fits staging).

    Requires ``out.fits`` plus the usual cutout/sigma/mask/proto inputs when
    present. Missing production sidecars (psfex.xml, pipeline_summary, …) are
    tolerated so staged refits can be inspected the same way as production.
    """
    hdir = os.path.abspath(hdir)
    if frb is None:
        base = os.path.basename(hdir.rstrip("\\/"))
        frb = base[:-4] if base.endswith("_all") else base
        if frb.startswith("galfit_"):
            # staged as .../Re-fits/<FRB>/galfit_sky_<label>
            parent = os.path.basename(os.path.dirname(hdir))
            frb = parent or frb

    out_path = os.path.join(hdir, "out.fits")
    if not os.path.isfile(out_path):
        raise FileNotFoundError(f"{frb}: out.fits missing in {hdir}")

    with fits.open(out_path, memmap=False) as hdul:
        data = np.asarray(hdul[1].data, dtype=float)
        model = np.asarray(hdul[2].data, dtype=float)
        resid = np.asarray(hdul[3].data, dtype=float)

    closure = float(np.nanmax(np.abs(data - model - resid)))

    sigma = _read_image(os.path.join(hdir, "host_sigma.fits"))
    if sigma is None:
        # fall back to production host dir if this is a staged refit
        prod = host_dir(frb)
        sigma = _read_image(os.path.join(prod, "host_sigma.fits"))
    if sigma is None:
        sigma = np.full_like(data, np.nan)
    mask_img = _read_image(os.path.join(hdir, "host_mask.fits"))
    if mask_img is None:
        mask_img = _read_image(os.path.join(host_dir(frb), "host_mask.fits"))
    mask = np.zeros(data.shape, dtype=bool) if mask_img is None else (mask_img != 0)
    prod_mask_img = _read_image(os.path.join(host_dir(frb), "host_mask.fits"))
    metric_mask = mask.copy()
    if prod_mask_img is not None and np.asarray(prod_mask_img).shape == mask.shape:
        metric_mask = metric_mask | (np.asarray(prod_mask_img) != 0)

    hdr = parse_out_header(out_path)
    comps = hdr["components"]
    if not comps:
        raise ValueError(f"{frb}: no Sersic component in out.fits header")
    host = comps[0]

    log_params, _strategy = parse_fitlog_file(
        os.path.join(hdir, "fit.log"), sersic_component_index=0
    )
    recon = {}
    for hkey, lkey in (("q", "b_a"), ("pa", "pa"), ("re", "re"), ("n", "n"), ("mag", "mag")):
        lv = log_params.get(lkey)
        hv = host.get(hkey)
        if lv is not None and hv is not None and math.isfinite(float(hv)):
            recon[f"d_{hkey}"] = float(lv) - float(hv)
    recon["chi2nu_log"] = float(log_params.get("chi2nu", float("nan")))

    feed_path = os.path.join(hdir, "galfit.feedme")
    if not os.path.isfile(feed_path):
        feed_path = os.path.join(host_dir(frb), "galfit.feedme")
    feed = parse_feedme(feed_path) if os.path.isfile(feed_path) else {}
    ps = feed.get("plate_scale") or (0.262, 0.262)

    psfex_path = os.path.join(hdir, "psfex.xml")
    if not os.path.isfile(psfex_path):
        psfex_path = os.path.join(host_dir(frb), "psfex.xml")
    psfex = parse_psfex_xml(psfex_path) if os.path.isfile(psfex_path) else {}

    psf = _read_image(os.path.join(hdir, "proto_image.fits"))
    if psf is None:
        psf = _read_image(os.path.join(host_dir(frb), "proto_image.fits"))
    psf_fwhm = psfex.get("FWHM_FromFluxRadius_Mean", float("nan"))
    if (not math.isfinite(psf_fwhm) or psf_fwhm <= 0) and psf is not None:
        psf_fwhm = psf_second_moments(psf).get("fwhm_px", float("nan"))

    notes: list[str] = []
    if closure > 1e-4 * max(1.0, float(np.nanmax(np.abs(data)))):
        notes.append(f"residual closure {closure:.3g} larger than expected")
    for key, val in recon.items():
        if key.startswith("d_") and math.isfinite(val) and abs(val) > 0.02:
            notes.append(f"fit.log vs out.fits header disagree on {key[2:]} by {val:.3g}")

    prod_dir = host_dir(frb)
    return HostData(
        frb=frb,
        dir=hdir,
        data=data,
        model=model,
        resid=resid,
        sigma=sigma,
        mask=mask,
        # GALFIT reports 1-based FITS pixel coordinates.
        xc=float(host["xc"]) - 1.0,
        yc=float(host["yc"]) - 1.0,
        q=float(host["q"]),
        pa=float(host["pa"]),
        re=float(host["re"]),
        n=float(host["n"]),
        mag=float(host["mag"]),
        q_err=float(host["q_err"]),
        pa_err=float(host["pa_err"]),
        re_err=float(host["re_err"]),
        mag_err=float(host["mag_err"]),
        chi2nu_log=float(recon.get("chi2nu_log", float("nan"))),
        sky_level=float(hdr["sky"].get("level", float("nan"))),
        neighbours=[
            {**c, "xc": float(c["xc"]) - 1.0, "yc": float(c["yc"]) - 1.0}
            for c in comps[1:]
            if str(c.get("type", "sersic")).lower() == "sersic"
        ],
        meta=_safe_json(os.path.join(prod_dir, "cutout_meta.json")),
        summary=_safe_json(os.path.join(prod_dir, "pipeline_summary.json")),
        sky_audit=_safe_json(os.path.join(prod_dir, "sky_fit_audit.json")),
        plate_scale=float(ps[0]),
        magzpt=float(feed.get("magzpt", hdr.get("magzpt", float("nan")))),
        psf=psf if psf is not None else np.zeros((1, 1)),
        psf_fwhm=float(psf_fwhm),
        residual_closure=closure,
        log_vs_header=recon,
        notes=notes,
        metric_mask=metric_mask,
        other_components=[
            {**c, "xc": float(c["xc"]) - 1.0, "yc": float(c["yc"]) - 1.0}
            for c in comps[1:]
        ],
    )


def load_host(frb: str) -> HostData:
    return load_host_from_dir(host_dir(frb), frb=frb)


# --------------------------------------------------------------------------
# geometry
# --------------------------------------------------------------------------

def elliptical_coords(
    shape: tuple[int, int], xc: float, yc: float, q: float, pa_deg: float
) -> tuple[np.ndarray, np.ndarray]:
    """Semi-major-axis radius ``a`` and in-ellipse azimuth ``theta``.

    ``x' = a cos(theta)`` along the major axis, ``y' = a q sin(theta)`` along the
    minor axis, so ``a`` is the semi-major axis of the ellipse through a pixel.
    GALFIT measures PA counter-clockwise from ``+y``.
    """
    ny, nx = shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    dx = xx - float(xc)
    dy = yy - float(yc)
    phi = math.radians(float(pa_deg))
    cos_p, sin_p = math.cos(phi), math.sin(phi)
    x_maj = -dx * sin_p + dy * cos_p
    y_min = dx * cos_p + dy * sin_p
    qs = max(float(q), 1e-3)
    a = np.hypot(x_maj, y_min / qs)
    theta = np.arctan2(y_min / qs, x_maj)
    return a, theta


def valid_mask(host: HostData, neighbour_re_factor: float = 1.0) -> np.ndarray:
    """Pixels a residual metric may use.

    Finite data/model, positive sigma, and not in the metric mask. The metric
    mask is the GALFIT mask plus the production neighbour mask, so a companion
    that was unmasked only so GALFIT could fit it is still excluded from RFF
    and χ². Fitted extra components are punched out as well (1 Re for a
    Sérsic; 2 × FWHM for a PSF/Moffat).
    """
    ignore = host.metric_mask if host.metric_mask is not None else host.mask
    ok = (
        np.isfinite(host.data)
        & np.isfinite(host.model)
        & np.isfinite(host.sigma)
        & (host.sigma > 0)
        & ~ignore
    )
    extras = list(host.other_components or []) or list(host.neighbours or [])
    ny, nx = host.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    for comp in extras:
        ctype = str(comp.get("type", "sersic")).lower()
        xc = float(comp.get("xc", float("nan")))
        yc = float(comp.get("yc", float("nan")))
        if not (math.isfinite(xc) and math.isfinite(yc)):
            continue
        re_n = float(comp.get("re", float("nan")))
        if ctype == "sersic" and math.isfinite(re_n) and re_n > 0:
            q_n = float(comp.get("q", 1.0))
            pa_n = float(comp.get("pa", 0.0))
            a_n, _ = elliptical_coords(
                host.shape, xc, yc, q_n if math.isfinite(q_n) else 1.0,
                pa_n if math.isfinite(pa_n) else 0.0,
            )
            ok &= a_n > neighbour_re_factor * re_n
            continue
        rad = 2.0 * float(host.psf_fwhm) if math.isfinite(host.psf_fwhm) else float("nan")
        fwhm = float(comp.get("fwhm", float("nan")))
        if math.isfinite(fwhm) and fwhm > 0:
            rad = max(rad, 2.0 * fwhm) if math.isfinite(rad) else 2.0 * fwhm
        if not math.isfinite(rad) or rad <= 0:
            rad = 8.0
        ok &= np.hypot(xx - xc, yy - yc) > rad
    return ok


def sky_annulus_mask(host: HostData, inner_re: float = 3.0) -> np.ndarray:
    """Unmasked pixels beyond ``inner_re`` Re — the blank-sky reference region."""
    a, _ = elliptical_coords(host.shape, host.xc, host.yc, host.q, host.pa)
    return valid_mask(host) & (a > inner_re * host.re)


def sigma_calibration_ratio(host: HostData, inner_re: float = 3.0) -> dict:
    """Empirical noise divided by the sigma map GALFIT was given.

    Post-fit audit of the sigma map GALFIT used. Combined with
    ``chi2nu_local / r_sigma^2`` this puts absolute chi2nu back on an
    empirically calibrated noise scale; see FIT_VERIFICATION_CHECKS.md §3.4.
    """
    sky = sky_annulus_mask(host, inner_re)
    npix = int(np.count_nonzero(sky))
    if npix < 50:
        # Compact stamp: fall back to whatever unmasked sky exists.
        sky = sky_annulus_mask(host, 1.5)
        npix = int(np.count_nonzero(sky))
    if npix < 20:
        return {"sigma_calibration_ratio": float("nan"), "sigma_cal_npix": npix}
    res = host.resid[sky]
    emp = 1.4826 * float(np.median(np.abs(res - np.median(res))))
    med_sigma = float(np.median(host.sigma[sky]))
    ratio = emp / med_sigma if med_sigma > 0 else float("nan")
    return {
        "sigma_calibration_ratio": ratio,
        "sigma_cal_npix": npix,
        "sky_mad_adu": emp,
        "sigma_map_median_adu": med_sigma,
        "sky_annulus_inner_re": inner_re if npix >= 50 else 1.5,
    }


def azimuthal_annuli(
    a_map: np.ndarray, re: float, fwhm: float, a_min_re: float = 0.2,
    a_max_re: float = 3.0,
) -> list[tuple[float, float]]:
    """Elliptical annulus edges from ``a_min_re`` Re to ``a_max_re`` Re.

    Width is ``max(1 px, FWHM/2)`` so annuli are never narrower than the
    resolution element, which would only correlate neighbouring rings.
    """
    width = max(1.0, (float(fwhm) / 2.0) if math.isfinite(fwhm) and fwhm > 0 else 1.0)
    lo = max(a_min_re * re, 0.5 * width)
    hi = a_max_re * re
    if not math.isfinite(hi) or hi <= lo:
        return []
    edges = np.arange(lo, hi + width, width)
    return [(float(edges[i]), float(edges[i + 1])) for i in range(len(edges) - 1)]


# --------------------------------------------------------------------------
# analytic model reconstruction
# --------------------------------------------------------------------------

def sersic_bn(n: float) -> float:
    """Ciotti & Bertin approximation; accurate to ~1e-6 for n > 0.36."""
    n = float(n)
    return 2.0 * n - 1.0 / 3.0 + 4.0 / (405.0 * n) + 46.0 / (25515.0 * n * n)


# GALFIT (Peng et al. 2002) point-samples pixels far from a component centre and
# only subdivides the inner region where the profile curvature is large. The
# polar-grid integrator is for steep Nuker cusps; for Sersic, square sub-grids
# are what GALFIT itself uses. We mirror that: native sampling everywhere,
# k x k block-average only inside ``CORE_RADIUS_PX`` of the centre.
CORE_RADIUS_PX = 5.0

# Dense enough to catch the non-monotonic match-to-GALFIT surface. Even k is
# allowed: block-average integration does not require odd factors, and hosts
# such as 20230626A prefer k=6 over the old odd-only grid.
OVERSAMPLE_GRID = (1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 15, 21)


def _oversample_factor(re: float, n: float) -> int:
    """Fallback when no per-host choice has been stored yet."""
    if not math.isfinite(re) or re <= 0:
        return 3
    if re >= 15.0 and n <= 2.0:
        return 1
    if re >= 8.0:
        return 3
    if re >= 3.0:
        return 5
    return 9


def sersic_stamp(
    shape: tuple[int, int], xc: float, yc: float, q: float, pa: float,
    re: float, n: float, total_flux: float, oversample: int | None = None,
    core_radius: float = CORE_RADIUS_PX,
) -> np.ndarray:
    """Sersic profile normalized analytically to ``total_flux``.

    Pixel integration follows GALFIT's distance-dependent scheme (Peng et al.
    2002): evaluate at the pixel centre everywhere, and replace only the inner
    ``core_radius`` with a ``k x k`` square sub-grid average. Analytic
    normalization (rather than the stamp sum) keeps the total flux equal to the
    magnitude even when the wing falls off the cutout.
    """
    from scipy.special import gamma  # local: only needed here

    bn = sersic_bn(n)
    # F_tot = 2 pi n Re^2 q e^bn bn^(-2n) Gamma(2n) for I(Re) = 1
    norm = (2.0 * math.pi * n * re * re * q * math.exp(bn)
            * bn ** (-2.0 * n) * float(gamma(2.0 * n)))

    def _point(sub_shape, sxc, syc):
        a, _ = elliptical_coords(sub_shape, sxc, syc, q, pa)
        return np.exp(-bn * ((np.maximum(a, 1e-6) / re) ** (1.0 / n) - 1.0))

    def _block(sub_shape, sxc, syc, k):
        sy, sx = sub_shape
        # Pixel i spans fine cells [i*k, i*k+k); its geometric centre maps to
        # (i + 0.5)*k - 0.5 on the fine grid. Fine-grid radii are in fine-pixel
        # units, so Re is scaled by k to match.
        a, _ = elliptical_coords(
            (sy * k, sx * k), (sxc + 0.5) * k - 0.5, (syc + 0.5) * k - 0.5, q, pa
        )
        fine = np.exp(-bn * ((np.maximum(a, 1e-6) / (re * k)) ** (1.0 / n) - 1.0))
        return fine.reshape(sy, k, sx, k).mean(axis=(1, 3))

    k = _oversample_factor(re, n) if oversample is None else int(oversample)
    prof = _point(shape, xc, yc)
    if k > 1 and math.isfinite(core_radius) and core_radius > 0:
        yy, xx = np.indices(shape)
        core = (xx - xc) ** 2 + (yy - yc) ** 2 <= core_radius ** 2
        if np.any(core):
            ys, xs = np.where(core)
            y0, y1 = int(ys.min()), int(ys.max()) + 1
            x0, x1 = int(xs.min()), int(xs.max()) + 1
            patch = _block((y1 - y0, x1 - x0), xc - x0, yc - y0, k)
            local = core[y0:y1, x0:x1]
            prof[y0:y1, x0:x1][local] = patch[local]
    if not math.isfinite(norm) or norm <= 0:
        norm = float(prof.sum())
    return prof * (total_flux / norm)


def convolve_psf(image: np.ndarray, psf: np.ndarray) -> np.ndarray:
    from scipy.signal import fftconvolve

    kernel = np.clip(np.nan_to_num(np.asarray(psf, dtype=float)), 0.0, None)
    total = kernel.sum()
    if total <= 0:
        return image
    return fftconvolve(image, kernel / total, mode="same")


# Core radii searched when matching GALFIT. ``inf`` = oversample the whole
# stamp (old uniform behaviour); finite values mirror Peng's distance cut.
CORE_RADIUS_GRID = (5.0, 8.0, math.inf)


def build_model(
    host: "HostData",
    *,
    q: float | None = None,
    pa: float | None = None,
    re: float | None = None,
    n: float | None = None,
    mag: float | None = None,
    include_sky: bool = True,
    include_neighbours: bool = True,
    oversample: int | None = None,
    core_radius: float | None = None,
) -> np.ndarray:
    """Rebuild GALFIT's PSF-convolved model, optionally with the host perturbed.

    Reproduces the production model plane to a fraction of a percent of peak,
    which is what lets the Fourier estimator be calibrated by finite differences
    instead of an analytic gradient.
    """
    def _flux(m: float) -> float:
        return 10.0 ** ((host.magzpt - float(m)) / 2.5)

    if oversample is None:
        oversample = getattr(host, "oversample", None)
    if core_radius is None:
        core_radius = getattr(host, "core_radius", CORE_RADIUS_PX)

    img = sersic_stamp(
        host.shape, host.xc, host.yc,
        host.q if q is None else float(q),
        host.pa if pa is None else float(pa),
        host.re if re is None else float(re),
        host.n if n is None else float(n),
        _flux(host.mag if mag is None else mag),
        oversample=oversample,
        core_radius=core_radius,
    )
    if include_neighbours:
        for comp in host.neighbours:
            if not all(math.isfinite(float(comp.get(k, float("nan"))))
                       for k in ("xc", "yc", "q", "pa", "re", "n", "mag")):
                continue
            img = img + sersic_stamp(
                host.shape, comp["xc"], comp["yc"], comp["q"], comp["pa"],
                comp["re"], comp["n"], _flux(comp["mag"]),
                oversample=oversample,
                core_radius=core_radius,
            )
    out = convolve_psf(img, host.psf)
    if include_sky and math.isfinite(host.sky_level):
        out = out + host.sky_level
    return out


def _recon_error(
    host: "HostData", oversample: int | None, core_radius: float,
) -> tuple[float, float]:
    recon = build_model(host, oversample=oversample, core_radius=core_radius)
    peak = float(np.nanmax(np.abs(host.model - host.sky_level)))
    diff = float(np.nanmax(np.abs(recon - host.model)))
    fsum_g = float(np.nansum(host.model - host.sky_level))
    fsum_r = float(np.nansum(recon - host.sky_level))
    return (diff / peak if peak > 0 else float("nan"),
            (fsum_r / fsum_g - 1.0) if fsum_g != 0 else float("nan"))


def model_reconstruction_error(host: "HostData") -> dict:
    """Pick (k, core_radius) that best reproduces GALFIT's model plane.

    GALFIT's integrator is distance-dependent and not identical to a uniform
    high-order pixel integral, so the match-to-GALFIT error is non-monotonic in
    k: larger k is not always better. Both the subdivision factor and the core
    radius are therefore chosen by direct comparison against ``out.fits`` HDU 2.
    """
    best: tuple[float, int | None, float | None, float] = (
        float("inf"), None, None, float("nan"),
    )
    tried: dict[str, float] = {}
    for core in CORE_RADIUS_GRID:
        for k in OVERSAMPLE_GRID:
            try:
                err, flux = _recon_error(host, k, core)
            except Exception:
                continue
            label = f"k{k}_r{'full' if math.isinf(core) else f'{core:g}'}"
            tried[label] = err
            if math.isfinite(err) and err < best[0]:
                best = (err, k, core, flux)
    best_err, best_k, best_core, best_flux = best
    if best_k is None:
        return {"model_recon_status": "failed"}
    host.oversample = best_k
    host.core_radius = best_core
    return {
        "model_recon_status": "ok",
        "model_recon_max_frac": best_err,
        "model_recon_flux_frac": best_flux,
        "model_recon_oversample": best_k,
        "model_recon_core_radius": (
            None if best_core is None or math.isinf(best_core) else best_core
        ),
        "model_recon_by_oversample": tried,
    }


def counts_to_mu(flux_per_px: np.ndarray, magzpt: float, plate_scale: float) -> np.ndarray:
    """Surface brightness in mag/arcsec^2 from counts per pixel."""
    with np.errstate(divide="ignore", invalid="ignore"):
        return magzpt - 2.5 * np.log10(np.asarray(flux_per_px) / plate_scale**2)


def weighted_mean(values: np.ndarray, errors: np.ndarray) -> tuple[float, float]:
    """Inverse-variance weighted mean and its error; ``(nan, nan)`` if empty."""
    v = np.asarray(values, dtype=float)
    e = np.asarray(errors, dtype=float)
    ok = np.isfinite(v) & np.isfinite(e) & (e > 0)
    if not np.any(ok):
        return float("nan"), float("nan")
    w = 1.0 / e[ok] ** 2
    mean = float(np.sum(w * v[ok]) / np.sum(w))
    err = float(1.0 / math.sqrt(np.sum(w)))
    return mean, err
