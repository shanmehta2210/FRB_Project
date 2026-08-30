"""Stage a GALFIT re-fit and re-run the full verification suite on it.

A re-fit is not just a new ``out.fits`` + panel. The verification checks
(RFF, isophotes, Fourier, sky perturbation, AstroPhot, visual, …) are re-run
against the new model and written under ``Re-fits/<FRB>/<label>/``, with the
production panel copied alongside for comparison.

Examples
--------
    python run_refit.py 20181112A --sky-adu 6.317e-5 --fix-n 1 --label n1_sky
    python run_refit.py 20181112A --sky-from-protocol --fix-n 1
    python run_refit.py 20221101B --add-psf-at-masked-star --psf-dmag 2 --label psf
    python run_refit.py 20221101B --add-psf-at-masked-star --fix-n 1 --fix-psf-xy --label n1_psf_xy
    python run_refit.py 20221101B --add-sersic-at-masked-star
    python run_refit.py 20221101B --add-moffat-at-masked-star
"""

from __future__ import annotations

import argparse
import contextlib
import math
import os
import re
import shutil
import sys
import time
from typing import Any

import numpy as np
from astropy.io import fits

VER_DIR = os.path.dirname(os.path.abspath(__file__))
if VER_DIR not in sys.path:
    sys.path.insert(0, VER_DIR)

import vercommon as vc  # noqa: E402
from checks import sky_perturb as sp  # noqa: E402
import run_verification as rv  # noqa: E402

REFITS_ROOT = os.path.join(VER_DIR, "Re-fits")

_SPREAD_STAR_MAX = 0.005
_SPREAD_SIGMA = 3.0
_CLASS_STAR_MAX = 0.75


def _ldac_table(path: str):
    from astropy.table import Table
    if not os.path.isfile(path):
        return None
    try:
        return Table.read(path, hdu=2, format="fits")
    except Exception:
        try:
            return Table.read(path, format="ascii.sextractor")
        except Exception:
            return None


def _col(row, name, default=float("nan")):
    if name not in row.colnames:
        return default
    try:
        return float(row[name])
    except (TypeError, ValueError):
        return default


def masked_object_star_cut(prod_dir: str, seg_id: int | None) -> dict[str, Any]:
    """CLASS_STAR / SPREAD used by Phase 3a to mask the neighbor as a star."""
    out: dict[str, Any] = {
        "seg_id": seg_id,
        "class_star": None,
        "class_star_threshold": _CLASS_STAR_MAX,
        "stellar_by_class_star": None,
        "spread_model": None,
        "spreaderr_model": None,
        "spread_plus_3sig": None,
        "spread_star_max": _SPREAD_STAR_MAX,
        "point_like_by_spread": None,
        "psfex_fwhm_px": None,
        "fwhm_image_px": None,
        "flux_radius_px": None,
        "mag_auto": None,
        "psf_number": None,
        "spread_match_arcsec": None,
        "elongation": None,
        "theta_image": None,
        "flags": None,
        "mask_reason": None,
        "psfex_chi2": None,
        "psfex_ellipticity": None,
        "fwhm_psfcat_px": None,
    }
    xml = os.path.join(prod_dir, "psfex.xml")
    if os.path.isfile(xml):
        psfex = vc.parse_psfex_xml(xml)
        out["psfex_fwhm_px"] = psfex.get("FWHM_FromFluxRadius_Mean")
        out["psfex_chi2"] = psfex.get("Chi2_Mean")
        out["psfex_ellipticity"] = psfex.get("Ellipticity_Mean")
    if seg_id is None:
        return out
    cat = _ldac_table(os.path.join(prod_dir, "image.cat"))
    if cat is None or "NUMBER" not in cat.colnames:
        return out
    sel = cat["NUMBER"] == int(seg_id)
    if not np.any(sel):
        return out
    row = cat[sel][0]
    cs = _col(row, "CLASS_STAR")
    out["class_star"] = cs
    out["stellar_by_class_star"] = bool(cs >= _CLASS_STAR_MAX) if np.isfinite(cs) else False
    out["fwhm_image_px"] = _col(row, "FWHM_IMAGE")
    out["flux_radius_px"] = _col(row, "FLUX_RADIUS")
    out["mag_auto"] = _col(row, "MAG_AUTO")
    out["elongation"] = _col(row, "ELONGATION")
    out["theta_image"] = _col(row, "THETA_IMAGE")
    out["flags"] = _col(row, "FLAGS")
    out["mask_reason"] = "mask_star" if out["stellar_by_class_star"] else "would_fit"
    psf = _ldac_table(os.path.join(prod_dir, "image.psf.cat"))
    if psf is None or "SPREAD_MODEL" not in psf.colnames:
        return out
    from astropy.coordinates import SkyCoord
    ra = _col(row, "ALPHAWIN_J2000")
    dec = _col(row, "DELTAWIN_J2000")
    if not (np.isfinite(ra) and np.isfinite(dec)):
        return out
    coord = SkyCoord(ra=ra, dec=dec, unit="deg", frame="icrs")
    psf_c = SkyCoord(ra=np.asarray(psf["ALPHAWIN_J2000"]),
                     dec=np.asarray(psf["DELTAWIN_J2000"]),
                     unit="deg", frame="icrs")
    seps = coord.separation(psf_c).arcsec
    j = int(np.argmin(seps))
    out["spread_match_arcsec"] = float(seps[j])
    if float(seps[j]) > 0.5:
        return out
    sm = float(psf["SPREAD_MODEL"][j])
    se = float(psf["SPREADERR_MODEL"][j])
    out["psf_number"] = int(psf["NUMBER"][j])
    out["spread_model"] = sm
    out["spreaderr_model"] = se
    out["spread_plus_3sig"] = sm + _SPREAD_SIGMA * se
    out["point_like_by_spread"] = bool(out["spread_plus_3sig"] < _SPREAD_STAR_MAX)
    if "FWHM_IMAGE" in psf.colnames:
        out["fwhm_psfcat_px"] = float(psf["FWHM_IMAGE"][j])
    if "FLUX_RADIUS" in psf.colnames and not np.isfinite(out["flux_radius_px"] or float("nan")):
        out["flux_radius_px"] = float(psf["FLUX_RADIUS"][j])
    if out["stellar_by_class_star"]:
        out["mask_reason"] = "mask_star"
    elif out["point_like_by_spread"]:
        out["mask_reason"] = "mask_spread"
    else:
        out["mask_reason"] = "would_fit"
    return out


_COMP_HEADER = re.compile(r"^\s*#\s*Component number:\s*(\d+)", re.IGNORECASE)
_PARAM = re.compile(r"^(\s*)(\d+)\)\s*(.*?)(\s*#.*)?$")

# Sidecars needed so verification checks see the same supporting files as
# production (catalog sky MAD, PSFEx XML, constraints, …).
_SIDECARS = (
    "host_cutout.fits",
    "host_sigma.fits",
    "host_mask.fits",
    "proto_image.fits",
    "image.cat",
    "image.psf.cat",
    "psfex.xml",
    "constraints.txt",
    "cutout_meta.json",
    "sky_fit_audit.json",
    "pipeline_summary.json",
    "reference_photometry.json",
    "zero_points.json",
)


def _edit_feedme_fixed(
    lines: list[str],
    comps: list[dict],
    *,
    sky_value: float | None,
    sky_comp: int,
    fix_n: float | None,
    reseed: bool = True,
    sky_seed_free: float | None = None,
    fix_host_xy: bool = False,
) -> list[str]:
    """Reseed Sérsics at best-fit; optionally fix sky and/or host n.

    ``sky_value`` → sky held fixed (flag 0).
    ``sky_seed_free`` → sky free (flag 1) but reseeding at that ADU (used when
    only ``n`` is fixed so sky can still float from the production level).
    """
    out = list(lines)
    comp_id = 0
    comp_type = ""
    sersic_seen = 0
    params: dict | None = None

    for i, line in enumerate(lines):
        head = _COMP_HEADER.match(line)
        if head:
            comp_id = int(head.group(1))
            comp_type = ""
            params = None
            continue
        m = _PARAM.match(line)
        if not m:
            continue
        indent, num, body, comment = m.group(1), m.group(2), m.group(3), m.group(4) or ""
        if num == "0":
            comp_type = body.split()[0].lower() if body.split() else ""
            if comp_type == "sersic":
                params = comps[sersic_seen] if sersic_seen < len(comps) else None
                sersic_seen += 1
            continue

        if comp_type == "sky" and comp_id == sky_comp and num == "1":
            if sky_value is not None:
                out[i] = f"{indent}1) {sky_value:.8f} 0{comment}"
            elif sky_seed_free is not None:
                out[i] = f"{indent}1) {sky_seed_free:.8f} 1{comment}"
            continue

        if comp_type != "sersic" or params is None:
            continue

        # Host = first Sérsic. Fix n only on the host component.
        is_host = sersic_seen == 1
        if num == "5" and is_host and fix_n is not None:
            out[i] = f"{indent}5) {float(fix_n):.4f} 0{comment}"
            continue
        if num == "1" and is_host and fix_host_xy:
            out[i] = f"{indent}1) {params['xc']:.4f} {params['yc']:.4f} 0 0{comment}"
            continue

        if not reseed:
            continue
        new = {
            "1": f"{params['xc']:.4f} {params['yc']:.4f} 1 1",
            "3": f"{params['mag']:.4f} 1",
            "4": f"{params['re']:.4f} 1",
            "5": f"{params['n']:.4f} 1",
            "9": f"{params['q']:.4f} 1",
            "10": f"{params['pa']:.4f} 1",
        }.get(num)
        # Don't overwrite a just-fixed n with the free reseed.
        if num == "5" and is_host and fix_n is not None:
            continue
        if new is not None:
            out[i] = f"{indent}{num}) {new}{comment}"
    return out


def _strip_param_constraints(src: str, dst: str, drop: set[str]) -> None:
    """Drop constraint rows for named params (e.g. {'n'} when n is fixed)."""
    try:
        with open(src, encoding="utf-8", errors="replace") as f:
            lines = f.read().splitlines()
    except OSError:
        lines = []
    kept = []
    for ln in lines:
        toks = ln.split()
        if len(toks) >= 2 and toks[1].lower() in drop:
            continue
        kept.append(ln)
    with open(dst, "w", newline="\n", encoding="utf-8") as f:
        f.write("\n".join(kept) + ("\n" if kept else ""))


def _cutout_bounds(meta: dict) -> tuple[int, int, int, int] | None:
    b = meta.get("cutout_bounds")
    if not isinstance(b, (list, tuple)) or len(b) != 4:
        return None
    xmin, xmax, ymin, ymax = (int(v) for v in b)
    return xmin, xmax, ymin, ymax


def _star_from_segmap(prod_dir: str, meta: dict) -> dict | None:
    """Centroid of the masked (non-host) seg island inside the stamp."""
    seg_path = os.path.join(prod_dir, "segmentation_map.fits")
    if not os.path.isfile(seg_path):
        return None
    bounds = _cutout_bounds(meta)
    if bounds is None:
        return None
    xmin, xmax, ymin, ymax = bounds
    host_id = int(meta.get("host_number") or 0)
    seg = np.asarray(fits.getdata(seg_path))
    while seg.ndim > 2:
        seg = seg[0]
    sub = seg[ymin:ymax, xmin:xmax]
    ids = [int(v) for v in np.unique(sub) if int(v) > 0 and int(v) != host_id]
    if not ids:
        return None
    cut = np.asarray(fits.getdata(os.path.join(prod_dir, "host_cutout.fits")), float)
    best_id, best_peak = ids[0], -np.inf
    for uid in ids:
        pix = cut[sub == uid]
        peak = float(np.nanmax(pix)) if pix.size else -np.inf
        if peak > best_peak:
            best_id, best_peak = uid, peak
    yy, xx = np.where(sub == best_id)
    if yy.size == 0:
        return None
    flux = np.clip(cut[yy, xx], 0, None)
    w = flux if float(np.sum(flux)) > 0 else np.ones(yy.size)
    peak_i = int(np.nanargmax(cut[yy, xx]))
    py, px = int(yy[peak_i]), int(xx[peak_i])
    peak = float(cut[py, px])
    core = (np.abs(yy - py) <= 3) & (np.abs(xx - px) <= 3) & (cut[yy, xx] >= 0.3 * peak)
    if not np.any(core):
        core = np.ones(yy.size, dtype=bool)
    cw = np.clip(cut[yy[core], xx[core]], 0, None)
    if float(np.sum(cw)) <= 0:
        cw = np.ones(int(np.count_nonzero(core)))
    return {
        "seg_id": int(best_id),
        "xc": float(np.average(xx, weights=w)) + 1.0,
        "yc": float(np.average(yy, weights=w)) + 1.0,
        "xc_core": float(np.average(xx[core], weights=cw)) + 1.0,
        "yc_core": float(np.average(yy[core], weights=cw)) + 1.0,
        "npix": int(yy.size),
    }


def _star_from_mask(wkdir: str) -> dict | None:
    """Fallback: brightest compact blob in host_mask (ignore huge no-data regions)."""
    from collections import deque

    mask_path = os.path.join(wkdir, "host_mask.fits")
    cut_path = os.path.join(wkdir, "host_cutout.fits")
    if not (os.path.isfile(mask_path) and os.path.isfile(cut_path)):
        return None
    mask = np.asarray(fits.getdata(mask_path)) > 0
    cut = np.asarray(fits.getdata(cut_path), float)
    ny, nx = mask.shape
    labeled = np.zeros(mask.shape, dtype=int)
    lab = 0
    for y in range(ny):
        for x in range(nx):
            if not mask[y, x] or labeled[y, x]:
                continue
            lab += 1
            q = deque([(y, x)])
            labeled[y, x] = lab
            while q:
                cy, cx = q.popleft()
                for dy, dx in ((-1, 0), (1, 0), (0, -1), (0, 1)):
                    ny2, nx2 = cy + dy, cx + dx
                    if (0 <= ny2 < ny and 0 <= nx2 < nx
                            and mask[ny2, nx2] and not labeled[ny2, nx2]):
                        labeled[ny2, nx2] = lab
                        q.append((ny2, nx2))
    n_pix = mask.size
    best_id, best_peak = 0, -np.inf
    for uid in range(1, lab + 1):
        npix = int(np.count_nonzero(labeled == uid))
        if npix < 4 or npix > 0.25 * n_pix:
            continue
        pix = cut[labeled == uid]
        peak = float(np.nanmax(pix)) if pix.size else -np.inf
        if peak > best_peak:
            best_id, best_peak = uid, peak
    if best_id == 0:
        return None
    yy, xx = np.where(labeled == best_id)
    flux = np.clip(cut[yy, xx], 0, None)
    w = flux if float(np.sum(flux)) > 0 else np.ones(yy.size)
    return {
        "seg_id": None,
        "xc": float(np.average(xx, weights=w)) + 1.0,
        "yc": float(np.average(yy, weights=w)) + 1.0,
        "npix": int(yy.size),
        "mask_label": int(best_id),
        "labeled": labeled,
    }


def _unmask_star(wkdir: str, star: dict, prod_dir: str, meta: dict) -> None:
    mask_path = os.path.join(wkdir, "host_mask.fits")
    with fits.open(mask_path, mode="update") as hdul:
        mask = np.asarray(hdul[0].data)
        if star.get("seg_id") is not None:
            bounds = _cutout_bounds(meta)
            seg = np.asarray(fits.getdata(os.path.join(prod_dir, "segmentation_map.fits")))
            while seg.ndim > 2:
                seg = seg[0]
            xmin, xmax, ymin, ymax = bounds
            sub = seg[ymin:ymax, xmin:xmax]
            sel = sub == int(star["seg_id"])
        else:
            sel = star["labeled"] == int(star["mask_label"])
        sigma_path = os.path.join(wkdir, "host_sigma.fits")
        if os.path.isfile(sigma_path):
            sigma = np.asarray(fits.getdata(sigma_path), float)
            finite = np.isfinite(sigma)
            huge = finite & (sigma >= 1.0e20)
            sel = sel & ~huge
        mask = np.asarray(mask, dtype=np.int16)
        mask[sel] = 0
        hdul[0].data = mask
        hdul.flush()


def _sky_component_index(lines: list[str]) -> int:
    for i, line in enumerate(lines):
        head = _COMP_HEADER.match(line)
        if not head:
            continue
        window = "\n".join(lines[i:i + 4]).lower()
        if "sky" in window:
            return i
    raise ValueError("no sky component in feedme")


def _renumber_components_from(lines: list[str], sky_i: int) -> list[str]:
    out = list(lines)
    for i, line in enumerate(out):
        head = _COMP_HEADER.match(line)
        if not head:
            continue
        cid = int(head.group(1))
        if i >= sky_i:
            out[i] = re.sub(
                r"(#\s*Component number:\s*)\d+",
                rf"\g<1>{cid + 1}",
                line,
                count=1,
            )
    return out


def _insert_psf_before_sky(lines: list[str], xc: float, yc: float,
                           mag: float, *, fix_xy: bool = False) -> list[str]:
    """Add a PSF component immediately before the sky block."""
    sky_i = _sky_component_index(lines)
    out = _renumber_components_from(lines, sky_i)
    psf_id = int(_COMP_HEADER.match(lines[sky_i]).group(1))
    xy_flag = "0 0" if fix_xy else "1 1"
    block = [
        f"# Component number: {psf_id}",
        " 0) psf  # Component type",
        f" 1) {xc:.4f} {yc:.4f} {xy_flag}  # position x y",
        f" 3) {mag:.4f} 1  # Integrated magnitude",
        " Z) 0  # Skip this model",
        "",
    ]
    return out[:sky_i] + block + out[sky_i:]


def _insert_sersic_before_sky(
    lines: list[str],
    xc: float,
    yc: float,
    mag: float,
    re_px: float,
    n: float = 1.0,
    q: float = 0.9,
    pa: float = 0.0,
) -> list[str]:
    """Add a free Sérsic immediately before the sky block (same coords as PSF)."""
    sky_i = _sky_component_index(lines)
    out = _renumber_components_from(lines, sky_i)
    cid = int(_COMP_HEADER.match(lines[sky_i]).group(1))
    block = [
        f"# Component number: {cid}",
        " 0) sersic  # Component type",
        f" 1) {xc:.4f} {yc:.4f} 1 1  # position x y",
        f" 3) {mag:.4f} 1  # Integrated magnitude",
        f" 4) {re_px:.4f} 1  # effective radius (pix)",
        f" 5) {n:.4f} 1  # sersic index",
        " 6) 0.0000 0  # ----",
        " 7) 0.0000 0  # ----",
        " 8) 0.0000 0  # ----",
        f" 9) {q:.4f} 1  # Axis ratio (b/a)",
        f"10) {pa:.4f} 1  # Position angle",
        " Z) 0  # Skip this model",
        "",
    ]
    return out[:sky_i] + block + out[sky_i:]


def _insert_moffat_before_sky(
    lines: list[str],
    xc: float,
    yc: float,
    mag: float,
    fwhm: float,
    c: float = 3.5,
    q: float = 0.95,
    pa: float = 0.0,
) -> list[str]:
    """Add a free Moffat immediately before the sky block.

    GALFIT convolves this with the empirical PSF, so ``fwhm`` is extra width
    on top of PSFEx, not a replacement for it.
    """
    sky_i = _sky_component_index(lines)
    out = _renumber_components_from(lines, sky_i)
    cid = int(_COMP_HEADER.match(lines[sky_i]).group(1))
    block = [
        f"# Component number: {cid}",
        " 0) moffat  # Component type",
        f" 1) {xc:.4f} {yc:.4f} 1 1  # position x y",
        f" 3) {mag:.4f} 1  # Integrated magnitude",
        f" 4) {fwhm:.4f} 1  # FWHM (pix)",
        f" 5) {c:.4f} 1  # power-law / concentration",
        " 6) 0.0000 0  # ----",
        " 7) 0.0000 0  # ----",
        " 8) 0.0000 0  # ----",
        f" 9) {q:.4f} 1  # Axis ratio (b/a)",
        f"10) {pa:.4f} 1  # Position angle",
        " Z) 0  # Skip this model",
        "",
    ]
    return out[:sky_i] + block + out[sky_i:]


def _companion_mag_seed(cut: dict, host_mag: float, dmag: float,
                        prod_dir: str) -> float:
    """Calibrated MAG_AUTO if present, else host minus dmag."""
    mag = cut.get("mag_auto")
    if isinstance(mag, (int, float)) and math.isfinite(float(mag)):
        mag = float(mag)
        if 8.0 <= mag <= 40.0:
            return min(max(mag, 8.5), 39.0)
        if mag < 0:
            zp = None
            zpath = os.path.join(prod_dir, "zero_points.json")
            zpj = vc.read_json(zpath) if os.path.isfile(zpath) else None
            if isinstance(zpj, dict):
                zp = (zpj.get("zp_auto") or zpj.get("zp_aper") or zpj.get("zp")
                      or (zpj.get("aper") or {}).get("zp"))
            if zp is None:
                summ = vc.read_json(os.path.join(prod_dir, "pipeline_summary.json")) or {}
                zp = (summ.get("photometry") or {}).get("zp_aper")
            if zp is not None and math.isfinite(float(zp)):
                cal = mag + float(zp)
                if 8.0 <= cal <= 40.0:
                    return min(max(cal, 8.5), 39.0)
    return min(max(float(host_mag) - float(dmag), 8.5), 39.0)


def _sersic_seed_from_star_cut(cut: dict, host_mag: float, dmag: float) -> dict:
    mag = cut.get("mag_auto")
    if (mag is None or not math.isfinite(float(mag))
            or float(mag) < 8.0 or float(mag) > 40.0):
        mag = float(host_mag) - float(dmag)
    mag = min(max(float(mag), 8.5), 39.0)
    re_px = cut.get("flux_radius_px")
    if re_px is None or not math.isfinite(float(re_px)) or float(re_px) <= 0:
        re_px = 2.0
    re_px = min(max(float(re_px), 1.5), 20.0)
    elong = cut.get("elongation")
    if elong is not None and math.isfinite(float(elong)) and float(elong) > 0:
        q = 1.0 / float(elong)
    else:
        q = 0.9
    q = min(max(q, 0.1), 1.0)
    theta = cut.get("theta_image")
    pa = (float(theta) - 90.0) if theta is not None and math.isfinite(float(theta)) else 0.0
    return {"mag": float(mag), "re": float(re_px), "n": 1.0, "q": float(q), "pa": float(pa)}


def _fmt_num(val, fmt: str = "{:.4f}") -> str:
    if val is None:
        return "-"
    try:
        x = float(val)
    except (TypeError, ValueError):
        return str(val)
    if not math.isfinite(x):
        return "-"
    return fmt.format(x)


def print_star_cut(cut: dict) -> None:
    """Stdout dump of the Phase 3a star/spread decision for the unmasked island."""
    print(
        f"  star-cut  seg={cut.get('seg_id')}  reason={cut.get('mask_reason')}  "
        f"CLASS_STAR={_fmt_num(cut.get('class_star'))}  "
        f"(mask if >= {_fmt_num(cut.get('class_star_threshold'), '{:.2f}')})  "
        f"stellar={cut.get('stellar_by_class_star')}"
    )
    print(
        f"  SPREAD_MODEL={_fmt_num(cut.get('spread_model'), '{:.5f}')} +/- "
        f"{_fmt_num(cut.get('spreaderr_model'), '{:.5f}')}  "
        f"SPREAD+3sig={_fmt_num(cut.get('spread_plus_3sig'), '{:.5f}')}  "
        f"(mask if < {_fmt_num(cut.get('spread_star_max'), '{:.3f}')})  "
        f"point_like={cut.get('point_like_by_spread')}  "
        f"match={_fmt_num(cut.get('spread_match_arcsec'), '{:.3f}')}\""
    )
    print(
        f"  PSFEx FWHM={_fmt_num(cut.get('psfex_fwhm_px'), '{:.3f}')} px  "
        f"FWHM_IMAGE={_fmt_num(cut.get('fwhm_image_px'), '{:.3f}')} px  "
        f"FLUX_RADIUS={_fmt_num(cut.get('flux_radius_px'), '{:.3f}')} px  "
        f"MAG_AUTO={_fmt_num(cut.get('mag_auto'), '{:.3f}')}  "
        f"PSFEx chi2={_fmt_num(cut.get('psfex_chi2'), '{:.2f}')}"
    )


def stage_refit(
    frb: str,
    label: str,
    *,
    sky_adu: float | None = None,
    fix_n: float | None = None,
    add_psf_masked_star: bool = False,
    add_sersic_masked_star: bool = False,
    add_moffat_masked_star: bool = False,
    psf_dmag: float = 2.0,
    fix_psf_xy: bool = False,
    fix_host_xy: bool = False,
) -> str:
    """Build ``Re-fits/<FRB>/<label>/`` ready for GALFIT + verification."""
    prod = vc.load_host(frb)
    wkdir = os.path.join(REFITS_ROOT, frb, label)
    if os.path.isdir(wkdir):
        shutil.rmtree(wkdir)
    os.makedirs(wkdir, exist_ok=True)

    for name in _SIDECARS:
        src = os.path.join(prod.dir, name)
        if os.path.isfile(src):
            shutil.copy2(src, os.path.join(wkdir, name))

    hdr = vc.parse_out_header(os.path.join(prod.dir, "out.fits"))
    comps = hdr["components"]
    sky_comp = int(hdr["sky"].get("comp", len(comps) + 1))

    feed = vc.parse_feedme(os.path.join(prod.dir, "galfit.feedme"))
    # If sky is left free, reseed it at the production best-fit (still flag 1).
    prod_sky = float(hdr["sky"].get("level", float("nan")))
    sky_seed = prod_sky if (sky_adu is None and math.isfinite(prod_sky)) else None
    edited = _edit_feedme_fixed(
        feed["lines"], comps,
        sky_value=sky_adu, sky_comp=sky_comp, fix_n=fix_n, reseed=True,
        sky_seed_free=sky_seed, fix_host_xy=fix_host_xy,
    )
    with open(os.path.join(wkdir, "galfit.feedme"), "w", newline="\n",
              encoding="utf-8") as f:
        f.write("\n".join(edited) + "\n")

    # Sky fixed → drop sky constraints; n fixed → drop n bounds row.
    drop: set[str] = set()
    if sky_adu is not None:
        # sky_perturb helper drops by component id; here rewrite constraints.
        csrc = os.path.join(prod.dir, "constraints.txt")
        try:
            with open(csrc, encoding="utf-8", errors="replace") as f:
                clines = f.read().splitlines()
        except OSError:
            clines = []
        kept = []
        for ln in clines:
            toks = ln.split()
            if not toks:
                kept.append(ln)
                continue
            if sky_adu is not None and toks[0] == str(sky_comp):
                continue
            if fix_n is not None and len(toks) >= 2 and toks[1].lower() == "n":
                continue
            kept.append(ln)
        with open(os.path.join(wkdir, "constraints.txt"), "w", newline="\n",
                  encoding="utf-8") as f:
            f.write("\n".join(kept) + ("\n" if kept else ""))
    elif fix_n is not None:
        _strip_param_constraints(
            os.path.join(prod.dir, "constraints.txt"),
            os.path.join(wkdir, "constraints.txt"),
            {"n"},
        )

    n_comp = int(add_psf_masked_star) + int(add_sersic_masked_star) + int(add_moffat_masked_star)
    if n_comp > 1:
        raise SystemExit("choose one companion: PSF, Sérsic, or Moffat")

    psf_info: dict[str, Any] | None = None
    sersic_info: dict[str, Any] | None = None
    moffat_info: dict[str, Any] | None = None
    star_cut: dict[str, Any] | None = None
    if n_comp:
        meta_cut = vc.read_json(os.path.join(wkdir, "cutout_meta.json")) or {}
        star = _star_from_segmap(prod.dir, meta_cut) or _star_from_mask(wkdir)
        if star is None:
            raise SystemExit(f"{frb}: no masked star island found to add a companion")
        host_mag = float(prod.mag)
        star_cut = masked_object_star_cut(prod.dir, star.get("seg_id"))
        print_star_cut(star_cut)
        _unmask_star(wkdir, star, prod.dir, meta_cut)
        feed_path = os.path.join(wkdir, "galfit.feedme")
        with open(feed_path, encoding="utf-8", errors="replace") as f:
            flines = f.read().splitlines()
        mag_seed = _companion_mag_seed(star_cut, host_mag, psf_dmag, prod.dir)
        xc_star = float(star.get("xc_core") or star["xc"])
        yc_star = float(star.get("yc_core") or star["yc"])
        if not fix_psf_xy:
            xc_star, yc_star = float(star["xc"]), float(star["yc"])
        print(
            f"  companion xy=({xc_star:.3f},{yc_star:.3f})  "
            f"island=({star['xc']:.3f},{star['yc']:.3f})  "
            f"xy_fixed={bool(fix_psf_xy)}  mag_seed={mag_seed:.3f}"
        )
        elong = star_cut.get("elongation")
        q_seed = 0.95
        if elong is not None and math.isfinite(float(elong)) and float(elong) > 0:
            q_seed = min(max(1.0 / float(elong), 0.1), 1.0)
        theta = star_cut.get("theta_image")
        pa_seed = ((float(theta) - 90.0)
                   if theta is not None and math.isfinite(float(theta)) else 0.0)
        pa_seed = vc.wrap_pa(pa_seed)
        cpath = os.path.join(wkdir, "constraints.txt")
        if add_sersic_masked_star:
            seed = _sersic_seed_from_star_cut(star_cut, host_mag, psf_dmag)
            flines = _insert_sersic_before_sky(
                flines, xc_star, yc_star,
                seed["mag"], seed["re"], seed["n"], seed["q"], seed["pa"],
            )
            with open(cpath, "a", encoding="utf-8") as f:
                f.write("2 n 0.5 to 6.0\n")
                f.write("2 re 1.5 to 100.0\n")
                f.write("2 mag 8.0 to 40.0\n")
            sersic_info = {
                "seg_id": star.get("seg_id"),
                "xc": star["xc"], "yc": star["yc"], "npix": star["npix"],
                "host_mag": host_mag, **seed, "star_cut": star_cut,
            }
        elif add_moffat_masked_star:
            fwhm_seed = 1.5
            c_seed = 3.5
            flines = _insert_moffat_before_sky(
                flines, xc_star, yc_star, mag_seed,
                fwhm_seed, c_seed, q_seed, pa_seed,
            )
            with open(cpath, "a", encoding="utf-8") as f:
                f.write("2 mag 8.0 to 40.0\n")
            moffat_info = {
                "seg_id": star.get("seg_id"),
                "xc": star["xc"], "yc": star["yc"], "npix": star["npix"],
                "host_mag": host_mag, "mag": mag_seed,
                "fwhm": fwhm_seed, "c": c_seed, "q": q_seed, "pa": pa_seed,
                "star_cut": star_cut,
            }
        else:
            flines = _insert_psf_before_sky(
                flines, xc_star, yc_star, mag_seed, fix_xy=fix_psf_xy,
            )
            with open(cpath, "a", encoding="utf-8") as f:
                f.write("2 mag 8.0 to 40.0\n")
            psf_info = {
                "seg_id": star.get("seg_id"),
                "xc": xc_star, "yc": yc_star, "npix": star["npix"],
                "xc_island": star["xc"], "yc_island": star["yc"],
                "xy_fixed": bool(fix_psf_xy),
                "mag_seed": mag_seed, "dmag": float(psf_dmag),
                "host_mag": host_mag, "star_cut": star_cut,
            }
        with open(feed_path, "w", newline="\n", encoding="utf-8") as f:
            f.write("\n".join(flines) + "\n")

    meta = {
        "frb": frb,
        "label": label,
        "sky_fixed_adu": sky_adu,
        "n_fixed": fix_n,
        "host_xy_fixed": bool(fix_host_xy),
        "psf_added": bool(psf_info),
        "psf_xy_fixed": bool(fix_psf_xy),
        "sersic_added": bool(sersic_info),
        "moffat_added": bool(moffat_info),
        "psf": psf_info,
        "sersic": sersic_info,
        "moffat": moffat_info,
        "star_cut": star_cut,
        "production_sky_adu": prod.sky_level,
        "production": {
            "q": prod.q, "n": prod.n, "re": prod.re,
            "mag": prod.mag, "pa": prod.pa,
        },
    }
    vc.write_json(os.path.join(wkdir, "refit_meta.json"), meta)
    return wkdir


def run_galfit_refit(wkdir: str) -> dict[str, Any]:
    run_galfit = sp._run_galfit()
    for stale in ("fit.log", "out.fits"):
        path = os.path.join(wkdir, stale)
        if os.path.isfile(path):
            os.remove(path)
    ok = False
    try:
        with open(os.path.join(wkdir, "galfit_stdout.log"), "w",
                  encoding="utf-8") as fh, contextlib.redirect_stdout(fh):
            ok = bool(run_galfit(wkdir))
    except Exception as exc:
        return {"status": f"launch_failed: {type(exc).__name__}: {exc}",
                "converged": False}
    fit = sp._read_result(wkdir)
    fit["converged"] = ok
    return fit


def production_panel_src(frb: str) -> str | None:
    """Canonical production panel — ``outputs/panels/<FRB>.png`` only."""
    path = os.path.join(vc.OUT_ROOT, "panels", f"{frb}.png")
    return path if os.path.isfile(path) else None


def copy_production_panel(frb: str, dest_path: str) -> str | None:
    """Byte-copy the production panel. Never regenerate it."""
    src = production_panel_src(frb)
    if src is None:
        return None
    os.makedirs(os.path.dirname(os.path.abspath(dest_path)) or ".", exist_ok=True)
    shutil.copy2(src, dest_path)
    return dest_path


def run_refit(
    frb: str,
    *,
    label: str,
    sky_adu: float | None = None,
    fix_n: float | None = None,
    checks: list[str] | None = None,
    force: bool = True,
    add_psf_masked_star: bool = False,
    add_sersic_masked_star: bool = False,
    add_moffat_masked_star: bool = False,
    psf_dmag: float = 2.0,
    fix_psf_xy: bool = False,
    fix_host_xy: bool = False,
) -> dict[str, Any]:
    """Stage → GALFIT → full verification suite for one re-fit label."""
    checks = checks or list(rv.CHECK_ORDER)
    wkdir = stage_refit(
        frb, label, sky_adu=sky_adu, fix_n=fix_n,
        add_psf_masked_star=add_psf_masked_star,
        add_sersic_masked_star=add_sersic_masked_star,
        add_moffat_masked_star=add_moffat_masked_star,
        psf_dmag=psf_dmag,
        fix_psf_xy=fix_psf_xy,
        fix_host_xy=fix_host_xy,
    )
    fit = run_galfit_refit(wkdir)
    out_fits = os.path.join(wkdir, "out.fits")
    components: list[dict] = []
    sky: dict = {}
    if os.path.isfile(out_fits):
        try:
            hdr = vc.parse_out_header(out_fits)
            components = hdr.get("components") or []
            sky = hdr.get("sky") or {}
        except Exception:
            pass
    if fit.get("status") != "ok":
        return {"frb": frb, "label": label, "wkdir": wkdir, "fit": fit,
                "components": components, "sky": sky, "verification": None}

    t0 = time.time()
    ver = rv.run_host_dir(frb, wkdir, wkdir, checks, force=force)

    summary = {
        "frb": frb,
        "label": label,
        "wkdir": wkdir,
        "sky_fixed_adu": sky_adu,
        "n_fixed": fix_n,
        "host_xy_fixed": fix_host_xy,
        "psf_added": add_psf_masked_star,
        "psf_xy_fixed": fix_psf_xy,
        "sersic_added": add_sersic_masked_star,
        "moffat_added": add_moffat_masked_star,
        "fit": fit,
        "components": components,
        "sky": sky,
        "verification": ver,
        "verification_runtime_s": round(time.time() - t0, 3),
        "panel_png": os.path.join(wkdir, "panel.png"),
    }
    panel_src = os.path.join(wkdir, "panel.png")
    if ((add_psf_masked_star or add_sersic_masked_star or add_moffat_masked_star)
            and os.path.isfile(panel_src)):
        pub = os.path.join(vc.OUT_ROOT, "panels", f"{frb}_{label}.png")
        os.makedirs(os.path.dirname(pub), exist_ok=True)
        shutil.copy2(panel_src, pub)
        summary["published_panel"] = pub
    vc.write_json(os.path.join(wkdir, "refit_summary.json"), summary)
    return summary


def _sky_from_protocol(frb: str) -> float:
    path = os.path.join(REFITS_ROOT, frb, "sky_protocol.json")
    data = vc.read_json(path)
    sky = (data.get("consensus") or {}).get("sky_adu")
    if sky is None or not math.isfinite(float(sky)):
        raise SystemExit(f"no consensus sky in {path}; run sky_protocol.py first")
    return float(sky)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("frb")
    ap.add_argument("--label", default=None,
                    help="subdir under Re-fits/<FRB>/ (default derived from flags)")
    ap.add_argument("--sky-adu", type=float, default=None)
    ap.add_argument("--sky-from-protocol", action="store_true",
                    help="use Re-fits/<FRB>/sky_protocol.json consensus sky")
    ap.add_argument("--fix-n", type=float, default=None,
                    help="hold host Sérsic index fixed at this value")
    ap.add_argument("--add-psf-at-masked-star", action="store_true",
                    help="unmask the Phase 3a star island and add a free PSF")
    ap.add_argument("--fix-psf-xy", action="store_true",
                    help="hold the PSF (or companion) at the star-core photocenter")
    ap.add_argument("--fix-host-xy", action="store_true",
                    help="hold the host Sérsic at the production photocenter")
    ap.add_argument("--add-sersic-at-masked-star", action="store_true",
                    help="unmask the same island and add a free Sérsic instead of a PSF")
    ap.add_argument("--add-moffat-at-masked-star", action="store_true",
                    help="unmask the same island and add a free Moffat (PSF-convolved)")
    ap.add_argument("--psf-dmag", type=float, default=2.0,
                    help="companion mag seed = host mag minus this if MAG_AUTO missing (default 2)")
    ap.add_argument("--checks", nargs="+", default=["all"],
                    choices=["all", *rv.CHECK_ORDER])
    ap.add_argument("--force", action="store_true", default=True)
    args = ap.parse_args(argv)

    sky = args.sky_adu
    if args.sky_from_protocol:
        sky = _sky_from_protocol(args.frb)

    parts = []
    if sky is not None:
        parts.append("sky")
    if args.fix_n is not None:
        parts.append(f"n{args.fix_n:g}")
    if args.add_psf_at_masked_star:
        parts.append("psf")
    if args.fix_psf_xy:
        parts.append("xy")
    if args.fix_host_xy:
        parts.append("hostxy")
    if args.add_sersic_at_masked_star:
        parts.append("sersic")
    if args.add_moffat_at_masked_star:
        parts.append("moffat")
    label = args.label or ("_".join(parts) if parts else "refit")

    checks = rv.CHECK_ORDER if "all" in args.checks else [
        c for c in rv.CHECK_ORDER if c in args.checks
    ]
    print(f"Re-fit {args.frb}  label={label}  sky={sky}  fix_n={args.fix_n}  "
          f"psf_star={args.add_psf_at_masked_star}  "
          f"fix_psf_xy={args.fix_psf_xy}  "
          f"sersic_star={args.add_sersic_at_masked_star}  "
          f"moffat_star={args.add_moffat_at_masked_star}")
    print(f"  checks: {', '.join(checks)}")
    summary = run_refit(
        args.frb, label=label, sky_adu=sky, fix_n=args.fix_n,
        checks=checks, force=True,
        add_psf_masked_star=args.add_psf_at_masked_star,
        add_sersic_masked_star=args.add_sersic_at_masked_star,
        add_moffat_masked_star=args.add_moffat_at_masked_star,
        psf_dmag=args.psf_dmag,
        fix_psf_xy=args.fix_psf_xy,
        fix_host_xy=args.fix_host_xy,
    )
    fit = summary.get("fit") or {}
    print(f"  GALFIT: status={fit.get('status')}  q={fit.get('q')}  "
          f"n={fit.get('n')}  Re={fit.get('re')}  m={fit.get('mag')}  "
          f"sky={fit.get('sky')}")
    for c in summary.get("components") or []:
        print(
            f"  comp {c.get('comp')} {c.get('type')}: "
            f"xy=({_fmt_num(c.get('xc'))},{_fmt_num(c.get('yc'))})  "
            f"m={_fmt_num(c.get('mag'))}  Re={_fmt_num(c.get('re'))}  "
            f"n={_fmt_num(c.get('n'))}  FWHM={_fmt_num(c.get('fwhm'))}  "
            f"c={_fmt_num(c.get('c'))}  q={_fmt_num(c.get('q'))}"
        )
    ver = summary.get("verification") or {}
    print(f"  verification: {ver.get('status')}")
    print(f"  panel: {summary.get('panel_png')}")
    if summary.get("published_panel"):
        print(f"  published: {summary['published_panel']}")
    return 0 if fit.get("status") == "ok" else 1


if __name__ == "__main__":
    raise SystemExit(main())
