import os
import sys
import argparse
import json

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.dirname(_HERE))
from pipeline_shared import get_logger  # noqa: E402
from sersic_init import effective_re_px  # noqa: E402

log = get_logger("phase3a")

def get_table_from_ldac(filename, frame=1):
    if frame > 0: frame = frame * 2
    return Table.read(filename, hdu=frame, format='fits')

def _pad_bbox(xmin, xmax, ymin, ymax, pad, shape):
    """Pad an [xmin,xmax,ymin,ymax) bbox by `pad` px, clipped to image shape."""
    ny, nx = shape
    return (
        int(max(xmin - pad, 0)),
        int(min(xmax + pad, nx)),
        int(max(ymin - pad, 0)),
        int(min(ymax + pad, ny)),
    )


def _clamp_bbox_around_center(xmin, xmax, ymin, ymax, cx, cy, max_side, shape):
    """If either side exceeds max_side, replace with a max_side window centred on (cx, cy)."""
    if max_side <= 0:
        return xmin, xmax, ymin, ymax
    ny_img, nx_img = shape
    nx, ny = xmax - xmin, ymax - ymin
    if nx <= max_side and ny <= max_side:
        return xmin, xmax, ymin, ymax
    half = max_side // 2
    icx, icy = int(round(float(cx))), int(round(float(cy)))
    xmin = max(0, icx - half)
    ymin = max(0, icy - half)
    xmax = min(nx_img, xmin + max_side)
    ymax = min(ny_img, ymin + max_side)
    if xmax - xmin < max_side:
        xmin = max(0, xmax - max_side)
    if ymax - ymin < max_side:
        ymin = max(0, ymax - max_side)
    return xmin, xmax, ymin, ymax


def _bbox_of_objids(seg_map, objids):
    """Tight (xmin, xmax, ymin, ymax) covering every pixel of the given IDs."""
    if not objids:
        return None
    y, x = np.where(np.isin(seg_map, list(objids)))
    if len(x) == 0 or len(y) == 0:
        return None
    return int(x.min()), int(x.max() + 1), int(y.min()), int(y.max() + 1)


def _ids_touching_roi(seg_map, xmin, xmax, ymin, ymax, drop):
    """Seg IDs that have at least one pixel inside ROI [xmin:xmax, ymin:ymax)."""
    sub = seg_map[ymin:ymax, xmin:xmax]
    ids = set(int(v) for v in np.unique(sub))
    ids.discard(0)
    for d in drop:
        ids.discard(int(d))
    return ids


def _class_star_for_number(cat, number: int):
    """Return CLASS_STAR for this SExtractor NUMBER, or None if missing."""
    sel = cat["NUMBER"] == number
    if not np.any(sel) or "CLASS_STAR" not in cat.colnames:
        return None
    return float(cat["CLASS_STAR"][sel][0])


def _is_stellar(cat, uid: int, class_star_max: float) -> bool:
    """True -> treat as point source: mask, do not add Sérsic."""
    cs = _class_star_for_number(cat, uid)
    if cs is None:
        return False
    return cs >= class_star_max


def _pixel_dist_from_host(cat, host_uid: int, uid: int) -> float:
    """Separation in pixels between two SExtractor NUMBER entries."""
    hsel = cat["NUMBER"] == int(host_uid)
    usel = cat["NUMBER"] == int(uid)
    if not np.any(hsel) or not np.any(usel):
        return float("inf")
    hx = float(cat["X_IMAGE"][hsel][0])
    hy = float(cat["Y_IMAGE"][hsel][0])
    ux = float(cat["X_IMAGE"][usel][0])
    uy = float(cat["Y_IMAGE"][usel][0])
    return float(np.hypot(ux - hx, uy - hy))


def _catalog_row_for_number(cat, number: int):
    """Return the first astropy table row for SExtractor NUMBER, or None."""
    sel = cat["NUMBER"] == int(number)
    if not np.any(sel):
        return None
    return cat[sel][0]


def _grow_roi_for_fit_ids(seg_map, fit_ids, host_pad: int, max_side: int, shape):
    """Pad union bbox of fit_ids; return new bounds or None if growth blocked."""
    grow_bbox = _bbox_of_objids(seg_map, list(fit_ids))
    if grow_bbox is None:
        return None
    gxmin, gxmax, gymin, gymax = _pad_bbox(
        *grow_bbox, pad=int(host_pad), shape=shape
    )
    nx, ny = gxmax - gxmin, gymax - gymin
    if max_side > 0 and (nx > max_side or ny > max_side):
        return None
    return gxmin, gxmax, gymin, gymax


def resolve_neighbor_re_roi(
    seg_map,
    cat,
    host_uid: int,
    *,
    host_pad: int = 20,
    re_sep_factor: float = 3.0,
    max_roi_iterations: int = 8,
    max_cutout_side: int = 512,
    neighbor_class_star_max: float = 0.75,
    spread_by_number=None,
    catalog_path: str = "",
    psf_match_arcsec: float = 0.5,
):
    """Build ROI + fit/mask sets via host_pad and re_sep_factor × Re_neighbor.

    For every non-host seg island that clips the ROI:
      - stars / point-like SPREAD → mask (never expand)
      - sep(host, neighbor) > re_sep_factor × Re_neighbor → mask
      - else → jointly fit; grow ROI so all fit members have host_pad around them

    Returns
    -------
    xmin, xmax, ymin, ymax, fit_objids, mask_objids, decisions
        decisions maps uid -> (kind, detail) for logging/tests.
        kind in {"host", "fit", "mask_far", "mask_star", "mask_spread", "mask_nocat"}.
    """
    host_uid = int(host_uid)
    host_bbox = _bbox_of_objids(seg_map, [host_uid])
    if host_bbox is None:
        raise ValueError(f"Host NUMBER={host_uid} not in segmentation map")

    xmin, xmax, ymin, ymax = _pad_bbox(
        *host_bbox, pad=int(host_pad), shape=seg_map.shape
    )
    fit_objids = {host_uid}
    mask_objids: set = set()
    decisions: dict = {host_uid: ("host", 1.0)}
    if spread_by_number is None:
        spread_by_number = {}

    max_iters = max(1, int(max_roi_iterations))
    max_side = int(max_cutout_side)
    factor = float(re_sep_factor)

    for it in range(1, max_iters + 1):
        touching = _ids_touching_roi(
            seg_map, xmin, xmax, ymin, ymax, drop=[]
        )
        new_expanded: set = set()
        for uid in sorted(touching):
            uid = int(uid)
            if uid == host_uid:
                continue
            if uid in fit_objids or uid in mask_objids:
                continue

            stellar = _is_stellar(cat, uid, float(neighbor_class_star_max))
            spread = _spread_for_seg_number(
                uid, cat, spread_by_number, catalog_path, float(psf_match_arcsec)
            )
            point_like = (
                spread is not None
                and _spread_is_point_source(spread[0], spread[1])
            )
            if stellar or point_like:
                mask_objids.add(uid)
                decisions[uid] = (
                    "mask_star" if stellar else "mask_spread",
                    0.0,
                )
                continue

            row = _catalog_row_for_number(cat, uid)
            if row is None:
                mask_objids.add(uid)
                decisions[uid] = ("mask_nocat", 0.0)
                continue

            re_n = effective_re_px(row)
            sep = _pixel_dist_from_host(cat, host_uid, uid)
            thresh = factor * re_n
            if sep > thresh:
                mask_objids.add(uid)
                decisions[uid] = ("mask_far", sep / re_n if re_n > 0 else float("inf"))
            else:
                fit_objids.add(uid)
                decisions[uid] = ("fit", sep / re_n if re_n > 0 else 0.0)
                new_expanded.add(uid)

        if not new_expanded:
            break

        grown = _grow_roi_for_fit_ids(
            seg_map, fit_objids, int(host_pad), max_side, seg_map.shape
        )
        if grown is None:
            # Cap / missing bbox: keep current ROI; demote brand-new expanders
            # that still aren't fully coverable under the side limit.
            log.warning(
                f"    [iter {it}] ROI growth blocked (max_cutout_side={max_side} "
                f"or empty bbox); new expanders {sorted(new_expanded)} stay "
                "in fit set with current bounds (overlapping pixels only)."
            )
            break
        xmin, xmax, ymin, ymax = grown
        log.info(
            f"    [iter {it}] expanded ROI for {sorted(new_expanded)} "
            f"-> bounds X=[{xmin}:{xmax}] Y=[{ymin}:{ymax}]"
        )

    # Any leftover touching IDs not yet classified → mask (safety net).
    for uid in _ids_touching_roi(seg_map, xmin, xmax, ymin, ymax, drop=[]):
        uid = int(uid)
        if uid in fit_objids or uid in mask_objids:
            continue
        mask_objids.add(uid)
        decisions[uid] = ("mask_far", -1.0)

    return xmin, xmax, ymin, ymax, fit_objids, mask_objids, decisions


# Same star/galaxy cut as Phase 2 (run_photometry_astropath.py AstroPath candidates).
_SPREAD_STAR_MAX = 0.005
_SPREAD_SIGMA = 3.0


def _spread_is_point_source(spread_model: float, spread_err: float) -> bool:
    return (spread_model + _SPREAD_SIGMA * spread_err) < _SPREAD_STAR_MAX


def _psf_catalog_path(catalog_path: str) -> str:
    return os.path.join(os.path.dirname(os.path.abspath(catalog_path)), "image.psf.cat")


def _load_psf_catalog(catalog_path: str):
    """Phase 2 SExtractor catalog (SPREAD, AstroPath sex_number)."""
    psf_path = _psf_catalog_path(catalog_path)
    if not os.path.isfile(psf_path):
        return None
    return get_table_from_ldac(psf_path)


def _load_spread_by_number(catalog_path: str) -> dict:
    """NUMBER -> (SPREAD_MODEL, SPREADERR_MODEL) from image.psf.cat (Phase 2)."""
    psf = _load_psf_catalog(catalog_path)
    if psf is None or "SPREAD_MODEL" not in psf.colnames or "SPREADERR_MODEL" not in psf.colnames:
        return {}
    out = {}
    for row in psf:
        out[int(row["NUMBER"])] = (float(row["SPREAD_MODEL"]), float(row["SPREADERR_MODEL"]))
    return out


def _spread_for_catalog_index(idx: int, cat, spread_by_number: dict, catalog_path: str,
                              match_arcsec: float = 0.5):
    """SPREAD for a Phase-1 image.cat row from Phase-2 image.psf.cat via sky position.

    Phase-1 NUMBER and Phase-2 NUMBER are independent SExtractor indices; never
    key solely by NUMBER (see _seg_number_for_psf_number).
    """
    psf = _load_psf_catalog(catalog_path)
    if psf is None:
        return None
    coord = SkyCoord(
        ra=cat["ALPHAWIN_J2000"][idx],
        dec=cat["DELTAWIN_J2000"][idx],
        unit="deg",
        frame="icrs",
    )
    psf_coords = SkyCoord(
        ra=psf["ALPHAWIN_J2000"], dec=psf["DELTAWIN_J2000"], unit="deg", frame="icrs"
    )
    seps_arcsec = coord.separation(psf_coords).arcsec
    j = int(np.argmin(seps_arcsec))
    if float(seps_arcsec[j]) > match_arcsec:
        return None
    return spread_by_number.get(int(psf["NUMBER"][j]))


def _spread_for_seg_number(
    seg_number: int,
    cat,
    spread_by_number: dict,
    catalog_path: str,
    match_arcsec: float = 0.5,
):
    """SPREAD for a Phase-1 image.cat / segmentation NUMBER (sky-match to Phase 2)."""
    sel = cat["NUMBER"] == int(seg_number)
    if not np.any(sel):
        return None
    idx = int(np.where(sel)[0][0])
    return _spread_for_catalog_index(
        idx, cat, spread_by_number, catalog_path, match_arcsec
    )


def _seg_number_for_psf_number(psf_number: int, cat, catalog_path: str, max_sep_arcsec: float = 1.0):
    """Map AstroPath / Phase-2 NUMBER to Phase-1 segmentation NUMBER via sky position.

    image.psf.cat (Phase 2) and image.cat (Phase 1) NUMBER fields are independent
    SExtractor run indices — never assume they refer to the same source.
    """
    psf = _load_psf_catalog(catalog_path)
    if psf is None:
        raise SystemExit(
            f"[Phase 3a] AstroPath sex_number={psf_number} not in image.cat and "
            "image.psf.cat is missing."
        )
    sel = psf["NUMBER"] == int(psf_number)
    if not np.any(sel):
        raise SystemExit(
            f"[Phase 3a] AstroPath sex_number={psf_number} not found in image.psf.cat."
        )
    host = SkyCoord(
        ra=psf["ALPHAWIN_J2000"][sel][0],
        dec=psf["DELTAWIN_J2000"][sel][0],
        unit="deg",
        frame="icrs",
    )
    cat_coords = SkyCoord(
        ra=cat["ALPHAWIN_J2000"], dec=cat["DELTAWIN_J2000"], unit="deg", frame="icrs"
    )
    seps_arcsec = host.separation(cat_coords).arcsec
    idx = int(np.argmin(seps_arcsec))
    sep = float(seps_arcsec[idx])
    if sep > max_sep_arcsec:
        raise SystemExit(
            f"[Phase 3a] AstroPath sex_number={psf_number} has no image.cat match within "
            f"{max_sep_arcsec}\" (nearest {sep:.3f}\")."
        )
    seg_number = int(cat["NUMBER"][idx])
    if seg_number != int(psf_number):
        log.info(
            f"AstroPath sex_number={psf_number} (image.psf.cat) -> "
            f"seg NUMBER={seg_number} (image.cat; {sep:.3f}\")"
        )
    return seg_number


def _pick_host_at_target_pixel(
    w: WCS,
    target_coord: SkyCoord,
    seg_map,
    cat,
    spread_by_number: dict,
    catalog_path: str,
    psf_match_arcsec: float = 0.5,
):
    """Return image.cat NUMBER under the target sky position (preferred for CSV hosts)."""
    xp, yp = w.world_to_pixel(target_coord)
    ix = int(round(float(np.squeeze(xp))))
    iy = int(round(float(np.squeeze(yp))))
    ny, nx = seg_map.shape
    if not (0 <= ix < nx and 0 <= iy < ny):
        return None
    seg_id = int(seg_map[iy, ix])
    if seg_id <= 0:
        return None
    sel = cat["NUMBER"] == seg_id
    if not np.any(sel):
        return None
    idx = int(np.where(sel)[0][0])
    spread = _spread_for_catalog_index(
        idx, cat, spread_by_number, catalog_path, psf_match_arcsec
    )
    if spread is None:
        return None
    if _spread_is_point_source(*spread):
        return None
    host_coord = SkyCoord(
        ra=cat["ALPHAWIN_J2000"][idx],
        dec=cat["DELTAWIN_J2000"][idx],
        unit="deg",
        frame="icrs",
    )
    sep_arcsec = float(target_coord.separation(host_coord).arcsec)
    return int(idx), int(seg_id), sep_arcsec


def _pick_nearest_galaxy_host(
    target_coord: SkyCoord,
    cat,
    spread_by_number: dict,
    catalog_path: str,
    max_sep_arcsec: float,
    psf_match_arcsec: float = 0.5,
):
    """Nearest SExtractor source that passes the Phase 2 SPREAD star cut.

    Returns (catalog_index, objid, sep_arcsec) or raises SystemExit if none qualify.
    """
    if not spread_by_number:
        raise SystemExit(
            "[Phase 3a] Cannot apply galaxy SPREAD cut: image.psf.cat missing or lacks "
            "SPREAD_MODEL / SPREADERR_MODEL (run Phase 2 first)."
        )
    cat_coords = SkyCoord(
        ra=cat["ALPHAWIN_J2000"], dec=cat["DELTAWIN_J2000"], unit="deg", frame="icrs"
    )
    seps_arcsec = target_coord.separation(cat_coords).arcsec
    for idx in np.argsort(seps_arcsec):
        sep = float(seps_arcsec[idx])
        if sep > max_sep_arcsec:
            break
        num = int(cat["NUMBER"][idx])
        spread = _spread_for_catalog_index(idx, cat, spread_by_number, catalog_path, psf_match_arcsec)
        if spread is None:
            log.info(f"    skip #{num}: no SPREAD in image.psf.cat ({sep:.2f}\" from target)")
            continue
        sm, se = spread
        if _spread_is_point_source(sm, se):
            log.info(
                f"    skip #{num}: point source SPREAD={sm:.4f}+/-{se:.4f} "
                f"({sep:.2f}\" from target)"
            )
            continue
        return int(idx), num, sep
    raise SystemExit(
        f"[Phase 3a] No galaxy (SPREAD+3*SPREADERR >= {_SPREAD_STAR_MAX}) within "
        f"{max_sep_arcsec}\" of the target position. Cannot associate FRB to a star."
    )


def _resolve_target_from_astropath(posteriors_csv: str, frb_ra: float,
                                   frb_dec: float, min_posterior: float):
    """If Phase 2 produced AstroPath posteriors, return its best-host sky position.

    Returns (ra, dec, label, sex_number). ``sex_number`` is the SExtractor
    ``NUMBER`` of the AstroPath candidate (Phase 2 galaxy cut already applied).
    When ``sex_number`` is missing (older posteriors files), Phase 3a falls
    back to nearest-galaxy matching at the AstroPath RA/Dec.

    Falls back to the raw FRB coords when the file is absent, unreadable, or
    no candidate clears ``min_posterior``.
    """
    if not os.path.exists(posteriors_csv):
        return frb_ra, frb_dec, "FRB position (no astropath_posteriors.csv next to catalog)", None
    try:
        df = pd.read_csv(posteriors_csv)
    except Exception as exc:  # noqa: BLE001
        return frb_ra, frb_dec, f"FRB position (posteriors read failed: {exc})", None
    if "posterior_O" not in df.columns or len(df) == 0:
        return frb_ra, frb_dec, "FRB position (posteriors empty / wrong schema)", None
    best = df.sort_values("posterior_O", ascending=False).iloc[0]
    p_best = float(best["posterior_O"])
    if not np.isfinite(p_best) or p_best < min_posterior:
        return frb_ra, frb_dec, (
            f"FRB position (best AstroPath P(O)={p_best:.4f} < {min_posterior})"
        ), None
    sex_number = None
    if "sex_number" in best.index and pd.notna(best["sex_number"]):
        sex_number = int(best["sex_number"])
    label = (
        f"AstroPath host objid={best.get('objid', '?')}"
        + (f", NUMBER={sex_number}" if sex_number is not None else "")
        + f" (P(O)={p_best:.4f}, mag={float(best.get('mag', np.nan)):.2f})"
    )
    return float(best["ra"]), float(best["dec"]), label, sex_number


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--image", required=True, help="Input raw FITS image")
    parser.add_argument("--segmap", required=True, help="Input SEGMENTATION FITS image")
    parser.add_argument("--catalog", required=True, help="Input SExtractor Catalog (.cat LDAC)")
    parser.add_argument("--ra", type=float, required=True, help="FRB RA (deg). Overridden by AstroPath host if posteriors are present.")
    parser.add_argument("--dec", type=float, required=True, help="FRB Dec (deg). Overridden by AstroPath host if posteriors are present.")
    parser.add_argument("--outdir", default=".", help="Output directory")
    parser.add_argument("--astropath-posteriors", default=None,
                        help="Path to astropath_posteriors.csv. Defaults to "
                             "<catalog_dir>/astropath_posteriors.csv. If found "
                             "and the best P(O) clears --min-astropath-posterior, "
                             "that host's RA/Dec is used instead of --ra/--dec.")
    parser.add_argument("--min-astropath-posterior", type=float, default=0.05,
                        help="Minimum AstroPath posterior_O to trust over the "
                             "raw FRB position (default 0.05).")
    parser.add_argument("--no-astropath-override", action="store_true",
                        help="Ignore astropath_posteriors.csv and always use --ra/--dec.")
    parser.add_argument(
        "--max-host-sep-arcsec",
        type=float,
        default=5.0,
        help="Abort if the nearest SExtractor source to the target RA/Dec is farther "
             "than this (arcsec). Prevents fitting a noise peak when the catalog "
             "has no detection at the secure host position.",
    )
    parser.add_argument(
        "--host-pad",
        type=int,
        default=20,
        help="Padding (px) on each side of every jointly-fitted galaxy's "
             "segmentation bbox when building / growing the ROI (default 20).",
    )
    parser.add_argument(
        "--re-sep-factor",
        type=float,
        default=3.0,
        help="Mask a clipping neighbor galaxy when the host–neighbor centroid "
             "separation exceeds this factor times the neighbor's FLUX_RADIUS "
             "Re seed (same recipe as the GALFIT re initial guess). Otherwise "
             "expand the ROI to give both a --host-pad boundary and jointly "
             "fit (default 3.0).",
    )
    parser.add_argument(
        "--neighbor-class-star-max",
        type=float,
        default=0.75,
        help="SExtractor CLASS_STAR >= this: source is treated as a point "
             "source and goes to the bad-pixel mask instead of being fitted as "
             "a Sérsic (default 0.75).",
    )
    parser.add_argument(
        "--max-roi-iterations",
        type=int,
        default=8,
        help="Safety cap on the expand-and-recategorise loop (default 8).",
    )
    parser.add_argument(
        "--max-fit-components",
        type=int,
        default=25,
        help="Maximum Sérsic components in the cutout (host always kept). "
             "Excess neighbor galaxies are masked instead of fitted (default 25). "
             "Set 0 to disable the cap.",
    )
    parser.add_argument(
        "--max-cutout-side",
        type=int,
        default=512,
        help="Maximum cutout width or height in pixels; stops ROI expansion when "
             "the bbox would exceed this (default 512). Set 0 to disable.",
    )
    parser.add_argument(
        "--host-only-min-bbox-side",
        type=int,
        default=400,
        help="When the host seg bbox max side exceeds this (px), fit only the host "
             "and mask all neighbors (0 = disable; default 400).",
    )
    parser.add_argument(
        "--host-only-min-elongation",
        type=float,
        default=3.0,
        help="When the host ELONGATION exceeds this, fit only the host and mask "
             "neighbors (default 3.0).",
    )
    parser.add_argument(
        "--mag-aper-index",
        type=int,
        default=-1,
        help="Column index into MAG_APER for the production magnitude (MAG_40PX). "
             "Negative -> use the last (largest) aperture column. master_run sets "
             "this from the resolved aperture ladder.",
    )
    parser.add_argument(
        "--psf-match-arcsec",
        type=float,
        default=0.5,
        help="Sky-match tolerance (arcsec) when keying a Phase-1 image.cat row to "
             "the Phase-2 image.psf.cat SPREAD value (default 0.5).",
    )
    parser.add_argument(
        "--no-data-sigma",
        type=float,
        default=1.0e30,
        help="Sigma value assigned to invvar<=0 / non-finite pixels so GALFIT "
             "ignores them (default 1e30).",
    )
    parser.add_argument(
        "--sigma-rescale-min",
        type=float,
        default=0.5,
        help="Lower bound on k=sky_MAD/sigma_invvar; below this the invvar sigma "
             "map is globally rescaled to the empirical sky noise (default 0.5).",
    )
    parser.add_argument(
        "--sigma-rescale-max",
        type=float,
        default=2.0,
        help="Upper bound on k=sky_MAD/sigma_invvar; above this the invvar sigma "
             "map is globally rescaled to the empirical sky noise (default 2.0).",
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    img_path = os.path.abspath(args.image)
    img_dir = os.path.dirname(img_path)
    
    invvar_path = os.path.join(img_dir, "invvar.fits")
    if not os.path.exists(invvar_path):
         invvar_path = os.path.join(img_dir, os.path.basename(img_path).replace("_flux.fits", "_invvar.fits"))
    has_invvar = os.path.exists(invvar_path)

    # Load outputs
    seg_path = os.path.abspath(args.segmap)
    cat_path = os.path.abspath(args.catalog)
    
    with fits.open(seg_path) as hdul:
        seg_map = hdul[0].data

    cat = get_table_from_ldac(cat_path)

    # Load true image for WCS math
    with fits.open(img_path) as hdul:
        w = WCS(hdul[0].header).celestial
        img_data = np.squeeze(hdul[0].data)
        while img_data.ndim > 2: img_data = img_data[0]

    # Decide which RA/Dec to actually centre the cutout on.  Default rule:
    # prefer AstroPath's chosen host (if Phase 2 wrote a posteriors file
    # alongside the catalog) so GALFIT always fits the same source AstroPath
    # identified; otherwise fall back to the FRB position and let the
    # nearest-source heuristic below pick a target.
    astropath_sex_number = None
    if args.no_astropath_override:
        target_ra, target_dec, target_src = args.ra, args.dec, "localization host (--ra/--dec)"
    else:
        posteriors_path = args.astropath_posteriors or os.path.join(
            os.path.dirname(cat_path), "astropath_posteriors.csv")
        target_ra, target_dec, target_src, astropath_sex_number = _resolve_target_from_astropath(
            posteriors_path, args.ra, args.dec, args.min_astropath_posterior)
    log.info(f"Phase 3a target: {target_src} -> RA={target_ra:.6f}, Dec={target_dec:.6f}")

    spread_by_number = _load_spread_by_number(cat_path)
    target_coord = SkyCoord(ra=target_ra, dec=target_dec, unit="deg", frame="icrs")

    max_astropath_host_sep = min(float(args.max_host_sep_arcsec), 2.0)
    tried_astropath_fallback = False
    while True:
        if astropath_sex_number is not None:
            psf_number = int(astropath_sex_number)
            spread = spread_by_number.get(psf_number)
            if spread is None:
                raise SystemExit(
                    f"[Phase 3a] AstroPath host sex_number={psf_number} missing from image.psf.cat."
                )
            sm, se = spread
            if _spread_is_point_source(sm, se):
                raise SystemExit(
                    f"[Phase 3a] AstroPath host sex_number={psf_number} fails galaxy SPREAD cut "
                    f"(SPREAD={sm:.4f}+/-{se:.4f})."
                )
            target_objid = _seg_number_for_psf_number(psf_number, cat, cat_path)
            sel = cat["NUMBER"] == target_objid
            host_coord = SkyCoord(
                ra=cat["ALPHAWIN_J2000"][sel][0],
                dec=cat["DELTAWIN_J2000"][sel][0],
                unit="deg",
                frame="icrs",
            )
            sep_arcsec = float(target_coord.separation(host_coord).arcsec)
            if sep_arcsec > max_astropath_host_sep:
                log.warning(
                    f"AstroPath sex_number={psf_number} -> seg NUMBER={target_objid} is "
                    f"{sep_arcsec:.2f}\" from the association centre "
                    f"(>{max_astropath_host_sep:.1f}\") — picking nearest galaxy at AstroPath RA/Dec"
                )
                idx, target_objid, sep_arcsec = _pick_nearest_galaxy_host(
                    target_coord, cat, spread_by_number, cat_path,
                    float(args.max_host_sep_arcsec), float(args.psf_match_arcsec),
                )
                log.info(
                    f"Host pick: nearest galaxy at AstroPath centre -> seg NUMBER={target_objid} "
                    f"({sep_arcsec:.2f}\" from {target_src})"
                )
            else:
                log.info(
                    f"Host pick: AstroPath sex_number={psf_number} -> seg NUMBER={target_objid} "
                    f"(SPREAD={sm:.4f}+/-{se:.4f}; {sep_arcsec:.2f}\" from association centre)"
                )
            log.info(f"Assured Target Host NUMBER={target_objid}")
            break

        log.info(
            f"Host pick: nearest galaxy within {args.max_host_sep_arcsec}\" "
            f"(SPREAD+{_SPREAD_SIGMA}*SPREADERR >= {_SPREAD_STAR_MAX}; "
            f"{len(spread_by_number)} sources with SPREAD in image.psf.cat)"
        )
        seg_pick = None
        if args.no_astropath_override:
            seg_pick = _pick_host_at_target_pixel(
                w, target_coord, seg_map, cat, spread_by_number, cat_path,
                float(args.psf_match_arcsec),
            )
            if seg_pick is not None:
                idx, target_objid, sep_arcsec = seg_pick
                log.info(
                    f"Host pick: seg map at target pixel -> image.cat NUMBER={target_objid} "
                    f"({sep_arcsec:.2f}\" from {target_src})"
                )
        if seg_pick is None:
            try:
                idx, target_objid, sep_arcsec = _pick_nearest_galaxy_host(
                    target_coord, cat, spread_by_number, cat_path,
                    float(args.max_host_sep_arcsec), float(args.psf_match_arcsec),
                )
            except SystemExit:
                if (
                    not args.no_astropath_override
                    or tried_astropath_fallback
                ):
                    raise
                tried_astropath_fallback = True
                posteriors_path = args.astropath_posteriors or os.path.join(
                    os.path.dirname(cat_path), "astropath_posteriors.csv")
                ap_ra, ap_dec, ap_src, ap_sex = _resolve_target_from_astropath(
                    posteriors_path, args.ra, args.dec, args.min_astropath_posterior)
                if ap_sex is None:
                    raise
                log.warning(
                    "Localization host (--ra/--dec): no galaxy within "
                    f"{args.max_host_sep_arcsec}\" — falling back to AstroPath host ({ap_src})"
                )
                target_ra, target_dec, target_src = ap_ra, ap_dec, ap_src
                target_coord = SkyCoord(ra=target_ra, dec=target_dec, unit="deg", frame="icrs")
                astropath_sex_number = ap_sex
                continue

        if sep_arcsec > 1.0:
            log.warning(
                f"Host galaxy #{target_objid} is {sep_arcsec:.2f} arcsec "
                f"from the target position."
            )
        log.info(
            f"Assured Target Host NUMBER={target_objid} "
            f"({sep_arcsec:.2f}\" from {target_src})"
        )
        break

    log.info(
        "Neighbor policy (Re-separation): start from host seg bbox + "
        f"{int(args.host_pad)} px pad; for every seg island that clips the ROI:\n"
        f"      star / point-like SPREAD          -> mask (never expand)\n"
        f"      sep > {float(args.re_sep_factor):.2f} × Re_neighbor  -> mask "
        "(Re = FLUX_RADIUS GALFIT seed)\n"
        f"      sep <= {float(args.re_sep_factor):.2f} × Re_neighbor -> jointly fit; "
        f"grow ROI so every fit member has a {int(args.host_pad)} px pad\n"
        "      Host is always GALFIT component 1"
    )

    host_bbox = _bbox_of_objids(seg_map, [int(target_objid)])
    if host_bbox is None:
        log.error("Target not in segmentation map.")
        return
    target_px, target_py = w.world_to_pixel(target_coord)
    target_px = float(np.squeeze(target_px))
    target_py = float(np.squeeze(target_py))

    try:
        xmin, xmax, ymin, ymax, fit_objids, mask_objids, decisions = (
            resolve_neighbor_re_roi(
                seg_map,
                cat,
                int(target_objid),
                host_pad=int(args.host_pad),
                re_sep_factor=float(args.re_sep_factor),
                max_roi_iterations=int(args.max_roi_iterations),
                max_cutout_side=int(args.max_cutout_side),
                neighbor_class_star_max=float(args.neighbor_class_star_max),
                spread_by_number=spread_by_number,
                catalog_path=cat_path,
                psf_match_arcsec=float(args.psf_match_arcsec),
            )
        )
    except ValueError as exc:
        log.error(str(exc))
        return

    # Host is ALWAYS fitted as a Sérsic — even if CLASS_STAR / SPREAD disagree.
    fit_objids.add(int(target_objid))
    mask_objids.discard(int(target_objid))

    kind_labels = {
        "host": "host -> fit",
        "fit": "near (sep<=factor*Re_n) -> fit",
        "mask_far": "far (sep>factor*Re_n) -> mask",
        "mask_star": "stellar -> mask",
        "mask_spread": "point SPREAD -> mask",
        "mask_nocat": "no catalog row -> mask",
    }
    for uid in sorted(decisions.keys()):
        kind, detail = decisions[uid]
        lbl = kind_labels.get(kind, kind)
        if kind in ("fit", "mask_far") and detail is not None and detail >= 0:
            log.info(f"    objid {uid}: sep/Re_n={float(detail):.2f}  {lbl}")
        else:
            log.info(f"    objid {uid}: {lbl}")

    max_fit = int(args.max_fit_components)
    if max_fit > 0 and len(fit_objids) > max_fit:
        host_id = int(target_objid)
        extras = sorted(
            fit_objids - {host_id},
            key=lambda u: _pixel_dist_from_host(cat, host_id, u),
        )
        keep_extras = extras[: max(0, max_fit - 1)]
        drop = set(extras) - set(keep_extras)
        fit_objids = {host_id} | set(keep_extras)
        mask_objids |= drop
        log.warning(
            f"Capped fit components at {max_fit} (host + {len(keep_extras)} neighbors); "
            f"masked {len(drop)} galaxy(ies) that would otherwise be fitted."
        )

    host_sel = cat["NUMBER"] == int(target_objid)
    host_elong = 1.0
    host_bbox_w = host_bbox[1] - host_bbox[0]
    host_bbox_h = host_bbox[3] - host_bbox[2]
    if np.any(host_sel) and "ELONGATION" in cat.colnames:
        host_elong = float(cat[host_sel][0]["ELONGATION"])
    min_bbox_side = int(args.host_only_min_bbox_side)
    min_elong = float(args.host_only_min_elongation)
    extended_host = (
        (min_bbox_side > 0 and max(host_bbox_w, host_bbox_h) >= min_bbox_side)
        or (min_elong > 0 and host_elong >= min_elong)
    )
    host_only_galfit = False
    if extended_host:
        all_touching = _ids_touching_roi(seg_map, xmin, xmax, ymin, ymax, drop=[])
        n_fit_before = len(fit_objids)
        fit_objids = {int(target_objid)}
        mask_objids = (mask_objids | all_touching) - {int(target_objid)}
        host_only_galfit = True
        log.warning(
            f"Extended host (bbox {host_bbox_w}x{host_bbox_h}px, elong={host_elong:.2f}): "
            f"host-only GALFIT — dropped {n_fit_before - 1} neighbor component(s), "
            f"masking {len(mask_objids)} seg object(s)."
        )

    log.info(
        f"Final cutout: {len(fit_objids)} Sérsic component(s), "
        f"{len(mask_objids)} mask-only object(s); "
        f"bounds X=[{xmin}:{xmax}] Y=[{ymin}:{ymax}]"
    )

    # Output the FITS cutout
    cutout_data = img_data[ymin:ymax, xmin:xmax]
    cutout_header = w[ymin:ymax, xmin:xmax].to_header() if w.naxis > 0 else hdul[0].header.copy()
    
    fits.PrimaryHDU(data=cutout_data, header=cutout_header).writeto(os.path.join(args.outdir, "host_cutout.fits"), overwrite=True)

    # 4a. Bad-pixel mask from segmentation (needed early so we can isolate sky pixels
    # when sanity-checking the absolute scale of the inverse-variance map).
    final_sub_map = seg_map[ymin:ymax, xmin:xmax]
    bad_pixel_mask = np.zeros_like(cutout_data, dtype=np.int16)

    # 4b. Output Sigma Cutout. Pixels with invvar <= 0 or non-finite (Legacy / DECaLS
    # convention for "no data" / masked) get a huge sigma so GALFIT effectively
    # ignores them; those same pixels are also marked bad in host_mask.fits below.
    #
    # Some Legacy Surveys frames are delivered with flux and invvar on inconsistent
    # unit conventions (e.g. flux in ~image counts but invvar in 1/nanomaggy^2).
    # In that case sigma = 1/sqrt(invvar) underestimates the real pixel noise by
    # many orders of magnitude and GALFIT's chi^2/nu blows up to ~1e8-1e9 even
    # when the structural fit itself is fine (observed for 20210807D, 20211127I,
    # 20220207C, 20220307B, 20220319D, 20220825A, 20220912A). We anchor the
    # cutout sigma to the empirically measured sky noise (robust MAD of unmasked,
    # un-segmented pixels) by applying a single global multiplicative factor,
    # which preserves the spatial structure of the invvar map (depth / coverage
    # variations) while pinning its absolute scale to the data.
    invvar_cutout = None
    if has_invvar:
        with fits.open(invvar_path) as hdul_inv:
            inv_data = np.squeeze(hdul_inv[0].data)
            while inv_data.ndim > 2: inv_data = inv_data[0]
        no_data_sigma = float(args.no_data_sigma)
        valid_invvar_full = np.isfinite(inv_data) & (inv_data > 0)
        sigma_full = np.full(inv_data.shape, no_data_sigma, dtype=np.float32)
        sigma_full[valid_invvar_full] = (1.0 / np.sqrt(inv_data[valid_invvar_full])).astype(np.float32)
        cutout_sigma = sigma_full[ymin:ymax, xmin:xmax].copy()
        invvar_cutout = inv_data[ymin:ymax, xmin:xmax]

        # Empirical sky noise from cutout (pixels with seg==0 and valid invvar)
        sky_ok = (final_sub_map == 0) & (invvar_cutout > 0) & np.isfinite(invvar_cutout) & np.isfinite(cutout_data)
        sky_pix = cutout_data[sky_ok]
        if sky_pix.size >= 50:
            sky_med = float(np.median(sky_pix))
            sky_mad = float(np.median(np.abs(sky_pix - sky_med)))
            sigma_emp = 1.4826 * sky_mad
            sigma_iv_sky = cutout_sigma[sky_ok]
            sigma_iv_med = float(np.median(sigma_iv_sky)) if sigma_iv_sky.size else 0.0
            if sigma_iv_med > 0 and sigma_emp > 0:
                k = sigma_emp / sigma_iv_med
                if (k < float(args.sigma_rescale_min)) or (k > float(args.sigma_rescale_max)):
                    keep = cutout_sigma < (no_data_sigma * 0.99)  # don't touch flagged "no data" pixels
                    cutout_sigma[keep] = (cutout_sigma[keep] * np.float32(k)).astype(np.float32)
                    log.warning(
                        f"host_sigma scale mismatch: sigma_invvar_med={sigma_iv_med:.3g}, "
                        f"sky_MAD*1.4826={sigma_emp:.3g}, k={k:.3g} -> rescaling sigma by k "
                        f"(preserves spatial structure; only absolute scale changes)."
                    )
                else:
                    log.info(
                        f"host_sigma scale OK (k={k:.3g}; "
                        f"sigma_invvar={sigma_iv_med:.3g}, sky_MAD*1.4826={sigma_emp:.3g})."
                    )
            else:
                log.warning("host_sigma scale check skipped: sigma_invvar or sigma_emp <= 0.")
        else:
            log.warning(f"host_sigma scale check skipped: only {sky_pix.size} sky pixels available.")

        fits.PrimaryHDU(data=cutout_sigma, header=cutout_header).writeto(
            os.path.join(args.outdir, "host_sigma.fits"), overwrite=True)

    # Mask every detection that is not explicitly fitted as a Sérsic (includes
    # mask_objids and any other stray labels inside the stamp).
    mask_condition = (final_sub_map > 0) & (~np.isin(final_sub_map, list(fit_objids)))
    bad_pixel_mask[mask_condition] = 1  # 1 is GALFIT standard bad pixel

    # Also flag any pixel whose invvar is zero / non-finite (no usable data there)
    if invvar_cutout is not None:
        invvar_bad = ~(np.isfinite(invvar_cutout) & (invvar_cutout > 0))
        bad_pixel_mask[invvar_bad] = 1

    fits.PrimaryHDU(data=bad_pixel_mask).writeto(os.path.join(args.outdir, "host_mask.fits"), overwrite=True)

    # 4. Output the Component Manifest for GALFIT
    comp_astropy = cat[np.isin(cat["NUMBER"], list(fit_objids))]
    
    # Export the production-aperture magnitude. Column index comes from
    # --mag-aper-index (master_run forwards the resolved ladder index); a
    # negative or out-of-range value falls back to the last/largest aperture.
    # The output column keeps the name MAG_40PX for downstream compatibility
    # even when the production aperture is not literally 40 px.
    if 'MAG_APER' in comp_astropy.colnames and len(comp_astropy['MAG_APER'].shape) > 1:
        _n_aper = comp_astropy['MAG_APER'].shape[1]
        _ai = args.mag_aper_index if 0 <= args.mag_aper_index < _n_aper else _n_aper - 1
        comp_astropy['MAG_40PX'] = comp_astropy['MAG_APER'][:, _ai]
        
    valid_cols = [col for col in comp_astropy.colnames if len(comp_astropy[col].shape) <= 1]
    comp_df = comp_astropy[valid_cols].to_pandas()
    # GALFIT feedme lists Sersic components in CSV row order (run_galfit_fitting.py).
    # The FRB host must be component 1 so fit.log's first `sersic` line matches
    # the host; neighbor galaxies follow.  Catalog row order is arbitrary.
    _tn = int(target_objid)
    comp_df = pd.concat(
        [comp_df[comp_df["NUMBER"] == _tn], comp_df[comp_df["NUMBER"] != _tn]],
        ignore_index=True,
    )
    # Shift X/Y coordinates to be relative to the cutout window (GALFIT indexing is 1-based, Python 0-based)
    # X_IMAGE is 1-based from SExtractor. xmin is 0-based index. 
    # The new X_IMAGE in cutout will be just original X_IMAGE - xmin
    comp_df['XC_CUTOUT'] = comp_df['X_IMAGE'] - xmin
    comp_df['YC_CUTOUT'] = comp_df['Y_IMAGE'] - ymin
    
    comp_df.to_csv(os.path.join(args.outdir, "host_components.csv"), index=False)
    cutout_meta = {
        "host_number": int(target_objid),
        "extended_host": bool(extended_host),
        "host_only_galfit": bool(host_only_galfit),
        "host_elongation": float(host_elong),
        "host_bbox_px": [int(host_bbox_w), int(host_bbox_h)],
        "cutout_bounds": [int(xmin), int(xmax), int(ymin), int(ymax)],
        "n_fit_components": int(len(fit_objids)),
        "n_mask_objects": int(len(mask_objids)),
        "host_pad": int(args.host_pad),
        "re_sep_factor": float(args.re_sep_factor),
        "neighbor_policy": "re_separation",
    }
    with open(os.path.join(args.outdir, "cutout_meta.json"), "w", encoding="utf-8") as f:
        json.dump(cutout_meta, f, indent=2)
    log.info(
        f"Complete. {len(fit_objids)} fitted + {len(mask_objids)} masked. "
        f"Saved to {args.outdir}"
    )

    # Generate Visual QA Plot
    try:
        import matplotlib.pyplot as plt
        from astropy.visualization import ZScaleInterval
        fig, axes = plt.subplots(1, 2, figsize=(14, 7))
        interval = ZScaleInterval()
        vmin, vmax = interval.get_limits(cutout_data)
        
        # Subplot 1: Clean Cutout with Centers marked
        axes[0].imshow(cutout_data, origin='lower', cmap='bone', vmin=vmin, vmax=vmax)
        axes[0].set_title("Generated GALFIT Cutout (Green=Target, Blue=Neighbors)")
        
        target_df = comp_df[comp_df['NUMBER'] == target_objid]
        if not target_df.empty:
             axes[0].scatter(target_df['XC_CUTOUT'], target_df['YC_CUTOUT'], s=200, facecolors='none', edgecolors='lime', lw=2)
             axes[0].text(target_df['XC_CUTOUT'].iloc[0], target_df['YC_CUTOUT'].iloc[0]-5, "Host", color='lime', ha='center', va='top')
        # CSV / association centre (may differ from SExtractor centroid on extended hosts)
        tx_cut = target_px - xmin
        ty_cut = target_py - ymin
        axes[0].scatter([tx_cut], [ty_cut], s=80, marker='+', c='yellow', lw=2)
        axes[0].text(tx_cut, ty_cut + 8, "Target", color='yellow', ha='center', va='bottom', fontsize=9)

        neighbor_df = comp_df[comp_df['NUMBER'] != target_objid]
        if not neighbor_df.empty:
             axes[0].scatter(neighbor_df['XC_CUTOUT'], neighbor_df['YC_CUTOUT'], s=100, facecolors='none', edgecolors='dodgerblue', lw=2)

        # Subplot 2: Bad Pixel Mask overlay
        axes[1].imshow(cutout_data, origin='lower', cmap='bone', vmin=vmin, vmax=vmax)
        axes[1].set_title("Bad pixel mask (red): point sources + fringe detections")
        
        mask_overlay = np.ma.masked_where(bad_pixel_mask == 0, bad_pixel_mask)
        axes[1].imshow(mask_overlay, origin='lower', cmap='autumn', alpha=0.5)
        
        qa_path = os.path.join(args.outdir, "qa_cutout_mask.png")
        plt.savefig(qa_path, bbox_inches='tight', dpi=150)
        plt.close(fig)
        log.info(f"Visual QA Map saved to {qa_path}")
    except Exception as e:
        log.warning(f"QA plotting failed: {e}")

if __name__ == "__main__":
    main()
