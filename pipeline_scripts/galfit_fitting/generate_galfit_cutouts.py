import os
import argparse

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u

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


def _spread_for_catalog_index(idx: int, cat, spread_by_number: dict, catalog_path: str):
    """SPREAD for a Phase-1 image.cat row, keyed by Phase-2 NUMBER when they differ."""
    num = int(cat["NUMBER"][idx])
    spread = spread_by_number.get(num)
    if spread is not None:
        return spread
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
    if float(seps_arcsec[j]) > 0.5:
        return None
    return spread_by_number.get(int(psf["NUMBER"][j]))


def _seg_number_for_psf_number(psf_number: int, cat, catalog_path: str, max_sep_arcsec: float = 1.0):
    """Map AstroPath / Phase-2 NUMBER to Phase-1 segmentation NUMBER via sky position."""
    if int(psf_number) in set(int(n) for n in cat["NUMBER"]):
        return int(psf_number)
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
        print(
            f"[*] AstroPath sex_number={psf_number} (image.psf.cat) -> "
            f"seg NUMBER={seg_number} (image.cat; {sep:.3f}\")"
        )
    return seg_number


def _pick_nearest_galaxy_host(
    target_coord: SkyCoord,
    cat,
    spread_by_number: dict,
    catalog_path: str,
    max_sep_arcsec: float,
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
        spread = _spread_for_catalog_index(idx, cat, spread_by_number, catalog_path)
        if spread is None:
            print(
                f"    skip #{num}: no SPREAD in image.psf.cat ({sep:.2f}\" from target)"
            )
            continue
        sm, se = spread
        if _spread_is_point_source(sm, se):
            print(
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
        help="Initial padding (px) added on each side of the host segmentation "
             "bbox to build the starting ROI (default 20).",
    )
    parser.add_argument(
        "--contain-thresh",
        type=float,
        default=0.95,
        help="A source whose seg pixels are at least this fraction inside the "
             "current ROI is treated as fully contained and gets fitted as a "
             "Sérsic (default 0.95).",
    )
    parser.add_argument(
        "--expand-thresh",
        type=float,
        default=0.50,
        help="A source whose seg pixels are at least this fraction inside the "
             "current ROI (but below --contain-thresh) is 'largely filled': the "
             "ROI grows to contain it fully and it is then fitted. Below this, "
             "the source is treated as a fringe and only the in-frame pixels "
             "are masked (default 0.50).",
    )
    parser.add_argument(
        "--neighbor-class-star-max",
        type=float,
        default=0.90,
        help="SExtractor CLASS_STAR >= this: source is treated as a point "
             "source and goes to the bad-pixel mask instead of being fitted as "
             "a Sérsic (default 0.90). Affects both contained and largely-filled "
             "categories.",
    )
    parser.add_argument(
        "--max-roi-iterations",
        type=int,
        default=6,
        help="Safety cap on the expand-and-recategorise loop (default 6).",
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
    print(f"[*] Phase 3a target: {target_src} -> RA={target_ra:.6f}, Dec={target_dec:.6f}")

    spread_by_number = _load_spread_by_number(cat_path)
    target_coord = SkyCoord(ra=target_ra, dec=target_dec, unit="deg", frame="icrs")

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
        print(
            f"[*] Host pick: AstroPath sex_number={psf_number} -> seg NUMBER={target_objid} "
            f"(SPREAD={sm:.4f}+/-{se:.4f}; {sep_arcsec:.2f}\" from association centre)"
        )
        print(f"[*] Assured Target Host NUMBER={target_objid}")
    else:
        print(
            f"[*] Host pick: nearest galaxy within {args.max_host_sep_arcsec}\" "
            f"(SPREAD+{_SPREAD_SIGMA}*SPREADERR >= {_SPREAD_STAR_MAX}; "
            f"{len(spread_by_number)} sources with SPREAD in image.psf.cat)"
        )
        idx, target_objid, sep_arcsec = _pick_nearest_galaxy_host(
            target_coord, cat, spread_by_number, cat_path, float(args.max_host_sep_arcsec)
        )
        if sep_arcsec > 1.0:
            print(
                f"WARNING: Host galaxy #{target_objid} is {sep_arcsec:.2f} arcsec "
                f"from the target position."
            )
        print(
            f"[*] Assured Target Host NUMBER={target_objid} "
            f"({sep_arcsec:.2f}\" from {target_src})"
        )

    print(
        "[*] Neighbor policy (per-source containment): for every seg island that "
        "touches the ROI we compute frac = pixels_in_ROI / total_pixels and "
        f"decide:\n"
        f"      frac >= {args.contain_thresh:.2f}  -> fully contained: fit (or mask if stellar)\n"
        f"      frac >= {args.expand_thresh:.2f}  -> largely filled : expand ROI to contain it, then fit/mask\n"
        f"      frac >  0.0                       -> small fringe   : mask in-frame pixels only (no ROI grow)\n"
        f"      no overlap                        -> ignored entirely"
    )

    # --- Starting ROI = host seg bbox + small pad ---
    host_bbox = _bbox_of_objids(seg_map, [int(target_objid)])
    if host_bbox is None:
        print("ERROR: Target not in segmentation map.")
        return
    xmin, xmax, ymin, ymax = _pad_bbox(*host_bbox, pad=int(args.host_pad), shape=seg_map.shape)

    contain_thresh = float(args.contain_thresh)
    expand_thresh = float(args.expand_thresh)
    cstar_lim = float(args.neighbor_class_star_max)

    # State maps so we have a clean record per source.
    decisions: dict = {}   # uid -> ("contained" | "expanded" | "fringe", frac)
    # Host is always fitted by construction; flag as contained 1.0.
    decisions[int(target_objid)] = ("contained", 1.0)

    # Expand-and-recategorise loop. We only grow the ROI; we never shrink it.
    # Termination: when no new island gets promoted from fringe -> expanded.
    max_iters = max(1, int(args.max_roi_iterations))
    for it in range(1, max_iters + 1):
        touching = _ids_touching_roi(seg_map, xmin, xmax, ymin, ymax, drop=[int(target_objid)])
        new_expanded: set = set()
        for uid in sorted(touching):
            total = int(np.sum(seg_map == uid))
            if total <= 0:
                continue
            in_roi = int(np.sum(seg_map[ymin:ymax, xmin:xmax] == uid))
            frac = in_roi / float(total)
            prev = decisions.get(uid)
            if frac >= contain_thresh:
                decisions[uid] = ("contained", frac)
            elif frac >= expand_thresh:
                decisions[uid] = ("expanded", frac)
                if prev is None or prev[0] == "fringe":
                    new_expanded.add(uid)
            else:
                # Only register as fringe if we haven't already classified it
                # as something stronger in a previous iteration.
                if prev is None or prev[0] == "fringe":
                    decisions[uid] = ("fringe", frac)

        if not new_expanded:
            break

        # Grow ROI so every "expanded" source is fully inside, then re-pad with
        # the host pad so we leave room for the model wings.
        expanded_ids = [u for u, (k, _) in decisions.items() if k == "expanded"]
        grow_ids = set(expanded_ids) | {int(target_objid)}
        grow_bbox = _bbox_of_objids(seg_map, list(grow_ids))
        if grow_bbox is None:
            break
        xmin, xmax, ymin, ymax = _pad_bbox(*grow_bbox, pad=int(args.host_pad), shape=seg_map.shape)
        print(f"    [iter {it}] expanded ROI for {sorted(new_expanded)} -> bounds X=[{xmin}:{xmax}] Y=[{ymin}:{ymax}]")

    # Final recompute of fracs in the converged ROI, so the logged number
    # matches what's actually on disk.
    converged_decisions: dict = {int(target_objid): ("contained", 1.0)}
    for uid in _ids_touching_roi(seg_map, xmin, xmax, ymin, ymax, drop=[int(target_objid)]):
        total = int(np.sum(seg_map == uid))
        if total <= 0:
            continue
        in_roi = int(np.sum(seg_map[ymin:ymax, xmin:xmax] == uid))
        frac = in_roi / float(total)
        if frac >= contain_thresh:
            converged_decisions[uid] = ("contained", frac)
        elif frac >= expand_thresh:
            # Should not normally happen post-convergence (we expanded for these),
            # but if the host pad + image edge clipped the bbox we may still see
            # one. Treat the same way: try to fit, else mask depending on type.
            converged_decisions[uid] = ("contained", frac)
        else:
            converged_decisions[uid] = ("fringe", frac)

    # Split into action groups. The host is ALWAYS fitted as a Sérsic — even if
    # SExtractor's CLASS_STAR happens to be >= cstar_lim. The whole purpose of
    # Phase 3a is to model the FRB host, so we never drop it into the mask just
    # because the star/galaxy classifier disagrees.
    fit_objids: set = {int(target_objid)}
    mask_objids: set = set()
    host_frac = converged_decisions.get(int(target_objid), ("contained", 1.0))[1]
    print(f"    objid {int(target_objid)}: frac_in_ROI={host_frac*100:.1f}%  host -> fit")
    for uid, (kind, frac) in converged_decisions.items():
        if int(uid) == int(target_objid):
            continue
        stellar = _is_stellar(cat, int(uid), cstar_lim)
        if kind == "fringe":
            mask_objids.add(int(uid))
            kind_lbl = "fringe -> mask"
        else:
            if stellar:
                mask_objids.add(int(uid))
                kind_lbl = "contained -> mask (stellar)"
            else:
                fit_objids.add(int(uid))
                kind_lbl = "contained -> fit (galaxy)"
        print(f"    objid {uid}: frac_in_ROI={frac*100:.1f}%  {kind_lbl}")

    print(
        f"[*] Final cutout: {len(fit_objids)} Sérsic component(s), "
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
        valid_invvar_full = np.isfinite(inv_data) & (inv_data > 0)
        sigma_full = np.full(inv_data.shape, 1.0e30, dtype=np.float32)
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
                if (k < 0.5) or (k > 2.0):
                    keep = cutout_sigma < 1.0e29  # don't touch flagged "no data" pixels
                    cutout_sigma[keep] = (cutout_sigma[keep] * np.float32(k)).astype(np.float32)
                    print(
                        f"[!] host_sigma scale mismatch: sigma_invvar_med={sigma_iv_med:.3g}, "
                        f"sky_MAD*1.4826={sigma_emp:.3g}, k={k:.3g} -> rescaling sigma by k "
                        f"(preserves spatial structure; only absolute scale changes)."
                    )
                else:
                    print(
                        f"[*] host_sigma scale OK (k={k:.3g}; "
                        f"sigma_invvar={sigma_iv_med:.3g}, sky_MAD*1.4826={sigma_emp:.3g})."
                    )
            else:
                print("[!] host_sigma scale check skipped: sigma_invvar or sigma_emp <= 0.")
        else:
            print(f"[!] host_sigma scale check skipped: only {sky_pix.size} sky pixels available.")

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
    
    # Export 40px magnitude specifically (Index 14 of the standard 15-aperture array)
    if 'MAG_APER' in comp_astropy.colnames and len(comp_astropy['MAG_APER'].shape) > 1:
        comp_astropy['MAG_40PX'] = comp_astropy['MAG_APER'][:, 14]
        
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
    print(
        f"[*] Complete. {len(fit_objids)} fitted + {len(mask_objids)} masked. "
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
        print(f"[*] Visual QA Map saved to {qa_path}")
    except Exception as e:
        print(f"[!] Warning: QA plotting failed: {e}")

if __name__ == "__main__":
    main()
