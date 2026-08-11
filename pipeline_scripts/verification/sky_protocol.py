"""Independent sky determination for GALFIT re-fits (B / E / F).

Uses the multi-band Legacy Survey stamps in ``large_cutouts/{FRB}_flux.fits``
(r-band plane) plus the production segmentation map. Analysis is restricted to
a physically reasonable box around the host — not the full ~10' stamp.

Methods
-------
B  Isophotal / annular flux-growth plateau (GALAPAGOS-style): annular median
   vs radius; sky = median on the outer plateau where the annulus curve is flat.
E  Aggressive mask-growth (Ji-style practical): dilate the source mask and take
   the robust mean of remaining pixels; report the stable dilation range.
F  Empty-patch sampling: medians of many small patches free of segmented
   sources inside the analysis box; optional local plane evaluated at the host.

The consensus sky is the median of methods that return a finite value, with the
MAD of those estimates as the uncertainty. Production ``sky_final_adu`` is
recorded for comparison only — it does not enter the consensus.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from typing import Any

import numpy as np
from astropy.io import fits
from scipy import ndimage

import vercommon as vc

LARGE_CUTOUTS = os.path.join(vc.REPO, "large_cutouts")
REFITS_ROOT = os.path.join(vc.VER_DIR, "Re-fits")

PIXSCALE_ARCSEC = 0.262  # Legacy Survey native; matches CD2_2 ≈ 0.262"/px
DEFAULT_BAND = "r"
BAND_ORDER = "griz"

# Restrict analysis to this half-width around the host (arcsec).
# ~90" is large vs host sizes here but << 10' stamp; avoids distant gradients.
DEFAULT_HALFBOX_ARCSEC = 90.0


# --------------------------------------------------------------------------
# IO
# --------------------------------------------------------------------------

def _band_index(header, band: str = DEFAULT_BAND) -> int:
    bands = str(header.get("BANDS") or "").strip().lower()
    if bands and band.lower() in bands:
        return bands.index(band.lower())
    key = f"BAND{BAND_ORDER.index(band.lower())}" if band.lower() in BAND_ORDER else None
    if key and str(header.get(key, "")).lower() == band.lower():
        return BAND_ORDER.index(band.lower())
    if band.lower() in BAND_ORDER:
        return BAND_ORDER.index(band.lower())
    raise ValueError(f"cannot locate band {band!r} in {bands!r}")


def load_large_rband(frb: str, band: str = DEFAULT_BAND) -> dict[str, Any]:
    flux_path = os.path.join(LARGE_CUTOUTS, f"{frb}_flux.fits")
    inv_path = os.path.join(LARGE_CUTOUTS, f"{frb}_invvar.fits")
    if not os.path.isfile(flux_path):
        raise FileNotFoundError(flux_path)
    with fits.open(flux_path, memmap=False) as hdul:
        data = np.asarray(hdul[0].data, dtype=float)
        hdr = hdul[0].header
    bi = _band_index(hdr, band)
    if data.ndim == 3:
        image = data[bi]
    elif data.ndim == 2:
        image = data
    else:
        raise ValueError(f"unexpected flux shape {data.shape}")

    invvar = None
    if os.path.isfile(inv_path):
        inv = np.asarray(fits.getdata(inv_path), dtype=float)
        invvar = inv[bi] if inv.ndim == 3 else inv

    cd = abs(float(hdr.get("CD2_2") or hdr.get("CDELT2") or 0.0))
    pixscale = cd * 3600.0 if cd > 0 else PIXSCALE_ARCSEC
    return {
        "image": image,
        "invvar": invvar,
        "header": hdr,
        "band": band,
        "band_index": bi,
        "pixscale_arcsec": pixscale,
        "flux_path": flux_path,
        "invvar_path": inv_path if os.path.isfile(inv_path) else None,
    }


def load_host_geometry(frb: str) -> dict[str, Any]:
    hdir = vc.host_dir(frb)
    meta = vc.read_json(os.path.join(hdir, "cutout_meta.json"))
    audit = vc.read_json(os.path.join(hdir, "sky_fit_audit.json"))
    seg_path = os.path.join(hdir, "segmentation_map.fits")
    if not os.path.isfile(seg_path):
        raise FileNotFoundError(seg_path)
    seg = np.asarray(fits.getdata(seg_path), dtype=np.int32)
    host_number = int(meta.get("host_number") or audit.get("host_number") or 0)
    if host_number <= 0:
        raise ValueError(f"{frb}: missing host_number")

    host_mask = seg == host_number
    if not np.any(host_mask):
        raise ValueError(f"{frb}: host {host_number} absent from segmentation map")
    ys, xs = np.nonzero(host_mask)
    yc = float(np.mean(ys))
    xc = float(np.mean(xs))
    # equivalent radius of host footprint (px)
    area = float(np.count_nonzero(host_mask))
    r_host = math.sqrt(area / math.pi)

    bounds = meta.get("cutout_bounds") or []
    return {
        "host_number": host_number,
        "seg": seg,
        "host_mask": host_mask,
        "xc": xc,
        "yc": yc,
        "r_host_px": r_host,
        "cutout_bounds": bounds,
        "sky_final_adu": audit.get("sky_final_adu"),
        "sky_ref_adu": audit.get("sky_ref_adu"),
        "host_dir": hdir,
    }


def analysis_box(
    shape: tuple[int, int],
    xc: float,
    yc: float,
    halfbox_arcsec: float,
    pixscale: float,
) -> tuple[int, int, int, int]:
    """Return y0,y1,x0,x1 (half-open) clipped to the image."""
    half = max(32, int(round(halfbox_arcsec / pixscale)))
    ny, nx = shape
    y0 = max(0, int(math.floor(yc - half)))
    y1 = min(ny, int(math.ceil(yc + half)) + 1)
    x0 = max(0, int(math.floor(xc - half)))
    x1 = min(nx, int(math.ceil(xc + half)) + 1)
    return y0, y1, x0, x1


# --------------------------------------------------------------------------
# robust stats
# --------------------------------------------------------------------------

def _finite(a: np.ndarray, mask: np.ndarray | None = None) -> np.ndarray:
    m = np.isfinite(a)
    if mask is not None:
        m &= mask
    return a[m]


def robust_mean_std(vals: np.ndarray) -> tuple[float, float, int]:
    v = np.asarray(vals, dtype=float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return float("nan"), float("nan"), 0
    med = float(np.median(v))
    mad = float(np.median(np.abs(v - med)))
    sigma = 1.4826 * mad if mad > 0 else float(np.std(v))
    if sigma > 0 and v.size >= 8:
        keep = np.abs(v - med) < 3.0 * sigma
        v = v[keep]
        if v.size:
            med = float(np.median(v))
            mad = float(np.median(np.abs(v - med)))
            sigma = 1.4826 * mad if mad > 0 else float(np.std(v))
    return med, sigma, int(v.size)


# --------------------------------------------------------------------------
# B — annular flux-growth / plateau
# --------------------------------------------------------------------------

def method_B_growth(
    image: np.ndarray,
    source_mask: np.ndarray,
    xc: float,
    yc: float,
    pixscale: float,
    r_max_arcsec: float | None = None,
) -> dict[str, Any]:
    """Annular median vs radius; sky from the outer flat plateau."""
    ny, nx = image.shape
    yy, xx = np.indices(image.shape)
    rr = np.hypot(xx - xc, yy - yc)

    r_host = max(2.0, math.sqrt(float(np.count_nonzero(source_mask & (rr < 50))) / math.pi))
    r_min = max(3.0 * r_host, 8.0 / pixscale)  # start well outside the host
    if r_max_arcsec is None:
        r_max = 0.9 * min(xc, yc, nx - 1 - xc, ny - 1 - yc)
    else:
        r_max = min(r_max_arcsec / pixscale, 0.9 * min(xc, yc, nx - 1 - xc, ny - 1 - yc))
    if r_max <= r_min + 5:
        return {"status": "region_too_small", "sky_adu": float("nan")}

    # exclude all segmented sources from annuli (neighbors contaminate growth)
    clean = ~source_mask & np.isfinite(image)
    n_ann = 24
    edges = np.linspace(r_min, r_max, n_ann + 1)
    radii, ann_med, ann_n = [], [], []
    for a, b in zip(edges[:-1], edges[1:]):
        sel = clean & (rr >= a) & (rr < b)
        vals = image[sel]
        if vals.size < 30:
            continue
        m, _, n = robust_mean_std(vals)
        radii.append(0.5 * (a + b))
        ann_med.append(m)
        ann_n.append(n)

    if len(ann_med) < 4:
        return {"status": "too_few_annuli", "sky_adu": float("nan")}

    radii = np.asarray(radii)
    ann_med = np.asarray(ann_med)
    # plateau = flattest contiguous window of annular medians (min std),
    # preferring outer windows on ties. Walking in from the edge alone fails
    # when large-scale residuals spike the outermost bins.
    n_bins = len(ann_med)
    win = max(4, n_bins // 4)
    best = None  # (std, -start, start, end)
    for start in range(0, n_bins - win + 1):
        end = start + win
        seg = ann_med[start:end]
        if not np.all(np.isfinite(seg)):
            continue
        score = (float(np.std(seg)), -start)
        if best is None or score < best[0:2]:
            best = (score[0], score[1], start, end)
    if best is None:
        start, end = n_bins // 2, n_bins
    else:
        start, end = best[2], best[3]
    plateau = ann_med[start:end]

    sky, sky_sig, n_plat = robust_mean_std(plateau)
    # also report growth-curve slope sky: dF/dA on outer annuli (should ≈ sky)
    # F cumulative on clean pixels is awkward with holes; skip formal F(R).
    return {
        "status": "ok",
        "sky_adu": sky,
        "sky_sigma_adu": sky_sig / max(math.sqrt(n_plat), 1.0),
        "n_plateau_bins": n_plat,
        "r_min_arcsec": float(r_min * pixscale),
        "r_max_arcsec": float(r_max * pixscale),
        "plateau_r_start_arcsec": float(radii[start] * pixscale),
        "annulus_median_adu": [float(x) for x in ann_med],
        "annulus_r_arcsec": [float(x * pixscale) for x in radii],
        "annulus_n": [int(x) for x in ann_n],
    }


# --------------------------------------------------------------------------
# E — aggressive mask growth
# --------------------------------------------------------------------------

def method_E_mask_growth(
    image: np.ndarray,
    source_mask: np.ndarray,
    dilations: tuple[int, ...] = (0, 2, 4, 6, 8, 12, 16, 24, 32, 48),
) -> dict[str, Any]:
    """Sky vs source-mask dilation; take the asymptotic stable end of the curve.

    Undilated medians are biased high by undetected wings; the useful answer is
    where further dilation stops moving the sky.
    """
    curve = []
    for d in dilations:
        if d <= 0:
            m = source_mask
        else:
            m = ndimage.binary_dilation(source_mask, iterations=int(d))
        sky_pix = _finite(image, ~m)
        sky, sig, n = robust_mean_std(sky_pix)
        curve.append({
            "dilation_px": int(d),
            "sky_adu": sky,
            "sky_sigma_adu": sig / max(math.sqrt(n), 1.0) if n else float("nan"),
            "n_pix": n,
            "masked_frac": float(np.count_nonzero(m)) / float(m.size),
        })

    skies = np.asarray([c["sky_adu"] for c in curve], dtype=float)
    ns = np.asarray([c["n_pix"] for c in curve], dtype=float)
    fracs = np.asarray([c["masked_frac"] for c in curve], dtype=float)
    ok = np.isfinite(skies) & (ns >= 500) & (fracs < 0.92)
    if np.count_nonzero(ok) < 2:
        return {"status": "unstable", "sky_adu": float("nan"), "curve": curve}

    idx = np.where(ok)[0]
    # asymptotic: among accepted points, pick the run of smallest |Δsky|
    # starting from the largest dilation and walking inward
    end = int(idx[-1])
    start = end
    diffs = np.abs(np.diff(skies))
    # typical step among high-dilatation pairs
    hi = diffs[max(0, end - 3):end] if end > 0 else diffs
    thr = max(float(np.median(hi)) * 2.5, 1e-9) if hi.size else 1e-9
    for i in range(end - 1, int(idx[0]) - 1, -1):
        if not ok[i]:
            break
        if abs(skies[i + 1] - skies[i]) <= thr:
            start = i
        else:
            break
    used = list(range(start, end + 1))
    sky, sig, n = robust_mean_std(skies[used])
    span = float(np.nanmax(skies[ok]) - np.nanmin(skies[ok]))
    return {
        "status": "ok",
        "sky_adu": sky,
        "sky_sigma_adu": sig / max(math.sqrt(n), 1.0),
        "stability_span_adu": span,
        "dilation_px_used": [int(dilations[i]) for i in used],
        "curve": curve,
    }


# --------------------------------------------------------------------------
# F — empty patches (+ optional local plane)
# --------------------------------------------------------------------------

def method_F_empty_patches(
    image: np.ndarray,
    source_mask: np.ndarray,
    xc: float,
    yc: float,
    patch: int = 16,
    n_want: int = 120,
    seed: int = 0,
) -> dict[str, Any]:
    """Median of medians from empty patches; plane sky at host position."""
    rng = np.random.default_rng(seed)
    ny, nx = image.shape
    # keep patches clear of segmented sources (Ji-style aggressive buffer)
    keepout = ndimage.binary_dilation(source_mask, iterations=max(6, patch // 2))
    free = ~keepout & np.isfinite(image)

    half = patch // 2
    meds = []
    tries = 0
    max_tries = n_want * 60
    while len(meds) < n_want and tries < max_tries:
        tries += 1
        cy = int(rng.integers(half, ny - half))
        cx = int(rng.integers(half, nx - half))
        sl = (slice(cy - half, cy - half + patch), slice(cx - half, cx - half + patch))
        if not np.all(free[sl]):
            continue
        meds.append(float(np.median(image[sl])))

    if len(meds) < 10:
        return {"status": "too_few_patches", "sky_adu": float("nan"), "n_patches": len(meds)}

    meds_arr = np.asarray(meds, dtype=float)
    # plain median-of-medians (no 3σ clip — clip pulls sky low on near-zero fields)
    sky = float(np.median(meds_arr))
    mad = float(np.median(np.abs(meds_arr - sky)))
    sig = 1.4826 * mad
    n = int(meds_arr.size)

    ys, xs = np.nonzero(free)
    if ys.size > 20000:
        pick = rng.choice(ys.size, size=20000, replace=False)
        ys, xs = ys[pick], xs[pick]
    z = image[ys, xs]
    A = np.column_stack([np.ones(ys.size), xs.astype(float), ys.astype(float)])
    try:
        coef, *_ = np.linalg.lstsq(A, z, rcond=None)
        plane_at_host = float(coef[0] + coef[1] * xc + coef[2] * yc)
        grad = math.hypot(float(coef[1]), float(coef[2]))
        plane_ok = grad * min(nx, ny) < 20.0 * (sig if math.isfinite(sig) and sig > 0 else abs(sky) + 1e-9)
    except Exception:
        plane_at_host = float("nan")
        plane_ok = False
        coef = (float("nan"), float("nan"), float("nan"))

    return {
        "status": "ok",
        "sky_adu": sky,
        "sky_sigma_adu": sig / max(math.sqrt(n), 1.0),
        "n_patches": n,
        "patch_px": patch,
        "plane_sky_adu": plane_at_host if plane_ok else float("nan"),
        "plane_coefs": {
            "a0": float(coef[0]),
            "ax": float(coef[1]),
            "ay": float(coef[2]),
        },
        "plane_used": bool(plane_ok),
    }


# --------------------------------------------------------------------------
# consensus + driver
# --------------------------------------------------------------------------

def _consensus(estimates: dict[str, float]) -> dict[str, Any]:
    vals = {k: float(v) for k, v in estimates.items() if math.isfinite(float(v))}
    if not vals:
        return {"sky_adu": float("nan"), "sky_sigma_adu": float("nan"), "n_methods": 0,
                "methods_used": []}
    arr = np.asarray(list(vals.values()), dtype=float)
    med = float(np.median(arr))
    mad = float(np.median(np.abs(arr - med)))
    sig = 1.4826 * mad if mad > 0 else float(np.std(arr)) if arr.size > 1 else 0.0
    return {
        "sky_adu": med,
        "sky_sigma_adu": sig,
        "n_methods": int(arr.size),
        "methods_used": sorted(vals.keys()),
        "span_adu": float(np.max(arr) - np.min(arr)),
    }


def determine_sky(
    frb: str,
    halfbox_arcsec: float = DEFAULT_HALFBOX_ARCSEC,
    band: str = DEFAULT_BAND,
) -> dict[str, Any]:
    cut = load_large_rband(frb, band=band)
    geo = load_host_geometry(frb)
    image_full = cut["image"]
    seg_full = geo["seg"]
    if seg_full.shape != image_full.shape:
        raise ValueError(
            f"{frb}: seg {seg_full.shape} != image {image_full.shape}"
        )

    y0, y1, x0, x1 = analysis_box(
        image_full.shape, geo["xc"], geo["yc"], halfbox_arcsec, cut["pixscale_arcsec"]
    )
    image = image_full[y0:y1, x0:x1]
    seg = seg_full[y0:y1, x0:x1]
    source_mask = seg > 0
    xc = geo["xc"] - x0
    yc = geo["yc"] - y0

    B = method_B_growth(
        image, source_mask, xc, yc, cut["pixscale_arcsec"],
        r_max_arcsec=0.95 * halfbox_arcsec,
    )
    E = method_E_mask_growth(image, source_mask)
    F = method_F_empty_patches(image, source_mask, xc, yc)

    estimates = {
        "B_growth": B.get("sky_adu"),
        "E_mask_growth": E.get("sky_adu"),
        "F_empty_patches": F.get("sky_adu"),
    }
    # include plane only when it agrees with empty-patch sky (else residual bias)
    if F.get("plane_used") and math.isfinite(float(F.get("plane_sky_adu", float("nan")))):
        fp = float(F["plane_sky_adu"])
        f0 = float(F.get("sky_adu", float("nan")))
        fsig = float(F.get("sky_sigma_adu", float("nan")))
        if math.isfinite(f0) and math.isfinite(fsig) and abs(fp - f0) <= max(3.0 * fsig, 1e-9):
            estimates["F_plane"] = fp

    cons = _consensus(estimates)
    # flag when methods disagree by more than ~the catalog sky MAD scale
    if cons["n_methods"] >= 2 and math.isfinite(cons.get("span_adu", float("nan"))):
        cons["agree"] = bool(cons["span_adu"] <= 3.0 * max(cons["sky_sigma_adu"], 1e-12))
    else:
        cons["agree"] = cons["n_methods"] >= 1
    prod = geo.get("sky_final_adu")
    out = {
        "frb": frb,
        "band": cut["band"],
        "band_index": cut["band_index"],
        "flux_path": cut["flux_path"],
        "pixscale_arcsec": cut["pixscale_arcsec"],
        "halfbox_arcsec": halfbox_arcsec,
        "analysis_box_yx": [y0, y1, x0, x1],
        "analysis_box_arcsec": [
            (y1 - y0) * cut["pixscale_arcsec"],
            (x1 - x0) * cut["pixscale_arcsec"],
        ],
        "host_number": geo["host_number"],
        "host_xy_full": [geo["xc"], geo["yc"]],
        "production_sky_final_adu": prod,
        "production_sky_ref_adu": geo.get("sky_ref_adu"),
        "methods": {"B": B, "E": E, "F": F},
        "estimates_adu": estimates,
        "consensus": cons,
        "delta_vs_production_adu": (
            cons["sky_adu"] - float(prod)
            if prod is not None and math.isfinite(cons["sky_adu"]) and math.isfinite(float(prod))
            else None
        ),
    }
    return out


def write_diagnostic_plot(result: dict[str, Any], path: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    B = result["methods"]["B"]
    E = result["methods"]["E"]
    cons = result["consensus"]["sky_adu"]
    prod = result.get("production_sky_final_adu")

    fig, axes = plt.subplots(1, 2, figsize=(9.5, 3.8))

    ax = axes[0]
    if B.get("status") == "ok":
        ax.plot(B["annulus_r_arcsec"], B["annulus_median_adu"], "o-", ms=3, label="annulus median")
        if B.get("plateau_r_start_arcsec") is not None:
            ax.axvline(B["plateau_r_start_arcsec"], color="0.5", ls=":", label="plateau start")
    if math.isfinite(cons):
        ax.axhline(cons, color="C2", ls="--", label=f"consensus {cons:.3e}")
    if prod is not None and math.isfinite(float(prod)):
        ax.axhline(float(prod), color="C3", ls="--", label=f"production {float(prod):.3e}")
    ax.set_xlabel("radius (arcsec)")
    ax.set_ylabel("annular median (ADU)")
    ax.set_title(f"{result['frb']} — B growth")
    ax.legend(fontsize=8, loc="best")

    ax = axes[1]
    if E.get("curve"):
        xs = [c["dilation_px"] for c in E["curve"]]
        ys = [c["sky_adu"] for c in E["curve"]]
        ax.plot(xs, ys, "o-", ms=3, label="E mask-growth")
    est = result.get("estimates_adu") or {}
    if est.get("F_empty_patches") is not None and math.isfinite(float(est["F_empty_patches"])):
        ax.axhline(float(est["F_empty_patches"]), color="C0", ls=":", label="F patches")
    if est.get("F_plane") is not None and math.isfinite(float(est["F_plane"])):
        ax.axhline(float(est["F_plane"]), color="C1", ls=":", label="F plane")
    if math.isfinite(cons):
        ax.axhline(cons, color="C2", ls="--", label="consensus")
    if prod is not None and math.isfinite(float(prod)):
        ax.axhline(float(prod), color="C3", ls="--", label="production")
    ax.set_xlabel("source-mask dilation (px)")
    ax.set_ylabel("sky (ADU)")
    ax.set_title(f"{result['frb']} — E / F")
    ax.legend(fontsize=8, loc="best")

    fig.suptitle(
        f"sky protocol  band={result['band']}  "
        f"box≈{result['halfbox_arcsec']:.0f}\"  "
        f"consensus={cons:.6g}" if math.isfinite(cons) else f"sky protocol {result['frb']}",
        fontsize=10,
    )
    fig.tight_layout()
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    fig.savefig(path, dpi=120)
    plt.close(fig)


def run_one(
    frb: str,
    halfbox_arcsec: float = DEFAULT_HALFBOX_ARCSEC,
    stability_boxes: tuple[float, ...] = (60.0, 90.0, 120.0),
) -> dict[str, Any]:
    result = determine_sky(frb, halfbox_arcsec=halfbox_arcsec)
    # multi-box stability (same methods; sky should not jump across 1–2')
    stab = []
    for hb in stability_boxes:
        if abs(hb - halfbox_arcsec) < 1e-6:
            stab.append({
                "halfbox_arcsec": hb,
                "consensus_sky_adu": result["consensus"]["sky_adu"],
                "span_adu": result["consensus"].get("span_adu"),
                "agree": result["consensus"].get("agree"),
                "estimates_adu": result["estimates_adu"],
            })
            continue
        alt = determine_sky(frb, halfbox_arcsec=hb)
        stab.append({
            "halfbox_arcsec": hb,
            "consensus_sky_adu": alt["consensus"]["sky_adu"],
            "span_adu": alt["consensus"].get("span_adu"),
            "agree": alt["consensus"].get("agree"),
            "estimates_adu": alt["estimates_adu"],
        })
    skies = [s["consensus_sky_adu"] for s in stab
             if s["consensus_sky_adu"] is not None
             and math.isfinite(float(s["consensus_sky_adu"]))]
    result["box_stability"] = {
        "boxes": stab,
        "sky_median_adu": float(np.median(skies)) if skies else float("nan"),
        "sky_span_adu": float(np.max(skies) - np.min(skies)) if len(skies) >= 2 else 0.0,
    }

    out_dir = os.path.join(REFITS_ROOT, frb)
    os.makedirs(out_dir, exist_ok=True)
    json_path = os.path.join(out_dir, "sky_protocol.json")
    plot_path = os.path.join(out_dir, "sky_protocol.png")
    vc.write_json(json_path, result)
    try:
        write_diagnostic_plot(result, plot_path)
        result["plot_path"] = plot_path
    except Exception as exc:
        result["plot_error"] = f"{type(exc).__name__}: {exc}"
    result["json_path"] = json_path
    return result


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("frb", help="FRB name, e.g. 20181112A")
    p.add_argument("--halfbox-arcsec", type=float, default=DEFAULT_HALFBOX_ARCSEC,
                   help="half-width of analysis box around host (default 90\")")
    args = p.parse_args(argv)
    result = run_one(args.frb, halfbox_arcsec=args.halfbox_arcsec)
    cons = result["consensus"]
    print(f"{args.frb}  band={result['band']}  box={args.halfbox_arcsec:.0f}\"")
    for k, v in result["estimates_adu"].items():
        print(f"  {k:16s}  {v:.6e}" if v is not None and math.isfinite(float(v)) else f"  {k:16s}  —")
    print(f"  {'consensus':16s}  {cons['sky_adu']:.6e}  "
          f"(sig={cons['sky_sigma_adu']:.2e}, span={cons.get('span_adu')}, "
          f"agree={cons.get('agree')})")
    prod = result.get("production_sky_final_adu")
    if prod is not None:
        dlt = result.get("delta_vs_production_adu")
        print(f"  {'production':16s}  {float(prod):.6e}  (delta={dlt})")
    stab = result.get("box_stability") or {}
    if stab.get("boxes"):
        print(f"  box_stability     median={stab.get('sky_median_adu'):.6e}  "
              f"span={stab.get('sky_span_adu')}")
    print(f"  wrote {result.get('json_path')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
