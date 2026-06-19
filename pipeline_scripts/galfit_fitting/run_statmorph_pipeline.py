"""Phase Statmorph — non-parametric morphology for the FRB host cutout.

Runs statmorph on the Phase 3a cutout (host_cutout.fits + host_sigma.fits) and
writes ``statmorph_results.json`` to the output directory.  This phase sits
between Phase 3a (cutout generation) and Phase 3b (GALFIT fitting) so its
metrics (CAS, Gini-M₂₀) can eventually be used to inform GALFIT initial
conditions or flag problematic fits.

Usage (standalone):
    python run_statmorph_pipeline.py \\
        --cutout host_cutout.fits --sigma host_sigma.fits \\
        --psf proto_image.fits   --outdir .

When invoked by master_run.py the arguments are wired automatically.
"""

import argparse
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pipeline_shared import get_logger  # noqa: E402

log = get_logger("statmorph")


def _load_fits_data(path: str) -> np.ndarray:
    from astropy.io import fits
    with fits.open(path) as hdul:
        data = np.array(hdul[0].data, dtype=np.float64)
        while data.ndim > 2:
            data = data[0]
    return data


def _clamp_weightmap(wmap: np.ndarray) -> np.ndarray:
    finite = np.isfinite(wmap) & (wmap > 0)
    if not np.any(finite):
        return np.ones_like(wmap)
    floor = max(float(np.nanpercentile(wmap[finite], 5)), 1e-4)
    clean = np.array(wmap, dtype=np.float64)
    clean[~finite] = floor
    return clean


def _make_segmap(image_sub: np.ndarray, std_bg: float) -> np.ndarray:
    """Build a segmentation map centred on the brightest / nearest-centre source."""
    from astropy.convolution import Gaussian2DKernel, convolve
    from photutils.segmentation import detect_sources

    kernel = Gaussian2DKernel(x_stddev=1.5)
    convolved = convolve(image_sub, kernel)
    ny, nx = image_sub.shape
    npixels = max(5, int(0.01 * nx * ny))
    segmap = detect_sources(convolved, 1.5 * std_bg, npixels=npixels)
    if segmap is None:
        segmap = detect_sources(convolved, std_bg, npixels=max(3, npixels // 2))
    if segmap is None:
        raise RuntimeError("no sources detected in segmentation")

    seg_array = segmap.data
    cy, cx = ny // 2, nx // 2
    center_label = seg_array[cy, cx]
    if center_label == 0:
        labels = np.unique(seg_array)
        labels = labels[labels > 0]
        if len(labels) == 0:
            raise RuntimeError("no labeled sources")
        best_label, best_dist = labels[0], math.inf
        for lab in labels:
            ys, xs = np.where(seg_array == lab)
            dist = math.hypot(float(np.mean(ys)) - cy, float(np.mean(xs)) - cx)
            if dist < best_dist:
                best_dist, best_label = dist, lab
        center_label = best_label
    return np.where(seg_array == center_label, center_label, 0).astype(np.int32)


def run_statmorph_on_cutout(
    cutout_path: str,
    sigma_path: str | None = None,
    psf_path: str | None = None,
) -> dict:
    """Run statmorph and return a dict of non-parametric morphology metrics."""
    import statmorph
    from astropy.stats import sigma_clipped_stats

    image = _load_fits_data(cutout_path)
    _, median_bg, std_bg = sigma_clipped_stats(image, sigma=3.0)
    image_sub = image - median_bg

    segmap = _make_segmap(image_sub, std_bg)

    weightmap = None
    if sigma_path and os.path.isfile(sigma_path):
        weightmap = _clamp_weightmap(_load_fits_data(sigma_path))

    psf = None
    if psf_path and os.path.isfile(psf_path):
        psf = _load_fits_data(psf_path)

    kwargs: dict = {"weightmap": weightmap} if weightmap is not None else {}
    if psf is not None:
        kwargs["psf"] = psf

    source_morphs = statmorph.source_morphology(image_sub, segmap, **kwargs)
    morph = source_morphs[0]

    def _safe(val):
        try:
            v = float(val)
            return round(v, 6) if math.isfinite(v) else None
        except (TypeError, ValueError):
            return None

    return {
        "gini": _safe(morph.gini),
        "m20": _safe(morph.m20),
        "concentration": _safe(morph.concentration),
        "asymmetry": _safe(morph.asymmetry),
        "smoothness": _safe(morph.smoothness),
        "rpetro_circ_px": _safe(morph.rpetro_circ),
        "r20_px": _safe(morph.r20),
        "r80_px": _safe(morph.r80),
        "sn_per_pixel": _safe(morph.sn_per_pixel),
        "gini_m20_merger": _safe(morph.gini_m20_merger),
        "gini_m20_bulge": _safe(morph.gini_m20_bulge),
        "flag": int(morph.flag),
        "flag_sersic": int(morph.flag_sersic),
    }


def main():
    parser = argparse.ArgumentParser(
        description="Run statmorph non-parametric morphology on a host cutout."
    )
    parser.add_argument("--cutout", required=True, help="Path to host_cutout.fits")
    parser.add_argument("--sigma", default=None, help="Path to host_sigma.fits (optional)")
    parser.add_argument("--psf", default=None, help="Path to PSF image (e.g. proto_image.fits)")
    parser.add_argument("--outdir", default=".", help="Output directory for statmorph_results.json")
    args = parser.parse_args()

    if not os.path.isfile(args.cutout):
        log.error(f"Cutout not found: {args.cutout}")
        sys.exit(1)

    try:
        results = run_statmorph_on_cutout(args.cutout, args.sigma, args.psf)
    except Exception as exc:
        log.error(f"Statmorph failed: {exc}")
        sys.exit(1)

    out_path = os.path.join(args.outdir, "statmorph_results.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)
    log.info(
        f"Statmorph OK: Gini={results.get('gini')}, M20={results.get('m20')}, "
        f"C={results.get('concentration')}, A={results.get('asymmetry')}, "
        f"S={results.get('smoothness')}"
    )
    log.info(f"Results written to {out_path}")


if __name__ == "__main__":
    main()
