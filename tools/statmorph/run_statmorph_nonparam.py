"""Run statmorph non-parametric morphology on the FRB host sample."""

from __future__ import annotations

import os
import math

import numpy as np
import pandas as pd
import statmorph
from astropy.convolution import Gaussian2DKernel, convolve
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.segmentation import detect_sources


def load_image(path: str) -> np.ndarray:
    with fits.open(path) as hdu:
        return np.array(hdu[0].data, dtype=np.float64)


def clamp_weightmap(weightmap: np.ndarray) -> np.ndarray:
    finite = np.isfinite(weightmap) & (weightmap > 0)
    if not np.any(finite):
        raise ValueError("no valid sigma values for weightmap")
    floor = np.nanpercentile(weightmap[finite], 5)
    floor = max(floor, 1e-4)
    clean = np.array(weightmap, dtype=np.float64)
    bad = ~finite
    clean[bad] = floor
    return clean


def make_segmap(image_sub: np.ndarray, std_bg: float) -> tuple[np.ndarray, int]:
    kernel = Gaussian2DKernel(x_stddev=1.5)
    convolved = convolve(image_sub, kernel)
    ny, nx = image_sub.shape
    npixels = max(5, int(0.01 * nx * ny))
    threshold = 1.5 * std_bg
    segmap = detect_sources(convolved, threshold, npixels=npixels)
    if segmap is None:
        threshold = std_bg
        segmap = detect_sources(convolved, threshold, npixels=max(3, npixels // 2))
    if segmap is None:
        raise RuntimeError("no sources detected in segmentation")

    segmap_array = segmap.data
    cy, cx = ny // 2, nx // 2
    center_label = segmap_array[cy, cx]
    if center_label == 0:
        labels = np.unique(segmap_array)
        labels = labels[labels > 0]
        if len(labels) == 0:
            raise RuntimeError("no labeled sources in segmentation")
        best_label = labels[0]
        best_dist = math.inf
        for lab in labels:
            ys, xs = np.where(segmap_array == lab)
            dist = math.hypot(np.mean(ys) - cy, np.mean(xs) - cx)
            if dist < best_dist:
                best_dist = dist
                best_label = lab
        center_label = best_label
    segmap_central = np.where(segmap_array == center_label, center_label, 0).astype(np.int32)
    return segmap_central, np.sum(segmap_central > 0)


def run() -> pd.DataFrame:
    psf_dir = "psfs/downsampled_psfs"
    out_root = os.path.join("pipeline_scripts", "Output")
    frb_names = []
    if os.path.isdir(out_root):
        for entry in sorted(os.listdir(out_root)):
            if not entry.endswith("_all"):
                continue
            cutout = os.path.join(out_root, entry, "host_cutout.fits")
            if os.path.isfile(cutout):
                frb_names.append(entry.replace("_all", ""))
    print(f"Found {len(frb_names)} galaxies")

    results = []
    for frb_name in frb_names:
        print("\n", "=" * 60)
        print(f"Processing {frb_name}")
        flux_path = os.path.join(out_root, f"{frb_name}_all", "host_cutout.fits")
        sigma_path = os.path.join(out_root, f"{frb_name}_all", "host_sigma.fits")
        psf_path = os.path.join(psf_dir, f"{frb_name}_1x_psf.fits")

        image = load_image(flux_path)
        sigma_image = load_image(sigma_path)
        psf_image = load_image(psf_path)

        mean_bg, median_bg, std_bg = sigma_clipped_stats(image, sigma=3.0)
        image_sub = image - median_bg
        segmap_central, seg_pixels = make_segmap(image_sub, std_bg)
        weightmap = clamp_weightmap(sigma_image)

        try:
            source_morphs = statmorph.source_morphology(
                image_sub,
                segmap_central,
                weightmap=weightmap,
                psf=psf_image,
            )
        except Exception as exc:
            print(f"  statmorph error: {exc}")
            continue

        morph = source_morphs[0]
        print(f"  Segmap center label pixels: {seg_pixels}")
        print(f"  Gini={morph.gini:.4f}, M20={morph.m20:.4f}, C={morph.concentration:.4f}")

        results.append(
            {
                "frb_name": frb_name,
                "gini": morph.gini,
                "m20": morph.m20,
                "concentration": morph.concentration,
                "asymmetry": morph.asymmetry,
                "smoothness": morph.smoothness,
                "rpetro": morph.rpetro_circ,
                "r20": morph.r20,
                "r80": morph.r80,
                "sn_per_pixel": morph.sn_per_pixel,
                "gini_m20_merger": morph.gini_m20_merger,
                "gini_m20_bulge": morph.gini_m20_bulge,
            }
        )

    df = pd.DataFrame(results)
    output_path = "statmorph_nonparam_results.csv"
    df.to_csv(output_path, index=False)
    print(f"\nSaved {len(df)} rows to {output_path}")
    return df


if __name__ == "__main__":
    run()