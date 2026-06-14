"""Run AstroPhot Sersic fits on all FRB hosts with or without PSF convolution."""

from __future__ import annotations

import argparse
import os
from typing import Any

import numpy as np
import pandas as pd
import astrophot as ap
from astropy.io import fits


def get_pixelscale(header: Any) -> float:
    try:
        return abs(header["CD2_2"]) * 3600
    except KeyError:
        try:
            return header["CDELT2"] * 3600
        except KeyError:
            return 0.258


def load_flux_and_variance(frb_name: str, root: str | None = None) -> tuple[np.ndarray, np.ndarray, float]:
    flux_path = os.path.join("pipeline_scripts", "Output", f"{frb_name}_all", "host_cutout.fits")
    sigma_path = os.path.join("pipeline_scripts", "Output", f"{frb_name}_all", "host_sigma.fits")
    if root:
        flux_path = os.path.join(root, flux_path)
        sigma_path = os.path.join(root, sigma_path)

    with fits.open(flux_path) as hdu:
        flux_data = np.array(hdu[0].data, dtype=np.float64)
        pixelscale = get_pixelscale(hdu[0].header)

    with fits.open(sigma_path) as hdu:
        sigma_data = np.array(hdu[0].data, dtype=np.float64)

    variance_data = sigma_data ** 2
    bad_mask = ~np.isfinite(variance_data) | (variance_data <= 0)
    if np.any(bad_mask):
        good_var = variance_data[np.isfinite(variance_data) & (variance_data > 0)]
        fallback = np.nanmax(good_var) * 10 if len(good_var) > 0 else 1.0
        variance_data[bad_mask] = fallback

    return flux_data, variance_data, pixelscale


def load_psf(frb_name: str, pixelscale: float, psf_dir: str) -> ap.image.PSF_Image:
    psf_path = os.path.join(psf_dir, f"{frb_name}_1x_psf.fits")
    with fits.open(psf_path) as hdu:
        psf_data = np.array(hdu[0].data, dtype=np.float64)
    return ap.image.PSF_Image(data=psf_data, pixelscale=pixelscale)


def _pipeline_frb_names() -> list[str]:
    out_root = os.path.join("pipeline_scripts", "Output")
    if not os.path.isdir(out_root):
        return []
    names = []
    for entry in sorted(os.listdir(out_root)):
        if not entry.endswith("_all"):
            continue
        cutout = os.path.join(out_root, entry, "host_cutout.fits")
        if os.path.isfile(cutout):
            names.append(entry.replace("_all", ""))
    return names


def run_fit(use_psf: bool, output_csv: str, fixed_n: float | None = None) -> pd.DataFrame:
    psf_dir = "psfs/downsampled_psfs"
    frb_names = _pipeline_frb_names()

    print(f"Found {len(frb_names)} galaxies")
    print(f"Fit mode: {'with PSF' if use_psf else 'without PSF'}")
    if fixed_n is not None:
        print(f"Sersic index mode: fixed n = {fixed_n}")
    else:
        print("Sersic index mode: free n")

    results = []

    for frb_name in frb_names:
        print(f"\n{'=' * 60}")
        print(f"Processing {frb_name}")

        flux_data, variance_data, pixelscale = load_flux_and_variance(frb_name)

        target = ap.image.Target_Image(
            data=flux_data,
            pixelscale=pixelscale,
            variance=variance_data,
        )

        model_kwargs = {
            "name": "model with target",
            "model_type": "sersic galaxy model",
            "target": target,
        }
        if use_psf:
            model_kwargs["psf"] = load_psf(frb_name, pixelscale, psf_dir)

        model = ap.models.AstroPhot_Model(**model_kwargs)
        # AstroPhot defaults to psf_mode='none'; set to 'full' to apply convolution.
        model.psf_mode = "full" if use_psf else "none"
        model.initialize()

        if fixed_n is not None:
            model["n"].value = float(fixed_n)
            model["n"].locked = True

        fitter = ap.fit.LM(model, verbose=1)
        result = fitter.fit()

        q = model["q"].value.item()
        inclination_angle = np.degrees(np.arccos(np.clip(q, -1, 1)))
        error_q = model["q"].uncertainty.item()
        # Keep the same uncertainty expression used in the notebook script.
        error_i = 1 / np.sqrt(1 - q**2 * error_q)
        n_sersic = model["n"].value.item()
        pa = model["PA"].value.item()
        re = model["Re"].value.item()

        chi2 = fitter.res_loss()
        ndf = fitter.ndf
        chi2_nu = chi2 / ndf if ndf else np.nan

        print(
            f"  q={q:.4f}, i={inclination_angle:.2f} deg, n={n_sersic:.4f}, "
            f"Re={re:.4f}, chi2={chi2:.4f}, ndf={ndf}, chi2_nu={chi2_nu:.4f}"
        )

        results.append(
            {
                "frb_name": frb_name,
                "fit_mode": "with_psf" if use_psf else "no_psf",
                "n_fixed": float(fixed_n) if fixed_n is not None else np.nan,
                "q": q,
                "inclination_angle": inclination_angle,
                "error_q": error_q,
                "error_i": error_i,
                "n": n_sersic,
                "PA": pa,
                "Re": re,
                "chi2": chi2,
                "ndf": ndf,
                "chi2_nu": chi2_nu,
                "fit_message": result.message,
            }
        )

    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    print(f"\nSaved {len(df)} rows to {output_csv}")
    return df


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["with_psf", "no_psf"], required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--fixed-n", type=float, default=None)
    args = parser.parse_args()

    run_fit(use_psf=args.mode == "with_psf", output_csv=args.output, fixed_n=args.fixed_n)


if __name__ == "__main__":
    main()