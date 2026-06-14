import os
import numpy as np
import pandas as pd
import astrophot as ap
from astropy.io import fits


def get_pixelscale(header):
    try:
        return abs(header["CD2_2"]) * 3600
    except KeyError:
        try:
            return header["CDELT2"] * 3600
        except KeyError:
            return 0.258


def load_flux_sigma_pixelscale(frb_name):
    flux_path = os.path.join("pipeline_scripts", "Output", f"{frb_name}_all", "host_cutout.fits")
    sigma_path = os.path.join("pipeline_scripts", "Output", f"{frb_name}_all", "host_sigma.fits")

    with fits.open(flux_path) as hdu:
        flux_data = np.array(hdu[0].data, dtype=np.float64)
        flux_header = hdu[0].header.copy()
    with fits.open(sigma_path) as hdu:
        sigma_data = np.array(hdu[0].data, dtype=np.float64)
        sigma_header = hdu[0].header.copy()

    pixelscale = get_pixelscale(flux_header)
    return flux_data, sigma_data, pixelscale, flux_header, sigma_header


def load_psf_image(frb_name, pixelscale):
    psf_path = os.path.join("psfs", "downsampled_psfs", f"{frb_name}_1x_psf.fits")
    with fits.open(psf_path) as hdu:
        psf_data = np.array(hdu[0].data, dtype=np.float64)
    return ap.image.PSF_Image(data=psf_data, pixelscale=pixelscale)


def sigma_to_variance(sigma):
    variance = sigma**2
    bad = ~np.isfinite(variance) | (variance <= 0)
    if np.any(bad):
        good = variance[np.isfinite(variance) & (variance > 0)]
        fallback = np.nanmax(good) * 10 if len(good) > 0 else 1.0
        variance[bad] = fallback
    return variance


def main():
    baseline = pd.read_csv("AstroPhot_Analysis/results/astrophot_psf_sigma_inclination_angles.csv")
    scales = pd.read_csv("AstroPhot_Analysis/results/highchi2_sigma_scale_recommendations.csv")

    scale_map = dict(zip(scales["frb_name"], scales["sigma_scale_per_object"]))
    baseline_map = dict(zip(baseline["frb_name"], baseline["chi2_nu"]))

    out_rows = []
    scaling_rows = []
    scaled_sigma_dir = "tools/astrophot/results/host_sigma_scaled_highchi2"
    os.makedirs(scaled_sigma_dir, exist_ok=True)

    for frb_name, sigma_scale in scale_map.items():
        print(f"\n{'='*60}")
        print(f"Processing {frb_name} with sigma scale {sigma_scale:.6f}")

        flux, sigma, pixelscale, flux_header, sigma_header = load_flux_sigma_pixelscale(frb_name)
        sigma_scaled = sigma * float(sigma_scale)

        # Save scaled sigma map for reproducibility.
        out_sigma_path = os.path.join(scaled_sigma_dir, f"{frb_name}_sigma_scaled.fits")
        fits.PrimaryHDU(data=sigma_scaled, header=sigma_header).writeto(out_sigma_path, overwrite=True)

        variance_scaled = sigma_to_variance(sigma_scaled)
        target = ap.image.Target_Image(data=flux, pixelscale=pixelscale, variance=variance_scaled)

        model = ap.models.AstroPhot_Model(
            name="model with target",
            model_type="sersic galaxy model",
            target=target,
            psf=load_psf_image(frb_name, pixelscale),
        )
        model.psf_mode = "full"
        model.initialize()

        fitter = ap.fit.LM(model, verbose=1)
        result = fitter.fit()

        q = model["q"].value.item()
        incl = np.degrees(np.arccos(np.clip(q, -1.0, 1.0)))
        nval = model["n"].value.item()
        re = model["Re"].value.item()
        pa = model["PA"].value.item()

        chi2 = float(fitter.res_loss())
        ndf = int(fitter.ndf)
        chi2_nu = chi2 / ndf if ndf else np.nan
        chi2_old = float(baseline_map.get(frb_name, np.nan))

        print(
            f"  chi2_nu old={chi2_old:.6f} new={chi2_nu:.6f} "
            f"(ratio {chi2_nu/chi2_old if np.isfinite(chi2_old) and chi2_old!=0 else np.nan:.6e})"
        )

        out_rows.append(
            {
                "frb_name": frb_name,
                "sigma_scale": float(sigma_scale),
                "chi2_nu_old": chi2_old,
                "chi2_nu_new": chi2_nu,
                "chi2_nu_ratio_new_over_old": chi2_nu / chi2_old if np.isfinite(chi2_old) and chi2_old != 0 else np.nan,
                "q": q,
                "inclination_angle": incl,
                "n": nval,
                "PA": pa,
                "Re": re,
                "chi2": chi2,
                "ndf": ndf,
                "fit_message": result.message,
            }
        )

        scaling_rows.append(
            {
                "frb_name": frb_name,
                "sigma_scale": float(sigma_scale),
                "sigma_median_original": float(np.nanmedian(sigma[np.isfinite(sigma) & (sigma > 0)])),
                "sigma_median_scaled": float(np.nanmedian(sigma_scaled[np.isfinite(sigma_scaled) & (sigma_scaled > 0)])),
                "scaled_sigma_path": out_sigma_path,
            }
        )

    results_df = pd.DataFrame(out_rows).sort_values("chi2_nu_new", ascending=False)
    scaling_df = pd.DataFrame(scaling_rows)

    results_path = "AstroPhot_Analysis/results/highchi2_refit_scaled_sigma_results.csv"
    scaling_path = "AstroPhot_Analysis/results/highchi2_sigma_scalings_applied.csv"
    summary_path = "AstroPhot_Analysis/results/highchi2_refit_scaled_sigma_summary.csv"

    results_df.to_csv(results_path, index=False)
    scaling_df.to_csv(scaling_path, index=False)

    summary = pd.DataFrame(
        [
            {
                "n_objects": len(results_df),
                "median_chi2_nu_old": float(results_df["chi2_nu_old"].median()),
                "median_chi2_nu_new": float(results_df["chi2_nu_new"].median()),
                "max_chi2_nu_new": float(results_df["chi2_nu_new"].max()),
                "min_chi2_nu_new": float(results_df["chi2_nu_new"].min()),
                "mean_abs_log10_ratio": float(np.mean(np.abs(np.log10(results_df["chi2_nu_ratio_new_over_old"])))),
                "n_with_chi2_nu_new_below_5": int((results_df["chi2_nu_new"] < 5).sum()),
            }
        ]
    )
    summary.to_csv(summary_path, index=False)

    print("\nSaved", results_path)
    print("Saved", scaling_path)
    print("Saved", summary_path)
    print("\nBefore vs after chi2_nu:")
    print(results_df[["frb_name", "sigma_scale", "chi2_nu_old", "chi2_nu_new", "chi2_nu_ratio_new_over_old"]].to_string(index=False))


if __name__ == "__main__":
    main()
