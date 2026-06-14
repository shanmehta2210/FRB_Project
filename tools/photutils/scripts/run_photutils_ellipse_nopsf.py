import argparse
import os
import warnings

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.convolution import convolve
from photutils.isophote import Ellipse, EllipseGeometry
from photutils.morphology import data_properties
from photutils.isophote import build_ellipse_model
from photutils.segmentation import make_2dgaussian_kernel, detect_sources

warnings.filterwarnings("ignore")


def hubble_inclination_from_q(q, q0=0.2):
    if not np.isfinite(q) or not np.isfinite(q0):
        return np.nan
    if q <= q0:
        return 90.0
    if q >= 1.0:
        return 0.0
    val = (q * q - q0 * q0) / (1.0 - q0 * q0)
    val = np.clip(val, 0.0, 1.0)
    return float(np.degrees(np.arccos(np.sqrt(val))))


def hubble_inclination_err_from_qerr(q, q_err, q0=0.2):
    if not np.isfinite(q) or not np.isfinite(q_err) or not np.isfinite(q0):
        return np.nan
    if q_err <= 0:
        return 0.0
    q_lo = max(0.0, q - q_err)
    q_hi = min(1.0, q + q_err)
    i_lo = hubble_inclination_from_q(q_lo, q0=q0)
    i_hi = hubble_inclination_from_q(q_hi, q0=q0)
    if not np.isfinite(i_lo) or not np.isfinite(i_hi):
        return np.nan
    return abs(i_hi - i_lo) / 2.0


def choose_half_light_isophote(isolist):
    # Archive notebook idea: use half-light isophote from cumulative isophotal flux.
    tflux = np.asarray(isolist.tflux_e, dtype=float)
    sma = np.asarray(isolist.sma, dtype=float)
    eps = np.asarray(isolist.eps, dtype=float)
    eps_err = np.asarray(isolist.ellip_err, dtype=float)
    pa = np.asarray(isolist.pa, dtype=float)
    pa_err = np.asarray(isolist.pa_err, dtype=float)

    valid = np.isfinite(tflux) & np.isfinite(sma) & np.isfinite(eps) & (eps >= 0) & (eps < 0.95)
    if not np.any(valid):
        return None

    idx = np.where(valid)[0]
    t = tflux[idx]
    # Ensure monotonic behavior for stable half-light selection.
    t = np.maximum.accumulate(t)
    half = 0.5 * t[-1]
    k = int(np.searchsorted(t, half, side="left"))
    k = min(max(k, 0), len(idx) - 1)
    best = idx[k]

    best_sma = sma[best]
    
    # Implementing median isophotal averaging to wash out clump-induced spikes
    mask = valid & (sma >= 0.5 * best_sma) & (sma <= 1.5 * best_sma)
    if np.count_nonzero(mask) >= 3:
        # Use median of the annulus to ignore localized clumps/knots
        med_eps = np.median(eps[mask])
        med_pa = np.median(pa[mask])
        med_eps_err = np.median(eps_err[mask])
        med_pa_err = np.median(pa_err[mask])
    else:
        # Fallback to single half-light proxy if annulus is too narrow
        med_eps = eps[best]
        med_pa = pa[best]
        med_eps_err = eps_err[best]
        med_pa_err = pa_err[best]

    return {
        "sma_pix": float(best_sma),
        "eps": float(med_eps),
        "eps_err": float(med_eps_err) if np.isfinite(med_eps_err) else np.nan,
        "pa_rad": float(med_pa),
        "pa_err_rad": float(med_pa_err) if np.isfinite(med_pa_err) else np.nan,
        "n_isophotes": int(np.count_nonzero(valid)),
        "rep_index": int(best),
        "half_light_flux_proxy": float(half),
    }


def fit_single_frb(frb_name, q0_value, galaxy_dir, minsma, step, maxsma_frac):
    if galaxy_dir == "pipeline":
        run_dir = os.path.join("pipeline_scripts", "Output", f"{frb_name}_all")
        flux_path = os.path.join(run_dir, "host_cutout.fits")
        sigma_path = os.path.join(run_dir, "host_sigma.fits")
    else:
        flux_path = os.path.join(galaxy_dir, f"{frb_name}_flux.fits")
        sigma_path = os.path.join(galaxy_dir, f"{frb_name}_sigma.fits")

    if not os.path.exists(flux_path):
        return {"FRB": frb_name, "status": "missing_flux_file"}
    if not os.path.exists(sigma_path):
        return {"FRB": frb_name, "status": "missing_sigma_file"}

    data = np.squeeze(np.asarray(fits.getdata(flux_path), dtype=float))
    sigma = np.squeeze(np.asarray(fits.getdata(sigma_path), dtype=float))

    if data.ndim != 2:
        return {"FRB": frb_name, "status": f"invalid_flux_shape_{data.shape}"}
    if sigma.shape != data.shape:
        return {"FRB": frb_name, "status": f"shape_mismatch_flux_{data.shape}_sigma_{sigma.shape}"}

    finite = np.isfinite(data) & np.isfinite(sigma) & (sigma > 0)
    if np.count_nonzero(finite) < 50:
        return {"FRB": frb_name, "status": "too_few_valid_pixels"}

    _, bkg_med, bkg_std = sigma_clipped_stats(data[finite], sigma=3.0, maxiters=10)
    image = data - bkg_med
    image_fit = np.where(finite, image, 0.0)

    # Segmentation masking to exclude non-target sources (like stars or companions)
    bad_pixel_mask = np.zeros_like(image_fit, dtype=bool)
    try:
        kernel = make_2dgaussian_kernel(3.0, size=3)
        smoothed = convolve(image_fit, kernel)
        threshold = bkg_std * 2.5
        segm = detect_sources(smoothed, threshold, npixels=5)
        if segm is not None:
            # Find label closest to center
            ny, nx = image_fit.shape
            cx, cy = int(nx / 2), int(ny / 2)
            central_labels = segm.data[max(0, cy-5):min(ny, cy+6), max(0, cx-5):min(nx, cx+6)]
            vals = central_labels[central_labels > 0]
            if len(vals) > 0:
                target_label = np.bincount(vals.flatten()).argmax()
                bad_pixel_mask = (segm.data > 0) & (segm.data != target_label)
                # Hide companions by replacing with background level (0.0 after subtraction)
                image_fit[bad_pixel_mask] = 0.0
                print(f"Masked {np.count_nonzero(bad_pixel_mask)} pixels for {frb_name}")
    except Exception as e:
        print(f"Segmentation masking failed for {frb_name}: {e}")

    # Robust initial values inspired by archive notebook manual-geometry workflow.
    props = data_properties(np.nan_to_num(image_fit, nan=0.0))
    x0 = float(getattr(props.xcentroid, "value", props.xcentroid))
    y0 = float(getattr(props.ycentroid, "value", props.ycentroid))
    a_sig = float(getattr(props.semimajor_sigma, "value", props.semimajor_sigma))
    b_sig = float(getattr(props.semiminor_sigma, "value", props.semiminor_sigma))
    if a_sig <= 0 or b_sig <= 0:
        return {"FRB": frb_name, "status": "invalid_moments"}

    eps0 = float(np.clip(1.0 - (b_sig / a_sig), 0.02, 0.90))
    pa0 = float(getattr(props.orientation, "value", props.orientation))
    sma0 = float(max(minsma + 0.5, a_sig))

    ny, nx = image.shape
    maxsma = float(max(minsma + 2.0, min(nx, ny) * maxsma_frac))

    # Primary adaptive strategy bank for most galaxies.
    strategies = [
        {"sma": sma0, "eps": eps0, "pa": pa0, "step": step, "fix_center": True, "fix_pa": False, "fix_eps": False},
        {"sma": max(minsma + 0.5, 1.5 * sma0), "eps": eps0, "pa": pa0, "step": 0.06, "fix_center": True, "fix_pa": False, "fix_eps": False},
        {"sma": max(minsma + 0.5, 0.7 * sma0), "eps": min(0.8, max(0.05, eps0)), "pa": pa0, "step": 0.12, "fix_center": True, "fix_pa": True, "fix_eps": False},
        {"sma": max(minsma + 1.0, min(image.shape) * 0.12), "eps": 0.35, "pa": 0.0, "step": 0.08, "fix_center": True, "fix_pa": False, "fix_eps": False},
        {"sma": max(minsma + 1.0, min(image.shape) * 0.18), "eps": 0.5, "pa": np.pi / 4.0, "step": 0.08, "fix_center": False, "fix_pa": False, "fix_eps": False},
    ]

    isolist = None
    used_strategy = -1
    last_exc = ""
    for i, st in enumerate(strategies):
        geom = EllipseGeometry(x0=x0, y0=y0, sma=st["sma"], eps=st["eps"], pa=st["pa"])
        ell = Ellipse(image_fit, geometry=geom)
        try:
            cur = ell.fit_image(
                minsma=minsma,
                maxsma=maxsma,
                step=st["step"],
                nclip=3,
                fix_center=st["fix_center"],
                fix_pa=st["fix_pa"],
                fix_eps=st["fix_eps"],
            )
            if len(cur) > 0:
                isolist = cur
                used_strategy = i
                break
        except Exception as exc:
            last_exc = str(exc)

    # Bounded fallback grid for stubborn objects that still fail the primary bank.
    if isolist is None:
        fallback_sma = [
            max(minsma + 0.1, sma0),
            max(minsma + 0.1, 1.5 * sma0),
            max(minsma + 0.1, 0.7 * sma0),
            max(minsma + 0.5, min(image.shape) * 0.12),
            max(minsma + 1.0, min(image.shape) * 0.18),
            max(minsma + 0.1, 4.0),
            max(minsma + 0.1, 6.0),
            max(minsma + 0.1, 8.0),
        ]
        fallback_eps = [eps0, 0.05, 0.15, 0.25, 0.35, 0.5]
        fallback_pa = [pa0, 0.0, np.pi / 4.0, np.pi / 2.0]
        fallback_steps = [0.03, 0.05]
        fallback_flags = [
            (False, False, False),
            (True, True, False),
            (True, False, False),
        ]

        fallback_idx = 0
        for sma_try in fallback_sma:
            for eps_try in fallback_eps:
                eps_try = float(np.clip(eps_try, 0.02, 0.90))
                for pa_try in fallback_pa:
                    for step_try in fallback_steps:
                        for fix_center, fix_pa, fix_eps in fallback_flags:
                            geom = EllipseGeometry(
                                x0=x0,
                                y0=y0,
                                sma=float(sma_try),
                                eps=eps_try,
                                pa=float(pa_try),
                            )
                            ell = Ellipse(image_fit, geometry=geom)
                            try:
                                cur = ell.fit_image(
                                    minsma=minsma,
                                    maxsma=maxsma,
                                    step=step_try,
                                    nclip=3,
                                    fix_center=fix_center,
                                    fix_pa=fix_pa,
                                    fix_eps=fix_eps,
                                )
                                if len(cur) > 0:
                                    isolist = cur
                                    used_strategy = len(strategies) + fallback_idx
                                    break
                            except Exception as exc:
                                last_exc = str(exc)
                            fallback_idx += 1
                        if isolist is not None:
                            break
                    if isolist is not None:
                        break
                if isolist is not None:
                    break
            if isolist is not None:
                break

    if isolist is None:
        if last_exc:
            return {"FRB": frb_name, "status": f"ellipse_fit_failed: {last_exc}"}
        return {"FRB": frb_name, "status": "no_isophotes"}

    rep = choose_half_light_isophote(isolist)
    if rep is None:
        return {"FRB": frb_name, "status": "no_valid_half_light_isophote"}

    q = 1.0 - rep["eps"]
    q_err = rep["eps_err"] if np.isfinite(rep["eps_err"]) else np.nan
    inc = hubble_inclination_from_q(q, q0=q0_value)
    inc_err = hubble_inclination_err_from_qerr(q, q_err, q0=q0_value)

    return {
        "FRB": frb_name,
        "status": "ok",
        "source_flux": flux_path,
        "source_sigma": sigma_path,
        "used_psf_convolution": False,
        "q0": float(q0_value),
        "b_over_a": float(q),
        "b_over_a_err": float(q_err) if np.isfinite(q_err) else np.nan,
        "inclination_deg": float(inc) if np.isfinite(inc) else np.nan,
        "inclination_err_deg": float(inc_err) if np.isfinite(inc_err) else np.nan,
        "sma_rep_pix": float(rep["sma_pix"]),
        "ellipticity_rep": float(rep["eps"]),
        "pa_deg": float(np.degrees(rep["pa_rad"])),
        "pa_err_deg": float(np.degrees(rep["pa_err_rad"])) if np.isfinite(rep["pa_err_rad"]) else np.nan,
        "x0_pix": float(x0),
        "y0_pix": float(y0),
        "n_isophotes": int(rep["n_isophotes"]),
        "isophote_index": int(rep["rep_index"]),
        "half_light_flux_proxy": float(rep["half_light_flux_proxy"]),
        "bkg_median": float(bkg_med),
        "bkg_std": float(bkg_std),
        "minsma": float(minsma),
        "maxsma": float(maxsma),
        "step": float(step),
        "fit_strategy_index": int(used_strategy),
    }


def main():
    parser = argparse.ArgumentParser(description="Run no-PSF photutils ellipse inclination analysis on pipeline host cutouts")
    parser.add_argument("--master-csv", default="master_frb_summary.csv")
    parser.add_argument(
        "--galaxy-dir",
        default="pipeline",
        help="Use 'pipeline' for pipeline_scripts/Output/<FRB>_all/host_cutout.fits (default), "
        "or a flat directory with legacy <FRB>_flux.fits naming",
    )
    parser.add_argument("--out-csv", default="tools/photutils/results/photutils_ellipse_inclination_angles_nopsf.csv")
    parser.add_argument("--minsma", type=float, default=2.0)
    parser.add_argument("--step", type=float, default=0.10)
    parser.add_argument("--maxsma-frac", type=float, default=0.45)
    args = parser.parse_args()

    master = pd.read_csv(args.master_csv)
    results = []

    print(f"FRBs in master: {len(master)}")
    print(f"Using host flux/sigma from: {args.galaxy_dir}")

    for _, row in master.iterrows():
        frb = str(row["FRB"])
        q0_val = float(row["q0"]) if ("q0" in row and np.isfinite(row["q0"])) else 0.2
        out = fit_single_frb(
            frb_name=frb,
            q0_value=q0_val,
            galaxy_dir=args.galaxy_dir,
            minsma=args.minsma,
            step=args.step,
            maxsma_frac=args.maxsma_frac,
        )
        results.append(out)
        print(f"{frb}: {out.get('status', 'unknown')}")

    df = pd.DataFrame(results)
    os.makedirs(os.path.dirname(args.out_csv), exist_ok=True)
    df.to_csv(args.out_csv, index=False)

    print("--- Summary ---")
    print(df["status"].value_counts(dropna=False))
    print(f"Saved {len(df)} rows to {args.out_csv}")


if __name__ == "__main__":
    main()
