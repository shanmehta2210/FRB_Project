#!/usr/bin/env python3
"""Generate GalSim mock exponential disks on a rigid parameter grid."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import galsim
import numpy as np
from astropy.io import fits

from sim_utils import (
    TOOL_DIR,
    ensure_output_layout,
    flux_e_from_mag,
    galfit_center_1based,
    iter_grid_points,
    load_config,
    outputs_dir,
    re_pix_from_arcsec,
    realization_seed,
    resolve_coadd_exptime,
    resolve_zeropoint_e,
    sky_e_per_pix,
)


def build_psf(cfg: dict) -> galsim.Moffat:
    p = cfg["physics"]
    return galsim.Moffat(
        beta=float(p["psf_moffat_beta"]),
        fwhm=float(p["psf_fwhm_arcsec"]),
    )


def build_galaxy(
    *,
    re_arcsec: float,
    flux_e: float,
    ba: float,
    pa_deg: float,
) -> galsim.GSObject:
    gal = galsim.Exponential(half_light_radius=re_arcsec, flux=flux_e)
    return gal.shear(q=ba, beta=pa_deg * galsim.degrees)


def sample_snr_at_re(
    *,
    re_arcsec: float,
    mag: float,
    zeropoint_e: float,
    pixel_scale: float,
    stamp_px: int,
    sky_mag_arcsec2: float,
    psf: galsim.Moffat,
    exptime: float = 1.0,
    ba: float = 0.9,
    pa_deg: float = 30.0,
) -> float:
    """SNR at Re along the major axis after coadd exposure (Poisson)."""
    flux_e = flux_e_from_mag(mag, zeropoint_e)
    gal = build_galaxy(re_arcsec=re_arcsec, flux_e=flux_e, ba=ba, pa_deg=pa_deg)
    obs = galsim.Convolve(gal, psf)
    img = obs.drawImage(
        nx=stamp_px,
        ny=stamp_px,
        scale=pixel_scale,
        method="no_pixel",
    )
    arr = np.asarray(img.array, dtype=float)
    sky_rate = sky_e_per_pix(sky_mag_arcsec2, zeropoint_e, pixel_scale)

    xc, yc = galfit_center_1based(stamp_px)
    cx0 = xc - 1.0
    cy0 = yc - 1.0
    re_pix = re_arcsec / pixel_scale
    pa_rad = math.radians(pa_deg)
    dx = re_pix * math.cos(pa_rad)
    dy = re_pix * math.sin(pa_rad)
    ix = int(round(cx0 + dx))
    iy = int(round(cy0 + dy))
    ix = min(max(ix, 0), stamp_px - 1)
    iy = min(max(iy, 0), stamp_px - 1)

    obj_rate = max(arr[iy, ix], 0.0)
    signal = exptime * obj_rate
    noise = math.sqrt(max(exptime * (obj_rate + sky_rate), 1e-30))
    return signal / noise


def calibrate_coadd_exptime(cfg: dict, zeropoint_e: float) -> float:
    """Scale coadd exposure (not ZP) until SNR@Re matches target at calibrate_mag."""
    p = cfg["physics"]
    cal_mag = float(p.get("calibrate_mag", 20.0))
    cal_re = float(p.get("calibrate_re_arcsec", 1.0))
    target = float(p.get("target_snr_at_re", 1.0))
    psf = build_psf(cfg)

    lo, hi = 1.0, 1.0e6
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        snr = sample_snr_at_re(
            re_arcsec=cal_re,
            mag=cal_mag,
            zeropoint_e=zeropoint_e,
            pixel_scale=float(p["pixel_scale"]),
            stamp_px=int(p["stamp_px"]),
            sky_mag_arcsec2=float(p["sky_mag_arcsec2"]),
            psf=psf,
            exptime=mid,
        )
        if snr > target:
            hi = mid
        else:
            lo = mid
    t = 0.5 * (lo + hi)
    sky_pix = sky_e_per_pix(float(p["sky_mag_arcsec2"]), zeropoint_e, float(p["pixel_scale"]))
    print(
        f"Calibrated coadd_exptime_sec={t:.2f} "
        f"(SNR@Re={sample_snr_at_re(re_arcsec=cal_re, mag=cal_mag, zeropoint_e=zeropoint_e, pixel_scale=float(p['pixel_scale']), stamp_px=int(p['stamp_px']), sky_mag_arcsec2=float(p['sky_mag_arcsec2']), psf=psf, exptime=t):.3f} "
        f"at mag={cal_mag}, Re={cal_re}\"; sky={sky_pix:.4f} e-/pix/s, sky*{t:.0f}={sky_pix*t:.2f} e-/pix)"
    )
    return t


def write_psf_fits(psf: galsim.Moffat, path: Path, pixel_scale: float) -> None:
    psf_img = psf.drawImage(scale=pixel_scale, method="auto")
    arr = np.asarray(psf_img.array, dtype=np.float64)
    total = arr.sum()
    if total > 0:
        arr /= total
    hdu = fits.PrimaryHDU(arr)
    hdu.header["CDELT1"] = pixel_scale
    hdu.header["CDELT2"] = pixel_scale
    hdu.writeto(path, overwrite=True)


def write_mock_fits(
    path: Path,
    data: np.ndarray,
    pixel_scale: float,
    exptime: float,
    zeropoint_e: float,
) -> None:
    hdu = fits.PrimaryHDU(np.asarray(data, dtype=np.float32))
    hdu.header["GAIN"] = 1.0
    hdu.header["EXPTIME"] = exptime
    hdu.header["ZP"] = zeropoint_e
    hdu.header["CDELT1"] = pixel_scale
    hdu.header["CDELT2"] = pixel_scale
    hdu.writeto(path, overwrite=True)


def generate_one(
    *,
    point: dict,
    cfg: dict,
    zeropoint_e: float,
    coadd_exptime: float,
    psf: galsim.Moffat,
    mock_dir: Path,
) -> dict:
    p = cfg["physics"]
    pixel_scale = float(p["pixel_scale"])
    stamp_px = int(p["stamp_px"])
    sky_rate = sky_e_per_pix(float(p["sky_mag_arcsec2"]), zeropoint_e, pixel_scale)
    sky_e = sky_rate * coadd_exptime

    flux_e = flux_e_from_mag(point["mag_true"], zeropoint_e)
    gal = build_galaxy(
        re_arcsec=point["re_arcsec_true"],
        flux_e=flux_e,
        ba=point["ba_true"],
        pa_deg=point["pa_true"],
    )
    obs = galsim.Convolve(gal, psf)

    base_seed = int(p["rng_seed"])
    seed = realization_seed(base_seed, point["galaxy_id"], point["realization"])
    rng = galsim.BaseDeviate(seed)

    noiseless = obs.drawImage(nx=stamp_px, ny=stamp_px, scale=pixel_scale, method="auto")
    obj_arr = np.asarray(noiseless.array, dtype=float)
    obj_plus_sky = coadd_exptime * (obj_arr + sky_rate)

    noise = galsim.PoissonNoise(rng)
    noisy = galsim.Image(obj_plus_sky, wcs=noiseless.wcs)
    noisy.addNoise(noise)
    noisy_arr = np.asarray(noisy.array, dtype=float)

    # GALFIT cannot handle free sky + free mag when sky ~ 26 e-/pix (numerical
    # blow-up, mag snaps to ~25). Write sky-subtracted stamps for fitting; sky
    # component is seeded at 0 but left free (over-estimation trap on faint wings).
    fit_arr = noisy_arr - sky_e
    sigma_arr = np.sqrt(np.maximum(coadd_exptime * (obj_arr + sky_rate), 0.0))

    gid = point["galaxy_id"]
    out_sub = mock_dir / gid
    out_sub.mkdir(parents=True, exist_ok=True)
    write_mock_fits(out_sub / "mock.fits", fit_arr, pixel_scale, coadd_exptime, zeropoint_e)
    write_mock_fits(out_sub / "sigma.fits", sigma_arr, pixel_scale, coadd_exptime, zeropoint_e)

    xc, yc = galfit_center_1based(stamp_px)
    re_pix = re_pix_from_arcsec(point["re_arcsec_true"], pixel_scale)
    snr = sample_snr_at_re(
        re_arcsec=point["re_arcsec_true"],
        mag=point["mag_true"],
        zeropoint_e=zeropoint_e,
        pixel_scale=pixel_scale,
        stamp_px=stamp_px,
        sky_mag_arcsec2=float(p["sky_mag_arcsec2"]),
        psf=psf,
        exptime=coadd_exptime,
        ba=point["ba_true"],
        pa_deg=point["pa_true"],
    )

    return {
        "galaxy_id": gid,
        "realization": point["realization"],
        "ba_true": point["ba_true"],
        "re_arcsec_true": point["re_arcsec_true"],
        "re_pix_true": round(re_pix, 4),
        "mag_true": point["mag_true"],
        "pa_true": point["pa_true"],
        "flux_e": round(flux_e, 6),
        "sky_e_per_pix": round(sky_e, 8),
        "sky_rate_per_pix": round(sky_rate, 8),
        "sky_fit_seed": 0.0,
        "mock_sky_subtracted": True,
        "coadd_exptime_sec": round(coadd_exptime, 4),
        "zeropoint_e": round(zeropoint_e, 4),
        "xc": round(xc, 4),
        "yc": round(yc, 4),
        "snr_at_re": round(snr, 4),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=TOOL_DIR / "config.yaml",
        help="Path to config.yaml",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Generate at most N grid points (0 = all)",
    )
    args = parser.parse_args()

    cfg = load_config(args.config)
    layout = ensure_output_layout(cfg)
    zeropoint_e = resolve_zeropoint_e(cfg)
    coadd_exptime = resolve_coadd_exptime(cfg, lambda c: calibrate_coadd_exptime(c, zeropoint_e))
    psf = build_psf(cfg)
    pixel_scale = float(cfg["physics"]["pixel_scale"])

    psf_path = layout["mocks"] / "psf.fits"
    if not psf_path.is_file():
        write_psf_fits(psf, psf_path, pixel_scale)
        print(f"Wrote shared PSF -> {psf_path}")

    rows: list[dict] = []
    for i, point in enumerate(iter_grid_points(cfg)):
        if args.limit and i >= args.limit:
            break
        row = generate_one(
            point=point,
            cfg=cfg,
            zeropoint_e=zeropoint_e,
            coadd_exptime=coadd_exptime,
            psf=psf,
            mock_dir=layout["mocks"],
        )
        rows.append(row)
        if (i + 1) % 25 == 0 or args.limit == 1:
            print(f"Generated {i + 1}: {point['galaxy_id']}")

    catalog_path = layout["catalogs"] / "truth_catalog.csv"
    fieldnames = [
        "galaxy_id",
        "realization",
        "ba_true",
        "re_arcsec_true",
        "re_pix_true",
        "mag_true",
        "pa_true",
        "flux_e",
        "sky_e_per_pix",
        "sky_rate_per_pix",
        "sky_fit_seed",
        "mock_sky_subtracted",
        "coadd_exptime_sec",
        "zeropoint_e",
        "xc",
        "yc",
        "snr_at_re",
    ]
    with open(catalog_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} entries -> {catalog_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
