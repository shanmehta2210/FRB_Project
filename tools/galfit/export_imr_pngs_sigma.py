from __future__ import annotations

from pathlib import Path
from typing import Dict

import numpy as np
from astropy.io import fits
from PIL import Image


UPSCALE_FACTOR = 4
DEFAULT_STRETCH = "linear"
ASINH_ALPHA = 10.0


def sanitize_for_png(data: np.ndarray, vmin: float, vmax: float) -> np.ndarray:
    """Replace non-finite values so PNG export remains valid."""
    return np.nan_to_num(data, nan=vmin, posinf=vmax, neginf=vmin)


def read_galfit_imr(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read GALFIT image/model/residual arrays from extensions 1/2/3."""
    with fits.open(path) as hdul:
        if len(hdul) < 4:
            raise ValueError(f"Expected >=4 HDUs in {path}, found {len(hdul)}")
        image = np.asarray(hdul[1].data, dtype=float)
        model = np.asarray(hdul[2].data, dtype=float)
        residual = np.asarray(hdul[3].data, dtype=float)
    return image, model, residual


def save_png(data: np.ndarray, out_path: Path, vmin: float, vmax: float) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    safe = sanitize_for_png(data, vmin=vmin, vmax=vmax)
    safe = np.flipud(safe)

    # Keep linear scaling but write a higher-resolution PNG for cleaner PDF rendering.
    denom = vmax - vmin
    if denom <= 0:
        norm = np.zeros_like(safe, dtype=np.uint8)
    else:
        scaled = (safe - vmin) / denom
        if DEFAULT_STRETCH == "asinh":
            scaled = np.arcsinh(ASINH_ALPHA * scaled) / np.arcsinh(ASINH_ALPHA)
        scaled = np.clip(scaled, 0.0, 1.0)
        norm = (scaled * 255.0).astype(np.uint8)

    image = Image.fromarray(norm, mode="L")
    # Use nearest-neighbor so panel details remain crisp in the PDF.
    upscaled = image.resize(
        (image.width * UPSCALE_FACTOR, image.height * UPSCALE_FACTOR),
        resample=Image.Resampling.NEAREST,
    )
    upscaled.save(out_path, optimize=True)


def main() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    galfit_root = repo_root / "tools" / "galfit"
    runs_root = galfit_root / "runs"
    out_root = galfit_root / "imr_pngs"
    out_root.mkdir(parents=True, exist_ok=True)

    sigma_outfits = list(runs_root.glob("*/with_psf_sigma/out.fits")) + list(
        runs_root.glob("*/no_psf_sigma/out.fits")
    )

    by_frb: Dict[str, Dict[str, Path]] = {}
    for out_fits in sigma_outfits:
        frb = out_fits.parents[1].name
        run_type = out_fits.parent.name
        by_frb.setdefault(frb, {})[run_type] = out_fits

    exported = 0
    for frb in sorted(by_frb):
        runs = by_frb[frb]
        frb_out_dir = out_root / frb

        # Use the sigma image scale (same image for with/without PSF runs).
        image_source = runs.get("with_psf_sigma") or runs.get("no_psf_sigma")
        if image_source is None:
            continue

        image, _, _ = read_galfit_imr(image_source)
        finite = np.isfinite(image)
        if not finite.any():
            print(f"[skip] {frb}: image has no finite pixels")
            continue

        vmin = float(np.nanmin(image[finite]))
        vmax = float(np.nanmax(image[finite]))
        if not np.isfinite(vmin) or not np.isfinite(vmax):
            print(f"[skip] {frb}: invalid min/max")
            continue
        if vmin == vmax:
            vmax = vmin + 1.0

        # Save image once per FRB.
        save_png(image, frb_out_dir / f"{frb}_image.png", vmin=vmin, vmax=vmax)

        if "with_psf_sigma" in runs:
            _, model, residual = read_galfit_imr(runs["with_psf_sigma"])
            save_png(
                model,
                frb_out_dir / f"{frb}_with_psf_model.png",
                vmin=vmin,
                vmax=vmax,
            )
            save_png(
                residual,
                frb_out_dir / f"{frb}_with_psf_residual.png",
                vmin=vmin,
                vmax=vmax,
            )

        if "no_psf_sigma" in runs:
            _, model, residual = read_galfit_imr(runs["no_psf_sigma"])
            save_png(
                model,
                frb_out_dir / f"{frb}_no_psf_model.png",
                vmin=vmin,
                vmax=vmax,
            )
            save_png(
                residual,
                frb_out_dir / f"{frb}_no_psf_residual.png",
                vmin=vmin,
                vmax=vmax,
            )

        exported += 1
        print(f"[ok] {frb}: exported PNG set to {frb_out_dir}")

    print(f"Completed. FRBs exported: {exported}")
    print(f"Output directory: {out_root}")


if __name__ == "__main__":
    main()
