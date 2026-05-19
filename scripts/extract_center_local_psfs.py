import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import psfex_local as psfex
from astropy.io import fits


def extract_local_psf(psf_binary_path, frb_x, frb_y, output_fits_path, size=25):
    """
    Reads a spatially varying .psf file and extracts the exact 2D PSF
    rendered at the FRB host's specific coordinates.
    """
    # 1. Load the spatially varying mathematical model
    psf_model = psfex.PSFEx(str(psf_binary_path))

    # 2. Render the PSF at the exact pixel location of the FRB, directly at
    # requested output size (native resampling path, no manual padding).
    # get_rec(y, x, out_size) returns the reconstructed 2D array.
    # Note: psfex typically expects (row, column) which is (Y, X)
    local_psf_array = psf_model.get_rec(float(frb_y), float(frb_x), out_size=int(size))

    if local_psf_array is None:
        raise RuntimeError(f"PSFEx returned None for ({frb_x}, {frb_y})")

    # 3. Validate fixed output shape.
    ny, nx = local_psf_array.shape
    out_size = int(size)
    if (ny, nx) != (out_size, out_size):
        raise ValueError(f"Expected reconstructed shape {(out_size, out_size)}, got {(ny, nx)}")

    fixed_psf = local_psf_array.astype(np.float64)

    # 4. Normalize it (Crucial for GALFIT/AstroPhot)
    total = float(np.sum(fixed_psf))
    if not np.isfinite(total) or total == 0.0:
        raise ValueError("PSF sum is non-finite or zero; cannot normalize")
    fixed_psf = fixed_psf / total

    # 5. Save as a standard FITS file
    output_fits_path = Path(output_fits_path)
    output_fits_path.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(fixed_psf).writeto(output_fits_path, overwrite=True)

    return fixed_psf, f"render_size={out_size};native_recon={psf_model.recon_ny}x{psf_model.recon_nx}"


def save_png(array2d, png_path):
    png_path = Path(png_path)
    vmin = float(np.nanmin(array2d))
    vmax = float(np.nanmax(array2d))
    if not np.isfinite(vmin) or not np.isfinite(vmax):
        raise ValueError("PSF contains non-finite-only values")
    if vmax == vmin:
        vmax = vmin + 1e-12

    plt.figure(figsize=(4, 4))
    plt.imshow(array2d, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)
    plt.axis("off")
    plt.tight_layout(pad=0)
    plt.savefig(png_path, dpi=180, bbox_inches="tight", pad_inches=0)
    plt.close()


def main():
    repo_root = Path(__file__).resolve().parents[1]
    runs_root = repo_root / "psfs" / "PSFEx + SExtractor" / "runs"
    out_root = repo_root / "psfs" / "PSFEx + SExtractor" / "final_center_psfs"
    out_root.mkdir(parents=True, exist_ok=True)

    rows = []

    for run_dir in sorted([p for p in runs_root.iterdir() if p.is_dir()]):
        frb = run_dir.name
        psf_path = run_dir / "input" / "psfex" / "sextractor_catalog.psf"
        image_path = run_dir / "input" / "sextractor" / "image.fits"

        if not psf_path.exists() or not image_path.exists():
            rows.append(
                {
                    "frb": frb,
                    "status": "skipped_missing_input",
                    "message": "Missing .psf or image.fits",
                    "center_x": "",
                    "center_y": "",
                    "fits_out": "",
                    "png_out": "",
                }
            )
            continue

        try:
            with fits.open(image_path) as hdul:
                image_data = hdul[0].data
            if image_data is None or image_data.ndim < 2:
                raise ValueError("image.fits has no 2D data")

            height, width = int(image_data.shape[-2]), int(image_data.shape[-1])
            center_x = width / 2.0
            center_y = height / 2.0

            frb_out_dir = out_root / frb
            frb_out_dir.mkdir(parents=True, exist_ok=True)
            fits_out = frb_out_dir / f"{frb}_center_psf_25.fits"
            png_out = frb_out_dir / f"{frb}_center_psf_25.png"

            psf_arr, size_note = extract_local_psf(
                psf_binary_path=psf_path,
                frb_x=center_x,
                frb_y=center_y,
                output_fits_path=fits_out,
                size=25,
            )
            save_png(psf_arr, png_out)

            rows.append(
                {
                    "frb": frb,
                    "status": "ok",
                    "message": size_note,
                    "center_x": f"{center_x:.3f}",
                    "center_y": f"{center_y:.3f}",
                    "fits_out": str(fits_out.relative_to(repo_root)).replace("\\", "/"),
                    "png_out": str(png_out.relative_to(repo_root)).replace("\\", "/"),
                }
            )
            print(f"[ok] {frb}: center=({center_x:.3f},{center_y:.3f})")
        except Exception as exc:
            rows.append(
                {
                    "frb": frb,
                    "status": "error",
                    "message": str(exc),
                    "center_x": "",
                    "center_y": "",
                    "fits_out": "",
                    "png_out": "",
                }
            )
            print(f"[error] {frb}: {exc}")

    manifest_path = out_root / "final_center_psf_manifest.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "frb",
                "status",
                "message",
                "center_x",
                "center_y",
                "fits_out",
                "png_out",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    ok_count = sum(1 for r in rows if r["status"] == "ok")
    print(f"Wrote manifest: {manifest_path}")
    print(f"Success: {ok_count}/{len(rows)}")


if __name__ == "__main__":
    main()
