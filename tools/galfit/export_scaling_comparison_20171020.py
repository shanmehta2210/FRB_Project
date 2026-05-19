from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits


def _normalize_with_image_scale(arr: np.ndarray, vmin: float, vmax: float) -> np.ndarray:
    arr = np.nan_to_num(arr, nan=vmin, posinf=vmax, neginf=vmin)
    arr = np.flipud(arr)  # Match report orientation used in the current workflow.
    denom = vmax - vmin
    if denom <= 0:
        return np.zeros_like(arr, dtype=float)
    return np.clip((arr - vmin) / denom, 0.0, 1.0)


def _stretch(x: np.ndarray, kind: str) -> np.ndarray:
    if kind == "linear":
        y = x
    elif kind == "sqrt":
        y = np.sqrt(x)
    elif kind == "log":
        y = np.log1p(1000.0 * x) / np.log1p(1000.0)
    elif kind == "asinh":
        y = np.arcsinh(10.0 * x) / np.arcsinh(10.0)
    elif kind == "sinh":
        y = np.sinh(3.0 * x) / np.sinh(3.0)
    elif kind == "power2":
        y = x**2
    else:
        raise ValueError(f"Unsupported scaling: {kind}")
    return np.clip(y, 0.0, 1.0)


def _save_gray(path: Path, arr: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.imsave(path, arr, cmap="gray", vmin=0.0, vmax=1.0)


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    out_fits = root / "tools" / "galfit" / "runs" / "20171020A" / "with_psf_sigma" / "out.fits"
    out_dir = (
        root
        / "Reports"
        / "00 Galfit (AS) verification"
        / "scaling_tests"
        / "20171020A_with_psf_sigma"
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    with fits.open(out_fits) as hdul:
        image = np.asarray(hdul[1].data, dtype=float)
        model = np.asarray(hdul[2].data, dtype=float)
        residual = np.asarray(hdul[3].data, dtype=float)

    finite = np.isfinite(image)
    if not finite.any():
        raise RuntimeError("Image has no finite values; cannot derive fixed scale.")

    vmin = float(np.nanmin(image[finite]))
    vmax = float(np.nanmax(image[finite]))

    image_norm = _normalize_with_image_scale(image, vmin=vmin, vmax=vmax)
    model_norm = _normalize_with_image_scale(model, vmin=vmin, vmax=vmax)
    residual_norm = _normalize_with_image_scale(residual, vmin=vmin, vmax=vmax)

    scalings = ["linear", "sqrt", "log", "asinh", "sinh", "power2"]

    # Per-scaling IMR triplets for close inspection.
    for s in scalings:
        img_s = _stretch(image_norm, s)
        mod_s = _stretch(model_norm, s)
        res_s = _stretch(residual_norm, s)

        _save_gray(out_dir / f"{s}_image.png", img_s)
        _save_gray(out_dir / f"{s}_model.png", mod_s)
        _save_gray(out_dir / f"{s}_residual.png", res_s)

        fig, axes = plt.subplots(1, 3, figsize=(9, 3))
        for ax, arr, title in zip(
            axes,
            [img_s, mod_s, res_s],
            ["Image", "Model", "Residual"],
            strict=True,
        ):
            ax.imshow(arr, cmap="gray", vmin=0.0, vmax=1.0)
            ax.set_title(title, fontsize=10)
            ax.axis("off")
        fig.suptitle(f"20171020A with_psf_sigma - {s}", fontsize=11)
        fig.tight_layout()
        fig.savefig(out_dir / f"{s}_imr_panel.png", dpi=200, bbox_inches="tight")
        plt.close(fig)

    # Overview sheet across all scalings.
    fig, axes = plt.subplots(len(scalings), 3, figsize=(9, 2.4 * len(scalings)))
    for i, s in enumerate(scalings):
        img_s = _stretch(image_norm, s)
        mod_s = _stretch(model_norm, s)
        res_s = _stretch(residual_norm, s)

        for j, arr in enumerate([img_s, mod_s, res_s]):
            axes[i, j].imshow(arr, cmap="gray", vmin=0.0, vmax=1.0)
            axes[i, j].axis("off")
            if i == 0:
                axes[i, j].set_title(["Image", "Model", "Residual"][j], fontsize=10)
            if j == 0:
                axes[i, j].set_ylabel(s, rotation=0, labelpad=25, fontsize=10, va="center")

    fig.suptitle("20171020A with_psf_sigma - Scaling Comparison", fontsize=12)
    fig.tight_layout()
    fig.savefig(out_dir / "20171020A_with_psf_sigma_scaling_overview.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved scaling test outputs to: {out_dir}")
    print("Scalings:", ", ".join(scalings))


if __name__ == "__main__":
    main()
