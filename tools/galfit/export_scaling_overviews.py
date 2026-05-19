from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits


FRB_IDS = ["20171020A", "20191001A", "20210807D"]
RUN_CASE = "with_psf_sigma"
SCALINGS = ["linear", "sqrt", "log", "asinh", "sinh", "power2"]


def _normalize_with_image_scale(arr: np.ndarray, vmin: float, vmax: float) -> np.ndarray:
    arr = np.nan_to_num(arr, nan=vmin, posinf=vmax, neginf=vmin)
    arr = np.flipud(arr)  # Keep report orientation convention.
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


def _build_overview(out_fits: Path, out_png: Path, title: str) -> None:
    with fits.open(out_fits) as hdul:
        image = np.asarray(hdul[1].data, dtype=float)
        model = np.asarray(hdul[2].data, dtype=float)
        residual = np.asarray(hdul[3].data, dtype=float)

    finite = np.isfinite(image)
    if not finite.any():
        raise RuntimeError(f"Image has no finite values: {out_fits}")

    vmin = float(np.nanmin(image[finite]))
    vmax = float(np.nanmax(image[finite]))

    image_norm = _normalize_with_image_scale(image, vmin=vmin, vmax=vmax)
    model_norm = _normalize_with_image_scale(model, vmin=vmin, vmax=vmax)
    residual_norm = _normalize_with_image_scale(residual, vmin=vmin, vmax=vmax)

    fig, axes = plt.subplots(len(SCALINGS), 3, figsize=(9, 2.4 * len(SCALINGS)))
    for i, s in enumerate(SCALINGS):
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

    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    runs_root = root / "tools" / "galfit" / "runs"
    tests_root = root / "Reports" / "00 Galfit (AS) verification" / "scaling_tests"

    for frb in FRB_IDS:
        out_fits = runs_root / frb / RUN_CASE / "out.fits"
        if not out_fits.exists():
            print(f"[skip] {frb}: missing {out_fits}")
            continue

        out_dir = tests_root / f"{frb}_{RUN_CASE}"
        out_dir.mkdir(parents=True, exist_ok=True)

        # Keep only the final overview output in each folder.
        for existing in out_dir.iterdir():
            if existing.is_file():
                existing.unlink()

        out_png = out_dir / f"{frb}_{RUN_CASE}_scaling_overview.png"
        _build_overview(
            out_fits=out_fits,
            out_png=out_png,
            title=f"{frb} {RUN_CASE} - Scaling Comparison",
        )
        print(f"[ok] {frb}: {out_png}")


if __name__ == "__main__":
    main()
