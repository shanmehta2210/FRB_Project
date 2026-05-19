import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits


def load_fits_image(path: Path) -> np.ndarray:
    with fits.open(path) as hdul:
        arr = np.array(hdul[0].data, dtype=float)
    if arr.ndim == 3:
        # Legacy cutouts may contain multiple planes; use the first plane for diagnostics.
        arr = arr[0]
    elif arr.ndim > 3:
        arr = np.squeeze(arr)
        while arr.ndim > 2:
            arr = arr[0]
    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    return arr


def load_galfit_hdu_image(path: Path, ext: int) -> np.ndarray:
    with fits.open(path) as hdul:
        arr = np.array(hdul[ext].data, dtype=float)
    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    return arr


def robust_limits(arr: np.ndarray) -> tuple[float, float]:
    lo = float(np.nanpercentile(arr, 2.0))
    hi = float(np.nanpercentile(arr, 98.0))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo = float(np.nanmin(arr))
        hi = float(np.nanmax(arr))
        if hi <= lo:
            hi = lo + 1.0
    return lo, hi


def apply_arcsin_stretch(arr: np.ndarray, ref_min: float, ref_max: float) -> np.ndarray:
    span = ref_max - ref_min
    if (not np.isfinite(span)) or span <= 0:
        span = 1.0
    scaled = 2.0 * ((arr - ref_min) / span) - 1.0
    scaled = np.clip(scaled, -1.0, 1.0)
    return np.arcsin(scaled)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--incl-csv", default="legacy_vs_galfit_inclination_comparison.csv")
    parser.add_argument("--legacy-imr-dir", default="tools/legacy/imr_fits")
    parser.add_argument("--galfit-imr-dir", default="tools/galfit/imr_pngs")
    parser.add_argument("--galfit-runs-dir", default="tools/galfit/runs")
    parser.add_argument("--galfit-run-subdir", default="with_psf_sigma", help="GALFIT run subdirectory to use (e.g., with_psf or with_psf_sigma)")
    parser.add_argument("--exclude-types", default="REX")
    parser.add_argument("--top-n", type=int, default=6)
    parser.add_argument("--out-dir", default="plots/plots_legacy_cdf")
    parser.add_argument("--tag", default="top_delta_i")
    args = parser.parse_args()

    df = pd.read_csv(args.incl_csv)
    df = df[df["ls_inc_deg"].notna() & df["galfit_inc_psf_deg"].notna()].copy()

    exclude_types = [t.strip().upper() for t in str(args.exclude_types).split(",") if t.strip()]
    if exclude_types and "type_ls" in df.columns:
        df = df[~df["type_ls"].astype(str).str.upper().isin(exclude_types)].copy()

    if "delta_deg_ls_minus_galfit" in df.columns:
        df["abs_delta_i_deg"] = pd.to_numeric(df["delta_deg_ls_minus_galfit"], errors="coerce").abs()
    else:
        di = pd.to_numeric(df["ls_inc_deg"], errors="coerce") - pd.to_numeric(df["galfit_inc_psf_deg"], errors="coerce")
        df["abs_delta_i_deg"] = di.abs()
        df["delta_deg_ls_minus_galfit"] = di

    def has_complete_imr(frb_name: str) -> bool:
        legacy_dir = Path(args.legacy_imr_dir) / frb_name
        legacy_ok = all(
            (legacy_dir / f"{frb_name}_{suffix}.fits").exists()
            for suffix in ["image", "model", "resid"]
        )
        galfit_out = Path(args.galfit_runs_dir) / frb_name / args.galfit_run_subdir / "out.fits"
        return legacy_ok and galfit_out.exists()

    df = df[df["FRB"].astype(str).map(has_complete_imr)].copy()

    df = df.sort_values("abs_delta_i_deg", ascending=False).head(args.top_n).copy()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    list_csv = out_dir / f"top_delta_i_frbs_{args.tag}.csv"
    df.to_csv(list_csv, index=False)

    n_rows = len(df)
    if n_rows == 0:
        print("No rows available for panel generation.")
        return

    fig, axes = plt.subplots(n_rows, 6, figsize=(18, 3.2 * n_rows))
    if n_rows == 1:
        axes = np.array([axes])

    panel_cmap = "gray"

    col_titles = ["Legacy Image", "Legacy Model", "Legacy Residual", "GALFIT Image", "GALFIT Model", "GALFIT Residual"]
    for j, title in enumerate(col_titles):
        axes[0, j].set_title(title, fontsize=10)

    # Keep the panel uncluttered; per-row FRB labels are rendered inside each row.

    for i, row in enumerate(df.itertuples(index=False)):
        frb = str(row.FRB)

        # Legacy IMR FITS
        legacy_dir = Path(args.legacy_imr_dir) / frb
        legacy_paths = {
            "image": legacy_dir / f"{frb}_image.fits",
            "model": legacy_dir / f"{frb}_model.fits",
            "resid": legacy_dir / f"{frb}_resid.fits",
        }

        legacy_data = {}
        for key in ["image", "model", "resid"]:
            p = legacy_paths[key]
            if p.exists():
                legacy_data[key] = load_fits_image(p)
            else:
                legacy_data[key] = None

        # GALFIT IMR from requested run subdirectory.
        galfit_out = Path(args.galfit_runs_dir) / frb / args.galfit_run_subdir / "out.fits"

        galfit_data = {"image": None, "model": None, "resid": None}
        if galfit_out.exists():
            try:
                galfit_data["image"] = load_galfit_hdu_image(galfit_out, 1)
                galfit_data["model"] = load_galfit_hdu_image(galfit_out, 2)
                galfit_data["resid"] = load_galfit_hdu_image(galfit_out, 3)
            except Exception:
                pass

        # Use GALFIT image limits as row reference, then apply the same arcsin stretch to all six panels.
        ref_img = galfit_data["image"]
        if ref_img is not None:
            disp_vmin = float(np.nanmin(ref_img))
            disp_vmax = float(np.nanmax(ref_img))
        else:
            row_arrays = [a for a in [
                legacy_data["model"], legacy_data["resid"],
                galfit_data["image"], galfit_data["model"], galfit_data["resid"],
            ] if a is not None]
            if row_arrays:
                combo = np.concatenate([a.ravel() for a in row_arrays])
                disp_vmin = float(np.nanmin(combo))
                disp_vmax = float(np.nanmax(combo))
            else:
                disp_vmin, disp_vmax = 0.0, 1.0
        if (not np.isfinite(disp_vmin)) or (not np.isfinite(disp_vmax)) or (disp_vmax <= disp_vmin):
            disp_vmin, disp_vmax = 0.0, 1.0

        for k, key in enumerate(["image", "model", "resid"]):
            ax = axes[i, k]
            arr = legacy_data[key]
            if arr is not None:
                arr_disp = apply_arcsin_stretch(arr, disp_vmin, disp_vmax)
                ax.imshow(arr_disp, origin="lower", cmap=panel_cmap, vmin=-np.pi / 2.0, vmax=np.pi / 2.0)
            else:
                ax.text(0.5, 0.5, "missing", ha="center", va="center")
            ax.set_xticks([])
            ax.set_yticks([])

        for k, key in enumerate(["image", "model", "resid"], start=3):
            ax = axes[i, k]
            arr = galfit_data[key]
            if arr is not None:
                arr_disp = apply_arcsin_stretch(arr, disp_vmin, disp_vmax)
                ax.imshow(arr_disp, origin="lower", cmap=panel_cmap, vmin=-np.pi / 2.0, vmax=np.pi / 2.0)
            else:
                ax.text(0.5, 0.5, "missing", ha="center", va="center")
            ax.set_xticks([])
            ax.set_yticks([])

        # High-contrast row label inside the first panel to avoid clipping/overlap.
        axes[i, 0].text(
            0.02,
            0.98,
            f"FRB {frb}",
            transform=axes[i, 0].transAxes,
            ha="left",
            va="top",
            fontsize=12,
            fontweight="bold",
            color="white",
            bbox={"facecolor": "black", "alpha": 0.65, "pad": 2, "edgecolor": "none"},
        )

        # Intentionally omit per-row left-side metadata text to prevent overlap.

    plt.tight_layout(rect=[0.0, 0.0, 1.0, 0.99], h_pad=1.0)
    out_png = out_dir / f"panel_top_delta_i_{args.tag}.png"
    out_pdf = out_dir / f"panel_top_delta_i_{args.tag}.pdf"
    plt.savefig(out_png, dpi=250, bbox_inches="tight")
    plt.savefig(out_pdf, dpi=250, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved {list_csv}")
    print(f"Saved {out_png}")
    print(f"Saved {out_pdf}")


if __name__ == "__main__":
    main()
