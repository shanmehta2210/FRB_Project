from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.modeling import fitting, models
from astropy.nddata import block_reduce


ANGLES = [0, 45, 90, 135]


def safe_ratio(num: float, den: float) -> float:
    if not np.isfinite(num) or not np.isfinite(den) or den == 0.0:
        return np.nan
    return float(num / den)


def safe_delta(a: float, b: float) -> float:
    if not np.isfinite(a) or not np.isfinite(b):
        return np.nan
    return float(a - b)


def better_is_lower(a: float, b: float) -> str:
    if not np.isfinite(a) and not np.isfinite(b):
        return "unknown"
    if not np.isfinite(a):
        return "PSFEx"
    if not np.isfinite(b):
        return "ePSF"
    if a < b:
        return "ePSF"
    if b < a:
        return "PSFEx"
    return "tie"


def compute_fwhm(profile: np.ndarray, r_vals: np.ndarray) -> float:
    if profile.size == 0 or r_vals.size == 0:
        return np.nan

    peak = float(np.max(profile))
    if not np.isfinite(peak) or peak <= 0.0:
        return np.nan

    half_max = peak / 2.0
    idx = np.where(profile < half_max)[0]
    if len(idx) == 0:
        return np.nan

    cross_idx = int(idx[0])
    if cross_idx > 0:
        r1, r2 = float(r_vals[cross_idx - 1]), float(r_vals[cross_idx])
        p1, p2 = float(profile[cross_idx - 1]), float(profile[cross_idx])
        if p2 == p1:
            return np.nan
        r_hm = r1 + (half_max - p1) * (r2 - r1) / (p2 - p1)
    else:
        r_hm = float(r_vals[0])

    return 2.0 * r_hm


def extract_radial_profile(image: np.ndarray, center: tuple[int, int], angle_deg: float, max_radius: int):
    cx, cy = center
    r_vals = np.arange(0, max_radius, 1.0)

    angle_rad = np.deg2rad(angle_deg)
    x_vals = cx + r_vals * np.cos(angle_rad)
    y_vals = cy + r_vals * np.sin(angle_rad)

    x_idx = np.round(x_vals).astype(int)
    y_idx = np.round(y_vals).astype(int)

    valid = (
        (x_idx >= 0)
        & (x_idx < image.shape[1])
        & (y_idx >= 0)
        & (y_idx < image.shape[0])
    )

    return r_vals[valid], image[y_idx[valid], x_idx[valid]]


def fwhm_four_angles(image: np.ndarray) -> tuple[list[float], float]:
    cy, cx = np.unravel_index(np.argmax(image), image.shape)
    center = (int(cx), int(cy))

    max_r = min(cx, cy, image.shape[1] - cx, image.shape[0] - cy)
    if max_r < 2:
        max_r = min(image.shape[0], image.shape[1]) // 2

    vals = []
    for ang in ANGLES:
        r_vals, profile = extract_radial_profile(image, center, ang, max_r)
        vals.append(compute_fwhm(profile, r_vals))

    return vals, float(np.nanmean(vals))


def fit_moffat(image: np.ndarray) -> dict:
    y_shape, x_shape = image.shape
    y, x = np.mgrid[:y_shape, :x_shape]

    center_y, center_x = y_shape // 2, x_shape // 2
    max_amp = float(np.nanmax(image))
    if not np.isfinite(max_amp) or max_amp <= 0.0:
        return {
            "alpha": np.nan,
            "gamma": np.nan,
            "fwhm": np.nan,
            "max_frac_resid_pct": np.nan,
        }

    init = models.Moffat2D(
        amplitude=max_amp,
        x_0=float(center_x),
        y_0=float(center_y),
        gamma=2.0,
        alpha=3.0,
    )

    fitter = fitting.LevMarLSQFitter()
    try:
        fit = fitter(init, x, y, image)
        model_data = fit(x, y)
        resid = image - model_data
        max_frac_resid = float(np.nanmax(np.abs(resid)) / max_amp)
        alpha = float(fit.alpha.value)
        gamma = float(fit.gamma.value)
        if alpha <= 0.0 or gamma <= 0.0:
            fwhm = np.nan
        else:
            fwhm = float(2.0 * gamma * np.sqrt(2.0 ** (1.0 / alpha) - 1.0))
    except Exception:
        alpha = np.nan
        gamma = np.nan
        fwhm = np.nan
        max_frac_resid = np.nan

    return {
        "alpha": alpha,
        "gamma": gamma,
        "fwhm": fwhm,
        "max_frac_resid_pct": max_frac_resid * 100.0 if np.isfinite(max_frac_resid) else np.nan,
    }


def read_fits_2d(path: Path) -> np.ndarray:
    with fits.open(path) as hdul:
        for hdu in hdul:
            if isinstance(hdu.data, np.ndarray) and hdu.data.ndim == 2:
                arr = np.asarray(hdu.data, dtype=np.float64)
                return arr
    raise ValueError(f"No 2D image found in {path}")


def main():
    repo_root = Path(__file__).resolve().parents[1]
    summary_path = repo_root / "psf_fwhm_summary.csv"
    compact_path = repo_root / "psf_constant_comparison.csv"
    epsf_root = repo_root / "psfs"
    psfex_root = repo_root / "psfs" / "PSFEx + SExtractor" / "final_center_psfs"

    df = pd.read_csv(summary_path)
    if "FRB" not in df.columns:
        raise ValueError("psf_fwhm_summary.csv must contain FRB column")

    generated_cols = [
        "ePSF_moffat_alpha",
        "ePSF_moffat_gamma",
        "ePSF_moffat_fwhm",
        "ePSF_max_frac_resid_pct",
        "PSFEx_FWHM_0",
        "PSFEx_FWHM_45",
        "PSFEx_FWHM_90",
        "PSFEx_FWHM_135",
        "PSFEx_Avg_FWHM",
        "PSFEx_moffat_alpha",
        "PSFEx_moffat_gamma",
        "PSFEx_moffat_fwhm",
        "PSFEx_max_frac_resid_pct",
        "Delta_Avg_FWHM_PSFEx_minus_ePSF",
        "Delta_moffat_fwhm_PSFEx_minus_ePSF",
        "ePSF_exp_minus_moffat_fwhm",
        "PSFEx_exp_minus_moffat_fwhm",
        "ePSF_exp_over_moffat_fwhm_ratio",
        "PSFEx_exp_over_moffat_fwhm_ratio",
        "Abs_ePSF_exp_minus_moffat_fwhm",
        "Abs_PSFEx_exp_minus_moffat_fwhm",
        "Delta_fit_resid_pct_PSFEx_minus_ePSF",
        "Better_Moffat_Fit_by_MaxFracResid",
        "Better_ExpVsMoffat_Agreement",
    ]

    # Keep this updater idempotent: remove previously generated columns and
    # stale merge artifacts like *_x/*_y before adding fresh values.
    drop_cols = []
    for col in df.columns:
        if col in generated_cols:
            drop_cols.append(col)
            continue
        if col.endswith("_x") or col.endswith("_y"):
            base = col[:-2]
            if base in generated_cols:
                drop_cols.append(col)
    if drop_cols:
        df = df.drop(columns=drop_cols)

    records = []
    for frb in df["FRB"].astype(str):
        epsf_path = epsf_root / f"{frb}_flux_psf.fits"
        psfex_path = psfex_root / frb / f"{frb}_center_psf_25.fits"

        row = {"FRB": frb}

        epsf_moffat_fit = {"fwhm": np.nan}
        if epsf_path.exists():
            epsf_img = read_fits_2d(epsf_path)
            # Keep historical convention: ePSF products are oversampled by 4x.
            if epsf_img.shape[0] >= 4 and epsf_img.shape[1] >= 4:
                epsf_img_native = block_reduce(epsf_img, 4, func=np.mean)
            else:
                epsf_img_native = epsf_img

            epsf_moffat_fit = fit_moffat(epsf_img_native)
            row["ePSF_moffat_alpha"] = epsf_moffat_fit["alpha"]
            row["ePSF_moffat_gamma"] = epsf_moffat_fit["gamma"]
            row["ePSF_moffat_fwhm"] = epsf_moffat_fit["fwhm"]
            row["ePSF_max_frac_resid_pct"] = epsf_moffat_fit["max_frac_resid_pct"]
        else:
            row["ePSF_moffat_alpha"] = np.nan
            row["ePSF_moffat_gamma"] = np.nan
            row["ePSF_moffat_fwhm"] = np.nan
            row["ePSF_max_frac_resid_pct"] = np.nan

        if psfex_path.exists():
            psfex_img = read_fits_2d(psfex_path)
            fwhms, avg_fwhm = fwhm_four_angles(psfex_img)
            row["PSFEx_FWHM_0"] = fwhms[0]
            row["PSFEx_FWHM_45"] = fwhms[1]
            row["PSFEx_FWHM_90"] = fwhms[2]
            row["PSFEx_FWHM_135"] = fwhms[3]
            row["PSFEx_Avg_FWHM"] = avg_fwhm

            mfit = fit_moffat(psfex_img)
            row["PSFEx_moffat_alpha"] = mfit["alpha"]
            row["PSFEx_moffat_gamma"] = mfit["gamma"]
            row["PSFEx_moffat_fwhm"] = mfit["fwhm"]
            row["PSFEx_max_frac_resid_pct"] = mfit["max_frac_resid_pct"]
        else:
            row["PSFEx_FWHM_0"] = np.nan
            row["PSFEx_FWHM_45"] = np.nan
            row["PSFEx_FWHM_90"] = np.nan
            row["PSFEx_FWHM_135"] = np.nan
            row["PSFEx_Avg_FWHM"] = np.nan
            row["PSFEx_moffat_alpha"] = np.nan
            row["PSFEx_moffat_gamma"] = np.nan
            row["PSFEx_moffat_fwhm"] = np.nan
            row["PSFEx_max_frac_resid_pct"] = np.nan

        row["Delta_Avg_FWHM_PSFEx_minus_ePSF"] = (
            row["PSFEx_Avg_FWHM"] - float(df.loc[df["FRB"] == frb, "Avg_FWHM"].iloc[0])
            if pd.notna(row["PSFEx_Avg_FWHM"])
            else np.nan
        )

        row["Delta_moffat_fwhm_PSFEx_minus_ePSF"] = (
            row["PSFEx_moffat_fwhm"] - epsf_moffat_fit["fwhm"]
            if np.isfinite(row["PSFEx_moffat_fwhm"]) and np.isfinite(epsf_moffat_fit["fwhm"])
            else np.nan
        )

        # Constant/consistent comparison: experimental-vs-theoretical agreement per method.
        epsf_avg = float(df.loc[df["FRB"] == frb, "Avg_FWHM"].iloc[0])
        psfex_avg = row["PSFEx_Avg_FWHM"]

        row["ePSF_exp_minus_moffat_fwhm"] = safe_delta(epsf_avg, row["ePSF_moffat_fwhm"])
        row["PSFEx_exp_minus_moffat_fwhm"] = safe_delta(psfex_avg, row["PSFEx_moffat_fwhm"])

        row["ePSF_exp_over_moffat_fwhm_ratio"] = safe_ratio(epsf_avg, row["ePSF_moffat_fwhm"])
        row["PSFEx_exp_over_moffat_fwhm_ratio"] = safe_ratio(psfex_avg, row["PSFEx_moffat_fwhm"])

        row["Abs_ePSF_exp_minus_moffat_fwhm"] = (
            abs(row["ePSF_exp_minus_moffat_fwhm"])
            if np.isfinite(row["ePSF_exp_minus_moffat_fwhm"])
            else np.nan
        )
        row["Abs_PSFEx_exp_minus_moffat_fwhm"] = (
            abs(row["PSFEx_exp_minus_moffat_fwhm"])
            if np.isfinite(row["PSFEx_exp_minus_moffat_fwhm"])
            else np.nan
        )

        row["Delta_fit_resid_pct_PSFEx_minus_ePSF"] = safe_delta(
            row["PSFEx_max_frac_resid_pct"], row["ePSF_max_frac_resid_pct"]
        )
        row["Better_Moffat_Fit_by_MaxFracResid"] = better_is_lower(
            row["ePSF_max_frac_resid_pct"], row["PSFEx_max_frac_resid_pct"]
        )

        row["Better_ExpVsMoffat_Agreement"] = better_is_lower(
            row["Abs_ePSF_exp_minus_moffat_fwhm"], row["Abs_PSFEx_exp_minus_moffat_fwhm"]
        )

        records.append(row)

    cmp_df = pd.DataFrame(records)
    merged = pd.merge(df, cmp_df, on="FRB", how="left")

    # Preserve legacy columns while adding explicit ePSF aliases for clarity.
    for legacy, explicit in [
        ("moffat_alpha", "ePSF_moffat_alpha"),
        ("moffat_gamma", "ePSF_moffat_gamma"),
        ("moffat_fwhm", "ePSF_moffat_fwhm"),
        ("max_frac_resid_pct", "ePSF_max_frac_resid_pct"),
    ]:
        if legacy in merged.columns and explicit in merged.columns:
            merged[explicit] = merged[explicit].where(merged[explicit].notna(), merged[legacy])

    merged.to_csv(summary_path, index=False)

    compact_cols = [
        "FRB",
        "Avg_FWHM",
        "PSFEx_Avg_FWHM",
        "Delta_Avg_FWHM_PSFEx_minus_ePSF",
        "ePSF_moffat_fwhm",
        "PSFEx_moffat_fwhm",
        "Delta_moffat_fwhm_PSFEx_minus_ePSF",
        "ePSF_exp_minus_moffat_fwhm",
        "PSFEx_exp_minus_moffat_fwhm",
        "ePSF_exp_over_moffat_fwhm_ratio",
        "PSFEx_exp_over_moffat_fwhm_ratio",
        "Abs_ePSF_exp_minus_moffat_fwhm",
        "Abs_PSFEx_exp_minus_moffat_fwhm",
        "ePSF_max_frac_resid_pct",
        "PSFEx_max_frac_resid_pct",
        "Delta_fit_resid_pct_PSFEx_minus_ePSF",
        "Better_Moffat_Fit_by_MaxFracResid",
        "Better_ExpVsMoffat_Agreement",
    ]
    compact_existing = [c for c in compact_cols if c in merged.columns]
    compact = merged[compact_existing].copy()
    compact.to_csv(compact_path, index=False)

    with np.errstate(invalid="ignore"):
        avg_delta = float(np.nanmean(merged["Delta_Avg_FWHM_PSFEx_minus_ePSF"]))
        moffat_delta = float(np.nanmean(merged["Delta_moffat_fwhm_PSFEx_minus_ePSF"]))
    print(f"Updated {summary_path}")
    print(f"Wrote compact comparison {compact_path}")
    print(f"Rows: {len(merged)}")
    print(f"Mean Delta_Avg_FWHM_PSFEx_minus_ePSF: {avg_delta:.4f}")
    print(f"Mean Delta_moffat_fwhm_PSFEx_minus_ePSF: {moffat_delta:.4f}")


if __name__ == "__main__":
    main()
