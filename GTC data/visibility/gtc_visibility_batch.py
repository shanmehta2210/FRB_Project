#!/usr/bin/env python3
"""Batch GTC rigorous visibility filter for FRB targets.

Applies five deterministic gates for inclination-safe imaging at GTC
(Roque de los Muchachos) and writes per-target CSV reports under
``GTC data/visibility/``.

Examples:
  python "GTC data/visibility/gtc_visibility_batch.py" --date 2026-06-24
  python "GTC data/visibility/gtc_visibility_batch.py" --date 2026-06-24 --end-date 2026-07-24
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timedelta
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import AltAz, SkyCoord, get_body
from astropy.time import Time
from astroplan import FixedTarget, Observer
from astroplan.constraints import (
    AirmassConstraint,
    AltitudeConstraint,
    AtNightConstraint,
    MoonSeparationConstraint,
)
from astroplan.moon import moon_illumination
from astroplan.utils import time_grid_from_range

VIS_DIR = Path(__file__).resolve().parent
GTC_DATA = VIS_DIR.parent
REPO = GTC_DATA.parent
DEFAULT_CSV = REPO / "master_frb_localization.csv"
NIGHTLY_DIR = VIS_DIR / "nightly"
SUMMARY_DIR = VIS_DIR / "summaries"

GTC_MIN_ELEVATION_DEG = 25.0
GTC_MAX_ELEVATION_DEG = 72.0
GTC_LAT_DEG = 28.0 + 45.0 / 60.0 + 24.0 / 3600.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Batch GTC visibility filter with rigorous science gates.",
    )
    parser.add_argument(
        "--date",
        required=True,
        help="First observing night (YYYY-MM-DD) at La Palma.",
    )
    parser.add_argument(
        "--end-date",
        default=None,
        help="Last night inclusive (YYYY-MM-DD). Writes nightly CSVs plus a summary rollup.",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=DEFAULT_CSV,
        help=f"Target catalog CSV (default: {DEFAULT_CSV.name}).",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output CSV for single-night run (default: visibility/nightly/gtc_visibility_DATE.csv).",
    )
    parser.add_argument(
        "--time-resolution-min",
        type=float,
        default=5.0,
        help="Time grid step in minutes (default: 5).",
    )
    parser.add_argument(
        "--min-duration-min",
        type=float,
        default=30.0,
        help="Minimum contiguous valid window in minutes (default: 30).",
    )
    parser.add_argument(
        "--max-airmass",
        type=float,
        default=1.5,
        help="Gate 2: maximum airmass (default: 1.5).",
    )
    parser.add_argument(
        "--min-moon-sep",
        type=float,
        default=30.0,
        help="Gate 3: minimum moon separation in degrees (default: 30).",
    )
    parser.add_argument(
        "--min-elevation",
        type=float,
        default=GTC_MIN_ELEVATION_DEG,
        help=f"Gate 1 mechanical minimum elevation (default: {GTC_MIN_ELEVATION_DEG}).",
    )
    parser.add_argument(
        "--max-elevation",
        type=float,
        default=GTC_MAX_ELEVATION_DEG,
        help=f"Gate 5 dome vignetting maximum elevation (default: {GTC_MAX_ELEVATION_DEG}).",
    )
    parser.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress per-night summaries when scanning a date range.",
    )
    return parser.parse_args()


def dec_threshold_for_min_elevation(lat_deg: float, h_min_deg: float) -> float:
    return h_min_deg - 90.0 + lat_deg


def iter_nights(start: str, end: str) -> list[str]:
    d0 = datetime.strptime(start, "%Y-%m-%d").date()
    d1 = datetime.strptime(end, "%Y-%m-%d").date()
    if d1 < d0:
        raise ValueError(f"--end-date {end} is before --date {start}")
    nights: list[str] = []
    d = d0
    while d <= d1:
        nights.append(d.strftime("%Y-%m-%d"))
        d += timedelta(days=1)
    return nights


def load_targets(csv_path: Path) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    required = {"frb", "ra_deg", "dec_deg"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"CSV missing columns: {sorted(missing)}")
    df = df.dropna(subset=["ra_deg", "dec_deg"]).copy()
    if "coord_semantics" not in df.columns:
        df["coord_semantics"] = ""
    return df


def _targets_from_df(df: pd.DataFrame) -> list[FixedTarget]:
    return [
        FixedTarget(
            coord=SkyCoord(ra=row["ra_deg"] * u.deg, dec=row["dec_deg"] * u.deg),
            name=str(row["frb"]),
        )
        for _, row in df.iterrows()
    ]


def _constraint_mask(
    constraints: list,
    observer: Observer,
    targets: list[FixedTarget],
    times: Time,
) -> np.ndarray:
    n_times = len(times)
    n_targets = len(targets)
    out = np.ones((n_times, n_targets), dtype=bool)
    for j, target in enumerate(targets):
        if not constraints:
            continue
        masks = [c.compute_constraint(times, observer, target) for c in constraints]
        out[:, j] = np.logical_and.reduce(masks)
    return out


def _longest_contiguous_run(mask_1d: np.ndarray, min_len: int) -> tuple[int, int] | None:
    if mask_1d.size == 0:
        return None
    best: tuple[int, int] | None = None
    best_len = 0
    start = 0
    for i, val in enumerate(mask_1d):
        if not val:
            start = i + 1
            continue
        run_len = i - start + 1
        if run_len >= min_len and run_len > best_len:
            best_len = run_len
            best = (start, i + 1)
    return best


def _format_utc(t: Time) -> str:
    return t.isot.replace("T", " ")[:19]


def _parallactic_angle_deg(observer: Observer, target: FixedTarget, t: Time) -> float:
    lst = observer.local_sidereal_time(t)
    ha = (lst - target.coord.ra).to(u.rad).value
    lat = observer.location.lat.rad
    dec = target.coord.dec.rad
    eta = np.arctan2(np.sin(ha), np.tan(lat) * np.cos(dec) - np.sin(dec) * np.cos(ha))
    return float(np.degrees(eta))


def _diagnostics_at_time(
    observer: Observer,
    target: FixedTarget,
    t: Time,
) -> dict[str, float]:
    altaz_frame = AltAz(obstime=t, location=observer.location)
    target_az = target.coord.transform_to(altaz_frame)
    moon = get_body("moon", t, observer.location)
    sep = target.coord.separation(moon).deg
    alt = target_az.alt.deg
    airmass = float(target_az.secz) if target_az.alt.deg > 0 else np.nan
    pa = _parallactic_angle_deg(observer, target, t)
    return {
        "min_airmass": float(airmass),
        "min_moon_sep_deg": float(sep),
        "max_target_alt_deg": float(alt),
        "parallactic_angle_deg": float(pa),
        "moon_illumination_frac": float(moon_illumination(t)),
    }


def evaluate_visibility(
    df: pd.DataFrame,
    obs_date: str,
    *,
    time_resolution_min: float,
    min_duration_min: float,
    max_airmass: float,
    min_moon_sep: float,
    min_elevation: float,
    max_elevation: float,
) -> pd.DataFrame:
    observer = Observer.at_site("Roque de los Muchachos")
    dec_min = dec_threshold_for_min_elevation(GTC_LAT_DEG, min_elevation)

    obs_noon = Time(f"{obs_date} 12:00:00")
    window_start = observer.twilight_evening_astronomical(obs_noon, which="next")
    window_end = observer.twilight_morning_astronomical(obs_noon, which="next")

    times = time_grid_from_range(
        [window_start, window_end],
        time_resolution=time_resolution_min * u.min,
    )
    min_samples = max(1, int(np.ceil(min_duration_min / time_resolution_min)))

    targets = _targets_from_df(df)
    night = AtNightConstraint.twilight_astronomical()

    mask_g2 = _constraint_mask(
        [night, AirmassConstraint(max=max_airmass)],
        observer,
        targets,
        times,
    )
    mask_g3 = _constraint_mask(
        [night, MoonSeparationConstraint(min=min_moon_sep * u.deg)],
        observer,
        targets,
        times,
    )
    mask_g4 = _constraint_mask([night], observer, targets, times)
    mask_g5 = _constraint_mask(
        [night, AltitudeConstraint(max=max_elevation * u.deg)],
        observer,
        targets,
        times,
    )
    mask_all = _constraint_mask(
        [
            night,
            AirmassConstraint(max=max_airmass),
            MoonSeparationConstraint(min=min_moon_sep * u.deg),
            AltitudeConstraint(max=max_elevation * u.deg),
        ],
        observer,
        targets,
        times,
    )

    gate1 = df["dec_deg"].to_numpy() >= dec_min
    gate2 = np.any(mask_g2, axis=0)
    gate3 = np.any(mask_g3, axis=0)
    gate4 = np.any(mask_g4, axis=0)
    gate5 = np.any(mask_g5, axis=0)

    mid_jd = (window_start.jd + window_end.jd) / 2.0
    moon_illum_at_mid = float(moon_illumination(Time(mid_jd, format="jd")))

    rows: list[dict] = []
    for i, (_, src) in enumerate(df.iterrows()):
        combined = mask_all[:, i]
        run = _longest_contiguous_run(combined, min_samples)
        duration_ok = run is not None

        rigorous = bool(
            gate1[i] and gate2[i] and gate3[i] and gate4[i] and gate5[i] and duration_ok
        )

        best_start = ""
        best_end = ""
        diag = {
            "min_airmass": np.nan,
            "min_moon_sep_deg": np.nan,
            "max_target_alt_deg": np.nan,
            "parallactic_angle_deg": np.nan,
            "moon_illumination_frac": moon_illum_at_mid,
        }

        if run is not None:
            s, e = run
            best_start = _format_utc(times[s])
            best_end = _format_utc(times[e - 1])
            mid_idx = s + (e - s) // 2
            diag = _diagnostics_at_time(observer, targets[i], times[mid_idx])

        rows.append(
            {
                "obs_date": obs_date,
                "frb": src["frb"],
                "ra_deg": round(float(src["ra_deg"]), 6),
                "dec_deg": round(float(src["dec_deg"]), 6),
                "coord_semantics": src.get("coord_semantics", ""),
                "gate1_mechanical_pass": bool(gate1[i]),
                "gate2_airmass_pass": bool(gate2[i]),
                "gate3_moon_pass": bool(gate3[i]),
                "gate4_dark_pass": bool(gate4[i]),
                "gate5_dome_pass": bool(gate5[i]),
                "duration_pass": duration_ok,
                "rigorous_science_pass": rigorous,
                "obs_window_start_utc": _format_utc(window_start),
                "obs_window_end_utc": _format_utc(window_end),
                "best_start_utc": best_start,
                "best_end_utc": best_end,
                "min_airmass": round(diag["min_airmass"], 3)
                if np.isfinite(diag["min_airmass"])
                else "",
                "min_moon_sep_deg": round(diag["min_moon_sep_deg"], 2)
                if np.isfinite(diag["min_moon_sep_deg"])
                else "",
                "max_target_alt_deg": round(diag["max_target_alt_deg"], 2)
                if np.isfinite(diag["max_target_alt_deg"])
                else "",
                "moon_illumination_frac": round(diag["moon_illumination_frac"], 3),
                "parallactic_angle_deg": round(diag["parallactic_angle_deg"], 2)
                if np.isfinite(diag["parallactic_angle_deg"])
                else "",
            }
        )

    return pd.DataFrame(rows)


def print_summary(out_df: pd.DataFrame, *, label: str = "") -> None:
    n = len(out_df)
    header = f"Targets evaluated: {n}"
    if label:
        header = f"{label} — {header}"
    print(header)
    for col, gate_label in [
        ("gate1_mechanical_pass", "Gate 1 (mechanical)"),
        ("gate2_airmass_pass", "Gate 2 (airmass)"),
        ("gate3_moon_pass", "Gate 3 (moon)"),
        ("gate4_dark_pass", "Gate 4 (dark)"),
        ("gate5_dome_pass", "Gate 5 (dome)"),
        ("duration_pass", "Duration"),
        ("rigorous_science_pass", "All gates"),
    ]:
        count = int(out_df[col].sum())
        print(f"  {gate_label}: {count}/{n}")


def _eval_kwargs(args: argparse.Namespace) -> dict:
    return {
        "time_resolution_min": args.time_resolution_min,
        "min_duration_min": args.min_duration_min,
        "max_airmass": args.max_airmass,
        "min_moon_sep": args.min_moon_sep,
        "min_elevation": args.min_elevation,
        "max_elevation": args.max_elevation,
    }


def build_range_summary(long_df: pd.DataFrame, n_nights: int) -> pd.DataFrame:
    passed = long_df[long_df["rigorous_science_pass"]].copy()
    grouped = passed.groupby("frb", as_index=False).agg(
        ra_deg=("ra_deg", "first"),
        dec_deg=("dec_deg", "first"),
        coord_semantics=("coord_semantics", "first"),
        n_pass_nights=("obs_date", "count"),
        pass_dates=("obs_date", lambda s: ";".join(sorted(s))),
        first_pass_date=("obs_date", "min"),
        last_pass_date=("obs_date", "max"),
    )
    all_frbs = long_df[["frb", "ra_deg", "dec_deg", "coord_semantics"]].drop_duplicates("frb")
    summary = all_frbs.merge(grouped, on=["frb", "ra_deg", "dec_deg", "coord_semantics"], how="left")
    summary["n_pass_nights"] = summary["n_pass_nights"].fillna(0).astype(int)
    summary["pass_dates"] = summary["pass_dates"].fillna("")
    summary["first_pass_date"] = summary["first_pass_date"].fillna("")
    summary["last_pass_date"] = summary["last_pass_date"].fillna("")
    summary["n_nights_scanned"] = n_nights
    summary["pass_fraction"] = (summary["n_pass_nights"] / n_nights).round(3)
    summary = summary.sort_values(["n_pass_nights", "frb"], ascending=[False, True])
    return summary


def build_nightly_counts(long_df: pd.DataFrame) -> pd.DataFrame:
    counts = (
        long_df.groupby("obs_date", as_index=False)
        .agg(
            n_pass=("rigorous_science_pass", "sum"),
            n_targets=("frb", "count"),
            moon_illumination_frac=("moon_illumination_frac", "first"),
        )
        .sort_values("obs_date")
    )
    return counts


def run_date_range(df: pd.DataFrame, args: argparse.Namespace) -> int:
    nights = iter_nights(args.date, args.end_date)
    NIGHTLY_DIR.mkdir(parents=True, exist_ok=True)
    SUMMARY_DIR.mkdir(parents=True, exist_ok=True)

    nightly_frames: list[pd.DataFrame] = []
    eval_kw = _eval_kwargs(args)

    for i, night in enumerate(nights, 1):
        out_df = evaluate_visibility(df, night, **eval_kw)
        nightly_path = NIGHTLY_DIR / f"gtc_visibility_{night}.csv"
        out_df.to_csv(nightly_path, index=False)
        nightly_frames.append(out_df)
        if not args.quiet:
            print(f"[{i}/{len(nights)}] Wrote {nightly_path}")
            print_summary(out_df, label=night)

    long_df = pd.concat(nightly_frames, ignore_index=True)
    long_path = SUMMARY_DIR / f"gtc_visibility_long_{args.date}_to_{args.end_date}.csv"
    long_df.to_csv(long_path, index=False)

    frb_summary = build_range_summary(long_df, len(nights))
    frb_path = SUMMARY_DIR / f"gtc_availability_by_frb_{args.date}_to_{args.end_date}.csv"
    frb_summary.to_csv(frb_path, index=False)

    night_summary = build_nightly_counts(long_df)
    night_path = SUMMARY_DIR / f"gtc_availability_by_night_{args.date}_to_{args.end_date}.csv"
    night_summary.to_csv(night_path, index=False)

    print(f"\nScanned {len(nights)} nights ({args.date} → {args.end_date})")
    print(f"  Long format: {long_path}")
    print(f"  Per-FRB rollup: {frb_path}")
    print(f"  Per-night counts: {night_path}")
    print(f"  FRBs passing ≥1 night: {(frb_summary['n_pass_nights'] > 0).sum()}/{len(frb_summary)}")
    print(f"  FRBs passing all nights: {(frb_summary['n_pass_nights'] == len(nights)).sum()}")
    return 0


def main() -> int:
    args = parse_args()
    csv_path = args.csv if args.csv.is_absolute() else REPO / args.csv
    df = load_targets(csv_path)
    eval_kw = _eval_kwargs(args)

    if args.end_date:
        return run_date_range(df, args)

    out_df = evaluate_visibility(df, args.date, **eval_kw)

    if args.out is not None:
        out_path = args.out if args.out.is_absolute() else REPO / args.out
    else:
        NIGHTLY_DIR.mkdir(parents=True, exist_ok=True)
        out_path = NIGHTLY_DIR / f"gtc_visibility_{args.date}.csv"

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, index=False)
    print(f"Wrote {out_path}")
    print_summary(out_df)
    return 0


if __name__ == "__main__":
    sys.exit(main())
