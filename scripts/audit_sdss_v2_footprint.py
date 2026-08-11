#!/usr/bin/env python3
"""
Audit SDSS v2 null catalog: footprint, deduplication, strict mag20 pool.

Writes outputs under plots/plots_null/v2/sdss_audit/.

Run from repo root::

    python scripts/audit_sdss_v2_footprint.py
    python scripts/audit_sdss_v2_footprint.py --compare-v1
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (
    JOINT_DEC_MAX,
    JOINT_DEC_MIN,
    MAG_LIMIT,
    Q0,
    SDSS_MAG20_LIMIT,
    SDSS_MIN_STRICT_MAG20_POOL,
    SDSS_UR_MAX_CDF,
    count_strict_mag20_pool,
    cut_funnel_rows,
    footprint_summary,
)
from pipeline_null_plot_utils import DEFAULT_SDSS, DEFAULT_SDSS_V2, REPO_ROOT

OUT_ROOT = REPO_ROOT / "plots" / "plots_null" / "v2" / "sdss_audit"


def _load_batch_log(v2_csv: Path) -> pd.DataFrame | None:
    log_path = v2_csv.with_suffix(".batch_log.csv")
    if log_path.is_file():
        return pd.read_csv(log_path)
    return None


def plot_ra_dec(
    v2: pd.DataFrame,
    v1: pd.DataFrame | None,
    out_stem: Path,
    *,
    max_points: int = 80_000,
) -> None:
    fig, axes = plt.subplots(1, 2 if v1 is not None else 1, figsize=(12 if v1 is not None else 7, 5))
    if v1 is None:
        axes = [axes]

    for ax, df, label, color in zip(
        axes,
        [v2, v1] if v1 is not None else [v2],
        ["SDSS v2", "SDSS v1"] if v1 is not None else ["SDSS v2"],
        ["#4daf4a", "#377eb8"] if v1 is not None else ["#4daf4a"],
    ):
        sub = df
        if len(sub) > max_points:
            sub = sub.sample(n=max_points, random_state=42)
        ra = np.radians(pd.to_numeric(sub["RA_ICRS"], errors="coerce"))
        dec = np.radians(pd.to_numeric(sub["DE_ICRS"], errors="coerce"))
        ok = np.isfinite(ra) & np.isfinite(dec)
        ax.scatter(ra[ok], dec[ok], s=0.3, c=color, alpha=0.15, rasterized=True)
        ax.set_xlabel("RA (rad)")
        ax.set_ylabel("Dec (rad)")
        ax.set_title(f"{label} (N={len(df):,}, subsample shown)")
        ax.axhline(np.radians(JOINT_DEC_MIN), color="0.5", ls="--", lw=0.8, label="joint Dec clip")
        ax.axhline(np.radians(JOINT_DEC_MAX), color="0.5", ls="--", lw=0.8)
        ax.legend(fontsize=7)

    fig.suptitle("SDSS null catalog sky distribution", fontsize=12)
    plt.tight_layout()
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_stem.with_suffix(".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_stem.with_suffix('.png')}")


def plot_dec_histogram(v2: pd.DataFrame, v1: pd.DataFrame | None, out_stem: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4))
    dec_v2 = pd.to_numeric(v2["DE_ICRS"], errors="coerce")
    ax.hist(dec_v2, bins=60, density=True, alpha=0.55, color="#4daf4a", label=f"v2 (N={len(v2):,})")
    if v1 is not None:
        dec_v1 = pd.to_numeric(v1["DE_ICRS"], errors="coerce")
        ax.hist(dec_v1, bins=60, density=True, alpha=0.45, color="#377eb8", label=f"v1 (N={len(v1):,})")
    ax.axvline(JOINT_DEC_MIN, color="0.4", ls="--", label="joint footprint edges")
    ax.axvline(JOINT_DEC_MAX, color="0.4", ls="--")
    ax.set_xlabel("Dec (deg)")
    ax.set_ylabel("Normalized density")
    ax.set_title("Declination distribution")
    ax.legend(fontsize=8)
    plt.tight_layout()
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_stem.with_suffix(".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_stem.with_suffix('.png')}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--v2-csv", type=Path, default=DEFAULT_SDSS_V2)
    parser.add_argument("--v1-csv", type=Path, default=DEFAULT_SDSS)
    parser.add_argument("--compare-v1", action="store_true")
    parser.add_argument("--out-root", type=Path, default=OUT_ROOT)
    args = parser.parse_args()

    if not args.v2_csv.is_file():
        raise SystemExit(f"Missing v2 catalog: {args.v2_csv}")

    args.out_root.mkdir(parents=True, exist_ok=True)
    v2 = pd.read_csv(args.v2_csv)
    v1 = pd.read_csv(args.v1_csv) if args.compare_v1 and args.v1_csv.is_file() else None

    footprint_summary(v2, "SDSS v2 audit")
    pool_mag20 = count_strict_mag20_pool(v2)
    pool_mag21 = count_strict_mag20_pool(v2, mag_limit=MAG_LIMIT)

    summary_rows = [
        {"metric": "n_rows", "value": len(v2)},
        {"metric": "strict_mag20_pool", "value": pool_mag20},
        {"metric": "min_strict_mag20_required", "value": SDSS_MIN_STRICT_MAG20_POOL},
        {
            "metric": "lnL_coverage_frac",
            "value": float(pd.to_numeric(v2.get("lnLExp_r"), errors="coerce").notna().mean()),
        },
        {"metric": "n_duplicate_objID", "value": int(v2["objID"].duplicated().sum()) if "objID" in v2.columns else -1},
        {"metric": "dec_min", "value": float(pd.to_numeric(v2["DE_ICRS"], errors="coerce").min())},
        {"metric": "dec_max", "value": float(pd.to_numeric(v2["DE_ICRS"], errors="coerce").max())},
        {"metric": "ra_min", "value": float(pd.to_numeric(v2["RA_ICRS"], errors="coerce").min())},
        {"metric": "ra_max", "value": float(pd.to_numeric(v2["RA_ICRS"], errors="coerce").max())},
    ]
    summary_df = pd.DataFrame(summary_rows)
    summary_path = args.out_root / "audit_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"Wrote {summary_path}")

    funnel = cut_funnel_rows(
        v2,
        survey="sdss_v2",
        mag_limit=SDSS_MAG20_LIMIT,
        sample_mode="strict",
        mag_column="modelMag_r",
        q_column="expAB_r",
        q0=Q0,
        is_legacy=False,
        sdss_ur_max=SDSS_UR_MAX_CDF,
        sdss_exp_winner_only=True,
    )
    funnel_mag21 = cut_funnel_rows(
        v2,
        survey="sdss_v2",
        mag_limit=MAG_LIMIT,
        sample_mode="strict",
        mag_column="modelMag_r",
        q_column="expAB_r",
        q0=Q0,
        is_legacy=False,
        sdss_ur_max=SDSS_UR_MAX_CDF,
        sdss_exp_winner_only=True,
    )
    funnel_df = pd.concat([pd.DataFrame(funnel), pd.DataFrame(funnel_mag21)], ignore_index=True)
    funnel_path = args.out_root / "cut_funnel_v2.csv"
    funnel_df.to_csv(funnel_path, index=False)
    print(f"Wrote {funnel_path}")
    print(f"[*] strict mag20 pool N={pool_mag20}, mag21 N={pool_mag21}")

    batch_log = _load_batch_log(args.v2_csv)
    if batch_log is not None:
        batch_path = args.out_root / "batch_log_copy.csv"
        batch_log.to_csv(batch_path, index=False)
        print(f"Wrote {batch_path}")

    plot_ra_dec(v2, v1, args.out_root / "sky_ra_dec")
    plot_dec_histogram(v2, v1, args.out_root / "dec_histogram")

    readme = args.out_root / "README.md"
    readme.write_text(
        f"""# SDSS v2 catalog audit

Generated by `scripts/audit_sdss_v2_footprint.py`.

| Metric | Value |
|--------|-------|
| Rows | {len(v2):,} |
| Strict mag20 pool | {pool_mag20:,} |
| Required minimum | {SDSS_MIN_STRICT_MAG20_POOL:,} |
| lnL coverage | {summary_rows[3]['value']:.1%} |
| Duplicate objID | {summary_rows[4]['value']} |
| Dec range | {summary_rows[5]['value']:.2f} … {summary_rows[6]['value']:.2f} |

See `cut_funnel_v2.csv` for per-stage counts at mag&lt;20 and mag&lt;21.
""",
        encoding="utf-8",
    )
    print(f"Wrote {readme}")

    if pool_mag20 < SDSS_MIN_STRICT_MAG20_POOL:
        raise SystemExit(
            f"[FAIL] strict mag20 pool {pool_mag20} < {SDSS_MIN_STRICT_MAG20_POOL}"
        )
    if "objID" in v2.columns and v2["objID"].duplicated().any():
        raise SystemExit("[FAIL] duplicate objID in v2 catalog")


if __name__ == "__main__":
    main()
