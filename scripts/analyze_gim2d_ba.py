#!/usr/bin/env python3
"""GIM2D b/a diagnostics: mag scatter, cos(i) vs isotropic, morphology breakdown."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPTS))

from cosmos_ba_audit_utils import COSMOS_PLOTS, HST_DISK_ENTIRE_CSV, SDSS_DISK_ENTIRE_CSV  # noqa: E402
from null_catalog_utils import Q0, hubble_cosi_from_ba  # noqa: E402

SENTINEL = -999999.0


def cosi_from_ba(ba: np.ndarray, q0: float = Q0) -> np.ndarray:
    return np.array([hubble_cosi_from_ba(float(v), q0=q0) for v in ba], dtype=float)


def isotropic_ba(n: int, q0: float = Q0, seed: int = 42) -> tuple[np.ndarray, np.ndarray]:
    """Random orientations with P(i) proportional to sin(i), i in [0, pi/2]."""
    rng = np.random.default_rng(seed)
    cos_i = 1.0 - rng.uniform(0.0, 1.0, n)
    ba = np.sqrt(q0**2 + cos_i**2 * (1.0 - q0**2))
    return ba, cos_i


def print_pool_stats(label: str, ba: np.ndarray, cosi: np.ndarray, mag: np.ndarray | None = None) -> None:
    print(f"--- {label} ---")
    print(f"  N={len(ba):,}")
    print(f"  median b/a={np.median(ba):.3f}, mean b/a={np.mean(ba):.3f}")
    print(f"  median cos(i)={np.median(cosi):.3f}, mean cos(i)={np.mean(cosi):.3f}")
    print(f"  frac cos(i)>0.8 (face-on): {(cosi > 0.8).mean() * 100:.1f}%")
    print(f"  frac cos(i)<0.3 (edge-on): {(cosi < 0.3).mean() * 100:.1f}%")
    if mag is not None:
        for cut in (20, 21, 22):
            sel = mag <= cut
            if sel.sum() >= 100:
                print(
                    f"  mag<={cut}: median cos(i)={np.median(cosi[sel]):.3f}, "
                    f"median b/a={np.median(ba[sel]):.3f}, N={int(sel.sum()):,}"
                )


def plot_mag_ba_scatter(
    mag: np.ndarray,
    ba: np.ndarray,
    cosi: np.ndarray,
    *,
    out_png: Path,
    title: str,
) -> None:
    rng = np.random.default_rng(42)
    idx = rng.choice(len(mag), size=min(15000, len(mag)), replace=False)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    ax = axes[0]
    ax.scatter(mag[idx], ba[idx], s=3, alpha=0.15, c="0.35", rasterized=True)
    bins = np.arange(15, 26.5, 0.5)
    xc, med = [], []
    for lo in bins:
        mask = (mag >= lo) & (mag < lo + 0.5)
        if mask.sum() >= 30:
            med.append(float(np.median(ba[mask])))
            xc.append(lo + 0.25)
    ax.plot(xc, med, "r-o", lw=2, ms=4, label="median b/a")
    ax.set_xlabel("ACS_MAG_AUTO")
    ax.set_ylabel("GIM2D b/a = 1 - ELL_GIM2D")
    ax.set_title(title)
    ax.set_ylim(0, 1.02)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    ax = axes[1]
    ax.scatter(mag[idx], cosi[idx], s=3, alpha=0.15, c="0.35", rasterized=True)
    medc = []
    for lo in bins:
        mask = (mag >= lo) & (mag < lo + 0.5)
        if mask.sum() >= 30:
            medc.append(float(np.median(cosi[mask])))
    ax.plot(xc[: len(medc)], medc, "b-o", lw=2, ms=4, label="median cos(i)")
    ax.axhline(0.5, color="k", ls="--", lw=1, label="isotropic median = 0.5")
    ax.set_xlabel("ACS_MAG_AUTO")
    ax.set_ylabel(f"cos(i) Hubble q0={Q0}")
    ax.set_title("Inclination proxy vs magnitude")
    ax.set_ylim(0, 1.02)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[*] Wrote {out_png}")


def plot_cosi_histogram(cosi: np.ndarray, cosi_iso: np.ndarray, out_png: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(cosi, bins=50, density=True, alpha=0.65, label="GIM2D", color="#e41a1c")
    ax.hist(cosi_iso, bins=50, density=True, alpha=0.35, label="Isotropic sim", color="0.45")
    ax.axvline(float(np.median(cosi)), color="#e41a1c", lw=2, label=f"GIM2D median={np.median(cosi):.2f}")
    ax.axvline(0.5, color="k", ls="--", lw=1.5, label="isotropic median=0.5")
    ax.set_xlabel("cos(i)")
    ax.set_ylabel("density")
    ax.set_title("GIM2D disk cos(i) vs isotropic random orientations")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[*] Wrote {out_png}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hst-csv", type=Path, default=HST_DISK_ENTIRE_CSV)
    parser.add_argument("--out-dir", type=Path, default=COSMOS_PLOTS)
    args = parser.parse_args()

    df = pd.read_csv(args.hst_csv)
    mag = pd.to_numeric(df["ACS_MAG_AUTO"], errors="coerce").to_numpy(dtype=float)
    ba = pd.to_numeric(df["b_a"], errors="coerce").to_numpy(dtype=float)
    ok = np.isfinite(mag) & np.isfinite(ba) & (ba > 0) & (ba <= 1)
    mag, ba = mag[ok], ba[ok]
    cosi = cosi_from_ba(ba)

    ba_iso, cosi_iso = isotropic_ba(500_000)
    print_pool_stats("GIM2D disk entire (TYPE=2)", ba, cosi, mag)
    print_pool_stats("Isotropic simulation", ba_iso, cosi_iso)

    if "TYPE" in df.columns:
        sub = df.loc[ok].copy()
        sub["cosi"] = cosi
        print("--- GIM2D TYPE (median cos i) ---")
        type_names = {1: "bulge-dom", 2: "disk/spiral", 3: "irregular", 9: "unknown/fit-fail"}
        for t in sorted(sub["TYPE"].dropna().unique()):
            s = sub.loc[sub["TYPE"] == t]
            if len(s) < 50:
                continue
            name = type_names.get(int(t) if t == int(t) else t, str(t))
            print(
                f"  TYPE {int(t) if t == int(t) else t} ({name}): "
                f"N={len(s):,}, median cos(i)={s['cosi'].median():.3f}, "
                f"median b/a={s['b_a'].median():.3f}"
            )

    plot_mag_ba_scatter(
        mag,
        ba,
        cosi,
        out_png=args.out_dir / "gim2d_mag_ba_scatter.png",
        title=f"Zurich GIM2D disk (TYPE=2): magnitude vs b/a (N={len(ba):,})",
    )
    plot_cosi_histogram(cosi, cosi_iso, args.out_dir / "gim2d_cosi_vs_isotropic.png")

    if SDSS_DISK_ENTIRE_CSV.is_file():
        sdss = pd.read_csv(SDSS_DISK_ENTIRE_CSV)
        sm = pd.to_numeric(sdss["modelMag_r"], errors="coerce").to_numpy(dtype=float)
        sb = pd.to_numeric(sdss["expAB_r"], errors="coerce").to_numpy(dtype=float)
        ok_s = np.isfinite(sm) & np.isfinite(sb) & (sb > 0) & (sb <= 1)
        sc = cosi_from_ba(sb[ok_s])
        print_pool_stats("SDSS disk entire (all mags)", sb[ok_s], sc, sm[ok_s])
        print_pool_stats("SDSS disk entire mag<=22", sb[ok_s & (sm <= 22)], sc[sm[ok_s] <= 22])


if __name__ == "__main__":
    main()
