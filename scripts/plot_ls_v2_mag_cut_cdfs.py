"""
Legacy v2 EXP: cos(i) ECDFs at modelMag_r <= 20, 21, 22.

b/a from Tractor e1,e2; Hubble cos(i) with q0=0.2.

Modes:
  strict         — pool b/a > q0  (original); writes under cdfs/strict/
  strict_scaled  — pool q0 < b/a <= 0.8; CDF of cos(i)/cos(i)|_{b/a=0.8}
                   so the truncated face-on edge (b/a=0.8) maps to 1.
                   Writes under cdfs/strict_scaled/

Run from repo root::

    python scripts/plot_ls_v2_mag_cut_cdfs.py
    python scripts/plot_ls_v2_mag_cut_cdfs.py --mode strict_scaled
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

from null_catalog_utils import (  # noqa: E402
    LS_CATALOG_V2_EXP_DEFAULT,
    Q0,
    hubble_cosi_from_ba,
)
from pipeline_null_plot_utils import REPO_ROOT  # noqa: E402

CDFS_ROOT = REPO_ROOT / "plots" / "plots_null" / "v2" / "ls_audit" / "cdfs"
DEFAULT_CSV = REPO_ROOT / LS_CATALOG_V2_EXP_DEFAULT
MAG_CUTS = (20.0, 21.0, 22.0)
USECOLS = ("modelMag_r", "shape_e1", "shape_e2")
BA_FACE_CAP = 0.8  # strict_scaled: drop b/a above this (REX-like face-on vicinity)


def q_from_e1e2(e1: np.ndarray, e2: np.ndarray) -> np.ndarray:
    eabs = np.hypot(np.asarray(e1, dtype=float), np.asarray(e2, dtype=float))
    q = np.full(eabs.shape, np.nan, dtype=float)
    good = np.isfinite(eabs) & (eabs < 1.0)
    q[good] = (1.0 - eabs[good]) / (1.0 + eabs[good])
    return q


def empirical_cdf(vals: np.ndarray, x: np.ndarray) -> np.ndarray:
    s = np.sort(vals)
    return np.searchsorted(s, x, side="right") / max(1, len(s))


def load_ba_mag(path: Path) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(path, usecols=list(USECOLS))
    mag = pd.to_numeric(df["modelMag_r"], errors="coerce").to_numpy(dtype=float)
    ba = q_from_e1e2(
        pd.to_numeric(df["shape_e1"], errors="coerce").to_numpy(dtype=float),
        pd.to_numeric(df["shape_e2"], errors="coerce").to_numpy(dtype=float),
    )
    ok = np.isfinite(mag) & np.isfinite(ba) & (ba > Q0) & (ba <= 1.0)
    return mag[ok], ba[ok]


def cosi_array(ba: np.ndarray) -> np.ndarray:
    return np.array([hubble_cosi_from_ba(float(q), q0=Q0) for q in ba], dtype=float)


def plot_one(
    cosi: np.ndarray,
    mag_limit: float,
    out_path: Path,
    *,
    mode: str,
    cosi_scale: float,
) -> dict:
    x = np.linspace(0.0, 1.0, 401)
    y = empirical_cdf(cosi, x)
    med = float(np.median(cosi))
    n = len(cosi)

    if mode == "strict_scaled":
        xlab = r"$\cos(i)\,/\,\cos(i)|_{b/a=0.8}$"
        title = (
            f"Legacy v2 EXP (scaled)  |  modelMag$_r$ ≤ {mag_limit:g}\n"
            f"N={n:,}  |  median scaled cos(i)={med:.3f}  |  "
            rf"${Q0:g}<b/a\leq{BA_FACE_CAP:g}$; "
            rf"scale cos$(i)|_{{b/a={BA_FACE_CAP:g}}}={cosi_scale:.3f}$"
        )
        label = "Legacy EXP (ba≤0.8, scaled)"
    else:
        xlab = r"$\cos(i)$"
        title = (
            f"Legacy v2 EXP  |  modelMag$_r$ ≤ {mag_limit:g}\n"
            f"N={n:,}  |  median cos(i)={med:.3f}  |  "
            r"$b/a$ from $e_1,e_2$, $b/a>q_0$"
        )
        label = "Legacy EXP null"

    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    ax.plot(x, y, color="#377eb8", linewidth=2.0, label=label)
    ax.plot((0, 1), (0, 1), "k--", linewidth=1.2, label="Uniform")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(xlab)
    ax.set_ylabel("Cumulative distribution")
    ax.set_title(title)
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return {"mag_limit": mag_limit, "n": n, "median_cosi": round(med, 4)}


def plot_overlay(
    panels: list[tuple[float, np.ndarray]],
    out_path: Path,
    *,
    mode: str,
    cosi_scale: float,
) -> None:
    x = np.linspace(0.0, 1.0, 401)
    colors = {20.0: "#1b9e77", 21.0: "#d95f02", 22.0: "#7570b3"}
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot((0, 1), (0, 1), "k--", linewidth=1.2, label="Uniform", zorder=1)
    for mag_limit, cosi in panels:
        y = empirical_cdf(cosi, x)
        med = float(np.median(cosi))
        ax.plot(
            x,
            y,
            color=colors.get(mag_limit, "#377eb8"),
            linewidth=2.0,
            label=f"≤ {mag_limit:g}  (N={len(cosi):,}, med={med:.3f})",
        )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    if mode == "strict_scaled":
        ax.set_xlabel(r"$\cos(i)\,/\,\cos(i)|_{b/a=0.8}$")
        ax.set_title(
            "Legacy v2 EXP null cos(i) ECDFs (scaled)\n"
            rf"${Q0:g}<b/a\leq{BA_FACE_CAP:g}$; "
            rf"divided by $\cos(i)|_{{b/a={BA_FACE_CAP:g}}}$ $={cosi_scale:.3f}$"
        )
    else:
        ax.set_xlabel(r"$\cos(i)$")
        ax.set_title(
            "Legacy v2 EXP null cos(i) ECDFs\n"
            r"$b/a=(1-|e|)/(1+|e|)$ from $e_1,e_2$; Hubble $q_0=0.2$; $b/a>q_0$"
        )
    ax.set_ylabel("Cumulative distribution")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument(
        "--mode",
        choices=("strict", "strict_scaled"),
        default="strict",
        help="strict: b/a>q0; strict_scaled: q0<b/a<=0.8 with cos(i) rescaled to 1 at ba=0.8",
    )
    parser.add_argument(
        "--out-root",
        type=Path,
        default=None,
        help="Default: plots/.../cdfs/<mode>/",
    )
    args = parser.parse_args()
    out_root = args.out_root if args.out_root is not None else CDFS_ROOT / args.mode

    print(f"[*] Loading {args.csv} ...", flush=True)
    mag, ba = load_ba_mag(args.csv)
    print(f"[*] Hubble-valid pool (b/a > {Q0:g}): N={len(mag):,}", flush=True)

    cosi_scale = 1.0
    if args.mode == "strict_scaled":
        keep = ba <= BA_FACE_CAP
        mag, ba = mag[keep], ba[keep]
        cosi_scale = float(hubble_cosi_from_ba(BA_FACE_CAP, q0=Q0))
        print(
            f"[*] strict_scaled: kept b/a ≤ {BA_FACE_CAP:g}: N={len(mag):,}  "
            f"cos(i)|_ba={BA_FACE_CAP:g} = {cosi_scale:.4f}",
            flush=True,
        )

    rows: list[dict] = []
    overlay: list[tuple[float, np.ndarray]] = []
    for lim in MAG_CUTS:
        mask = mag <= lim
        cosi = cosi_array(ba[mask])
        if args.mode == "strict_scaled":
            cosi = cosi / cosi_scale
            cosi = np.clip(cosi, 0.0, 1.0)
        out = out_root / f"mag{int(lim)}.png"
        rows.append(plot_one(cosi, lim, out, mode=args.mode, cosi_scale=cosi_scale))
        overlay.append((lim, cosi))
        print(
            f"  mag≤{lim:g}: N={len(cosi):,}  median={rows[-1]['median_cosi']}",
            flush=True,
        )

    plot_overlay(overlay, out_root / "overlay.png", mode=args.mode, cosi_scale=cosi_scale)
    summary = pd.DataFrame(rows)
    summary.to_csv(out_root / "summary.csv", index=False)
    print(f"[*] Wrote {out_root}", flush=True)


if __name__ == "__main__":
    main()
