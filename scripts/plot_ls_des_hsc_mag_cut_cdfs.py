"""
cos(i) ECDFs at mag <= 20, 21, 22 for LS (scaled), DES, and HSC — plus overlays.

LS only uses the face-on cap + cos(i) rescaling (REX / near-round issue):
  - b/a > q0, b/a <= 0.8
  - cos(i)' = cos(i) / cos(i)|_{b/a=0.8}

DES / HSC (no REX gate — do NOT scale):
  - b/a > q0
  - plain Hubble cos(i)

Catalogs:
  LS  — Tractor type=EXP (v2 fullsky)
  DES — Y1 morph EXP analogue (0.4 < n < 1.5), calibrated ba/mag
  HSC — Kawinwanichakij EXP analogue (0.4 < n < 1.5)

Outputs:
  plots/plots_null/v2/ls_audit/cdfs/strict_scaled/
  plots/plots_null/v2/des_audit/cdfs/strict/
  plots/plots_null/v2/hsc_audit/cdfs/strict/
  plots/plots_null/v2/comparisons/

Run from repo root::

    python scripts/plot_ls_des_hsc_mag_cut_cdfs.py
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
    DES_Y1_MORPH_EXP_DEFAULT,
    HSC_KAWIN_EXP_DEFAULT,
    LS_CATALOG_V2_EXP_DEFAULT,
    Q0,
    hubble_cosi_from_ba,
)
from pipeline_null_plot_utils import REPO_ROOT  # noqa: E402

MAG_CUTS = (20.0, 21.0, 22.0)
BA_FACE_CAP = 0.8
PLOTS_NULL = REPO_ROOT / "plots" / "plots_null"

SURVEYS = {
    "ls": {
        "label": "LS EXP (scaled)",
        "short": "LS",
        "csv": REPO_ROOT / LS_CATALOG_V2_EXP_DEFAULT,
        "out": PLOTS_NULL / "v2" / "ls_audit" / "cdfs" / "strict_scaled",
        "color": "#377eb8",
        "scaled": True,
        "mag_col": None,
        "ba_col": None,
    },
    "des": {
        "label": "DES EXP-analogue",
        "short": "DES",
        "csv": REPO_ROOT / DES_Y1_MORPH_EXP_DEFAULT,
        "out": PLOTS_NULL / "v2" / "des_audit" / "cdfs" / "strict",
        "color": "#e41a1c",
        "scaled": False,
        "mag_col": "mag_r",
        "ba_col": "ba_r",
    },
    "hsc": {
        "label": "HSC EXP-analogue",
        "short": "HSC",
        "csv": REPO_ROOT / HSC_KAWIN_EXP_DEFAULT,
        "out": PLOTS_NULL / "v2" / "hsc_audit" / "cdfs" / "strict",
        "color": "#4daf4a",
        "scaled": False,
        "mag_col": "mag",
        "ba_col": "ba",
    },
}
COMPARE_OUT = PLOTS_NULL / "v2" / "comparisons"


def q_from_e1e2(e1: np.ndarray, e2: np.ndarray) -> np.ndarray:
    eabs = np.hypot(np.asarray(e1, dtype=float), np.asarray(e2, dtype=float))
    q = np.full(eabs.shape, np.nan, dtype=float)
    good = np.isfinite(eabs) & (eabs < 1.0)
    q[good] = (1.0 - eabs[good]) / (1.0 + eabs[good])
    return q


def empirical_cdf(vals: np.ndarray, x: np.ndarray) -> np.ndarray:
    s = np.sort(vals)
    return np.searchsorted(s, x, side="right") / max(1, len(s))


def cosi_array(ba: np.ndarray) -> np.ndarray:
    ba = np.asarray(ba, dtype=float)
    val = (ba * ba - Q0 * Q0) / (1.0 - Q0 * Q0)
    out = np.zeros_like(val)
    ok = np.isfinite(val) & (val >= 0.0)
    out[ok] = np.sqrt(np.clip(val[ok], 0.0, 1.0))
    return out


def load_ls(path: Path) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(path, usecols=["modelMag_r", "shape_e1", "shape_e2"])
    mag = pd.to_numeric(df["modelMag_r"], errors="coerce").to_numpy(dtype=float)
    ba = q_from_e1e2(
        pd.to_numeric(df["shape_e1"], errors="coerce").to_numpy(dtype=float),
        pd.to_numeric(df["shape_e2"], errors="coerce").to_numpy(dtype=float),
    )
    return mag, ba


def load_ba_mag(path: Path, mag_col: str, ba_col: str) -> tuple[np.ndarray, np.ndarray]:
    df = pd.read_csv(path, usecols=[mag_col, ba_col])
    mag = pd.to_numeric(df[mag_col], errors="coerce").to_numpy(dtype=float)
    ba = pd.to_numeric(df[ba_col], errors="coerce").to_numpy(dtype=float)
    return mag, ba


def apply_pool(
    mag: np.ndarray, ba: np.ndarray, *, scaled: bool
) -> tuple[np.ndarray, np.ndarray, float]:
    """Strict ba>q0; LS-only also caps ba<=0.8 and returns cosi_scale."""
    if scaled:
        ok = np.isfinite(mag) & np.isfinite(ba) & (ba > Q0) & (ba <= BA_FACE_CAP)
        cosi_scale = float(hubble_cosi_from_ba(BA_FACE_CAP, q0=Q0))
    else:
        ok = np.isfinite(mag) & np.isfinite(ba) & (ba > Q0) & (ba <= 1.0)
        cosi_scale = 1.0
    return mag[ok], ba[ok], cosi_scale


def xlab_for(scaled: bool) -> str:
    if scaled:
        return r"$\cos(i)\,/\,\cos(i)|_{b/a=0.8}$"
    return r"$\cos(i)$"


def plot_one(
    cosi: np.ndarray,
    mag_limit: float,
    out_path: Path,
    *,
    label: str,
    color: str,
    scaled: bool,
) -> dict:
    x = np.linspace(0.0, 1.0, 401)
    y = empirical_cdf(cosi, x)
    med = float(np.median(cosi)) if len(cosi) else float("nan")
    n = len(cosi)

    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    ax.plot(x, y, color=color, linewidth=2.0, label=label)
    ax.plot((0, 1), (0, 1), "k--", linewidth=1.2, label="Uniform")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(xlab_for(scaled))
    ax.set_ylabel("Cumulative distribution")
    ax.set_title(f"{label}  |  mag <= {mag_limit:g}\nN={n:,}  |  median={med:.3f}")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return {"mag_limit": mag_limit, "n": n, "median_cosi": round(med, 4) if n else np.nan}


def plot_survey_overlay(
    panels: list[tuple[float, np.ndarray]],
    out_path: Path,
    *,
    label: str,
    color: str,
    scaled: bool,
) -> None:
    x = np.linspace(0.0, 1.0, 401)
    mag_colors = {20.0: "#1b9e77", 21.0: "#d95f02", 22.0: "#7570b3"}
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot((0, 1), (0, 1), "k--", linewidth=1.2, label="Uniform", zorder=1)
    for mag_limit, cosi in panels:
        if len(cosi) == 0:
            continue
        y = empirical_cdf(cosi, x)
        med = float(np.median(cosi))
        ax.plot(
            x,
            y,
            color=mag_colors.get(mag_limit, color),
            linewidth=2.0,
            label=f"<= {mag_limit:g}  (N={len(cosi):,}, med={med:.3f})",
        )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(xlab_for(scaled))
    ax.set_ylabel("Cumulative distribution")
    ax.set_title(label)
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_compare(
    series: dict[str, np.ndarray],
    mag_limit: float,
    out_path: Path,
) -> None:
    x = np.linspace(0.0, 1.0, 401)
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot((0, 1), (0, 1), "k--", linewidth=1.2, label="Uniform", zorder=1)
    for key, cosi in series.items():
        cfg = SURVEYS[key]
        if len(cosi) == 0:
            continue
        y = empirical_cdf(cosi, x)
        med = float(np.median(cosi))
        tag = " scaled" if cfg["scaled"] else ""
        ax.plot(
            x,
            y,
            color=cfg["color"],
            linewidth=2.0,
            label=f"{cfg['short']}{tag}  (N={len(cosi):,}, med={med:.3f})",
        )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"$\cos(i)$  (LS only: scaled to 1 at $b/a=0.8$)")
    ax.set_ylabel("Cumulative distribution")
    ax.set_title(
        f"LS (scaled) / DES / HSC  |  mag <= {mag_limit:g}\n"
        r"DES+HSC: plain Hubble $\cos(i)$, $b/a>q_0$; LS only: $b/a\leq0.8$ + rescale"
    )
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def run_survey(key: str) -> dict[float, np.ndarray]:
    cfg = SURVEYS[key]
    print(f"[*] {cfg['short']}: loading {cfg['csv'].name} ...", flush=True)
    if key == "ls":
        mag, ba = load_ls(cfg["csv"])
    else:
        mag, ba = load_ba_mag(cfg["csv"], cfg["mag_col"], cfg["ba_col"])

    mag, ba, cosi_scale = apply_pool(mag, ba, scaled=cfg["scaled"])
    if cfg["scaled"]:
        print(
            f"    LS scaled pool {Q0:g}<ba<={BA_FACE_CAP:g}: N={len(mag):,}  "
            f"cosi_scale={cosi_scale:.4f}",
            flush=True,
        )
    else:
        print(f"    strict pool ba>{Q0:g}: N={len(mag):,}  (no ba cap / no scale)", flush=True)

    out_root = cfg["out"]
    out_root.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    overlay: list[tuple[float, np.ndarray]] = []
    by_mag: dict[float, np.ndarray] = {}

    for lim in MAG_CUTS:
        mask = mag <= lim
        cosi = cosi_array(ba[mask])
        if cfg["scaled"]:
            cosi = np.clip(cosi / cosi_scale, 0.0, 1.0)
        by_mag[lim] = cosi
        row = plot_one(
            cosi,
            lim,
            out_root / f"mag{int(lim)}.png",
            label=cfg["label"],
            color=cfg["color"],
            scaled=cfg["scaled"],
        )
        row["survey"] = cfg["short"]
        rows.append(row)
        overlay.append((lim, cosi))
        print(
            f"    mag<={lim:g}: N={len(cosi):,}  median={row['median_cosi']}",
            flush=True,
        )

    plot_survey_overlay(
        overlay,
        out_root / "overlay.png",
        label=cfg["label"],
        color=cfg["color"],
        scaled=cfg["scaled"],
    )
    pd.DataFrame(rows).to_csv(out_root / "summary.csv", index=False)
    print(f"    wrote {out_root}", flush=True)
    return by_mag


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--survey", choices=("ls", "des", "hsc", "all"), default="all")
    args = parser.parse_args()
    keys = ("ls", "des", "hsc") if args.survey == "all" else (args.survey,)

    all_by_mag: dict[str, dict[float, np.ndarray]] = {}
    for k in keys:
        all_by_mag[k] = run_survey(k)

    if args.survey == "all":
        COMPARE_OUT.mkdir(parents=True, exist_ok=True)
        compare_rows: list[dict] = []
        for lim in MAG_CUTS:
            series = {k: all_by_mag[k][lim] for k in ("ls", "des", "hsc")}
            plot_compare(series, lim, COMPARE_OUT / f"mag{int(lim)}.png")
            for k, cosi in series.items():
                compare_rows.append(
                    {
                        "survey": SURVEYS[k]["short"],
                        "scaled": SURVEYS[k]["scaled"],
                        "mag_limit": lim,
                        "n": len(cosi),
                        "median_cosi": round(float(np.median(cosi)), 4) if len(cosi) else np.nan,
                    }
                )
            print(f"[*] compare mag<={lim:g} -> {COMPARE_OUT / f'mag{int(lim)}.png'}", flush=True)

        pd.DataFrame(compare_rows).to_csv(COMPARE_OUT / "summary.csv", index=False)

        fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
        x = np.linspace(0.0, 1.0, 401)
        for ax, lim in zip(axes, MAG_CUTS):
            ax.plot((0, 1), (0, 1), "k--", linewidth=1.0, label="Uniform")
            for k in ("ls", "des", "hsc"):
                cosi = all_by_mag[k][lim]
                if len(cosi) == 0:
                    continue
                cfg = SURVEYS[k]
                y = empirical_cdf(cosi, x)
                med = float(np.median(cosi))
                tag = " scaled" if cfg["scaled"] else ""
                ax.plot(
                    x,
                    y,
                    color=cfg["color"],
                    linewidth=2.0,
                    label=f"{cfg['short']}{tag} (N={len(cosi):,}, med={med:.3f})",
                )
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_xlabel(r"$\cos(i)$")
            ax.set_title(f"mag <= {lim:g}")
            ax.grid(True, alpha=0.3)
            ax.legend(loc="upper left", fontsize=8)
        axes[0].set_ylabel("Cumulative distribution")
        fig.suptitle(
            "LS scaled vs DES/HSC plain Hubble cos(i)  |  EXP / EXP-analogue",
            fontsize=12,
        )
        fig.tight_layout()
        fig.savefig(COMPARE_OUT / "overlay_all.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"[*] Wrote {COMPARE_OUT / 'overlay_all.png'}", flush=True)


if __name__ == "__main__":
    main()
