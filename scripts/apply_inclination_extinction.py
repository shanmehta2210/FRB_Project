"""
Dust-selection CDFs: raw vs Unterborn A1 vs survival trim (strict ba > q0).

All cos(i) pools require ``b/a > q0`` (q0=0.2) so Hubble cos(i) is defined.
LS also caps ``b/a <= 0.8`` and rescales cos(i) (REX); that is separate from dust.

Laws (Unterborn & Ryden 2008): Δm = 1.27 (log10 q)^2
  A1:        keep m^f = m - Δm(q)  <= lim   (adds edge-ons)
  survival:  keep m_edge = m^f + Δm(q0) <= lim  (drops near-limit face-ons)

Outputs under ``plots/plots_null/v2/extinction/``.

Run from repo root::

    python scripts/apply_inclination_extinction.py
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
    delta_m_unterborn,
    edge_on_mag,
    face_on_mag,
    hubble_cosi_from_ba,
)
from pipeline_null_plot_utils import REPO_ROOT  # noqa: E402

MAG_CUTS = (20.0, 21.0, 22.0)
BA_FACE_CAP = 0.8
OUT_ROOT = REPO_ROOT / "plots" / "plots_null" / "v2" / "extinction"
# Edge-on reference for survival trim = Hubble floor (same q0 as cos(i)).
Q_EDGE = Q0
DM_EDGE = float(delta_m_unterborn(Q_EDGE))

SURVEYS = {
    "ls": {
        "short": "LS",
        "label": "LS EXP (scaled, strict)",
        "csv": REPO_ROOT / LS_CATALOG_V2_EXP_DEFAULT,
        "color": "#377eb8",
        "scaled": True,
        "mag_col": None,
        "ba_col": None,
    },
    "des": {
        "short": "DES",
        "label": "DES EXP-analogue (strict)",
        "csv": REPO_ROOT / DES_Y1_MORPH_EXP_DEFAULT,
        "color": "#e41a1c",
        "scaled": False,
        "mag_col": "mag_r",
        "ba_col": "ba_r",
    },
    "hsc": {
        "short": "HSC",
        "label": "HSC EXP-analogue (strict)",
        "csv": REPO_ROOT / HSC_KAWIN_EXP_DEFAULT,
        "color": "#4daf4a",
        "scaled": False,
        "mag_col": "mag",
        "ba_col": "ba",
    },
}


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
    """Hubble cos(i); caller must already enforce ba > q0."""
    ba = np.asarray(ba, dtype=float)
    if len(ba) and np.any(ba <= Q0):
        raise ValueError(f"strict pool violated: found ba <= q0={Q0:g}")
    val = (ba * ba - Q0 * Q0) / (1.0 - Q0 * Q0)
    return np.sqrt(np.clip(val, 0.0, 1.0))


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


def apply_strict_pool(
    mag: np.ndarray, ba: np.ndarray, *, scaled: bool
) -> tuple[np.ndarray, np.ndarray, float]:
    """Strict Hubble pool: ba > q0. LS also ba <= 0.8 for REX scaling."""
    if scaled:
        ok = np.isfinite(mag) & np.isfinite(ba) & (ba > Q0) & (ba <= BA_FACE_CAP)
        cosi_scale = float(hubble_cosi_from_ba(BA_FACE_CAP, q0=Q0))
    else:
        ok = np.isfinite(mag) & np.isfinite(ba) & (ba > Q0) & (ba <= 1.0)
        cosi_scale = 1.0
    mag_o, ba_o = mag[ok], ba[ok]
    if len(ba_o) == 0 or float(np.min(ba_o)) <= Q0:
        raise RuntimeError("strict ba > q0 pool empty or contaminated")
    return mag_o, ba_o, cosi_scale


def make_cosi(ba: np.ndarray, *, cosi_scale: float) -> np.ndarray:
    cosi = cosi_array(ba)
    if cosi_scale != 1.0:
        cosi = np.clip(cosi / cosi_scale, 0.0, 1.0)
    return cosi


def xlab_for(scaled: bool) -> str:
    if scaled:
        return rf"$\cos(i)\,/\,\cos(i)|_{{b/a={BA_FACE_CAP:g}}}$  (strict $b/a>q_0$)"
    return r"$\cos(i)$  (strict $b/a>q_0$)"


def plot_three_way(
    modes: list[tuple[str, np.ndarray, str, float]],
    mag_limit: float,
    out_path: Path,
    *,
    label: str,
    color: str,
    scaled: bool,
) -> None:
    """modes: (name, cosi, linestyle, alpha)"""
    x = np.linspace(0.0, 1.0, 401)
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot((0, 1), (0, 1), "k--", linewidth=1.2, label="Uniform", zorder=1)
    for name, cosi, ls, alpha in modes:
        if len(cosi) == 0:
            continue
        med = float(np.median(cosi))
        ax.plot(
            x,
            empirical_cdf(cosi, x),
            color=color,
            linestyle=ls,
            linewidth=2.0,
            alpha=alpha,
            label=f"{name}  (N={len(cosi):,}, med={med:.3f})",
        )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(xlab_for(scaled))
    ax.set_ylabel("Cumulative distribution")
    ax.set_title(
        f"{label}  |  limit = {mag_limit:g}\n"
        rf"strict $b/a > q_0={Q0:g}$  |  raw / A1 / survival"
    )
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_compare_mode(
    series: dict[str, np.ndarray],
    mag_limit: float,
    out_path: Path,
    *,
    mode_title: str,
) -> None:
    x = np.linspace(0.0, 1.0, 401)
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.plot((0, 1), (0, 1), "k--", linewidth=1.2, label="Uniform", zorder=1)
    for key, cosi in series.items():
        cfg = SURVEYS[key]
        if len(cosi) == 0:
            continue
        med = float(np.median(cosi))
        tag = " scaled" if cfg["scaled"] else ""
        ax.plot(
            x,
            empirical_cdf(cosi, x),
            color=cfg["color"],
            linewidth=2.0,
            label=f"{cfg['short']}{tag}  (N={len(cosi):,}, med={med:.3f})",
        )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"$\cos(i)$  (LS scaled; all strict $b/a>q_0$)")
    ax.set_ylabel("Cumulative distribution")
    ax.set_title(f"{mode_title}  |  limit = {mag_limit:g}")
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def run_survey(key: str) -> dict[float, dict[str, np.ndarray]]:
    cfg = SURVEYS[key]
    print(f"[*] {cfg['short']}: loading {cfg['csv'].name} ...", flush=True)
    if key == "ls":
        mag, ba = load_ls(cfg["csv"])
    else:
        mag, ba = load_ba_mag(cfg["csv"], cfg["mag_col"], cfg["ba_col"])

    mag, ba, cosi_scale = apply_strict_pool(mag, ba, scaled=cfg["scaled"])
    m_face = face_on_mag(mag, ba)
    m_edge = edge_on_mag(mag, ba, q_edge=Q_EDGE)
    print(
        f"    STRICT pool ba>{Q0:g}"
        + (f" & ba<={BA_FACE_CAP:g}" if cfg["scaled"] else "")
        + f": N={len(mag):,}  min_ba={float(np.min(ba)):.4f}  "
        f"med_ba={float(np.median(ba)):.3f}  med_dm={float(np.median(mag - m_face)):.3f}",
        flush=True,
    )

    survey_out = OUT_ROOT / "cdfs" / key
    survey_out.mkdir(parents=True, exist_ok=True)
    by_lim: dict[float, dict[str, np.ndarray]] = {}

    for lim in MAG_CUTS:
        raw_mask = mag <= lim
        a1_mask = m_face <= lim
        # Survival: would still pass lim if viewed at q_edge (= q0).
        surv_mask = m_edge <= lim

        packs = {
            "raw": (ba[raw_mask], mag[raw_mask]),
            "a1": (ba[a1_mask], mag[a1_mask]),
            "survival": (ba[surv_mask], mag[surv_mask]),
        }
        out: dict[str, np.ndarray] = {}
        for mode, (ba_m, _) in packs.items():
            out[f"ba_{mode}"] = ba_m
            out[mode] = make_cosi(ba_m, cosi_scale=cosi_scale)
        by_lim[lim] = out

        plot_three_way(
            [
                ("raw $m$", out["raw"], "-", 0.45),
                (r"A1 $m^f$", out["a1"], "-", 1.0),
                (rf"survival $m_{{\rm edge}}$", out["survival"], "--", 0.9),
            ],
            lim,
            survey_out / f"mag{int(lim)}_before_after.png",
            label=cfg["label"],
            color=cfg["color"],
            scaled=cfg["scaled"],
        )
        # Also keep a dedicated survival vs raw panel name for clarity
        plot_three_way(
            [
                ("raw $m$", out["raw"], "-", 0.45),
                (rf"survival $m_{{\rm edge}}$", out["survival"], "-", 1.0),
            ],
            lim,
            survey_out / f"mag{int(lim)}_survival.png",
            label=cfg["label"],
            color=cfg["color"],
            scaled=cfg["scaled"],
        )

        print(
            f"    lim={lim:g}: raw N={len(out['raw']):,} med={np.median(out['raw']):.4f}  |  "
            f"A1 N={len(out['a1']):,} med={np.median(out['a1']):.4f}  |  "
            f"surv N={len(out['survival']):,} med={np.median(out['survival']):.4f}  "
            f"(surv/raw={len(out['survival']) / max(1, len(out['raw'])):.2f})",
            flush=True,
        )

    return by_lim


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--survey", choices=("ls", "des", "hsc", "all"), default="all")
    args = parser.parse_args()
    keys = ("ls", "des", "hsc") if args.survey == "all" else (args.survey,)

    print(
        f"[*] strict ba > q0={Q0:g}; survival q_edge={Q_EDGE:g} "
        f"(dm_edge={DM_EDGE:.3f} mag)",
        flush=True,
    )

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    all_data: dict[str, dict[float, dict[str, np.ndarray]]] = {}
    for k in keys:
        all_data[k] = run_survey(k)

    summary_rows: list[dict] = []
    funnel_rows: list[dict] = []

    for k in keys:
        cfg = SURVEYS[k]
        for lim in MAG_CUTS:
            d = all_data[k][lim]
            for mode in ("raw", "a1", "survival"):
                cosi = d[mode]
                ba = d[f"ba_{mode}"]
                summary_rows.append(
                    {
                        "survey": cfg["short"],
                        "scaled": cfg["scaled"],
                        "strict_ba_gt_q0": True,
                        "q0": Q0,
                        "mag_limit": lim,
                        "mode": mode,
                        "n": len(cosi),
                        "median_cosi": round(float(np.median(cosi)), 4) if len(cosi) else np.nan,
                        "median_ba": round(float(np.median(ba)), 4) if len(ba) else np.nan,
                        "min_ba": round(float(np.min(ba)), 4) if len(ba) else np.nan,
                    }
                )
            funnel_rows.append(
                {
                    "survey": cfg["short"],
                    "mag_limit": lim,
                    "n_raw": len(d["raw"]),
                    "n_a1": len(d["a1"]),
                    "n_survival": len(d["survival"]),
                    "n_a1_minus_raw": len(d["a1"]) - len(d["raw"]),
                    "n_surv_minus_raw": len(d["survival"]) - len(d["raw"]),
                    "frac_surv_of_raw": round(
                        len(d["survival"]) / max(1, len(d["raw"])), 4
                    ),
                    "median_cosi_raw": round(float(np.median(d["raw"])), 4)
                    if len(d["raw"])
                    else np.nan,
                    "median_cosi_a1": round(float(np.median(d["a1"])), 4)
                    if len(d["a1"])
                    else np.nan,
                    "median_cosi_survival": round(float(np.median(d["survival"])), 4)
                    if len(d["survival"])
                    else np.nan,
                }
            )

    pd.DataFrame(summary_rows).to_csv(OUT_ROOT / "summary_before_after.csv", index=False)
    pd.DataFrame(funnel_rows).to_csv(OUT_ROOT / "funnel.csv", index=False)

    if args.survey == "all":
        compare_out = OUT_ROOT / "compare"
        compare_out.mkdir(parents=True, exist_ok=True)
        for lim in MAG_CUTS:
            for mode, title, fname in (
                ("a1", rf"A1 $m^f$ (strict $b/a>q_0$)", f"mag{int(lim)}_a1.png"),
                (
                    "survival",
                    rf"survival $m_{{\rm edge}}$ (strict $b/a>q_0$)",
                    f"mag{int(lim)}_survival.png",
                ),
            ):
                series = {k: all_data[k][lim][mode] for k in ("ls", "des", "hsc")}
                plot_compare_mode(series, lim, compare_out / fname, mode_title=title)
            print(f"[*] compare mag<={lim:g} (A1 + survival)", flush=True)

        for mode, fname, supt in (
            ("a1", "overlay_all_a1.png", r"A1 $m^f$ re-cut"),
            ("survival", "overlay_all_survival.png", r"survival $m_{\rm edge}$ trim"),
        ):
            fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
            x = np.linspace(0.0, 1.0, 401)
            for ax, lim in zip(axes, MAG_CUTS):
                ax.plot((0, 1), (0, 1), "k--", linewidth=1.0, label="Uniform")
                for k in ("ls", "des", "hsc"):
                    cosi = all_data[k][lim][mode]
                    if len(cosi) == 0:
                        continue
                    cfg = SURVEYS[k]
                    med = float(np.median(cosi))
                    tag = " scaled" if cfg["scaled"] else ""
                    ax.plot(
                        x,
                        empirical_cdf(cosi, x),
                        color=cfg["color"],
                        linewidth=2.0,
                        label=f"{cfg['short']}{tag} (N={len(cosi):,}, med={med:.3f})",
                    )
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
                ax.set_xlabel(r"$\cos(i)$")
                ax.set_title(f"limit = {lim:g}")
                ax.grid(True, alpha=0.3)
                ax.legend(loc="upper left", fontsize=8)
            axes[0].set_ylabel("Cumulative distribution")
            fig.suptitle(
                f"{supt}  |  strict $b/a>q_0={Q0:g}$  |  LS scaled / DES / HSC",
                fontsize=12,
            )
            fig.tight_layout()
            fig.savefig(OUT_ROOT / fname, dpi=300, bbox_inches="tight")
            plt.close(fig)

    q = np.linspace(0.05, 1.0, 200)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(q, 1.27 * (np.log10(q) ** 2), color="#333333", lw=2, label=r"$\Delta m(q)$")
    ax.axvline(Q0, color="#e41a1c", ls="--", lw=1.2, label=rf"$q_0={Q0:g}$ (survival edge)")
    ax.axhline(DM_EDGE, color="#e41a1c", ls=":", lw=1.0)
    ax.set_xlabel(r"$b/a$")
    ax.set_ylabel(r"$\Delta m_r = 1.27(\log_{10} q)^2$")
    ax.set_title("Unterborn & Ryden dimming law")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT_ROOT / "delta_m_vs_ba.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"[*] Wrote {OUT_ROOT}", flush=True)
    print(pd.DataFrame(funnel_rows).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
