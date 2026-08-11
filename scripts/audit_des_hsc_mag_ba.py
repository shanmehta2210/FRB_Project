"""
DES Y1 morph + HSC Kawinwanichakij — magnitude vs axis-ratio audit plots.

Mirrors ``scripts/audit_ls_v2_mag_ba.py`` outputs:
  - mag_histogram_ba.png / .csv
  - median_ba_vs_mag.png / .csv
  - ba_vs_mag_scatter.png
  - README.md

Default: full 500k samples under ``plots/plots_null/v2/des_audit/`` and
``plots/plots_null/v2/hsc_audit/``. Pass ``--exp-analogue`` for the n-window
EXP-analogue catalogs (written to ``v2/des_audit/exp_analogue/`` etc.).

Run from repo root::

    python scripts/audit_des_hsc_mag_ba.py
    python scripts/audit_des_hsc_mag_ba.py --exp-analogue
    python scripts/audit_des_hsc_mag_ba.py --survey des
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
    DES_Y1_MORPH_SAMPLE_DEFAULT,
    EXP_ANALOGUE_N_MAX,
    EXP_ANALOGUE_N_MIN,
    HSC_KAWIN_EXP_DEFAULT,
    HSC_KAWIN_SAMPLE_DEFAULT,
)
from pipeline_null_plot_utils import REPO_ROOT  # noqa: E402

MAG_CLIP = 15.0
BIN_WIDTH = 0.5
PLOT_MAG_MAX = 26.0
MIN_BIN_N = 1
SCATTER_OVERLAY_N = 40_000

SURVEYS = {
    "des": {
        "label": "DES Y1 morph (Tarsitano+2018)",
        "short": "DES Y1",
        "sample_csv": REPO_ROOT / DES_Y1_MORPH_SAMPLE_DEFAULT,
        "exp_csv": REPO_ROOT / DES_Y1_MORPH_EXP_DEFAULT,
        "out_sample": REPO_ROOT / "plots" / "plots_null" / "v2" / "des_audit",
        "out_exp": REPO_ROOT / "plots" / "plots_null" / "v2" / "des_audit" / "exp_analogue",
        "mag_col": "mag_r",
        "ba_col": "ba_r",
        "n_col": "n_r",
        "ra_col": "RA_ICRS",
        "dec_col": "DE_ICRS",
        "id_col": "COADD_OBJECTS_ID",
        "mag_xlabel": r"MAG_SERSIC$_r$ (calibrated)",
        "ba_ylabel": r"$b/a = 1 - \varepsilon_{\mathrm{Sersic},r}$",
        "ba_note": r"$b/a=1-\varepsilon$ (already calibrated; App. B)",
        "color": "#e41a1c",
    },
    "hsc": {
        "label": "HSC Kawinwanichakij+2021",
        "short": "HSC",
        "sample_csv": REPO_ROOT / HSC_KAWIN_SAMPLE_DEFAULT,
        "exp_csv": REPO_ROOT / HSC_KAWIN_EXP_DEFAULT,
        "out_sample": REPO_ROOT / "plots" / "plots_null" / "v2" / "hsc_audit",
        "out_exp": REPO_ROOT / "plots" / "plots_null" / "v2" / "hsc_audit" / "exp_analogue",
        "mag_col": "mag",
        "ba_col": "ba",
        "n_col": "n_sersic",
        "ra_col": "RA_ICRS",
        "dec_col": "DE_ICRS",
        "id_col": "object_id",
        "mag_xlabel": r"fitted mag$_r$",
        "ba_ylabel": r"$b/a$ (fitted_q)",
        "ba_note": r"$b/a=$ fitted_q",
        "color": "#4daf4a",
    },
}


def mag_bin_edges(mag: np.ndarray, *, clip: float, width: float) -> np.ndarray:
    finite = mag[np.isfinite(mag)]
    hi = float(np.ceil(finite.max() / width) * width) if len(finite) else clip + width
    edges: list[float] = [-np.inf, clip]
    m = clip
    while m < hi:
        m += width
        edges.append(m)
    return np.asarray(edges, dtype=float)


def mag_bin_labels(edges: np.ndarray) -> list[str]:
    labels: list[str] = []
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        if not np.isfinite(lo):
            labels.append(f"< {hi:g}")
        else:
            labels.append(f"{lo:g}–{hi:g}")
    return labels


def bin_mask(mag: np.ndarray, lo: float, hi: float, *, last: bool) -> np.ndarray:
    del last
    if not np.isfinite(lo):
        return np.isfinite(mag) & (mag <= hi)
    return np.isfinite(mag) & (mag > lo) & (mag <= hi)


def load_pool(path: Path, cfg: dict) -> pd.DataFrame:
    print(f"[*] Loading {path} ...", flush=True)
    usecols = [cfg["mag_col"], cfg["ba_col"], cfg["ra_col"], cfg["dec_col"], cfg["id_col"]]
    n_col = cfg["n_col"]
    # n is optional (present on both catalogs)
    peek = pd.read_csv(path, nrows=0)
    if n_col in peek.columns:
        usecols.append(n_col)
    df = pd.read_csv(path, usecols=usecols)
    mag = pd.to_numeric(df[cfg["mag_col"]], errors="coerce").to_numpy(dtype=float)
    ba = pd.to_numeric(df[cfg["ba_col"]], errors="coerce").to_numpy(dtype=float)
    ok = np.isfinite(mag) & np.isfinite(ba) & (ba >= 0.0) & (ba <= 1.0)
    out = {
        "id": df.loc[ok, cfg["id_col"]].to_numpy(),
        "mag": mag[ok],
        "ba": ba[ok],
        "RA_ICRS": pd.to_numeric(df.loc[ok, cfg["ra_col"]], errors="coerce").to_numpy(),
        "DE_ICRS": pd.to_numeric(df.loc[ok, cfg["dec_col"]], errors="coerce").to_numpy(),
    }
    if n_col in df.columns:
        out["n"] = pd.to_numeric(df.loc[ok, n_col], errors="coerce").to_numpy(dtype=float)
    print(
        f"[*] Finite mag + b/a: N={ok.sum():,}  (dropped {(~ok).sum():,})",
        flush=True,
    )
    return pd.DataFrame(out)


def make_bin_summary(mag: np.ndarray, ba: np.ndarray, edges: np.ndarray) -> pd.DataFrame:
    labels = mag_bin_labels(edges)
    rows: list[dict] = []
    n_tot = len(mag)
    for b in range(len(edges) - 1):
        lo, hi = edges[b], edges[b + 1]
        mask = bin_mask(mag, lo, hi, last=(b == len(edges) - 2))
        n = int(np.sum(mask))
        rows.append(
            {
                "bin_label": labels[b],
                "mag_lo": float(lo) if np.isfinite(lo) else np.nan,
                "mag_hi": float(hi),
                "n_galaxies": n,
                "frac_of_pool_pct": round(100.0 * n / max(1, n_tot), 4),
                "median_ba": round(float(np.median(ba[mask])), 4) if n else np.nan,
                "mean_ba": round(float(np.mean(ba[mask])), 4) if n else np.nan,
            }
        )
    return pd.DataFrame(rows)


def plot_display_summary(summary: pd.DataFrame, *, mag_max: float, min_n: int) -> pd.DataFrame:
    hi_ok = summary["mag_hi"] <= mag_max + 1e-9
    n_ok = summary["n_galaxies"] >= min_n
    return summary.loc[hi_ok & n_ok].copy()


def plot_mag_histogram(
    summary: pd.DataFrame,
    n_pool: int,
    out_png: Path,
    *,
    cfg: dict,
    subtitle: str,
) -> None:
    fig, ax = plt.subplots(figsize=(12, 5))
    x = np.arange(len(summary))
    counts = summary["n_galaxies"].to_numpy()
    bars = ax.bar(x, counts, color=cfg["color"], edgecolor="white", linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(summary["bin_label"], rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Number of galaxies")
    ax.set_xlabel(f"{cfg['mag_xlabel']} bin")
    ax.set_title(f"{cfg['short']} catalog vs magnitude\nN={n_pool:,}  |  {subtitle}")
    for bar, n in zip(bars, counts):
        if n > 0:
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height(),
                f"{n:,}",
                ha="center",
                va="bottom",
                fontsize=5.5,
            )
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_median_ba_vs_mag(
    summary: pd.DataFrame,
    n_pool: int,
    out_png: Path,
    *,
    cfg: dict,
    subtitle: str,
) -> None:
    plot = summary.loc[summary["n_galaxies"] > 0].copy()
    centers = []
    for _, r in plot.iterrows():
        if not np.isfinite(r["mag_lo"]):
            centers.append(r["mag_hi"] - 0.25)
        else:
            centers.append(0.5 * (r["mag_lo"] + r["mag_hi"]))
    centers = np.asarray(centers, dtype=float)

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(
        centers,
        plot["median_ba"],
        "o-",
        color=cfg["color"],
        markersize=5,
        linewidth=1.5,
        label=f"median b/a (N={n_pool:,})",
    )
    ax.set_xlabel(f"{cfg['mag_xlabel']} (bin center)")
    ax.set_ylabel(f"median {cfg['ba_ylabel']}")
    ax.set_title(
        f"{cfg['short']}: median axis ratio vs magnitude\n"
        f"{cfg['ba_note']}  |  {subtitle}"
    )
    ax.set_ylim(0.0, 1.0)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_ba_vs_mag_scatter(
    mag: np.ndarray,
    ba: np.ndarray,
    n_pool: int,
    out_png: Path,
    *,
    cfg: dict,
    subtitle: str,
    overlay_n: int,
    mag_max: float,
    seed: int = 42,
) -> None:
    in_range = np.isfinite(mag) & (mag <= mag_max)
    mag_p = mag[in_range]
    ba_p = ba[in_range]

    fig, ax = plt.subplots(figsize=(10, 6))
    hb = ax.hexbin(
        mag_p,
        ba_p,
        gridsize=80,
        bins="log",
        cmap="viridis",
        mincnt=1,
        extent=(
            float(np.nanpercentile(mag_p, 0.1)),
            float(mag_max),
            0.0,
            1.0,
        ),
    )
    cb = fig.colorbar(hb, ax=ax, pad=0.02)
    cb.set_label(r"log$_{10}$(N) per hex")

    rng = np.random.default_rng(seed)
    if len(mag_p) > overlay_n:
        idx = rng.choice(len(mag_p), size=overlay_n, replace=False)
        ax.scatter(
            mag_p[idx],
            ba_p[idx],
            s=2,
            c="white",
            alpha=0.08,
            linewidths=0,
            rasterized=True,
            zorder=2,
        )

    ax.set_xlabel(cfg["mag_xlabel"])
    ax.set_ylabel(cfg["ba_ylabel"])
    ax.set_ylim(0.0, 1.0)
    ax.set_xlim(right=mag_max)
    ax.set_title(
        f"{cfg['short']}: $b/a$ vs magnitude (hexbin + {overlay_n:,} overlay pts)\n"
        f"N={n_pool:,} shown for mag ≤ {mag_max:g}  |  {subtitle}"
    )
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_readme(
    out_root: Path,
    *,
    cfg: dict,
    catalog_rel: str,
    subtitle: str,
    n_pool: int,
    ra_span: tuple[float, float],
    dec_span: tuple[float, float],
    median_ba: float,
    mean_ba: float,
    mag_med: float,
    regenerate_cmd: str,
) -> None:
    text = f"""# {cfg['label']} — mag / b/a audit

Catalog: `{catalog_rel}`

| Metric | Value |
|--------|-------|
| Rows (finite mag + b/a) | {n_pool:,} |
| RA span | {ra_span[0]:.3f} … {ra_span[1]:.3f} |
| Dec span | {dec_span[0]:.3f} … {dec_span[1]:.3f} |
| median mag | {mag_med:.3f} |
| median b/a | {median_ba:.4f} |
| mean b/a | {mean_ba:.4f} |

Selection: {subtitle}

## Axis ratio

{cfg['ba_note']}

## Plots

| File | Content |
|------|---------|
| `mag_histogram_ba.png` | Galaxy counts per 0.5 mag bin (bright open bin `<15`; display through mag ≤ 26) |
| `median_ba_vs_mag.png` | Median b/a per mag bin (same display cut) |
| `ba_vs_mag_scatter.png` | b/a vs mag (log hexbin + thin scatter overlay; mag ≤ 26) |

CSV companions retain the full magnitude range: `mag_histogram_ba.csv`, `median_ba_vs_mag.csv`.

## Regenerate

```bash
{regenerate_cmd}
```
"""
    (out_root / "README.md").write_text(text, encoding="utf-8")


def run_one(
    survey: str,
    *,
    exp_analogue: bool,
    mag_clip: float,
    bin_width: float,
    overlay_n: int,
) -> None:
    cfg = SURVEYS[survey]
    csv_path = cfg["exp_csv"] if exp_analogue else cfg["sample_csv"]
    out_root = cfg["out_exp"] if exp_analogue else cfg["out_sample"]
    if exp_analogue:
        subtitle = (
            f"EXP analogue  {EXP_ANALOGUE_N_MIN}<n<{EXP_ANALOGUE_N_MAX}  |  no color / q cuts"
        )
        regen = f"python scripts/audit_des_hsc_mag_ba.py --survey {survey} --exp-analogue"
    else:
        subtitle = "full sample  |  no color / q / n cuts"
        regen = f"python scripts/audit_des_hsc_mag_ba.py --survey {survey}"

    pool = load_pool(csv_path, cfg)
    mag = pool["mag"].to_numpy(dtype=float)
    ba = pool["ba"].to_numpy(dtype=float)
    n_pool = len(pool)

    edges = mag_bin_edges(mag, clip=mag_clip, width=bin_width)
    summary = make_bin_summary(mag, ba, edges)

    out_root.mkdir(parents=True, exist_ok=True)
    summary.to_csv(out_root / "mag_histogram_ba.csv", index=False)
    summary.to_csv(out_root / "median_ba_vs_mag.csv", index=False)

    display = plot_display_summary(summary, mag_max=PLOT_MAG_MAX, min_n=MIN_BIN_N)
    plot_mag_histogram(
        display, n_pool, out_root / "mag_histogram_ba.png", cfg=cfg, subtitle=subtitle
    )
    plot_median_ba_vs_mag(
        display, n_pool, out_root / "median_ba_vs_mag.png", cfg=cfg, subtitle=subtitle
    )
    plot_ba_vs_mag_scatter(
        mag,
        ba,
        n_pool,
        out_root / "ba_vs_mag_scatter.png",
        cfg=cfg,
        subtitle=subtitle,
        overlay_n=overlay_n,
        mag_max=PLOT_MAG_MAX,
    )

    try:
        catalog_rel = str(csv_path.relative_to(REPO_ROOT)).replace("\\", "/")
    except ValueError:
        catalog_rel = str(csv_path)

    write_readme(
        out_root,
        cfg=cfg,
        catalog_rel=catalog_rel,
        subtitle=subtitle,
        n_pool=n_pool,
        ra_span=(float(pool["RA_ICRS"].min()), float(pool["RA_ICRS"].max())),
        dec_span=(float(pool["DE_ICRS"].min()), float(pool["DE_ICRS"].max())),
        median_ba=round(float(np.median(ba)), 4),
        mean_ba=round(float(np.mean(ba)), 4),
        mag_med=round(float(np.median(mag)), 3),
        regenerate_cmd=regen,
    )
    print(f"[*] Wrote {out_root}", flush=True)
    print(summary.to_string(index=False), flush=True)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--survey", choices=("des", "hsc", "both"), default="both")
    p.add_argument(
        "--exp-analogue",
        action="store_true",
        help=f"Use EXP-analogue catalogs ({EXP_ANALOGUE_N_MIN}<n<{EXP_ANALOGUE_N_MAX})",
    )
    p.add_argument("--mag-clip", type=float, default=MAG_CLIP)
    p.add_argument("--bin-width", type=float, default=BIN_WIDTH)
    p.add_argument("--overlay-n", type=int, default=SCATTER_OVERLAY_N)
    args = p.parse_args()

    surveys = ("des", "hsc") if args.survey == "both" else (args.survey,)
    for s in surveys:
        run_one(
            s,
            exp_analogue=args.exp_analogue,
            mag_clip=args.mag_clip,
            bin_width=args.bin_width,
            overlay_n=args.overlay_n,
        )


if __name__ == "__main__":
    main()
