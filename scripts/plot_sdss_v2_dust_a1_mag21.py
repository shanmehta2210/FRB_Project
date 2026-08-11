"""
SDSS v2 production null + Unterborn A1 dust correction at modelMag_r <= 21.

Method (Unterborn & Ryden 2008): Δm = 1.27 (log10 q)^2
  A1: keep face-on mag m^f = m_obs - Δm(q) <= 21  (re-admits dust-faded edge-ons)

Pool: production SDSS v2 cuts (u-r < 2.3, lnLExp > lnLDeV, expAB_r > q0).
Outputs: plots/plots_null/v2/sdss_audit/sdss_cut_evolution/dust_extinction/

Run from repo root::

    python scripts/plot_sdss_v2_dust_a1_mag21.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (  # noqa: E402
    Q0,
    SDSS_CATALOG_V2_DEFAULT,
    SDSS_UR_MAX_CDF,
    UNTERBORN_DELTA_M_COEFF,
    cosi_array_from_df,
    delta_m_unterborn,
    face_on_mag,
    filter_sdss_drop_dev_winners,
    filter_sdss_ur,
)
from pipeline_null_plot_utils import REPO_ROOT  # noqa: E402

MAG_LIMIT = 21.0
OUT = (
    REPO_ROOT
    / "plots"
    / "plots_null"
    / "v2"
    / "sdss_audit"
    / "sdss_cut_evolution"
    / "dust_extinction"
)


def empirical_cdf(vals: np.ndarray, x: np.ndarray) -> np.ndarray:
    s = np.sort(vals[np.isfinite(vals)])
    return np.searchsorted(s, x, side="right") / max(1, len(s))


def load_production_base() -> pd.DataFrame:
    path = REPO_ROOT / SDSS_CATALOG_V2_DEFAULT
    print(f"[*] Loading {path.name} ...", flush=True)
    usecols = [
        "modelMag_r",
        "modelMag_u",
        "u_r",
        "expAB_r",
        "lnLExp_r",
        "lnLDeV_r",
    ]
    df = pd.read_csv(path, usecols=usecols)
    df = filter_sdss_ur(df, SDSS_UR_MAX_CDF)
    df = filter_sdss_drop_dev_winners(df)
    ba = pd.to_numeric(df["expAB_r"], errors="coerce")
    mag = pd.to_numeric(df["modelMag_r"], errors="coerce")
    ok = ba.notna() & mag.notna() & (ba > Q0) & (ba <= 1.0) & np.isfinite(mag)
    out = df.loc[ok].copy()
    out["expAB_r"] = ba.loc[ok].to_numpy()
    out["modelMag_r"] = mag.loc[ok].to_numpy()
    print(
        f"    production base (u-r<{SDSS_UR_MAX_CDF}, lnL exp, ba>{Q0:g}): N={len(out):,}",
        flush=True,
    )
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    base = load_production_base()
    mag = base["modelMag_r"].to_numpy(dtype=float)
    ba = base["expAB_r"].to_numpy(dtype=float)
    dm = delta_m_unterborn(ba)
    m_face = face_on_mag(mag, ba)

    raw_mask = mag <= MAG_LIMIT
    a1_mask = m_face <= MAG_LIMIT
    # galaxies added by A1: face-on mag passes but observed mag fails
    added_mask = a1_mask & ~raw_mask

    raw = base.loc[raw_mask]
    a1 = base.loc[a1_mask]
    added = base.loc[added_mask]

    cosi_raw = cosi_array_from_df(raw, q_col="expAB_r", q0=Q0)
    cosi_a1 = cosi_array_from_df(a1, q_col="expAB_r", q0=Q0)
    cosi_added = cosi_array_from_df(added, q_col="expAB_r", q0=Q0) if len(added) else np.array([])

    rows = [
        {
            "mode": "raw_m_obs",
            "cut": f"modelMag_r <= {MAG_LIMIT:g}",
            "n": len(cosi_raw),
            "median_cosi": round(float(np.median(cosi_raw)), 4),
            "median_ba": round(float(np.median(raw["expAB_r"])), 4),
            "median_delta_m": round(float(np.median(dm[raw_mask])), 4),
        },
        {
            "mode": "a1_m_face",
            "cut": f"m^f = m - {UNTERBORN_DELTA_M_COEFF:g}(log10 q)^2 <= {MAG_LIMIT:g}",
            "n": len(cosi_a1),
            "median_cosi": round(float(np.median(cosi_a1)), 4),
            "median_ba": round(float(np.median(a1["expAB_r"])), 4),
            "median_delta_m": round(float(np.median(dm[a1_mask])), 4),
        },
        {
            "mode": "a1_added_only",
            "cut": "m^f<=21 but m_obs>21 (dust-rescued edge-ons)",
            "n": len(cosi_added),
            "median_cosi": round(float(np.median(cosi_added)), 4) if len(cosi_added) else np.nan,
            "median_ba": round(float(np.median(added["expAB_r"])), 4) if len(added) else np.nan,
            "median_delta_m": round(float(np.median(dm[added_mask])), 4) if len(added) else np.nan,
        },
    ]
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "summary_mag21.csv", index=False)
    print(summary.to_string(index=False), flush=True)

    # CDF overlay
    x = np.linspace(0, 1, 401)
    fig, ax = plt.subplots(figsize=(7.2, 7.2))
    ax.plot((0, 1), (0, 1), "k--", lw=1.2, label="Uniform")
    ax.plot(
        x,
        empirical_cdf(cosi_raw, x),
        color="#e41a1c",
        lw=2.0,
        label=f"raw $m_{{\\rm obs}}\\leq{MAG_LIMIT:g}$  N={len(cosi_raw):,}, med={np.median(cosi_raw):.3f}",
    )
    ax.plot(
        x,
        empirical_cdf(cosi_a1, x),
        color="#377eb8",
        lw=2.0,
        label=f"A1 $m^f\\leq{MAG_LIMIT:g}$  N={len(cosi_a1):,}, med={np.median(cosi_a1):.3f}",
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(rf"$\cos(i)$ (Hubble, $q_0={Q0:g}$, expAB$_r$)")
    ax.set_ylabel("Cumulative distribution")
    ax.set_title(
        "SDSS v2 production (u−r<2.3, lnL exp) — Unterborn A1 dust\n"
        rf"$\Delta m = {UNTERBORN_DELTA_M_COEFF:g}\,(\log_{{10}} q)^2$; "
        f"adds {len(added):,} edge-ons with $m_{{\\rm obs}}>21$"
    )
    ax.legend(loc="upper left", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "cdf_raw_vs_a1_mag21.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    # ba of added vs raw
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    bins = np.linspace(Q0, 1.0, 25)
    ax.hist(raw["expAB_r"], bins=bins, density=True, histtype="step", lw=2.0, color="#e41a1c", label="raw")
    if len(added):
        ax.hist(
            added["expAB_r"],
            bins=bins,
            density=True,
            histtype="step",
            lw=2.0,
            color="#377eb8",
            label=f"A1-added (N={len(added):,})",
        )
    ax.set_xlabel(r"$\exp AB_r$")
    ax.set_ylabel("Density")
    ax.set_title("Axis ratios: raw mag≤21 vs dust-rescued (A1-added)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "ba_hist_raw_vs_a1_added.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    readme = f"""# SDSS v2 — Unterborn A1 dust correction (mag ≤ 21)

Production null pool: `u−r < {SDSS_UR_MAX_CDF}`, `lnLExp > lnLDeV`, `expAB_r > {Q0:g}`.

## Method

Unterborn & Ryden 2008 A1:

$$\\Delta m = {UNTERBORN_DELTA_M_COEFF:g}\\,(\\log_{{10}} q)^2,\\qquad m^f = m_{{\\rm obs}} - \\Delta m(q)$$

- **raw:** keep `modelMag_r ≤ {MAG_LIMIT:g}`
- **A1:** keep `m^f ≤ {MAG_LIMIT:g}` (re-admits edge-ons dust-faded past the observed limit)

## Results

| Mode | N | median cos(i) | median b/a |
|------|--:|--------------:|-----------:|
| raw m_obs≤21 | {len(cosi_raw):,} | {np.median(cosi_raw):.4f} | {np.median(raw['expAB_r']):.4f} |
| A1 m^f≤21 | {len(cosi_a1):,} | {np.median(cosi_a1):.4f} | {np.median(a1['expAB_r']):.4f} |
| A1-added only | {len(cosi_added):,} | {(float(np.median(cosi_added)) if len(cosi_added) else float('nan')):.4f} | {(float(np.median(added['expAB_r'])) if len(added) else float('nan')):.4f} |

Δ median cos(i) (A1 − raw) = **{float(np.median(cosi_a1) - np.median(cosi_raw)):+.4f}**

## Files

- `cdf_raw_vs_a1_mag21.png`
- `ba_hist_raw_vs_a1_added.png`
- `summary_mag21.csv`

```bash
python scripts/plot_sdss_v2_dust_a1_mag21.py
```
"""
    (OUT / "README.md").write_text(readme, encoding="utf-8")
    print(f"[*] Wrote -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
