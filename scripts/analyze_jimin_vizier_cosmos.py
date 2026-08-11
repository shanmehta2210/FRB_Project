"""
Cross-check VizieR SDSS_DR16_cosmos.txt vs Jimin SkyServer catalogs + Hubble CDFs.

VizieR file uses rdVell = deVAB_r (not expAB_r). Outputs under
plots/plots_null/v2/sdss_audit/Jimin/vizier_cosmos/.
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

from null_catalog_utils import Q0, hubble_cosi_from_ba  # noqa: E402
from pipeline_null_plot_utils import REPO_ROOT  # noqa: E402

JIMIN = REPO_ROOT / "plots" / "plots_null" / "v2" / "sdss_audit" / "Jimin"
TXT = JIMIN / "catalog" / "SDSS_DR16_cosmos.txt"
CSV = JIMIN / "catalog" / "SDSS_DR16_cosmos.csv"
OUT = JIMIN / "vizier_cosmos"
RA0, DEC0 = 150.1255, 2.2108
MATCH_TOL_DEG = 1.5 / 3600.0


def parse_vizier_txt(path: Path) -> pd.DataFrame:
    cols = None
    rows: list[list[str]] = []
    with path.open(encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            if line.startswith("RA_ICRS"):
                cols = line.split()
                continue
            if cols is None:
                continue
            parts = line.split()
            if len(parts) < len(cols):
                continue
            rows.append(parts[: len(cols)])
    df = pd.DataFrame(rows, columns=cols)
    for c in cols:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    return df


def empirical_cdf(vals: np.ndarray, x: np.ndarray) -> np.ndarray:
    s = np.sort(vals[np.isfinite(vals)])
    return np.searchsorted(s, x, side="right") / max(1, len(s))


def cosi_from_ba(ba: np.ndarray, q0: float = Q0) -> np.ndarray:
    return np.array([hubble_cosi_from_ba(float(v), q0=q0) for v in ba], dtype=float)


def strict_ba(ba: pd.Series, q0: float = Q0) -> np.ndarray:
    v = pd.to_numeric(ba, errors="coerce").to_numpy(dtype=float)
    ok = np.isfinite(v) & (v > q0) & (v <= 1.0)
    return v[ok]


def in_cone(ra: np.ndarray, dec: np.ndarray, ra0: float, dec0: float, r_deg: float = 1.0) -> np.ndarray:
    dra = (ra - ra0) * np.cos(np.deg2rad(dec0))
    ddec = dec - dec0
    return np.sqrt(dra**2 + ddec**2) <= r_deg


def match_nearest(a_ra, a_dec, b_ra, b_dec, tol: float = MATCH_TOL_DEG):
    order = np.argsort(b_ra)
    bra_s, bdec_s = b_ra[order], b_dec[order]
    matched = np.zeros(len(a_ra), dtype=bool)
    idx_b = np.full(len(a_ra), -1, dtype=int)
    for i, (ra, dec) in enumerate(zip(a_ra, a_dec)):
        lo = np.searchsorted(bra_s, ra - 2 * tol, side="left")
        hi = np.searchsorted(bra_s, ra + 2 * tol, side="right")
        if hi <= lo:
            continue
        d2 = (bra_s[lo:hi] - ra) ** 2 + (bdec_s[lo:hi] - dec) ** 2
        j = int(np.argmin(d2))
        if d2[j] <= tol**2:
            matched[i] = True
            idx_b[i] = int(order[lo + j])
    return matched, idx_b


def plot_cdf_series(
    series: list[tuple[str, np.ndarray, str]],
    *,
    title: str,
    out: Path,
) -> None:
    x = np.linspace(0, 1, 401)
    fig, ax = plt.subplots(figsize=(7.2, 7.2))
    ax.plot((0, 1), (0, 1), "k--", lw=1.2, label="Uniform")
    for lab, cosi, col in series:
        if len(cosi) == 0:
            continue
        med = float(np.median(cosi))
        ax.plot(
            x,
            empirical_cdf(cosi, x),
            color=col,
            lw=2.0,
            label=f"{lab}  N={len(cosi):,}, med={med:.3f}",
        )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(rf"$\cos(i)$ (Hubble, $q_0={Q0:g}$)")
    ax.set_ylabel("Cumulative distribution")
    ax.set_title(title)
    ax.legend(loc="upper left", fontsize=7.5)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    viz = parse_vizier_txt(TXT)
    viz.to_csv(CSV, index=False)

    v1 = pd.read_csv(JIMIN / "catalog" / "v1_fracDev0_and_lnL.csv")
    v2 = pd.read_csv(JIMIN / "catalog" / "v2_lnL_exp.csv")
    v1.columns = [c.lower() for c in v1.columns]
    v2.columns = [c.lower() for c in v2.columns]
    ecol = "expab_r" if "expab_r" in v2.columns else "expAB_r"
    dcol = "devab_r" if "devab_r" in v2.columns else "deVAB_r"

    # Cross-match stats
    m_v1, _ = match_nearest(
        viz["RA_ICRS"].to_numpy(),
        viz["DE_ICRS"].to_numpy(),
        v1["ra"].to_numpy(),
        v1["dec"].to_numpy(),
    )
    m_v2, i_v2 = match_nearest(
        viz["RA_ICRS"].to_numpy(),
        viz["DE_ICRS"].to_numpy(),
        v2["ra"].to_numpy(),
        v2["dec"].to_numpy(),
    )
    v1_cone = v1[in_cone(v1["ra"].to_numpy(), v1["dec"].to_numpy(), RA0, DEC0)]
    v2_cone = v2[in_cone(v2["ra"].to_numpy(), v2["dec"].to_numpy(), RA0, DEC0)]
    mv1_rev, _ = match_nearest(
        v1_cone["ra"].to_numpy(),
        v1_cone["dec"].to_numpy(),
        viz["RA_ICRS"].to_numpy(),
        viz["DE_ICRS"].to_numpy(),
    )
    mv2_rev, _ = match_nearest(
        v2_cone["ra"].to_numpy(),
        v2_cone["dec"].to_numpy(),
        viz["RA_ICRS"].to_numpy(),
        viz["DE_ICRS"].to_numpy(),
    )

    # Matched ba comparison
    sub = viz.loc[m_v2].copy()
    sub["expAB_r"] = v2.iloc[i_v2[m_v2]][ecol].to_numpy()
    sub["deVAB_jimin"] = v2.iloc[i_v2[m_v2]][dcol].to_numpy()
    ba_diff_dev = np.nanmedian(np.abs(sub["rdVell"] - sub["deVAB_jimin"]))
    ba_diff_exp = np.nanmedian(np.abs(sub["rdVell"] - sub["expAB_r"]))

    rows = [
        {"metric": "vizier_N", "value": len(viz)},
        {"metric": "vizier_area_deg2_approx", "value": round(np.pi, 4)},
        {"metric": "vizier_density_per_deg2", "value": round(len(viz) / np.pi, 2)},
        {"metric": "vizier_rmag_le_21", "value": int((viz["rmag"] <= 21).sum())},
        {"metric": "vizier_rPmag_le_21", "value": int((viz["rPmag"] <= 21).sum())},
        {"metric": "jimin_v1_N", "value": len(v1)},
        {"metric": "jimin_v2_N", "value": len(v2)},
        {"metric": "jimin_v1_in_1deg_cone", "value": len(v1_cone)},
        {"metric": "jimin_v2_in_1deg_cone", "value": len(v2_cone)},
        {"metric": "vizier_matched_v1_1p5asec", "value": int(m_v1.sum())},
        {"metric": "vizier_matched_v2_1p5asec", "value": int(m_v2.sum())},
        {"metric": "v1_cone_found_in_vizier_pct", "value": round(100 * mv1_rev.mean(), 2)},
        {"metric": "v2_cone_found_in_vizier_pct", "value": round(100 * mv2_rev.mean(), 2)},
        {"metric": "median_abs_rdVell_minus_deVAB", "value": round(float(ba_diff_dev), 6)},
        {"metric": "median_abs_rdVell_minus_expAB", "value": round(float(ba_diff_exp), 6)},
    ]
    pd.DataFrame(rows).to_csv(OUT / "crossmatch_summary.csv", index=False)

    # Hubble CDFs from VizieR rdVell (= deVAB_r)
    pools = {
        "all_rPmag_lt22": viz,
        "rPmag_le_21": viz[viz["rPmag"] <= 21],
        "rmag_le_21": viz[viz["rmag"] <= 21],
        "rmag_12_21": viz[(viz["rmag"] >= 12) & (viz["rmag"] <= 21)],
    }
    cdf_rows = []
    series_overlay = []
    colors = {
        "all_rPmag_lt22": "#984ea3",
        "rPmag_le_21": "#ff7f00",
        "rmag_le_21": "#4daf4a",
        "rmag_12_21": "#a65628",
    }
    for key, frame in pools.items():
        ba = strict_ba(frame["rdVell"])
        cosi = cosi_from_ba(ba)
        med = float(np.median(cosi)) if len(cosi) else float("nan")
        cdf_rows.append(
            {
                "sample": f"vizier_{key}_deVAB",
                "ba_column": "rdVell (=deVAB_r)",
                "n_sql": len(frame),
                "n_strict": len(cosi),
                "median_cosi": round(med, 4) if len(cosi) else np.nan,
            }
        )
        plot_cdf_series(
            [(f"VizieR {key} (deVAB)", cosi, colors[key])],
            title=(
                f"VizieR SDSS DR16 COSMOS 1° — Hubble from rdVell (deVAB)\n"
                f"{key} | ba > {Q0:g}"
            ),
            out=OUT / "plots" / f"cdf_{key}_deVAB.png",
        )
        series_overlay.append((f"viz {key} deVAB", cosi, colors[key]))

    # Jimin cone-restricted for fairer overlay
    ba_v1 = strict_ba(v1_cone[ecol])
    ba_v2 = strict_ba(v2_cone[ecol])
    cosi_v1 = cosi_from_ba(ba_v1)
    cosi_v2 = cosi_from_ba(ba_v2)
    cdf_rows.append(
        {
            "sample": "jimin_v1_in_cone_expAB",
            "ba_column": "expAB_r",
            "n_sql": len(v1_cone),
            "n_strict": len(cosi_v1),
            "median_cosi": round(float(np.median(cosi_v1)), 4),
        }
    )
    cdf_rows.append(
        {
            "sample": "jimin_v2_in_cone_expAB",
            "ba_column": "expAB_r",
            "n_sql": len(v2_cone),
            "n_strict": len(cosi_v2),
            "median_cosi": round(float(np.median(cosi_v2)), 4),
        }
    )

    # Matched objects: same galaxies, deVAB vs expAB
    ba_dev_m = strict_ba(sub["rdVell"])
    ba_exp_m = strict_ba(sub["expAB_r"])
    cosi_dev_m = cosi_from_ba(ba_dev_m)
    cosi_exp_m = cosi_from_ba(ba_exp_m)
    cdf_rows.append(
        {
            "sample": "matched_v2_deVAB",
            "ba_column": "rdVell",
            "n_sql": len(sub),
            "n_strict": len(cosi_dev_m),
            "median_cosi": round(float(np.median(cosi_dev_m)), 4),
        }
    )
    cdf_rows.append(
        {
            "sample": "matched_v2_expAB",
            "ba_column": "expAB_r",
            "n_sql": len(sub),
            "n_strict": len(cosi_exp_m),
            "median_cosi": round(float(np.median(cosi_exp_m)), 4),
        }
    )

    plot_cdf_series(
        [
            ("VizieR rmag≤21 (deVAB)", cosi_from_ba(strict_ba(pools["rmag_le_21"]["rdVell"])), "#4daf4a"),
            ("Jimin V2 in cone (expAB)", cosi_v2, "#377eb8"),
            ("Jimin V1 in cone (expAB)", cosi_v1, "#e41a1c"),
        ],
        title=(
            "VizieR COSMOS vs Jimin (1° cone)\n"
            "WARNING: VizieR ba = deVAB_r; Jimin ba = expAB_r"
        ),
        out=OUT / "plots" / "cdf_overlay_vizier_vs_jimin.png",
    )
    plot_cdf_series(
        [
            ("Same objs: deVAB (rdVell)", cosi_dev_m, "#984ea3"),
            ("Same objs: expAB", cosi_exp_m, "#377eb8"),
        ],
        title=(
            "Matched VizieR ∩ Jimin V2 — same galaxies\n"
            "Hubble cos(i): deVAB vs expAB"
        ),
        out=OUT / "plots" / "cdf_matched_deVAB_vs_expAB.png",
    )
    plot_cdf_series(
        series_overlay,
        title="VizieR COSMOS DR16 — Hubble from deVAB (rdVell), mag variants",
        out=OUT / "plots" / "cdf_overlay_vizier_mag_variants.png",
    )

    pd.DataFrame(cdf_rows).to_csv(OUT / "cdf_summary.csv", index=False)
    print("[*] Wrote", OUT)
    for r in rows:
        print(f"  {r['metric']}: {r['value']}")


if __name__ == "__main__":
    main()
