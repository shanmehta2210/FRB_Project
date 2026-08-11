"""
Decisive test: do Jimin's lnLStar / score / Photoz cuts bias cos(i) toward edge-on
(i.e. remove face-on/round galaxies), explaining why production v2 (mag<=21) looks
MORE face-on than Jimin?

Fetch the shared field pool ONCE with all diagnostic columns, then slice in pandas.
Field: RA[148,152] x Dec[0,4]. LEFT JOIN Photoz so we can measure the nnAvgZ cut too.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astroquery.sdss import SDSS

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import Q0, hubble_cosi_from_ba  # noqa: E402
from pipeline_null_plot_utils import REPO_ROOT  # noqa: E402

OUT = REPO_ROOT / "plots" / "plots_null" / "v2" / "sdss_audit" / "Jimin" / "cut_bias"
CACHE = OUT / "field_pool_raw.csv"
DR = 16

FETCH_SQL = """
SELECT
  p.objid, p.ra, p.dec,
  p.r AS modelmag_r, p.u AS modelmag_u,
  p.expAB_r, p.deVAB_r,
  p.lnLExp_r, p.lnLDeV_r, p.lnLStar_r,
  p.fracDeV_r, p.score,
  pz.nnAvgZ AS nnavgz
FROM PhotoObj AS p
LEFT OUTER JOIN Photoz AS pz ON pz.objid = p.objid
WHERE p.ra BETWEEN 148.0 AND 152.0
  AND p.dec BETWEEN 0.0 AND 4.0
  AND p.mode = 1 AND p.clean = 1 AND p.type_r = 3
  AND p.r BETWEEN 12 AND 21
""".strip()


def fetch() -> pd.DataFrame:
    if CACHE.exists():
        print(f"[*] Loading cache {CACHE}")
        return pd.read_csv(CACHE)
    print("[*] Fetching field pool from DR16 (may take a minute) ...")
    tbl = SDSS.query_sql(FETCH_SQL, data_release=DR, timeout=1200)
    df = tbl.to_pandas()
    df.columns = [c.lower() for c in df.columns]
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(CACHE, index=False)
    return df


def cosi(ba: np.ndarray, q0: float = Q0) -> np.ndarray:
    return np.array([hubble_cosi_from_ba(float(v), q0=q0) for v in ba], dtype=float)


def strict_cosi(ba_series: pd.Series, q0: float = Q0) -> np.ndarray:
    v = pd.to_numeric(ba_series, errors="coerce").to_numpy(dtype=float)
    ok = np.isfinite(v) & (v > q0) & (v <= 1.0)
    return cosi(v[ok], q0=q0)


def ecdf(vals: np.ndarray, x: np.ndarray) -> np.ndarray:
    s = np.sort(vals[np.isfinite(vals)])
    return np.searchsorted(s, x, side="right") / max(1, len(s))


def summ(ba_series: pd.Series) -> tuple[int, float]:
    c = strict_cosi(ba_series)
    return len(c), (float(np.median(c)) if len(c) else float("nan"))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = fetch()
    n0 = len(df)
    print(f"[*] Field pool (type_r=3, mode, clean, r 12-21): N={n0:,}")

    exp = pd.to_numeric(df["lnlexp_r"], errors="coerce")
    dev = pd.to_numeric(df["lnldev_r"], errors="coerce")
    star = pd.to_numeric(df["lnlstar_r"], errors="coerce")
    score = pd.to_numeric(df["score"], errors="coerce")
    nnz = pd.to_numeric(df["nnavgz"], errors="coerce")
    ur = pd.to_numeric(df["modelmag_u"], errors="coerce") - pd.to_numeric(
        df["modelmag_r"], errors="coerce"
    )
    frac = pd.to_numeric(df["fracdev_r"], errors="coerce")

    lnl_exp_wins = exp > dev

    # === Decisive: for each Jimin cut, cos(i) of KEPT vs REMOVED, on the lnL-exp base ===
    base = df[lnl_exp_wins].copy()
    b_star = pd.to_numeric(base["lnlstar_r"], errors="coerce")
    b_score = pd.to_numeric(base["score"], errors="coerce")
    b_nnz = pd.to_numeric(base["nnavgz"], errors="coerce")
    b_frac = pd.to_numeric(base["fracdev_r"], errors="coerce")

    rows = []

    def add(name: str, mask_keep: pd.Series, frame: pd.DataFrame) -> None:
        keep = frame[mask_keep]
        rem = frame[~mask_keep]
        nk, mk = summ(keep["expab_r"])
        nr, mr = summ(rem["expab_r"])
        rows.append(
            {
                "cut": name,
                "n_kept": nk,
                "median_cosi_kept": round(mk, 4),
                "n_removed": nr,
                "median_cosi_removed": round(mr, 4) if nr else np.nan,
                "shift_removed_minus_kept": round(mr - mk, 4) if nr else np.nan,
            }
        )

    add("lnLStar_r < -10", b_star < -10, base)
    add("score > 0.8", b_score > 0.8, base)
    add("nnAvgZ > 0 (has photoz)", b_nnz > 0, base)
    add("fracDeV_r = 0", b_frac == 0, base)
    add("u-r < 2.3 (ours)", (pd.to_numeric(base["modelmag_u"], errors="coerce")
                             - pd.to_numeric(base["modelmag_r"], errors="coerce")) < 2.3, base)
    bias = pd.DataFrame(rows)
    bias.to_csv(OUT / "cut_bias_kept_vs_removed.csv", index=False)
    print("\n=== cos(i): kept vs REMOVED (base = lnLExp wins) ===")
    print(bias.to_string(index=False))

    # === End-to-end CDFs: Production (mag<=21) vs Jimin, same box ===
    # Production strict: lnL exp wins + u-r<2.3 + expAB>0.2
    prod_mask = lnl_exp_wins & (ur < 2.3)
    prod = df[prod_mask]
    # Jimin V2 strict: lnLStar<-10 + nnz>0 + score>0.8 + lnL + expAB>0.2
    jv2_mask = lnl_exp_wins & (star < -10) & (nnz > 0) & (score > 0.8)
    jv2 = df[jv2_mask]
    jv1_mask = jv2_mask & (frac == 0)
    jv1 = df[jv1_mask]
    # "pool" = just galaxies+mag+lnL (no purity, no color)
    pool = df[lnl_exp_wins]

    series = [
        ("Pool: gal+mag+lnL (no extra cuts)", strict_cosi(pool["expab_r"]), "#999999"),
        ("Production: +u-r<2.3", strict_cosi(prod["expab_r"]), "#4daf4a"),
        ("Jimin V2: +lnLStar+photoz+score", strict_cosi(jv2["expab_r"]), "#377eb8"),
        ("Jimin V1: +fracDeV=0", strict_cosi(jv1["expab_r"]), "#e41a1c"),
    ]
    x = np.linspace(0, 1, 401)
    fig, ax = plt.subplots(figsize=(7.6, 7.6))
    ax.plot((0, 1), (0, 1), "k--", lw=1.1, label="Uniform")
    cdf_rows = []
    for lab, c, col in series:
        med = float(np.median(c)) if len(c) else float("nan")
        ax.plot(x, ecdf(c, x), color=col, lw=2.0, label=f"{lab}\n  N={len(c):,}, med={med:.3f}")
        cdf_rows.append({"selection": lab, "n_strict": len(c), "median_cosi": round(med, 4)})
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(rf"$\cos(i)$ (Hubble, $q_0={Q0:g}$)")
    ax.set_ylabel("Cumulative distribution")
    ax.set_title("Same field (RA148-152, Dec0-4), model r<=21\nProduction vs Jimin selections (expAB_r)")
    ax.legend(loc="lower right", fontsize=7.5)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "cdf_prod_vs_jimin_samefield.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(cdf_rows).to_csv(OUT / "cdf_prod_vs_jimin_summary.csv", index=False)
    print("\n=== End-to-end CDFs (same field) ===")
    print(pd.DataFrame(cdf_rows).to_string(index=False))

    # === Incremental: start production-like, add Jimin cuts one by one ===
    steps = [
        ("gal+mag+lnL", lnl_exp_wins),
        ("+lnLStar<-10", lnl_exp_wins & (star < -10)),
        ("+nnAvgZ>0", lnl_exp_wins & (star < -10) & (nnz > 0)),
        ("+score>0.8 (=Jimin V2)", lnl_exp_wins & (star < -10) & (nnz > 0) & (score > 0.8)),
        ("+fracDeV=0 (=Jimin V1)", lnl_exp_wins & (star < -10) & (nnz > 0) & (score > 0.8) & (frac == 0)),
    ]
    inc_rows = []
    prev = None
    for name, m in steps:
        c = strict_cosi(df[m]["expab_r"])
        med = float(np.median(c)) if len(c) else float("nan")
        inc_rows.append(
            {
                "stage": name,
                "n_strict": len(c),
                "median_cosi": round(med, 4),
                "delta_median_vs_prev": (round(med - prev, 4) if prev is not None else 0.0),
            }
        )
        prev = med
    pd.DataFrame(inc_rows).to_csv(OUT / "incremental_median_shift.csv", index=False)
    print("\n=== Incremental median cos(i) as Jimin cuts stack ===")
    print(pd.DataFrame(inc_rows).to_string(index=False))

    # === lnLStar vs roundness: is star-like-ness correlated with face-on/round? ===
    g = df[lnl_exp_wins].copy()
    g_ba = pd.to_numeric(g["expab_r"], errors="coerce")
    g_star = pd.to_numeric(g["lnlstar_r"], errors="coerce")
    # bin by expAB, show fraction failing lnLStar<-10 and median lnLStar
    bins = np.linspace(0.05, 1.0, 20)
    g["ba_bin"] = pd.cut(g_ba, bins)
    corr_rows = []
    for b, sub in g.groupby("ba_bin", observed=True):
        st = pd.to_numeric(sub["lnlstar_r"], errors="coerce")
        sc = pd.to_numeric(sub["score"], errors="coerce")
        corr_rows.append(
            {
                "ba_bin": str(b),
                "n": len(sub),
                "frac_fail_lnLStar_ge_-10": round(float((st >= -10).mean()), 4),
                "median_lnLStar": round(float(st.median()), 2),
                "frac_fail_score_le_0.8": round(float((sc <= 0.8).mean()), 4),
            }
        )
    pd.DataFrame(corr_rows).to_csv(OUT / "roundness_vs_cut_failure.csv", index=False)
    print("\n=== Do rounder (face-on) galaxies fail lnLStar/score more? ===")
    print(pd.DataFrame(corr_rows).to_string(index=False))

    print(f"\n[*] Wrote outputs -> {OUT}")


if __name__ == "__main__":
    main()
