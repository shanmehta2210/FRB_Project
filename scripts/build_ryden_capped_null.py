"""
Capped-Ryden null: bake the REX b/a<=0.8 ceiling into the intrinsic shape law.

Idea (user proposal): instead of the ad-hoc ``scaled`` renormalisation cos i -> cos i/H(0.8),
build a *genuine* Ryden/elliptical-disk distribution whose intrinsic FACE-ON b/a is
uniform below 0.8 and falls off exponentially above it (scale ``lam``). Since a disk's
projected b/a never exceeds its face-on b/a, this enforces the observed b/a<=~0.8 cap
physically, while keeping a real cos(i) spread (unlike the degenerate delta that
``scaled`` corresponds to; see SCALED_IS_DEGENERATE_RYDEN.md).

Pipeline:
  1. Fit (mu_g, sig_g, lam) of CappedShapeParams (cap=0.8 fixed) to the LS EXP b/a
     histogram via differential_evolution + chi2.
  2. Sample P(cos i | b/a) from the fitted capped model for LS at each mag cut.
  3. Overlay the resulting cos(i) CDF on ``scaled`` (Hubble/0.8) to see how close it is.

Outputs under plots/plots_null/v2/ls_audit/scaled_ryden_capped/.

Run from repo root::

    python scripts/build_ryden_capped_null.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from elliptical_disk_model import (  # noqa: E402
    CappedShapeParams,
    ConditionalCosiSampler,
    chi2_hist,
    generate_q_capped,
    histogram_q,
    sample_faceon_ba_capped,
)
from fit_ls_scaled_elliptical import load_ls  # noqa: E402
from null_catalog_utils import (  # noqa: E402
    LS_CATALOG_V2_EXP_DEFAULT,
    Q0,
    hubble_cosi_from_ba,
)
from pipeline_null_plot_utils import REPO_ROOT  # noqa: E402

OUT_DIR = REPO_ROOT / "plots" / "plots_null" / "v2" / "ls_audit" / "scaled_ryden_capped"
MAG_CUTS = (20.0, 21.0, 22.0)
BA_FACE_CAP = 0.8
BIN_WIDTH = 0.05
N_MODEL_FIT = 120_000
N_MODEL = 400_000
N_POOL = 3_000_000
H08 = float(hubble_cosi_from_ba(BA_FACE_CAP, q0=Q0))


def scaled_cosi(q: np.ndarray) -> np.ndarray:
    val = (q * q - Q0 * Q0) / (1.0 - Q0 * Q0)
    return np.clip(np.sqrt(np.clip(val, 0.0, 1.0)) / H08, 0.0, 1.0)


def ecdf(v: np.ndarray, x: np.ndarray) -> np.ndarray:
    s = np.sort(v[np.isfinite(v)])
    return np.searchsorted(s, x, side="right") / max(1, len(s))


def med(v: np.ndarray) -> float:
    v = v[np.isfinite(v)]
    return float(np.median(v)) if v.size else float("nan")


def model_hist_capped(
    params: CappedShapeParams,
    n: int,
    data_counts: np.ndarray,
    edges: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    q = generate_q_capped(params, n, rng)
    m, _ = np.histogram(q, bins=edges)
    m = m.astype(float)
    n_data, n_mod = float(data_counts.sum()), float(m.sum())
    if n_mod > 0:
        m *= n_data / n_mod
    return m


def fit_capped(data_counts: np.ndarray, edges: np.ndarray, *, seed: int = 42) -> CappedShapeParams:
    # params: mu_g (thickness), sig_g, lam (fall-off); cap fixed at 0.8
    bounds = [(0.05, 0.5), (0.01, 0.2), (1e-3, 0.5)]

    def objective(x: np.ndarray, n_draw: int = N_MODEL_FIT) -> float:
        p = CappedShapeParams(mu_g=x[0], sig_g=x[1], cap=BA_FACE_CAP, lam=x[2])
        rng = np.random.default_rng(seed + int(1e6 * abs(x[0] + x[2])) % 100000)
        return chi2_hist(data_counts, model_hist_capped(p, n_draw, data_counts, edges, rng))

    print("[*] Fitting capped model (mu_g, sig_g, lam), cap=0.8 fixed ...", flush=True)
    res = differential_evolution(
        objective, bounds=bounds, seed=seed, maxiter=20, popsize=12,
        mutation=(0.5, 1.0), recombination=0.7, polish=True, tol=0.02, atol=1.0,
        updating="immediate", workers=1,
    )
    best = CappedShapeParams(mu_g=res.x[0], sig_g=res.x[1], cap=BA_FACE_CAP, lam=res.x[2])
    chi2 = objective(res.x, n_draw=N_MODEL)
    print(f"    best chi2={chi2:.1f}  {best}", flush=True)
    best_chi2 = chi2
    globals()["_LAST_CHI2"] = best_chi2
    return best


def plot_intrinsic_law(params: CappedShapeParams, out: Path) -> None:
    rng = np.random.default_rng(0)
    f = sample_faceon_ba_capped(500_000, params.cap, params.lam, rng)
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.hist(f, bins=100, range=(0, 1), density=True, color="#7fbf7b", edgecolor="white")
    ax.axvline(params.cap, color="#d7191c", ls="--", lw=1.5, label=f"cap={params.cap:g}")
    ax.set_xlabel(r"intrinsic face-on $b/a = 1-\varepsilon$")
    ax.set_ylabel("density")
    ax.set_title(f"Capped face-on b/a law (uniform below cap, exp fall-off; lam={params.lam:.3f})")
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_hist_fit(data_counts, model_counts, centers, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    w = centers[1] - centers[0]
    ax.bar(centers, data_counts, width=w * 0.9, color="#9ecae1", edgecolor="white",
           label="LS EXP data")
    ax.plot(centers, model_counts, "o-", color="#e41a1c", ms=4, lw=1.5, label="capped model")
    ax.axvline(BA_FACE_CAP, color="k", ls=":", lw=1, label="cap=0.8")
    ax.set_xlabel(r"projected $b/a$")
    ax.set_ylabel("counts (model scaled to data N)")
    ax.set_title("LS EXP b/a vs capped-Ryden fit")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "cdfs").mkdir(exist_ok=True)

    mag, ba = load_ls(REPO_ROOT / LS_CATALOG_V2_EXP_DEFAULT)
    data_counts, edges, centers = histogram_q(ba[ba > 0], bin_width=BIN_WIDTH)

    best = fit_capped(data_counts, edges)
    chi2 = float(globals().get("_LAST_CHI2", float("nan")))

    rng = np.random.default_rng(1)
    model_counts = model_hist_capped(best, N_MODEL, data_counts, edges, rng)
    plot_hist_fit(data_counts, model_counts, centers, OUT_DIR / "ba_hist_data_vs_model.png")
    plot_intrinsic_law(best, OUT_DIR / "intrinsic_faceon_law.png")

    # conditional cos(i) sampler from the fitted capped model
    pool_rng = np.random.default_rng(2)
    qm, cm, _ = generate_q_capped(best, N_POOL, pool_rng, return_angles=True)
    sampler = ConditionalCosiSampler.from_pool(qm, cm)
    srng = np.random.default_rng(3)

    x = np.linspace(0, 1, 401)
    rows: list[dict] = []
    fig_all, axes = plt.subplots(1, 3, figsize=(16, 5.2), sharey=True)
    for ax, lim in zip(axes, MAG_CUTS):
        m_cap = (mag <= lim) & (ba > Q0)                       # capped model keeps b/a<=~0.8 naturally
        m_scaled = (mag <= lim) & (ba > Q0) & (ba <= BA_FACE_CAP)
        c_cap = sampler.sample(ba[m_cap], srng)
        c_scaled = scaled_cosi(ba[m_scaled])
        ks = float(np.max(np.abs(ecdf(c_cap, x) - ecdf(c_scaled, x))))

        fig, ax1 = plt.subplots(figsize=(6.8, 6.8))
        for a in (ax, ax1):
            a.plot((0, 1), (0, 1), "k--", lw=1.0, label="Uniform")
            a.plot(x, ecdf(c_scaled, x), color="#4daf4a", lw=3.0,
                   label=f"scaled (Hubble/0.8)  med={med(c_scaled):.3f}")
            a.plot(x, ecdf(c_cap, x), color="#e41a1c", lw=1.6, ls="--",
                   label=f"capped-Ryden  med={med(c_cap):.3f}")
            a.set_xlim(0, 1)
            a.set_ylim(0, 1)
            a.grid(True, alpha=0.3)
            a.set_xlabel(r"$\cos(i)$  (strict $b/a>q_0$)")
        ax.set_title(f"mag limit={lim:g}   max CDF diff={ks:.4f}")
        ax1.set_ylabel("Cumulative distribution")
        ax1.set_title(f"LS capped-Ryden vs scaled — mag limit={lim:g}")
        ax1.legend(loc="upper left", fontsize=8)
        fig.tight_layout()
        fig.savefig(OUT_DIR / "cdfs" / f"mag{int(lim)}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

        rows.append({
            "mag_limit": lim,
            "n": int(m_cap.sum()),
            "median_cosi_capped": round(med(c_cap), 4),
            "median_cosi_scaled": round(med(c_scaled), 4),
            "max_cdf_diff_vs_scaled": round(ks, 4),
        })
        print(f"    mag<={lim:g}: capped med={med(c_cap):.3f}  scaled med={med(c_scaled):.3f}  "
              f"max CDF diff={ks:.4f}", flush=True)

    axes[0].set_ylabel("Cumulative distribution")
    axes[0].legend(loc="upper left", fontsize=7)
    fig_all.suptitle(
        "LS EXP: capped-Ryden (uniform face-on b/a, exp fall-off past 0.8) vs ad-hoc scaled",
        fontsize=12)
    fig_all.tight_layout()
    fig_all.savefig(OUT_DIR / "cdf_compare.png", dpi=300, bbox_inches="tight")
    plt.close(fig_all)

    pd.DataFrame(rows).to_csv(OUT_DIR / "summary.csv", index=False)
    (OUT_DIR / "fit_params.json").write_text(json.dumps({
        "track": "scaled_ryden_capped",
        "model": "CappedShapeParams: face-on b/a uniform below cap, exp fall-off above",
        "params": {k: float(v) for k, v in best.to_dict().items()},
        "chi2_bahist": chi2,
        "n_model_pool": N_POOL,
        "note": "Apply the identical fitted capped sampler to FRB host b/a for comparison.",
    }, indent=2), encoding="utf-8")
    print(f"[*] Wrote outputs -> {OUT_DIR}")


if __name__ == "__main__":
    main()
