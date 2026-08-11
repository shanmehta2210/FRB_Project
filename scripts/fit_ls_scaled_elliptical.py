"""
Fit Ryden / Padilla elliptical-disk models to LS EXP b/a; write scaled_ryden / scaled_padilla.

Ryden track: 4-param shape fit, then Unterborn A1 mag re-cut CDFs.
Padilla track: 5-param shape+E0 fit with psi(theta) ~ 10^{-0.4 E}.

Run from repo root::

    python scripts/fit_ls_scaled_elliptical.py --mode both
"""

from __future__ import annotations

import argparse
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
    PADILLA_SEED,
    RYDEN_SEED,
    ConditionalCosiSampler,
    ShapeParams,
    chi2_hist,
    histogram_q,
    inclination_vs_q,
    model_histogram,
)
from null_catalog_utils import (  # noqa: E402
    LS_CATALOG_V2_EXP_DEFAULT,
    Q0,
    face_on_mag,
    hubble_cosi_from_ba,
)
from pipeline_null_plot_utils import REPO_ROOT  # noqa: E402

OUT_ROOT = REPO_ROOT / "plots" / "plots_null" / "v2" / "ls_audit"
MAG_CUTS = (20.0, 21.0, 22.0)
BIN_WIDTH = 0.05
N_MODEL = 200_000
N_MODEL_FIT = 80_000  # fewer draws inside optimizer
N_MODEL_COSI = 2_000_000  # pool for conditional P(cos i | q) sampling
RNG_SEED = 42
BA_FACE_CAP = 0.8  # ad-hoc scaled reference only


def q_from_e1e2(e1: np.ndarray, e2: np.ndarray) -> np.ndarray:
    eabs = np.hypot(np.asarray(e1, dtype=float), np.asarray(e2, dtype=float))
    q = np.full(eabs.shape, np.nan, dtype=float)
    good = np.isfinite(eabs) & (eabs < 1.0)
    q[good] = (1.0 - eabs[good]) / (1.0 + eabs[good])
    return q


def load_ls(path: Path) -> tuple[np.ndarray, np.ndarray]:
    print(f"[*] Loading {path.name} ...", flush=True)
    df = pd.read_csv(path, usecols=["modelMag_r", "shape_e1", "shape_e2"])
    mag = pd.to_numeric(df["modelMag_r"], errors="coerce").to_numpy(dtype=float)
    ba = q_from_e1e2(
        pd.to_numeric(df["shape_e1"], errors="coerce").to_numpy(dtype=float),
        pd.to_numeric(df["shape_e2"], errors="coerce").to_numpy(dtype=float),
    )
    ok = np.isfinite(mag) & np.isfinite(ba) & (ba > 0.0) & (ba <= 1.0)
    print(f"    finite mag+ba: N={ok.sum():,}", flush=True)
    return mag[ok], ba[ok]


def cosi_hubble(ba: np.ndarray) -> np.ndarray:
    ba = np.asarray(ba, dtype=float)
    val = (ba * ba - Q0 * Q0) / (1.0 - Q0 * Q0)
    out = np.zeros_like(val)
    ok = np.isfinite(val) & (val >= 0.0)
    out[ok] = np.sqrt(np.clip(val[ok], 0.0, 1.0))
    return out


def params_jsonable(p: ShapeParams) -> dict:
    d = p.to_dict()
    return {k: float(v) for k, v in d.items()}


def empirical_cdf(vals: np.ndarray, x: np.ndarray) -> np.ndarray:
    s = np.sort(vals)
    return np.searchsorted(s, x, side="right") / max(1, len(s))


def fit_params(
    data_counts: np.ndarray,
    edges: np.ndarray,
    *,
    with_E0: bool,
    n_model: int,
    seed: int,
) -> tuple[ShapeParams, float]:
    if with_E0:
        # mu_g, sig_g, mu, sig, E0
        bounds = [
            (0.05, 0.55),
            (0.01, 0.25),
            (-3.5, -0.5),
            (0.3, 1.5),
            (0.0, 2.0),
        ]
        x0_seed = [
            PADILLA_SEED.mu_g,
            PADILLA_SEED.sig_g,
            PADILLA_SEED.mu,
            PADILLA_SEED.sig,
            PADILLA_SEED.E0,
        ]
    else:
        bounds = [
            (0.05, 0.55),
            (0.01, 0.25),
            (-3.5, -0.5),
            (0.3, 1.5),
        ]
        x0_seed = [RYDEN_SEED.mu_g, RYDEN_SEED.sig_g, RYDEN_SEED.mu, RYDEN_SEED.sig]

    def objective(x: np.ndarray, n_draw: int = N_MODEL_FIT) -> float:
        if with_E0:
            p = ShapeParams(mu_g=x[0], sig_g=x[1], mu=x[2], sig=x[3], E0=x[4])
        else:
            p = ShapeParams(mu_g=x[0], sig_g=x[1], mu=x[2], sig=x[3], E0=0.0)
        rng = np.random.default_rng(seed + int(1e6 * abs(x[0] + x[2])) % 100000)
        model = model_histogram(p, n_draw, data_counts, edges, rng)
        return chi2_hist(data_counts, model)

    print(
        f"[*] Fitting {'Padilla (5-param)' if with_E0 else 'Ryden (4-param)'} "
        f"with differential_evolution ...",
        flush=True,
    )
    seed_chi2 = objective(np.asarray(x0_seed, dtype=float))
    print(f"    literature-seed chi2={seed_chi2:.2f}", flush=True)

    result = differential_evolution(
        objective,
        bounds=bounds,
        seed=seed,
        maxiter=18,
        popsize=10,
        mutation=(0.5, 1.0),
        recombination=0.7,
        polish=True,
        workers=1,
        updating="immediate",
        atol=1.0,
        tol=0.02,
    )
    x = result.x
    if with_E0:
        best = ShapeParams(mu_g=x[0], sig_g=x[1], mu=x[2], sig=x[3], E0=x[4])
    else:
        best = ShapeParams(mu_g=x[0], sig_g=x[1], mu=x[2], sig=x[3], E0=0.0)
    # Final chi2 with full n_model draws
    chi2 = objective(x, n_draw=n_model)
    print(f"    best chi2={chi2:.2f}  params={best}", flush=True)
    return best, chi2


def plot_hist_compare(
    data_counts: np.ndarray,
    model_counts: np.ndarray,
    centers: np.ndarray,
    out: Path,
    *,
    title: str,
) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    width = centers[1] - centers[0] if len(centers) > 1 else 0.05
    ax.bar(
        centers,
        data_counts,
        width=width * 0.9,
        color="#9ecae1",
        edgecolor="white",
        label="LS EXP data",
    )
    ax.plot(centers, model_counts, "o-", color="#e41a1c", ms=4, lw=1.5, label="Model")
    ax.set_xlabel(r"projected $b/a$")
    ax.set_ylabel("Counts (model scaled to data N)")
    ax.set_title(title)
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_param_slices(
    data_counts: np.ndarray,
    edges: np.ndarray,
    best: ShapeParams,
    out: Path,
    *,
    with_E0: bool,
    n_model: int,
    seed: int,
) -> None:
    """Local chi2 slices around the best fit (mu-sig and mu_g-sig_g)."""
    rng = np.random.default_rng(seed)
    mu_grid = np.linspace(best.mu - 0.8, best.mu + 0.8, 9)
    sig_grid = np.linspace(max(0.15, best.sig - 0.4), best.sig + 0.4, 9)
    zg = np.full((len(sig_grid), len(mu_grid)), np.nan)
    for i, s in enumerate(sig_grid):
        for j, m in enumerate(mu_grid):
            p = ShapeParams(best.mu_g, best.sig_g, float(m), float(s), best.E0 if with_E0 else 0.0)
            mod = model_histogram(p, n_model // 2, data_counts, edges, rng)
            zg[i, j] = chi2_hist(data_counts, mod)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    ax = axes[0]
    im = ax.pcolormesh(mu_grid, sig_grid, zg, shading="auto", cmap="viridis_r")
    ax.plot(best.mu, best.sig, "r*", ms=12, label="best")
    ax.set_xlabel(r"$\mu=\langle\ln\varepsilon\rangle$")
    ax.set_ylabel(r"$\sigma(\ln\varepsilon)$")
    ax.set_title(r"$\chi^2$ slice (eps)")
    ax.legend(loc="best")
    fig.colorbar(im, ax=ax, fraction=0.046)

    mug = np.linspace(max(0.05, best.mu_g - 0.15), min(0.55, best.mu_g + 0.15), 9)
    sigg = np.linspace(max(0.01, best.sig_g - 0.08), best.sig_g + 0.08, 9)
    zg2 = np.full((len(sigg), len(mug)), np.nan)
    for i, s in enumerate(sigg):
        for j, m in enumerate(mug):
            p = ShapeParams(float(m), float(s), best.mu, best.sig, best.E0 if with_E0 else 0.0)
            mod = model_histogram(p, n_model // 2, data_counts, edges, rng)
            zg2[i, j] = chi2_hist(data_counts, mod)
    ax = axes[1]
    im2 = ax.pcolormesh(mug, sigg, zg2, shading="auto", cmap="viridis_r")
    ax.plot(best.mu_g, best.sig_g, "r*", ms=12, label="best")
    ax.set_xlabel(r"$\mu_\gamma$ (thickness)")
    ax.set_ylabel(r"$\sigma_\gamma$")
    ax.set_title(r"$\chi^2$ slice (thickness)")
    ax.legend(loc="best")
    fig.colorbar(im2, ax=ax, fraction=0.046)
    fig.suptitle("Local parameter slices around best fit")
    fig.tight_layout()
    fig.savefig(out, dpi=250, bbox_inches="tight")
    plt.close(fig)


def plot_inclination(params: ShapeParams, out: Path, *, title: str) -> None:
    tab = inclination_vs_q(params, n=400_000, bin_width=BIN_WIDTH, rng=np.random.default_rng(0))
    q = tab["q"]
    # Hubble reference
    hub = np.array([hubble_cosi_from_ba(float(qq), q0=Q0) if qq > Q0 else np.nan for qq in q])

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.fill_between(q, tab["cosi_lo"], tab["cosi_hi"], color="#377eb8", alpha=0.25, label="16-84%")
    ax.plot(q, tab["cosi_med"], "o-", color="#377eb8", ms=4, label="elliptical-disk median")
    ax.plot(q, hub, "s--", color="#e41a1c", ms=3, label=rf"Hubble $q_0={Q0:g}$")
    ax.axvline(BA_FACE_CAP, color="#4daf4a", ls=":", lw=1.2, label=r"ad-hoc $b/a=0.8$")
    ax.set_xlabel(r"projected $b/a$")
    ax.set_ylabel(r"$\cos i$")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_title(title)
    ax.legend(loc="best", fontsize=8)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(tab).to_csv(out.with_suffix(".csv"), index=False)


def plot_cdf(cosi: np.ndarray, mag_limit: float, out: Path, *, label: str, color: str) -> dict:
    x = np.linspace(0, 1, 401)
    med = float(np.median(cosi)) if len(cosi) else float("nan")
    fig, ax = plt.subplots(figsize=(6.5, 6.5))
    ax.plot((0, 1), (0, 1), "k--", lw=1.2, label="Uniform")
    if len(cosi):
        ax.plot(x, empirical_cdf(cosi, x), color=color, lw=2.0, label=label)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"$\cos(i)$ (elliptical-disk model; strict $b/a>q_0$)")
    ax.set_ylabel("Cumulative distribution")
    ax.set_title(f"{label}\nmag limit={mag_limit:g}  N={len(cosi):,}  med={med:.3f}")
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return {"mag_limit": mag_limit, "n": len(cosi), "median_cosi": round(med, 4) if len(cosi) else np.nan}


def write_summary(
    out_dir: Path,
    *,
    track: str,
    params: ShapeParams,
    chi2: float,
    n_data: int,
    cdf_rows: list[dict],
    notes: str,
) -> None:
    med_eps = float(np.exp(params.mu))
    text = f"""# {track}

## Fit

| Parameter | Value |
|-----------|------:|
| N (fit pool) | {n_data:,} |
| chi2 | {chi2:.2f} |
| mu_g (mean thickness C/A) | {params.mu_g:.4f} |
| sig_g | {params.sig_g:.4f} |
| mu (mean ln eps) | {params.mu:.4f} |
| sig | {params.sig:.4f} |
| median eps = exp(mu) | {med_eps:.4f} |
| E0 (Padilla dust) | {params.E0:.4f} |

Literature seeds: Ryden (mu_g~0.22, ln eps~-1.85); Padilla (ln e~-2.33, E0~0.45).

## Mag-cut model cos(i) CDFs (strict b/a > q0={Q0:g})

| mag_limit | N | median cos(i) model | median cos(i) Hubble |
|----------:|--:|--------------------:|---------------------:|
"""
    for r in cdf_rows:
        mh = r.get("median_cosi_hubble", "")
        text += f"| {r['mag_limit']:g} | {r['n']:,} | {r['median_cosi']} | {mh} |\n"
    text += f"\n{notes}\n"
    (out_dir / "summary.md").write_text(text, encoding="utf-8")


def run_ryden(
    mag: np.ndarray, ba: np.ndarray, out_dir: Path, *, n_model: int = N_MODEL
) -> ShapeParams:
    out_dir.mkdir(parents=True, exist_ok=True)
    data_counts, edges, centers = histogram_q(ba, bin_width=BIN_WIDTH)
    best, chi2 = fit_params(
        data_counts, edges, with_E0=False, n_model=n_model, seed=RNG_SEED
    )
    rng = np.random.default_rng(RNG_SEED)
    model_c = model_histogram(best, n_model, data_counts, edges, rng)

    plot_hist_compare(
        data_counts,
        model_c,
        centers,
        out_dir / "ba_hist_data_vs_model.png",
        title="LS EXP vs Ryden elliptical-disk model",
    )
    plot_param_slices(
        data_counts,
        edges,
        best,
        out_dir / "param_contours.png",
        with_E0=False,
        n_model=n_model,
        seed=RNG_SEED,
    )
    plot_inclination(
        best,
        out_dir / "inclination_vs_ba.png",
        title="Ryden elliptical disk: cos i vs b/a",
    )
    sampler = ConditionalCosiSampler(
        best, np.random.default_rng(RNG_SEED), n_model=N_MODEL_COSI
    )

    params_out = {
        "track": "scaled_ryden",
        "params": params_jsonable(best),
        "chi2": float(chi2),
        "n_data": int(len(ba)),
        "bin_width": BIN_WIDTH,
        "n_model": n_model,
        "dust": "Unterborn A1 after shape fit",
        "literature": {"ryden_2004": params_jsonable(RYDEN_SEED)},
    }
    (out_dir / "fit_params.json").write_text(json.dumps(params_out, indent=2), encoding="utf-8")

    # Unterborn A1 + model cos(i) CDFs
    m_face = face_on_mag(mag, ba)
    cdf_rng = np.random.default_rng(RNG_SEED + 100)
    cdf_rows: list[dict] = []
    cdfs = out_dir / "cdfs"
    for lim in MAG_CUTS:
        mask = (m_face <= lim) & (ba > Q0)
        cosi = sampler.sample(ba[mask], cdf_rng)
        row = plot_cdf(
            cosi,
            lim,
            cdfs / f"mag{int(lim)}.png",
            label="LS Ryden+A1 (model cos i)",
            color="#377eb8",
        )
        row["median_cosi_hubble"] = (
            round(float(np.median(cosi_hubble(ba[mask]))), 4) if mask.sum() else np.nan
        )
        cdf_rows.append(row)
        print(
            f"    A1 mag<={lim:g}: N={row['n']:,} med_model={row['median_cosi']} "
            f"med_hubble={row['median_cosi_hubble']}",
            flush=True,
        )
    pd.DataFrame(cdf_rows).to_csv(cdfs / "summary.csv", index=False)

    write_summary(
        out_dir,
        track="scaled_ryden",
        params=best,
        chi2=chi2,
        n_data=len(ba),
        cdf_rows=cdf_rows,
        notes=(
            "Shape fit is Ryden photometry-only. Dust handled separately via Unterborn "
            "face-on mag re-cut for the CDFs in `cdfs/`. Per-galaxy cos(i) is DRAWN from "
            "the model posterior P(cos i | b/a) (one sample each), not the per-bin median: "
            "a given b/a maps to a distribution of inclinations, so using the median "
            "collapses every galaxy onto a narrow band and produces a degenerate "
            "near-vertical CDF.\n\n"
            "CAVEAT (see CIRCULARITY_CHECK.md): the model ASSUMES isotropic orientation, "
            "so the recovered cos(i) CDF is ~uniform by construction whenever the b/a fit "
            "is decent. It is NOT independent evidence that LS is isotropic. The useful "
            "products are the fitted shape params and the per-galaxy P(cos i | b/a) "
            "posterior (for individual FRB hosts).\n\n"
            "Note: best-fit face-on ellipticity is much larger than Ryden's SDSS value "
            "(median eps~0.58 vs ~0.16). That is expected for Tractor EXP-only: REX removes "
            "near-round disks, so the apparent q distribution requires strong intrinsic "
            "ellipticity (and/or selection) to suppress q~1.\n"
        ),
    )
    return best


def run_padilla(
    mag: np.ndarray, ba: np.ndarray, out_dir: Path, *, n_model: int = N_MODEL
) -> ShapeParams:
    out_dir.mkdir(parents=True, exist_ok=True)
    data_counts, edges, centers = histogram_q(ba, bin_width=BIN_WIDTH)
    best, chi2 = fit_params(
        data_counts, edges, with_E0=True, n_model=n_model, seed=RNG_SEED + 1
    )
    rng = np.random.default_rng(RNG_SEED + 1)
    model_c = model_histogram(best, n_model, data_counts, edges, rng)

    plot_hist_compare(
        data_counts,
        model_c,
        centers,
        out_dir / "ba_hist_data_vs_model.png",
        title="LS EXP vs Padilla shape+dust model",
    )
    plot_param_slices(
        data_counts,
        edges,
        best,
        out_dir / "param_contours.png",
        with_E0=True,
        n_model=n_model,
        seed=RNG_SEED + 1,
    )
    plot_inclination(
        best,
        out_dir / "inclination_vs_ba.png",
        title=rf"Padilla elliptical disk ($E_0={best.E0:.2f}$): cos i vs b/a",
    )
    sampler = ConditionalCosiSampler(
        best, np.random.default_rng(RNG_SEED + 1), n_model=N_MODEL_COSI
    )

    params_out = {
        "track": "scaled_padilla",
        "params": params_jsonable(best),
        "chi2": float(chi2),
        "n_data": int(len(ba)),
        "bin_width": BIN_WIDTH,
        "n_model": n_model,
        "dust": "Padilla E0 in psi(theta); photometric adaptation (no 1/Vmax)",
        "literature": {"padilla_strauss_2008": params_jsonable(PADILLA_SEED)},
    }
    (out_dir / "fit_params.json").write_text(json.dumps(params_out, indent=2), encoding="utf-8")

    # CDFs: observed mag cut; cos(i) from fitted elliptical+dust model
    cdf_rng = np.random.default_rng(RNG_SEED + 101)
    cdf_rows: list[dict] = []
    cdfs = out_dir / "cdfs"
    for lim in MAG_CUTS:
        mask = (mag <= lim) & (ba > Q0)
        cosi = sampler.sample(ba[mask], cdf_rng)
        row = plot_cdf(
            cosi,
            lim,
            cdfs / f"mag{int(lim)}.png",
            label="LS Padilla (model cos i)",
            color="#984ea3",
        )
        row["median_cosi_hubble"] = (
            round(float(np.median(cosi_hubble(ba[mask]))), 4) if mask.sum() else np.nan
        )
        cdf_rows.append(row)
        print(
            f"    padilla mag<={lim:g}: N={row['n']:,} med_model={row['median_cosi']} "
            f"med_hubble={row['median_cosi_hubble']}",
            flush=True,
        )
    pd.DataFrame(cdf_rows).to_csv(cdfs / "summary.csv", index=False)

    write_summary(
        out_dir,
        track="scaled_padilla",
        params=best,
        chi2=chi2,
        n_data=len(ba),
        cdf_rows=cdf_rows,
        notes=(
            "Joint Padilla shape+E0 fit. LS has no redshifts; psi(theta) uses "
            "10^{-0.4 E(theta)} only (not full 1/Vmax+LF). CDFs use observed mag cuts, "
            "strict b/a>q0. Per-galaxy cos(i) is DRAWN from the model posterior "
            "P(cos i | b/a) (one sample each), not the per-bin median.\n\n"
            "Large E0 and large face-on ellipticity partly compensate for the EXP/REX "
            "selection (deficit of round systems) in this photometric adaptation.\n"
        ),
    )
    return best


def plot_compare_all(
    ryden_params: ShapeParams,
    padilla_params: ShapeParams,
    mag: np.ndarray,
    ba: np.ndarray,
    out: Path,
) -> None:
    """Overlay CDFs: ad-hoc scaled (Hubble/0.8), Ryden+A1 model, Padilla model."""
    m_face = face_on_mag(mag, ba)
    cosi_scale = float(hubble_cosi_from_ba(BA_FACE_CAP, q0=Q0))
    samp_r = ConditionalCosiSampler(
        ryden_params, np.random.default_rng(RNG_SEED), n_model=N_MODEL_COSI
    )
    samp_p = ConditionalCosiSampler(
        padilla_params, np.random.default_rng(RNG_SEED + 1), n_model=N_MODEL_COSI
    )
    cmp_rng = np.random.default_rng(RNG_SEED + 200)

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    x = np.linspace(0, 1, 401)
    for ax, lim in zip(axes, MAG_CUTS):
        ax.plot((0, 1), (0, 1), "k--", lw=1.0, label="Uniform")

        m_ad = (mag <= lim) & (ba > Q0) & (ba <= BA_FACE_CAP)
        c_ad = np.clip(cosi_hubble(ba[m_ad]) / cosi_scale, 0, 1)
        m_ry = (m_face <= lim) & (ba > Q0)
        c_ry = samp_r.sample(ba[m_ry], cmp_rng)
        m_pa = (mag <= lim) & (ba > Q0)
        c_pa = samp_p.sample(ba[m_pa], cmp_rng)

        for cosi, lab, col in (
            (c_ad, f"ad-hoc scaled (N={len(c_ad):,})", "#4daf4a"),
            (c_ry, f"Ryden+A1 (N={len(c_ry):,})", "#377eb8"),
            (c_pa, f"Padilla (N={len(c_pa):,})", "#984ea3"),
        ):
            if len(cosi) == 0:
                continue
            med = float(np.median(cosi))
            ax.plot(
                x,
                empirical_cdf(cosi, x),
                color=col,
                lw=2.0,
                label=f"{lab}, med={med:.3f}",
            )
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel(r"$\cos(i)$")
        ax.set_title(f"mag limit = {lim:g}")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper left", fontsize=7)
    axes[0].set_ylabel("Cumulative distribution")
    fig.suptitle("LS EXP: ad-hoc scaled vs Ryden+Unterborn vs Padilla", fontsize=12)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode", choices=("ryden", "padilla", "both"), default="both")
    p.add_argument("--n-model", type=int, default=N_MODEL)
    args = p.parse_args()

    n_model = args.n_model
    mag, ba = load_ls(REPO_ROOT / LS_CATALOG_V2_EXP_DEFAULT)
    ryden_dir = OUT_ROOT / "scaled_ryden"
    padilla_dir = OUT_ROOT / "scaled_padilla"

    ryden_p = padilla_p = None
    if args.mode in ("ryden", "both"):
        print("[*] === scaled_ryden ===", flush=True)
        ryden_p = run_ryden(mag, ba, ryden_dir, n_model=n_model)
    if args.mode in ("padilla", "both"):
        print("[*] === scaled_padilla ===", flush=True)
        padilla_p = run_padilla(mag, ba, padilla_dir, n_model=n_model)
    if args.mode == "both" and ryden_p is not None and padilla_p is not None:
        plot_compare_all(
            ryden_p,
            padilla_p,
            mag,
            ba,
            OUT_ROOT / "scaled_ryden_padilla_cdf_compare.png",
        )
        print(
            f"[*] Wrote compare -> {OUT_ROOT / 'scaled_ryden_padilla_cdf_compare.png'}",
            flush=True,
        )

    print("[*] Done.", flush=True)


if __name__ == "__main__":
    main()
