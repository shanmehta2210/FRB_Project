"""
Prove that our ad-hoc `scaled` is the Ryden/elliptical-disk model in its fully
degenerate limit, and quantify the physical approximations required.

Analytic reduction
------------------
Ryden projected axis ratio (Binney 1985 / Ryden 2004 eqs 12-15). Put the intrinsic
face-on ellipticity to zero, eps -> 0 (a PERFECTLY CIRCULAR disk). Then e2=eps(2-eps)=0
so B=0 and the azimuthal (phi) dependence vanishes, leaving

    A = cos^2(theta) + gamma^2 sin^2(theta),  C = 1,
    q^2 = A/C = cos^2(theta) (1 - gamma^2) + gamma^2.

Invert:  cos(theta) = sqrt( (q^2 - gamma^2) / (1 - gamma^2) ).

That is EXACTLY the Hubble formula with q0 = gamma, and it is deterministic (no phi
scatter). So the elliptical-disk model collapses to the Hubble deprojection iff:

  (1) eps = 0            intrinsic disks are perfectly circular  (mu -> -inf)
  (2) sigma_gamma = 0    a single universal thickness            (delta function)
  (3) gamma = q0 = 0.2   thickness equals the edge-on axis ratio

`scaled` adds one non-shape step: the REX renormalisation cos i -> cos i / H(0.8)
(equivalently relocating the face-on edge to the empirically-observed b/a=0.8 cap).

This script (a) verifies the sampled degenerate-Ryden cos(i) equals the Hubble value,
and (b) overlays the degenerate-Ryden+renorm CDF on `scaled` for LS at each mag cut.

Run from repo root::

    python scripts/ryden_reduces_to_scaled.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from elliptical_disk_model import ConditionalCosiSampler, ShapeParams, generate_q  # noqa: E402
from fit_ls_scaled_elliptical import load_ls  # noqa: E402
from null_catalog_utils import LS_CATALOG_V2_EXP_DEFAULT, Q0, face_on_mag, hubble_cosi_from_ba  # noqa: E402
from pipeline_null_plot_utils import REPO_ROOT  # noqa: E402

OUT_DIR = REPO_ROOT / "plots" / "plots_null" / "v2" / "ls_audit" / "scaled_ryden_fixed"
MAG_CUTS = (20.0, 21.0, 22.0)
BA_FACE_CAP = 0.8
N_POOL = 3_000_000
H08 = float(hubble_cosi_from_ba(BA_FACE_CAP, q0=Q0))  # H(0.8) ~ 0.7906

# Degenerate elliptical-disk params == Hubble deprojection with q0=0.2
DEGENERATE = ShapeParams(mu_g=Q0, sig_g=1e-4, mu=-12.0, sig=1e-3, E0=0.0)


def hubble(q: np.ndarray) -> np.ndarray:
    val = (q * q - Q0 * Q0) / (1.0 - Q0 * Q0)
    return np.sqrt(np.clip(val, 0.0, 1.0))


def scaled_cosi(q: np.ndarray) -> np.ndarray:
    return np.clip(hubble(q) / H08, 0.0, 1.0)


def ecdf(v: np.ndarray, x: np.ndarray) -> np.ndarray:
    s = np.sort(v[np.isfinite(v)])
    return np.searchsorted(s, x, side="right") / max(1, len(s))


def main() -> None:
    # (a) verify degenerate model reproduces Hubble deterministically
    rng = np.random.default_rng(0)
    q, cost, _ = generate_q(DEGENERATE, 500_000, rng, return_angles=True)
    ok = q > Q0
    resid = np.abs(cost[ok] - hubble(q[ok]))
    print(f"[a] degenerate-Ryden cos(theta) vs Hubble(q): "
          f"max|diff|={resid.max():.2e}, median={np.median(resid):.2e}")

    # (b) LS CDFs: scaled vs degenerate-Ryden(sampled)+renorm
    sampler = ConditionalCosiSampler(DEGENERATE, np.random.default_rng(1), n_model=N_POOL)
    mag, ba = load_ls(REPO_ROOT / LS_CATALOG_V2_EXP_DEFAULT)
    m_face = face_on_mag(mag, ba)
    srng = np.random.default_rng(5)

    x = np.linspace(0, 1, 401)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2), sharey=True)
    for ax, lim in zip(axes, MAG_CUTS):
        mask = (m_face <= lim) & (ba > Q0) & (ba <= BA_FACE_CAP)
        c_scaled = scaled_cosi(ba[mask])
        c_ryden = np.clip(sampler.sample(ba[mask], srng) / H08, 0.0, 1.0)
        ks = float(np.max(np.abs(ecdf(c_scaled, x) - ecdf(c_ryden, x))))
        ax.plot((0, 1), (0, 1), "k--", lw=1.0, label="Uniform")
        ax.plot(x, ecdf(c_scaled, x), color="#4daf4a", lw=3.2,
                label=f"scaled (Hubble/0.8)  med={np.median(c_scaled):.3f}")
        ax.plot(x, ecdf(c_ryden, x), color="#e41a1c", lw=1.5, ls="--",
                label=f"degenerate-Ryden + renorm  med={np.median(c_ryden):.3f}")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel(r"$\cos(i)$")
        ax.set_title(f"mag limit={lim:g}   max CDF diff = {ks:.4f}")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper left", fontsize=8)
        print(f"[b] mag<={lim:g}: max CDF diff (scaled vs degenerate-Ryden) = {ks:.4f}")
    axes[0].set_ylabel("Cumulative distribution")
    fig.suptitle(
        "scaled == elliptical-disk model in the degenerate limit "
        r"($\varepsilon=0,\ \sigma_\gamma=0,\ \gamma=q_0=0.2$) + REX renorm /0.8",
        fontsize=12)
    fig.tight_layout()
    out = OUT_DIR / "ryden_equals_scaled.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[*] Wrote {out}")


if __name__ == "__main__":
    main()
