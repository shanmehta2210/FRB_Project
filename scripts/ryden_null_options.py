"""
Which b/a -> cos(i) conversion should be the null for FRB comparison?

Compares three ways to turn LS EXP b/a into a cos(i) CDF, at each mag cut:

  (A) ad-hoc "scaled"      : deterministic Hubble thin-disk, cos i = f(b/a), q0.
                             = elliptical-disk model with a DELTA-function shape.
  (B) Ryden REFIT-to-LS    : shape params fit to LS b/a -> ~uniform BY CONSTRUCTION
                             (isotropy assumed + good b/a fit). Circular; no power.
  (C) Ryden FIXED-lit      : shape params FROZEN at Ryden (2004) SDSS values, NOT
                             refit. LS b/a is then NOT fully explained by those
                             shapes, so cos(i) departs from uniform -> curved,
                             like (A) but with realistic shape scatter.

(C) is the physically-motivated generalisation of (A) and is the recommended null:
freeze shapes from the literature, apply the SAME frozen model to FRB hosts.

Run from repo root::

    python scripts/ryden_null_options.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kstest

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from elliptical_disk_model import (  # noqa: E402
    PADILLA_SEED,
    RYDEN_SEED,
    ConditionalCosiSampler,
    ShapeParams,
)
from fit_ls_scaled_elliptical import cosi_hubble, load_ls  # noqa: E402
from null_catalog_utils import LS_CATALOG_V2_EXP_DEFAULT, Q0, face_on_mag  # noqa: E402
from pipeline_null_plot_utils import REPO_ROOT  # noqa: E402

OUT_DIR = REPO_ROOT / "plots" / "plots_null" / "v2" / "ls_audit" / "scaled_ryden"
MAG_CUTS = (20.0, 21.0, 22.0)
N_POOL = 2_000_000


def ecdf(vals: np.ndarray, x: np.ndarray) -> np.ndarray:
    s = np.sort(vals[np.isfinite(vals)])
    return np.searchsorted(s, x, side="right") / max(1, len(s))


def ks_uniform(v: np.ndarray) -> float:
    v = v[np.isfinite(v)]
    return float(kstest(v, "uniform").statistic) if v.size else float("nan")


def main() -> None:
    refit = ShapeParams(**json.loads((OUT_DIR / "fit_params.json").read_text())["params"])
    print(f"[*] refit-to-LS params : {refit}")
    print(f"[*] Ryden-2004 lit     : {RYDEN_SEED}")
    print(f"[*] Padilla-2008 lit   : {PADILLA_SEED}")

    rng = np.random.default_rng(7)
    samp_refit = ConditionalCosiSampler(refit, np.random.default_rng(1), n_model=N_POOL)
    samp_lit = ConditionalCosiSampler(RYDEN_SEED, np.random.default_rng(2), n_model=N_POOL)
    samp_pad = ConditionalCosiSampler(PADILLA_SEED, np.random.default_rng(3), n_model=N_POOL)

    mag, ba = load_ls(REPO_ROOT / LS_CATALOG_V2_EXP_DEFAULT)
    m_face = face_on_mag(mag, ba)

    x = np.linspace(0, 1, 401)
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2), sharey=True)
    for ax, lim in zip(axes, MAG_CUTS):
        m_obs = (mag <= lim) & (ba > Q0)
        m_a1 = (m_face <= lim) & (ba > Q0)
        q_obs = ba[m_obs]

        series = [
            ("A: scaled (Hubble thin-disk)", cosi_hubble(q_obs), "#4daf4a"),
            ("B: Ryden refit-to-LS", samp_refit.sample(ba[m_a1], rng), "#377eb8"),
            ("C: Ryden fixed-lit (2004)", samp_lit.sample(q_obs, rng), "#e41a1c"),
            ("C': Padilla fixed-lit (2008)", samp_pad.sample(q_obs, rng), "#984ea3"),
        ]
        ax.plot((0, 1), (0, 1), "k--", lw=1.1, label="Uniform")
        for lab, cosi, col in series:
            med = float(np.median(cosi[np.isfinite(cosi)]))
            ax.plot(x, ecdf(cosi, x), color=col, lw=2.0,
                    label=f"{lab}\n  med={med:.3f}, KS={ks_uniform(cosi):.3f}")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel(r"$\cos(i)$")
        ax.set_title(f"mag limit = {lim:g}  (N={int(m_obs.sum()):,})")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper left", fontsize=7)
    axes[0].set_ylabel("Cumulative distribution")
    fig.suptitle(
        "LS EXP cos(i) null: refit-to-sample -> uniform (circular); "
        "fixed-literature shapes -> curved, like scaled",
        fontsize=12,
    )
    fig.tight_layout()
    out = OUT_DIR / "null_options_compare.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[*] Wrote {out}")


if __name__ == "__main__":
    main()
