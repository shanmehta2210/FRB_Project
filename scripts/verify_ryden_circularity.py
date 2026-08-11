"""
Is the ~uniform cos(i) CDF from scaled_ryden a real result or a tautology?

Ryden (2004) / Padilla & Strauss (2008) *assume* random orientation (uniform in
cos theta) as an INPUT, then fit the intrinsic shape distribution to the observed
b/a histogram. So the model marginal P(cos i) is uniform BY CONSTRUCTION.

When we sample cos(i) ~ P_model(cos i | b/a) for the observed galaxies we recover

    P_rec(cos i) = INT P_model(cos i | q) P_obs(q) dq.

If the shape fit is good, P_obs(q) ~= P_model(q), and then P_rec -> P_model(cos i)
= uniform. So the diagonal CDF is guaranteed whenever the b/a fit is decent; it is
NOT an independent test that LS orientations are isotropic (shape/orientation are
degenerate given b/a alone).

This script demonstrates that empirically with three controls. Run from repo root::

    python scripts/verify_ryden_circularity.py
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

from elliptical_disk_model import ConditionalCosiSampler, ShapeParams, generate_q  # noqa: E402
from fit_ls_scaled_elliptical import load_ls  # noqa: E402
from null_catalog_utils import LS_CATALOG_V2_EXP_DEFAULT, Q0  # noqa: E402
from pipeline_null_plot_utils import REPO_ROOT  # noqa: E402

OUT_DIR = REPO_ROOT / "plots" / "plots_null" / "v2" / "ls_audit" / "scaled_ryden"
N_POOL = 2_000_000


def ks_uniform(cosi: np.ndarray) -> float:
    cosi = cosi[np.isfinite(cosi)]
    if cosi.size == 0:
        return float("nan")
    return float(kstest(cosi, "uniform").statistic)


def ecdf(vals: np.ndarray, x: np.ndarray) -> np.ndarray:
    s = np.sort(vals[np.isfinite(vals)])
    return np.searchsorted(s, x, side="right") / max(1, len(s))


def main() -> None:
    params = ShapeParams(**json.loads((OUT_DIR / "fit_params.json").read_text())["params"])
    print(f"[*] Fitted Ryden params: {params}")

    rng = np.random.default_rng(0)
    sampler = ConditionalCosiSampler(params, np.random.default_rng(1), n_model=N_POOL)

    # (0) model marginal cos(theta) is uniform BY CONSTRUCTION
    _, cos_model, _ = generate_q(params, 1_000_000, rng, return_angles=True)
    print(f"[0] model marginal cos(theta): KS vs uniform = {ks_uniform(cos_model):.4f}")

    # (1) real observed b/a -> recovered cos(i)
    mag, ba = load_ls(REPO_ROOT / LS_CATALOG_V2_EXP_DEFAULT)
    obs = ba[ba > Q0]
    cosi_obs = sampler.sample(obs, rng)
    ks_obs = ks_uniform(cosi_obs)
    print(f"[1] observed b/a -> recovered cos(i): KS vs uniform = {ks_obs:.4f}")

    # (2) CONTROL: synthetic b/a drawn FROM the model itself (perfect fit case)
    q_syn = generate_q(params, obs.size, rng)
    q_syn = q_syn[q_syn > Q0]
    cosi_syn = sampler.sample(q_syn, rng)
    ks_syn = ks_uniform(cosi_syn)
    print(f"[2] model-drawn b/a -> recovered cos(i): KS vs uniform = {ks_syn:.4f}  (tautology floor)")

    # (3) COUNTER-CONTROL: deliberately distort the observed b/a distribution
    #     (upweight round galaxies) WITHOUT refitting. If the CDF still went uniform
    #     it would prove no data dependence; instead it departs -> uniformity is
    #     inherited from P_obs ~= P_model, i.e. from the fit, not from isotropy.
    w = obs**4  # heavily favour round (face-on-looking) systems
    w = w / w.sum()
    idx = rng.choice(obs.size, size=obs.size, p=w)
    cosi_dist = sampler.sample(obs[idx], rng)
    ks_dist = ks_uniform(cosi_dist)
    print(f"[3] round-biased b/a -> recovered cos(i): KS vs uniform = {ks_dist:.4f}  (data DOES matter)")

    # KS between observed and model b/a histograms (the actual fit residual)
    ks_ba = float(kstest(obs, lambda x: ecdf(q_syn, x)).statistic)
    print(f"[*] KS(observed b/a, model b/a) = {ks_ba:.4f}  (drives the [1] residual)")

    # figure
    x = np.linspace(0, 1, 401)
    fig, ax = plt.subplots(figsize=(7.5, 7))
    ax.plot((0, 1), (0, 1), "k--", lw=1.3, label="Uniform (assumed input)")
    ax.plot(x, ecdf(cosi_syn, x), color="#999999", lw=2.0,
            label=f"[2] model-drawn b/a (KS={ks_syn:.3f})")
    ax.plot(x, ecdf(cosi_obs, x), color="#377eb8", lw=2.2,
            label=f"[1] observed b/a (KS={ks_obs:.3f})")
    ax.plot(x, ecdf(cosi_dist, x), color="#e41a1c", lw=2.2,
            label=f"[3] round-biased b/a (KS={ks_dist:.3f})")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"$\cos(i)$ sampled from $P_{\rm model}(\cos i\,|\,b/a)$")
    ax.set_ylabel("Cumulative distribution")
    ax.set_title(
        "Ryden cos(i) CDF is inherited from the b/a fit, not a test of isotropy\n"
        "[1] and [2] hug the diagonal by construction; [3] shows real data-dependence"
    )
    ax.legend(loc="upper left", fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    out = OUT_DIR / "circularity_check.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[*] Wrote {out}")


if __name__ == "__main__":
    main()
