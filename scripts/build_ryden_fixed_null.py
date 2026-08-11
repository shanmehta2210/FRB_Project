"""
Ryden FIXED-literature shape model (+ Unterborn A1 dust) as the LS cos(i) null.

The REX trilemma (see summary.md):
  1. Fit our own shapes to LS -> the shape distribution absorbs ALL the b/a
     structure and cos(i) collapses to ~uniform (circular; isotropy assumed in,
     isotropy out). See scaled_ryden/CIRCULARITY_CHECK.md.
  2. Use SDSS Ryden (2004) shapes DIRECTLY -> LS is REX-truncated (round/face-on
     disks removed), which the SDSS-calibrated model reads as real inclination,
     so the median cos(i) is pushed too edge-on (~0.41 after dust).
  3. Use SDSS Ryden shapes but ASSUME b/a=0.8 is the face-on ceiling (de-REX the
     axis ratio, b/a -> b/a/0.8) -> curved, median ~0.5, matching ad-hoc scaled.
     Imposing the 0.8 ceiling is the price of directly using Ryden on REX data.

Two tracks emitted, each raw + Unterborn A1 dust, strict b/a > q0:
  pure_lit    : cos(i) ~ P_Ryden2004(cos i | b/a)                 (option 2)
  ceiling08   : cos(i) ~ P_Ryden2004(cos i | b/a/0.8), b/a<=0.8   (option 3)
Reference overlay: ad-hoc "scaled" (Hubble cos i / cos i|_0.8) + A1.

Outputs under plots/plots_null/v2/ls_audit/scaled_ryden_fixed/. FRB hosts get the
IDENTICAL frozen sampler + the same b/a treatment for a like-for-like comparison.

Run from repo root::

    python scripts/build_ryden_fixed_null.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from elliptical_disk_model import RYDEN_SEED, ConditionalCosiSampler  # noqa: E402
from fit_ls_scaled_elliptical import load_ls  # noqa: E402
from null_catalog_utils import (  # noqa: E402
    LS_CATALOG_V2_EXP_DEFAULT,
    Q0,
    face_on_mag,
    hubble_cosi_from_ba,
)
from pipeline_null_plot_utils import REPO_ROOT  # noqa: E402

OUT_DIR = REPO_ROOT / "plots" / "plots_null" / "v2" / "ls_audit" / "scaled_ryden_fixed"
MAG_CUTS = (20.0, 21.0, 22.0)
BA_FACE_CAP = 0.8
N_POOL = 2_000_000
COSI_SCALE = float(hubble_cosi_from_ba(BA_FACE_CAP, q0=Q0))  # ~0.791


def scaled_cosi(ba: np.ndarray) -> np.ndarray:
    """Ad-hoc scaled reference: Hubble cos(i) normalised so b/a=0.8 -> 1."""
    val = (ba * ba - Q0 * Q0) / (1.0 - Q0 * Q0)
    cosi = np.sqrt(np.clip(val, 0.0, 1.0))
    return np.clip(cosi / COSI_SCALE, 0.0, 1.0)


def empirical_cdf(vals: np.ndarray, x: np.ndarray) -> np.ndarray:
    s = np.sort(vals[np.isfinite(vals)])
    return np.searchsorted(s, x, side="right") / max(1, len(s))


def med(v: np.ndarray) -> float:
    v = v[np.isfinite(v)]
    return float(np.median(v)) if v.size else float("nan")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / "cdfs").mkdir(exist_ok=True)
    print(f"[*] Frozen Ryden-2004 shape params: {RYDEN_SEED}")
    sampler = ConditionalCosiSampler(RYDEN_SEED, np.random.default_rng(2), n_model=N_POOL)
    rng = np.random.default_rng(11)

    mag, ba = load_ls(REPO_ROOT / LS_CATALOG_V2_EXP_DEFAULT)
    m_face = face_on_mag(mag, ba)

    x = np.linspace(0, 1, 401)
    rows: list[dict] = []
    fig_all, axes = plt.subplots(1, 3, figsize=(16, 5.2), sharey=True)

    for ax, lim in zip(axes, MAG_CUTS):
        # pool masks
        m_pure = (m_face <= lim) & (ba > Q0)                      # option 2, dusted
        m_ceil = (m_face <= lim) & (ba > Q0) & (ba <= BA_FACE_CAP)  # option 3, dusted

        # option 2: SDSS Ryden applied directly to observed b/a
        c_pure = sampler.sample(ba[m_pure], rng)
        # option 3: de-REX (b/a -> b/a/0.8), then SDSS Ryden
        c_ceil = sampler.sample(np.clip(ba[m_ceil] / BA_FACE_CAP, 0.0, 1.0), rng)
        # ad-hoc scaled reference
        c_scaled = scaled_cosi(ba[m_ceil])

        series = [
            ("ceiling-0.8 (Ryden + de-REX) + A1", c_ceil, "#e41a1c", 2.4, "-"),
            ("pure-lit (SDSS Ryden) + A1", c_pure, "#377eb8", 2.0, "--"),
            ("ad-hoc scaled + A1 (ref)", c_scaled, "#4daf4a", 1.8, ":"),
        ]
        fig, ax1 = plt.subplots(figsize=(6.8, 6.8))
        for a in (ax, ax1):
            a.plot((0, 1), (0, 1), "k--", lw=1.1, label="Uniform")
            for lab, cosi, col, lw, ls in series:
                a.plot(x, empirical_cdf(cosi, x), color=col, lw=lw, ls=ls,
                       label=f"{lab}\n  N={int(np.isfinite(cosi).sum()):,}, med={med(cosi):.3f}")
            a.set_xlim(0, 1)
            a.set_ylim(0, 1)
            a.grid(True, alpha=0.3)
            a.set_xlabel(r"$\cos(i)$  (strict $b/a>q_0$)")
        ax.set_title(f"mag limit = {lim:g}")
        ax1.set_ylabel("Cumulative distribution")
        ax1.set_title(f"LS Ryden fixed-lit (2004) + dust — mag limit={lim:g}")
        ax1.legend(loc="upper left", fontsize=8)
        fig.tight_layout()
        fig.savefig(OUT_DIR / "cdfs" / f"mag{int(lim)}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

        rows.append({
            "mag_limit": lim,
            "n_pure_a1": int(m_pure.sum()),
            "median_cosi_pure_a1": round(med(c_pure), 4),
            "n_ceil_a1": int(m_ceil.sum()),
            "median_cosi_ceiling_a1": round(med(c_ceil), 4),
            "median_cosi_scaled_a1_ref": round(med(c_scaled), 4),
        })
        print(f"    mag<={lim:g}: pure-lit+A1 med={med(c_pure):.3f} (N={int(m_pure.sum()):,})  "
              f"ceiling08+A1 med={med(c_ceil):.3f} (N={int(m_ceil.sum()):,})  "
              f"scaled+A1 ref med={med(c_scaled):.3f}")

    axes[0].set_ylabel("Cumulative distribution")
    axes[0].legend(loc="upper left", fontsize=7)
    fig_all.suptitle(
        "LS EXP null: Ryden fixed-lit (2004) — pure-lit vs ceiling-0.8 vs ad-hoc scaled "
        "(all + Unterborn A1 dust)", fontsize=12)
    fig_all.tight_layout()
    fig_all.savefig(OUT_DIR / "cdf_compare.png", dpi=300, bbox_inches="tight")
    plt.close(fig_all)

    pd.DataFrame(rows).to_csv(OUT_DIR / "summary.csv", index=False)

    (OUT_DIR / "fit_params.json").write_text(json.dumps({
        "track": "scaled_ryden_fixed",
        "shape_source": "Ryden 2004 SDSS DR1 (frozen, NOT refit to LS)",
        "params": {k: float(v) for k, v in RYDEN_SEED.to_dict().items()},
        "dust": "Unterborn A1 face-on mag cut, coeff=1.27",
        "n_model_pool": N_POOL,
        "tracks": {
            "pure_lit": "cos i ~ P_Ryden2004(cos i | b/a); no de-REX",
            "ceiling08": "cos i ~ P_Ryden2004(cos i | b/a/0.8), b/a<=0.8 (assume 0.8 face-on)",
        },
        "note": "Apply the IDENTICAL frozen sampler + same b/a treatment to FRB host b/a.",
    }, indent=2), encoding="utf-8")

    tbl = "| mag | pure-lit + A1 | ceiling-0.8 + A1 | scaled + A1 (ref) |\n"
    tbl += "|----:|--------------:|-----------------:|------------------:|\n"
    for r in rows:
        tbl += (f"| {r['mag_limit']:g} | {r['median_cosi_pure_a1']} | "
                f"{r['median_cosi_ceiling_a1']} | {r['median_cosi_scaled_a1_ref']} |\n")
    (OUT_DIR / "summary.md").write_text(
        "# scaled_ryden_fixed — Ryden (2004) frozen shapes + Unterborn A1 dust\n\n"
        "## The REX trilemma\n\n"
        "REX (round-object excision) removes near-round disks from LS Tractor EXP, so LS\n"
        "has a deficit of face-on/round galaxies. This forces a three-way choice:\n\n"
        "1. **Refit shapes to LS** -> the shape distribution absorbs all the b/a\n"
        "   structure and cos(i) collapses to ~uniform (circular; isotropy in = out).\n"
        "   See `scaled_ryden/CIRCULARITY_CHECK.md`.\n"
        "2. **Use SDSS Ryden (2004) shapes directly** (`pure_lit`) -> the REX round-\n"
        "   deficit is misread as real inclination, pushing the median too edge-on\n"
        "   (~0.41 after dust).\n"
        "3. **Use SDSS Ryden shapes + assume b/a=0.8 is face-on** (`ceiling08`;\n"
        "   de-REX via b/a -> b/a/0.8) -> curved, median ~0.5, matching ad-hoc scaled.\n\n"
        "Assuming the 0.8 ceiling is the price of directly using literature Ryden shapes\n"
        "on REX-truncated data. Per-galaxy cos(i) is drawn from P(cos i | b/a) of the\n"
        "FROZEN model (never refit); dust via the Unterborn face-on mag re-cut.\n\n"
        "## Median cos(i) (Unterborn A1 dust)\n\n" + tbl +
        "\nApply the identical frozen sampler + same b/a treatment (fit_params.json) to\n"
        "FRB host b/a and compare FRB-vs-null (two-sample KS/AD), not FRB-vs-uniform.\n",
        encoding="utf-8")
    print(f"[*] Wrote outputs -> {OUT_DIR}")


if __name__ == "__main__":
    main()
