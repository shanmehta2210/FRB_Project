"""Figure 3: the DM_host--inclination analysis of Marnoch et al. (2026) repeated on our hosts.

(a) Rest-frame DM_host (their definition) versus axis ratio q with an LTSFIT line, 1 and 2.6 sigma bands
    and clipped points (their Figure 8); their best-fitting relation dashed, shared hosts outlined.
(b) Spearman r_S with 95% Fisher intervals: their published values, DM_host for our sample and their
    subsamples, and the propagation-related burst properties plus old-star extinction.
(c) Scatter about the Macquart relation before and after subtracting a host term (their Figure 13),
    for their Table 1 and for our hosts.

Inputs: scripts/marnoch_replication.py outputs and plots/plots_null/v2/host_correlations/correlations.csv
Output: Publication/materials/figures/fig3_correlations.{pdf,png}
"""
from __future__ import annotations

import contextlib
import io
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_paper_inclination_figures as P  # noqa: E402  (style)
import marnoch_replication as M  # noqa: E402
from ltsfit.lts_linefit import lts_linefit  # noqa: E402

OUT = P.REPO_ROOT / "Publication" / "materials" / "figures"
GREY = "0.45"

PROPS = [("log_abs_rm", r"$|{\rm RM}|$"), ("log_tau", "Scattering time"),
         ("lin_pol", "Linear polarization"), ("av_old", r"$A_V$, old stars")]


def panel_a(ax, h: pd.DataFrame, shared: set[str]) -> None:
    x, y = h.q.to_numpy(float), h.dm_mac_rest.to_numpy(float)
    sy = (h.sig_dm * (1 + h.z)).to_numpy(float)
    p, scale = M.lts(x, y, h.sig_q.to_numpy(float), sy, 0.5, full=True)
    print("panel a LTS: intercept, slope", p.ab, "+/-", p.ab_err, "sig_int", p.sig_int, "clipped", (~p.mask).sum(), "sy scale", scale)
    a, b = p.ab
    xx = np.linspace(0.2, 1.0, 50)
    yy = a + b * (xx - 0.5)
    for k, alpha in ((2.6, 0.12), (1.0, 0.22)):
        ax.fill_between(xx, yy - k * p.sig_int, yy + k * p.sig_int, color=P.FRB_C, alpha=alpha, lw=0, zorder=1)
    ax.plot(xx, yy, color=P.FRB_C, lw=1.0, zorder=2, label="LTS fit, this work")
    ax.plot(xx, -309.3 * (xx - 0.5) + 155.9, "k--", lw=0.9, zorder=2, label="Marnoch et al. (2026)")
    sh = h.frb.isin(shared).to_numpy()
    keep = p.mask
    for sel, kw in ((keep & ~sh, dict(color=P.FRB_C, lw=0)), (keep & sh, dict(color=P.FRB_C, edgecolor="k", lw=0.7)),
                    (~keep, dict(facecolor="white", edgecolor=P.FRB_C, lw=0.8))):
        ax.scatter(x[sel], y[sel], s=14, zorder=4, **kw)
    rp, pp = stats.pearsonr(x, y)
    rs, ps = stats.spearmanr(x, y)
    ax.text(0.97, 0.97, f"$N = {len(x)}$\n$r_P = {rp:.2f}$, $p = {pp:.2f}$\n$r_S = {rs:.2f}$, $p = {ps:.2f}$",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.5)
    ax.set_xlim(0.2, 1.0)
    ax.set_ylim(-230, 980)
    ax.set_xlabel(r"Axis ratio $q$")
    ax.set_ylabel(r"DM$_{\rm host}$ (pc cm$^{-3}$, rest frame)")
    ax.legend(loc="upper left", frameon=False, fontsize=6, bbox_to_anchor=(0.0, 0.86))
    ax.text(0.03, 0.03, "(a)", transform=ax.transAxes, fontsize=9)


def panel_b(ax, R: pd.DataFrame, corr: pd.DataFrame) -> None:
    def fr(r, n):
        lo, hi = M.fisher_ci(r, n)
        return r, lo, hi, n
    rows = [("Marnoch et al., all", fr(-0.437, 34), GREY, True),
            (r"Marnoch et al., $q > 0.2$", fr(-0.358, 30), GREY, False),
            ("Marnoch et al., shared", fr(-0.19, 20), GREY, False)]
    for sub, lab in (("all", "This work, all"), ("m_r < 21", r"$m_r < 21$"), ("n < 2", r"$n < 2$"),
                     ("bulge fits removed", "Bulge fits removed"), ("shared with Marnoch", "Shared")):
        r = R[(R.subsample == sub) & (R.y == "dm_mac_rest") & (R.x == "q")].iloc[0]
        rows.append((lab, (r.r_S, r.rS_lo, r.rS_hi, int(r.N)), P.FRB_C, r.p_S < 0.05))
    c = corr.set_index("property")
    for k, lab in PROPS:
        r, n = float(c.loc[k, "spearman"]), int(c.loc[k, "N"])
        # cos i rises as q rises, so r_S with cos i is on the same scale as r_S with q
        rows.append((lab, fr(r, n), P.HSC_C, float(c.loc[k, "spearman_p"]) < 0.05))
    y = np.arange(len(rows))[::-1].astype(float)
    y[3:] -= 0.6
    y[8:] -= 0.6
    for yy, (lab, (r, lo, hi, n), col, filled) in zip(y, rows):
        ax.errorbar(r, yy, xerr=[[r - lo], [hi - r]], fmt="o", ms=3.5, color=col, mfc=col if filled else "white",
                    lw=0.9, capsize=1.5, zorder=3)
        ax.text(1.03, yy, f"{n}", transform=ax.get_yaxis_transform(), va="center", fontsize=6)
    ax.axvline(0, color="0.7", lw=0.5, zorder=1)
    ax.axvline(-0.44, color="k", ls="--", lw=0.8, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels([r[0] for r in rows], fontsize=6.5)
    ax.text(1.03, y.max() + 0.9, "$N$", transform=ax.get_yaxis_transform(), va="bottom", fontsize=6)
    ax.set_ylim(y.min() - 0.8, y.max() + 1.6)
    ax.set_xlim(-1, 1)
    ax.set_xlabel(r"Spearman $r_S$ with $q$ (or $\cos i$)")
    ax.text(-0.98, y[3] + 0.55, r"DM$_{\rm host}$", fontsize=6, color="0.3", va="center")
    ax.text(-0.98, y[8] + 0.55, "Other properties", fontsize=6, color="0.3", va="center")
    ax.tick_params(axis="y", length=0)
    ax.text(0.03, 0.03, "(b)", transform=ax.transAxes, fontsize=9)


def panel_c(ax, mac: pd.DataFrame) -> None:
    corrs = [("none", "No correction"), ("constant, best-fitting mean (no q)", r"Constant DM$_{\rm host}$"),
             ("Marnoch model DM(q)", r"Marnoch et al. DM$_{\rm host}(q)$")]
    samples = [(mac["sample"].unique()[0], "Marnoch\net al."), (mac["sample"].unique()[1], "This\nwork")]
    cols = ["0.75", P.HSC_C, "k"]
    w = 0.26
    for j, (s, slab) in enumerate(samples):
        for k, ((key, lab), col) in enumerate(zip(corrs, cols)):
            v = float(mac[(mac["sample"] == s) & (mac.correction == key)].rmse.iloc[0])
            ax.bar(j + (k - 1) * w, v, w * 0.92, color=col, label=lab if j == 0 else None)
            ax.text(j + (k - 1) * w, v + 3, f"{v:.0f}", ha="center", va="bottom", fontsize=6)
    ax.set_xticks([0, 1])
    ax.set_xticklabels([f"{lab}\n($N = {int(mac[mac['sample'] == s].N.iloc[0])}$)" for s, lab in samples], fontsize=7)
    ax.set_ylabel(r"RMSE about $\langle{\rm DM}_{\rm cosmic}\rangle(z)$ (pc cm$^{-3}$)")
    ax.set_ylim(0, 285)
    ax.legend(loc="upper left", frameon=False, fontsize=6, ncol=1, handlelength=1.0)
    ax.text(0.03, 0.03, "(c)", transform=ax.transAxes, fontsize=9, color="w")


def main() -> None:
    h = M.load()
    R = pd.read_csv(M.OUT / "replication_correlations.csv")
    mac = pd.read_csv(M.OUT / "replication_macquart.csv")
    corr = pd.read_csv(M.CORR / "correlations.csv")
    shared = set(pd.read_csv(M.OUT / "marnoch2026_table1.csv", dtype={"frb": str}).frb) & set(h.frb)
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.9), gridspec_kw={"width_ratios": [1.0, 1.05, 0.75], "wspace": 0.62})
    panel_a(axes[0], h, shared)
    panel_b(axes[1], R, corr)
    panel_c(axes[2], mac)
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"fig3_correlations.{ext}", dpi=250, bbox_inches="tight")
    plt.close(fig)
    print("wrote fig3_correlations; shared hosts:", len(shared))


if __name__ == "__main__":
    main()
