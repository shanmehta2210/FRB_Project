"""Appendix figures, styled like Figure 1 (plot_results_figures.py / plot_paper_inclination_figures.py).

figA_control_hist      control-pool cos i histogram vs uniform          (3.4 x 2.8 in)
figB_gini_m20          Gini-M20 of the hosts   , Lotz+2008 boundaries   (3.4 x 2.8 in)
figA_control_variants  CDFs of the control variants of Table A1          (7.1 x 2.6 in, 3 panels)
figC_q0                CDFs for q0 = 0.10, 0.15, 0.20, 0.25 with 68% bands (2 x 2)

p_m annotations are read from Publication/materials/appendix/robustness/robustness_results.json
(written by scripts/appendix_robustness.py), so the figures always carry the latest run.

Outputs: Publication/materials/appendix/robustness/{pdf,png} + APPENDIX_FIGURES.md
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import plot_paper_inclination_figures as P  # noqa: E402  (rcParams, colours, helpers)
import appendix_robustness as R              # noqa: E402  (pool loaders)
from null_catalog_utils import hubble_cosi_from_ba  # noqa: E402

OUT = R.OUT
RES = json.loads((OUT / "robustness_results.json").read_text(encoding="utf-8"))
PANEL = (3.4, 2.8)
GRID = np.linspace(0.0, 1.0, 500)

PANEL_MARGINS = dict(left=0.165, right=0.975, bottom=0.165, top=0.845)  # identical box for the two side-by-side panels


def save(fig: plt.Figure, stem: str, tight: bool = True) -> None:
    kw = dict(bbox_inches="tight", pad_inches=0.02) if tight else {}
    fig.savefig(OUT / f"{stem}.pdf", **kw)
    fig.savefig(OUT / f"{stem}.png", dpi=200, **kw)
    plt.close(fig)
    print("wrote", stem)


def band(pool: np.ndarray, n_sample: int, rng: np.random.Generator):
    cdf = np.empty((P.N_DRAWS, len(GRID)))
    for k in range(P.N_DRAWS):
        cdf[k] = P.ecdf_at(pool[rng.choice(len(pool), size=n_sample, replace=False)], GRID)
    return np.percentile(cdf, 16, axis=0), np.percentile(cdf, 84, axis=0)


def pm(block: str, key: str, value) -> float:
    for row in RES[block]:
        if row.get(key) == value:
            return row["ad_pm"]
    raise KeyError((block, key, value))


# ----------------------------------------------------------------------------
def fig_control_hist(pool: np.ndarray) -> dict:
    fig, ax = plt.subplots(figsize=PANEL)
    counts, edges, _ = ax.hist(pool, bins=20, range=(0, 1), density=True, color=P.HSC_C, alpha=0.5, lw=0,
                               label="Control sample")
    ax.hist(pool, bins=20, range=(0, 1), density=True, histtype="step", color="#0B4F7E", lw=0.8)
    ax.axhline(1.0, color="k", ls="--", lw=0.8, label="Random orientation")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1.15 * counts.max())
    ax.set_xlabel(r"$\cos i$")
    ax.set_ylabel("Probability density")
    ax.tick_params(right=True)
    P.top_inclination(ax)
    ax.legend(loc="upper left", frameon=False, fontsize=7)
    fig.subplots_adjust(**PANEL_MARGINS)
    save(fig, "figA_control_hist", tight=False)
    return {"N": int(len(pool)), "bins": 20, "max_density": float(counts.max()),
            "median_i": float(np.degrees(np.arccos(np.median(pool)))), "frac_i_gt60": float(np.mean(pool < 0.5))}


# ----------------------------------------------------------------------------
def fig_gini_m20(n: int) -> dict:
    d = pd.read_csv(OUT / f"gini_m20_{n}.csv", dtype={"frb": str})
    assert len(d) == n
    rel, unr = d[d.reliable], d[~d.reliable]
    xlo, xhi = -0.5, -3.0
    ylo, yhi = min(0.30, d.sm_gini.min() - 0.02), max(0.70, d.sm_gini.max() + 0.03)
    m_meet = (0.33 - 0.80) / (0.14 + 0.14)        # -1.68
    g_meet = -0.14 * m_meet + 0.33                 # 0.565

    fig, ax = plt.subplots(figsize=PANEL)
    m = np.linspace(xhi, xlo, 200)
    ax.plot(m, -0.14 * m + 0.33, "k-", lw=0.8, zorder=1)
    m2 = np.linspace(xhi, m_meet, 100)
    ax.plot(m2, 0.14 * m2 + 0.80, "k-", lw=0.8, zorder=1)
    ax.text(-0.6, 0.66, "Mergers", color="0.45", fontsize=7, style="italic", ha="left")
    ax.text(-2.92, 0.45, "E/S0/Sa", color="0.45", fontsize=7, style="italic", ha="right")
    ax.text(-0.6, 0.315, "Sb–Irr", color="0.45", fontsize=7, style="italic", ha="left")
    ax.scatter(rel.sm_m20, rel.sm_gini, s=18, color=P.FRB_C, edgecolor=P.FRB_C, lw=0.6, label="Reliable", zorder=3)
    ax.scatter(unr.sm_m20, unr.sm_gini, s=18, facecolor="white", edgecolor=P.FRB_C, lw=0.7, label="Unreliable", zorder=3)
    offsets = {"20190608B": (4, 4), "20220319D": (4, -8), "20231120A": (-5, -9), "20220509G": (4, -8)}
    for r in d[d.lotz2008 != "late"].itertuples():
        dx, dy = offsets.get(r.frb, (4, 3))
        ax.annotate(r.frb, (r.sm_m20, r.sm_gini), xytext=(dx, dy), textcoords="offset points", fontsize=6,
                    ha="left" if dx > 0 else "right")
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ylo, yhi)
    ax.set_xlabel(r"$M_{20}$")
    ax.set_ylabel("Gini")
    ax.tick_params(top=True, right=True)
    ax.legend(loc="lower right", frameon=False, fontsize=7)
    fig.subplots_adjust(**PANEL_MARGINS)
    save(fig, "figB_gini_m20", tight=False)
    return {"n_reliable": int(len(rel)), "n_unreliable": int(len(unr)), "meet_M20": m_meet, "meet_G": g_meet,
            "labelled": d.loc[d.lotz2008 != "late", "frb"].tolist(), "xlim": [xlo, xhi], "ylim": [ylo, yhi]}


# ----------------------------------------------------------------------------
def _cdf_panel(ax, frb: np.ndarray, pool: np.ndarray, p: float, title: str, seed: int, *, top: bool, title_pad: float) -> dict:
    """One control-variant panel: band from this pool (size = len(frb)), control CDF, FRB step, AD p_m."""
    lo, hi = band(pool, len(frb), np.random.default_rng(seed))
    ax.plot([0, 1], [0, 1], "k--", lw=0.6, zorder=1)
    ax.fill_between(GRID, lo, hi, color=P.HSC_C, alpha=0.28, lw=0, label=f"68% band, samples of {len(frb)}", zorder=2)
    ax.plot(GRID, P.ecdf_at(pool, GRID), color=P.HSC_C, lw=1.0, label=f"Control sample (N = {len(pool):,})", zorder=3)
    x, y = P.step_xy(frb)
    ax.step(x, y, where="post", color=P.FRB_C, lw=1.5, label=f"FRB hosts (N = {len(frb)})", zorder=5)
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(0.97, 0.05, rf"AD $p_{{\rm m}} = {p:.2f}$", transform=ax.transAxes, ha="right", va="bottom", fontsize=8)
    ax.tick_params(right=True)
    if top:  # title sits just above the inclination-axis label
        P.top_inclination(ax).set_title(title, fontsize=8, pad=title_pad)
    else:
        ax.tick_params(top=True)
        ax.set_title(title, fontsize=8, pad=title_pad)
    ax.legend(loc="upper left", frameon=False, fontsize=6.5)
    return {"panel": title, "N_frb": int(len(frb)), "N_pool": int(len(pool)), "ad_pm": p}


def _finish_grid(fig, axes, stem: str) -> None:
    axes = np.atleast_2d(axes)
    for ax in axes[-1]:
        ax.set_xlabel(r"$\cos i$")
    for ax in axes[:, 0]:
        ax.set_ylabel("Cumulative fraction")
    fig.tight_layout(h_pad=0.8, w_pad=0.6)
    save(fig, stem)


def fig_ctrl_sersic(parent: pd.DataFrame, frb: np.ndarray) -> list:
    items = [("$0.4 < n < 2.0$", dict(n_lo=0.4, n_hi=2.0), pm("A2", "n_hi", 2.0)),
             ("$0.4 < n < 2.5$", dict(n_lo=0.4, n_hi=2.5), pm("A2", "label", "$0.4 < n < 2.5$")),
             ("$n < 2.5$", dict(n_lo=0.0, n_hi=2.5), pm("A2", "label", "$n < 2.5$")),
             ("No morphology cut", dict(n_lo=0.0, n_hi=10.0), pm("A2", "n_hi", 10.0))]
    fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.6), sharex=True, sharey=True)
    log = []
    for k, (ax, (title, kw, p)) in enumerate(zip(axes.ravel(), items)):
        pool = R.cosi_from_q(R.select(parent, **kw)["ba"].to_numpy(), P.Q0)
        log.append(_cdf_panel(ax, frb, pool, p, title, P.SEED + 30 + k, top=k < 2, title_pad=4))
    _finish_grid(fig, axes, "figA_ctrl_sersic")
    return log


def fig_ctrl_mag(parent: pd.DataFrame, hosts: pd.DataFrame) -> list:
    frb = hosts["cosi_point"].to_numpy(dtype=float)
    mag = hosts["mag"].to_numpy(dtype=float)
    items = []
    for m in (21.0, 20.0):
        items.append((rf"$r \leq {m:.0f}$ (control only)", m, frb, pm("A3", "label", f"control $r \\le {m:.0f}$, all {len(hosts)} hosts")))
        items.append((rf"$r \leq {m:.0f}$ (both samples)", m, frb[mag <= m], pm("A3", "label", f"both samples $r \\le {m:.0f}$")))
    fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.6), sharex=True, sharey=True)
    log = []
    for k, (ax, (title, m, f, p)) in enumerate(zip(axes.ravel(), items)):
        pool = R.cosi_from_q(R.select(parent, mag=m)["ba"].to_numpy(), P.Q0)
        log.append(_cdf_panel(ax, f, pool, p, title, P.SEED + 40 + k, top=k < 2, title_pad=4))
    _finish_grid(fig, axes, "figA_ctrl_mag")
    return log


def fig_ctrl_color(parent: pd.DataFrame, frb: np.ndarray) -> list:
    gi = parent["gmag"] - parent["imag"]
    items = [("Star-forming", parent["quiescent_flag"] == 0, pm("A5", "label", "star-forming only (`quiescent_flag = 0`)")),
             ("$g - i < 1.2$", gi < 1.2, pm("A5", "label", "blue, $g - i < 1.2$")),
             ("$g - i \\geq 1.2$", gi >= 1.2, pm("A5", "label", "red, $g - i \\ge 1.2$"))]
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.6), sharex=True, sharey=True)
    log = []
    for k, (ax, (title, extra, p)) in enumerate(zip(axes, items)):
        pool = R.cosi_from_q(R.select(parent, extra=extra)["ba"].to_numpy(), P.Q0)
        log.append(_cdf_panel(ax, frb, pool, p, title, P.SEED + 50 + k, top=True, title_pad=4))
    _finish_grid(fig, axes, "figA_ctrl_color")
    return log


# ----------------------------------------------------------------------------
def fig_q0(hosts: pd.DataFrame, parent: pd.DataFrame) -> list:
    par_q = parent.loc[(parent["n_sersic"] > 0.4) & (parent["n_sersic"] < 1.5) & (parent["rmag"] <= 22.0), "ba"].to_numpy(dtype=float)
    qser = P.magnitude_q(hosts)
    q = qser.to_numpy(dtype=float)
    fig, axes = plt.subplots(2, 2, figsize=(7.1, 6.2), sharex=True, sharey=True)
    log = []
    for k, (ax, q0) in enumerate(zip(axes.ravel(), (0.10, 0.15, 0.20, 0.25))):
        frb_c = np.array([hubble_cosi_from_ba(float(v), q0=q0) for v in q[q > q0]])
        hsc_c = np.array([hubble_cosi_from_ba(float(v), q0=q0) for v in par_q[par_q > q0]])
        lo, hi = band(hsc_c, len(frb_c), np.random.default_rng(P.SEED + 20 + k))
        p = pm("C1", "q0", q0)
        ax.plot([0, 1], [0, 1], "k--", lw=0.6, zorder=1)
        ax.fill_between(GRID, lo, hi, color=P.HSC_C, alpha=0.28, lw=0, label=f"68% band, samples of {len(frb_c)}", zorder=2)
        ax.plot(GRID, P.ecdf_at(hsc_c, GRID), color=P.HSC_C, lw=1.0, label=f"Control sample (N = {len(hsc_c):,})", zorder=3)
        x, y = P.step_xy(frb_c)
        ax.step(x, y, where="post", color=P.FRB_C, lw=1.5, label=f"FRB hosts (N = {len(frb_c)})", zorder=5)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        ax.text(0.97, 0.05, rf"$q_0 = {q0:.2f}$" + "\n" + rf"AD $p_{{\rm m}} = {p:.2f}$", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=8)
        ax.tick_params(right=True)
        P.top_inclination(ax)
        ax.legend(loc="upper left", frameon=False, fontsize=6.5)
        log.append({"q0": q0, "N_frb": int(len(frb_c)), "N_pool": int(len(hsc_c)), "ad_pm": p})
    for ax in axes[1]:
        ax.set_xlabel(r"$\cos i$")
    for ax in axes[:, 0]:
        ax.set_ylabel("Cumulative fraction")
    fig.tight_layout(h_pad=1.2)
    save(fig, "figC_q0")
    return log


# ----------------------------------------------------------------------------
def main() -> None:
    hosts = P.load_hosts(); N = len(hosts)
    frb = hosts["cosi_point"].to_numpy(dtype=float)
    parent = R.load_parent()
    base = R.cosi_from_q(R.select(parent)["ba"].to_numpy(), P.Q0)
    assert len(base) == 24_450

    a = fig_control_hist(base)
    b = fig_gini_m20(N)
    vs = fig_ctrl_sersic(parent, frb)
    vm = fig_ctrl_mag(parent, hosts)
    vc = fig_ctrl_color(parent, frb)
    c = fig_q0(hosts, parent)
    for old in ("figA_control_variants.pdf", "figA_control_variants.png"):  # replaced by the three figures above
        (OUT / old).unlink(missing_ok=True)

    md = ["# Appendix figures", "",
          "Source: `scripts/appendix_figures.py`. Style inherited from Figure 1 (`plot_paper_inclination_figures.py` rcParams: serif/STIX, "
          "9 pt axis labels, 8 pt ticks, inward ticks; control `#1F77B4`, FRB hosts `#D62728`, reference lines black). "
          "The two side-by-side panels are saved at exactly 3.4 × 2.8 in with identical axes margins (`subplots_adjust`, no bbox cropping) so they line up in one row; "
          "the multi-panel figures use `bbox_inches='tight', pad_inches=0.02`. PNGs at 200 dpi. "
          f"$p_{{\\rm m}}$ annotations are read from `robustness_results.json` (run of {RES['n_draws']:,} draws).", "",
          "## figA_control_hist (3.4 × 2.8 in, left panel)", "",
          f"Histogram of $\\cos i$ for the {a['N']:,} control galaxies (adopted pool: `goodfits_flag = 1`, $r \\le 22$, $0.4 < n < 1.5$, $q > 0.2$, $q_0 = 0.2$), "
          f"{a['bins']} bins of width 0.05, `density=True` (uniform = 1). Filled blue ($\\alpha = 0.5$) with a darker blue step outline; dashed black line at 1 labelled "
          f"\"Random orientation\". Top axis: inclination at 90, 75, 60, 45, 30, 0 deg via `top_inclination`. y from 0 to 1.15 × tallest bar ({a['max_density']:.2f}). "
          f"No median line, quartile edges or FRB hosts. Median $i = {a['median_i']:.1f}^\\circ$, $f(i > 60^\\circ) = {100*a['frac_i_gt60']:.1f}\\%$ go in the caption.", "",
          "## figB_gini_m20 (3.4 × 2.8 in, right panel)", "",
          f"$G$ and $M_{{20}}$ from `gini_m20_{N}.csv`. Reliable = `flag = 0`, $\\langle S/N\\rangle_{{\\rm pix}} \\ge 2.5$, $R_e/{{\\rm FWHM}} \\ge 1.0$: "
          f"{b['n_reliable']} filled red circles; {b['n_unreliable']} unreliable as open red circles (s = 18). $x$ reversed, {b['xlim'][0]} to {b['xlim'][1]}; "
          f"$y$ {b['ylim'][0]:.2f} to {b['ylim'][1]:.2f}. Lotz et al. (2008) boundaries (eqs. 4–6 of that paper, checked): merger line $G = -0.14\\,M_{{20}} + 0.33$ across the full range; "
          f"E/S0/Sa–Sb/Irr line $G = 0.14\\,M_{{20}} + 0.80$ drawn only for $M_{{20}} \\le {b['meet_M20']:.2f}$ (the lines meet at $G = {b['meet_G']:.3f}$). "
          f"Grey italic region labels. Only the four hosts outside Sb–Irr are labelled: {', '.join(b['labelled'])}. Legend lower right (fewest points).", "",
          "## Control-variant CDFs: figA_ctrl_sersic, figA_ctrl_mag, figA_ctrl_color", "",
          "These replace `figA_control_variants`. One control variant per panel, in the style of `figC_q0`: uniform diagonal (dashed black, lw 0.6); "
          "the 68% band of 10,000 samples drawn from *that panel's* control pool, each of that panel's FRB sample size, on the 500-point $\\cos i$ grid (light blue); "
          "that variant's control CDF (solid blue, lw 1.0); the FRB step CDF (red, lw 1.5, `where='post'`, legend gives the panel's $N$); "
          f"\"AD $p_{{\\rm m}}$\" lower right from `robustness_results.json` (tables A2, A3, A5; the adopted configuration reproduces the main-text run, seed 42, AD $p_{{\\rm m}} = {pm('A2', 'n_hi', 1.5):.2f}$); "
          "short title above. `sharex`/`sharey`; $y$ label on the left column, $x$ label on the bottom row, top inclination axis on the top row only.", "",
          "`figA_ctrl_sersic` (7.1 × 5.6 in, 2 × 2): $0.4 < n < 2.0$ | $0.4 < n < 2.5$ / $n < 2.5$ | no morphology cut.", "",
          f"`figA_ctrl_mag` (7.1 × 5.6 in, 2 × 2): left column the limit on the control only (all {N} hosts), right column on both samples "
          f"(the band then uses samples of {vm[1]['N_frb']} and {vm[3]['N_frb']}); top row $r \\le 21$, bottom row $r \\le 20$.", "",
          "`figA_ctrl_color` (7.1 × 2.6 in, 1 × 3): star-forming (`quiescent_flag = 0`) | $g - i < 1.2$ | $g - i \\ge 1.2$.", "",
          "| figure | panel | $N_{\\rm FRB}$ | $N_{\\rm pool}$ | AD $p_{\\rm m}$ |", "|---|---|---|---|---|"]
    for stem, rows in (("figA_ctrl_sersic", vs), ("figA_ctrl_mag", vm), ("figA_ctrl_color", vc)):
        md += [f"| {stem} | {r['panel']} | {r['N_frb']} | {r['N_pool']:,} | {r['ad_pm']:.2f} |" for r in rows]
    md += ["", "## figC_q0 (7.1 × 6.2 in, 2 × 2)", "",
           "Revision of the old `figA_q0`: per panel the $q > q_0$ cut is applied to both samples and $\\cos i$ recomputed with that $q_0$. Added (1) the 68% band of 10,000 draws from "
           "that panel's control pool, each of that panel's FRB sample size; (2) the AD $p_{\\rm m}$ from the C1 table, lower right; (3) the top inclination axis; legend entry renamed "
           "\"Control sample\". Panel title replaced by the in-panel $q_0$ label.", "",
           "| $q_0$ | $N_{\\rm FRB}$ | $N_{\\rm pool}$ | AD $p_{\\rm m}$ |", "|---|---|---|---|"]
    md += [f"| {r['q0']:.2f} | {r['N_frb']} | {r['N_pool']:,} | {r['ad_pm']:.2f} |" for r in c]
    md += ["", "LaTeX placement (two panels in one row):", "", "```latex",
           "\\begin{figure*}", "\\centering",
           "\\includegraphics[width=3.4in]{figA_control_hist.pdf}\\hfill", "\\includegraphics[width=3.4in]{figB_gini_m20.pdf}",
           "\\caption{...}", "\\end{figure*}", "```", "",
           "Control-variant figures (full width, one per float):", "", "```latex",
           "\\begin{figure*}\\centering\\includegraphics[width=7.1in]{figA_ctrl_sersic.pdf}\\caption{...}\\end{figure*}",
           "\\begin{figure*}\\centering\\includegraphics[width=7.1in]{figA_ctrl_mag.pdf}\\caption{...}\\end{figure*}",
           "\\begin{figure*}\\centering\\includegraphics[width=7.1in]{figA_ctrl_color.pdf}\\caption{...}\\end{figure*}",
           "```", ""]
    (OUT / "APPENDIX_FIGURES.md").write_text("\n".join(md), encoding="utf-8")
    print("wrote APPENDIX_FIGURES.md")


if __name__ == "__main__":
    main()
