"""Results-section figures, built exactly to the manuscript captions.

fig1_cdf   (a) 51 FRB hosts vs the HSC control with the 68% band of 10,000
           random samples of 51; (b) the 23 Bhardwaj et al. (2024) hosts with
           their published inclinations (grey), our GALFIT remeasurement (red)
           and our independent AstroPhot remeasurement (orange), with the 68%
           band for samples of 23.  Dashed line = uniform.
fig2_bins  Counts of the 51 hosts in four bins defined by the quartiles of the
           control sample's cos i; dashed line at 12.75; shaded 68% band from
           random samples of 51 control galaxies.

Also writes RESULTS_FIGURES.md with the AstroPhot - GALFIT inclination RMS for
the 23 and the 51, and the numbers plotted.

Outputs: Publication/materials/figures/fig1_cdf.{pdf,png}, fig2_bins.{pdf,png},
         RESULTS_FIGURES.md, results_fig_data.csv
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import plot_paper_inclination_figures as P  # noqa: E402  (style, loaders, helpers)

OUT = P.REPO_ROOT / "Publication" / "materials" / "figures"
SEED51 = P.REPO_ROOT / "Publication" / "materials" / "host_properties" / "seed_51.csv"
AP_C = P.ORANGE


def save(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=250, bbox_inches="tight")
    plt.close(fig)
    print("wrote", stem)


def band(pool: np.ndarray, n_sample: int, grid: np.ndarray, rng: np.random.Generator):
    cdf = np.empty((P.N_DRAWS, len(grid)))
    for k in range(P.N_DRAWS):
        cdf[k] = P.ecdf_at(pool[rng.choice(len(pool), size=n_sample, replace=False)], grid)
    return np.percentile(cdf, 16, axis=0), np.percentile(cdf, 84, axis=0)


def fig1(hosts: pd.DataFrame, hsc: np.ndarray) -> None:
    grid = np.linspace(0, 1, 500)
    frb51 = hosts["cosi_point"].to_numpy(dtype=float)
    n23 = hosts.loc[hosts["i_2024"].notna()]
    lo51, hi51 = band(hsc, 51, grid, np.random.default_rng(P.SEED))
    lo23, hi23 = band(hsc, 23, grid, np.random.default_rng(P.SEED + 1))
    hsc_cdf = P.ecdf_at(hsc, grid)

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.2), sharey=True)
    for ax, lo, hi, nlab in ((axes[0], lo51, hi51, 51), (axes[1], lo23, hi23, 23)):
        ax.plot([0, 1], [0, 1], "k--", lw=0.6, label="Uniform", zorder=1)
        ax.fill_between(grid, lo, hi, color=P.HSC_C, alpha=0.28, lw=0,
                        label=f"68% band, samples of {nlab}", zorder=2)
        ax.plot(grid, hsc_cdf, color=P.HSC_C, lw=1.0, label="Control sample", zorder=3)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel(r"$\cos i$")
        P.top_inclination(ax)
    x, y = P.step_xy(frb51)
    axes[0].step(x, y, where="post", color=P.FRB_C, lw=1.5, label="FRB hosts (N = 51)", zorder=5)
    x, y = P.step_xy(np.cos(np.radians(n23["i_2024"].to_numpy(dtype=float))))
    axes[1].step(x, y, where="post", color=P.GREY, lw=1.2, label="Bhardwaj et al. (2024), published", zorder=4)
    x, y = P.step_xy(n23["cosi_point"].to_numpy(dtype=float))
    axes[1].step(x, y, where="post", color=P.FRB_C, lw=1.5, label="Same 23, GALFIT (this work)", zorder=5)
    x, y = P.step_xy(np.cos(np.radians(n23["i_astrophot"].to_numpy(dtype=float))))
    axes[1].step(x, y, where="post", color=AP_C, lw=1.2, ls="-", label="Same 23, AstroPhot (this work)", zorder=6)
    axes[0].set_ylabel("Cumulative fraction")
    for ax, letter in zip(axes, "ab"):
        ax.legend(loc="upper left", frameon=False, fontsize=6.5)
        ax.text(0.97, 0.05, f"({letter})", transform=ax.transAxes, ha="right", va="bottom", fontsize=9)
    fig.tight_layout()
    save(fig, "fig1_cdf")


def fig2(hosts: pd.DataFrame, hsc: np.ndarray) -> dict:
    rng = np.random.default_rng(P.SEED + 4)
    edges = np.quantile(hsc, np.linspace(0, 1, 5))
    edges[0], edges[-1] = -np.inf, np.inf
    frb = hosts["cosi_point"].to_numpy(dtype=float)
    counts = np.histogram(frb, bins=edges)[0]
    draws = np.empty((P.N_DRAWS, 4), dtype=int)
    for k in range(P.N_DRAWS):
        draws[k] = np.histogram(hsc[rng.choice(len(hsc), size=51, replace=False)], bins=edges)[0]
    lo, hi = np.percentile(draws, 16, axis=0), np.percentile(draws, 84, axis=0)
    finite = edges.copy()
    finite[0], finite[-1] = 0.0, 1.0
    i_hi, i_lo = P.cosi_to_i(finite[:-1]), P.cosi_to_i(finite[1:])
    labels = [f"{b:.0f}°–{a:.0f}°" for a, b in zip(i_hi, i_lo)]  # low–high in degrees

    fig, ax = plt.subplots(figsize=(3.5, 3.0))
    x = np.arange(4)
    for k in range(4):
        ax.fill_between([k - 0.5, k + 0.5], lo[k], hi[k], color=P.HSC_C, alpha=0.25, lw=0,
                        label="68% range, samples of 51" if k == 0 else None, zorder=1)
    ax.axhline(12.75, color="k", ls="--", lw=0.7, label="Expected 12.75", zorder=2)
    ax.bar(x, counts, color=P.FRB_C, width=0.6, label="FRB hosts (N = 51)", zorder=3)
    for k, c in enumerate(counts):
        ax.text(k, c + 0.4, str(int(c)), ha="center", va="bottom", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_xlabel(r"Inclination bin (control-sample quartiles of $\cos i$)")
    ax.set_ylabel("Number of FRB hosts")
    ax.set_ylim(0, max(24, hi.max() + 4))
    ax.set_xlim(-0.6, 3.6)
    ax.legend(frameon=False, fontsize=6, loc="upper left")
    ax.tick_params(top=False)
    fig.tight_layout()
    save(fig, "fig2_bins")
    return {"edges_cosi": [float(e) for e in finite], "labels": labels, "counts": counts.tolist(),
            "lo68": lo.tolist(), "hi68": hi.tolist()}


def rms_block(name: str, d: np.ndarray) -> str:
    d = d[np.isfinite(d)]
    return (f"| {name} | {len(d)} | {np.median(d):+.2f} | {np.mean(d):+.2f} | "
            f"{np.sqrt(np.mean(d ** 2)):.2f} | {1.4826 * np.median(np.abs(d - np.median(d))):.2f} | "
            f"{int((np.abs(d) > 10).sum())} |")


def main() -> None:
    hosts = P.load_hosts()
    assert len(hosts) == 51
    seed = pd.read_csv(SEED51, dtype={"frb": str}).set_index("frb")
    hosts["i_astrophot"] = hosts["frb"].map(seed["i_astrophot_deg"])
    hosts["q_astrophot"] = hosts["frb"].map(seed["q_astrophot"])
    hsc, n_hsc = P.load_hsc_pool(mag_limit=P.MAG_LIMIT, q0=P.Q0)

    fig1(hosts, hsc)
    bins = fig2(hosts, hsc)

    n23 = hosts.loc[hosts["i_2024"].notna()].copy()
    d23 = (n23["i_astrophot"] - n23["i_deg"]).to_numpy(dtype=float)
    d51 = (hosts["i_astrophot"] - hosts["i_deg"]).to_numpy(dtype=float)
    d23_pub_g = (n23["i_deg"] - n23["i_2024"]).to_numpy(dtype=float)
    d23_pub_a = (n23["i_astrophot"] - n23["i_2024"]).to_numpy(dtype=float)

    hosts[["frb", "q", "i_deg", "q_astrophot", "i_astrophot", "i_2024", "cosi_point", "re_over_fwhm"]].to_csv(
        OUT / "results_fig_data.csv", index=False)

    md = ["# Results-section figures", "",
          "Built by `scripts/plot_results_figures.py` (seed 20240924; 10,000 control draws). Inputs: winning-leg "
          "GALFIT $q$ (`statistics/confirmed50_q.csv`), AstroPhot $q$ from the verification run "
          "(`host_properties/seed_51.csv`, column `q_astrophot`), HSC exponential-disk control pool "
          f"($r \\le 22$, $q > 0.2$; $N = {n_hsc:,}$). $\\cos^2 i = (q^2 - q_0^2)/(1 - q_0^2)$, $q_0 = 0.2$.", "",
          "## fig1_cdf", "",
          "(a) 51 hosts vs control with the 68% band of samples of 51. (b) The 23 hosts of Bhardwaj et al. (2024): "
          "published inclinations (grey), GALFIT remeasurement (red), AstroPhot remeasurement (orange), band for "
          "samples of 23. The two remeasurements lie on top of each other and inside the band; the published curve "
          "does not.", "",
          "## AstroPhot − GALFIT inclination differences", "",
          "| sample | $N$ | median | mean | RMS | $1.4826\\,$MAD | $|\\Delta i| > 10^\\circ$ |",
          "|---|---|---|---|---|---|---|",
          rms_block("AstroPhot − GALFIT, the 23 of Bhardwaj+2024", d23),
          rms_block("AstroPhot − GALFIT, all 51", d51),
          rms_block("GALFIT (this work) − published 2024", d23_pub_g),
          rms_block("AstroPhot (this work) − published 2024", d23_pub_a), ""]
    big = n23.loc[np.abs(d23) > 5, ["frb", "i_2024", "i_deg", "i_astrophot", "re_over_fwhm"]]
    if len(big):
        md += ["Hosts of the 23 with $|i_{\\rm AstroPhot} - i_{\\rm GALFIT}| > 5^\\circ$:", "",
               "| FRB | published | GALFIT | AstroPhot | $R_e$/FWHM |", "|---|---|---|---|---|"]
        for _, r in big.iterrows():
            md.append(f"| {r.frb} | {r.i_2024:.0f} | {r.i_deg:.1f} | {r.i_astrophot:.1f} | {r.re_over_fwhm:.2f} |")
        md.append("")
    md += ["## fig2_bins", "",
           "Bins are the quartiles of the control sample's $\\cos i$ (edge-on to face-on). Expected 12.75 per bin.", "",
           "| bin | $i$ range | FRB count | 68% range of control samples of 51 |", "|---|---|---|---|"]
    for k in range(4):
        md.append(f"| {k + 1} | {bins['labels'][k]} | {bins['counts'][k]} | {bins['lo68'][k]:.0f}–{bins['hi68'][k]:.0f} |")
    md.append("")
    rms23 = float(np.sqrt(np.mean(d23[np.isfinite(d23)] ** 2)))
    tex = rf"""\begin{{figure*}}[t]
\centering
\includegraphics[width=\textwidth]{{fig1_cdf.pdf}}
\caption{{CDFs of $\cos i$. (a) The 51 FRB hosts (red) and the control sample
(blue), with the 68\% bound of 10{{,}}000 random samples of 51 control galaxies
(shaded). (b) The 23 hosts of \citet{{Bhardwaj2024}}, using their published
inclinations (grey), our \textsc{{galfit}} remeasurement (red) and an independent
\textsc{{AstroPhot}} remeasurement of the same images (orange; RMS difference
between the two codes ${rms23:.1f}^\circ$), with the 68\% bound for samples of 23.
The dashed line shows a uniform distribution.}}
\label{{fig:cdf}}
\end{{figure*}}

\begin{{figure}}[t]
\centering
\includegraphics[width=\columnwidth]{{fig2_bins.pdf}}
\caption{{Number of FRB hosts in each of four inclination bins, defined by the
quartiles of the control sample's $\cos i$ distribution; the bin ranges are
marked on the horizontal axis. The dashed line marks the expected 12.75 hosts
per bin, and the shaded band the 68\% range from random samples of 51 control
galaxies.}}
\label{{fig:bins}}
\end{{figure}}
"""
    (OUT / "results_figures.tex").write_text(tex, encoding="utf-8")
    md += ["## LaTeX", "", "Figure environments with the updated panel-(b) caption: `results_figures.tex`.", ""]
    (OUT / "RESULTS_FIGURES.md").write_text("\n".join(md), encoding="utf-8")
    print("\n".join(md))


if __name__ == "__main__":
    main()
