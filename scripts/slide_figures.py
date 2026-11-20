"""Slide-ready figures.

disk_gallery_{1,2,3}.png   colour images of the 51 disk hosts, 18 per 16:9 page, with survey and disk evidence
inclination_vs_control.png FRB host cos i CDF against the HSC disk control
dm_host_vs_inclination.png, rm_vs_inclination.png, scattering_vs_inclination.png, dm_host_vs_stellar_mass.png

Output: Publication/materials/slides/
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import appendix_robustness as A  # noqa: E402
import host_gallery as G  # noqa: E402
import plot_paper_inclination_figures as P  # noqa: E402
import robust_correlations as R  # noqa: E402

OUT = ROOT / "Publication" / "materials" / "slides"
MOR = ROOT / "Publication" / "materials" / "appendix" / "morphology"
PER_PAGE, NCOL = 18, 6
RED, BLUE, GREY = "#C0392B", "#2E86C1", "0.55"
CLASS_ORDER = {"reliable": 0, "lit": 1, "profile": 2}

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"], "mathtext.fontset": "dejavusans",
    "font.size": 15, "axes.labelsize": 17, "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 13,
    "axes.linewidth": 1.1, "xtick.direction": "out", "ytick.direction": "out", "xtick.top": False,
    "ytick.right": False, "xtick.major.size": 5, "ytick.major.size": 5, "axes.spines.top": False,
    "axes.spines.right": False,
})


# ---------------------------------------------------------------- disk gallery
def lit_label(note: str) -> str:
    m = re.search(r"\(lit\. n=([\d.]+)\)", note)
    note = f"Sérsic n = {m.group(1)}" if m else note.replace("late type", "late-type")
    return f"Literature: {note}"


def disk_hosts() -> pd.DataFrame:
    seed = pd.read_csv(ROOT / "Publication/materials/host_properties/seed_56.csv", dtype={"frb": str}).set_index("frb")
    ev = pd.read_csv(MOR / "disk_evidence.csv", dtype={"frb": str}).set_index("frb")
    wise = pd.read_csv(MOR / "wise_colours.csv", dtype={"frb": str}).set_index("frb").wise_class
    rows = []
    for frb in seed.index:
        cls = ev.cls.get(frb, "reliable")
        if cls not in CLASS_ORDER:
            continue
        if cls == "reliable":
            why = "Gini–M20: Sb–Irr"
        elif cls == "lit":
            why = lit_label(ev.loc[frb, "class"])
        else:
            why = f"Sérsic n = {seed.loc[frb, 'n_galfit']:.1f} (disk-like)"
        w = wise.get(frb, "")
        rows.append(dict(frb=frb, z=seed.loc[frb, "z"], cls=cls, why=why,
                         wise="WISE: disk colours" if w in ("star-forming disk", "intermediate disk") else ""))
    return pd.DataFrame(rows).assign(o=lambda d: d.cls.map(CLASS_ORDER)).sort_values(["o", "z"]).reset_index(drop=True)


def colour_tile(frb: str) -> tuple[np.ndarray | None, str]:
    win = pd.read_csv(G.WIN, dtype={"frb": str}).set_index("frb")
    ra, dec, _ = G.host_radec(frb)
    (data, _, _), _, _ = G.stamps(frb, win.loc[frb, "wkdir"])
    hdr = G.fits.getheader(ROOT / "large_cutouts" / f"{frb}_flux.fits")
    ps = abs(float(hdr.get("CD1_1") or hdr.get("CDELT1"))) * 3600
    return G.colour_image(frb, ra, dec, max(data.shape) * ps)


def gallery(d: pd.DataFrame) -> None:
    for page, k in enumerate(range(0, len(d), PER_PAGE), 1):
        sub = d.iloc[k:k + PER_PAGE]
        nrow = int(np.ceil(len(sub) / NCOL))
        fig, axes = plt.subplots(nrow, NCOL, figsize=(13.33, 2.2 * nrow), facecolor="black")
        for ax in axes.ravel():
            ax.set_facecolor("black")
            ax.axis("off")
        for ax, r in zip(axes.ravel(), sub.itertuples()):
            img, survey = colour_tile(r.frb)
            if img is not None:
                ax.imshow(img, interpolation="lanczos")
            ax.text(0.04, 0.96, f"FRB {r.frb}", transform=ax.transAxes, color="w", fontsize=10.5, weight="bold", va="top")
            if np.isfinite(r.z):
                ax.text(0.04, 0.84, f"z = {r.z:.3f}", transform=ax.transAxes, color="w", fontsize=9, va="top")
            ax.text(0.96, 0.96, survey, transform=ax.transAxes, color="0.85", fontsize=8.5, va="top", ha="right")
            ax.text(0.04, 0.04, r.why + (f"\n{r.wise}" if r.wise else ""), transform=ax.transAxes, color="#FFD54F",
                    fontsize=8.5, weight="bold", va="bottom", linespacing=1.15)
        fig.subplots_adjust(0.003, 0.003, 0.997, 0.997, 0.02, 0.02)
        fig.savefig(OUT / f"disk_gallery_{page}.png", dpi=220, facecolor="black")
        plt.close(fig)
        print(f"disk_gallery_{page}: {len(sub)} hosts", flush=True)


# ---------------------------------------------------------------- inclination vs control
def inclination_vs_control() -> None:
    hosts = P.load_hosts()
    frb = hosts["cosi_point"].to_numpy(float)
    pool = A.cosi_from_q(A.select(A.load_parent())["ba"].to_numpy(), A.Q0)
    grid = np.linspace(0, 1, 500)
    lo, hi, _ = P.hsc_band(pool, len(frb), grid, np.random.default_rng(1))
    pm = pd.read_csv(MOR / "disk_subsample_tests.csv").set_index("label").loc["All hosts"]
    fig, ax = plt.subplots(figsize=(8, 5.6))
    ax.plot([0, 1], [0, 1], color=GREY, ls="--", lw=1.2, label="Random orientation")
    ax.fill_between(grid, lo, hi, color=BLUE, alpha=0.25, lw=0, label=f"Control: draws of {len(frb)} (68%)")
    ax.plot(grid, P.ecdf_at(pool, grid), color=BLUE, lw=2, label=f"Control disks (HSC, N = {len(pool):,})")
    x, y = P.step_xy(frb)
    ax.step(x, y, where="post", color=RED, lw=2.6, label=f"FRB hosts (N = {len(frb)})")
    ax.set(xlim=(0, 1), ylim=(0, 1), xlabel=r"$\cos i$", ylabel="Cumulative fraction")
    top = ax.secondary_xaxis("top", functions=(lambda c: np.degrees(np.arccos(np.clip(c, 0, 1))),
                                                lambda i: np.cos(np.radians(i))))
    top.set_xticks([90, 75, 60, 45, 30, 0])
    top.set_xlabel("Inclination (deg)", fontsize=15)
    ax.spines["top"].set_visible(True)
    ax.legend(loc="upper left", frameon=False)
    ax.text(0.97, 0.05, f"Anderson–Darling $p$ = {pm.ad_pm:.2f}\nMann–Whitney $p$ = {pm.mw2_pm:.2f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=14)
    finish(fig, "inclination_vs_control")


# ---------------------------------------------------------------- scatter plots
def finish(fig, stem: str) -> None:
    fig.tight_layout()
    fig.savefig(OUT / f"{stem}.png", dpi=300, facecolor="white")
    plt.close(fig)
    print(stem, flush=True)


def scatter(x, y, xerr, xlabel, ylabel, stem, logy=False, xlim=None) -> None:
    ok = np.isfinite(x) & np.isfinite(y) & ((y > 0) if logy else True)
    x, y = x[ok], y[ok]
    r, p = spearmanr(x, y)
    fig, ax = plt.subplots(figsize=(8, 5.6))
    if xerr is not None:
        ax.errorbar(x, y, xerr=xerr[ok], fmt="none", ecolor="0.75", elinewidth=1, zorder=1)
    ax.scatter(x, y, s=70, color=RED, edgecolor="white", linewidth=0.8, zorder=2)
    if logy:
        ax.set_yscale("log")
    ax.set(xlabel=xlabel, ylabel=ylabel)
    if xlim:
        ax.set_xlim(*xlim)
    ax.text(0.97, 0.95, f"Spearman $r_S$ = {r:+.2f}\n$p$ = {p:.2f}, N = {len(x)}", transform=ax.transAxes,
            ha="right", va="top", fontsize=14, bbox=dict(fc="white", ec="0.8", boxstyle="round,pad=0.4"))
    finish(fig, stem)


def correlations() -> None:
    d = R.load()
    sig_i = P.sigma_i_deg(d.q.to_numpy(float), d.sigma_q_B.to_numpy(float), np.random.default_rng(3))
    i = d.i_deg.to_numpy(float)
    inc = "Host inclination $i$ (deg)"
    scatter(i, d.dm_host_rest.to_numpy(float), sig_i, inc, r"DM$_{\rm host}$ (pc cm$^{-3}$, rest frame)",
            "dm_host_vs_inclination", xlim=(0, 90))
    scatter(i, d.abs_rm.to_numpy(float), sig_i, inc, r"|RM| (rad m$^{-2}$)", "rm_vs_inclination", logy=True, xlim=(0, 90))
    scatter(i, d.tau_scatt_ms.to_numpy(float), sig_i, inc, r"Scattering time $\tau$ (ms)", "scattering_vs_inclination",
            logy=True, xlim=(0, 90))
    scatter(d.logmstar.to_numpy(float), d.dm_host_rest.to_numpy(float), None, r"$\log_{10}(M_*/M_\odot)$",
            r"DM$_{\rm host}$ (pc cm$^{-3}$, rest frame)", "dm_host_vs_stellar_mass")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    d = disk_hosts()
    print(len(d), d.cls.value_counts().to_dict())
    gallery(d)
    inclination_vs_control()
    correlations()


if __name__ == "__main__":
    main()
