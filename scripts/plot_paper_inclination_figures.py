"""Paper figures for the FRB host inclination comparison.

Seed 20240924. Outputs PDF and PNG under plots/paper_figures/.
"""
from __future__ import annotations

import json
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from scipy import stats

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import REPO_ROOT, hubble_cosi_from_ba  # noqa: E402
from plot_frb_inflated_error_cdf import _hubble_inc_deg  # noqa: E402
from plot_frb_vs_hsc_confirmed50 import OUT_DIR, load_hsc_pool  # noqa: E402

SEED = 20240924
Q0 = 0.2
N_DRAWS = 10_000
N_ERR = 2_000
FLOOR_DEG = 5.0
P_2024 = 4.0 / 23.0
MAG_LIMIT = 22.0
OUT = REPO_ROOT / "plots" / "paper_figures"
NATURE = _SCRIPTS.parent / "Archive" / "original paper" / "frb_sample.txt"
HSC_EXP = REPO_ROOT / "catalog" / "HSC_kawinwanichakij_exp_analogue.csv"
HSC_PARENT = REPO_ROOT / "catalog" / "HSC_kawinwanichakij_sample_500k.csv"
RESULTS = REPO_ROOT / "pipeline_galfit_results.csv"
MORPH_REJECT = {"20190711A", "20221101B", "20230930A", "20221116A"}

FRB_C = "#D62728"
HSC_C = "#1F77B4"
GREY = "#7F7F7F"
ORANGE = "#FF7F0E"
GREEN = "#2CA02C"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIXGeneral", "Times New Roman", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 8,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.linewidth": 0.6,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def save(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{stem}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {stem}")


def cosi_to_i(cosi: np.ndarray) -> np.ndarray:
    return np.degrees(np.arccos(np.clip(cosi, 0.0, 1.0)))


def ecdf_at(sample: np.ndarray, grid: np.ndarray) -> np.ndarray:
    s = np.sort(np.asarray(sample, dtype=float))
    return np.searchsorted(s, grid, side="right") / len(s)


def step_xy(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    v = np.sort(np.asarray(values, dtype=float))
    n = len(v)
    x = np.concatenate([[0.0], v, [1.0]])
    y = np.concatenate([[0.0], np.arange(1, n + 1) / n, [1.0]])
    return x, y


def hsc_band(pool: np.ndarray, n_sample: int, grid: np.ndarray, rng: np.random.Generator):
    n = len(pool)
    cdf = np.empty((N_DRAWS, len(grid)), dtype=float)
    counts60 = np.empty(N_DRAWS, dtype=int)
    # quartile counts use the caller
    cut = np.cos(np.radians(60.0))
    for i in range(N_DRAWS):
        draw = pool[rng.choice(n, size=n_sample, replace=False)]
        cdf[i] = ecdf_at(draw, grid)
        counts60[i] = int(np.sum(draw < cut))
    lo = np.percentile(cdf, 16.0, axis=0)
    hi = np.percentile(cdf, 84.0, axis=0)
    return lo, hi, counts60


def top_inclination(ax) -> None:
    ax.tick_params(top=False)
    top = ax.twiny()
    top.set_xlim(ax.get_xlim())
    degs = [90, 75, 60, 45, 30, 0]
    top.set_xticks(np.cos(np.radians(degs)))
    top.set_xticklabels([str(d) for d in degs])
    top.set_xlabel(r"Inclination $i$ (deg)")
    top.tick_params(direction="in", labelsize=8, length=3)
    return top


def fwhm_arcsec(frb: str) -> float:
    path = REPO_ROOT / "pipeline_scripts" / "Output" / f"{frb}_all" / "psfex.xml"
    if not path.is_file():
        return np.nan
    root = ET.parse(path).getroot()
    fields = root.findall(".//{*}FIELD")
    tds = root.findall(".//{*}TR/{*}TD")
    for fld, td in zip(fields, tds):
        if fld.attrib.get("name") == "FWHM_WCS_Mean":
            return float(td.text)
    return np.nan


def re_arcsec(frb: str, wkdir: str) -> float:
    prod = REPO_ROOT / "pipeline_scripts" / "Output" / f"{frb}_all" / "pipeline_summary.json"
    scale = np.nan
    if prod.is_file():
        summary = json.loads(prod.read_text(encoding="utf-8"))
        scale = float(summary.get("plate_scale_arcsec_px", np.nan))
        host = summary.get("galfit_host") or {}
        if "Re-fits" not in wkdir.replace("\\", "/") and host.get("re_arcsec") is not None:
            return float(host["re_arcsec"])
    sky = REPO_ROOT / wkdir / "sky.json"
    if not sky.is_file():
        sky = REPO_ROOT / "pipeline_scripts" / "verification" / "outputs" / "per_host" / frb / "sky.json"
    if sky.is_file() and np.isfinite(scale):
        data = json.loads(sky.read_text(encoding="utf-8"))
        if data.get("re_galfit_px") is not None:
            return float(data["re_galfit_px"]) * scale
    if prod.is_file():
        host = json.loads(prod.read_text(encoding="utf-8")).get("galfit_host") or {}
        if host.get("re_arcsec") is not None:
            return float(host["re_arcsec"])
    return np.nan


def sigma_i_deg(q: np.ndarray, sig: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    out = np.empty(len(q))
    for i, (qi, si) in enumerate(zip(q, sig)):
        qq = np.clip(rng.normal(qi, si, size=4000), Q0 + 1e-6, 1.0)
        inc = _hubble_inc_deg(qq, q0=Q0)
        inc = np.clip(inc + rng.normal(0.0, FLOOR_DEG, size=len(qq)), 0.0, 90.0)
        out[i] = 0.5 * (np.percentile(inc, 84) - np.percentile(inc, 16))
    return out


def load_hosts(path: Path | None = None) -> pd.DataFrame:
    hosts = pd.read_csv(path or OUT_DIR / "confirmed50_q.csv", dtype={"frb": str})
    paper = []
    for line in NATURE.read_text(encoding="utf-8").splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        paper.append((parts[1].strip(), int(parts[2].replace('"', "").split(",")[0])))
    paper = pd.DataFrame(paper, columns=["frb", "i_2024"])
    hosts = hosts.merge(paper, on="frb", how="left")
    hosts["i_deg"] = cosi_to_i(hosts["cosi_point"].to_numpy(dtype=float))
    hosts["fwhm"] = [fwhm_arcsec(f) for f in hosts["frb"]]
    hosts["re"] = [re_arcsec(f, w) for f, w in zip(hosts["frb"], hosts["wkdir"])]
    hosts["re_over_fwhm"] = hosts["re"] / hosts["fwhm"]
    rng = np.random.default_rng(SEED + 3)
    hosts["sigma_i"] = sigma_i_deg(
        hosts["q"].to_numpy(dtype=float),
        hosts["sigma_q_B"].to_numpy(dtype=float),
        rng,
    )
    return hosts


def load_hsc_ba(n_min: float, n_max: float, ba_cut: bool) -> pd.DataFrame:
    cols = ["rmag", "ba", "n_sersic", "goodfits_flag"]
    df = pd.read_csv(HSC_PARENT, usecols=cols)
    r = pd.to_numeric(df["rmag"], errors="coerce")
    ba = pd.to_numeric(df["ba"], errors="coerce")
    n = pd.to_numeric(df["n_sersic"], errors="coerce")
    ok = (
        (df["goodfits_flag"] == 1)
        & np.isfinite(r) & (r <= MAG_LIMIT)
        & np.isfinite(ba) & (ba > 0) & (ba <= 1)
        & np.isfinite(n) & (n > n_min) & (n < n_max)
    )
    if ba_cut:
        ok &= ba > Q0
    return pd.DataFrame({"rmag": r[ok], "ba": ba[ok], "n": n[ok]})


def fig1(hosts: pd.DataFrame, hsc: np.ndarray, rng: np.random.Generator) -> None:
    grid = np.linspace(0.0, 1.0, 500)
    frb = hosts["cosi_point"].to_numpy(dtype=float)
    nature = hosts.loc[hosts["i_2024"].notna()]
    lo51, hi51, _ = hsc_band(hsc, 51, grid, rng)
    lo23, hi23, _ = hsc_band(hsc, 23, grid, np.random.default_rng(SEED + 1))
    hsc_cdf = ecdf_at(hsc, grid)

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.15), sharey=True)
    panels = [
        (axes[0], frb, lo51, hi51, "HSC draws of 51 (68%)", "FRB hosts (N = 51)", r"AD $p_m = 0.61$", "a"),
        (axes[1], nature["cosi_point"].to_numpy(dtype=float), lo23, hi23,
         "HSC draws of 23 (68%)", "This work", r"AD $p_m = 0.93$ (this work)", "b"),
    ]
    for ax, sample, lo, hi, band_label, frb_label, note, letter in panels:
        ax.plot([0, 1], [0, 1], "k--", lw=0.6, label="Uniform", zorder=1)
        ax.fill_between(grid, lo, hi, color=HSC_C, alpha=0.3, lw=0, label=band_label, zorder=2)
        ax.plot(grid, hsc_cdf, color=HSC_C, lw=1.0, label="HSC disks", zorder=3)
        x, y = step_xy(sample)
        ax.step(x, y, where="post", color=FRB_C, lw=1.5, label=frb_label, zorder=4)
        if letter == "b":
            x0, y0 = step_xy(np.cos(np.radians(nature["i_2024"].to_numpy(dtype=float))))
            ax.step(x0, y0, where="post", color=GREY, lw=1.2, label="2024 values", zorder=3)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel(r"$\cos i$")
        ax.legend(loc="upper left", frameon=False)
        ax.text(0.97, 0.06, note, transform=ax.transAxes, ha="right", va="bottom", fontsize=8)
        ax.text(0.03, 0.97, f"({letter})", transform=ax.transAxes, ha="left", va="top", fontsize=9)
        top_inclination(ax)
    axes[0].set_ylabel("Cumulative fraction")
    fig.tight_layout()
    save(fig, "fig1_cdf")


def fig2(hosts: pd.DataFrame) -> None:
    sub = hosts.loc[hosts["i_2024"].notna()].copy()
    x = sub["i_2024"].to_numpy(dtype=float)
    y = sub["i_deg"].to_numpy(dtype=float)
    yerr = sub["sigma_i"].to_numpy(dtype=float)
    c = sub["re_over_fwhm"].to_numpy(dtype=float)
    fig, axes = plt.subplots(
        2, 1, figsize=(3.5, 6.0), sharex=True,
        gridspec_kw={"height_ratios": [1.35, 1.0], "hspace": 0.08},
    )
    ax = axes[0]
    ax.errorbar(x, y, yerr=yerr, fmt="none", ecolor="0.45", elinewidth=0.6, capsize=0, zorder=2)
    sc = ax.scatter(x, y, c=c, cmap="viridis", s=16, zorder=3, edgecolors="none")
    ax.plot([0, 90], [0, 90], "k--", lw=0.6, zorder=1)
    ax.axvline(60, color="0.4", ls=":", lw=0.6)
    ax.axhline(60, color="0.4", ls=":", lw=0.6)
    cross = (x < 60) & (y > 60)
    offsets = {
        "20210410D": (-4, 8),
        "20220914A": (-34, 2),
        "20190102C": (5, 4),
        "20190714A": (5, -8),
    }
    for _, row in sub.loc[cross].iterrows():
        ax.annotate(
            row["frb"], (row["i_2024"], row["i_deg"]),
            textcoords="offset points", xytext=offsets.get(row["frb"], (3, 3)),
            fontsize=5, color="0.15",
        )
    ax.set_xlim(0, 90)
    ax.set_ylim(0, 90)
    ax.set_ylabel(r"$i$ (this work) [deg]")
    ax.set_aspect("equal", adjustable="box")
    ax.tick_params(top=True, right=True)
    cb = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label(r"$R_e$ / FWHM")
    cb.ax.tick_params(labelsize=7)

    ax2 = axes[1]
    ax2.scatter(x, y - x, c=c, cmap="viridis", s=16, edgecolors="none", norm=sc.norm)
    ax2.axhline(0, color="k", lw=0.6)
    ax2.set_xlabel(r"$i$ (Bhardwaj et al. 2024) [deg]")
    ax2.set_ylabel(r"$\Delta i$ [deg]")
    ax2.tick_params(top=True, right=True)
    fig.tight_layout()
    save(fig, "fig2_old_vs_new")
    print("crossed 60:", ", ".join(sub.loc[cross, "frb"].astype(str)))


def fig3(hsc: np.ndarray, rng: np.random.Generator) -> None:
    n = len(hsc)
    cut = np.cos(np.radians(60.0))
    counts = np.empty(N_DRAWS, dtype=int)
    for i in range(N_DRAWS):
        draw = hsc[rng.choice(n, size=51, replace=False)]
        counts[i] = int(np.sum(draw < cut))
    bins = np.arange(-0.5, 35.5, 1.0)
    xs = np.arange(0, 36)
    fig, ax = plt.subplots(figsize=(3.5, 3.0))
    ax.hist(
        counts, bins=bins, density=True, histtype="stepfilled",
        color=HSC_C, alpha=0.5, label="HSC draws",
    )
    p_hsc = float(np.mean(hsc < cut))
    p_scaled = (P_2024 / 0.5) * p_hsc
    ax.plot(xs, stats.binom.pmf(xs, 51, P_2024), color=ORANGE, lw=1.0, marker="o", ms=2.0, label="2024 rate (0.174)")
    ax.plot(xs, stats.binom.pmf(xs, 51, p_scaled), color=GREEN, lw=1.0, marker="o", ms=2.0, label="2024 rate scaled to HSC (0.138)")
    ax.axvline(15, color=FRB_C, lw=1.0, label="Observed (15)")
    ax.set_xlim(0, 35)
    ax.set_xlabel(r"$N(i > 60^\circ)$ out of 51")
    ax.set_ylabel("Probability")
    ax.legend(frameon=False, fontsize=6)
    ax.tick_params(top=True, right=True)
    fig.tight_layout()
    save(fig, "fig3_edgeon_count")


def fig4(hosts: pd.DataFrame, hsc: np.ndarray, rng: np.random.Generator) -> None:
    edges = np.quantile(hsc, np.linspace(0, 1, 5))
    edges[0], edges[-1] = -np.inf, np.inf
    frb = hosts["cosi_point"].to_numpy(dtype=float)
    counts = np.histogram(frb, bins=edges)[0]
    draw_counts = np.empty((N_DRAWS, 4), dtype=int)
    n = len(hsc)
    for i in range(N_DRAWS):
        draw = hsc[rng.choice(n, size=51, replace=False)]
        draw_counts[i] = np.histogram(draw, bins=edges)[0]
    lo = np.percentile(draw_counts, 16, axis=0)
    hi = np.percentile(draw_counts, 84, axis=0)
    finite = edges.copy()
    finite[0], finite[-1] = hsc.min(), hsc.max()
    i_hi = cosi_to_i(finite[:-1])
    i_lo = cosi_to_i(finite[1:])
    labels = [f"Bin {k}" for k in range(1, 5)]
    top_labels = [f"{a:.0f}–{b:.0f}" for a, b in zip(i_hi, i_lo)]

    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.15))
    x = np.arange(4)
    ax = axes[0]
    ax.bar(x, counts, color=FRB_C, width=0.7, label="FRB hosts")
    yerr = np.vstack([np.maximum(12.75 - lo, 0.0), np.maximum(hi - 12.75, 0.0)])
    ax.errorbar(x, np.full(4, 12.75), yerr=yerr, fmt="none", ecolor="k", elinewidth=0.7, capsize=2, label="HSC draws")
    ax.axhline(12.75, color="k", ls="--", lw=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Count")
    ax.set_xlabel("HSC quartile, edge-on to face-on")
    ax.set_ylim(0, 24)
    ax.legend(frameon=False, fontsize=6, loc="upper center", ncol=2)
    ax.tick_params(top=False)
    top = ax.twiny()
    top.set_xlim(ax.get_xlim())
    top.set_xticks(x)
    top.set_xticklabels(top_labels, fontsize=6)
    top.set_xlabel(r"$i$ range (deg)")
    top.tick_params(direction="in", length=3)
    ax.text(0.03, 0.95, "(a)", transform=ax.transAxes, va="top")

    ax = axes[1]
    counts_2024 = np.array([1, 3, 8, 11])
    ranges = ["76–90", "60–76", "42–60", "0–42"]
    ax.bar(x, counts_2024, color=GREY, width=0.7, label="2024 values")
    ax.axhline(5.75, color="k", ls="--", lw=0.6, label="Expected 5.75")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel("SDSS quartile, edge-on to face-on")
    ax.set_ylim(0, 16)
    ax.legend(frameon=False, fontsize=6)
    ax.tick_params(top=False)
    top = ax.twiny()
    top.set_xlim(ax.get_xlim())
    top.set_xticks(x)
    top.set_xticklabels(ranges, fontsize=6)
    top.set_xlabel(r"$i$ range (deg)")
    top.tick_params(direction="in", length=3)
    ax.text(0.03, 0.95, "(b)", transform=ax.transAxes, va="top")
    fig.tight_layout()
    save(fig, "fig4_quartiles")
    print("quartile counts", counts.tolist())


def fig5(hosts: pd.DataFrame) -> None:
    prod = pd.read_csv(RESULTS, dtype={"frb": str})
    mag = pd.to_numeric(prod["mag"], errors="coerce")
    ba = pd.to_numeric(prod["b_a"], errors="coerce")
    keep = mag.notna() & (mag <= MAG_LIMIT) & ba.notna() & ~prod["frb"].isin(MORPH_REJECT)
    tab = prod.loc[keep, ["frb"]].copy()
    tab["ba"] = ba[keep].to_numpy()
    win = hosts.set_index("frb")["q"]
    tab["ba"] = tab["frb"].map(win).fillna(tab["ba"])
    hsc = load_hsc_ba(0.4, 1.5, ba_cut=False)
    bins = np.arange(0.0, 1.0 + 1e-9, 0.05)
    fig, ax = plt.subplots(figsize=(3.5, 3.0))
    ax.hist(tab["ba"], bins=bins, color=FRB_C, histtype="step", lw=1.2, label=f"FRB hosts (N = {len(tab)})")
    ax.hist(
        hsc["ba"], bins=bins, weights=np.full(len(hsc), len(tab) / len(hsc)),
        histtype="step", color=HSC_C, lw=1.2, label="HSC disks",
    )
    ax.axvline(0.2, color="k", ls="--", lw=0.6)
    n_low = int(np.sum(tab["ba"] <= 0.2))
    hsc_low = float(np.mean(hsc["ba"] <= 0.2))
    ax.text(
        0.98, 0.28,
        f"FRB: {n_low}/{len(tab)} = {100 * n_low / len(tab):.1f}%\n"
        f"HSC: {100 * hsc_low:.1f}%",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=7,
    )
    ax.set_xlim(0, 1)
    ax.set_xlabel(r"$b/a$")
    ax.set_ylabel("Count")
    ax.legend(frameon=False, fontsize=7, loc="upper left")
    ax.tick_params(top=True, right=True)
    fig.tight_layout()
    save(fig, "fig5_ba_hist")
    print(f"fig5 N={len(tab)} below 0.2: {n_low}")


def fig_pvalues(hosts: pd.DataFrame, hsc: np.ndarray) -> None:
    q = hosts["q"].to_numpy(dtype=float)
    sig = hosts["sigma_q_B"].to_numpy(dtype=float)
    point = hosts["cosi_point"].to_numpy(dtype=float)
    with np.errstate(invalid="ignore"):
        ad0 = float(stats.anderson_ksamp([hsc, point]).pvalue)
    ks0 = float(stats.ks_2samp(point, hsc, alternative="less").pvalue)
    mw0 = float(stats.mannwhitneyu(point, hsc, alternative="greater").pvalue)
    rng = np.random.default_rng(SEED)
    ad_p = np.empty(N_ERR)
    ks_p = np.empty(N_ERR)
    mw_p = np.empty(N_ERR)
    for i in range(N_ERR):
        qq = np.clip(rng.normal(q, sig), Q0 + 1e-6, 1.0)
        inc = np.clip(_hubble_inc_deg(qq, q0=Q0) + rng.normal(0.0, FLOOR_DEG, size=len(q)), 0.0, 90.0)
        cosi = np.cos(np.radians(inc))
        ad_p[i] = float(stats.anderson_ksamp([hsc, cosi]).pvalue)
        ks_p[i] = float(stats.ks_2samp(cosi, hsc, alternative="less").pvalue)
        mw_p[i] = float(stats.mannwhitneyu(cosi, hsc, alternative="greater").pvalue)
        if i % 500 == 0:
            print("pvalue draw", i)
    fig, axes = plt.subplots(1, 3, figsize=(7.1, 2.6))
    specs = [
        (ad_p, ad0, "Anderson–Darling"),
        (ks_p, ks0, "KS, face-on"),
        (mw_p, mw0, "Mann–Whitney, face-on"),
    ]
    for ax, (vals, mark, title) in zip(axes, specs):
        ax.hist(vals, bins=20, color=HSC_C, alpha=0.7)
        ax.axvline(mark, color=FRB_C, lw=1.0)
        ax.set_xlabel(r"$p$")
        ax.set_title(title, fontsize=8)
        ax.tick_params(top=True, right=True)
    axes[0].set_ylabel("Draws")
    fig.tight_layout()
    save(fig, "figA_pvalues")


def magnitude_q(hosts: pd.DataFrame) -> pd.Series:
    """Confirmed hosts with m_r <= 22. Winning-leg q, else the production q.

    Hosts rejected as bad fits stay out. A confirmed host with q <= 0.2
    is kept so a lower q0 can admit it.
    """
    prod = pd.read_csv(RESULTS, dtype={"frb": str})
    hc = pd.read_csv(
        REPO_ROOT / "pipeline_scripts" / "verification" / "host_confirmation.csv",
        dtype={"frb": str},
    )
    confirmed = set(hc.loc[hc["confirmed"].astype(str).str.lower().eq("true"), "frb"])
    mag = pd.to_numeric(prod["mag"], errors="coerce")
    ba = pd.to_numeric(prod["b_a"], errors="coerce")
    keep = mag.notna() & (mag <= MAG_LIMIT) & ba.notna() & prod["frb"].isin(confirmed)
    tab = prod.loc[keep, ["frb"]].copy()
    tab["q"] = ba[keep].to_numpy()
    win = hosts.set_index("frb")["q"]
    tab["q"] = tab["frb"].map(win).fillna(tab["q"])
    return tab.set_index("frb")["q"]


def fig_q0(hosts: pd.DataFrame) -> None:
    parent = load_hsc_ba(0.4, 1.5, ba_cut=False)
    q0s = [0.10, 0.15, 0.20, 0.25]
    fig, axes = plt.subplots(2, 2, figsize=(7.1, 6.2), sharex=True, sharey=True)
    qseries = magnitude_q(hosts)
    q = qseries.to_numpy(dtype=float)
    print("q0 pool", len(q), "min", float(np.min(q)))
    for ax, q0 in zip(axes.ravel(), q0s):
        added = qseries[(qseries > q0) & (qseries <= 0.2)]
        if len(added):
            print(f"q0={q0:.2f} adds", ", ".join(f"{i} ({v:.2f})" for i, v in added.items()))
        frb_q = q[q > q0]
        hsc_q = parent.loc[parent["ba"] > q0, "ba"].to_numpy(dtype=float)
        frb_c = np.array([hubble_cosi_from_ba(float(v), q0=q0) for v in frb_q])
        hsc_c = np.array([hubble_cosi_from_ba(float(v), q0=q0) for v in hsc_q])
        grid = np.linspace(0, 1, 400)
        ax.plot(grid, ecdf_at(hsc_c, grid), color=HSC_C, lw=1.0, label="HSC")
        x, y = step_xy(frb_c)
        ax.step(x, y, where="post", color=FRB_C, lw=1.2, label=f"FRB (N = {len(frb_c)})")
        ax.plot([0, 1], [0, 1], "k--", lw=0.5)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_title(rf"$q_0 = {q0:.2f}$", fontsize=8)
        ax.tick_params(top=True, right=True)
        ax.legend(frameon=False, fontsize=6)
    for ax in axes[1]:
        ax.set_xlabel(r"$\cos i$")
    for ax in axes[:, 0]:
        ax.set_ylabel("Cumulative fraction")
    fig.tight_layout()
    save(fig, "figA_q0")


def fig_nrange(hosts: pd.DataFrame) -> None:
    windows = [(0.4, 1.5), (0.5, 1.2), (0.8, 1.2), (0.4, 2.0)]
    frb = np.sort(hosts["cosi_point"].to_numpy(dtype=float))
    fig, axes = plt.subplots(2, 2, figsize=(7.1, 6.2), sharex=True, sharey=True)
    for ax, (lo, hi) in zip(axes.ravel(), windows):
        hsc = load_hsc_ba(lo, hi, ba_cut=True)
        cosi = np.array([hubble_cosi_from_ba(float(v), q0=Q0) for v in hsc["ba"]])
        grid = np.linspace(0, 1, 400)
        ax.plot(grid, ecdf_at(cosi, grid), color=HSC_C, lw=1.0, label=f"HSC (N = {len(cosi):,})")
        x, y = step_xy(frb)
        ax.step(x, y, where="post", color=FRB_C, lw=1.2, label="FRB (N = 51)")
        ax.plot([0, 1], [0, 1], "k--", lw=0.5)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_title(rf"${lo:.1f} < n < {hi:.1f}$", fontsize=8)
        ax.tick_params(top=True, right=True)
        ax.legend(frameon=False, fontsize=6)
    for ax in axes[1]:
        ax.set_xlabel(r"$\cos i$")
    for ax in axes[:, 0]:
        ax.set_ylabel("Cumulative fraction")
    fig.tight_layout()
    save(fig, "figA_nrange")


def main() -> None:
    print(f"seed {SEED}")
    OUT.mkdir(parents=True, exist_ok=True)
    hosts = load_hosts()
    if len(hosts) != 51:
        raise RuntimeError(f"expected 51 hosts, got {len(hosts)}")
    hosts.to_csv(OUT / "fig_host_inputs.csv", index=False)
    print(
        "Re/FWHM missing",
        int((~np.isfinite(hosts["re_over_fwhm"])).sum()),
        "sigma_i median",
        float(np.median(hosts["sigma_i"])),
    )
    hsc, n_hsc = load_hsc_pool(mag_limit=MAG_LIMIT, q0=Q0)
    print(f"HSC {n_hsc}")
    rng = np.random.default_rng(SEED)
    fig1(hosts, hsc, rng)
    fig2(hosts)
    fig3(hsc, np.random.default_rng(SEED + 2))
    fig4(hosts, hsc, np.random.default_rng(SEED + 4))
    fig5(hosts)
    fig_pvalues(hosts, hsc)
    fig_q0(hosts)
    fig_nrange(hosts)


if __name__ == "__main__":
    main()
