"""Slide-ready figures.

disk_gallery_{1,2,3}.png   colour images of the 51 disk hosts, 18 per 16:9 page, with survey and disk evidence
inclination_vs_control.png FRB host cos i CDF against the HSC disk control
dm_host_vs_inclination.png, rm_vs_inclination.png, scattering_vs_inclination.png, dm_host_vs_stellar_mass.png

Error bars are 16th-84th percentiles of Monte Carlo draws: inclination from sigma_q (formal + sky) through the
Hubble relation plus a 5 deg floor; DM_host from DM_cosmic drawn from the Macquart et al. (2020) distribution and
DM_MW (NE2001) with a 20% error; M* from the quoted literature errors (split normal), |RM| and tau likewise in log10. Spearman r_S
is quoted with the 16th-84th percentile range over the same draws.

Output: Publication/materials/slides/
"""
from __future__ import annotations

import re
import sys
from dataclasses import dataclass
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
import build_host_property_table as B  # noqa: E402
import dm_host_physics as DP  # noqa: E402
import host_gallery as G  # noqa: E402
import plot_paper_inclination_figures as P  # noqa: E402
import robust_correlations as R  # noqa: E402

OUT = ROOT / "Publication" / "materials" / "slides"
MOR = ROOT / "Publication" / "materials" / "appendix" / "morphology"
PER_PAGE, NCOL = 18, 6
N_MC = 2000
MSTAR_ASSUMED_DEX = 0.2
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
    rng = np.random.default_rng(5)
    q, sq = hosts["q"].to_numpy(float), np.nan_to_num(hosts["sigma_q_B"].to_numpy(float))
    qq = np.clip(rng.normal(q, sq, (N_MC, len(q))), P.Q0 + 1e-6, 1.0)
    inc = np.clip(P._hubble_inc_deg(qq, q0=P.Q0) + rng.normal(0, P.FLOOR_DEG, qq.shape), 0, 90)
    cdf = np.array([P.ecdf_at(c, grid) for c in np.cos(np.radians(inc))])
    flo, fhi = np.percentile(cdf, [16, 84], axis=0)
    ax.fill_between(grid, flo, fhi, color=RED, alpha=0.2, lw=0, label="FRB hosts: inclination errors (68%)")
    x, y = P.step_xy(frb)
    ax.step(x, y, where="post", color=RED, lw=2.6, label=f"FRB hosts (N = {len(frb)})")
    ax.set(xlim=(0, 1), ylim=(0, 1), xlabel=r"$\cos i$", ylabel="Cumulative fraction")
    top = ax.secondary_xaxis("top", functions=(lambda c: np.degrees(np.arccos(np.clip(c, 0, 1))),
                                                lambda i: np.cos(np.radians(i))))
    top.set_xticks([90, 75, 60, 45, 30, 0])
    top.set_xlabel("Inclination (deg)", fontsize=15)
    ax.spines["top"].set_visible(True)
    ax.legend(loc="upper left", frameon=False)
    ax.text(0.97, 0.05, f"Median over inclination-error draws:\nAnderson–Darling $p$ = {pm.ad_pm:.2f}\n"
                        f"Mann–Whitney $p$ = {pm.mw2_pm:.2f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=14)
    finish(fig, "inclination_vs_control")


# ---------------------------------------------------------------- scatter plots
def finish(fig, stem: str) -> None:
    fig.tight_layout()
    fig.savefig(OUT / f"{stem}.png", dpi=300, facecolor="white")
    plt.close(fig)
    print(stem, flush=True)


@dataclass
class Quantity:
    value: np.ndarray
    draws: np.ndarray
    label: str
    assumed: np.ndarray | None = None

    def bars(self, ok):
        lo, hi = np.nanpercentile(self.draws[:, ok], [16, 84], axis=0)
        v = self.value[ok]
        return np.vstack([np.clip(v - lo, 0, None), np.clip(hi - v, 0, None)])


def split_normal(v, ep, em, rng, n=N_MC):
    u = rng.standard_normal((n, len(v)))
    return v + np.where(u > 0, u * ep, u * em)


def split_lognormal(v, ep, em, rng):
    """Positive quantities: split normal in log10, with v + ep and v - em at +-1 sigma (lower side mirrors the
    upper when em >= v)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        sp = np.log10((v + ep) / v)
        sm = np.where(em < v, np.log10(v / (v - em)), sp)
    return 10 ** split_normal(np.log10(v), sp, sm, rng)


def literature_errors(d: pd.DataFrame, q: str) -> tuple[np.ndarray, np.ndarray]:
    """err_plus, err_minus of the adopted literature value; if that row has none, the error another source
    quotes for the same burst."""
    L = B.load_long()
    ep, em = np.full(len(d), np.nan), np.full(len(d), np.nan)
    for k, f in enumerate(d.frb):
        rows = L[(L.frb == f) & (L.quantity == q)]
        got = B.pick(L, f, q)
        cand = ([got[1:3]] if got else []) + list(zip(rows.err_plus, rows.err_minus))
        for p, m in cand:
            p, m = abs(B.to_num(p)), abs(B.to_num(m))
            if np.isfinite(p):
                ep[k], em[k] = p, m if np.isfinite(m) else p
                break
    return ep, em


def quantities(d: pd.DataFrame) -> dict[str, Quantity]:
    rng = np.random.default_rng(7)
    q, sq = d.q.to_numpy(float), np.nan_to_num(d.sigma_q_B.to_numpy(float))
    qq = np.clip(rng.normal(q, sq, (N_MC, len(q))), P.Q0 + 1e-6, 1.0)
    inc = np.clip(P._hubble_inc_deg(qq, q0=P.Q0) + rng.normal(0, P.FLOOR_DEG, qq.shape), 0, 90)

    DP.N_MC = N_MC
    dm = DP.dm_host_draws(d, rng)

    ms = d.logmstar.to_numpy(float)
    p57 = pd.read_csv(ROOT / "Publication/materials/host_properties/host_properties_57.csv",
                      dtype={"frb": str}).set_index("frb")
    ep = pd.to_numeric(d.frb.map(p57.logmstar_ep), errors="coerce").abs().to_numpy()
    em = pd.to_numeric(d.frb.map(p57.logmstar_em), errors="coerce").abs().to_numpy()
    lp, lm = literature_errors(d, "logmstar")
    ep, em = np.where(np.isfinite(ep), ep, lp), np.where(np.isfinite(em), em, lm)
    assumed = np.isfinite(ms) & ~np.isfinite(ep)
    ep, em = np.where(assumed, MSTAR_ASSUMED_DEX, ep), np.where(assumed, MSTAR_ASSUMED_DEX, em)

    rm = d.abs_rm.to_numpy(float)
    rp, rmn = literature_errors(d, "rm_rad_m2")
    tau = d.tau_scatt_ms.to_numpy(float)
    tp, tm = literature_errors(d, "tau_scatt_ms")
    for name, v, e in (("|RM|", rm, rp), ("tau", tau, tp)):
        miss = np.isfinite(v) & ~np.isfinite(e)
        if miss.any():
            raise ValueError(f"{name} without an error: {list(d.frb[miss])}")

    return dict(
        inc=Quantity(d.i_deg.to_numpy(float), inc, "Host inclination $i$ (deg)"),
        dm=Quantity(d.dm_host_rest.to_numpy(float), dm, r"DM$_{\rm host}$ (pc cm$^{-3}$, rest frame)"),
        ms=Quantity(ms, split_normal(ms, ep, em, rng), r"$\log_{10}(M_*/M_\odot)$", assumed),
        rm=Quantity(rm, split_lognormal(rm, rp, rmn, rng), r"|RM| (rad m$^{-2}$)"),
        tau=Quantity(tau, split_lognormal(tau, tp, tm, rng), r"Scattering time $\tau$ (ms)"),
    )


def scatter(X: Quantity, Y: Quantity, stem: str, logy=False, xlim=None) -> None:
    ok = np.isfinite(X.value) & np.isfinite(Y.value)
    x, y = X.value[ok], Y.value[ok]
    r, p = spearmanr(x, y)
    r_mc = np.array([spearmanr(a, b)[0] for a, b in zip(X.draws[:, ok], Y.draws[:, ok])])
    r_lo, r_hi = np.percentile(r_mc, [16, 84])
    flagged = X.assumed[ok] if X.assumed is not None else np.zeros(len(x), bool)

    fig, ax = plt.subplots(figsize=(8, 5.6))
    ax.errorbar(x, y, xerr=X.bars(ok), yerr=Y.bars(ok), fmt="none", ecolor="0.7", elinewidth=1, zorder=1)
    ax.scatter(x[~flagged], y[~flagged], s=70, color=RED, edgecolor="white", linewidth=0.8, zorder=2)
    if flagged.any():
        ax.scatter(x[flagged], y[flagged], s=70, facecolor="white", edgecolor=RED, linewidth=1.6, zorder=2,
                   label=fr"No quoted error ($\pm${MSTAR_ASSUMED_DEX} dex assumed)")
        ax.legend(loc="lower left", frameon=False, fontsize=12)
    if logy:
        ax.set_yscale("log")
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi * (hi / lo) ** 0.22)
    else:
        lo, hi = ax.get_ylim()
        ax.set_ylim(lo, hi + 0.22 * (hi - lo))
    ax.set(xlabel=X.label, ylabel=Y.label)
    if xlim:
        ax.set_xlim(*xlim)
    ax.text(0.97, 0.95, f"Spearman $r_S$ = {r:+.2f} ({r_lo:+.2f}, {r_hi:+.2f})\n$p$ = {p:.2f}, N = {len(x)}",
            transform=ax.transAxes, ha="right", va="top", fontsize=14,
            bbox=dict(fc="white", ec="0.8", boxstyle="round,pad=0.4"))
    finish(fig, stem)
    print(f"  r_S = {r:+.3f} [{r_lo:+.3f}, {r_hi:+.3f}] p = {p:.3f} N = {len(x)}", flush=True)


def correlations() -> None:
    Q = quantities(R.load())
    scatter(Q["inc"], Q["dm"], "dm_host_vs_inclination", xlim=(0, 90))
    scatter(Q["inc"], Q["rm"], "rm_vs_inclination", logy=True, xlim=(0, 90))
    scatter(Q["inc"], Q["tau"], "scattering_vs_inclination", logy=True, xlim=(0, 90))
    scatter(Q["ms"], Q["dm"], "dm_host_vs_stellar_mass")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    d = disk_hosts()
    print(len(d), d.cls.value_counts().to_dict())
    gallery(d)
    inclination_vs_control()
    correlations()


if __name__ == "__main__":
    main()
