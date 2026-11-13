"""What sets the host DM? Star formation, stellar mass, the host halo, redshift selection, and the burst's
propagation (scattering, RM, polarization).

Every relation gets the same outlier and robustness checks:
  Spearman and Kendall, bootstrap 95% interval, delete-1/2/3 jackknife over all subsets (fraction of 3-host
  subsets with p < 0.05), the most influential host, Theil-Sen slope, Kendall's partial tau given z
  (Akritas & Siebert 1996; p from permutations within redshift strata), and, for DM_host, 2,000 Monte Carlo
  draws of DM_cosmic from the Macquart et al. (2020) distribution and of DM_MW (NE2001, 20%).
Benjamini-Hochberg q-values are computed over all hypothesis tests; the sanity checks are kept separate.

Output: Publication/materials/appendix/dm_host_physics/{tests.csv, halo_budget.csv, DM_HOST_PHYSICS.md, figure}
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.cosmology import Planck18 as cosmo
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_host_property_table as BT  # noqa: E402
import foreground_halo_dm as F  # noqa: E402
import redshift_bins as ZB  # noqa: E402
import robust_correlations as RC  # noqa: E402
import robust_extras as E  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Publication" / "materials" / "appendix" / "dm_host_physics"
PROPS_57 = ROOT / "Publication" / "materials" / "host_properties" / "host_properties_57.csv"
SEED = ROOT / "Publication" / "materials" / "host_properties" / "seed_56.csv"
N_MC = 2_000
N_PERM = 2_000
MACQUART_F = 0.2
DM_MW_FRAC_ERR = 0.2
SHMR_SCATTER_DEX = 0.15
DM = "dm_host_rest"

LABELS = {
    "z": "z", DM: "DM_host", "dm_host_obs": "DM_host (observer frame)", "logmstar": "log M*",
    "log_sfr": "log SFR", "log_ssfr": "log sSFR", "log_sigma_sfr": "log Sigma_SFR", "delta_ms": "Delta MS",
    "log_sigma_star_best": "log Sigma*", "re_kpc_best": "R_e", "offset_kpc": "offset (kpc)",
    "offset_norm": "offset / R_e", "av_old": "A_V (old stars)", "av": "A_V (young stars)",
    "dm_halo_pred": "predicted host-halo DM", "log_tau": "log tau", "log_abs_rm": "log |RM|",
    "lin_pol": "linear polarization", "metallicity_logz": "stellar metallicity", "tm_gyr": "mass-weighted age",
}

SAMPLES = {
    "All": lambda h: h,
    "z < 0.2": lambda h: h[h.z < 0.2],
    "z >= 0.2": lambda h: h[h.z >= 0.2],
    "Prospector SFRs": lambda h: h[h.sfr_method.str.startswith("prospector", na=False)],
    "Prospector, z >= 0.2": lambda h: h[h.sfr_method.str.startswith("prospector", na=False) & (h.z >= 0.2)],
    "literature R_e": lambda h: h[h.re_kpc.notna()],
}

# (theme, x, y, sample, what it tests)
TESTS = [
    ("Redshift selection", "z", "logmstar", "All", "flux limit: more massive hosts at higher z"),
    ("Redshift selection", "z", "log_sfr", "All", "flux limit: more star-forming hosts at higher z"),
    ("Redshift selection", "z", "log_ssfr", "All", "cosmic sSFR rises with z"),
    ("Redshift selection", "z", "re_kpc_best", "All", "larger hosts at higher z"),
    ("Redshift selection", "z", "offset_kpc", "All", "localization in kpc degrades with z"),
    ("Redshift selection", "z", "dm_host_obs", "All", "observer-frame DM_host diluted by 1/(1+z)"),
    ("Redshift selection", "z", DM, "All", "Verdi+26: no evolution"),
    ("Stars and star formation", "logmstar", DM, "All", "Leung+25: anticorrelation"),
    ("Stars and star formation", "log_sfr", DM, "All", "Li+25: positive"),
    ("Stars and star formation", "log_sfr", DM, "z < 0.2", ""),
    ("Stars and star formation", "log_sfr", DM, "z >= 0.2", "earlier lead: r_S = +0.55"),
    ("Stars and star formation", "log_sfr", DM, "Prospector SFRs", "one SFR method"),
    ("Stars and star formation", "log_sfr", DM, "Prospector, z >= 0.2", "one SFR method"),
    ("Stars and star formation", "log_ssfr", DM, "All", "Li+25: tight positive"),
    ("Stars and star formation", "log_ssfr", DM, "z >= 0.2", ""),
    ("Stars and star formation", "log_sigma_sfr", DM, "All", "H-alpha emission measure: DM ~ Sigma_SFR^0.5"),
    ("Stars and star formation", "log_sigma_sfr", DM, "z >= 0.2", ""),
    ("Stars and star formation", "log_sigma_sfr", DM, "Prospector SFRs", "one SFR method"),
    ("Stars and star formation", "log_sigma_sfr", DM, "literature R_e", "no R_e from our fits"),
    ("Stars and star formation", "delta_ms", DM, "All", "gas-rich hosts above the main sequence"),
    ("Stars and star formation", "log_sigma_star_best", DM, "All", "denser discs"),
    ("Stars and star formation", "re_kpc_best", DM, "All", "longer path through larger discs"),
    ("Stars and star formation", "offset_norm", DM, "All", "Mannings+21, Marnoch+26: less gas at large offset"),
    ("Stars and star formation", "av_old", DM, "All", "dust column traces gas column"),
    ("Stars and star formation", "av", DM, "All", "dust column traces gas column"),
    ("Host halo", "dm_halo_pred", DM, "All", "halo gas sets the floor of DM_host"),
    ("Propagation", DM, "log_tau", "All", "Cordes+22: scattering grows with DM"),
    ("Propagation", DM, "log_abs_rm", "All", "Verdi+26: RM grows with DM"),
    ("Propagation", "log_sigma_sfr", "log_tau", "All", "turbulent star-forming gas scatters"),
    ("Propagation", "log_sigma_star_best", "log_tau", "All", "Glowacki+25: positive"),
    ("Propagation", "offset_norm", "log_tau", "All", "less scattering at large offset"),
    ("Propagation", "log_sigma_sfr", "log_abs_rm", "All", "stronger fields in star-forming gas"),
    ("Propagation", "offset_norm", "log_abs_rm", "All", "less RM at large offset"),
    ("Propagation", "log_abs_rm", "lin_pol", "All", "depolarization by Faraday rotation"),
    ("Sanity checks", "logmstar", "log_sfr", "All", "star-forming main sequence"),
    ("Sanity checks", "logmstar", "metallicity_logz", "All", "mass-metallicity relation"),
    ("Sanity checks", "logmstar", "re_kpc_best", "All", "mass-size relation"),
    ("Sanity checks", "logmstar", "tm_gyr", "All", "downsizing: massive hosts older"),
    ("Sanity checks", "logmstar", "av", "All", "dustier massive hosts"),
]


# ---------------------------------------------------------------- derived quantities
def speagle14_log_sfr(logmstar, z):
    t = cosmo.age(z).value
    return (0.84 - 0.026 * t) * logmstar - (6.51 - 0.11 * t)


def host_halo_dm(logms, logms_err, z, offset_kpc, rng, n=500):
    """Rest-frame DM from the FRB to the edge of its own halo: half the mNFW chord at the projected offset."""
    logmh = F.logmh_from_logms(logms + rng.normal(0, logms_err, n), np.full(n, z))
    logmh = logmh + rng.normal(0, SHMR_SCATTER_DEX, n)
    theta = offset_kpc / F.kpc_per_arcsec(z)
    dm = F.dm_mnfw(logmh, np.full(n, z), np.full(n, theta)) * (1 + z) / 2
    return np.percentile(dm, [16, 50, 84])


def load() -> pd.DataFrame:
    h = ZB.hosts()
    extra = pd.read_csv(PROPS_57, dtype={"frb": str})[["frb", "logmstar_ep", "logmstar_em", "sfr_method"]]
    seed = pd.read_csv(SEED, dtype={"frb": str})[["frb", "re_arcsec_galfit"]]
    h = h.merge(extra, on="frb", how="left").merge(seed, on="frb", how="left")

    h["log_sfr"] = np.log10(h.sfr.where(h.sfr > 0))
    z = h.z.fillna(0.1).to_numpy()
    h["re_kpc_ours"] = np.where(h.z.notna(), h.re_arcsec_galfit * F.kpc_per_arcsec(z), np.nan)
    h["re_kpc_best"] = h.re_kpc.fillna(h.re_kpc_ours)
    area = 2 * np.pi * h.re_kpc_best ** 2
    h["log_sigma_sfr"] = h.log_sfr - np.log10(area)
    h["log_sigma_star_best"] = h.logmstar - np.log10(area)
    h["offset_norm"] = h.offset_kpc / h.re_kpc_best
    h["delta_ms"] = np.where(h.z.notna(), h.log_sfr - speagle14_log_sfr(h.logmstar, z), np.nan)

    rng = np.random.default_rng(11)
    halo = np.full((len(h), 3), np.nan)
    for k, r in enumerate(h.itertuples()):
        if np.isfinite(r.logmstar) and np.isfinite(r.offset_kpc) and np.isfinite(r.z):
            err = max(np.nanmean([abs(r.logmstar_ep), abs(r.logmstar_em)]), 0.1) if np.isfinite(r.logmstar_ep) else 0.1
            halo[k] = host_halo_dm(r.logmstar, err, r.z, r.offset_kpc, rng)
    h["dm_halo_p16"], h["dm_halo_pred"], h["dm_halo_p84"] = halo.T
    h["dm_ism"] = h[DM] - h.dm_halo_pred
    return h


def macquart_delta_draws(z, n, rng):
    """Draws of Delta = DM_cosmic / <DM_cosmic> from the Macquart et al. (2020) distribution."""
    sigma = MACQUART_F / np.sqrt(z)
    d = np.linspace(1e-3, 6.0, 20000)

    def pdf(c0):
        p = d ** -3 * np.exp(-((d ** -3 - c0) ** 2) / (18 * sigma ** 2))
        return p / np.trapezoid(p, d)

    lo, hi = -5.0, 5.0
    for _ in range(60):
        c0 = 0.5 * (lo + hi)
        lo, hi = (c0, hi) if np.trapezoid(d * pdf(c0), d) > 1 else (lo, c0)
    cdf = np.cumsum(pdf(c0))
    return np.interp(rng.random(n), cdf / cdf[-1], d)


def dm_host_draws(h, rng) -> np.ndarray:
    """(N_MC, N_hosts) draws of rest-frame DM_host."""
    out = np.full((N_MC, len(h)), np.nan)
    for k, r in enumerate(h.itertuples()):
        if not (np.isfinite(r.z) and np.isfinite(r.dm_obs) and np.isfinite(r.dm_mw)):
            continue
        cosmic = BT.mean_dm_cosmic(r.z) * macquart_delta_draws(r.z, N_MC, rng)
        mw = r.dm_mw * np.clip(rng.normal(1, DM_MW_FRAC_ERR, N_MC), 0, None)
        out[:, k] = (r.dm_obs - mw - cosmic) * (1 + r.z)
    return out


# ---------------------------------------------------------------- statistics
def partial_tau(x, y, z):
    """Kendall's partial tau of x and y given z (Akritas & Siebert 1996)."""
    txy, txz, tyz = (stats.kendalltau(a, b)[0] for a, b in ((x, y), (x, z), (y, z)))
    return (txy - txz * tyz) / np.sqrt((1 - txz ** 2) * (1 - tyz ** 2))


def partial_tau_test(x, y, z, rng, block=5):
    """Partial tau, its bootstrap 95% interval, and a p-value from shuffling y within redshift strata."""
    t = partial_tau(x, y, z)
    strata = np.array_split(np.argsort(z), max(len(z) // block, 1))
    null = np.empty(N_PERM)
    for k in range(N_PERM):
        yy = y.copy()
        for s in strata:
            yy[s] = y[rng.permutation(s)]
        null[k] = partial_tau(x, yy, z)
    boot = []
    for _ in range(1000):
        i = rng.integers(0, len(x), len(x))
        boot.append(partial_tau(x[i], y[i], z[i]))
    lo, hi = np.nanpercentile(boot, [2.5, 97.5])
    return t, lo, hi, float((np.sum(np.abs(null) >= abs(t)) + 1) / (N_PERM + 1))


def spearman_rows(x, Y):
    """Spearman r_S and p of x against every row of Y."""
    rx = stats.rankdata(x) - (len(x) + 1) / 2
    ry = stats.rankdata(Y, axis=1) - (len(x) + 1) / 2
    r = (ry @ rx) / np.sqrt((rx ** 2).sum() * (ry ** 2).sum(1))
    t = r * np.sqrt((len(x) - 2) / np.clip(1 - r ** 2, 1e-12, None))
    return r, 2 * stats.t.sf(np.abs(t), len(x) - 2)


def run_test(h, draws, theme, xcol, ycol, sample, claim, rng) -> dict:
    sub = SAMPLES[sample](h)
    sub = sub[np.isfinite(sub[xcol]) & np.isfinite(sub[ycol])]
    x, y = sub[xcol].to_numpy(float), sub[ycol].to_numpy(float)
    row = dict(theme=theme, x=xcol, y=ycol, sample=sample, claim=claim,
               **RC.robust_row(sub.frb.to_numpy(), x, y, rng))
    for d in (2, 3):
        lo, hi, pmin, pmax, fsig = E.jackknife(x, y, d)
        row.update({f"jk{d}_min": lo, f"jk{d}_max": hi, f"jk{d}_pmax": pmax, f"jk{d}_fsig": fsig})
    if "z" not in (xcol, ycol):
        z = sub.z.to_numpy(float)
        row["ptau"], row["ptau_lo"], row["ptau_hi"], row["p_ptau"] = partial_tau_test(x, y, z, rng)
    if DM in (xcol, ycol):
        other = x if ycol == DM else y
        dm = draws[:, sub.index.to_numpy()]
        r, p = spearman_rows(other, dm)
        row["mc_r16"], row["mc_r50"], row["mc_r84"] = np.percentile(r, [16, 50, 84])
        row["mc_fsig"] = float(np.mean(p < 0.05))
        p_ad = [RC.ad_split(*((other, dd) if ycol == DM else (dd, other)), exact=False) for dd in dm]
        row["mc_ad_fsig"] = float(np.mean(np.array(p_ad) < 0.05))
        fg = sub[np.isfinite(sub.dm_host_fg_excess)]
        row["rS_fg"], row["p_fg"] = stats.spearmanr(fg[xcol if xcol != DM else "dm_host_fg_excess"],
                                                    fg[ycol if ycol != DM else "dm_host_fg_excess"])
        row["N_fg"] = len(fg)
    return row


def sigma_sfr_slope(h, rng):
    """Theil-Sen slope of log DM_host against log Sigma_SFR (positive DM_host only); H-alpha EM predicts 0.5."""
    d = h[(h[DM] > 0) & np.isfinite(h.log_sigma_sfr)]
    x, y = d.log_sigma_sfr.to_numpy(), np.log10(d[DM].to_numpy())
    boot = [stats.theilslopes(y[i], x[i])[0] for i in (rng.integers(0, len(x), len(x)) for _ in range(2000))]
    return len(d), stats.theilslopes(y, x)[0], *np.percentile(boot, [2.5, 97.5])


# ---------------------------------------------------------------- report
def fmt_p(p):
    return "" if not np.isfinite(p) else (f"{p:.2f}" if p >= 0.1 else f"{p:.3f}" if p >= 0.001 else "<0.001")


def p_text(p):
    return "p < 0.001" if p < 0.001 else f"p = {fmt_p(p)}"


def table(d: pd.DataFrame) -> list[str]:
    lines = ["| x | y | sample | claim | N | r_S (p) | AD p (BH q) | 95% CI | tau (p) | JK1 | JK3: p<0.05 (max p) | "
             "most influential (r_S without) | partial tau given z (p) | DM_cosmic MC: r_S, p<0.05 (AD) | "
             "foreground-corrected | BH q |",
             "|" + "---|" * 16]
    for r in d.itertuples():
        ptau = f"{r.ptau:+.2f} ({fmt_p(r.p_ptau)})" if np.isfinite(r.ptau) else ""
        mc = (f"{r.mc_r50:+.2f} [{r.mc_r16:+.2f}, {r.mc_r84:+.2f}], {r.mc_fsig:.2f} ({r.mc_ad_fsig:.2f})"
              if np.isfinite(r.mc_r50) else "")
        fg = f"{r.rS_fg:+.2f} ({fmt_p(r.p_fg)})" if np.isfinite(r.rS_fg) else ""
        q = f"{r.q_BH:.2f}" if np.isfinite(r.q_BH) else ""
        ad = f"{fmt_p(r.p_AD)} ({r.q_BH_AD:.2f})" if np.isfinite(r.q_BH_AD) else fmt_p(r.p_AD)
        lines.append(f"| {LABELS[r.x]} | {LABELS[r.y]} | {r.sample} | {r.claim} | {r.N} | {r.r_S:+.2f} ({fmt_p(r.p_S)}) | "
                     f"{ad} | "
                     f"[{r.rS_boot_lo:+.2f}, {r.rS_boot_hi:+.2f}] | {r.tau:+.2f} ({fmt_p(r.p_tau)}) | "
                     f"{r.rS_jack_min:+.2f} to {r.rS_jack_max:+.2f} | {r.jk3_fsig:.2f} ({fmt_p(r.jk3_pmax)}) | "
                     f"{r.most_influential} ({r.rS_without_it:+.2f}) | {ptau} | {mc} | {fg} | {q} |")
    return lines


def figure(h, d):
    panels = [("dm_halo_pred", DM, "Predicted host-halo DM (pc cm$^{-3}$)"),
              ("log_sfr", DM, r"log SFR (M$_\odot$ yr$^{-1}$)"),
              ("log_sigma_sfr", DM, r"log $\Sigma_{\rm SFR}$ (M$_\odot$ yr$^{-1}$ kpc$^{-2}$)"),
              (DM, "log_tau", None), (DM, "log_abs_rm", None), ("z", "logmstar", "Redshift")]
    ylab = {DM: r"DM$_{\rm host}$ (rest frame, pc cm$^{-3}$)", "log_tau": r"log $\tau$ (ms)",
            "log_abs_rm": r"log |RM| (rad m$^{-2}$)", "logmstar": r"log $M_\star$ (M$_\odot$)"}
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.4))
    for ax, (xc, yc, xl), letter in zip(axes.flat, panels, "abcdef"):
        s = h[np.isfinite(h[xc]) & np.isfinite(h[yc])]
        sc = ax.scatter(s[xc], s[yc], c=s.z, cmap="viridis", vmin=0, vmax=0.8, s=18, edgecolor="k", lw=0.3)
        r = d[(d.x == xc) & (d.y == yc) & (d["sample"] == "All")].iloc[0]
        ax.text(0.03, 0.97, f"N = {r.N}\n$r_S$ = {r.r_S:+.2f} ({p_text(r.p_S)})\n"
                f"$\\tau$ = {r.tau:+.2f} ({p_text(r.p_tau)})\nAD {p_text(r.p_AD)}",
                transform=ax.transAxes, va="top", fontsize=7.5)
        ax.text(0.97, 0.03, f"({letter})", transform=ax.transAxes, ha="right", fontsize=9)
        ax.set_xlabel(xl or ylab[xc], fontsize=8.5)
        ax.set_ylabel(ylab[yc], fontsize=8.5)
        ax.tick_params(labelsize=7.5)
    lim = [0, max(h.dm_halo_pred.max(), 1) * 1.1]
    axes[0, 0].plot(lim, lim, "k--", lw=0.8)
    fig.colorbar(sc, ax=axes, shrink=0.6, label="Redshift")
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"dm_host_physics.{ext}", dpi=200, bbox_inches="tight")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    h = load().reset_index(drop=True)
    rng = np.random.default_rng(3)
    draws = dm_host_draws(h, rng)
    rows = []
    for t in TESTS:
        rows.append(run_test(h, draws, *t, rng))
        r = rows[-1]
        print(f"{t[1]:>20} {t[2]:>14} {t[3]:>16} N={r['N']:2d} r_S={r['r_S']:+.2f} p={r['p_S']:.3f}", flush=True)
    d = pd.DataFrame(rows)
    family = d.theme != "Sanity checks"
    d["q_BH"], d["q_BH_AD"] = np.nan, np.nan
    d.loc[family, "q_BH"] = stats.false_discovery_control(d.loc[family, "p_S"])
    d.loc[family, "q_BH_AD"] = stats.false_discovery_control(d.loc[family, "p_AD"])
    d.to_csv(OUT / "tests.csv", index=False)

    budget = h.loc[np.isfinite(h.dm_halo_pred), ["frb", "z", "logmstar", "offset_kpc", DM, "dm_halo_p16",
                                                 "dm_halo_pred", "dm_halo_p84", "dm_ism"]]
    budget.to_csv(OUT / "halo_budget.csv", index=False)
    n_s, slope, s_lo, s_hi = sigma_sfr_slope(h, rng)

    md = ["# What sets DM_host?", "",
          f"{len(h)} hosts; {np.isfinite(h[DM]).sum()} with DM_host. DM_host = (DM - DM_MW - median DM_cosmic)(1 + z). "
          "Sigma = quantity / (2 pi R_e^2), with the literature R_e where published and our GALFIT R_e otherwise; "
          "offset / R_e uses the literature offset; "
          "Delta MS is log SFR minus the Speagle+14 main sequence at the host redshift. "
          f"BH q over the {family.sum()} hypothesis tests (sanity checks excluded), separately for Spearman and AD. "
          "AD: two-sample Anderson-Darling p that y has the same distribution below and above the median of x "
          "(permutation p; inside the DM_cosmic MC, the fraction of draws with asymptotic AD p < 0.05).", ""]
    for theme in d.theme.unique():
        md += [f"## {theme}", ""] + table(d[d.theme == theme]) + [""]
    halo = budget.dropna(subset=[DM])
    md += ["## Host-halo budget", "",
           f"Predicted rest-frame DM from the FRB to the edge of its own halo (Moster+13 halo mass, mNFW gas, half chord "
           f"at the projected offset), N = {len(halo)}: median {halo.dm_halo_pred.median():.0f} pc cm^-3 "
           f"(16-84%: {halo.dm_halo_pred.quantile(0.16):.0f}-{halo.dm_halo_pred.quantile(0.84):.0f}); median DM_host "
           f"{halo[DM].median():.0f}; median halo share of DM_host (DM_host > 0) "
           f"{(halo.dm_halo_pred / halo[DM]).where(halo[DM] > 0).median():.2f}; hosts with DM_host below the halo "
           f"prediction: {(halo[DM] < halo.dm_halo_pred).sum()} / {len(halo)}.", "",
           "## Sigma_SFR slope", "",
           f"Theil-Sen slope of log DM_host vs log Sigma_SFR (DM_host > 0, N = {n_s}): {slope:.2f} "
           f"[{s_lo:.2f}, {s_hi:.2f}]; the H-alpha emission-measure argument predicts 0.5.", "",
           "## Main-sequence offsets", "",
           f"Delta MS (N = {np.isfinite(h.delta_ms).sum()}): median {h.delta_ms.median():+.2f} dex; within +-0.3 dex: "
           f"{(h.delta_ms.abs() <= 0.3).mean() / np.isfinite(h.delta_ms).mean():.2f}; more than 0.6 dex below: "
           f"{(h.delta_ms < -0.6).sum()}."]
    (OUT / "DM_HOST_PHYSICS.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    figure(h, d)
    print("\n".join(md))


if __name__ == "__main__":
    main()
