"""Collect the main evidence for the three results of the paper into Publication/materials/major_summary.

  1_inclination_null        FRB host inclinations match field disks; no property depends on inclination
  2_correlations_null       correlations among z, M*, SFR and DM_host within bins of each other
  3_late_type_hosts         the hosts are overwhelmingly disk (late-type) galaxies

Existing figures are copied; tables are rendered as figures from the CSV outputs of the analysis scripts;
a few summary figures are drawn here. Run after the analysis scripts (binned_correlations.py last).
"""
from __future__ import annotations

import json
import re
import shutil
import sys
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from dm_host_physics import LABELS  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
MAT = ROOT / "Publication" / "materials"
APP = MAT / "appendix"
OUT = MAT / "major_summary"
INC, COR, MOR = OUT / "1_inclination_null", OUT / "2_correlations_null", OUT / "3_late_type_hosts"

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "stix", "font.serif": ["STIXGeneral"]})
RED, BLUE, GREY = "#D62728", "#1F77B4", "0.55"


# ---------------------------------------------------------------- helpers
def copy(src: Path, dst: Path):
    for ext in (".pdf", ".png"):
        if src.with_suffix(ext).exists():
            shutil.copy(src.with_suffix(ext), dst.with_suffix(ext))


def save(fig, path: Path):
    for ext in (".pdf", ".png"):
        fig.savefig(path.with_suffix(ext), dpi=200, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


def clean(s) -> str:
    s = "" if pd.isna(s) else str(s)
    s = re.sub(r"\\citep?\{([^}]*)\}", lambda m: "(" + m.group(1).replace(",", ";") + ")", s)
    s = re.sub(r"\\(le|ge)(?![a-z])", r"\\\1q", s)
    return (s.replace("--", "–").replace("`", "").replace("\\nodata", "–").replace(">=", "≥").replace("<=", "≤"))


def fmt_p(p) -> str:
    if pd.isna(p):
        return "–"
    return r"$\leq 10^{-4}$" if p <= 1e-4 else "< 0.001" if p < 0.001 else f"{p:.3f}" if p < 0.1 else f"{p:.2f}"


def table_figure(rows, header, path, title, note="", widths=None, fontsize=8, row_h=0.26):
    """Render rows (lists of strings) as a figure. A row with one cell is a section heading. widths in inches."""
    widths = np.asarray(widths if widths is not None else [1.2] * len(header), float)
    width_in = widths.sum() * 0.72 * fontsize / 8 + 0.2
    widths = widths / widths.sum()
    n = len(rows) + 1
    fig = plt.figure(figsize=(width_in, row_h * n + 0.9 + 0.16 * note.count("\n")))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.axis("off")
    top, h = 1 - 0.5 / fig.get_figheight(), row_h / fig.get_figheight()
    edges = np.r_[0, np.cumsum(widths)] * 0.98 + 0.01
    ax.text(0.01, 1 - 0.18 / fig.get_figheight(), title, fontsize=fontsize + 2, weight="bold", va="center")

    def put(y, cells, bold=False, shade=None):
        if shade:
            ax.add_patch(plt.Rectangle((0.01, y - h / 2), 0.98, h, color=shade, lw=0))
        if len(cells) == 1:
            ax.text(0.015, y, cells[0], fontsize=fontsize, style="italic", weight="bold", va="center")
            return
        for k, c in enumerate(cells):
            ax.text(edges[k] + 0.005, y, c, fontsize=fontsize, va="center", weight="bold" if bold else "normal",
                    wrap=False)

    put(top, header, bold=True, shade="#d9d9d9")
    for j, r in enumerate(rows):
        put(top - (j + 1) * h, r, shade="#f2f2f2" if (j % 2 and len(r) > 1) else ("#e8eef7" if len(r) == 1 else None))
    ax.plot([0.01, 0.99], [top - (n - 0.5) * h] * 2, color="k", lw=0.6)
    if note:
        ax.text(0.01, top - (n - 0.2) * h, note, fontsize=fontsize - 1, va="top", color="0.25")
    save(fig, path)


def inclination_tests(d: pd.DataFrame, label_col="label", n_col="N_frb"):
    return [[clean(getattr(r, label_col)), str(int(getattr(r, n_col))),
             f"{int(r.N_pool):,}" if np.isfinite(r.N_pool) else "–", clean(r.quartile_counts) or "–",
             fmt_p(r.ad_pm), fmt_p(r.mw2_pm), fmt_p(r.mw1_pm)] for r in d.itertuples()]


INC_HEADER = ["Variant", r"$N_{\rm FRB}$", r"$N_{\rm control}$", "Quartile counts (face-on → edge-on)",
              r"AD $p_{\rm m}$", r"MW two-sided $p_{\rm m}$", r"MW one-sided $p_{\rm m}$"]
INC_WIDTHS = [3.6, 0.8, 1.0, 2.2, 1.0, 1.4, 1.4]
PM_NOTE = (r"$p_{\rm m}$: fraction of 10,000 random control samples of the same size that differ from the control "
           "at least as much as the FRB hosts.\nThe one-sided MW test is run in the direction of an edge-on "
           r"deficit. With 10,000 draws the smallest possible $p_{\rm m}$ is $10^{-4}$.")
AD_NOTE = ("AD p: two-sample Anderson–Darling test that y (in 'y vs x'; the property, against $q$) has the same "
           "distribution for hosts below and above the median of x (permutation p, about 5,000 shuffles).")


# ---------------------------------------------------------------- 1: inclination
def pvalue_summary():
    rob = APP / "robustness"
    groups = []
    fv = pd.read_csv(APP / "fit_variants" / "fit_variants.csv")
    groups.append(("Adopted (GALFIT, PSF-convolved)", [("56 hosts vs HSC disks", fv.iloc[0])]))
    b24 = json.loads((MAT / "statistics" / "nature23_tests.json").read_text())
    groups.append(("Hosts of Bhardwaj et al. (2024), remeasured", [
        ("23 hosts, this work", dict(ad_pm=b24["ad_p_matched"], mw2_pm=np.nan, mw1_pm=b24["mw_p_matched"])),
        ("22 shared hosts, this work", fv.iloc[4])]))
    err = pd.read_csv(rob / "error_mc_realizations.csv")
    groups.append(("Measurement errors (2,000 realizations; median, 16–84%)",
                   [("q and 5° systematic redrawn", {k: err[k] for k in ("ad_pm", "mw2_pm", "mw1_pm")})]))
    ctrl = []
    for f in ("A2_sersic_window.csv", "A3_magnitude_cuts.csv", "A5_colour_cut.csv"):
        d = pd.read_csv(rob / f)
        d = d[~d.label.str.contains("adopted|r \\\\le 22", regex=True)]
        ctrl += [(clean(r.label), r) for _, r in d.iterrows()]
    groups.append(("Control-sample variants", ctrl))
    q0 = pd.read_csv(rob / "C1_q0.csv")
    groups.append((r"Intrinsic thickness $q_0$", [(rf"$q_0$ = {r.q0:.2f}", r) for _, r in q0.iterrows()]))
    ds = pd.read_csv(APP / "morphology" / "disk_subsample_tests.csv").iloc[1:]
    groups.append(("Host morphology subsets", [(clean(r.label), r) for _, r in ds.iterrows()]))
    zb = pd.read_csv(APP / "redshift_bins" / "inclination_zbins.csv").dropna(subset=["ad_pm"])
    groups.append(("Redshift bins (full control; control matched in z_phot)",
                   [(f"{r.bin} ({'full' if r.control == 'full' else 'matched'})", r) for _, r in zb.iterrows()]))
    groups.append(("For contrast: fits without the PSF, and the published inclinations",
                   [("GALFIT, no PSF", fv.iloc[2]), ("AstroPhot, no PSF", fv.iloc[3]),
                    ("Bhardwaj et al. (2024) published, 22 shared", fv.iloc[5])]))

    flat = []
    for g, items in groups:
        flat.append((g, None))
        flat += items
    fig, ax = plt.subplots(figsize=(7.6, 0.24 * len(flat) + 1.2))
    marks = (("ad_pm", "o", "AD"), ("mw2_pm", "s", "MW two-sided"), ("mw1_pm", "^", "MW one-sided (edge-on deficit)"))
    labels, contrast = [], False
    for j, (lab, r) in enumerate(flat):
        if r is None:
            labels.append(lab)
            contrast = "contrast" in lab
            ax.axhspan(j - 0.5, j + 0.5, color="#e8eef7", lw=0)
            continue
        labels.append("   " + lab)
        for k, (col, m, name) in enumerate(marks):
            v = r[col] if col in r else np.nan
            off = (k - 1) * 0.22
            if isinstance(v, pd.Series):
                lo, mid, hi = np.percentile(v, [16, 50, 84])
                ax.errorbar(max(mid, 1e-4), j + off, xerr=[[mid - lo], [hi - mid]], fmt=m, ms=4,
                            color=[RED, BLUE, "#2CA02C"][k], lw=0.8)
                continue
            if pd.isna(v):
                continue
            ax.plot(max(float(v), 1e-4), j + off, m, ms=4, color=GREY if contrast else [RED, BLUE, "#2CA02C"][k],
                    label=name if j == 1 else None)
    ax.axvline(0.05, color="k", ls="--", lw=0.8)
    ax.text(0.05, len(flat) - 0.6, " p = 0.05", fontsize=7, va="bottom")
    ax.set_xscale("log")
    ax.set_xlim(8e-5, 1.3)
    ax.set_ylim(len(flat) - 0.5, -1.5)
    ax.set_yticks(range(len(flat)), labels, fontsize=7)
    for t, (lab, r) in zip(ax.get_yticklabels(), flat):
        if r is None:
            t.set_weight("bold")
    ax.set_xlabel(r"$p_{\rm m}$ (fraction of random control samples at least as different)", fontsize=8)
    ax.tick_params(axis="x", labelsize=7)
    ax.legend(fontsize=7, loc="upper left", frameon=True)
    ax.set_title("Every inclination test: FRB hosts vs field disks (grey: no-PSF fits and published values, "
                 "for contrast)", fontsize=9)
    save(fig, INC / "03_pvalue_summary_all_tests")


def error_mc_hist():
    err = pd.read_csv(APP / "robustness" / "error_mc_realizations.csv")
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    bins = np.linspace(0, 1, 41)
    for col, lab, c in (("ad_pm", "AD", RED), ("mw2_pm", "MW two-sided", BLUE), ("mw1_pm", "MW one-sided", "#2CA02C")):
        ax.hist(err[col], bins=bins, histtype="step", lw=1.3, color=c,
                label=f"{lab}: median {err[col].median():.2f}; p < 0.05 in {np.mean(err[col] < 0.05):.1%}")
    ax.axvline(0.05, color="k", ls="--", lw=0.8)
    ax.set_xlabel(r"$p_{\rm m}$", fontsize=9)
    ax.set_ylabel("Realizations", fontsize=9)
    ax.set_title("Inclination tests through the measurement errors (2,000 realizations)", fontsize=9)
    ax.legend(fontsize=7, frameon=False)
    save(fig, INC / "04_measurement_error_realizations")


def marnoch_q():
    d = pd.read_csv(MAT / "marnoch_recheck" / "shared_q_comparison.csv", dtype={"frb": str})
    fig, ax = plt.subplots(figsize=(3.8, 3.6))
    ax.plot([0.2, 1], [0.2, 1], "k--", lw=0.7)
    ax.scatter(d.q_marnoch, d.q_ours, s=22, color=RED, edgecolor="k", lw=0.4, zorder=3)
    r, p = stats.spearmanr(d.q_marnoch, d.q_ours)
    ax.text(0.04, 0.96, f"N = {len(d)}\n$r_S$ = {r:.2f}\nmedian $|\\Delta q|$ = {d.dq.abs().median():.3f}",
            transform=ax.transAxes, va="top", fontsize=8)
    ax.set_xlabel(r"$q$, Marnoch et al. (2026)", fontsize=9)
    ax.set_ylabel(r"$q$, this work", fontsize=9)
    ax.set_title("Axis ratios of the shared hosts", fontsize=9)
    save(fig, INC / "13_marnoch_shared_axis_ratios")


def inclination_tables():
    rob = APP / "robustness"
    fv = pd.read_csv(APP / "fit_variants" / "fit_variants.csv")
    b24 = json.loads((MAT / "statistics" / "nature23_tests.json").read_text())
    pw = json.loads((MAT / "statistics" / "power_bhardwaj_deficit.json").read_text())
    err = pd.read_csv(rob / "error_mc_realizations.csv")
    ad = fv.iloc[0]
    pct = lambda c: f"{err[c].median():.2f} ({err[c].quantile(.16):.2f}–{err[c].quantile(.84):.2f})"
    rows = [[r"Median $\cos i$ (median $i$)", f"{ad.median_cosi_frb:.3f} ({ad.median_i_frb:.1f}°)",
             f"{ad.median_cosi_pool:.3f} ({ad.median_i_pool:.1f}°)"],
            [r"Fraction with $i > 60°$", f"{ad.frac_i_gt60_frb:.2f}", f"{ad.frac_i_gt60_pool:.2f}"],
            ["Quartile counts, face-on → edge-on (expected 14 each)", ad.quartile_counts, "–"],
            [r"AD / MW two-sided / MW one-sided $p_{\rm m}$, 56 hosts",
             f"{ad.ad_pm:.2f} / {ad.mw2_pm:.2f} / {ad.mw1_pm:.2f}", "–"],
            [r"Same, with measurement errors (median, 16–84%)",
             f"{pct('ad_pm')} / {pct('mw2_pm')} / {pct('mw1_pm')}", "–"],
            [r"Realizations with $p_{\rm m} < 0.05$ (AD / MW two-sided / one-sided)",
             f"{np.mean(err.ad_pm < .05):.1%} / {np.mean(err.mw2_pm < .05):.1%} / {np.mean(err.mw1_pm < .05):.1%}", "–"],
            [r"AD / MW one-sided $p_{\rm m}$, the 23 hosts of Bhardwaj et al. (2024) remeasured",
             f"{b24['ad_p_matched']:.2f} / {b24['mw_p_matched']:.2f}", "–"],
            ["Power: AD (MW) detects the published deficit at p < 0.05, samples of 56",
             f"{pw['N56']['ad_detect']:.2%} ({pw['N56']['mw_detect']:.2%})", "–"],
            ["Power: same, samples of 23",
             f"{pw['N23']['ad_detect']:.2%} ({pw['N23']['mw_detect']:.2%})", "–"]]
    table_figure(rows, ["Quantity", "FRB hosts", "Control (HSC disks)"], INC / "02b_table_headline_tests",
                 "Headline inclination tests", PM_NOTE + "\nPower: 2,000 samples drawn with replacement from the 23 "
                 "published inclinations of Bhardwaj et al. (2024).", widths=[5.4, 4.6, 1.8])

    ctrl = []
    for f, sec in (("A2_sersic_window.csv", "Sérsic-index range"), ("A3_magnitude_cuts.csv", "Magnitude limit"),
                   ("A5_colour_cut.csv", "Colour and star formation")):
        ctrl += [[sec]] + inclination_tests(pd.read_csv(rob / f))
    table_figure(ctrl, INC_HEADER, INC / "05_table_control_variants", "Sensitivity to the control-sample selection",
                 PM_NOTE + "\nEach block varies one criterion; the others are those of the adopted control "
                 r"($0.4 < n < 1.5$, $r \leq 22$, $q > 0.2$).", widths=INC_WIDTHS)
    q0 = pd.read_csv(rob / "C1_q0.csv")
    q0["label"] = [rf"$q_0$ = {v:.2f}" for v in q0.q0]
    table_figure(inclination_tests(q0), INC_HEADER, INC / "07_table_q0", r"Sensitivity to the intrinsic thickness $q_0$",
                 PM_NOTE + "\nThe cut $q > q_0$ is applied to both samples.", widths=INC_WIDTHS)
    fvt = [[clean(r.label), str(int(r.N_frb)), f"{r.median_i_frb:.1f}°", r.quartile_counts, fmt_p(r.ad_pm),
            fmt_p(r.mw2_pm), fmt_p(r.mw1_pm)] for r in fv.itertuples()]
    fvt.insert(4, ["Hosts shared with Bhardwaj et al. (2024)"])
    table_figure(fvt, ["Axis ratios", r"$N_{\rm FRB}$", "Median $i$", "Quartile counts", r"AD $p_{\rm m}$",
                       r"MW two-sided $p_{\rm m}$", r"MW one-sided $p_{\rm m}$"], INC / "09_table_psf_fit_variants",
                 "Fits with and without the PSF", PM_NOTE + "\nWithout the PSF, seeing circularizes marginally resolved "
                 "hosts and makes them look face-on; this reproduces the published deficit.", widths=INC_WIDTHS)
    zb = pd.read_csv(APP / "redshift_bins" / "inclination_zbins.csv")
    zb["label"] = [f"{b}, control {'full' if c == 'full' else c if isinstance(c, str) else '–'}"
                   for b, c in zip(zb.bin, zb.control)]
    table_figure(inclination_tests(zb), INC_HEADER, INC / "10_table_redshift_bins",
                 "Inclination tests in redshift bins", PM_NOTE + "\nMatched: control restricted to the same photometric-"
                 "redshift range. The z ≥ 0.4 bin (7 hosts) is too small to test.", widths=INC_WIDTHS)

    rc = pd.read_csv(APP / "robust_correlations" / "robust_correlations.csv")
    rct = [[clean(r.label), str(int(r.N)), f"{r.r_S:+.2f} ({fmt_p(r.p_S)})", fmt_p(r.p_AD),
            f"[{r.rS_boot_lo:+.2f}, {r.rS_boot_hi:+.2f}]",
            f"{r.tau:+.2f} ({fmt_p(r.p_tau)})", f"{r.rS_jack_min:+.2f} to {r.rS_jack_max:+.2f}",
            f"{r.theil_sen_slope:.3g} [{r.ts_lo:.3g}, {r.ts_hi:.3g}]"] for r in rc.itertuples()]
    table_figure(rct, ["Property vs axis ratio $q$", "N", r"$r_S$ (p)", "AD p", "bootstrap 95%", r"Kendall $\tau$ (p)",
                       "delete-1 range", "Theil–Sen slope [95%]"], INC / "12_table_correlations_with_axis_ratio",
                 "Host and burst properties vs axis ratio: outlier-robust statistics",
                 "Bootstrap: 10,000 resamples. Since every host has $q > q_0$, rank correlations with $q$ and with "
                 r"$\cos i$ are identical (sign flipped for $i$)." "\n" + AD_NOTE + "\nA$_V$ of the old stars is the only "
                 "Spearman p < 0.05 case (dust in inclined disks).",
                 widths=[4.4, 0.5, 1.2, 0.7, 1.3, 1.3, 1.3, 1.9], fontsize=7.5)
    cc = pd.read_csv(MAT / "marnoch_recheck" / "combined_correlations.csv")
    cct = [[clean(r.sample), str(int(r.N)), f"{r.r_S:+.2f} ({fmt_p(r.p_S)})", fmt_p(r.p_AD),
            f"[{r.boot_lo:+.2f}, {r.boot_hi:+.2f}]",
            f"{r.tau:+.2f} ({fmt_p(r.p_tau)})", f"{r.jack_min:+.2f} to {r.jack_max:+.2f}"] for r in cc.itertuples()]
    table_figure(cct, ["Sample", "N", r"$r_S$ (p)", "AD p", "bootstrap 95%", r"Kendall $\tau$ (p)", "delete-1 range"],
                 INC / "14_table_marnoch_combined", "DM$_{\\rm host}$ vs $q$: this work, Marnoch et al. (2026), and the union",
                 "Union: 69 unique hosts, Marnoch DM$_{\\rm host}$ shifted onto our definition; the 20 shared hosts take "
                 "either our values or theirs.\n" + AD_NOTE, widths=[3.6, 0.5, 1.3, 0.8, 1.4, 1.4, 1.4])
    zc = pd.read_csv(APP / "redshift_bins" / "correlations_zbins.csv")
    bins = list(dict.fromkeys(zc.bin))
    def cell(s):
        if not len(s) or not np.isfinite(s.r_S.iloc[0]):
            return "–" if not len(s) else f"– (N = {int(s.N.iloc[0])})"
        return (f"{s.r_S.iloc[0]:+.2f} ({fmt_p(s.p_S.iloc[0])}; {int(s.N.iloc[0])})\n"
                f"AD {fmt_p(s.p_AD.iloc[0])}")

    zct = [[clean(p)] + [cell(zc[(zc.property == p) & (zc.bin == b)]) for b in bins]
           for p in dict.fromkeys(zc.property)]
    table_figure(zct, ["Property vs $q$"] + [clean(b) for b in bins], INC / "15_table_correlations_with_q_in_redshift_bins",
                 "Correlations with axis ratio within redshift bins",
                 "Cells: Spearman $r_S$ (p; N), and below it the AD p (–: fewer than 4 hosts in one half).\n" + AD_NOTE,
                 widths=[2.6] + [1.6] * len(bins), fontsize=7.5, row_h=0.36)


def build_inclination():
    INC.mkdir(parents=True, exist_ok=True)
    fig, rob, bc = MAT / "figures", APP / "robustness", APP / "binned_correlations"
    copy(fig / "fig1_cdf.pdf", INC / "01_cdf_hosts_vs_control")
    copy(fig / "fig2_bins.pdf", INC / "02_quartile_bins")
    inclination_tables()
    pvalue_summary()
    error_mc_hist()
    for k, name in enumerate(("sersic", "mag", "color")):
        copy(rob / f"figA_ctrl_{name}.pdf", INC / f"06{'abc'[k]}_control_variants_{name}")
    copy(rob / "figC_q0.pdf", INC / "08_q0_cdfs")
    copy(rob / "figA_control_hist.pdf", INC / "11_control_cosi_histogram")
    copy(fig / "fig3_correlations.pdf", INC / "12a_dm_host_vs_axis_ratio")
    marnoch_q()
    copy(bc / "grid_inclination.pdf", INC / "16_binned_properties_vs_inclination_grid")
    copy(bc / "heatmap_inclination.pdf", INC / "17_binned_properties_vs_inclination_heatmap")
    copy(bc / "global_tests.pdf", INC / "18_global_tests_both_grids")


# ---------------------------------------------------------------- 2: correlations
def correlation_tables():
    t = pd.read_csv(APP / "dm_host_physics" / "tests.csv")
    rows, theme = [], None
    for r in t.itertuples():
        if r.theme != theme:
            theme = r.theme
            rows.append([clean(theme)])
        partial = f"{r.ptau:+.2f} ({fmt_p(r.p_ptau)})" if np.isfinite(r.ptau) else "–"
        mc = f"{r.mc_r50:+.2f} ({r.mc_fsig:.2f}; {r.mc_ad_fsig:.2f})" if np.isfinite(r.mc_r50) else "–"
        rows.append([f"{LABELS[r.y]} vs {LABELS[r.x]}", clean(r.sample), str(int(r.N)), f"{r.r_S:+.2f} ({fmt_p(r.p_S)})",
                     fmt_p(r.p_AD), f"[{r.rS_boot_lo:+.2f}, {r.rS_boot_hi:+.2f}]",
                     f"{r.jk3_fsig:.2f}" if np.isfinite(r.jk3_fsig) else "–", partial, mc,
                     f"{r.q_BH:.2f}" if np.isfinite(r.q_BH) else "–",
                     f"{r.q_BH_AD:.2f}" if np.isfinite(r.q_BH_AD) else "–"])
    table_figure(rows, ["Pair", "Sample", "N", r"$r_S$ (p)", "AD p", "bootstrap 95%", "delete-3: frac. p<0.05",
                        r"partial $\tau$ | z (p)", "DM$_{\\rm cosmic}$ MC $r_S$ (frac.; AD)", "BH q", "BH q (AD)"],
                 COR / "04_table_dm_host_tests", "What DM$_{\\rm host}$ and the host properties correlate with (whole sample)",
                 "BH q: Benjamini–Hochberg over the 34 hypothesis tests (sanity checks excluded), for Spearman and for AD. "
                 "MC: median $r_S$ over 2,000 draws of DM$_{\\rm cosmic}$ and DM$_{\\rm MW}$, and the fractions with "
                 "Spearman p < 0.05 and AD p < 0.05.\n" + AD_NOTE,
                 widths=[3.4, 1.7, 0.45, 1.25, 0.7, 1.3, 1.6, 1.4, 2.0, 0.6, 0.8], fontsize=7, row_h=0.22)


def build_correlations():
    COR.mkdir(parents=True, exist_ok=True)
    bc = APP / "binned_correlations"
    copy(bc / "grid.pdf", COR / "01_binned_pairs_grid")
    copy(bc / "heatmap.pdf", COR / "02_binned_pairs_heatmap")
    copy(bc / "binned_by_z__dm_host_rest_vs_log_sfr.png", COR / "03_dm_host_vs_sfr_binned_by_z")
    correlation_tables()
    copy(APP / "dm_host_physics" / "dm_host_physics.pdf", COR / "05_dm_host_physics_panels")
    copy(MAT / "foreground_halos" / "foreground_diagnostic.pdf", COR / "06_foreground_halo_dm_per_frb")
    copy(MAT / "foreground_halos" / "audit" / "audit.pdf", COR / "08_foreground_halo_audit")
    copy(bc / "global_tests.pdf", COR / "07_global_tests_both_grids")


# ---------------------------------------------------------------- 3: morphology
CLASS_ORDER = [("reliable", "Sb–Irr in Gini–$M_{20}$ (reliable)", "#2a9d8f"),
               ("lit", "Disk or spiral in the literature", "#4c78a8"),
               ("profile", "Disk-like profile ($n_{\\rm free}$ < 2.5)", "#9ecae9"),
               ("ambig", "Ambiguous", "#bbbbbb"), ("early", "Early-type", RED)]


def host_classes():
    seed = pd.read_csv(MAT / "host_properties" / "seed_56.csv", dtype={"frb": str})
    ev = pd.read_csv(APP / "morphology" / "disk_evidence.csv", dtype={"frb": str})
    seed["cls"] = seed.frb.map(ev.set_index("frb").cls).fillna("reliable")
    seed["sersic_n_free"] = seed.n_galfit
    return seed


def morphology_overview():
    h = host_classes()
    gm = pd.read_csv(APP / "robustness" / "gini_m20_56.csv", dtype={"frb": str})
    h = h.merge(gm[["frb", "re_over_fwhm"]], on="frb", how="left")
    counts = h.cls.value_counts()
    fig, (a, b) = plt.subplots(1, 2, figsize=(9.2, 3.4), gridspec_kw=dict(width_ratios=[1.05, 1]))
    left = 0
    for key, lab, col in CLASS_ORDER:
        n = int(counts.get(key, 0))
        a.barh(0, n, left=left, color=col, edgecolor="white", label=f"{lab}: {n}")
        if n:
            a.text(left + n / 2, 0, str(n), ha="center", va="center", fontsize=9, color="k")
        left += n
    disk = int(sum(counts.get(k, 0) for k in ("reliable", "lit", "profile")))
    a.axvline(disk, color="k", lw=1)
    a.text(disk / 2, 0.55, f"disk / late-type: {disk} of {len(h)}", ha="center", fontsize=9, weight="bold")
    a.set_xlim(0, len(h))
    a.set_ylim(-0.6, 0.8)
    a.set_yticks([])
    a.set_xlabel("Number of hosts", fontsize=9)
    a.legend(fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=2, frameon=False)
    a.set_title("(a) Morphological class of the 56 hosts", fontsize=9)
    bins = np.arange(0, 6.51, 0.5)
    res = h.re_over_fwhm >= 1
    b.hist([h.sersic_n_free[res], h.sersic_n_free[~res]], bins=bins, stacked=True, color=["#4c78a8", "#9ecae9"],
           edgecolor="white", label=[r"resolved ($R_e$/FWHM ≥ 1)", r"marginally resolved ($R_e$/FWHM < 1)"])
    b.axvline(2.5, color="k", ls="--", lw=0.8)
    frac = np.mean(h.sersic_n_free < 2.5)
    b.text(2.6, b.get_ylim()[1] * 0.92, f"n < 2.5: {frac:.0%}\nmedian n = {h.sersic_n_free.median():.2f}",
           fontsize=8, va="top")
    b.set_xlabel("Sérsic index $n$ (PSF-convolved GALFIT, $n$ free)", fontsize=9)
    b.set_ylabel("Hosts", fontsize=9)
    b.legend(fontsize=7, frameon=False, loc="upper right", bbox_to_anchor=(1, 0.75))
    b.set_title("(b) Sérsic index", fontsize=9)
    save(fig, MOR / "01_morphology_overview")


def morphology_tables():
    ev = pd.read_csv(APP / "morphology" / "disk_evidence.csv", dtype={"frb": str})
    order = {k: j for j, (k, _, _) in enumerate(CLASS_ORDER)}
    ev = ev.assign(o=ev.cls.map(order)).sort_values(["o", "frb"])
    gm = {"late": "Sb–Irr", "early": "E/S0/Sa", "merger": "Merger"}
    rows = [[r.frb, f"{r.n_free:.1f}", f"{r.n_statmorph:.1f}", f"{r.re_over_fwhm:.2f}", gm.get(r.gm20, r.gm20),
             clean(r.label), textwrap.shorten(clean(r.literature), 100, placeholder=" …")] for r in ev.itertuples()]
    table_figure(rows, ["FRB", r"$n_{\rm free}$", r"$n_{\rm sm}$", r"$R_e$/FWHM", "G–M$_{20}$", "Class", "Literature"],
                 MOR / "03_table_disk_evidence",
                 "Host-by-host evidence for the hosts with unreliable or non-Sb–Irr Gini–M$_{20}$",
                 "$n_{\\rm free}$: PSF-convolved GALFIT with $n$ free; $n_{\\rm sm}$: statmorph Sérsic fit. The other 24 "
                 "hosts are reliable Sb–Irr in Gini–M$_{20}$.", widths=[0.9, 0.5, 0.5, 0.75, 0.8, 1.4, 6.0],
                 fontsize=7, row_h=0.22)
    ds = pd.read_csv(APP / "morphology" / "disk_subsample_tests.csv")
    table_figure(inclination_tests(ds), INC_HEADER, MOR / "05_table_inclination_tests_disk_subsets",
                 "Inclination tests restricted to the most secure disks", PM_NOTE, widths=INC_WIDTHS)
    mz = pd.read_csv(APP / "redshift_bins" / "morphology_zbins.csv")
    rows = [[r.bin, str(int(r.N)), str(int(r.reliable)), str(int(r.lit)), str(int(r.profile)), str(int(r.ambig)),
             str(int(r.early)), f"{r.disk_frac:.2f}", f"{r.frac_n_lt_2p5:.2f}", f"{r.median_re_over_fwhm:.2f}"]
            for r in mz.itertuples()]
    table_figure(rows, ["Redshift bin", "N", "Sb–Irr (rel.)", "Lit. disk", "Disk profile", "Ambiguous", "Early",
                        "Disk fraction", "n < 2.5", "median $R_e$/FWHM"], MOR / "06_table_morphology_by_redshift",
                 "Host morphology by redshift", "55 hosts with a redshift. Disk fraction: reliable Sb–Irr, literature "
                 "disks and disk-like profiles.", widths=[1.4, 0.5, 1, 0.8, 1, 0.9, 0.6, 1, 0.7, 1.3])


def build_morphology():
    MOR.mkdir(parents=True, exist_ok=True)
    morphology_overview()
    copy(APP / "robustness" / "figB_gini_m20.pdf", MOR / "02_gini_m20")
    morphology_tables()
    for k in range(1, 5):
        copy(ROOT / "Publication" / "manuscript" / "galaxy" / f"gallery_{k}.pdf", MOR / f"04{'abcd'[k - 1]}_gallery_{k}")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    build_inclination()
    build_correlations()
    build_morphology()
    print("\n".join(sorted(str(p.relative_to(OUT)) for p in OUT.rglob("*.png"))))


if __name__ == "__main__":
    main()
