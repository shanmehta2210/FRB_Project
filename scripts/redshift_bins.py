"""Inclination tests, correlations with axis ratio, and host morphology in redshift bins.

1. Matched-draw AD / MW tests of cos i against the HSC control in each bin, against the full control and,
   for bins at z >= 0.2, against controls restricted to the same photometric-redshift range. The HSC
   catalog covers 0.2 < z < 1.0, so bins below z = 0.2 are compared with the full control and with its
   lowest-redshift slice (0.2 <= z_phot < 0.3).
2. Spearman r_S (p) and the AD face-on / edge-on halves test of every property against q in each bin.
3. Morphology class (disk_evidence.py) and Sersic index in each bin.

Output: Publication/materials/appendix/redshift_bins/{inclination_zbins.csv, correlations_zbins.csv,
        morphology_zbins.csv, REDSHIFT_BINS.md}
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import appendix_robustness as A  # noqa: E402
import plot_paper_inclination_figures as P  # noqa: E402
import robust_correlations as RC  # noqa: E402

OUT = P.REPO_ROOT / "Publication" / "materials" / "appendix" / "redshift_bins"
MORPH = P.REPO_ROOT / "Publication" / "materials" / "appendix" / "morphology" / "disk_evidence.csv"
N_DRAWS = 10_000
N_PERM = 4_999
BINS = [("All", -np.inf, np.inf), ("z < 0.1", -np.inf, 0.1), ("z < 0.2", -np.inf, 0.2),
        ("z >= 0.2", 0.2, np.inf), ("0.2 <= z < 0.4", 0.2, 0.4), ("z >= 0.4", 0.4, np.inf)]
PROPS = [("dm_host_rest", "DM_host"), ("dm_host_fg_excess", "DM_host, foreground excess removed"),
         ("log_abs_rm", "log|RM|"), ("log_tau", "log tau"), ("lin_pol", "Linear polarization"),
         ("av_old", "A_V, old stars"), ("logmstar", "log M*"), ("sfr", "SFR"), ("log_ssfr", "log sSFR")]
DISK = ["reliable", "lit", "profile"]


def hosts() -> pd.DataFrame:
    h = RC.load()
    lh = P.load_hosts()[["frb", "cosi_point", "re_over_fwhm"]]
    h = h.merge(lh, on="frb", how="inner")
    ev = pd.read_csv(MORPH, dtype={"frb": str})[["frb", "cls"]]
    h = h.merge(ev, on="frb", how="left")
    h["cls"] = h.cls.fillna("reliable")
    return h


def control() -> pd.DataFrame:
    cols = ["rmag", "ba", "n_sersic", "goodfits_flag", "photoz_best"]
    df = pd.read_csv(A.HSC_PARENT, usecols=cols).apply(pd.to_numeric, errors="coerce")
    ok = (df.goodfits_flag == 1) & np.isfinite(df.rmag) & np.isfinite(df.ba) & (df.ba > 0) & (df.ba <= 1) \
        & np.isfinite(df.n_sersic)
    pool = A.select(df.loc[ok].reset_index(drop=True))
    pool["cosi"] = A.cosi_from_q(pool.ba.to_numpy(), A.Q0)
    return pool


def ad_halves(x, y, rng):
    med = np.median(x)
    a, b = y[x > med], y[x <= med]
    if min(len(a), len(b)) < 4:
        return np.nan
    return float(stats.anderson_ksamp([a, b], method=stats.PermutationMethod(n_resamples=N_PERM,
                                                                            random_state=rng)).pvalue)


def inclination(h: pd.DataFrame, pool: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for lab, lo, hi in BINS:
        s = h[(h.z >= lo) & (h.z < hi)] if np.isfinite(lo) or np.isfinite(hi) else h
        if len(s) < 8:
            rows.append(dict(bin=lab, control="", N_frb=len(s)))
            continue
        ctrls = [("full", pool)]
        if lo >= 0.2:
            ctrls.append((f"z_phot in [{lo}, {hi})", pool[(pool.photoz_best >= lo) & (pool.photoz_best < hi)]))
        elif hi <= 0.2:
            ctrls.append(("z_phot in [0.2, 0.3)", pool[(pool.photoz_best >= 0.2) & (pool.photoz_best < 0.3)]))
        for clab, c in ctrls:
            r = A.matched_tests(s.cosi_point.to_numpy(float), c.cosi.to_numpy(float), None, N_DRAWS)
            rows.append(dict(bin=lab, control=clab, **r))
            print(lab, clab, r["N_frb"], r["N_pool"], r["quartile_counts"], round(r["ad_pm"], 3),
                  round(r["mw2_pm"], 3), round(r["mw1_pm"], 3), flush=True)
    return pd.DataFrame(rows)


def correlations(h: pd.DataFrame, rng) -> pd.DataFrame:
    rows = []
    for lab, lo, hi in BINS:
        s = h[(h.z >= lo) & (h.z < hi)]
        for col, plab in PROPS:
            d = s[np.isfinite(s.q) & np.isfinite(s[col])]
            row = dict(bin=lab, property=plab, N=len(d), r_S=np.nan, p_S=np.nan, p_AD=np.nan)
            if len(d) >= 6:
                row["r_S"], row["p_S"] = stats.spearmanr(d.q, d[col])
                row["p_AD"] = ad_halves(d.q.to_numpy(float), d[col].to_numpy(float), rng)
            rows.append(row)
    return pd.DataFrame(rows)


def morphology(h: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for lab, lo, hi in BINS:
        s = h[(h.z >= lo) & (h.z < hi)]
        c = s.cls.value_counts()
        rows.append(dict(bin=lab, N=len(s), **{k: int(c.get(k, 0)) for k in ["reliable", "lit", "profile", "ambig", "early"]},
                         disk_frac=float(s.cls.isin(DISK).mean()) if len(s) else np.nan,
                         median_n=float(s.sersic_n.median()), frac_n_lt_2p5=float((s.sersic_n < 2.5).mean()),
                         median_re_over_fwhm=float(s.re_over_fwhm.median()),
                         frac_resolved=float((s.re_over_fwhm >= 1).mean())))
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    h = hosts()
    hz = h[np.isfinite(h.z)]
    rng = np.random.default_rng(11)
    mo = morphology(hz)
    mo.to_csv(OUT / "morphology_zbins.csv", index=False)
    co = correlations(hz, rng)
    co.to_csv(OUT / "correlations_zbins.csv", index=False)
    inc = inclination(hz, control())
    inc.to_csv(OUT / "inclination_zbins.csv", index=False)

    disk = hz.cls.isin(DISK)
    trends = {
        "Sersic n vs z": stats.spearmanr(hz.z, hz.sersic_n),
        "R_e/FWHM vs z": stats.spearmanr(hz.z, hz.re_over_fwhm),
        "q vs z": stats.spearmanr(hz.z, hz.q),
    }
    mw = stats.mannwhitneyu(hz.z[disk], hz.z[~disk]) if (~disk).sum() else None

    md = ["# Redshift bins", "",
          f"{len(hz)} of {len(h)} hosts have a redshift. The HSC control covers 0.2 < z_phot < 1.0.", "",
          "## Host morphology", "",
          "Classes: reliable Sb--Irr (Gini--M20), literature disk, disk-like profile (n_free < 2.5), ambiguous, early-type.", "",
          "| bin | N | reliable | lit. disk | profile | ambiguous | early | disk fraction | median n | n < 2.5 | median R_e/FWHM | resolved |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    md += [f"| {r.bin} | {r.N} | {r.reliable} | {r.lit} | {r.profile} | {r.ambig} | {r.early} | {r.disk_frac:.2f} | "
           f"{r.median_n:.2f} | {r.frac_n_lt_2p5:.2f} | {r.median_re_over_fwhm:.2f} | {r.frac_resolved:.2f} |"
           for r in mo.itertuples()]
    md += ["", "Trends with redshift (Spearman): "
           + "; ".join(f"{k}: r_S = {v[0]:+.2f} (p = {v[1]:.3f})" for k, v in trends.items())
           + (f"; z of disks vs non-disks: MW p = {mw.pvalue:.3f}" if mw is not None else ""), "",
           "## Inclination vs control", "",
           "| bin | control | N_FRB | N_pool | quartile counts | AD p_m | MW 2-sided p_m | MW 1-sided p_m |",
           "|---|---|---|---|---|---|---|---|"]
    for r in inc.itertuples():
        if not r.control:
            md.append(f"| {r.bin} | too few hosts | {r.N_frb} | | | | | |")
            continue
        md.append(f"| {r.bin} | {r.control} | {r.N_frb} | {int(r.N_pool):,} | {r.quartile_counts} | {r.ad_pm:.3f} | "
                  f"{r.mw2_pm:.3f} | {r.mw1_pm:.3f} |")
    md += ["", "## Correlations with q", "", "r_S (p_S) [p_AD]; blank where N < 6.", "",
           "| property | " + " | ".join(b for b, _, _ in BINS) + " |", "|---|" + "---|" * len(BINS)]
    for col, plab in PROPS:
        cells = []
        for b, _, _ in BINS:
            r = co[(co.bin == b) & (co.property == plab)].iloc[0]
            cells.append("" if not np.isfinite(r.r_S) else
                         f"N={r.N}: {r.r_S:+.2f} ({r.p_S:.2f})" + (f" [{r.p_AD:.2f}]" if np.isfinite(r.p_AD) else ""))
        md.append(f"| {plab} | " + " | ".join(cells) + " |")
    (OUT / "REDSHIFT_BINS.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))


if __name__ == "__main__":
    main()
