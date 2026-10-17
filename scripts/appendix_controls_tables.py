"""Numbers behind the appendix control figures (figA_q0, figA_nrange).

Reuses the loaders of ``plot_paper_inclination_figures.py`` so the samples are
identical to the figures.  Writes:

  Publication/materials/appendix/controls_q0.csv
  Publication/materials/appendix/controls_nrange.csv
  Publication/materials/appendix/CONTROLS_Q0_NRANGE.md
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import plot_paper_inclination_figures as P  # noqa: E402
from null_catalog_utils import hubble_cosi_from_ba  # noqa: E402

OUT = P.REPO_ROOT / "Publication" / "materials" / "appendix"
Q0_LIST = [0.10, 0.15, 0.20, 0.25]
N_WINDOWS = [(0.4, 1.5), (0.5, 1.2), (0.8, 1.2), (0.4, 2.0)]
N_MATCHED = 10_000


def tests(frb: np.ndarray, hsc: np.ndarray, rng: np.random.Generator) -> dict:
    ks = stats.ks_2samp(frb, hsc, alternative="two-sided")
    ks1 = stats.ks_2samp(frb, hsc, alternative="less")  # FRB CDF below HSC => more edge-on
    mw = stats.mannwhitneyu(frb, hsc, alternative="two-sided")
    mw1 = stats.mannwhitneyu(frb, hsc, alternative="greater")
    # matched-N KS: p from the sampling distribution of D for N-sized HSC draws
    n = len(frb)
    d_obs = float(ks.statistic)
    d_null = np.empty(N_MATCHED)
    for k in range(N_MATCHED):
        sub = rng.choice(hsc, size=n, replace=False)
        d_null[k] = stats.ks_2samp(sub, hsc).statistic
    return {
        "N_frb": n, "N_hsc": len(hsc),
        "median_cosi_frb": float(np.median(frb)), "median_cosi_hsc": float(np.median(hsc)),
        "mean_cosi_frb": float(np.mean(frb)), "mean_cosi_hsc": float(np.mean(hsc)),
        "frac_i_gt_60_frb": float(np.mean(frb < 0.5)), "frac_i_gt_60_hsc": float(np.mean(hsc < 0.5)),
        "KS_D": d_obs, "KS_p_two_sided": float(ks.pvalue), "KS_p_one_sided_edgeon": float(ks1.pvalue),
        "KS_p_matchedN": float(np.mean(d_null >= d_obs)),
        "MWU_p_two_sided": float(mw.pvalue), "MWU_p_one_sided_edgeon": float(mw1.pvalue),
        "Welch_t_p": float(stats.ttest_ind(frb, hsc, equal_var=False).pvalue),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    hosts = P.load_hosts()
    n_hosts = len(hosts)
    rng = np.random.default_rng(P.SEED + 11)

    # ---- q0 ---------------------------------------------------------------
    parent = P.load_hsc_ba(0.4, 1.5, ba_cut=False)
    qseries = P.magnitude_q(hosts)
    rows = []
    for q0 in Q0_LIST:
        frb_q = qseries[qseries > q0]
        hsc_q = parent.loc[parent["ba"] > q0, "ba"].to_numpy(dtype=float)
        frb_c = np.array([hubble_cosi_from_ba(float(v), q0=q0) for v in frb_q])
        hsc_c = np.array([hubble_cosi_from_ba(float(v), q0=q0) for v in hsc_q])
        added = qseries[(qseries > q0) & (qseries <= 0.2)]
        removed = qseries[(qseries <= q0) & (qseries > 0.2)]
        r = {"q0": q0, "hosts_added_vs_0.2": ", ".join(f"{i} ({v:.2f})" for i, v in added.items()),
             "hosts_removed_vs_0.2": ", ".join(f"{i} ({v:.2f})" for i, v in removed.items())}
        r.update(tests(frb_c, hsc_c, rng))
        rows.append(r)
    q0df = pd.DataFrame(rows)
    q0df.to_csv(OUT / "controls_q0.csv", index=False)

    # ---- Sersic-n window of the HSC control -------------------------------
    frb = hosts["cosi_point"].to_numpy(dtype=float)
    rows = []
    for lo, hi in N_WINDOWS:
        hsc = P.load_hsc_ba(lo, hi, ba_cut=True)
        cosi = np.array([hubble_cosi_from_ba(float(v), q0=P.Q0) for v in hsc["ba"]])
        r = {"n_min": lo, "n_max": hi, "median_n_hsc": float(np.median(hsc["n"]))}
        r.update(tests(frb, cosi, rng))
        rows.append(r)
    ndf = pd.DataFrame(rows)
    ndf.to_csv(OUT / "controls_nrange.csv", index=False)

    # ---- markdown ---------------------------------------------------------
    def f3(x: float) -> str:
        return f"{x:.3f}"

    md = ["# Appendix controls: intrinsic thickness $q_0$ and the HSC Sérsic-index window", "",
          "Source: `scripts/appendix_controls_tables.py`, which reuses the loaders of "
          "`scripts/plot_paper_inclination_figures.py` (figures `figA_q0`, `figA_nrange`). "
          "Point estimates only (winning-leg $q$, no error smearing). HSC control: Kawinwanichakij et al. "
          "parent sample, `goodfits_flag = 1`, $r \\le 22$. One-sided $p$ tests the direction claimed in 2024 "
          "(FRB hosts more edge-on, i.e. smaller $\\cos i$). Matched-$N$ KS $p$ is the fraction of 10,000 "
          "$N$-sized HSC subsamples with $D \\ge D_{\\rm obs}$.", "",
          "## 1. Intrinsic axis ratio $q_0$", "",
          "$\\cos^2 i = (q^2 - q_0^2)/(1 - q_0^2)$ applied to both samples with the same $q_0$; the $q > q_0$ "
          "cut is re-applied to both, so the FRB sample size changes with $q_0$ (hosts with $q \\le 0.2$ that "
          "pass the magnitude and verification cuts re-enter at lower $q_0$).", "",
          "| $q_0$ | $N_{\\rm FRB}$ | $N_{\\rm HSC}$ | hosts added / removed vs 0.2 | median $\\cos i$ FRB / HSC | "
          "$f(i>60^\\circ)$ FRB / HSC | KS $D$ | KS $p$ (2-sided) | KS $p$ (edge-on) | KS $p$ matched-$N$ | "
          "MWU $p$ (2-sided) | MWU $p$ (edge-on) |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in q0df.iterrows():
        added, removed = r["hosts_added_vs_0.2"], r["hosts_removed_vs_0.2"]
        chg = ((f"+ {added}" if added else "") + (f" − {removed}" if removed else "")).strip() or "—"
        md.append(f"| {r.q0:.2f} | {r.N_frb} | {r.N_hsc:,} | {chg} | {f3(r.median_cosi_frb)} / {f3(r.median_cosi_hsc)} | "
                  f"{f3(r.frac_i_gt_60_frb)} / {f3(r.frac_i_gt_60_hsc)} | {f3(r.KS_D)} | {f3(r.KS_p_two_sided)} | "
                  f"{f3(r.KS_p_one_sided_edgeon)} | {f3(r.KS_p_matchedN)} | {f3(r.MWU_p_two_sided)} | {f3(r.MWU_p_one_sided_edgeon)} |")
    md += ["", "## 2. Sérsic-index window of the HSC control", "",
           f"FRB sample fixed ({n_hosts}, $q_0 = 0.2$). The HSC pool is re-selected with different $n$ windows "
           "(the adopted window is $0.4 < n < 1.5$).", "",
           "| $n$ window | $N_{\\rm HSC}$ | median $n$ | median $\\cos i$ HSC | $f(i>60^\\circ)$ HSC | KS $D$ | "
           "KS $p$ (2-sided) | KS $p$ (edge-on) | KS $p$ matched-$N$ | MWU $p$ (2-sided) | MWU $p$ (edge-on) |",
           "|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in ndf.iterrows():
        md.append(f"| ${r.n_min:.1f} < n < {r.n_max:.1f}$ | {int(r.N_hsc):,} | {r.median_n_hsc:.2f} | {f3(r.median_cosi_hsc)} | "
                  f"{f3(r.frac_i_gt_60_hsc)} | {f3(r.KS_D)} | {f3(r.KS_p_two_sided)} | {f3(r.KS_p_one_sided_edgeon)} | "
                  f"{f3(r.KS_p_matchedN)} | {f3(r.MWU_p_two_sided)} | {f3(r.MWU_p_one_sided_edgeon)} |")
    md += ["", f"FRB sample: median $\\cos i$ = {np.median(frb):.3f}, $f(i>60^\\circ)$ = {np.mean(frb < 0.5):.3f}.", ""]
    (OUT / "CONTROLS_Q0_NRANGE.md").write_text("\n".join(md), encoding="utf-8")
    sys.stdout.reconfigure(encoding="utf-8")
    print("\n".join(md))


if __name__ == "__main__":
    main()
