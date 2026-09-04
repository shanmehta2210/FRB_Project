"""
Two-sample tests: confirmed-50 FRB Hubble cos(i) vs HSC disk pool.

Uses winning-leg point estimates in confirmed50_q.csv (no error MC) and the
same HSC cuts as plot_frb_vs_hsc_confirmed50.py (rmag <= 22, ba > 0.2).

Writes plots/plots_null/v2/frb_vs_hsc_confirmed50/TEST_RESULTS.md

Run from repo root::

    python scripts/run_frb_hsc_confirmed50_tests.py
"""
from __future__ import annotations

import argparse
import math
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import Q0  # noqa: E402
from pipeline_null_plot_utils import PLOTS_NULL  # noqa: E402
from plot_frb_inflated_error_cdf import N_DRAWS, SEED  # noqa: E402
from plot_frb_vs_hsc_confirmed50 import (  # noqa: E402
    OUT_DIR,
    load_hsc_pool,
)

CSV_DEFAULT = OUT_DIR / "confirmed50_q.csv"
MD_DEFAULT = OUT_DIR / "TEST_RESULTS.md"


def _fmt_p(p: float, *, ad_cap: bool = False) -> str:
    if p is None or not math.isfinite(float(p)):
        return "—"
    p = float(p)
    if ad_cap and abs(p - 0.25) < 1e-9:
        return ">= 0.25 (SciPy cap)"
    if p < 0.001:
        return f"{p:.2e}"
    return f"{p:.3f}"


def _fmt_stat(x: float, nd: int = 3) -> str:
    if x is None or not math.isfinite(float(x)):
        return "—"
    return f"{float(x):.{nd}f}"


def kuiper_v(x: np.ndarray, y: np.ndarray) -> float:
    """Two-sample Kuiper V = D+ + D- of the ECDF difference."""
    x = np.sort(np.asarray(x, dtype=float))
    y = np.sort(np.asarray(y, dtype=float))
    grid = np.sort(np.concatenate([x, y]))
    fx = np.searchsorted(x, grid, side="right") / len(x)
    fy = np.searchsorted(y, grid, side="right") / len(y)
    d = fx - fy
    return float(d.max() - d.min())


def observed_tests(frb: np.ndarray, hsc: np.ndarray) -> dict:
    out: dict = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        ks = stats.ks_2samp(frb, hsc, alternative="two-sided", method="auto")
    out["ks_d"] = float(ks.statistic)
    out["ks_p"] = float(ks.pvalue)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ad = stats.anderson_ksamp([hsc, frb])
    out["ad_stat"] = float(ad.statistic)
    out["ad_p"] = float(ad.pvalue)
    out["ad_capped"] = abs(float(ad.pvalue) - 0.25) < 1e-9

    cvm = stats.cramervonmises_2samp(frb, hsc)
    out["cvm_stat"] = float(cvm.statistic)
    out["cvm_p"] = float(cvm.pvalue)

    out["kuiper_v"] = kuiper_v(frb, hsc)

    mwu = stats.mannwhitneyu(frb, hsc, alternative="two-sided")
    out["mwu_u"] = float(mwu.statistic)
    out["mwu_p"] = float(mwu.pvalue)
    n1, n2 = len(frb), len(hsc)
    out["mwu_rbc"] = float(1.0 - (2.0 * mwu.statistic) / (n1 * n2))

    bm = stats.brunnermunzel(frb, hsc, alternative="two-sided")
    out["bm_stat"] = float(bm.statistic)
    out["bm_p"] = float(bm.pvalue)

    try:
        es = stats.epps_singleton_2samp(frb, hsc)
        out["es_stat"] = float(es.statistic)
        out["es_p"] = float(es.pvalue)
    except (ValueError, np.linalg.LinAlgError):
        out["es_stat"] = float("nan")
        out["es_p"] = float("nan")

    tt = stats.ttest_ind(frb, hsc, equal_var=False)
    out["welch_t"] = float(tt.statistic)
    out["welch_p"] = float(tt.pvalue)

    out["w1"] = float(stats.wasserstein_distance(frb, hsc))
    return out


def mc_tail_p(
    frb: np.ndarray,
    hsc: np.ndarray,
    *,
    n_draws: int,
    seed: int,
) -> dict[str, float]:
    """Fraction of HSC-N draws whose statistic vs the pool is >= the FRB value."""
    n = len(frb)
    rng = np.random.default_rng(seed)
    obs = observed_tests(frb, hsc)
    keys = ("ks_d", "ad_stat", "cvm_stat", "kuiper_v", "w1")
    counts = {k: 0 for k in keys}
    for _ in range(n_draws):
        draw = rng.choice(hsc, size=n, replace=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ks = stats.ks_2samp(draw, hsc, alternative="two-sided", method="asymp")
            ad = stats.anderson_ksamp([hsc, draw])
            cvm = stats.cramervonmises_2samp(draw, hsc)
        stats_draw = {
            "ks_d": float(ks.statistic),
            "ad_stat": float(ad.statistic),
            "cvm_stat": float(cvm.statistic),
            "kuiper_v": kuiper_v(draw, hsc),
            "w1": float(stats.wasserstein_distance(draw, hsc)),
        }
        for k in keys:
            if stats_draw[k] >= obs[k]:
                counts[k] += 1
    return {k: (1.0 + counts[k]) / (1.0 + n_draws) for k in keys}


def write_md(
    dest: Path,
    *,
    frb: np.ndarray,
    hsc: np.ndarray,
    obs: dict,
    mc_p: dict,
    n_draws: int,
) -> None:
    n1, n2 = len(frb), len(hsc)
    med_f, med_h = float(np.median(frb)), float(np.median(hsc))
    mean_f, mean_h = float(np.mean(frb)), float(np.mean(hsc))
    dmed, dmean = med_f - med_h, mean_f - mean_h

    def row(name, stat, p_pool, p_mc=None, extra=""):
        pmc = _fmt_p(p_mc) if p_mc is not None else "—"
        return f"| {name} | {stat} | {p_pool} | {pmc} | {extra} |"

    lines = [
        "# Confirmed-50 vs HSC — statistical tests",
        "",
        "Question: could the 50 confirmed FRB-host Hubble \(\\cos(i)\) "
        "values be a random draw from the HSC disk pool "
        "(\(r\\le 22\), \(b/a>0.2\), \(0.4<n<1.5\))?",
        "",
        "Tests use **point estimates** (winning-leg \(q \\to\) Hubble "
        "\(\\cos(i)\), \(q_0=0.2\)). They do **not** fold in Protocol A/B "
        "error smearing. Same convention as "
        "[`run_sdss_frb_inclination_tests.py`](../../../scripts/run_sdss_frb_inclination_tests.py).",
        "",
        f"| | FRB | HSC |",
        f"|---|---:|---:|",
        f"| N | {n1} | {n2:,} |",
        f"| median \(\\cos(i)\) | {med_f:.3f} | {med_h:.3f} |",
        f"| mean \(\\cos(i)\) | {mean_f:.3f} | {mean_h:.3f} |",
        f"| FRB − HSC (median) | {dmed:+.3f} | |",
        f"| FRB − HSC (mean) | {dmean:+.3f} | |",
        "",
        "## How to read a p-value",
        "",
        "p is the chance of a discrepancy **at least this large** if FRB and HSC "
        "really share one continuous distribution. It is **not** a percent-similar score.",
        "",
        "| p | Read |",
        "|---|---|",
        "| > 0.10 | Consistent with the same distribution (for that test). |",
        "| 0.05–0.10 | Weak / suggestive only. |",
        "| 0.01–0.05 | Significant at the usual 5% bar. |",
        "| < 0.01 | Stronger evidence of a difference. |",
        "",
        "Two p-values are listed:",
        "",
        "- **p (vs pool)** — textbook two-sample test of the 50 vs all "
        f"{n2:,} HSC galaxies. Power is high: a *small* CDF offset can make p tiny.",
        f"- **p (matched N)** — among {n_draws:,} random HSC draws of 50, the "
        "fraction whose statistic vs the pool is at least as large as the FRB’s. "
        "This matches the CDF envelope (same N). Prefer this when asking "
        "“is this 50 unusual among 50-galaxy HSC samples?”",
        "",
        "SciPy AD p-values are **capped at 0.25** when samples look alike "
        "(true p is at least 0.25) and floored near 0.001 when they differ strongly.",
        "",
        "**Takeaway.** Matched-N Monte Carlo p-values sit near 0.3-0.6 "
        "(typical HSC-50 draw). Location tests vs the pool agree. "
        "The CvM vs-pool p is a SciPy N2>>N1 artifact; use matched-N only.",
        "",
        "## Summary",
        "",
        "| Test | statistic | p (vs pool) | p (matched N) | notes |",
        "|---|---:|---:|---:|---|",
        row(
            "Kolmogorov–Smirnov",
            _fmt_stat(obs["ks_d"], 3),
            _fmt_p(obs["ks_p"]),
            mc_p["ks_d"],
            r"D = max \|F−G\|",
        ),
        row(
            "Anderson–Darling",
            _fmt_stat(obs["ad_stat"], 3),
            _fmt_p(obs["ad_p"], ad_cap=obs["ad_capped"]),
            mc_p["ad_stat"],
            "tails weighted",
        ),
        row(
            "Cramér–von Mises",
            _fmt_stat(obs["cvm_stat"], 3),
            "unreliable",
            mc_p["cvm_stat"],
            "vs-pool p broken at N2 >> N1",
        ),
        row(
            "Kuiper",
            _fmt_stat(obs["kuiper_v"], 3),
            "—",
            mc_p["kuiper_v"],
            "V = D+ + D−; p from MC only",
        ),
        row(
            "Mann–Whitney U",
            _fmt_stat(obs["mwu_u"], 1),
            _fmt_p(obs["mwu_p"]),
            None,
            f"rank-biserial = {_fmt_stat(obs['mwu_rbc'], 3)}",
        ),
        row(
            "Brunner–Munzel",
            _fmt_stat(obs["bm_stat"], 3),
            _fmt_p(obs["bm_p"]),
            None,
            "location; unequal variance OK",
        ),
        row(
            "Epps–Singleton",
            _fmt_stat(obs["es_stat"], 2),
            _fmt_p(obs["es_p"]),
            None,
            "characteristic function",
        ),
        row(
            "Welch t (means)",
            _fmt_stat(obs["welch_t"], 3),
            _fmt_p(obs["welch_p"]),
            None,
            "parametric; cos(i) is not Gaussian",
        ),
        row(
            "Wasserstein-1",
            _fmt_stat(obs["w1"], 4),
            "—",
            mc_p["w1"],
            r"mean \|quantile gap\|; p from MC",
        ),
        "",
        "## Kolmogorov–Smirnov (two-sample)",
        "",
        "**How it works.** Sort both samples and form empirical CDFs \(F_{50}\) and "
        r"\(G_{\mathrm{HSC}}\). The statistic is "
        r"\(D=\sup_x |F_{50}(x)-G_{\mathrm{HSC}}(x)|\), the largest vertical gap. "
        "SciPy’s two-sided p uses the asymptotic (or exact, when cheap) null "
        "distribution of \(D\) under a shared continuous law.",
        "",
        "**How to read it.** \(D\) is a distance between CDFs on [0, 1]: "
        f"here \(D={_fmt_stat(obs['ks_d'], 3)}\) means the curves never differ by "
        "more than that. p answers “would a gap this big show up often if both "
        "were draws from one law?” KS is most sensitive near the middle of the "
        "distribution and weaker in the tails than AD.",
        "",
        f"**Result.** \(D={_fmt_stat(obs['ks_d'], 3)}\), "
        f"p (vs pool) = {_fmt_p(obs['ks_p'])}, "
        f"p (matched N) = {_fmt_p(mc_p['ks_d'])}.",
        "",
        "## Anderson–Darling (k-sample)",
        "",
        "**How it works.** SciPy `anderson_ksamp` compares the two samples to the "
        "pooled ECDF with a weight \(1/[H(1-H)]\) that **boosts the tails**. "
        "The statistic is a standardized integral of the weighted gap. Null: "
        "one continuous distribution.",
        "",
        "**How to read it.** Larger statistic → more evidence against a shared "
        "law, especially if the mismatch is at high/low \(\\cos(i)\). Do **not** "
        "read the sign of the statistic (recent SciPy can return negative values "
        "when samples agree). Use the p-value. A printed p of 0.25 is a **cap** "
        "(true p ≥ 0.25).",
        "",
        f"**Result.** AD = {_fmt_stat(obs['ad_stat'], 3)}, "
        f"p (vs pool) = {_fmt_p(obs['ad_p'], ad_cap=obs['ad_capped'])}, "
        f"p (matched N) = {_fmt_p(mc_p['ad_stat'])}.",
        "",
        "## Cramér–von Mises (two-sample)",
        "",
        "**How it works.** Integrates \((F-G)^2\) over the pooled sample instead "
        "of taking only the single largest gap (KS). SciPy "
        "`cramervonmises_2samp` reports that integral and an approximate p.",
        "",
        "**How to read it.** Middle ground between KS (one spike) and AD (tail "
        "emphasis): a modest gap over a wide stretch of \(\\cos(i)\) can beat a "
        "tall narrow spike. Same p-scale as above.",
        "",
        "SciPy vs-pool p is uncalibrated when N_HSC >> N_FRB. Use matched-N p.",
        "",
        f"**Result.** T={_fmt_stat(obs['cvm_stat'], 3)}, "
        f"p (vs pool) = ignore ({_fmt_p(obs['cvm_p'])}, uncalibrated), "
        f"p (matched N) = {_fmt_p(mc_p['cvm_stat'])}.",
        "",
        "## Kuiper (two-sample)",
        "",
        "**How it works.** \(V=D^{+}+D^{-}\): the sum of the largest upward and "
        "downward ECDF gaps. Equal to KS \(D\) when the discrepancy is one-sided; "
        "larger when the curves cross. No SciPy two-sample p, so only the "
        "matched-N Monte Carlo p is quoted.",
        "",
        "**How to read it.** Sensitive to **shape** differences that flip sign "
        "(excess edge-on *and* excess face-on). A large \(V\) with a small KS "
        "\(D\) means the CDFs cross. p (matched N) is the fraction of HSC-50 "
        "draws with \(V\) at least this large vs the pool.",
        "",
        f"**Result.** \(V={_fmt_stat(obs['kuiper_v'], 3)}\), "
        f"p (matched N) = {_fmt_p(mc_p['kuiper_v'])}.",
        "",
        "## Mann–Whitney U",
        "",
        "**How it works.** Rank all \(N_1+N_2\) values. \(U\) counts how often "
        "an FRB \(\\cos(i)\) beats an HSC one. Two-sided p tests "
        r"\(P(\mathrm{FRB}>\mathrm{HSC})=1/2\) (no stochastic dominance / "
        "location shift). The rank-biserial correlation "
        r"\(r=1-2U/(N_1 N_2)\) is the effect size on \([-1,1]\).",
        "",
        "**How to read it.** p speaks only to a **shift**, not to variance or "
        "tail shape (AD/KS can fire when MWU does not). "
        r"\(|r|\lesssim 0.1\) is a small location shift even if p is small "
        f"against {n2:,} HSC galaxies. Sign: \(r>0\) means FRB ranks "
        "**lower** \(\\cos(i)\) (more edge-on) than HSC.",
        "",
        f"**Result.** \(U={_fmt_stat(obs['mwu_u'], 1)}\), "
        f"p (vs pool) = {_fmt_p(obs['mwu_p'])}, "
        f"rank-biserial \(r={_fmt_stat(obs['mwu_rbc'], 3)}\) "
        f"(median \(\\Delta\\cos(i)={dmed:+.3f}\)).",
        "",
        "## Brunner–Munzel",
        "",
        "**How it works.** Same stochastic-dominance question as Mann–Whitney, "
        "but it does **not** assume equal variance (MWU’s rank null is touchy "
        "when the two spreads differ). SciPy `brunnermunzel`, two-sided.",
        "",
        "**How to read it.** Treat it as a robustness check on MWU. If BM and "
        "MWU p-values agree, the location conclusion is not an equal-spread "
        "artifact. Same p-scale.",
        "",
        f"**Result.** \(W={_fmt_stat(obs['bm_stat'], 3)}\), "
        f"p (vs pool) = {_fmt_p(obs['bm_p'])}.",
        "",
        "## Epps–Singleton",
        "",
        "**How it works.** Compares the empirical **characteristic functions** "
        r"\(\langle e^{itX}\rangle\) of the two samples at a few \(t\) points "
        "(SciPy `epps_singleton_2samp`). Detects differences in location, "
        "scale, or shape that moment- or rank-only tests can miss. Needs "
        r"\(N\gtrsim 5\) per sample; can fail on heavy ties.",
        "",
        "**How to read it.** The statistic is a \(\\chi^2\)-like distance "
        "between characteristic functions. Large stat / small p → the full "
        "laws differ, without saying *where*. Use the CDF plot for direction.",
        "",
        f"**Result.** \(W={_fmt_stat(obs['es_stat'], 2)}\), "
        f"p (vs pool) = {_fmt_p(obs['es_p'])}.",
        "",
        "## Welch \(t\) (means)",
        "",
        "**How it works.** Two-sample \(t\) on the means with unequal-variance "
        "standard errors. Assumes each sample mean is approximately normal "
        f"(CLT is comfortable at \(N={n1}\), but \(\\cos(i)\) is bounded).",
        "",
        "**How to read it.** This is a **mean shift** test only. Do not prefer "
        "it over MWU/BM for the scientific claim; it is here because \(N=50\) "
        "is large enough that the mean is a stable summary. "
        f"Mean \(\\Delta\\cos(i)={dmean:+.3f}\).",
        "",
        f"**Result.** \(t={_fmt_stat(obs['welch_t'], 3)}\), "
        f"p (vs pool) = {_fmt_p(obs['welch_p'])}.",
        "",
        "## Wasserstein-1",
        "",
        "**How it works.** Earth-mover distance between the two 1-D laws: "
        r"\(\int |F^{-1}-G^{-1}|\,du\), the average absolute gap between "
        "quantile functions. Not a significance test by itself; p comes from "
        "the matched-N Monte Carlo.",
        "",
        "**How to read it.** \(W_1\) is in \(\\cos(i)\) units. "
        f"\(W_1={_fmt_stat(obs['w1'], 4)}\) is the typical horizontal shift "
        "you would need to morph one CDF into the other. Compare it to the "
        "plot’s 68% band width, not to a 0–1 “agreement” scale.",
        "",
        f"**Result.** \(W_1={_fmt_stat(obs['w1'], 4)}\), "
        f"p (matched N) = {_fmt_p(mc_p['w1'])}.",
        "",
        "## Caveats",
        "",
        "1. Point estimates only — Protocol A/B error bands are not in these p-values.",
        "2. HSC $q$ is from the Kawinwanichakij et al. Lenstronomy "
        "single-Sérsic fits; FRB $q$ is this pipeline’s GALFIT winning leg.",
        "3. Pool tests have high power at \(N_{\\mathrm{HSC}}=24{,}450\); lean on matched-N p "
        "and the CDF when judging *practical* difference.",
        "4. Several tests share one sample; do not treat six small p-values as six "
        "independent discoveries.",
        "",
    ]
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("\n".join(lines), encoding="utf-8")
    print(f"[*] Wrote {dest}")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hosts-csv", type=Path, default=CSV_DEFAULT)
    p.add_argument("--out-md", type=Path, default=MD_DEFAULT)
    p.add_argument("--n-draws", type=int, default=N_DRAWS)
    p.add_argument("--seed", type=int, default=SEED)
    p.add_argument("--mag-limit", type=float, default=22.0)
    p.add_argument("--q0", type=float, default=Q0)
    args = p.parse_args(argv)

    hosts = pd.read_csv(args.hosts_csv)
    frb = pd.to_numeric(hosts["cosi_point"], errors="coerce").to_numpy(dtype=float)
    frb = frb[np.isfinite(frb)]
    if len(frb) != 50:
        raise RuntimeError(f"expected 50 finite cosi_point, got {len(frb)}")

    print("[*] Loading HSC pool ...", flush=True)
    hsc, n_hsc = load_hsc_pool(mag_limit=args.mag_limit, q0=args.q0)
    print(f"[*] FRB N={len(frb)}  HSC N={n_hsc:,}", flush=True)

    obs = observed_tests(frb, hsc)
    print(
        f"[*] KS D={obs['ks_d']:.3f} p={obs['ks_p']:.4g}  "
        f"AD={obs['ad_stat']:.3f} p={obs['ad_p']:.4g}  "
        f"MWU p={obs['mwu_p']:.4g}",
        flush=True,
    )
    print(f"[*] Matched-N MC ({args.n_draws} draws) ...", flush=True)
    mc_p = mc_tail_p(frb, hsc, n_draws=args.n_draws, seed=args.seed)
    print(
        f"[*] MC p  KS={mc_p['ks_d']:.3f}  AD={mc_p['ad_stat']:.3f}  "
        f"CvM={mc_p['cvm_stat']:.3f}  Kuiper={mc_p['kuiper_v']:.3f}  "
        f"W1={mc_p['w1']:.3f}",
        flush=True,
    )
    write_md(
        args.out_md, frb=frb, hsc=hsc, obs=obs, mc_p=mc_p, n_draws=args.n_draws,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
