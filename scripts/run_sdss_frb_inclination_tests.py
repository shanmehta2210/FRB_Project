#!/usr/bin/env python3
"""
Anderson–Darling and Mann–Whitney U tests: SDSS null vs FRB host inclination.

Scientific question
-------------------
Under the same catalog cuts used for mag-sliced inclination CDF plots, do FRB
host inclinations look like a random draw from the SDSS field-galaxy inclination
distribution, or is there evidence they differ?

This script does **not** refit galaxies or rebuild CDFs. It reuses the exact
selection logic from ``plot_null_mag_cut_cdfs.py`` and compares two 1D samples
per scenario with standard two-sample tests in ``scipy.stats``.

Comparison modes (default mag limits 20, 21, 22)
------------------------------------------------
**Mode A — strict cos(i)** [matches CDF x-axis]
  - SDSS: ``prepare_null_strict_color_base`` + ``slice_null_base_by_mag``
    (``u-r < 2.3``, ``best_model_ba_r > 0.2``, ``modelMag_r <= mag``).
  - FRB: ``frb_hosts_for_cdf(..., sample_mode="strict")`` — GALFIT ``mag`` cut,
    ``b/a > 0.2``; **no** host color cut.
  - SDSS inclination: Hubble ``cos(i)`` from ``best_model_ba_r`` (``q0=0.2``).
  - FRB inclination: ``cos(radians(inc))`` from GALFIT **point** ``inc``.

**Mode B — strict i (degrees)**
  - Same rows as Mode A.
  - SDSS: ``degrees(arccos(cos i))`` from the Hubble cos(i) array.
  - FRB: GALFIT ``inc`` column (degrees).
  - Tests are rank-equivalent to Mode A (monotonic transform), so statistics
    should match Mode A closely.

**Mode C — inclusive cos(i)** [relaxed axis-ratio cut]
  - SDSS: ``prepare_null_inclusive_color_base`` + mag slice — ``u-r < 2.3``,
    finite ``b/a in [0, 1]``, **no** ``b/a > 0.2`` requirement.
  - FRB: mag cut only (``sample_mode="inclusive"``).
  - Hubble mapping still used; when ``b/a <= q0`` the formula fails and
    ``cos(i)`` is clipped to 0 (face-on pile-up), as in inclusive null semantics.

Tests applied (per mode × mag limit)
------------------------------------
1. **Anderson–Darling k-sample** — ``scipy.stats.anderson_ksamp([sdss, frb])``
   - Null hypothesis: both samples come from the same continuous distribution.
   - Sensitive to differences in **shape** (tails, skew), not only the median.
   - Reports ``statistic`` and ``pvalue``. In SciPy ≥1.15, ``pvalue`` can be
     **capped at 0.25** when samples are very similar (not a bug).
   - Negative AD statistics can appear with the updated SciPy implementation;
     interpret using **p-value**, not the sign of the statistic alone.

2. **Mann–Whitney U** — ``scipy.stats.mannwhitneyu(..., alternative="two-sided")``
   - Null hypothesis: P(FRB value > SDSS value) = 0.5 (equal rank / location shift).
   - Rank-based; robust to outliers compared with a t-test.
   - ``U`` is **not** on a fixed 0–1 scale when ``N_SDSS >> N_FRB``; use **p-value**.
   - With ~10⁴–10⁵ SDSS galaxies vs ~30–60 FRBs, the test has **very high power**:
     small systematic shifts can yield small p-values even when CDFs look similar
     by eye.

Important differences vs CDF plots
------------------------------------
- CDF figures use **MC draws** for FRB (``inc`` ± ``inc_err``) and bootstrap
  subsampling of the null pool to size ``N_FRB``. These tests use:
  - full SDSS pool at each mag (not subsampled to ``N_FRB``), and
  - FRB **point estimates** only (no measurement-error propagation).
- SDSS inclinations come from catalog **b/a**; FRB inclinations from **GALFIT**.
  Tests compare the numbers used in plots, but pipelines differ by construction.

Inputs (defaults)
-------------------
- ``SDSS_catalog_v1_allsky_modelmr.csv`` — same as CDF driver.
- ``pipeline_galfit_results.csv`` — hosts with finite ``inc``.
- Constants: ``SDSS_UR_MAX_CDF=2.3``, ``q0=0.2``, columns ``modelMag_r``,
  ``best_model_ba_r``.

Output
------
- ``test_results.md`` at repo root (tables + interpretation guide).

Run from repo root::

    python scripts/run_sdss_frb_inclination_tests.py
    python scripts/run_sdss_frb_inclination_tests.py --mag-limits 20 21 22 --out-md test_results.md
"""

from __future__ import annotations

import argparse
import math
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (  # noqa: E402
    Q0,
    SDSS_UR_MAX_CDF,
    cosi_array_from_df,
    inc_deg_from_cosi,
    prepare_null_inclusive_color_base,
    prepare_null_strict_color_base,
    read_sdss_null_catalog,
    slice_null_base_by_mag,
)
from pipeline_null_plot_utils import (  # noqa: E402
    DEFAULT_PIPELINE,
    DEFAULT_SDSS,
    frb_hosts_for_cdf,
    load_pipeline_hosts,
)

MAG_LIMITS_DEFAULT = [20, 21, 22]
OUT_MD_DEFAULT = Path(__file__).resolve().parents[1] / "test_results.md"

# Interpretation block appended to test_results.md (kept in this file as source of truth).
INTERPRETATION_MD = """
## What these tests are asking

Each row tests whether the **FRB host inclination sample** could have been drawn from
the same distribution as the **SDSS field-galaxy null** sample, after applying the
listed magnitude and color (and, for Modes A/B, axis-ratio) cuts. A small p-value
means the two samples differ in a way that is unlikely under the test's null
hypothesis—not automatically that FRBs are "more edge-on" (check CDFs for direction).

These tests complement the CDF figures in
`plots/plots_null/v1_null_cdf_inclination/mag_cuts/`; they do not replace visual
inspection.

## How to read p-values (the scale)

**Direction:** Think of p as “how surprising would these data be if FRB and SDSS
inclinations really came from the same distribution?”

- **High p (close to 1)** → *not surprising* → data are **compatible with “similar”**
- **Low p (close to 0)** → *very surprising* → data **favour “different”**

**It is not a similarity percentage.** MWU p = 0.28 does **not** mean “28% different”
or “72% the same.” It is not on a 0–100% “how alike are the CDFs?” scale.

**Usual cutoffs (convention, not physics):**

| p-value | Plain-language read |
|---------|---------------------|
| **> 0.10** | No meaningful evidence of a difference (for this test). Treat as **consistent with the null**. |
| **0.05 – 0.10** | Weak / suggestive only; many fields still call this **not significant**. |
| **0.01 – 0.05** | **Significant** at the common 5% bar — difference is unlikely to be pure chance. |
| **< 0.01** | Stronger significance (1% bar) — still says nothing about *how big* the shift is. |

**Your numbers, on that scale:**

| MWU p | What it means |
|-------|----------------|
| **0.98** (mag < 20) | **Remarkably consistent with “no difference.”** If the true distributions were the same, you would often see a p-value this high. This is as “they look alike” as these tests get. |
| **0.28** (mag < 21) | **Not significantly different** at the usual 5% level. Nowhere near “really, really different”—it is the opposite: the data are **plausibly from the same distribution**, with no strong rank shift detected. |
| **0.009** (mag < 22) | **Significant** (~1% level): a systematic shift in inclination is unlikely to be chance alone. That is **real statistical evidence of a difference**, but it does **not** by itself mean the CDFs are wildly separated—you still check the plot for *size* and *direction*. |

**Effect size vs significance:** With ~100k SDSS galaxies, even a **small** CDF offset can
yield p ≈ 0.01. Significance answers “is there a detectable shift?” not “is it huge?”
Use the mag 22 CDF overlay for “how much” and “which way.”

**AD p = 0.25 (capped):** SciPy’s way of saying p is **at least 0.25**—even more
“everything looks fine / similar” than 0.28. Not a separate scale; same rule: high = similar.

## Column guide

| Column | Meaning |
|--------|---------|
| **mag limit** | Keep galaxies with `modelMag_r` (SDSS) or GALFIT `mag` (FRB) ≤ this value. |
| **N_SDSS** | Size of the SDSS null pool after all cuts for that mode. Should match `mag_cut_summary.csv` for strict modes. |
| **N_FRB** | Number of FRB hosts passing the matching FRB cuts. |
| **AD statistic** | Anderson–Darling k-sample statistic (SciPy `anderson_ksamp`). Larger values → more evidence against "same distribution." Can be negative in recent SciPy versions; rely on **AD p**. |
| **AD p** | Approximate p-value for equal distributions. Values reported as **0.25** are **capped** (SciPy: true p > 0.25; samples very similar). |
| **AD note** | Short readout of AD p (capped / significant). |
| **MWU U** | Mann–Whitney U statistic. Not normalized when sample sizes differ by orders of magnitude; use **MWU p** for inference. |
| **MWU p** | Two-sided p-value for equal rank distributions. Small p → systematic shift in inclination between samples. |

## How to read the three modes

**Mode A (strict cos i)** — Closest to the CDF plots (x-axis is cos i). Use this mode
when comparing test outcomes to `mag_cuts/magXX/sdss_strict/null_cdf_inclination.png`.

**Mode B (strict i deg)** — Same galaxies as Mode A, but inclination in degrees.
Because cos i is a monotonic function of i on [0°, 90°], **AD and MWU results
should match Mode A** (only tiny floating-point differences possible).

**Mode C (inclusive cos i)** — Drops the `b/a > 0.2` cut on both surveys (FRB: mag only).
The SDSS pool is larger and includes face-on systems piled up at cos i = 0 when Hubble's
formula fails. Use this to see whether conclusions depend on excluding very round galaxies.

## Interpreting your current results (strict modes A/B)

| mag | Rough takeaway |
|-----|----------------|
| **< 20** | No evidence for a different distribution (AD p capped at 0.25; MWU p ≈ 0.98). Bright, small FRB subsample (N=31). |
| **< 21** | Still consistent with the null (AD p capped; MWU p ≈ 0.28). |
| **< 22** | **Significant difference** at ~1% level (AD p ≈ 0.008; MWU p ≈ 0.009). Visually, check whether FRB CDF sits above/below the SDSS band in the mag22 plot—tests do not state the direction of the shift. |

Mode C at mag < 22 shows an even smaller p-value because the inclusive null includes
more face-on galaxies (cos i → 0), which can exaggerate FRB–null separation if FRBs
are less face-on on average.

## Limitations and caveats

1. **Sample size imbalance** — ~10⁴–10⁵ SDSS vs ~30–60 FRBs. Mann–Whitney has high
   power; mag 22 significance may reflect a modest visual offset, not a dramatic effect.
2. **No FRB measurement errors** — CDFs perturb `inc` using `inc_err`; these tests use
   point estimates only.
3. **Different measurement systems** — SDSS: Hubble formula on catalog b/a. FRB: GALFIT
   Sérsic fit inclinations. Systematic offsets can masquerade as distribution differences.
4. **Selection not identical** — SDSS null is color-trimmed (late-type proxy); FRB hosts
   are not color-cut. Matches the CDF methodology but is not a perfect physical match.
5. **Multiple comparisons** — Nine table rows (3 modes × 3 mag limits) without Bonferroni
   or similar correction; treat borderline cases cautiously.
6. **Mag limit is a hard cut** — Pools at mag 20, 21, 22 are nested (larger limits include
   fainter galaxies). Results are correlated across rows, not independent experiments.
7. **Not causal** — Rejecting "same distribution" does not identify astrophysics (e.g.
   host type, redshift, selection in FRB surveys) without further modeling.

## Method reference

Implementation: `scripts/run_sdss_frb_inclination_tests.py`. Shared cuts:
`scripts/null_catalog_utils.py`, `scripts/pipeline_null_plot_utils.py` (`frb_hosts_for_cdf`).
CDF driver: `scripts/plot_null_mag_cut_cdfs.py`.
"""


@dataclass
class TestRow:
    """One row of the results table (one mode × one mag limit)."""

    mode: str
    mag_limit: float
    variable: str
    n_sdss: int
    n_frb: int
    ad_statistic: float
    ad_pvalue: float
    ad_note: str
    mwu_u: float
    mwu_pvalue: float


def ad_pvalue_and_note(ad_result) -> tuple[float, str]:
    """
    Extract Anderson–Darling p-value and a short human-readable note.

    SciPy caps reported p-values at 0.25 when the true p-value is larger
    (samples are similar under the AD null).
    """
    p = float(ad_result.pvalue)
    if p >= 0.25 - 1e-12:
        return p, "p capped at 0.25 (not significant)"
    if p < 0.001:
        return p, "p < 0.001"
    return p, f"p = {p:.4f}"


def run_two_sample_tests(
    sdss_vals: np.ndarray, frb_vals: np.ndarray
) -> tuple[float, float, str, float, float]:
    """
    Run Anderson–Darling and Mann–Whitney U on two 1D samples.

    Parameters
    ----------
    sdss_vals, frb_vals
        Inclination-related scalars (cos i or i in degrees). Non-finite values
        are dropped before testing.

    Returns
    -------
    ad_statistic, ad_pvalue, ad_note, mwu_u, mwu_pvalue
    """
    sdss = np.asarray(sdss_vals, dtype=float)
    sdss = sdss[np.isfinite(sdss)]
    frb = np.asarray(frb_vals, dtype=float)
    frb = frb[np.isfinite(frb)]
    if len(sdss) == 0 or len(frb) == 0:
        raise RuntimeError("Empty sample after finite filter.")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        ad = stats.anderson_ksamp([sdss, frb])
    ad_p, ad_note = ad_pvalue_and_note(ad)
    mwu = stats.mannwhitneyu(sdss, frb, alternative="two-sided")
    return float(ad.statistic), ad_p, ad_note, float(mwu.statistic), float(mwu.pvalue)


def frb_cosi_from_inc(hosts: pd.DataFrame) -> np.ndarray:
    """FRB cos(i) from GALFIT inclination point estimates (degrees → radians)."""
    inc = pd.to_numeric(hosts["inc"], errors="coerce").to_numpy(dtype=float)
    inc = np.clip(inc, 0.0, 90.0)
    return np.cos(np.radians(inc))


def frb_inc_deg(hosts: pd.DataFrame) -> np.ndarray:
    """Finite FRB inclinations in degrees (GALFIT ``inc`` column)."""
    inc = pd.to_numeric(hosts["inc"], errors="coerce").to_numpy(dtype=float)
    return inc[np.isfinite(inc)]


def build_sdss_bases(
    sdss_raw: pd.DataFrame,
    *,
    sdss_mag_column: str,
    sdss_q_column: str,
    q0: float,
    sdss_ur_max: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build mag-agnostic SDSS pools for strict and inclusive tests.

    Mirrors the one-time catalog load in ``plot_null_mag_cut_cdfs.py``; magnitude
    limits are applied later via ``slice_null_base_by_mag``.
    """
    strict = prepare_null_strict_color_base(
        sdss_raw,
        mag_column=sdss_mag_column,
        q0=q0,
        q_column=sdss_q_column,
        is_legacy=False,
        sdss_ur_max=sdss_ur_max,
    )
    inclusive = prepare_null_inclusive_color_base(
        sdss_raw,
        mag_column=sdss_mag_column,
        q0=q0,
        q_column=sdss_q_column,
        is_legacy=False,
        sdss_ur_max=sdss_ur_max,
    )
    return strict, inclusive


def collect_rows(
    *,
    hosts: pd.DataFrame,
    sdss_strict_base: pd.DataFrame,
    sdss_inclusive_base: pd.DataFrame,
    mag_limits: list[float],
    q0: float,
    sdss_mag_column: str,
    sdss_q_column: str,
) -> list[TestRow]:
    """
    For each mag limit and comparison mode, build samples and run both tests.

    Sample construction is intentionally parallel to the CDF overlay logic in
    ``plot_null_mag_cut_cdfs.plot_one`` (same ``frb_hosts_for_cdf`` and
    ``cosi_array_from_df``), but without MC or null subsampling.
    """
    rows: list[TestRow] = []

    for mag_limit in sorted(mag_limits):
        sdss_strict = slice_null_base_by_mag(
            sdss_strict_base,
            mag_column=sdss_mag_column,
            mag_limit=mag_limit,
        )
        sdss_inclusive = slice_null_base_by_mag(
            sdss_inclusive_base,
            mag_column=sdss_mag_column,
            mag_limit=mag_limit,
        )
        frb_strict = frb_hosts_for_cdf(
            hosts,
            sample_mode="strict",
            q0=q0,
            mag_limit=mag_limit,
        )
        frb_inclusive = frb_hosts_for_cdf(
            hosts,
            sample_mode="inclusive",
            q0=q0,
            mag_limit=mag_limit,
        )

        # SDSS cos(i): Hubble formula on best_model_ba_r (see hubble_cosi_from_ba).
        cosi_sdss_strict = cosi_array_from_df(sdss_strict, q_col=sdss_q_column, q0=q0)
        cosi_frb_strict = frb_cosi_from_inc(frb_strict)
        inc_sdss_strict = inc_deg_from_cosi(cosi_sdss_strict)
        inc_frb_strict = frb_inc_deg(frb_strict)

        cosi_sdss_incl = cosi_array_from_df(sdss_inclusive, q_col=sdss_q_column, q0=q0)
        cosi_frb_incl = frb_cosi_from_inc(frb_inclusive)

        scenarios = [
            (
                "A_strict_cosi",
                "cos(i)",
                cosi_sdss_strict,
                cosi_frb_strict,
                len(sdss_strict),
                len(frb_strict),
            ),
            (
                "B_strict_i_deg",
                "i (deg)",
                inc_sdss_strict,
                inc_frb_strict,
                len(sdss_strict),
                len(frb_strict),
            ),
            (
                "C_inclusive_cosi",
                "cos(i)",
                cosi_sdss_incl,
                cosi_frb_incl,
                len(sdss_inclusive),
                len(frb_inclusive),
            ),
        ]

        for mode, variable, sdss_arr, frb_arr, n_sdss, n_frb in scenarios:
            ad_stat, ad_p, ad_note, mwu_u, mwu_p = run_two_sample_tests(sdss_arr, frb_arr)
            rows.append(
                TestRow(
                    mode=mode,
                    mag_limit=mag_limit,
                    variable=variable,
                    n_sdss=n_sdss,
                    n_frb=n_frb,
                    ad_statistic=ad_stat,
                    ad_pvalue=ad_p,
                    ad_note=ad_note,
                    mwu_u=mwu_u,
                    mwu_pvalue=mwu_p,
                )
            )
            print(
                f"[*] mag<{mag_limit:g} {mode}: N_sdss={n_sdss}, N_frb={n_frb}, "
                f"AD={ad_stat:.4g} ({ad_note}), MWU p={mwu_p:.4g}"
            )

    return rows


def format_pvalue(p: float) -> str:
    """Format p-values for markdown tables."""
    if not math.isfinite(p):
        return "nan"
    if p < 1e-4:
        return f"{p:.2e}"
    return f"{p:.4f}"


def write_results_md(rows: list[TestRow], out_path: Path) -> None:
    """Write results tables plus the interpretation guide to markdown."""
    lines = [
        "# SDSS vs FRB inclination statistical tests",
        "",
        "Anderson–Darling (`scipy.stats.anderson_ksamp`) and Mann–Whitney U "
        "(two-sided). SDSS: `u-r < 2.3`. FRB: GALFIT point estimates (no MC). "
        f"Hubble `q0={Q0}`. Driver: `scripts/run_sdss_frb_inclination_tests.py`.",
        "",
        "## Mode A — strict cos(i)",
        "",
        "SDSS: `modelMag_r` cut, `best_model_ba_r > 0.2`, Hubble cos(i). "
        "FRB: mag cut, `b/a > 0.2`, cos(GALFIT inc). Same cuts as CDF plots.",
        "",
        "| mag limit | N_SDSS | N_FRB | AD statistic | AD p | AD note | MWU U | MWU p |",
        "|-----------|--------|-------|--------------|--------|---------|-------|-------|",
    ]

    for r in sorted([x for x in rows if x.mode == "A_strict_cosi"], key=lambda x: x.mag_limit):
        lines.append(row_line(r))

    lines.extend(
        [
            "",
            "## Mode B — strict i (deg)",
            "",
            "Same sample selection as Mode A; SDSS i from arccos(cos i), FRB GALFIT `inc`.",
            "",
            "| mag limit | N_SDSS | N_FRB | AD statistic | AD p | AD note | MWU U | MWU p |",
            "|-----------|--------|-------|--------------|--------|---------|-------|-------|",
        ]
    )
    for r in sorted([x for x in rows if x.mode == "B_strict_i_deg"], key=lambda x: x.mag_limit):
        lines.append(row_line(r))

    lines.extend(
        [
            "",
            "## Mode C — inclusive cos(i)",
            "",
            "SDSS: mag + color, finite b/a in [0,1] (no b/a>0.2). "
            "FRB: mag cut only. Hubble cos(i) with cos i=0 when b/a ≤ q0.",
            "",
            "| mag limit | N_SDSS | N_FRB | AD statistic | AD p | AD note | MWU U | MWU p |",
            "|-----------|--------|-------|--------------|--------|---------|-------|-------|",
        ]
    )
    for r in sorted([x for x in rows if x.mode == "C_inclusive_cosi"], key=lambda x: x.mag_limit):
        lines.append(row_line(r))

    lines.append(INTERPRETATION_MD.strip())
    lines.append("")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def row_line(r: TestRow) -> str:
    return (
        f"| {r.mag_limit:g} | {r.n_sdss} | {r.n_frb} | {r.ad_statistic:.4g} | "
        f"{format_pvalue(r.ad_pvalue)} | {r.ad_note} | {r.mwu_u:.4g} | "
        f"{format_pvalue(r.mwu_pvalue)} |"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="SDSS vs FRB inclination two-sample tests (see module docstring).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Modes: A=strict cos(i) [CDF-matched], B=strict i(deg), C=inclusive cos(i). "
            "Output includes tables and an interpretation guide."
        ),
    )
    parser.add_argument(
        "--pipeline-csv",
        type=Path,
        default=DEFAULT_PIPELINE,
        help="GALFIT pipeline results (FRB hosts).",
    )
    parser.add_argument(
        "--sdss-csv",
        type=Path,
        default=DEFAULT_SDSS,
        help="SDSS v1 null catalog CSV.",
    )
    parser.add_argument(
        "--mag-limits",
        type=float,
        nargs="+",
        default=MAG_LIMITS_DEFAULT,
        help="modelMag_r / GALFIT mag upper limits (default: 20 21 22).",
    )
    parser.add_argument(
        "--q0",
        type=float,
        default=Q0,
        help="Intrinsic axis ratio in Hubble formula (default 0.2).",
    )
    parser.add_argument(
        "--sdss-mag-column",
        default="modelMag_r",
        help="SDSS magnitude column for mag cuts.",
    )
    parser.add_argument(
        "--sdss-q-column",
        default="expAB_r",
        help="SDSS axis-ratio column for Hubble cos(i).",
    )
    parser.add_argument(
        "--sdss-ur-max",
        type=float,
        default=SDSS_UR_MAX_CDF,
        help="SDSS u-r color cut (default 2.3).",
    )
    parser.add_argument(
        "--out-md",
        type=Path,
        default=OUT_MD_DEFAULT,
        help="Markdown output path (default: repo root test_results.md).",
    )
    args = parser.parse_args()

    print("[*] Loading SDSS null catalog...")
    sdss_raw = read_sdss_null_catalog(args.sdss_csv)
    sdss_strict_base, sdss_inclusive_base = build_sdss_bases(
        sdss_raw,
        sdss_mag_column=args.sdss_mag_column,
        sdss_q_column=args.sdss_q_column,
        q0=args.q0,
        sdss_ur_max=args.sdss_ur_max,
    )
    del sdss_raw
    print(
        f"[*] SDSS bases: strict N={len(sdss_strict_base)}, "
        f"inclusive N={len(sdss_inclusive_base)}"
    )

    hosts = load_pipeline_hosts(args.pipeline_csv)
    rows = collect_rows(
        hosts=hosts,
        sdss_strict_base=sdss_strict_base,
        sdss_inclusive_base=sdss_inclusive_base,
        mag_limits=args.mag_limits,
        q0=args.q0,
        sdss_mag_column=args.sdss_mag_column,
        sdss_q_column=args.sdss_q_column,
    )

    write_results_md(rows, args.out_md)
    print(f"[*] Wrote {len(rows)} test rows + interpretation to {args.out_md}")


if __name__ == "__main__":
    main()
