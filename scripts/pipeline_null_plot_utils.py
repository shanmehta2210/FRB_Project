"""
Shared utilities for pipeline null / diagnostic plots.
"""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate

from null_catalog_utils import (
    LEGACY_GR_MAX_CDF,
    LS_CATALOG_V1_DEFAULT,
    Q0,
    REPO_ROOT,
    SDSS_CATALOG_V1_DEFAULT,
    SDSS_CATALOG_V2_DEFAULT,
    SDSS_UR_MAX_CDF,
    assign_values_to_bin_edges,
    cosi_array_from_df,
    ensure_sdss_n_eff,
    equal_count_quantile_edges,
    filter_frb_hosts_mag,
    filter_frb_hosts_strict_ba,
    inc_deg_from_cosi,
    n_from_legacy_df,
    n_from_sdss_df,
    prepare_null_sample,
    read_legacy_null_catalog,
    read_sdss_null_catalog,
    re_arcsec_from_legacy_df,
    re_arcsec_from_sdss_df,
)

DEFAULT_PIPELINE = REPO_ROOT / "pipeline_galfit_results.csv"
DEFAULT_LOC = REPO_ROOT / "master_frb_localization.csv"
DEFAULT_SDSS = REPO_ROOT / SDSS_CATALOG_V1_DEFAULT
DEFAULT_SDSS_V2 = REPO_ROOT / SDSS_CATALOG_V2_DEFAULT
DEFAULT_LEGACY = REPO_ROOT / LS_CATALOG_V1_DEFAULT

PLATE_SCALE_DEFAULT = 0.262
PLATE_SCALE_BY_FRB: dict[str, float] = {
    "20171020A": 0.25,
    "20210807D": 0.25,
    "20211127I": 0.25,
    "20220207C": 0.25,
    "20220307B": 0.25,
    "20220319D": 0.25,
    "20220825A": 0.25,
    "20220912A": 0.25,
}

PLOTS_NULL = REPO_ROOT / "plots" / "plots_null"

np.random.seed(42)


def plate_scale_arcsec(frb: str) -> float:
    return PLATE_SCALE_BY_FRB.get(str(frb), PLATE_SCALE_DEFAULT)


def save_figure(fig: plt.Figure, out_stem: Path) -> None:
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    path = out_stem.with_suffix(".png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    print(f"Saved {path}")


def cdf_envelope(
    reference_vals: np.ndarray,
    n_sample: int,
    n_draws: int = 10000,
    rng_seed: int = 42,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(rng_seed)
    ref = np.asarray(reference_vals, dtype=float)
    ref = ref[np.isfinite(ref)]
    if len(ref) < n_sample:
        raise RuntimeError(f"Null pool too small: {len(ref)} < {n_sample}")

    total_samples = []
    for _ in range(n_draws):
        draw = np.sort(rng.choice(ref, size=n_sample, replace=False))
        idx_norm = [0.0] + [j / n_sample for j in range(1, n_sample + 1)] + [1.0]
        draw_ext = [0.0] + draw.tolist() + [1.0]
        total_samples.append(interpolate.interp1d(draw_ext, idx_norm))

    x = np.linspace(0, 1, 100)
    means, down, up = [], [], []
    lo = int(0.16 * n_draws)
    hi = int(0.84 * n_draws)
    for value in x:
        idx_sample = sorted(float(s(value)) for s in total_samples)
        means.append(float(np.mean(idx_sample)))
        down.append(float(idx_sample[lo]))
        up.append(float(idx_sample[hi]))
    return x, np.array(means), np.array(down), np.array(up)


def _x_grid_for_values(
    *arrays: np.ndarray,
    n_pts: int = 100,
    lo_percentile: float = 0.5,
    hi_pad_frac: float = 0.04,
) -> np.ndarray:
    """
    X grid for scalar CDF plots.

    Lower edge: low percentile (avoid log-like clutter near zero).
    Upper edge: **max** of all pooled values (FRB + null), not p99 — so a
    sparse null tail cannot truncate the host CDF when hosts have larger Re.
    """
    pooled = np.concatenate([a[np.isfinite(a)] for a in arrays if len(a)])
    if len(pooled) == 0:
        raise RuntimeError("No finite values for CDF grid.")
    lo = max(0.0, float(np.percentile(pooled, lo_percentile)))
    hi = float(np.max(pooled)) * (1.0 + hi_pad_frac)
    if hi <= lo:
        hi = lo + 1.0
    return np.linspace(lo, hi, n_pts)


def ecdf_on_grid(sorted_vals: np.ndarray, x_grid: np.ndarray) -> np.ndarray:
    n = len(sorted_vals)
    return np.searchsorted(sorted_vals, x_grid, side="right") / n


def cdf_envelope_scalar(
    reference_vals: np.ndarray,
    n_sample: int,
    x_grid: np.ndarray,
    n_draws: int = 10000,
    rng_seed: int = 42,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(rng_seed)
    ref = np.asarray(reference_vals, dtype=float)
    ref = ref[np.isfinite(ref)]
    if len(ref) < n_sample:
        raise RuntimeError(f"Null pool too small: {len(ref)} < {n_sample}")

    curves = []
    for _ in range(n_draws):
        draw = np.sort(rng.choice(ref, size=n_sample, replace=False))
        curves.append(ecdf_on_grid(draw, x_grid))

    arr = np.array(curves)
    lo_i = int(0.16 * n_draws)
    hi_i = int(0.84 * n_draws)
    mean = arr.mean(axis=0)
    down = np.sort(arr, axis=0)[lo_i]
    up = np.sort(arr, axis=0)[hi_i]
    return x_grid, mean, down, up


def mc_mean_cdf_from_draws(draws: list[list[float]]) -> tuple[np.ndarray, np.ndarray]:
    n_sample = len(draws[0])
    funcs = []
    for draw in draws:
        idx_norm = [0.0] + [j / n_sample for j in range(1, n_sample + 1)] + [1.0]
        draw_ext = [0.0] + list(draw) + [1.0]
        funcs.append(interpolate.interp1d(draw_ext, idx_norm))
    x = np.linspace(0, 1, 100)
    mean = [float(np.mean([float(f(value)) for f in funcs])) for value in x]
    return x, np.array(mean)


def build_frb_mc_draws_inc(
    hosts: pd.DataFrame,
    *,
    n_draws: int = 500,
    inc_col: str = "inc",
    inc_err_col: str = "inc_err",
) -> list[list[float]]:
    inc_vals: list[tuple[float, float]] = []
    for _, row in hosts.iterrows():
        mu = pd.to_numeric(row.get(inc_col), errors="coerce")
        sigma = pd.to_numeric(row.get(inc_err_col), errors="coerce")
        if not np.isfinite(mu):
            continue
        if not np.isfinite(sigma) or sigma < 0:
            sigma = 0.0
        inc_vals.append((min(90.0, max(0.0, float(mu))), float(sigma)))
    if not inc_vals:
        raise RuntimeError("No finite host inclinations.")
    draws: list[list[float]] = []
    for _ in range(n_draws):
        cosi = []
        for mu, sigma in inc_vals:
            sampled_inc = float(np.random.normal(mu, sigma)) if sigma > 0 else mu
            sampled_inc = min(90.0, max(0.0, sampled_inc))
            cosi.append(math.cos(math.radians(sampled_inc)))
        draws.append(sorted(cosi))
    return draws


def load_pipeline_hosts(path: Path = DEFAULT_PIPELINE) -> pd.DataFrame:
    df = pd.read_csv(path)
    inc = pd.to_numeric(df["inc"], errors="coerce")
    hosts = df.loc[inc.notna()].copy()
    if hosts.empty:
        raise RuntimeError(f"No hosts with finite inc in {path}")
    return hosts


def frb_hosts_for_cdf(
    hosts: pd.DataFrame,
    *,
    sample_mode: str,
    q0: float = Q0,
    ba_col: str = "b_a",
    mag_limit: float | None = None,
    frb_mag_column: str = "mag",
) -> pd.DataFrame:
    """
    FRB host rows used for inclination CDF overlays.

    If ``mag_limit`` is set, keep only GALFIT ``mag`` <= limit (matches null mag cut).
    strict: GALFIT b/a > q0 only (same rule as null strict pool).
    inclusive: all hosts passing any mag cut.
    """
    out = hosts
    if mag_limit is not None:
        out = filter_frb_hosts_mag(out, mag_limit=mag_limit, mag_column=frb_mag_column)
        if out.empty:
            raise RuntimeError(
                f"No FRB hosts with {frb_mag_column} <= {mag_limit}; cannot build CDF overlay."
            )
    if sample_mode == "strict":
        out = filter_frb_hosts_strict_ba(out, q0=q0, ba_col=ba_col)
        if out.empty:
            raise RuntimeError(
                f"No FRB hosts with {ba_col} > {q0} (strict); cannot build CDF overlay."
            )
        return out
    if sample_mode == "inclusive":
        return out
    raise ValueError(f"Unknown sample_mode: {sample_mode!r}")


def load_hosts_with_coords(
    pipeline_csv: Path = DEFAULT_PIPELINE,
    loc_csv: Path = DEFAULT_LOC,
) -> pd.DataFrame:
    hosts = load_pipeline_hosts(pipeline_csv)
    loc = pd.read_csv(loc_csv)
    merged = hosts.merge(
        loc[["frb", "ra_deg", "dec_deg", "survey"]],
        on="frb",
        how="left",
    )
    return merged


def frb_re_arcsec(hosts: pd.DataFrame) -> np.ndarray:
    re_px = pd.to_numeric(hosts["re"], errors="coerce")
    scales = hosts["frb"].map(plate_scale_arcsec)
    return (re_px * scales).to_numpy(dtype=float)


def frb_n(hosts: pd.DataFrame) -> np.ndarray:
    return pd.to_numeric(hosts["n"], errors="coerce").to_numpy(dtype=float)


def quartile_bin_edges(null_vals: np.ndarray) -> np.ndarray:
    v = np.asarray(null_vals, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) < 4:
        raise RuntimeError("Too few null values for quartile bins.")
    return np.quantile(v, [0.0, 0.25, 0.5, 0.75, 1.0])


def inclination_bin_labels() -> list[str]:
    return [
        r"$0° \leq i < 22.5°$",
        r"$22.5° \leq i < 45°$",
        r"$45° \leq i < 67.5°$",
        r"$67.5° \leq i \leq 90°$",
    ]


def fraction_in_angle_bins(angles_deg: np.ndarray, edges: np.ndarray | None = None) -> np.ndarray:
    """Fraction of sample in each equal-width inclination bin (sums to 1)."""
    if edges is None:
        edges = np.array([0.0, 22.5, 45.0, 67.5, 90.0])
    v = np.asarray(angles_deg, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return np.zeros(len(edges) - 1)
    counts, _ = count_in_bins(v, edges)
    total = counts.sum()
    return counts / total if total > 0 else counts


def plot_inclination_bin_fractions(
    *,
    frb_inc: np.ndarray,
    sdss_inc: np.ndarray,
    n_frb: int,
    n_sdss: int,
    out_stem: Path,
    sample_mode: str = "inclusive",
) -> None:
    """
    Grouped bar chart: fraction of objects per 22.5° inclination bin.

    Compares FRB hosts (GALFIT) to SDSS null (inclination from expAB_r, inclusive cuts).
    Uses **fractions** (not raw counts) so pools of different size are comparable.
    """
    edges = np.array([0.0, 22.5, 45.0, 67.5, 90.0])
    labels = inclination_bin_labels()
    frb_frac = fraction_in_angle_bins(frb_inc, edges)
    sdss_frac = fraction_in_angle_bins(sdss_inc, edges)

    x = np.arange(len(labels))
    width = 0.36
    fig, ax = plt.subplots(figsize=(9, 5))
    bars_frb = ax.bar(
        x - width / 2,
        frb_frac * 100,
        width,
        color="#d62728",
        label=f"FRB hosts ($N={n_frb}$)",
        edgecolor="white",
        linewidth=0.8,
    )
    bars_sdss = ax.bar(
        x + width / 2,
        sdss_frac * 100,
        width,
        color="#4daf4a",
        label=f"SDSS null ($N={n_sdss}$, {sample_mode})",
        edgecolor="white",
        linewidth=0.8,
    )
    uniform_pct = 100.0 / len(labels)
    ax.axhline(uniform_pct, color="0.35", linestyle="--", linewidth=1.2, label="Uniform in $i$ (25%)")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_ylabel("Fraction of sample (%)")
    ymax = max(100.0, float(np.max([frb_frac.max(), sdss_frac.max()]) * 100 * 1.15))
    ax.set_ylim(0, ymax)
    ax.set_title(
        f"Host inclination distribution in four equal-angle bins\n"
        f"($m_r<21$, SDSS null: {sample_mode}; FRB from pipeline GALFIT)"
    )
    ax.legend(loc="upper right", fontsize=9)

    def _label_bars(bars, fracs, n_pool: int):
        for bar, frac in zip(bars, fracs):
            h = bar.get_height()
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                h + 1.0,
                f"{frac * 100:.0f}%\n(n={int(round(frac * n_pool))})",
                ha="center",
                va="bottom",
                fontsize=8,
            )

    _label_bars(bars_frb, frb_frac, n_frb)
    _label_bars(bars_sdss, sdss_frac, n_sdss)

    fig.text(
        0.5,
        0.01,
        "SDSS: $i$ from expAB_r (lnL exp-wins) via Hubble formula ($q_0=0.2$). "
        "Earlier Legacy overlay used rescaled raw counts (N_frB/N_null); replaced by this fraction plot.",
        ha="center",
        fontsize=7.5,
        color="0.35",
    )
    plt.tight_layout(rect=[0, 0.04, 1, 1])
    save_figure(fig, out_stem)
    plt.close(fig)


def quantile_cosi_bin_table(
    edges: np.ndarray,
    *,
    sdss_cosi: np.ndarray,
    frb_hosts: pd.DataFrame,
    frb_cosi: np.ndarray,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build summary and per-FRB tables for equal-count cos(i) quantile bins.

    Returns (bin_summary_df, frb_assignments_df).
    """
    edges = np.asarray(edges, dtype=float)
    n_bins = len(edges) - 1
    sdss_idx = assign_values_to_bin_edges(sdss_cosi, edges)
    frb_idx = assign_values_to_bin_edges(frb_cosi, edges)

    summary_rows = []
    for b in range(n_bins):
        cos_lo, cos_hi = float(edges[b]), float(edges[b + 1])
        i_lo = float(inc_deg_from_cosi(np.array([cos_hi]))[0])
        i_hi = float(inc_deg_from_cosi(np.array([cos_lo]))[0])
        n_sdss = int(np.sum(sdss_idx == b))
        in_bin = frb_idx == b
        summary_rows.append(
            {
                "bin": b + 1,
                "cos_i_lo": cos_lo,
                "cos_i_hi": cos_hi,
                "i_deg_lo": min(i_lo, i_hi),
                "i_deg_hi": max(i_lo, i_hi),
                "n_sdss": n_sdss,
                "frac_sdss_pct": 100.0 * n_sdss / max(1, np.sum(sdss_idx >= 0)),
                "n_frb": int(np.sum(in_bin)),
                "frac_frb_pct": 100.0 * np.sum(in_bin) / max(1, len(frb_cosi)),
            }
        )
    summary_df = pd.DataFrame(summary_rows)

    frb_rows = []
    frb_col = "frb" if "frb" in frb_hosts.columns else None
    inc = pd.to_numeric(frb_hosts.get("inc"), errors="coerce")
    for j in range(len(frb_hosts)):
        name = frb_hosts.iloc[j][frb_col] if frb_col else f"host_{j}"
        frb_rows.append(
            {
                "frb": name,
                "inc_deg": float(inc.iloc[j]) if np.isfinite(inc.iloc[j]) else np.nan,
                "cos_i": float(frb_cosi[j]),
                "bin": int(frb_idx[j]) + 1 if frb_idx[j] >= 0 else np.nan,
            }
        )
    frb_df = pd.DataFrame(frb_rows)
    return summary_df, frb_df


def plot_frb_in_sdss_quantile_cosi_bins(
    *,
    sdss_cosi: np.ndarray,
    frb_hosts: pd.DataFrame,
    frb_cosi: np.ndarray,
    edges: np.ndarray,
    n_bins: int,
    mag_limit: float,
    out_stem: Path,
    sdss_label: str = "SDSS null",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Bar chart of FRB fraction per equal-count SDSS cos(i) quantile bin.

    Bin edges are defined from the SDSS pool so each bin holds ~100/n_bins % of null galaxies.
    """
    summary_df, frb_df = quantile_cosi_bin_table(
        edges,
        sdss_cosi=sdss_cosi,
        frb_hosts=frb_hosts,
        frb_cosi=frb_cosi,
    )
    n_bins_eff = len(edges) - 1
    expected_pct = 100.0 / n_bins_eff

    x = np.arange(n_bins_eff)
    frb_frac = summary_df["frac_frb_pct"].to_numpy(dtype=float)
    n_frb = len(frb_cosi)

    fig, ax = plt.subplots(figsize=(10, 5.5))

    bars = ax.bar(
        x,
        frb_frac,
        color="#d62728",
        edgecolor="white",
        linewidth=0.8,
        label=f"FRB hosts ($N={n_frb}$)",
    )
    ax.axhline(
        expected_pct,
        color="#4daf4a",
        linestyle="--",
        linewidth=1.5,
        label=f"SDSS equal-count target ({expected_pct:.1f}% per bin)",
    )
    for bar, row in zip(bars, summary_df.itertuples()):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 1.0,
            f"{row.n_frb}\n({row.frac_frb_pct:.0f}%)",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    ax.set_ylabel("FRB fraction (%)")
    ymax = max(100.0, float(frb_frac.max()) * 1.25, expected_pct * 2.5)
    ax.set_ylim(0, ymax)
    ax.set_title(
        f"FRB hosts in SDSS-defined inclination quantile bins\n"
        f"($m_r<{mag_limit:g}$, strict $b/a>0.2$, SDSS $u-r<2.3$; "
        f"{n_bins_eff} equal-count bins)"
    )
    ax.legend(loc="upper right", fontsize=9)

    tick_labels = []
    for row in summary_df.itertuples():
        hi_bracket = "]" if int(row.bin) == n_bins_eff else ")"
        tick_labels.append(
            f"Bin {row.bin}\n"
            rf"$i$: [{row.i_deg_lo:.0f}°, {row.i_deg_hi:.0f}°{hi_bracket}"
        )
    ax.set_xticks(x)
    ax.set_xticklabels(tick_labels, fontsize=9)
    ax.set_xlabel("Inclination bin (equal SDSS counts; edges from SDSS null)")

    fig.text(
        0.5,
        0.01,
        "Bins: equal SDSS counts (~12.5% each for 8 bins). "
        "Inclination ranges from SDSS Hubble cos(i) quantile edges.",
        ha="center",
        fontsize=8,
        color="0.35",
    )
    plt.tight_layout(rect=[0, 0.04, 1, 1])
    save_figure(fig, out_stem)
    plt.close(fig)
    return summary_df, frb_df


def count_in_bins(values: np.ndarray, edges: np.ndarray) -> tuple[np.ndarray, list[str]]:
    counts = []
    labels = []
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        if i < len(edges) - 2:
            mask = (values >= lo) & (values < hi)
            labels.append(f"[{lo:.3g}, {hi:.3g})")
        else:
            mask = (values >= lo) & (values <= hi)
            labels.append(f"[{lo:.3g}, {hi:.3g}]")
        counts.append(int(np.sum(mask)))
    return np.array(counts), labels


def default_font() -> font_manager.FontProperties:
    return font_manager.FontProperties(family="Arial", style="normal", size=8)


def add_inclination_top_axis(ax, font_prop) -> None:
    ax_top = ax.twiny()
    ax_top.set_xlim(ax.get_xlim())
    tick_vals = np.cos(np.radians([90, 78, 66, 53, 37, 0]))
    ax_top.set_xticks(tick_vals)
    ax_top.set_xticklabels(["90", "78", "66", "53", "37", "0"], fontproperties=font_prop)
    ax_top.set_xlabel("Inclination angle i (degrees)", fontproperties=font_prop, fontsize=10)


def plot_inclination_cdf_overlay(
    *,
    null_label: str,
    null_color: str,
    x_null: np.ndarray,
    mean_null: np.ndarray,
    lo_null: np.ndarray,
    hi_null: np.ndarray,
    frb_draws: list[list[float]],
    x_frb: np.ndarray,
    y_frb: np.ndarray,
    n_frb: int,
    title: str,
    out_stem: Path,
    mc_alpha: float = 0.03,
) -> None:
    font_prop = default_font()
    fig, ax = plt.subplots(figsize=(8, 8))
    y_steps = [0.0] + [i / n_frb for i in range(1, n_frb + 1)] + [1.0]
    for draw in frb_draws:
        x_draw = [0.0] + draw + [1.0]
        ax.step(x_draw, y_steps, where="mid", color="#d62728", linewidth=0.9, alpha=mc_alpha)
    ax.plot(x_frb, y_frb, color="#d62728", linewidth=2.0, label=f"FRB hosts (N={n_frb})")
    ax.plot(x_null, mean_null, color=null_color, linewidth=2.0, label=null_label)
    ax.fill_between(x_null, lo_null, hi_null, color=null_color, alpha=0.22, label="Null 68% CI")
    ax.plot((0, 1), (0, 1), "k--", linewidth=1.0, label="Uniform")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("cos(i)", fontproperties=font_prop, fontsize=11)
    ax.set_ylabel("Cumulative distribution", fontproperties=font_prop, fontsize=11)
    ax.set_title(title, fontsize=12)
    ax.legend(fontsize=8, loc="upper left")
    add_inclination_top_axis(ax, font_prop)
    plt.tight_layout()
    save_figure(fig, out_stem)
    plt.close(fig)


def plot_inclination_cdf_dual_null_overlay(
    *,
    nulls: tuple[
        tuple[str, str, np.ndarray, np.ndarray, np.ndarray, np.ndarray],
        ...,
    ],
    frb_draws: list[list[float]],
    x_frb: np.ndarray,
    y_frb: np.ndarray,
    n_frb: int,
    title: str,
    out_stem: Path,
    mc_alpha: float = 0.03,
) -> None:
    """Legacy + SDSS null CDFs and FRB overlay on one panel."""
    font_prop = default_font()
    fig, ax = plt.subplots(figsize=(8, 8))
    y_steps = [0.0] + [i / n_frb for i in range(1, n_frb + 1)] + [1.0]
    for draw in frb_draws:
        x_draw = [0.0] + draw + [1.0]
        ax.step(x_draw, y_steps, where="mid", color="#d62728", linewidth=0.9, alpha=mc_alpha)
    ax.plot(x_frb, y_frb, color="#d62728", linewidth=2.0, label=f"FRB hosts (N={n_frb})")
    for label, color, x_n, mn, lo, hi in nulls:
        ax.plot(x_n, mn, color=color, linewidth=2.0, label=label)
        ax.fill_between(x_n, lo, hi, color=color, alpha=0.18)
    ax.plot((0, 1), (0, 1), "k--", linewidth=1.0, label="Uniform")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("cos(i)", fontproperties=font_prop, fontsize=11)
    ax.set_ylabel("Cumulative distribution", fontproperties=font_prop, fontsize=11)
    ax.set_title(title, fontsize=12)
    ax.legend(fontsize=8, loc="upper left")
    add_inclination_top_axis(ax, font_prop)
    plt.tight_layout()
    save_figure(fig, out_stem)
    plt.close(fig)


def plot_scalar_cdf_overlay(
    *,
    null_label: str,
    null_color: str,
    null_vals: np.ndarray,
    frb_vals: np.ndarray,
    n_frb: int,
    xlabel: str,
    title: str,
    out_stem: Path,
    n_draws_null: int = 10000,
    frb_use_mc: bool = False,
    frb_draws: list[list[float]] | None = None,
    x_frb: np.ndarray | None = None,
    y_frb: np.ndarray | None = None,
    mc_alpha: float = 0.03,
) -> None:
    frb_finite = frb_vals[np.isfinite(frb_vals)]
    x_grid = _x_grid_for_values(null_vals, frb_finite)
    x_n, mean_n, lo_n, hi_n = cdf_envelope_scalar(
        null_vals, n_sample=n_frb, x_grid=x_grid, n_draws=n_draws_null
    )
    font_prop = default_font()
    fig, ax = plt.subplots(figsize=(8, 8))

    if frb_use_mc and frb_draws is not None and x_frb is not None and y_frb is not None:
        y_steps = [0.0] + [i / n_frb for i in range(1, n_frb + 1)] + [1.0]
        for draw in frb_draws:
            x_draw = [0.0] + list(draw) + [1.0]
            ax.step(x_draw, y_steps, where="mid", color="#d62728", linewidth=0.9, alpha=mc_alpha)
        ax.plot(x_frb, y_frb, color="#d62728", linewidth=2.0, label=f"FRB hosts (N={n_frb})")
    else:
        frb_sorted = np.sort(frb_finite)
        yf = ecdf_on_grid(frb_sorted, x_grid)
        ax.plot(x_grid, yf, color="#d62728", linewidth=2.0, label=f"FRB hosts (N={n_frb})")

    ax.plot(x_n, mean_n, color=null_color, linewidth=2.0, label=null_label)
    ax.fill_between(x_n, lo_n, hi_n, color=null_color, alpha=0.22, label="Null 68% CI")
    ax.set_xlim(float(x_grid[0]), float(x_grid[-1]))
    ax.set_ylim(0, 1)
    ax.set_xlabel(xlabel, fontproperties=font_prop, fontsize=11)
    ax.set_ylabel("Cumulative distribution", fontproperties=font_prop, fontsize=11)
    ax.set_title(title, fontsize=12)
    ax.legend(fontsize=8, loc="lower right")
    plt.tight_layout()
    save_figure(fig, out_stem)
    plt.close(fig)


def plot_histogram_four_bins(
    *,
    counts: np.ndarray,
    labels: list[str],
    title: str,
    ylabel: str,
    out_stem: Path,
    null_counts: np.ndarray | None = None,
    null_label: str | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))
    x = np.arange(len(counts))
    width = 0.35 if null_counts is not None else 0.6
    ax.bar(x - (width / 2 if null_counts is not None else 0), counts, width=width, color="#d62728", label="FRB hosts")
    if null_counts is not None and null_label:
        ax.bar(x + width / 2, null_counts, width=width, color="#377eb8", alpha=0.7, label=null_label)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=25, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    for i, c in enumerate(counts):
        ax.text(i - (width / 2 if null_counts is not None else 0), c + 0.3, str(int(c)), ha="center", fontsize=9)
    ax.legend()
    plt.tight_layout()
    save_figure(fig, out_stem)
    plt.close(fig)


def plot_random_host_inclination(
    *,
    null_cosi: np.ndarray,
    frb_draws: list[list[float]],
    x_frb: np.ndarray,
    y_frb: np.ndarray,
    n_frb: int,
    null_label: str,
    title: str,
    out_stem: Path,
    n_random: int = 200,
    n_draws_null: int = 10000,
) -> float:
    rng = np.random.default_rng(42)
    ref = null_cosi[np.isfinite(null_cosi)]
    x_n, mean_n, lo_n, hi_n = cdf_envelope(ref, n_sample=n_frb, n_draws=n_draws_null)

    fig, ax = plt.subplots(figsize=(8, 8))
    y_steps = [0.0] + [i / n_frb for i in range(1, n_frb + 1)] + [1.0]
    random_means = []
    for _ in range(n_random):
        draw = np.sort(rng.choice(ref, size=n_frb, replace=False))
        x_draw = [0.0] + draw.tolist() + [1.0]
        ax.step(x_draw, y_steps, where="mid", color="0.7", linewidth=0.6, alpha=0.25)
        random_means.append(float(np.mean(draw)))

    ax.plot(x_frb, y_frb, color="#d62728", linewidth=2.5, label=f"FRB hosts (N={n_frb})")
    ax.plot(x_n, mean_n, color="#377eb8", linewidth=2.0, label=null_label)
    ax.fill_between(x_n, lo_n, hi_n, color="#377eb8", alpha=0.2, label="Null 68% CI")
    frb_mean = float(np.mean([np.mean(d) for d in frb_draws]))
    p_close = float(np.mean([abs(m - frb_mean) <= 0.02 for m in random_means]))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("cos(i)")
    ax.set_ylabel("Cumulative distribution")
    ax.set_title(f"{title}\n(random mean within 0.02 of FRB: {p_close:.1%})")
    ax.legend(fontsize=8, loc="upper left")
    font_prop = default_font()
    add_inclination_top_axis(ax, font_prop)
    plt.tight_layout()
    save_figure(fig, out_stem)
    plt.close(fig)
    return p_close


def plot_delta_cleveland(
    diff: pd.DataFrame,
    delta_col: str,
    param_label: str,
    out_stem: Path,
) -> None:
    df = diff.copy()
    df["abs_delta"] = pd.to_numeric(df[delta_col], errors="coerce").abs()
    df = df.sort_values(delta_col)
    y = np.arange(len(df))
    colors = np.where(df["single_sersic"], "#377eb8", "#ff7f00")

    fig, ax = plt.subplots(figsize=(8, max(6, len(df) * 0.22)))
    ax.hlines(y, 0, df[delta_col], colors=colors, linewidth=1.2)
    ax.plot(df[delta_col], y, "o", color="black", markersize=4)
    ax.axvline(0, color="0.3", linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(df["frb"].astype(str), fontsize=7)
    ax.set_xlabel(f"Δ {param_label} (pipeline − master)")
    ax.set_title(f"Pipeline vs master: {param_label}")
    from matplotlib.lines import Line2D

    ax.legend(
        handles=[
            Line2D([0], [0], color="#377eb8", lw=2, label="single Sérsic"),
            Line2D([0], [0], color="#ff7f00", lw=2, label="multi Sérsic"),
        ],
        loc="lower right",
        fontsize=8,
    )
    plt.tight_layout()
    save_figure(fig, out_stem)
    plt.close(fig)


def plot_sky_map(
    hosts: pd.DataFrame,
    legacy_cut: pd.DataFrame,
    sdss_cut: pd.DataFrame,
    *,
    title: str,
    out_stem: Path,
    max_null_points: int = 15000,
    seed: int = 42,
) -> None:
    rng = np.random.default_rng(seed)
    fig = plt.figure(figsize=(12, 6))
    ax = fig.add_subplot(111, projection="mollweide")

    def subsample(df: pd.DataFrame) -> pd.DataFrame:
        if len(df) <= max_null_points:
            return df
        return df.sample(n=max_null_points, random_state=seed)

    for df, color, label in (
        (subsample(sdss_cut), "#4daf4a", "SDSS null"),
        (subsample(legacy_cut), "#377eb8", "Legacy null"),
    ):
        ra = np.radians(pd.to_numeric(df["RA_ICRS"], errors="coerce"))
        dec = np.radians(pd.to_numeric(df["DE_ICRS"], errors="coerce"))
        ok = np.isfinite(ra) & np.isfinite(dec)
        ax.scatter(
            ra[ok],
            dec[ok],
            s=1,
            c=color,
            alpha=0.08,
            linewidths=0,
            label=label,
            rasterized=True,
        )

    ra_h = np.radians(pd.to_numeric(hosts["ra_deg"], errors="coerce"))
    dec_h = np.radians(pd.to_numeric(hosts["dec_deg"], errors="coerce"))
    ax.scatter(ra_h, dec_h, s=80, c="#d62728", edgecolors="white", linewidths=0.5, label="FRB hosts", zorder=5)
    for r, d, name in zip(ra_h, dec_h, hosts["frb"].astype(str)):
        if np.isfinite(r) and np.isfinite(d):
            ax.text(r, d, f"  {name}", fontsize=5, color="#d62728", zorder=6)

    ax.set_title(title)
    ax.legend(loc="lower center", ncol=3, fontsize=8)
    ax.grid(True, alpha=0.3)
    note = ""
    if len(legacy_cut) > max_null_points or len(sdss_cut) > max_null_points:
        note = f" (null subsampled to {max_null_points}/survey)"
    fig.text(0.5, 0.02, note, ha="center", fontsize=8)
    plt.tight_layout()
    save_figure(fig, out_stem)
    plt.close(fig)


def load_null_cuts(
    legacy_csv: Path | None,
    sdss_csv: pd.DataFrame | Path | None,
    *,
    sample_mode: str,
    mag_limit: float,
    mag_column: str,
    q0: float,
    exclude_types: str,
    sdss_q_column: str,
    sdss_mag_column: str | None = None,
    exclude_sdss_ba_floor: bool = False,
    sdss_ur_max: float | None = SDSS_UR_MAX_CDF,
    legacy_gr_max: float | None = LEGACY_GR_MAX_CDF,
    legacy_df: pd.DataFrame | None = None,
    sdss_df: pd.DataFrame | None = None,
    extended_columns: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if legacy_df is None:
        if legacy_csv is None:
            raise ValueError("legacy_csv or legacy_df required")
        legacy_df = read_legacy_null_catalog(legacy_csv, extended=extended_columns)
    if sdss_df is None:
        if sdss_csv is None:
            raise ValueError("sdss_csv or sdss_df required")
        if isinstance(sdss_csv, Path):
            sdss_df = read_sdss_null_catalog(sdss_csv, extended=extended_columns)
        else:
            sdss_df = ensure_sdss_n_eff(sdss_csv)
    sdss_mag = sdss_mag_column or mag_column
    legacy_cut = prepare_null_sample(
        legacy_df,
        sample_mode=sample_mode,
        mag_column=mag_column,
        mag_limit=mag_limit,
        q0=q0,
        exclude_legacy_types=exclude_types,
        is_legacy=True,
        sdss_ur_max=None,
        legacy_gr_max=legacy_gr_max,
    )
    sdss_cut = prepare_null_sample(
        sdss_df,
        sample_mode=sample_mode,
        mag_column=sdss_mag,
        mag_limit=mag_limit,
        q0=q0,
        q_column=sdss_q_column,
        is_legacy=False,
        exclude_sdss_ba_floor=exclude_sdss_ba_floor,
        sdss_ur_max=sdss_ur_max,
        legacy_gr_max=None,
        sdss_exp_winner_only=True,
    )
    return legacy_cut, sdss_cut
