"""Aggregate per-host verification JSON into cohort tables, plots and flags.

The population-level parts of checks 1, 5, 6 and 7 live here rather than in the
checks themselves: a leakage test or a trend slope is meaningless one galaxy at
a time.
"""

from __future__ import annotations

import glob
import math
import os
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

VER_DIR = os.path.dirname(os.path.abspath(__file__))
if VER_DIR not in sys.path:
    sys.path.insert(0, VER_DIR)

import vercommon as vc  # noqa: E402

CHECKS = ("chi2", "rff", "fourier", "psf", "mag", "isophote", "sky", "astrophot")
PLOTS_ROOT = os.path.join(vc.OUT_ROOT, "plots")

# Provisional. Deliberately loose until they are calibrated on the real
# distribution of each metric; see FIT_VERIFICATION_CHECKS.md.
THRESHOLDS = {
    "rff_2re": 0.10,
    "fourier_dq": 0.05,
    "iso_dq_2re": 0.10,
    "dq_sky": 0.05,
    "dq_astrophot": 0.05,
    "dmag_ref": 0.5,
    "sigma_calibration_ratio": (0.5, 2.0),
}


CONFIRMATION_CSV = os.path.join(VER_DIR, "host_confirmation.csv")


def _load_confirmation() -> pd.DataFrame:
    """Per-host paper inclusion decisions (maintained by visual triage)."""
    if not os.path.isfile(CONFIRMATION_CSV):
        return pd.DataFrame(columns=["frb", "confirmed", "notes"])
    conf = pd.read_csv(CONFIRMATION_CSV)
    if conf.empty or "frb" not in conf.columns:
        return pd.DataFrame(columns=["frb", "confirmed", "notes"])
    conf = conf.copy()
    conf["frb"] = conf["frb"].astype(str)
    if "confirmed" not in conf.columns:
        conf["confirmed"] = pd.NA
    else:
        # Accept True/False/yes/no/1/0; blank -> not yet reviewed.
        conf["confirmed"] = conf["confirmed"].map(
            lambda v: (True if str(v).strip().lower() in ("1", "true", "yes", "y")
                       else False if str(v).strip().lower() in ("0", "false", "no", "n")
                       else pd.NA)
            if pd.notna(v) and str(v).strip() != "" else pd.NA
        )
    if "notes" not in conf.columns:
        conf["notes"] = ""
    return conf[["frb", "confirmed", "notes"]].drop_duplicates("frb", keep="last")


def collect(cohort: str = "all64") -> pd.DataFrame:
    """One row per host: every scalar every check wrote."""
    df = vc.cohort(cohort)
    base = df[["frb", "in_53", "snr_win", "mag", "b_a", "b_a_err", "re", "n",
               "pa", "chi2nu", "zp_ok", "ref_survey", "n_sersic_components"]].copy()
    # Checks re-report a few production columns for convenience; keep the
    # authoritative copy from the results CSV so the merge cannot suffix them.
    reserved = set(base.columns)
    rows = []
    for frb in base["frb"].astype(str):
        row: dict = {"frb": frb}
        hdir = vc.per_host_dir(frb, create=False)
        for check in CHECKS:
            payload = vc.read_json(os.path.join(hdir, f"{check}.json"))
            row[f"{check}_status"] = payload.get("status", "missing")
            for key, val in payload.items():
                if key.startswith("_") or key in ("status", "traceback", "error"):
                    continue
                if key in reserved and key != "frb":
                    continue
                if isinstance(val, (int, float, bool, str)) or val is None:
                    row.setdefault(key, val)
        rows.append(row)
    metrics = pd.DataFrame(rows)
    out = base.merge(metrics, on="frb", how="left")
    # Explicit triad for §4.4: older sky.json may omit q_sky_0.
    if "q_sky_0" not in out.columns:
        out["q_sky_0"] = out["b_a"]
    else:
        out["q_sky_0"] = out["q_sky_0"].fillna(out["b_a"])
    # Paper-inclusion triage (blank = not yet reviewed).
    conf = _load_confirmation()
    out = out.merge(conf, on="frb", how="left")
    return out


def _num(df: pd.DataFrame, col: str) -> np.ndarray:
    if col not in df.columns:
        return np.full(len(df), np.nan)
    return pd.to_numeric(df[col], errors="coerce").to_numpy(dtype=float)


def _spearman(x: np.ndarray, y: np.ndarray) -> tuple[float, float, int]:
    ok = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(ok) < 5:
        return float("nan"), float("nan"), int(np.count_nonzero(ok))
    r, p = stats.spearmanr(x[ok], y[ok])
    return float(r), float(p), int(np.count_nonzero(ok))


def _wls_slope(x: np.ndarray, y: np.ndarray) -> dict:
    ok = np.isfinite(x) & np.isfinite(y)
    n = int(np.count_nonzero(ok))
    if n < 5:
        return {"slope": float("nan"), "slope_err": float("nan"),
                "slope_sig": float("nan"), "n": n}
    res = stats.linregress(x[ok], y[ok])
    return {
        "slope": float(res.slope),
        "slope_err": float(res.stderr),
        "slope_sig": float(res.slope / res.stderr) if res.stderr > 0 else float("nan"),
        "r": float(res.rvalue),
        "p": float(res.pvalue),
        "n": n,
    }


def population(df: pd.DataFrame) -> dict:
    """The correlations and trends that only exist across the cohort."""
    sel = df["in_53"].fillna(False).to_numpy(dtype=bool)
    out: dict = {"n_hosts": int(len(df)), "n_in_53": int(sel.sum())}

    for tag, mask in (("all64", np.ones(len(df), bool)), ("in53", sel)):
        d = df[mask]
        # Check 6: is q tracking the instrument?
        r, p, n = _spearman(_num(d, "q_host"), _num(d, "psf_ellipticity"))
        out[f"psf_q_vs_epsf_spearman_{tag}"] = r
        out[f"psf_q_vs_epsf_p_{tag}"] = p
        out[f"psf_q_vs_epsf_n_{tag}"] = n

        r, p, n = _spearman(_num(d, "q_host"), _num(d, "fwhm_over_re"))
        out[f"psf_q_vs_fwhm_over_re_spearman_{tag}"] = r
        out[f"psf_q_vs_fwhm_over_re_p_{tag}"] = p

        # Alignment is the sharper test: an amplitude correlation dilutes
        # easily, but PA_host tracking PA_psf is unambiguous.
        dpa = _num(d, "dpa_host_psf_deg")
        dpa = dpa[np.isfinite(dpa)]
        if dpa.size >= 5:
            ks = stats.kstest(dpa / 90.0, "uniform")
            out[f"psf_dpa_ks_stat_{tag}"] = float(ks.statistic)
            out[f"psf_dpa_ks_p_{tag}"] = float(ks.pvalue)
            out[f"psf_dpa_median_deg_{tag}"] = float(np.median(dpa))
            out[f"psf_dpa_n_{tag}"] = int(dpa.size)

        # Check 7: magnitude leakage trends.
        dmag = _num(d, "dmag_ref")
        for name, xcol in (("re", "re_arcsec"), ("n", "n_sersic"),
                           ("sky", "sky_offset_sigma"), ("chi2nu", "chi2nu_global")):
            fit = _wls_slope(_num(d, xcol), dmag)
            for key, val in fit.items():
                out[f"dmag_vs_{name}_{key}_{tag}"] = val

        # Check 1: chi2nu should track SNR, not model failure.
        for name, xcol in (("snr", "snr_win"), ("re_over_fwhm", "re_over_fwhm"),
                           ("mag", "mag")):
            r, p, n = _spearman(_num(d, xcol), _num(d, "chi2nu_local_2re"))
            out[f"chi2nu_local_vs_{name}_spearman_{tag}"] = r
            out[f"chi2nu_local_vs_{name}_p_{tag}"] = p

        # Check 5 / 8: do the independent estimates of q agree?
        for name, col in (("isophote", "iso_q_at_1re_data"),
                          ("astrophot", "ap_q")):
            dq = _num(d, col) - _num(d, "b_a")
            dq = dq[np.isfinite(dq)]
            if dq.size:
                out[f"q_vs_{name}_median_offset_{tag}"] = float(np.median(dq))
                out[f"q_vs_{name}_scatter_{tag}"] = float(
                    1.4826 * np.median(np.abs(dq - np.median(dq)))
                )
                out[f"q_vs_{name}_n_{tag}"] = int(dq.size)

    for col in ("rff_2re", "fourier_dq", "dq_sky", "dq_astrophot",
                "iso_dq_2re", "sigma_calibration_ratio", "chi2nu_local_2re",
                "dmag_ref"):
        vals = _num(df[sel], col)
        vals = vals[np.isfinite(vals)]
        if vals.size:
            out[f"{col}_median_in53"] = float(np.median(vals))
            out[f"{col}_p16_in53"] = float(np.percentile(vals, 16))
            out[f"{col}_p84_in53"] = float(np.percentile(vals, 84))
    return out


def _bool_col_of(df: pd.DataFrame, name: str) -> np.ndarray:
    """Booleans round-tripped through JSON and CSV arrive as mixed types."""
    if name not in df.columns:
        return np.zeros(len(df), dtype=bool)
    return np.array([v is True or v is np.True_ or str(v) == "True"
                     for v in df[name].tolist()], dtype=bool)


def flags(df: pd.DataFrame) -> pd.DataFrame:
    """Per-check booleans and a provisional trust tier.

    Thresholds are provisional; the tier is a triage aid, not a verdict.
    """
    f = pd.DataFrame({"frb": df["frb"], "in_53": df["in_53"]})
    lo, hi = THRESHOLDS["sigma_calibration_ratio"]
    ratio = _num(df, "sigma_calibration_ratio")

    f["flag_sigma_miscalibrated"] = ~((ratio > lo) & (ratio < hi))
    f["flag_rff_high"] = np.abs(_num(df, "rff_2re")) > THRESHOLDS["rff_2re"]
    f["flag_rff_sky"] = np.abs(_num(df, "rff_outer_minus_inner")) > THRESHOLDS["rff_2re"]
    # Only count the Fourier dq where the estimator had leverage; otherwise a
    # barely-resolved host would be condemned by a number that means nothing.
    reliable = _bool_col_of(df, "fourier_reliable")
    f["fourier_reliable"] = reliable
    f["flag_fourier_dq"] = reliable & (np.abs(_num(df, "fourier_dq"))
                                       > THRESHOLDS["fourier_dq"])
    f["flag_fourier_dq_significant"] = reliable & (np.abs(_num(df, "fourier_dq_sig")) > 3.0)
    f["flag_fourier_unusable"] = ~reliable
    f["flag_iso_dq"] = np.abs(_num(df, "iso_dq_2re")) > THRESHOLDS["iso_dq_2re"]
    f["flag_sky_sensitive"] = _num(df, "dq_sky") > THRESHOLDS["dq_sky"]
    f["flag_astrophot_disagrees"] = np.abs(_num(df, "dq_astrophot")) > THRESHOLDS["dq_astrophot"]
    f["flag_dmag"] = np.abs(_num(df, "dmag_ref")) > THRESHOLDS["dmag_ref"]
    f["flag_param_at_bound"] = (_bool_col_of(df, "n_at_bound")
                                | _bool_col_of(df, "re_at_bound"))
    f["flag_q_near_floor"] = _num(df, "b_a") <= 0.25
    f["flag_unresolved"] = _num(df, "re_over_fwhm") < 1.0
    # "unresolved" and "stamp_too_small" are statements about the data, not
    # failures of the suite, so they must not masquerade as errors.
    benign = {"ok", "unresolved", "stamp_too_small"}
    f["flag_check_error"] = np.any(
        np.column_stack([~df[f"{c}_status"].fillna("missing").isin(benign).to_numpy()
                         for c in CHECKS if f"{c}_status" in df.columns]),
        axis=1,
    )

    # Only the flags that speak to whether q itself is trustworthy count toward
    # the tier. A high RFF on a bright spiral is astrophysics, not a bad fit.
    geometry = ["flag_fourier_dq", "flag_iso_dq", "flag_sky_sensitive",
                "flag_astrophot_disagrees"]
    f["n_geometry_flags"] = f[geometry].sum(axis=1)
    f["n_flags_total"] = f[[c for c in f.columns if c.startswith("flag_")]].sum(axis=1)
    # A host whose checks did not all complete cannot be graded: without that
    # distinction an unrun host would silently score as clean.
    f["trust_tier"] = np.where(
        f["flag_check_error"], "?",
        np.where(f["n_geometry_flags"] == 0, "A",
                 np.where(f["n_geometry_flags"] == 1, "B", "C")),
    )
    return f


def _scatter(ax, x, y, c=None, xlabel="", ylabel="", title="", one_to_one=False,
             cbar_label=""):
    ok = np.isfinite(x) & np.isfinite(y)
    kw = {"s": 26, "edgecolor": "k", "linewidth": 0.4}
    if c is not None and np.any(np.isfinite(c)):
        sc = ax.scatter(x[ok], y[ok], c=c[ok], cmap="viridis", **kw)
        cb = ax.figure.colorbar(sc, ax=ax, fraction=0.04, pad=0.02)
        cb.set_label(cbar_label, fontsize=8)
    else:
        ax.scatter(x[ok], y[ok], color="tab:blue", **kw)
    if one_to_one and np.any(ok):
        lim = [np.nanmin([x[ok].min(), y[ok].min()]),
               np.nanmax([x[ok].max(), y[ok].max()])]
        ax.plot(lim, lim, "k--", lw=0.9)
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(ylabel, fontsize=9)
    ax.set_title(title, fontsize=10)
    ax.tick_params(labelsize=8)


def population_plots(df: pd.DataFrame) -> list[str]:
    os.makedirs(PLOTS_ROOT, exist_ok=True)
    sel = df["in_53"].fillna(False).to_numpy(dtype=bool)
    d = df[sel]
    paths = []

    fig, axes = plt.subplots(2, 3, figsize=(15, 8.5))
    _scatter(axes[0, 0], _num(d, "snr_win"), _num(d, "chi2nu_local_2re"),
             xlabel="SNR (win)", ylabel=r"$\chi^2/\nu$ within $2R_e$",
             title="localized $\\chi^2/\\nu$ is an SNR meter")
    axes[0, 0].set_xscale("log"); axes[0, 0].set_yscale("log")

    _scatter(axes[0, 1], _num(d, "re_over_fwhm"), _num(d, "rff_2re"),
             xlabel=r"$R_e$ / FWHM", ylabel=r"RFF within $2R_e$",
             title="residual structure vs resolution")
    axes[0, 1].set_xscale("log")

    _scatter(axes[0, 2], _num(d, "psf_ellipticity"), _num(d, "q_host"),
             c=_num(d, "fwhm_over_re"), xlabel=r"$e_{\rm PSF}$", ylabel="$q$ (GALFIT)",
             title="PSF leakage: $q$ vs PSF ellipticity", cbar_label="FWHM/$R_e$")

    dpa = _num(d, "dpa_host_psf_deg")
    axes[1, 0].hist(dpa[np.isfinite(dpa)], bins=9, range=(0, 90),
                    color="tab:blue", edgecolor="k")
    axes[1, 0].axhline(np.count_nonzero(np.isfinite(dpa)) / 9.0, color="r",
                       ls="--", lw=1.0, label="uniform")
    axes[1, 0].set_xlabel(r"$|{\rm PA}_{\rm host} - {\rm PA}_{\rm PSF}|$ [deg]",
                          fontsize=9)
    axes[1, 0].set_ylabel("hosts", fontsize=9)
    axes[1, 0].set_title("alignment with the PSF (flat = no leakage)", fontsize=10)
    axes[1, 0].legend(fontsize=8, frameon=False)

    _scatter(axes[1, 1], _num(d, "b_a"), _num(d, "ap_q"), c=_num(d, "fwhm_over_re"),
             xlabel="$q$ (GALFIT)", ylabel="$q$ (AstroPhot)",
             title="independent refit", one_to_one=True, cbar_label="FWHM/$R_e$")

    _scatter(axes[1, 2], _num(d, "b_a"), _num(d, "iso_q_at_1re_data"),
             c=_num(d, "fwhm_over_re"), xlabel="$q$ (GALFIT)",
             ylabel=r"$q$ (isophote, $1R_e$)", title="isophotal comparison",
             one_to_one=True, cbar_label="FWHM/$R_e$")

    fig.suptitle(f"Fit verification — population diagnostics ({int(sel.sum())} hosts, "
                 "mag $\\leq$ 22 and $b/a$ > 0.2)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    path = os.path.join(PLOTS_ROOT, "population_diagnostics.png")
    fig.savefig(path, dpi=120); plt.close(fig)
    paths.append(path)

    fig, axes = plt.subplots(1, 4, figsize=(17, 4.0))
    for ax, (xcol, xlabel) in zip(axes, (
        ("re_arcsec", r"$R_e$ [arcsec]"), ("n_sersic", "Sersic $n$"),
        ("sky_offset_sigma", r"sky offset [$\sigma_{\rm sky}$]"),
        ("chi2nu_global", r"$\chi^2/\nu$"),
    )):
        _scatter(ax, _num(d, xcol), _num(d, "dmag_ref"), xlabel=xlabel,
                 ylabel=r"$\Delta$mag (GALFIT $-$ ref)")
        ax.axhline(0, color="k", lw=0.8)
    fig.suptitle("Magnitude leakage", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    path = os.path.join(PLOTS_ROOT, "mag_leakage.png")
    fig.savefig(path, dpi=120); plt.close(fig)
    paths.append(path)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    for ax, (col, label) in zip(axes, (
        ("fourier_dq", r"$\delta q$ from the $m{=}2$ residual"),
        ("dq_sky", r"$|\Delta q|$ from $\pm1\sigma$ sky"),
        ("dq_astrophot", r"$q_{\rm AstroPhot} - q_{\rm GALFIT}$"),
    )):
        vals = _num(d, col)
        vals = vals[np.isfinite(vals)]
        if vals.size:
            ax.hist(vals, bins=20, color="tab:purple", edgecolor="k")
            ax.axvline(0, color="k", lw=0.9)
            ax.axvline(float(np.median(vals)), color="r", ls="--", lw=1.0,
                       label=f"median {np.median(vals):+.4f}")
            ax.legend(fontsize=8, frameon=False)
        ax.set_xlabel(label, fontsize=9)
        ax.set_ylabel("hosts", fontsize=9)
    fig.suptitle("Three independent handles on the axis-ratio error", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    path = os.path.join(PLOTS_ROOT, "dq_comparison.png")
    fig.savefig(path, dpi=120); plt.close(fig)
    paths.append(path)
    return paths


def main(cohort: str = "all64") -> pd.DataFrame:
    os.makedirs(vc.TABLES_ROOT, exist_ok=True)
    df = collect(cohort)
    metrics_path = os.path.join(vc.TABLES_ROOT, "fit_verification_metrics.csv")
    df.to_csv(metrics_path, index=False)

    fl = flags(df)
    flags_path = os.path.join(vc.TABLES_ROOT, "fit_verification_flags.csv")
    fl.to_csv(flags_path, index=False)

    pop = population(df)
    vc.write_json(os.path.join(vc.TABLES_ROOT, "population_summary.json"), pop)

    try:
        plot_paths = population_plots(df)
    except Exception as exc:
        plot_paths = []
        print(f"[warn] population plots failed: {type(exc).__name__}: {exc}")

    have_panels = sorted(
        os.path.basename(os.path.dirname(p))
        for p in glob.glob(os.path.join(vc.PER_HOST_ROOT, "*", "panel.png"))
    )
    if have_panels:
        try:
            from checks.visual import contact_sheet

            sheet = contact_sheet(have_panels,
                                  os.path.join(PLOTS_ROOT, "contact_sheet.png"))
            plot_paths.append(sheet)
        except Exception as exc:
            print(f"[warn] contact sheet failed: {type(exc).__name__}: {exc}")

    print(f"\nWrote {metrics_path} ({len(df)} rows, {len(df.columns)} columns)")
    print(f"Wrote {flags_path}")
    for path in plot_paths:
        print(f"Wrote {path}")
    print("\nTrust tiers (53-host cut): "
          + ", ".join(f"{k}={v}" for k, v in
                      fl[fl["in_53"].fillna(False)]["trust_tier"]
                      .value_counts().sort_index().items()))
    return df


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "all64")
