"""Check 9 — per-host diagnostic panel.

Data and model use an asinh stretch on sky-subtracted flux over the data's
[1, 99] percentile window (grayscale). Residual stays in sigma units on
``RdBu_r`` clipped at +/-5.

Consumes the products of the earlier checks, so it runs last.
"""

from __future__ import annotations

import csv
import json
import math
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.ndimage import map_coordinates  # noqa: E402

import vercommon as vc  # noqa: E402

NAME = "visual"

RES_CLIP = 5.0
_ASINH_SOFT = 10.0
_ALT_PANEL = re.compile(
    r"outputs/panels/([A-Za-z0-9]+_(?:n1_sky|n1_psf|n1_moffat|n1|sky|psf|sersic|moffat))\.png"
)
_HC_CSV = os.path.join(os.path.dirname(os.path.abspath(vc.__file__)),
                       "host_confirmation.csv")
_REFITS = os.path.join(os.path.dirname(os.path.abspath(vc.__file__)), "Re-fits")


def _asinh_display(img: np.ndarray, lo: float, hi: float) -> np.ndarray:
    """Map ``img`` into [0, 1] with asinh stretch over ``[lo, hi]``."""
    x = np.clip((img - lo) / (hi - lo + 1e-30), 0.0, 1.0)
    return np.arcsinh(_ASINH_SOFT * x) / np.arcsinh(_ASINH_SOFT)


def _lotz2008_class(gini: float, m20: float) -> str | None:
    """Lotz et al. (2008) Gini–M20 morphology class (Sb–Irr / E–Sa / merger)."""
    if not (isinstance(gini, (int, float)) and isinstance(m20, (int, float))
            and math.isfinite(float(gini)) and math.isfinite(float(m20))):
        return None
    g, m = float(gini), float(m20)
    if g > -0.14 * m + 0.33:
        return "merger"
    if g > 0.14 * m + 0.80:
        return "early"  # E/S0/Sa
    return "late"  # Sb-Irr


def _lotz_label(host: vc.HostData) -> str:
    """Compact Lotz+2008 class from Phase ``statmorph_results.json`` (if present)."""
    path = os.path.join(vc.OUTPUT_ROOT, f"{host.frb}_all", "statmorph_results.json")
    if not os.path.isfile(path):
        return r"${\rm Lotz}=-$"
    try:
        with open(path, encoding="utf-8") as f:
            sm = json.load(f)
    except Exception:
        return r"${\rm Lotz}=-$"
    cls = _lotz2008_class(sm.get("gini"), sm.get("m20"))
    return rf"${{\rm Lotz}}={{\rm {cls}}}$" if cls else r"${\rm Lotz}=-$"


def _load_npz(outdir: str, name: str) -> dict | None:
    path = os.path.join(outdir, name)
    if not os.path.isfile(path):
        return None
    try:
        with np.load(path) as data:
            return {k: data[k] for k in data.files}
    except Exception:
        return None


def _ellipse_overlay(ax, host: vc.HostData, n_re: float, **kw) -> None:
    t = np.linspace(0, 2 * np.pi, 256)
    phi = math.radians(host.pa)
    a = n_re * host.re
    x_maj, y_min = a * np.cos(t), a * host.q * np.sin(t)
    ax.plot(host.xc - x_maj * np.sin(phi) + y_min * np.cos(phi),
            host.yc + x_maj * np.cos(phi) + y_min * np.sin(phi), **kw)


def _axis_cut_mu(image: np.ndarray, sigma: np.ndarray, host: vc.HostData,
                 sky: float, a_re: np.ndarray, *, minor: bool
                 ) -> tuple[np.ndarray, np.ndarray]:
    """Sky-subtracted SB along the major or minor axis vs elliptical a/Re.

    Both sides of the centre are averaged. Minor-axis samples sit at the same
    elliptical ``a`` as the major-axis ones (physical offset ``a*q``), so a
    correct elliptical model makes the two curves coincide.
    """
    phi = math.radians(host.pa)
    cos_p, sin_p = math.cos(phi), math.sin(phi)
    qs = max(float(host.q), 1e-3)
    ny, nx = image.shape
    mu = np.full(a_re.shape, np.nan)
    mu_err = np.full(a_re.shape, np.nan)
    img = np.asarray(image, dtype=float) - sky
    sig = np.asarray(sigma, dtype=float)
    for i, ar in enumerate(a_re):
        a = float(ar) * host.re
        fluxes, errs = [], []
        for sign in (+1.0, -1.0):
            if minor:
                x_maj, y_min = 0.0, sign * a * qs
            else:
                x_maj, y_min = sign * a, 0.0
            dx = -x_maj * sin_p + y_min * cos_p
            dy = x_maj * cos_p + y_min * sin_p
            x, y = host.xc + dx, host.yc + dy
            if not (0.5 <= x <= nx - 1.5 and 0.5 <= y <= ny - 1.5):
                continue
            # map_coordinates uses (row, col) = (y, x).
            f = float(map_coordinates(img, [[y], [x]], order=1, mode="nearest")[0])
            s = float(map_coordinates(sig, [[y], [x]], order=1, mode="nearest")[0])
            if math.isfinite(f) and math.isfinite(s) and s > 0 and f > 0:
                fluxes.append(f)
                errs.append(s)
        if not fluxes:
            continue
        fmean = float(np.mean(fluxes))
        # Average of two sides: variance adds, then /N^2.
        ferr = float(math.sqrt(np.sum(np.square(errs))) / len(fluxes))
        mu[i] = float(vc.counts_to_mu(np.array([fmean]), host.magzpt,
                                      host.plate_scale)[0])
        # dmu = (2.5/ln10) * dF/F
        mu_err[i] = (2.5 / math.log(10.0)) * (ferr / fmean) if fmean > 0 else np.nan
    return mu, mu_err


def _fmt(val, fmt: str, default: str = "—") -> str:
    if isinstance(val, (int, float)) and math.isfinite(float(val)):
        return fmt.format(float(val))
    return default


def _refit_fix_marker(outdir: str) -> str | None:
    """Compact tag for Re-fit panels: ``[n1]``, ``[sky]``, or ``[n1+sky]``."""
    meta = vc.read_json(os.path.join(outdir, "refit_meta.json"))
    if not meta:
        return None
    n_fixed = meta.get("n_fixed")
    sky_fixed = meta.get("sky_fixed_adu")
    has_n = isinstance(n_fixed, (int, float)) and math.isfinite(float(n_fixed))
    has_sky = isinstance(sky_fixed, (int, float)) and math.isfinite(float(sky_fixed))
    if has_n and has_sky:
        return "[n1+sky]" if abs(float(n_fixed) - 1.0) < 1e-6 else f"[n={float(n_fixed):g}+sky]"
    if has_n:
        return "[n1]" if abs(float(n_fixed) - 1.0) < 1e-6 else f"[n={float(n_fixed):g}]"
    if has_sky:
        return "[sky]"
    if meta.get("psf_added"):
        return "[psf]"
    if meta.get("sersic_added"):
        return "[sersic]"
    if meta.get("moffat_added"):
        return "[moffat]"
    return None


def _header_columns(host: vc.HostData, outdir: str) -> list[list[str]]:
    """Metric columns for the panel header band (no section titles)."""
    chi2 = vc.read_json(os.path.join(outdir, "chi2.json"))
    rff = vc.read_json(os.path.join(outdir, "rff.json"))
    fou = vc.read_json(os.path.join(outdir, "fourier.json"))
    skyj = vc.read_json(os.path.join(outdir, "sky.json"))
    ap = vc.read_json(os.path.join(outdir, "astrophot.json"))
    mag = vc.read_json(os.path.join(outdir, "mag.json"))

    ratio = chi2.get("sigma_calibration_ratio", float("nan"))
    gcorr = chi2.get("chi2nu_global_corrected", float("nan"))
    lcorr = chi2.get("chi2nu_local_2re_corrected", float("nan"))
    if not (isinstance(gcorr, (int, float)) and math.isfinite(gcorr)):
        graw = chi2.get("chi2nu_global", float("nan"))
        if (isinstance(ratio, (int, float)) and math.isfinite(ratio) and ratio > 0
                and isinstance(graw, (int, float)) and math.isfinite(graw)):
            gcorr = graw / ratio**2
    if not (isinstance(lcorr, (int, float)) and math.isfinite(lcorr)):
        lraw = chi2.get("chi2nu_local_2re", float("nan"))
        if (isinstance(ratio, (int, float)) and math.isfinite(ratio) and ratio > 0
                and isinstance(lraw, (int, float)) and math.isfinite(lraw)):
            lcorr = lraw / ratio**2

    m_gal = mag.get("mag_galfit", host.mag)
    m_ref = mag.get("ref_mag", float("nan"))
    # User convention: Δm = ref − m_GALFIT (opposite sign of dmag_ref in mag.json).
    dmag = (float(m_ref) - float(m_gal)
            if isinstance(m_ref, (int, float)) and isinstance(m_gal, (int, float))
            and math.isfinite(float(m_ref)) and math.isfinite(float(m_gal))
            else float("nan"))

    q_ap = ap.get("ap_q", float("nan"))
    dq_ap = (float(q_ap) - host.q
             if isinstance(q_ap, (int, float)) and math.isfinite(float(q_ap))
             else float("nan"))

    title = host.frb
    marker = _refit_fix_marker(outdir)
    if marker:
        title = f"{host.frb}  {marker}"

    return [
        # 1 — GALFIT geometry
        [
            title,
            f"$q={host.q:.3f}\\pm{_fmt(host.q_err, '{:.3f}')}$",
            f"PA$={vc.wrap_pa(host.pa):.1f}^\\circ$",
            f"$R_e={host.re:.1f}\\,$px $({host.re_arcsec:.2f}'')$",
        ],
        # 2 — n + photometry (+ mag uncertainties)
        [
            f"$n={host.n:.2f}$",
            f"$m={_fmt(m_gal, '{:.2f}')}$",
            f"$\\Delta m_{{\\rm ref-m}}={_fmt(dmag, '{:+.2f}')}$",
            f"$\\delta m_{{\\rm gal}}={_fmt(mag.get('mag_galfit_err'), '{:.3f}')}$",
            f"$\\delta m_{{\\rm ref}}={_fmt(mag.get('ref_mag_err'), '{:.3f}')}$",
        ],
        # 3 — fit quality + Lotz+2008 Gini–M20 class (most of the 53 cut = late)
        [
            f"$\\chi^2/\\nu_{{\\rm corr}}={_fmt(gcorr, '{:.2f}')}$",
            f"$\\chi^2/\\nu|_{{2R_e,{{\\rm corr}}}}={_fmt(lcorr, '{:.2f}')}$",
            (
                f"${{\\rm RFF}}_{{1R_e}}={_fmt(rff.get('rff_1re'), '{:+.3f}')}$"
            ),
            (
                f"${{\\rm RFF}}_{{2R_e}}={_fmt(rff.get('rff_2re'), '{:+.3f}')}"
                f"\\pm{_fmt(rff.get('rff_2re_err'), '{:.3f}')}$"
                if (isinstance(rff.get("rff_2re_err"), (int, float))
                    and math.isfinite(float(rff.get("rff_2re_err"))))
                else f"${{\\rm RFF}}_{{2R_e}}={_fmt(rff.get('rff_2re'), '{:+.3f}')}$"
            ),
            _lotz_label(host),
        ],
        # 4 — Fourier (δq + m=2 phase slope) + sky + AstroPhot
        [
            f"$\\delta q_{{\\rm Fou}}={_fmt(fou.get('fourier_dq'), '{:+.3f}')}$",
            f"$\\psi_2'={_fmt(fou.get('fourier_m2_phase_slope_deg_per_re'), '{:+.0f}')}"
            f"^\\circ/R_e$",
            # q_sky± = GALFIT q under sky ±1σ (not a free-q fit; q itself is free)
            f"$q_{{\\rm sky+}}={_fmt(skyj.get('q_sky_plus'), '{:.3f}')}$"
            f"   $q_{{\\rm sky-}}={_fmt(skyj.get('q_sky_minus'), '{:.3f}')}$",
            f"$\\Delta q_{{\\rm AP-G}}={_fmt(dq_ap, '{:+.3f}')}$",
        ],
    ]


def _draw_header(fig: plt.Figure, columns: list[list[str]]) -> None:
    """Multi-column metric strip above the plot grid — values only, no labels."""
    n = len(columns)
    for i, lines in enumerate(columns):
        x0 = (i + 0.5) / n
        fig.text(x0, 0.985, "\n".join(lines), ha="center", va="top",
                 fontsize=13, transform=fig.transFigure, linespacing=1.42)
        if i < n - 1:
            xd = (i + 1) / n
            fig.add_artist(plt.Line2D(
                [xd, xd], [0.855, 0.995], transform=fig.transFigure,
                color="0.75", lw=0.8, solid_capstyle="butt",
            ))


def run(host: vc.HostData, outdir: str, panel_path: str | None = None,
        metrics_outdir: str | None = None) -> dict:
    """Build the standard 2×3 verification panel for ``host``.

    ``outdir`` is where check products (``*.npz``, ``*.json``) are read from
    and, by default, where ``panel.png`` is written. Pass ``panel_path`` to
    write elsewhere (e.g. Re-fits). Pass ``metrics_outdir`` to read header /
    profile products from a different directory than the image host.
    """
    sky = host.sky_level if math.isfinite(host.sky_level) else 0.0
    mdir = metrics_outdir if metrics_outdir is not None else outdir
    fou = _load_npz(mdir, "fourier_profiles.npz")
    iso = _load_npz(mdir, "isophote_profiles.npz")

    fig = plt.figure(figsize=(16, 10.0))
    gs = fig.add_gridspec(2, 3, hspace=0.30, wspace=0.22,
                          top=0.80, bottom=0.06, left=0.05, right=0.98)

    # Data/model: asinh on sky-subtracted flux (shared 1–99% window from data).
    # Residual: sky-subtracted / σ, RdBu_r, ±5 (unchanged).
    data_flux = np.asarray(host.data, float) - sky
    model_flux = np.asarray(host.model, float) - sky
    data_flux = np.where(host.mask, np.nan, data_flux)
    model_flux = np.where(host.mask, np.nan, model_flux)
    finite = np.isfinite(data_flux)
    if finite.any():
        lo, hi = np.nanpercentile(data_flux[finite], [1.0, 99.0])
    else:
        lo, hi = 0.0, 1.0
    data_disp = _asinh_display(data_flux, float(lo), float(hi))
    model_disp = _asinh_display(model_flux, float(lo), float(hi))
    with np.errstate(divide="ignore", invalid="ignore"):
        resid_s = np.where(host.sigma > 0, host.resid / host.sigma, np.nan)
    resid_ignore = host.metric_mask if host.metric_mask is not None else host.mask
    resid_s = np.where(resid_ignore, np.nan, resid_s)

    for col, (img, title, cmap, vmin, vmax) in enumerate((
        (data_disp, r"data  (asinh, 1–99%)", "gray", 0.0, 1.0),
        (model_disp, r"model  (asinh, 1–99%)", "gray", 0.0, 1.0),
        (resid_s, r"residual / $\sigma$  (clipped $\pm5$)", "RdBu_r",
         -RES_CLIP, RES_CLIP),
    )):
        ax = fig.add_subplot(gs[0, col])
        im = ax.imshow(img, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax,
                       interpolation="nearest")
        _ellipse_overlay(ax, host, 1.0, color="k", lw=1.0)
        _ellipse_overlay(ax, host, 2.0, color="k", lw=0.7, ls="--")
        ax.set_title(title, fontsize=10)
        ax.set_xticks([]); ax.set_yticks([])
        if col == 2:
            fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)

    # Major- and minor-axis surface brightness vs elliptical a/Re.
    ax = fig.add_subplot(gs[1, 0])
    a_grid = np.linspace(0.4, 3.0, 27)
    mu_maj_d, e_maj_d = _axis_cut_mu(host.data, host.sigma, host, sky, a_grid,
                                     minor=False)
    mu_min_d, e_min_d = _axis_cut_mu(host.data, host.sigma, host, sky, a_grid,
                                     minor=True)
    mu_maj_m, _ = _axis_cut_mu(host.model, host.sigma, host, sky, a_grid,
                               minor=False)
    mu_min_m, _ = _axis_cut_mu(host.model, host.sigma, host, sky, a_grid,
                               minor=True)
    ok_maj = np.isfinite(mu_maj_d)
    ok_min = np.isfinite(mu_min_d)
    if np.any(ok_maj) or np.any(ok_min):
        if np.any(ok_maj):
            ax.errorbar(a_grid[ok_maj], mu_maj_d[ok_maj], yerr=e_maj_d[ok_maj],
                        fmt="o", ms=3, lw=0.8, color="k", label="data major")
            ax.plot(a_grid, mu_maj_m, "-", color="tab:red", lw=1.2,
                    label="model major")
        if np.any(ok_min):
            ax.errorbar(a_grid[ok_min], mu_min_d[ok_min], yerr=e_min_d[ok_min],
                        fmt="s", ms=3, lw=0.8, color="0.35",
                        label="data minor")
            ax.plot(a_grid, mu_min_m, "--", color="tab:orange", lw=1.2,
                    label="model minor")
        ax.invert_yaxis()
        ax.set_ylabel(r"$\mu$  [mag arcsec$^{-2}$]")
        ax.legend(fontsize=7, frameon=False)
    else:
        ax.text(0.5, 0.5, "no axis profile", ha="center", va="center",
                transform=ax.transAxes, fontsize=9, color="0.5")
    ax.set_xlabel(r"$a / R_e$  (elliptical)")
    ax.set_title("major- & minor-axis profiles", fontsize=10)

    # Isophotal q(a) only — PA twin removed (triage does not need it).
    ax = fig.add_subplot(gs[1, 1])
    if iso is not None:
        ax.errorbar(iso["a_re"], iso["q_data"], yerr=iso["q_data_err"], fmt="o",
                    ms=3, lw=0.8, color="k", label=r"$q$ data")
        ax.plot(iso["a_re"], iso["q_model"], "-", color="tab:red",
                label=r"$q$ model")
        ax.axhline(host.q, color="tab:blue", ls=":", lw=1.2,
                   label=r"GALFIT $q$ (intrinsic)")
        ax.set_ylim(0, 1)
        ax.legend(fontsize=8, frameon=False, loc="lower right")
    else:
        ax.text(0.5, 0.5, "no isophote profile", ha="center", va="center",
                transform=ax.transAxes, fontsize=9, color="0.5")
    ax.set_xlabel(r"$a / R_e$")
    ax.set_ylabel("$q$")
    ax.set_title("isophotal $q(a)$", fontsize=10)

    # dq(a) from the Fourier estimator: what the m=2 residual implies.
    ax = fig.add_subplot(gs[1, 2])
    if fou is not None:
        x, y, e = fou["a_mid_re"], fou["dq"], fou["dq_err"]
        ok = np.isfinite(x) & np.isfinite(y)
        ax.fill_between(x[ok], (y - e)[ok], (y + e)[ok], color="tab:purple",
                        alpha=0.25, lw=0)
        ax.plot(x[ok], y[ok], "-o", ms=3, lw=1.0, color="tab:purple")
        ax.axhline(0, color="k", lw=0.8)
        ax.axvline(2.0, color="0.6", ls="--", lw=0.8)
        lim = np.nanpercentile(np.abs(y[ok]), 95) if np.any(ok) else 0.1
        if math.isfinite(lim) and lim > 0:
            ax.set_ylim(-2.5 * lim, 2.5 * lim)
    else:
        ax.text(0.5, 0.5, "no Fourier profile", ha="center", va="center",
                transform=ax.transAxes, fontsize=9, color="0.5")
    ax.set_xlabel(r"$a / R_e$")
    ax.set_ylabel(r"$\delta q$")
    ax.set_title(r"Fourier $\delta q(a)$ from the $m{=}2$ residual", fontsize=10)

    _draw_header(fig, _header_columns(host, mdir))

    path = panel_path or os.path.join(outdir, "panel.png")
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    fig.savefig(path, dpi=110, bbox_inches="tight")
    plt.close(fig)
    return {"panel_png": os.path.relpath(path, vc.VER_DIR),
            "had_fourier": fou is not None, "had_isophotes": iso is not None}


def _stack_panels(top_path: str, bottom_path: str, out_path: str) -> str:
    """Stack two existing panel PNGs (production on top, refit below)."""
    top = plt.imread(top_path)
    bot = plt.imread(bottom_path)
    # match widths
    if top.shape[1] != bot.shape[1]:
        from PIL import Image
        tw = max(top.shape[1], bot.shape[1])
        def _resize(arr, w):
            im = Image.fromarray((np.clip(arr, 0, 1) * 255).astype(np.uint8)
                                 if arr.dtype != np.uint8 else arr)
            if arr.shape[1] == w:
                return np.asarray(im)
            h = int(round(arr.shape[0] * (w / arr.shape[1])))
            return np.asarray(im.resize((w, h), Image.Resampling.LANCZOS))
        top_u = _resize(top, tw)
        bot_u = _resize(bot, tw)
    else:
        def _to_u8(arr):
            if arr.dtype == np.uint8:
                return arr
            return (np.clip(arr, 0, 1) * 255).astype(np.uint8)
        top_u, bot_u = _to_u8(top), _to_u8(bot)
        # pad channel axis if needed
    if top_u.shape[-1] != bot_u.shape[-1]:
        # force RGB
        def _rgb(a):
            if a.shape[-1] == 4:
                return a[..., :3]
            return a
        top_u, bot_u = _rgb(top_u), _rgb(bot_u)
    gap = np.full((12, top_u.shape[1], top_u.shape[-1]), 255, dtype=np.uint8)
    stacked = np.concatenate([top_u[..., :3], gap[..., :3], bot_u[..., :3]], axis=0)
    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    plt.imsave(out_path, stacked)
    return out_path


def write_refit_panel(
    wkdir: str,
    out_path: str | None = None,
    *,
    frb: str | None = None,
    label: str = "refit",
    compare_production: bool = True,
) -> dict:
    """Write the standard verification panel for a staged refit.

    Also copies the production panel alongside (and stacks both into one PNG
    when the production panel exists). ``label`` is unused in the figure
    itself — kept for callers / filenames.
    """
    del label  # filename/label handled by caller; panel format matches production
    host = vc.load_host_from_dir(wkdir, frb=frb)
    frb = host.frb
    refit_root = os.path.dirname(os.path.abspath(wkdir))
    # If wkdir is .../Re-fits/FRB/galfit_sky_*, parent is the FRB folder.
    if os.path.basename(refit_root).startswith("galfit_"):
        refit_root = os.path.dirname(refit_root)
    if out_path is None:
        out_path = os.path.join(refit_root, "panel_refit.png")

    # Same panel format as verification; header/profiles from this fit only
    # (no production check JSONs — those describe the old model).
    result = run(host, outdir=wkdir, panel_path=out_path, metrics_outdir=wkdir)

    out = {
        "panel_refit_png": out_path,
        "panel_png": result["panel_png"],
    }
    if compare_production:
        prod_panel = os.path.join(vc.per_host_dir(frb, create=False), "panel.png")
        if not os.path.isfile(prod_panel):
            alt = os.path.join(vc.OUT_ROOT, "panels", f"{frb}.png")
            prod_panel = alt if os.path.isfile(alt) else prod_panel
        if os.path.isfile(prod_panel):
            import shutil
            prod_copy = os.path.join(os.path.dirname(out_path), "panel_production.png")
            shutil.copy2(prod_panel, prod_copy)
            out["panel_production_png"] = prod_copy
            stacked = os.path.join(os.path.dirname(out_path), "panel_production_and_refit.png")
            try:
                _stack_panels(prod_copy, out_path, stacked)
                out["panel_stacked_png"] = stacked
            except Exception as exc:
                out["stack_error"] = f"{type(exc).__name__}: {exc}"
    return out


def _confirmation_notes() -> dict[str, str]:
    out: dict[str, str] = {}
    if not os.path.isfile(_HC_CSV):
        return out
    with open(_HC_CSV, encoding="utf-8", newline="") as f:
        for row in csv.DictReader(f):
            frb = str(row.get("frb", "")).strip()
            notes = row.get("notes") or ""
            if frb:
                out[frb] = notes
    return out


def residual_workdir(frb: str, notes: str = "") -> tuple[str, str]:
    """Confirmed-leg (or PSF-only) workdir for residual tiles; else production."""
    refits = os.path.join(_REFITS, frb)
    if "star" in notes.lower():
        psf_only = os.path.join(refits, "psf_only")
        if os.path.isfile(os.path.join(psf_only, "out.fits")):
            return psf_only, "star"
    m = _ALT_PANEL.search(notes)
    if m:
        leg = m.group(1).split("_", 1)[1]
        cands = [os.path.join(refits, leg)]
        if leg == "psf":
            cands.append(os.path.join(refits, "sandbox"))
        for cand in cands:
            if os.path.isfile(os.path.join(cand, "out.fits")):
                return cand, leg
    return vc.host_dir(frb), "production"


def contact_sheet(frbs: list[str], out_path: str, ncols: int = 8) -> str:
    """Cohort triage sheet: residual in sigma units, one tile per host.

    Uses the confirmed-leg ``out.fits`` when ``host_confirmation.csv`` cites an
    alternate panel (or a star/PSF-only reject). Production otherwise.
    """
    frbs = list(frbs)
    notes_map = _confirmation_notes()
    nrows = max(1, math.ceil(len(frbs) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(2.0 * ncols, 2.15 * nrows))
    axes = np.atleast_1d(axes).ravel()
    for ax in axes:
        ax.set_xticks([]); ax.set_yticks([]); ax.axis("off")
    for ax, frb in zip(axes, frbs):
        ax.axis("on"); ax.set_xticks([]); ax.set_yticks([])
        notes = notes_map.get(frb, "")
        hdir, tag = residual_workdir(frb, notes)
        try:
            host = vc.load_host_from_dir(hdir, frb=frb)
        except Exception:
            ax.text(0.5, 0.5, f"{frb}\nunreadable", ha="center", va="center",
                    fontsize=6, transform=ax.transAxes)
            continue
        with np.errstate(divide="ignore", invalid="ignore"):
            snr = np.where(host.sigma > 0, host.resid / host.sigma, np.nan)
        snr = np.where(host.mask, np.nan, snr)
        re = host.re if math.isfinite(host.re) and host.re > 0 else 5.0
        half = max(int(3 * re), 15)
        y0, y1 = int(max(0, host.yc - half)), int(min(host.shape[0], host.yc + half))
        x0, x1 = int(max(0, host.xc - half)), int(min(host.shape[1], host.xc + half))
        ax.imshow(snr[y0:y1, x0:x1], origin="lower", cmap="RdBu_r",
                  vmin=-RES_CLIP, vmax=RES_CLIP)
        if tag == "star" or not math.isfinite(host.q):
            label = f"{frb} [{tag}]\nstar"
        elif tag != "production":
            label = f"{frb} [{tag}]\n$q$={host.q:.2f}"
        else:
            label = f"{frb}\n$q$={host.q:.2f}"
        ax.set_title(label, fontsize=6.5, pad=2)
    fig.suptitle(
        r"Residual / $\sigma$, clipped $\pm5$, cropped to $3R_e$"
        "  (confirmed-leg / PSF-only where gated)",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    fig.savefig(out_path, dpi=100)
    plt.close(fig)
    return out_path
