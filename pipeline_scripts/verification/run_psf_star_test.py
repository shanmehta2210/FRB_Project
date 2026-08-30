"""Single-PSF (star) vs Sérsic test for unresolved hosts.

Copies the confirmed-leg stamp/PSF/sky setup, replaces the host Sérsic with
one GALFIT ``psf`` component, and compares χ². Does not touch production
``Output/`` or confirmed panels.
"""
from __future__ import annotations

import contextlib
import math
import os
import re
import shutil
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits

VER = os.path.dirname(os.path.abspath(__file__))
if VER not in sys.path:
    sys.path.insert(0, VER)

import vercommon as vc  # noqa: E402
from checks import sky_perturb as sp  # noqa: E402
from run_sandbox import make_panel  # noqa: E402

REFITS = os.path.join(VER, "Re-fits")
OUT_CSV = os.path.join(VER, "outputs", "tables", "psf_star_test.csv")
OUT_PNG = os.path.join(VER, "outputs", "plots", "psf_star_test.png")

# Hosts originally noted as unresolved / star-like, later confirmed.
# Parent = confirmed-leg directory used for the science panel.
TARGETS: list[tuple[str, str]] = [
    ("20190711A", os.path.join(REFITS, "20190711A", "n1_sky")),
    ("20220725A", os.path.join(REFITS, "20220725A", "n1_sky")),
    ("20230526A", os.path.join(REFITS, "20230526A", "n1_sky")),
    ("20230626A", os.path.join(REFITS, "20230626A", "n1_sky")),
    ("20231220A", os.path.join(REFITS, "20231220A", "n1_sky")),
    ("20240229A", os.path.join(vc.OUTPUT_ROOT, "20240229A_all")),
]

_INPUTS = ("host_cutout.fits", "host_sigma.fits", "host_mask.fits", "proto_image.fits")
_CHI2 = re.compile(
    r"Chi\^2/nu\s*=\s*([0-9.eE+-]+).*?Chi\^2\s*=\s*([0-9.eE+-]+).*?Ndof\s*=\s*([0-9]+)"
)
_XY = re.compile(r"^\s*1\)\s+(\S+)\s+(\S+)", re.M)
_MAG = re.compile(r"^\s*3\)\s+(\S+)", re.M)
_RE = re.compile(r"^\s*4\)\s+(\S+)", re.M)
_N = re.compile(r"^\s*5\)\s+(\S+)", re.M)
_Q = re.compile(r"^\s*9\)\s+(\S+)", re.M)


def _first_sersic_block(text: str) -> str:
    parts = re.split(r"(?=# Component number:)", text)
    for p in parts:
        if re.search(r"^\s*0\)\s*sersic", p, re.M | re.I):
            return p
    raise ValueError("no sersic component")


def _parse_chi2(path: str) -> tuple[float, float, int]:
    text = open(path, encoding="utf-8", errors="replace").read()
    matches = list(_CHI2.finditer(text))
    if not matches:
        return float("nan"), float("nan"), 0
    m = matches[-1]
    return float(m.group(1)), float(m.group(2)), int(m.group(3))


def _host_params(g01: str) -> dict:
    block = _first_sersic_block(open(g01, encoding="utf-8", errors="replace").read())
    xy = _XY.search(block)
    mag = _MAG.search(block)
    re_px = _RE.search(block)
    n = _N.search(block)
    q = _Q.search(block)
    return {
        "x": float(xy.group(1)) if xy else float("nan"),
        "y": float(xy.group(2)) if xy else float("nan"),
        "mag": float(mag.group(1)) if mag else float("nan"),
        "re": float(re_px.group(1)) if re_px else float("nan"),
        "n": float(n.group(1)) if n else float("nan"),
        "q": float(q.group(1)) if q else float("nan"),
    }


def _psf_feedme_from_parent(parent_feedme: str, x: float, y: float, mag: float) -> str:
    """Keep control + sky; replace host Sérsic with a free PSF."""
    lines = open(parent_feedme, encoding="utf-8", errors="replace").read().splitlines()
    out: list[str] = []
    in_host = False
    host_done = False
    sky_buf: list[str] = []
    in_sky = False
    header: list[str] = []
    after_comps = False

    i = 0
    while i < len(lines):
        line = lines[i]
        if re.match(r"\s*#\s*Component number:\s*1\b", line) and not host_done:
            in_host = True
            i += 1
            continue
        if in_host:
            if re.match(r"\s*#\s*Component number:", line) or line.startswith("===="):
                in_host = False
                host_done = True
                continue  # reprocess this line
            i += 1
            continue
        if re.match(r"\s*#\s*Component number:", line) and "sky" in "".join(lines[i:i + 3]).lower():
            in_sky = True
        if in_sky:
            sky_buf.append(line)
            if line.startswith("===="):
                in_sky = False
                after_comps = True
            i += 1
            continue
        if after_comps:
            i += 1
            continue
        header.append(line)
        i += 1

    # Empty constraints — PSF has no Re/n to bound.
    header_txt = "\n".join(header)
    header_txt = re.sub(
        r"^G\).*$",
        "G) none  # no constraints (PSF: x,y,mag only)",
        header_txt,
        flags=re.M,
    )

    psf_block = f"""
# Component number: 1
 0) psf  # Component type
 1) {x:.4f} {y:.4f} 1 1  # position x y
 3) {mag:.4f} 1  # Integrated magnitude
 Z) 0  # Skip this model
"""
    sky = "\n".join(sky_buf).rstrip("=").rstrip()
    return header_txt.rstrip() + "\n" + psf_block + "\n" + sky + "\n" + ("=" * 80) + "\n"


def _core_chi2(out_fits: str, sigma_fits: str, mask_fits: str | None,
               xc: float, yc: float, radius: float) -> dict:
    with fits.open(out_fits) as hdul:
        data = np.asarray(hdul[1].data, float)
        model = np.asarray(hdul[2].data, float)
        resid = np.asarray(hdul[3].data, float)
    sigma = np.asarray(fits.getdata(sigma_fits), float)
    mask = np.zeros(data.shape, dtype=bool)
    if mask_fits and os.path.isfile(mask_fits):
        mask = np.asarray(fits.getdata(mask_fits), float) > 0
    yy, xx = np.indices(data.shape)
    # GALFIT 1-based xy → 0-based array
    rr = np.hypot(xx - (xc - 1.0), yy - (yc - 1.0))
    sel = (rr <= radius) & np.isfinite(resid) & np.isfinite(sigma) & (sigma > 0) & ~mask
    n = int(np.count_nonzero(sel))
    if n < 5:
        return {"core_npix": n, "core_chi2nu": float("nan"), "core_rms": float("nan")}
    chi2 = float(np.sum((resid[sel] / sigma[sel]) ** 2))
    return {
        "core_npix": n,
        "core_chi2nu": chi2 / n,
        "core_rms": float(np.sqrt(np.mean(resid[sel] ** 2))),
    }


def _asinh(img: np.ndarray, lo: float, hi: float) -> np.ndarray:
    x = np.clip((img - lo) / (hi - lo + 1e-30), 0.0, 1.0)
    return np.arcsinh(10.0 * x) / np.arcsinh(10.0)


def make_compare_panel(frb: str, sersic_out: str, psf_dir: str, dest: str) -> None:
    def load(path: str):
        with fits.open(path) as hdul:
            return (np.asarray(hdul[1].data, float),
                    np.asarray(hdul[2].data, float),
                    np.asarray(hdul[3].data, float))

    d_s, m_s, r_s = load(sersic_out)
    d_p, m_p, r_p = load(os.path.join(psf_dir, "out.fits"))
    finite = np.isfinite(d_s)
    lo, hi = np.nanpercentile(d_s[finite], [1, 99]) if finite.any() else (0.0, 1.0)
    rlim = np.nanpercentile(np.abs(r_s[np.isfinite(r_s)]), 98)
    rlim = max(float(rlim), 1e-8)

    fig, axes = plt.subplots(2, 3, figsize=(10.5, 7.0), constrained_layout=True)
    rows = [
        (d_s, m_s, r_s, "Sérsic (confirmed leg)"),
        (d_p, m_p, r_p, "PSF only (star)"),
    ]
    for row, (data, model, resid, label) in enumerate(rows):
        for col, (arr, title, is_res) in enumerate((
            (data, "data", False),
            (model, "model", False),
            (resid, "resid", True),
        )):
            ax = axes[row, col]
            if is_res:
                ax.imshow(arr, origin="lower", cmap="RdBu_r", vmin=-rlim, vmax=rlim,
                          interpolation="nearest")
            else:
                ax.imshow(_asinh(arr, lo, hi), origin="lower", cmap="gray",
                          vmin=0, vmax=1, interpolation="nearest")
            ax.set_title(f"{label} — {title}" if col == 0 else title, fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
    fig.suptitle(frb, fontsize=13)
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def setup_and_run(frb: str, parent: str) -> dict:
    dest = os.path.join(REFITS, frb, "psf_only")
    os.makedirs(dest, exist_ok=True)
    for name in _INPUTS:
        src = os.path.join(parent, name)
        if not os.path.isfile(src):
            raise FileNotFoundError(src)
        shutil.copy2(src, os.path.join(dest, name))

    g01 = os.path.join(parent, "galfit.01")
    feedme_src = os.path.join(parent, "galfit.feedme")
    if not os.path.isfile(g01):
        raise FileNotFoundError(g01)
    hp = _host_params(g01)
    chi2nu_s, chi2_s, ndof_s = _parse_chi2(g01)
    text = _psf_feedme_from_parent(feedme_src, hp["x"], hp["y"], hp["mag"])
    open(os.path.join(dest, "galfit.feedme"), "w", encoding="utf-8").write(text)

    for stale in ("fit.log", "out.fits"):
        p = os.path.join(dest, stale)
        if os.path.isfile(p):
            os.remove(p)

    run_galfit = sp._run_galfit()
    log_path = os.path.join(dest, "galfit_stdout.log")
    with open(log_path, "w", encoding="utf-8") as fh, contextlib.redirect_stdout(fh):
        ok = bool(run_galfit(dest))

    g01_p = os.path.join(dest, "galfit.01")
    chi2nu_p = chi2_p = float("nan")
    ndof_p = 0
    mag_p = float("nan")
    if os.path.isfile(g01_p):
        chi2nu_p, chi2_p, ndof_p = _parse_chi2(g01_p)
        try:
            mag_p = _host_params(g01_p.replace("galfit.01", "galfit.01"))["mag"]
        except Exception:
            mag_p = float("nan")
        # PSF block uses same mag regex on first component
        block = re.split(r"(?=# Component number:)",
                         open(g01_p, encoding="utf-8", errors="replace").read())
        for b in block:
            if re.search(r"^\s*0\)\s*psf", b, re.M | re.I):
                mm = _MAG.search(b)
                xy = _XY.search(b)
                mag_p = float(mm.group(1)) if mm else mag_p
                if xy:
                    hp["psf_x"] = float(xy.group(1))
                    hp["psf_y"] = float(xy.group(2))
                break

    core_s = _core_chi2(
        os.path.join(parent, "out.fits"),
        os.path.join(parent, "host_sigma.fits"),
        os.path.join(parent, "host_mask.fits"),
        hp["x"], hp["y"], radius=8.0,
    )
    core_p = {"core_npix": 0, "core_chi2nu": float("nan"), "core_rms": float("nan")}
    if os.path.isfile(os.path.join(dest, "out.fits")):
        core_p = _core_chi2(
            os.path.join(dest, "out.fits"),
            os.path.join(dest, "host_sigma.fits"),
            os.path.join(dest, "host_mask.fits"),
            hp.get("psf_x", hp["x"]), hp.get("psf_y", hp["y"]), radius=8.0,
        )
        make_panel(dest)
        make_compare_panel(
            frb,
            os.path.join(parent, "out.fits"),
            dest,
            os.path.join(dest, "compare_sersic_vs_psf.png"),
        )

    # Nested-model: Sérsic n=1 has 3 extra free params (Re, q, PA) vs PSF.
    dchi2 = chi2_p - chi2_s
    # Prefer PSF if it is not worse by more than the extra 3 params (AIC).
    # ΔAIC = (χ²_psf + 2*k_psf) - (χ²_sersic + 2*k_sersic); k_sersic = k_psf + 3
    # → ΔAIC = Δχ² - 6. Negative ΔAIC ⇒ PSF preferred.
    delta_aic = dchi2 - 6.0
    if not math.isfinite(dchi2):
        verdict = "failed"
    elif delta_aic < -2:
        verdict = "PSF preferred"
    elif abs(delta_aic) <= 2:
        verdict = "indistinguishable"
    else:
        verdict = "Sersic preferred"

    row = {
        "frb": frb,
        "parent": os.path.relpath(parent, VER),
        "galfit_ok": ok,
        "sersic_mag": hp["mag"],
        "sersic_re": hp["re"],
        "sersic_n": hp["n"],
        "sersic_q": hp["q"],
        "sersic_chi2nu": chi2nu_s,
        "sersic_chi2": chi2_s,
        "sersic_ndof": ndof_s,
        "psf_mag": mag_p,
        "psf_chi2nu": chi2nu_p,
        "psf_chi2": chi2_p,
        "psf_ndof": ndof_p,
        "delta_chi2_psf_minus_sersic": dchi2,
        "delta_aic_psf_minus_sersic": delta_aic,
        "sersic_core_chi2nu": core_s["core_chi2nu"],
        "psf_core_chi2nu": core_p["core_chi2nu"],
        "verdict": verdict,
    }
    print(
        f"[{frb}] {verdict}  chi2/nu Sersic={chi2nu_s:.3f}  PSF={chi2nu_p:.3f}  "
        f"dchi2={dchi2:+.1f}  dAIC={delta_aic:+.1f}  "
        f"core chi2/nu {core_s['core_chi2nu']:.2f}->{core_p['core_chi2nu']:.2f}"
    )
    return row


def _summary_figure(df: pd.DataFrame, path: str) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 4.6), constrained_layout=True)
    x = np.arange(len(df))
    w = 0.35
    ax.bar(x - w / 2, df["sersic_chi2nu"], w, label="Sérsic χ²/ν", color="#4a4a4a")
    ax.bar(x + w / 2, df["psf_chi2nu"], w, label="PSF χ²/ν", color="#b33")
    ax.set_xticks(x)
    ax.set_xticklabels(df["frb"], rotation=30, ha="right")
    ax.set_ylabel(r"GALFIT $\chi^2/\nu$")
    ax.set_title("Unresolved hosts: Sérsic vs single PSF")
    ax.legend(frameon=False)
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main() -> int:
    rows = []
    rc = 0
    for frb, parent in TARGETS:
        try:
            rows.append(setup_and_run(frb, parent))
        except Exception as exc:
            print(f"[{frb}] ERROR: {exc}", file=sys.stderr)
            rows.append({"frb": frb, "galfit_ok": False, "verdict": f"error: {exc}"})
            rc = 1
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT_CSV, index=False)
    if df["psf_chi2nu"].notna().any():
        _summary_figure(df.dropna(subset=["psf_chi2nu"]), OUT_PNG)
    print(f"\nwrote {OUT_CSV}")
    print(df[["frb", "sersic_chi2nu", "psf_chi2nu", "delta_chi2_psf_minus_sersic",
              "delta_aic_psf_minus_sersic", "verdict"]].to_string(index=False))
    return rc


if __name__ == "__main__":
    raise SystemExit(main())
