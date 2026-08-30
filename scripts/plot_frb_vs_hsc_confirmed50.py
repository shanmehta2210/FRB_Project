"""
Confirmed-50 FRB host cos(i) CDF vs HSC disk-analogue Monte Carlo.

FRB sample: in_53 ∩ confirmed (winning-leg q, same panel rule as the PPT).
HSC: Kawinwanichakij EXP analogue, photometric rmag ≤ 22, b/a > q0.
Each of 10_000 draws samples 50 HSC galaxies without replacement.

Protocol A: σ_q = 10 × q_err, then i = Hubble(q) + N(0, 5°)
    (existing inflated-error protocol).
Protocol B: σ_q = sqrt(q_err² + σ_q,sky² + σ_q,5°²) with
    σ_q,sky = |q_sky+ − q_sky−| / 2 and σ_q,5° the Hubble half-range
    at i ± 5°. No extra 5° after Hubble.

Run from repo root::

    python scripts/plot_frb_vs_hsc_confirmed50.py
"""
from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import interpolate

_SCRIPTS = Path(__file__).resolve().parent
REPO_ROOT = _SCRIPTS.parent
VER_DIR = REPO_ROOT / "pipeline_scripts" / "verification"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))
if str(VER_DIR) not in sys.path:
    sys.path.insert(0, str(VER_DIR))

import vercommon as vc  # noqa: E402
from null_catalog_utils import (  # noqa: E402
    HSC_KAWIN_EXP_DEFAULT,
    Q0,
    hubble_cosi_from_ba,
)
from pipeline_null_plot_utils import (  # noqa: E402
    PLOTS_NULL,
    add_inclination_top_axis,
    default_font,
    save_figure,
)
from plot_frb_inflated_error_cdf import (  # noqa: E402
    INC_FLOOR_DEG,
    N_DRAWS,
    SEED,
    X,
    _hubble_inc_deg,
    mc_inflated_cdf_envelope,
)
from plot_frb_vs_sdss_inflated_cdf import mc_subsample_envelope  # noqa: E402

OUT_DIR = PLOTS_NULL / "v2" / "frb_vs_hsc_confirmed50"
HC_CSV = VER_DIR / "host_confirmation.csv"
REFITS = VER_DIR / "Re-fits"
PER_HOST = Path(vc.OUT_ROOT) / "per_host"

_ALT = re.compile(
    r"outputs/panels/([A-Za-z0-9]+_(?:n1_sky|n1|sky|psf))\.png"
)
HSC_COL = "#4daf4a"
FRB_COL = "#d62728"


def _finite(x, default: float = float("nan")) -> float:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return default
    return v if math.isfinite(v) else default


def _leg_from_notes(notes: str) -> str:
    m = _ALT.search(notes or "")
    if not m:
        return "production"
    return m.group(1).split("_", 1)[1]


def _workdir(frb: str, leg: str) -> Path:
    if leg == "production":
        return Path(vc.host_dir(frb))
    candidates = [
        REFITS / frb / leg,
        REFITS / frb / "sandbox" if leg == "psf" else None,
        REFITS / frb / "sandbox" / leg,
    ]
    for path in candidates:
        if path is not None and (path / "out.fits").is_file():
            return path
    return Path(vc.host_dir(frb))


def _q_from_outfits(wkdir: Path) -> tuple[float, float]:
    hdr = vc.parse_out_header(str(wkdir / "out.fits"))
    host = (hdr.get("components") or [None])[0]
    if not host:
        return float("nan"), float("nan")
    return _finite(host.get("q")), _finite(host.get("q_err"), 0.0)


def _q_from_summary(wkdir: Path) -> tuple[float, float]:
    path = wkdir / "refit_summary.json"
    if not path.is_file():
        return float("nan"), float("nan")
    data = vc.read_json(str(path))
    fit = data.get("fit") or {}
    return _finite(fit.get("q")), _finite(fit.get("q_err"), 0.0)


def _sky_pm(wkdir: Path, frb: str) -> tuple[float, float, str]:
    for path in (wkdir / "sky.json", PER_HOST / frb / "sky.json"):
        if not path.is_file():
            continue
        data = vc.read_json(str(path))
        qp = _finite(data.get("q_sky_plus"))
        qm = _finite(data.get("q_sky_minus"))
        if math.isfinite(qp) and math.isfinite(qm):
            return qp, qm, str(path.relative_to(REPO_ROOT)).replace("\\", "/")
    return float("nan"), float("nan"), ""


def q_from_inc_deg(inc_deg: np.ndarray | float, q0: float = Q0) -> np.ndarray:
    inc = np.clip(np.asarray(inc_deg, dtype=float), 0.0, 90.0)
    c2 = np.cos(np.radians(inc)) ** 2
    return np.sqrt(q0**2 + (1.0 - q0**2) * c2)


def sigma_q_from_floor_deg(q: np.ndarray, floor_deg: float, q0: float = Q0) -> np.ndarray:
    inc = _hubble_inc_deg(np.asarray(q, dtype=float), q0=q0)
    q_hi = q_from_inc_deg(inc + floor_deg, q0=q0)
    q_lo = q_from_inc_deg(inc - floor_deg, q0=q0)
    return 0.5 * np.abs(q_hi - q_lo)


def sigma_q_protocol_b(
    q: np.ndarray,
    q_err: np.ndarray,
    sigma_q_sky: np.ndarray,
    *,
    floor_deg: float = INC_FLOOR_DEG,
    q0: float = Q0,
) -> np.ndarray:
    qe = np.nan_to_num(np.asarray(q_err, dtype=float), nan=0.0)
    qs = np.nan_to_num(np.asarray(sigma_q_sky, dtype=float), nan=0.0)
    q5 = sigma_q_from_floor_deg(q, floor_deg, q0=q0)
    return np.sqrt(np.maximum(qe, 0.0) ** 2 + np.maximum(qs, 0.0) ** 2 + q5**2)


def mc_sigma_q_envelope(
    ba: np.ndarray,
    sigma_q: np.ndarray,
    *,
    n_draws: int,
    q0: float,
    x: np.ndarray,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Median + 16/84 ECDF from q ~ N(ba, σ_q), Hubble, no extra i floor."""
    ba = np.asarray(ba, dtype=float)
    sigma_q = np.maximum(np.nan_to_num(np.asarray(sigma_q, dtype=float), nan=0.0), 0.0)
    n = len(ba)
    rng = np.random.default_rng(seed)
    y_corners = np.concatenate([[0.0], np.arange(1, n + 1) / n, [1.0]])
    ecdfs = np.empty((n_draws, len(x)), dtype=float)
    for i in range(n_draws):
        q = np.clip(rng.normal(ba, sigma_q), q0 + 1e-6, 1.0)
        inc = np.clip(_hubble_inc_deg(q, q0=q0), 0.0, 90.0)
        cosi = np.sort(np.cos(np.radians(inc)))
        x_corners = np.concatenate([[0.0], cosi, [1.0]])
        for j in range(1, len(x_corners)):
            if x_corners[j] <= x_corners[j - 1]:
                x_corners[j] = x_corners[j - 1] + 1e-12
        f = interpolate.interp1d(
            x_corners, y_corners, kind="linear",
            bounds_error=False, fill_value=(0.0, 1.0),
        )
        ecdfs[i] = f(x)
    return (
        np.median(ecdfs, axis=0),
        np.percentile(ecdfs, 16.0, axis=0),
        np.percentile(ecdfs, 84.0, axis=0),
    )


def confirmed50() -> pd.DataFrame:
    hc = pd.read_csv(HC_CSV, dtype={"frb": str})
    hc["confirmed"] = hc["confirmed"].astype(str).str.lower().eq("true")
    cohort = vc.cohort("all64")
    m = hc.merge(cohort, on="frb", suffixes=("", "_cohort"))
    m = m.loc[m["confirmed"] & m["in_53"]].sort_values("frb").reset_index(drop=True)
    if len(m) != 50:
        raise RuntimeError(f"expected 50 confirmed in-cut hosts, got {len(m)}")

    rows = []
    for _, row in m.iterrows():
        frb = str(row["frb"])
        notes = "" if pd.isna(row.get("notes")) else str(row["notes"])
        leg = _leg_from_notes(notes)
        wkdir = _workdir(frb, leg)
        used_prod = wkdir.resolve() == Path(vc.host_dir(frb)).resolve()
        if used_prod and leg != "production":
            # requested a re-fit that is missing; fall back and record it
            src = "production_fallback"
            q = _finite(row["b_a"])
            q_err = _finite(row["b_a_err"], 0.0)
        elif used_prod:
            src = "production"
            q = _finite(row["b_a"])
            q_err = _finite(row["b_a_err"], 0.0)
        else:
            src = f"refit:{leg}"
            q_s, qe_s = _q_from_summary(wkdir)
            q_o, qe_o = _q_from_outfits(wkdir)
            q = q_s if math.isfinite(q_s) else q_o
            q_err = qe_s if math.isfinite(qe_s) else qe_o
        q_plus, q_minus, sky_path = _sky_pm(wkdir, frb)
        if not (math.isfinite(q_plus) and math.isfinite(q_minus)):
            q_plus, q_minus, sky_path = _sky_pm(Path(vc.host_dir(frb)), frb)
        sky_ok = math.isfinite(q_plus) and math.isfinite(q_minus)
        sig_sky = 0.5 * abs(q_plus - q_minus) if sky_ok else 0.0
        rows.append(
            {
                "frb": frb,
                "leg": "production" if used_prod else leg,
                "q_source": src,
                "q": q,
                "q_err": q_err,
                "q_sky_plus": q_plus,
                "q_sky_minus": q_minus,
                "sigma_q_sky": sig_sky,
                "sky_missing": not sky_ok,
                "sky_path": sky_path,
                "mag": _finite(row.get("mag")),
                "wkdir": str(wkdir.relative_to(REPO_ROOT)).replace("\\", "/"),
            }
        )
    out = pd.DataFrame(rows)
    out["sigma_q_5deg"] = sigma_q_from_floor_deg(out["q"].to_numpy(), INC_FLOOR_DEG)
    out["sigma_q_A"] = 10.0 * np.nan_to_num(out["q_err"].to_numpy(), nan=0.0)
    out["sigma_q_B"] = sigma_q_protocol_b(
        out["q"].to_numpy(),
        out["q_err"].to_numpy(),
        out["sigma_q_sky"].to_numpy(),
    )
    out["cosi_point"] = [
        hubble_cosi_from_ba(float(q), q0=Q0) for q in out["q"].to_numpy()
    ]
    return out


def load_hsc_pool(*, mag_limit: float, q0: float) -> tuple[np.ndarray, int]:
    path = REPO_ROOT / HSC_KAWIN_EXP_DEFAULT
    df = pd.read_csv(path, usecols=["rmag", "ba"])
    rmag = pd.to_numeric(df["rmag"], errors="coerce")
    ba = pd.to_numeric(df["ba"], errors="coerce")
    ok = (
        np.isfinite(rmag) & (rmag <= mag_limit)
        & np.isfinite(ba) & (ba > q0) & (ba <= 1.0)
    )
    cosi = np.array(
        [hubble_cosi_from_ba(float(q), q0=q0) for q in ba[ok].to_numpy()],
        dtype=float,
    )
    return cosi, int(ok.sum())


def _overlay(
    *,
    x: np.ndarray,
    hsc_med, hsc_lo, hsc_hi,
    frb_med, frb_lo, frb_hi,
    n_frb: int,
    n_hsc: int,
    n_draws: int,
    q0: float,
    title: str,
    frb_label: str,
    stem: Path,
) -> None:
    font_prop = default_font()
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.plot((0, 1), (0, 1), "k--", lw=1.2, zorder=1)
    ax.fill_between(x, hsc_lo, hsc_hi, color=HSC_COL, alpha=0.22, linewidth=0, zorder=2)
    ax.plot(
        x, hsc_med, color=HSC_COL, lw=2.0, zorder=3,
        label=rf"HSC disks  ($N={n_frb}\times{n_draws // 1000}$k of {n_hsc:,}, 68% CI)",
    )
    ax.fill_between(x, frb_lo, frb_hi, color=FRB_COL, alpha=0.22, linewidth=0, zorder=4)
    ax.plot(x, frb_med, color=FRB_COL, lw=2.0, zorder=5, label=frb_label)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel(rf"$\cos(i)$ (Hubble, $q_0={q0:g}$)", fontproperties=font_prop, fontsize=11)
    ax.set_ylabel("Cumulative distribution", fontproperties=font_prop, fontsize=11)
    ax.set_title(title, fontsize=11)
    ax.legend(loc="upper left", fontsize=8, frameon=True)
    ax.grid(True, alpha=0.3)
    add_inclination_top_axis(ax, font_prop)
    fig.tight_layout()
    save_figure(fig, stem)
    plt.close(fig)


def write_readme(out_dir: Path, hosts: pd.DataFrame, n_hsc: int, n_draws: int) -> None:
    n = len(hosts)
    n_refit = int((hosts["leg"] != "production").sum())
    n_sky_miss = int(hosts["sky_missing"].sum())
    med_q = float(np.median(hosts["q"]))
    med_cosi = float(np.median(hosts["cosi_point"]))
    med_a = float(np.median(hosts["sigma_q_A"]))
    med_b = float(np.median(hosts["sigma_q_B"]))
    lines = [
        "# Confirmed-50 vs HSC MC cos(i) CDF",
        "",
        f"Paper sample: `in_53` intersect `confirmed` = **{n}**. Winning-leg q "
        f"({n_refit} re-fit, {n - n_refit} production). Hubble q0 = 0.2.",
        "",
        "HSC pool: Kawinwanichakij EXP analogue (`goodfits=1`, 0.4 < n < 1.5), "
        f"photometric r <= 22, b/a > 0.2. **N = {n_hsc:,}**. Each of {n_draws:,} "
        f"draws samples {n} galaxies without replacement.",
        "",
        "| | |",
        "|---|---|",
        f"| median q (FRB) | {med_q:.3f} |",
        f"| median cos(i) (point) | {med_cosi:.3f} |",
        f"| median sigma_q,A (10 x q_err) | {med_a:.3f} |",
        f"| median sigma_q,B (sky xor q_err xor 5 deg) | {med_b:.3f} |",
        f"| hosts missing sky +/- | {n_sky_miss} (those get sigma_q,sky = 0) |",
        "",
        "## Protocol A — `frb_vs_hsc_inflate10.png`",
        "",
        "sigma_q = 10 * q_err, then i = Hubble(q) + N(0, 5 deg).",
        "",
        "## Protocol B — `frb_vs_hsc_sky_quad.png`",
        "",
        "sigma_q,sky = |q+ - q-| / 2",
        "sigma_q,5deg = (1/2) |q(i+5 deg) - q(i-5 deg)|",
        "sigma_q = sqrt(q_err^2 + sigma_q,sky^2 + sigma_q,5deg^2)",
        "",
        "No extra 5 deg after Hubble. q_err is not multiplied by 10.",
        "",
        "Table: `confirmed50_q.csv`.",
        "",
    ]
    dest = out_dir / "README.md"
    dest.write_text("\n".join(lines), encoding="utf-8")
    print(f"[*] Wrote {dest}")


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--n-draws", type=int, default=N_DRAWS)
    p.add_argument("--q0", type=float, default=Q0)
    p.add_argument("--mag-limit", type=float, default=22.0)
    p.add_argument("--floor-deg", type=float, default=INC_FLOOR_DEG)
    p.add_argument("--out-dir", type=Path, default=OUT_DIR)
    args = p.parse_args(argv)

    hosts = confirmed50()
    n_frb = len(hosts)
    n_sky_miss = int(hosts["sky_missing"].sum())
    print(
        f"[*] FRB N={n_frb}  re-fit={int((hosts.leg != 'production').sum())}  "
        f"sky_missing={n_sky_miss}  med_q={float(np.median(hosts.q)):.3f}",
        flush=True,
    )
    if n_sky_miss:
        print("    " + ", ".join(hosts.loc[hosts.sky_missing, "frb"].astype(str)), flush=True)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.out_dir / "confirmed50_q.csv"
    hosts.to_csv(csv_path, index=False, float_format="%.6g")
    print(f"[*] Wrote {csv_path}", flush=True)

    print("[*] Loading HSC EXP analogue ...", flush=True)
    hsc_cosi, n_hsc = load_hsc_pool(mag_limit=args.mag_limit, q0=args.q0)
    if n_hsc < n_frb:
        raise RuntimeError(f"HSC pool N={n_hsc} < n_frb={n_frb}")
    print(
        f"[*] HSC pool N={n_hsc:,}  (rmag<={args.mag_limit:g}, ba>{args.q0:g})",
        flush=True,
    )

    print("[*] HSC subsample envelope ...", flush=True)
    hsc_med, hsc_lo, hsc_hi = mc_subsample_envelope(
        hsc_cosi, n_sample=n_frb, n_draws=args.n_draws, x=X, seed=SEED,
    )

    ba = hosts["q"].to_numpy(dtype=float)
    ba_err = hosts["q_err"].to_numpy(dtype=float)

    print("[*] Protocol A (10 x q_err + 5 deg) ...", flush=True)
    a_med, a_lo, a_hi = mc_inflated_cdf_envelope(
        ba, ba_err, n_draws=args.n_draws, inflate=10.0,
        floor_deg=args.floor_deg, q0=args.q0, x=X, seed=SEED,
    )
    _overlay(
        x=X, hsc_med=hsc_med, hsc_lo=hsc_lo, hsc_hi=hsc_hi,
        frb_med=a_med, frb_lo=a_lo, frb_hi=a_hi,
        n_frb=n_frb, n_hsc=n_hsc, n_draws=args.n_draws, q0=args.q0,
        title=rf"FRB vs HSC  $\cos(i)$  ($m_r\leq{args.mag_limit:g}$, "
              rf"$N={n_frb}$, $10\times\sigma_{{b/a}}+{args.floor_deg:g}^\circ$)",
        frb_label=rf"FRB confirmed  ($N={n_frb}$, $10\times\sigma_{{b/a}}+{args.floor_deg:g}^\circ$)",
        stem=args.out_dir / "frb_vs_hsc_inflate10",
    )

    print("[*] Protocol B (sky + q_err + 5 deg in quadrature) ...", flush=True)
    sig_b = hosts["sigma_q_B"].to_numpy(dtype=float)
    b_med, b_lo, b_hi = mc_sigma_q_envelope(
        ba, sig_b, n_draws=args.n_draws, q0=args.q0, x=X, seed=SEED,
    )
    _overlay(
        x=X, hsc_med=hsc_med, hsc_lo=hsc_lo, hsc_hi=hsc_hi,
        frb_med=b_med, frb_lo=b_lo, frb_hi=b_hi,
        n_frb=n_frb, n_hsc=n_hsc, n_draws=args.n_draws, q0=args.q0,
        title=rf"FRB vs HSC  $\cos(i)$  ($m_r\leq{args.mag_limit:g}$, "
              rf"$N={n_frb}$, sky$\oplus q_{{\rm err}}\oplus{args.floor_deg:g}^\circ$)",
        frb_label=rf"FRB confirmed  ($N={n_frb}$, sky$\oplus q_{{\rm err}}\oplus{args.floor_deg:g}^\circ$)",
        stem=args.out_dir / "frb_vs_hsc_sky_quad",
    )

    write_readme(args.out_dir, hosts, n_hsc, args.n_draws)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
