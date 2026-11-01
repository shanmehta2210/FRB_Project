"""
Replicate the DM_host--inclination analysis of Marnoch et al. (2026, Section 6 and 7.2) on our hosts.

Their definitions (eqs. 3-4): DM_host^Macquart = DM_FRB - DM_NE2001 - 40 - <DM_cosmic>(z) (mean Macquart
relation), rest frame x (1+z). DM variants of their Tables 3-6 except DM_YMW16 (no YMW16 values here).
Pearson and Spearman (p < 0.05), LTSFIT linear fits with 2.6 sigma clipping (Cappellari et al. 2013),
against q, i and cos i, for the full sample, m < 21, n < 2 and with their three bulge-fit hosts removed.
Figure 13 test: RMSE of DM_exgal about <DM_cosmic>(z) before and after subtracting their DM_host(q) model,
on their own Table 1 values and on ours, compared with subtracting a constant DM_host.
Also, for every property of correlations.csv: Fisher 95% CI on r_S, whether it excludes |r_S| = 0.44,
and the power to detect |r_S| = 0.44 at p < 0.05 for that N.

Writes Publication/materials/marnoch_recheck/{replication_correlations.csv, replication_macquart.csv,
property_power.csv, REPLICATION_MARNOCH.md}.
"""
from __future__ import annotations

import contextlib
import io
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_host_property_table as B  # noqa: E402
from ltsfit.lts_linefit import lts_linefit  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
CORR = ROOT / "plots" / "plots_null" / "v2" / "host_correlations"
OUT = ROOT / "Publication" / "materials" / "marnoch_recheck"
Q0, HALO, R_M = 0.2, 40.0, 0.44
BULGE_FITS = ["20190608B", "20220725A", "20231226A"]
MODEL = (-309.3, 0.5, 155.9)  # their eq. 5, rest frame, vs q


def hubble_i(q: np.ndarray) -> np.ndarray:
    return np.degrees(np.arccos(np.sqrt(np.clip((q ** 2 - Q0 ** 2) / (1 - Q0 ** 2), 0, 1))))


def fisher_ci(r: float, n: int) -> tuple[float, float]:
    z, s = np.arctanh(r), 1 / np.sqrt(n - 3)
    return float(np.tanh(z - 1.96 * s)), float(np.tanh(z + 1.96 * s))


def power(r: float, n: int) -> float:
    return float(stats.norm.cdf(abs(np.arctanh(r)) * np.sqrt(n - 3) - 1.96))


def load() -> pd.DataFrame:
    h = pd.read_csv(CORR / "host_props.csv", dtype={"frb": str})
    h = h[np.isfinite(h.z) & np.isfinite(h.dm_obs)].copy()
    mc = np.array([B.mean_dm_cosmic(z) for z in h.z])
    h["dm_frb"] = h.dm_obs
    h["dm_exgal"] = h.dm_obs - h.dm_mw - HALO
    h["dm_frb_minus_cosmic"] = h.dm_obs - mc
    h["dm_mac_obs"] = h.dm_exgal - mc
    h["dm_mac_rest"] = h.dm_mac_obs * (1 + h.z)
    h["i"] = hubble_i(h.q.to_numpy())
    h["cos_i"] = np.cos(np.radians(h.i))
    # NE2001 uncertainty taken as 20% of DM_MW, as the YMW16 difference used by Marnoch et al. is unavailable
    h["sig_dm"] = np.hypot(0.2 * h.dm_mw, 10.0)
    h["sig_q"] = h.sigma_q_B.fillna(h.sigma_q_B.median())
    return h


def lts(x, y, sx, sy, pivot, full=False):
    """LTSFIT; if the intrinsic-scatter root is not bracketed, retry with doubled DM errors."""
    for scale in (1, 2, 4):
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                p = lts_linefit(x, y, sx, sy * scale, clip=2.6, pivot=pivot, plot=False, text=False)
        except ValueError:
            continue
        if full:
            return p, scale
        return p.ab, p.ab_err, p.sig_int, int((~p.mask).sum())
    if full:
        return None, np.nan
    return (np.nan, np.nan), (np.nan, np.nan), np.nan, -1


def corr_row(d: pd.DataFrame, ycol: str, xcol: str, sub: str) -> dict:
    d = d[np.isfinite(d[ycol]) & np.isfinite(d[xcol])]
    x, y = d[xcol].to_numpy(float), d[ycol].to_numpy(float)
    rp, pp = stats.pearsonr(x, y)
    rs, ps = stats.spearmanr(x, y)
    lo, hi = fisher_ci(rs, len(d))
    row = dict(subsample=sub, y=ycol, x=xcol, N=len(d), r_P=rp, p_P=pp, r_S=rs, p_S=ps, rS_lo=lo, rS_hi=hi)
    if xcol in ("q", "cos_i", "i") and ycol.startswith("dm"):
        sq = d.sig_q.to_numpy()
        sx = {"q": sq,
              "cos_i": sq * d.q.to_numpy() / ((1 - Q0 ** 2) * np.maximum(d.cos_i.to_numpy(), 0.05)),
              "i": np.full(len(d), 5.0)}[xcol]
        sy = d.sig_dm.to_numpy() * ((1 + d.z.to_numpy()) if ycol == "dm_mac_rest" else 1)
        piv = round(float(np.median(x)), 1)
        ab, ab_err, sig_int, nclip = lts(x, y, sx, sy, piv)
        row.update(pivot=piv, slope=ab[1], slope_err=ab_err[1], intercept=ab[0], intercept_err=ab_err[0],
                   sig_int=sig_int, n_clipped=nclip)
    return row


def macquart_rmse(dm_host_obs: np.ndarray, q: np.ndarray, z: np.ndarray, label: str, fit=None) -> list[dict]:
    a, piv, b = MODEL
    out = [dict(sample=label, correction="none", N=len(q), rmse=float(np.sqrt(np.mean(dm_host_obs ** 2))))]
    corr = {"Marnoch model DM(q)": (a * (q - piv) + b) / (1 + z),
            "constant 155.9 (their intercept, no q)": np.full_like(q, b) / (1 + z),
            "constant, best-fitting mean (no q)": np.full_like(q, np.mean(dm_host_obs * (1 + z))) / (1 + z)}
    if fit is not None:
        corr["own LTS fit DM(q)"] = (fit[1] * (q - fit[2]) + fit[0]) / (1 + z)
    for k, c in corr.items():
        out.append(dict(sample=label, correction=k, N=len(q), rmse=float(np.sqrt(np.mean((dm_host_obs - c) ** 2)))))
    return out


def main() -> None:
    h = load()
    m = pd.read_csv(OUT / "marnoch2026_table1.csv", dtype={"frb": str})
    shared = sorted(set(m.frb) & set(h.frb))
    subs = {"all": h, "m_r < 21": h[h.m_r < 21], "n < 2": h[h.sersic_n < 2],
            "bulge fits removed": h[~h.frb.isin(BULGE_FITS)], "shared with Marnoch": h[h.frb.isin(shared)]}
    rows = []
    for sub, d in subs.items():
        for y in ("dm_mac_obs", "dm_mac_rest", "dm_exgal", "dm_frb", "dm_frb_minus_cosmic"):
            rows.append(corr_row(d, y, "q", sub))
        rows.append(corr_row(d, "dm_mac_rest", "i", sub))
        rows.append(corr_row(d, "dm_mac_rest", "cos_i", sub))
    for x in ("offset_kpc", "offset_normalized", "u_r_color", "g_r_color"):
        for sub in ("all", "m_r < 21"):
            rows.append(corr_row(subs[sub], "dm_mac_rest", x, sub))
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "replication_correlations.csv", index=False)

    full = R[(R.subsample == "all") & (R.y == "dm_mac_rest") & (R.x == "q")].iloc[0]
    mac = macquart_rmse(m.dm_host_obs.to_numpy(float), m.q.to_numpy(float), m.z.to_numpy(float), "Marnoch Table 1 (34)")
    mac += macquart_rmse(h.dm_mac_obs.to_numpy(float), h.q.to_numpy(float), h.z.to_numpy(float), f"this work ({len(h)})",
                         fit=(full.intercept, full.slope, full.pivot))
    M = pd.DataFrame(mac)
    M.to_csv(OUT / "replication_macquart.csv", index=False)

    c = pd.read_csv(CORR / "correlations.csv")
    P = []
    for r in c.itertuples():
        if not np.isfinite(r.spearman) or r.N < 5:
            continue
        lo, hi = fisher_ci(r.spearman, int(r.N))
        P.append(dict(property=r.property, label=r.label, N=int(r.N), r_S=r.spearman, p=r.spearman_p, lo=lo, hi=hi,
                      excludes_minus044=hi < -R_M or lo > -R_M, excludes_both=(lo > -R_M) and (hi < R_M),
                      power_044=power(R_M, int(r.N))))
    P = pd.DataFrame(P)
    P.to_csv(OUT / "property_power.csv", index=False)

    f = lambda v, n=2: "—" if pd.isna(v) else f"{v:.{n}f}"
    md = ["# Marnoch et al. (2026) analysis repeated on our hosts", "",
          "Source: `scripts/marnoch_replication.py`. DM definitions as their eqs. 3–4 (NE2001, 40 pc cm$^{-3}$ Milky Way halo, "
          "mean Macquart $\\langle$DM$_{\\rm cosmic}\\rangle$); LTSFIT, 2.6σ clip; $i$ from $q$ with $q_0 = 0.2$. "
          f"Hosts with $z$ and DM: {len(h)}. Their bulge-fit hosts removed in the † subsample: {', '.join(BULGE_FITS)}.", "",
          "## DM vs axis ratio, inclination and cos i (their Tables 3–6)", "",
          "| subsample | y | x | N | $r_P$ ($p$) | $r_S$ ($p$) | $r_S$ 95% CI | LTS slope | LTS intercept (pivot) | clipped |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for r in R.itertuples():
        fit = f"{f(r.slope,1)} ± {f(r.slope_err,1)}" if not pd.isna(getattr(r, "slope", np.nan)) else "—"
        icp = f"{f(r.intercept,1)} ({f(r.pivot,1)})" if not pd.isna(getattr(r, "intercept", np.nan)) else "—"
        md.append(f"| {r.subsample} | {r.y} | {r.x} | {r.N} | {f(r.r_P)} ({f(r.p_P,3)}) | {f(r.r_S)} ({f(r.p_S,3)}) | "
                  f"[{f(r.rS_lo)}, {f(r.rS_hi)}] | {fit} | {icp} | {'' if pd.isna(getattr(r,'n_clipped',np.nan)) else int(r.n_clipped)} |")
    md += ["", "Their full-sample relation: slope $-309.3 \\pm 69.5$, intercept $155.9 \\pm 17.4$ at $q = 0.5$, "
           "$r_P = -0.40$, $r_S = -0.44$ (N = 34).", "",
           "## Macquart-relation scatter (their Figure 13)", "",
           "RMSE of DM$_{\\rm exgal}$ about $\\langle$DM$_{\\rm cosmic}\\rangle(z)$, i.e. of observer-frame DM$_{\\rm host}^{\\rm Macquart}$, "
           "before and after subtracting a host term converted to the observer frame.", "",
           "| sample | correction | N | RMSE (pc cm$^{-3}$) |", "|---|---|---|---|"]
    md += [f"| {r.sample} | {r.correction} | {r.N} | {r.rmse:.0f} |" for r in M.itertuples()]
    md += ["", "## Which properties can exclude a correlation as strong as Marnoch's?", "",
           "Fisher 95% CI on $r_S$ with $\\cos i$; power = probability of $p < 0.05$ if the true $|r_S| = 0.44$.", "",
           "| property | N | $r_S$ ($p$) | 95% CI | excludes $-0.44$ | excludes $\\pm 0.44$ | power at $|r_S| = 0.44$ |",
           "|---|---|---|---|---|---|---|"]
    for r in P.sort_values("N", ascending=False).itertuples():
        md.append(f"| {r.label} | {r.N} | {f(r.r_S)} ({f(r.p,3)}) | [{f(r.lo)}, {f(r.hi)}] | {'yes' if r.excludes_minus044 else 'no'} | "
                  f"{'yes' if r.excludes_both else 'no'} | {r.power_044:.2f} |")
    (OUT / "REPLICATION_MARNOCH.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(M.to_string(index=False))
    print(R[R.x == "q"][["subsample", "y", "N", "r_P", "p_P", "r_S", "p_S", "slope", "slope_err", "intercept", "n_clipped"]].round(3).to_string(index=False))
    print(R[R.x != "q"][["subsample", "y", "x", "N", "r_P", "p_P", "r_S", "p_S"]].round(3).to_string(index=False))
    print(P[["property", "N", "r_S", "lo", "hi", "excludes_minus044", "power_044"]].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
