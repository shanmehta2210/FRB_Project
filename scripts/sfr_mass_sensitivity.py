"""Sensitivity checks for stellar mass, SFR and sSFR, as for DM_host.

For each property against q (and against rest-frame DM_host), in the full sample and in subsamples:
Spearman and Kendall, bootstrap 95% interval, delete-1/2/3 jackknife over all subsets, Theil-Sen slope,
the AD test between the face-on and edge-on halves, and (against q) 2,000 Monte Carlo draws of
q ~ N(q, sigma_q_B).

Output: Publication/materials/appendix/robust_correlations/{sfr_mass_sensitivity.csv, SFR_MASS_SENSITIVITY.md}
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import redshift_bins as ZB  # noqa: E402
import robust_correlations as RC  # noqa: E402
import robust_extras as E  # noqa: E402

N_MC = 2_000
Q0 = 0.2
PROPS = [("logmstar", "log M*"), ("log_sfr", "log SFR"), ("log_ssfr", "log sSFR")]


def subsamples(h: pd.DataFrame) -> list[tuple[str, pd.DataFrame]]:
    return [("All", h), ("m_r < 21", h[h.m_r < 21]), ("n < 2", h[h.sersic_n < 2]),
            ("bulge fits removed", h[~h.frb.isin(RC.BULGE_FITS)]), ("disks only", h[h.cls.isin(ZB.DISK)]),
            ("z < 0.2", h[h.z < 0.2]), ("z >= 0.2", h[h.z >= 0.2])]


def mc_q(q, sq, y, rng):
    r = np.empty(N_MC)
    p = np.empty(N_MC)
    for k in range(N_MC):
        qq = np.clip(q + rng.normal(0, 1, len(q)) * np.nan_to_num(sq), Q0 + 1e-3, 1.0)
        r[k], p[k] = stats.spearmanr(qq, y)
    return np.percentile(r, [16, 50, 84]), float(np.mean(p < 0.05))


def run(lab, sub, xcol, ycol, rng, mc=False) -> dict:
    d = sub[np.isfinite(sub[xcol]) & np.isfinite(sub[ycol])]
    x, y = d[xcol].to_numpy(float), d[ycol].to_numpy(float)
    row = dict(sample=lab, x=xcol, y=ycol, **RC.robust_row(d.frb.to_numpy(), x, y, rng))
    for k in (2, 3):
        if len(x) > k + 5:
            lo, hi, pmin, pmax, fsig = E.jackknife(x, y, k)
            row.update({f"jk{k}_min": lo, f"jk{k}_max": hi, f"jk{k}_pmin": pmin, f"jk{k}_fsig": fsig})
    row["p_AD"] = ZB.ad_halves(x, y, rng)
    if mc:
        (r16, r50, r84), fsig = mc_q(x, d.sigma_q_B.to_numpy(float), y, rng)
        row.update(mc_r16=r16, mc_r50=r50, mc_r84=r84, mc_fsig=fsig)
    return row


def main() -> None:
    h = ZB.hosts()
    h["log_sfr"] = np.log10(h.sfr.where(h.sfr > 0))
    rng = np.random.default_rng(5)
    rows = []
    for col, _ in PROPS:
        for lab, sub in subsamples(h):
            rows.append(run(lab, sub, "q", col, rng, mc=True))
            print("q", col, lab, rows[-1]["N"], round(rows[-1]["r_S"], 3), round(rows[-1]["p_S"], 3), flush=True)
    for col, _ in PROPS:
        for lab, sub in subsamples(h)[:1] + subsamples(h)[-2:]:
            rows.append(run(lab, sub, col, "dm_host_rest", rng))
            print("dm", col, lab, rows[-1]["N"], round(rows[-1]["r_S"], 3), round(rows[-1]["p_S"], 3), flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(RC.OUT / "sfr_mass_sensitivity.csv", index=False)

    lab = dict(PROPS, q="q", dm_host_rest="DM_host")
    md = ["# Sensitivity checks for stellar mass, SFR and sSFR", "",
          "Bootstrap 95% interval (10,000 resamples); JK1/JK3 = range of r_S with 1 or 3 hosts removed (all subsets), "
          "f = fraction of 3-host subsets with p < 0.05; AD = face-on vs edge-on halves; MC = r_S 16/50/84 and fraction "
          f"p < 0.05 over {N_MC:,} draws of q ~ N(q, sigma_q_B). SFR upper limits are excluded.", "",
          "| x | y | sample | N | r_S (p) | 95% CI | Kendall tau (p) | JK1 | JK3 (f) | Theil-Sen slope | p_AD | MC r_S | MC f |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in d.itertuples():
        jk3 = f"{r.jk3_min:+.2f} to {r.jk3_max:+.2f} ({r.jk3_fsig:.2f})" if np.isfinite(getattr(r, "jk3_min", np.nan)) else ""
        mc = (f"{r.mc_r16:+.2f} / {r.mc_r50:+.2f} / {r.mc_r84:+.2f}", f"{r.mc_fsig:.2f}") \
            if np.isfinite(getattr(r, "mc_r50", np.nan)) else ("", "")
        md.append(f"| {lab[r.x]} | {lab[r.y]} | {r.sample} | {r.N} | {r.r_S:+.2f} ({r.p_S:.3f}) | "
                  f"[{r.rS_boot_lo:+.2f}, {r.rS_boot_hi:+.2f}] | {r.tau:+.2f} ({r.p_tau:.3f}) | "
                  f"{r.rS_jack_min:+.2f} to {r.rS_jack_max:+.2f} | {jk3} | {r.theil_sen_slope:+.3g} "
                  f"[{r.ts_lo:+.3g}, {r.ts_hi:+.3g}] | {r.p_AD:.3f} | {mc[0]} | {mc[1]} |")
    (RC.OUT / "SFR_MASS_SENSITIVITY.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))


if __name__ == "__main__":
    main()
