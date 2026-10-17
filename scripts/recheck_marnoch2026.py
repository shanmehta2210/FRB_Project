"""
Re-run the Marnoch et al. (2026) DM_host vs axis-ratio tests from their Table 1.

Table 1 gives, per CRAFT host: z, PATH P, the observer-frame DM_host proxy
(DM_FRB - DM_MW - DM_halo - <DM_cosmic>(z), "DM_host^Macquart"), GALFIT q, i, n, m.
Rest frame = observer frame x (1+z). Their Section 6 quotes Pearson / Spearman for the
full sample (N = 34), m_galfit <= 21, n < 2, and the sample without the three
bulge-fit hosts (20190608B, 20220725A, 20231226A).

Variants: as published; without the four q < 0.2 hosts that they set to i = 90 deg
(our sample cut is q > 0.2); additionally without the repeater 20190711A.

Writes Publication/materials/marnoch_recheck/MARNOCH_RECHECK.md and .csv.
"""
from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Publication" / "materials" / "marnoch_recheck"
TEXT = ROOT / "Publication" / "materials" / "host_properties" / "papers" / "marnoch2026.txt"
Q0 = 0.2
BULGE = {"20190608B", "20220725A", "20231226A"}
REPEATERS = {"20190711A"}

ROW = re.compile(
    r"^(\d{8}[A-Z])[^\d]*(\d\.\d+)\s+([\d.]+)([gIR])(-?\d+)\((\d+)\)([\d.]+)\((\d+)\)(\d+)\s+"
    r"([\d.]+)\((\d+)\)([\d.]+)\((\d+)\)"
)


def parse_table(text: str) -> pd.DataFrame:
    rows = []
    for line in text.splitlines():
        line = line.replace("\u2212", "-").strip()
        m = ROW.match(line)
        if not m:
            continue
        g = m.groups()
        rows.append({"frb": g[0], "z": float(g[1]), "path_p": float(g[2]), "band": g[3],
                     "dm_host_obs": float(g[4]), "dm_host_obs_err": float(g[5]),
                     "q": float(g[6]), "i_table": float(g[8]), "n": float(g[9]), "m_galfit": float(g[11])})
    t = pd.DataFrame(rows).drop_duplicates("frb").reset_index(drop=True)
    t["dm_host_rest"] = t.dm_host_obs * (1 + t.z)
    c2 = ((t.q ** 2 - Q0 ** 2) / (1 - Q0 ** 2)).clip(lower=0, upper=1)
    t["i_hubble"] = np.degrees(np.arccos(np.sqrt(c2)))
    return t


def corr(x: np.ndarray, y: np.ndarray) -> dict:
    rp = stats.pearsonr(x, y)
    rs = stats.spearmanr(x, y)
    n = len(x)
    se = 1.06 / math.sqrt(n - 3)
    zr = math.atanh(rs.statistic)
    return {"N": n, "r_P": rp.statistic, "p_P": rp.pvalue, "r_S": rs.statistic, "p_S": rs.pvalue,
            "r_S_lo95": math.tanh(zr - 1.96 * se), "r_S_hi95": math.tanh(zr + 1.96 * se)}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    t = parse_table(TEXT.read_text(encoding="utf-8", errors="replace"))
    t.to_csv(OUT / "marnoch2026_table1.csv", index=False)
    lowq = set(t.loc[t.q < Q0, "frb"])

    variants = {
        "as published": set(),
        "drop q < 0.2": lowq,
        "drop q < 0.2 and repeater": lowq | REPEATERS,
    }
    subsamples = {
        "full": lambda d: d,
        "m_galfit <= 21": lambda d: d[d.m_galfit <= 21.0],
        "n < 2": lambda d: d[d.n < 2.0],
        "no bulge fits": lambda d: d[~d.frb.isin(BULGE)],
    }
    rows = []
    for vname, drop in variants.items():
        base = t[~t.frb.isin(drop)]
        for sname, f in subsamples.items():
            d = f(base)
            for dm in ("dm_host_obs", "dm_host_rest"):
                for x in ("q", "i_table", "i_hubble"):
                    if vname == "as published" and x == "i_hubble":
                        continue
                    if vname != "as published" and x == "i_table":
                        continue
                    r = corr(d[x].to_numpy(float), d[dm].to_numpy(float))
                    r.update({"variant": vname, "subsample": sname, "dm": dm, "x": x})
                    rows.append(r)
    R = pd.DataFrame(rows)
    R.to_csv(OUT / "marnoch_recheck.csv", index=False)

    claims = {
        ("full", "dm_host_obs", "q"): "r_P = -0.399 (p = 0.02), r_S = -0.425 (p = 0.012)",
        ("full", "dm_host_rest", "q"): "r_P = -0.402 (p = 0.018), r_S = -0.437 (p = 0.01)",
        ("full", "dm_host_rest", "i_table"): "p < 0.05",
        ("m_galfit <= 21", "dm_host_rest", "q"): "r_P = -0.538, r_S = -0.544 (p ~ 0.02)",
        ("n < 2", "dm_host_rest", "q"): "r_P = -0.43 (p = 0.022), r_S = -0.40 (p = 0.036)",
        ("no bulge fits", "dm_host_rest", "q"): "significant, moderately strong",
    }

    def f(x, nd=3):
        return f"{x:.{nd}f}"

    md = ["# Re-check of Marnoch et al. (2026) DM_host vs axis ratio", "",
          "Source: `scripts/recheck_marnoch2026.py`, using the per-host values of their Table 1 "
          f"(N = {len(t)}). DM_host is their observer-frame Macquart proxy; rest frame = observer × (1+z). "
          f"Hosts with $q < 0.2$ that they set to $i = 90^\\circ$: {', '.join(sorted(lowq))}. "
          "Our sample cut is $q > 0.2$, so these are dropped in the second variant; the third also drops the "
          "repeater 20190711A. $i$ is their tabulated value (as published) or Hubble with $q_0 = 0.2$. "
          "95% CI on $r_S$ from the Fisher transform.", "",
          "## 1. Reproduction of the published numbers", "",
          "| subsample | DM_host | vs | N | $r_P$ (p) | $r_S$ (p) | paper |", "|---|---|---|---|---|---|---|"]
    for _, r in R[R.variant == "as published"].iterrows():
        claim = claims.get((r.subsample, r.dm, r.x), "")
        md.append(f"| {r.subsample} | {r.dm} | {r.x} | {r.N} | {f(r.r_P)} ({f(r.p_P)}) | "
                  f"{f(r.r_S)} ({f(r.p_S)}) | {claim} |")
    for vname in list(variants)[1:]:
        md += ["", f"## {list(variants).index(vname) + 1}. {vname}", "",
               "| subsample | DM_host | vs | N | $r_P$ (p) | $r_S$ (p) | $r_S$ 95% CI |",
               "|---|---|---|---|---|---|---|"]
        for _, r in R[R.variant == vname].iterrows():
            md.append(f"| {r.subsample} | {r.dm} | {r.x} | {r.N} | {f(r.r_P)} ({f(r.p_P)}) | "
                      f"{f(r.r_S)} ({f(r.p_S)}) | [{f(r.r_S_lo95, 2)}, {f(r.r_S_hi95, 2)}] |")
    lo = t[t.frb.isin(lowq)][["frb", "q", "z", "dm_host_obs", "dm_host_rest", "m_galfit", "n"]]
    md += ["", "## Dropped low-$q$ hosts", "", "| FRB | q | z | DM_host obs | DM_host rest | m | n |",
           "|---|---|---|---|---|---|---|"]
    for _, r in lo.iterrows():
        md.append(f"| {r.frb} | {r.q:.4f} | {r.z:.4f} | {r.dm_host_obs:.0f} | {r.dm_host_rest:.0f} | "
                  f"{r.m_galfit:.2f} | {r.n:.2f} |")
    ours = pd.read_csv(ROOT / "plots" / "plots_null" / "v2" / "host_correlations" / "host_props.csv",
                       dtype={"frb": str})[["frb", "q", "dm_host_rest"]]
    J = t.merge(ours, on="frb", suffixes=("_M", "_O"))
    md += ["", f"## Hosts shared with our sample (N = {len(J)})", "",
           "Subscript M = Marnoch Table 1, O = ours (DM_host without halo term, so offset by a near-constant).", "",
           "| q | DM_host rest | N | $r_S$ (p) |", "|---|---|---|---|"]
    for qa in ("q_M", "q_O"):
        for da in ("dm_host_rest_M", "dm_host_rest_O"):
            rs = stats.spearmanr(J[qa], J[da])
            md.append(f"| {qa} | {da} | {len(J)} | {f(rs.statistic)} ({f(rs.pvalue)}) |")
    rq, rd = stats.spearmanr(J.q_M, J.q_O), stats.spearmanr(J.dm_host_rest_M, J.dm_host_rest_O)
    md += ["", f"Agreement: $q$ $r_S$ = {rq.statistic:.3f}; DM_host $r_S$ = {rd.statistic:.3f}.",
           f"Their hosts not in our sample: {', '.join(sorted(set(t.frb) - set(J.frb)))}.", ""]
    (OUT / "MARNOCH_RECHECK.md").write_text("\n".join(md), encoding="utf-8")
    sys.stdout.reconfigure(encoding="utf-8")
    print("\n".join(md))


if __name__ == "__main__":
    main()
