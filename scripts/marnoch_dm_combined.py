"""Our DM_host against Marnoch et al. (2026) for shared hosts, and DM_host-q correlations
in the union of both samples.

Output: Publication/materials/marnoch_recheck/{dm_host_shared.csv, combined_correlations.csv, COMBINED.md}
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import robust_correlations as RC  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
MR = ROOT / "Publication" / "materials" / "marnoch_recheck"


def rstats(x, y, rng):
    x, y = np.asarray(x, float), np.asarray(y, float)
    rs, ps = stats.spearmanr(x, y)
    tau, pt = stats.kendalltau(x, y)
    idx = rng.integers(0, len(x), (10_000, len(x)))
    boot = np.array([stats.spearmanr(x[i], y[i])[0] for i in idx])
    jack = np.array([stats.spearmanr(np.delete(x, k), np.delete(y, k))[0] for k in range(len(x))])
    return dict(N=len(x), r_S=rs, p_S=ps, boot_lo=np.nanpercentile(boot, 2.5), boot_hi=np.nanpercentile(boot, 97.5),
                tau=tau, p_tau=pt, jack_min=jack.min(), jack_max=jack.max())


def main() -> None:
    rng = np.random.default_rng(42)
    m = pd.read_csv(MR / "marnoch2026_table1.csv", dtype={"frb": str})
    h = RC.load()
    ours = h[["frb", "q", "dm_host_rest", "dm_host_fg", "dm_host_fg_excess"]].dropna(subset=["dm_host_rest"])

    sh = m.merge(ours, on="frb", suffixes=("_m", ""))
    sh = sh[["frb", "z", "q_m", "q", "dm_host_rest_m", "dm_host_rest", "dm_host_fg", "dm_host_fg_excess"]]
    sh.to_csv(MR / "dm_host_shared.csv", index=False)
    md = ["# DM_host: this work vs Marnoch et al. (2026)", "", f"Shared hosts with DM_host in both: {len(sh)}", ""]
    for col, lab in (("dm_host_rest", "ours, no foreground correction"), ("dm_host_fg", "ours, foreground removed"),
                     ("dm_host_fg_excess", "ours, foreground excess removed")):
        d = sh.dropna(subset=[col])
        diff = d[col] - d.dm_host_rest_m
        md.append(f"- {lab}: N={len(d)}, r_S={stats.spearmanr(d[col], d.dm_host_rest_m)[0]:.2f}, "
                  f"median diff={diff.median():+.0f}, median |diff|={diff.abs().median():.0f} pc cm^-3")

    rows = []
    mo = m[~m.frb.isin(ours.frb)]
    variants = {
        "Marnoch alone": (m.q, m.dm_host_rest),
        "This work alone": (ours.q, ours.dm_host_rest),
        "Union, shared hosts from this work": (pd.concat([ours.q, mo.q]), pd.concat([ours.dm_host_rest, mo.dm_host_rest])),
        "Union, shared hosts from Marnoch": (
            pd.concat([ours[~ours.frb.isin(m.frb)].q, m.q]),
            pd.concat([ours[~ours.frb.isin(m.frb)].dm_host_rest, m.dm_host_rest])),
        "Union, q > 0.2 only": None,
    }
    uq = pd.concat([ours.q, mo.q]).to_numpy()
    ud = pd.concat([ours.dm_host_rest, mo.dm_host_rest]).to_numpy()
    variants["Union, q > 0.2 only"] = (uq[uq > 0.2], ud[uq > 0.2])
    for lab, (x, y) in variants.items():
        r = rstats(x, y, rng)
        r["sample"] = lab
        rows.append(r)
    t = pd.DataFrame(rows)[["sample", "N", "r_S", "p_S", "boot_lo", "boot_hi", "tau", "p_tau", "jack_min", "jack_max"]]
    t.to_csv(MR / "combined_correlations.csv", index=False)
    md += ["", "## DM_host vs q", "", "| sample | N | r_S | p | 95% bootstrap | Kendall tau (p) | jackknife r_S |",
           "|---|---|---|---|---|---|---|"]
    for r in t.itertuples():
        md.append(f"| {r.sample} | {r.N} | {r.r_S:+.2f} | {r.p_S:.3f} | ({r.boot_lo:+.2f}, {r.boot_hi:+.2f}) | "
                  f"{r.tau:+.2f} ({r.p_tau:.3f}) | {r.jack_min:+.2f} to {r.jack_max:+.2f} |")
    (MR / "COMBINED.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print("\n".join(md))


if __name__ == "__main__":
    main()
