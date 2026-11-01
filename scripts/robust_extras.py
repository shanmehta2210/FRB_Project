"""Delete-d jackknife (d = 1, 2, 3; all subsets) and a rank-insensitive companion test for each
correlation with host axis ratio.

Companion test: two-sample Anderson-Darling between the property values of the more face-on half
(q above the median) and the more edge-on half of the hosts, with a permutation p-value. Spearman
tests for a monotonic trend; AD tests whether the two halves differ in distribution at all.

Rows: the properties of robust_correlations.py, the DM_host subsamples, and DM_host for the sample of
Marnoch et al. (2026) alone and combined with ours.

Output: Publication/materials/appendix/robust_correlations/{robust_extras.csv, ROBUST_EXTRAS.md}
"""
from __future__ import annotations

import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import marnoch_dm_combined as MDC  # noqa: E402
import robust_correlations as RC  # noqa: E402

ROOT = RC.ROOT
MR = ROOT / "Publication" / "materials" / "marnoch_recheck"
N_PERM = 9_999


def spearman_many(x: np.ndarray, y: np.ndarray, drop: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Spearman r_S and p for every row of `drop` (indices removed)."""
    n = len(x)
    keep = np.ones((len(drop), n), bool)
    np.put_along_axis(keep, drop, False, axis=1)
    m = n - drop.shape[1]
    xs = x[np.nonzero(keep)[1].reshape(len(drop), m)]
    ys = y[np.nonzero(keep)[1].reshape(len(drop), m)]
    rx = stats.rankdata(xs, axis=1)
    ry = stats.rankdata(ys, axis=1)
    rx -= rx.mean(1, keepdims=True)
    ry -= ry.mean(1, keepdims=True)
    r = (rx * ry).sum(1) / np.sqrt((rx ** 2).sum(1) * (ry ** 2).sum(1))
    t = r * np.sqrt((m - 2) / np.clip(1 - r ** 2, 1e-12, None))
    return r, 2 * stats.t.sf(np.abs(t), m - 2)


def jackknife(x, y, d):
    drop = np.array(list(combinations(range(len(x)), d)))
    r, p = spearman_many(x, y, drop)
    return r.min(), r.max(), p.min(), p.max(), float(np.mean(p < 0.05))


def ad_halves(x, y, rng):
    med = np.median(x)
    a, b = y[x > med], y[x <= med]
    res = stats.anderson_ksamp([a, b], method=stats.PermutationMethod(n_resamples=N_PERM, random_state=rng))
    return float(res.statistic), float(res.pvalue)


def samples() -> list[tuple[str, np.ndarray, np.ndarray]]:
    h = RC.load()
    out = []
    for col, lab in RC.PROPS:
        if col in ("dm_host_fg_cond", "dm_host_fg_spec"):
            continue
        out.append((lab, h.q.to_numpy(float), h[col].to_numpy(float)))
    for lab, sel in (("m_r < 21", h.m_r < 21), ("n < 2", h.sersic_n < 2), ("bulge fits removed", ~h.frb.isin(RC.BULGE_FITS))):
        s = h[sel]
        out.append((rf"DM$_{{\rm host}}$, {lab}", s.q.to_numpy(float), s.dm_host_rest.to_numpy(float)))
    m = pd.read_csv(MR / "marnoch2026_table1.csv", dtype={"frb": str})
    ours = h.dropna(subset=["dm_host_rest"])
    out.append(("Marnoch et al. (2026) sample", m.q.to_numpy(float), m.dm_host_rest.to_numpy(float)))
    m, _ = MDC.marnoch_on_our_scale(m, ours)
    mo = m[~m.frb.isin(ours.frb)]
    uq = np.r_[ours.q, mo.q]
    ud = np.r_[ours.dm_host_rest, mo.dm_host_rest]
    out.append(("Both samples (shared hosts from this work)", uq, ud))
    oo = ours[~ours.frb.isin(m.frb)]
    out.append(("Both samples (shared hosts from Marnoch)", np.r_[oo.q, m.q], np.r_[oo.dm_host_rest, m.dm_host_rest]))
    out.append(("Both samples, $q > 0.2$", uq[uq > 0.2], ud[uq > 0.2]))
    return out


def main() -> None:
    rng = np.random.default_rng(7)
    rows = []
    for lab, x, y in samples():
        ok = np.isfinite(x) & np.isfinite(y)
        x, y = x[ok], y[ok]
        rs, ps = stats.spearmanr(x, y)
        row = dict(label=lab, N=len(x), r_S=rs, p_S=ps)
        for d in (1, 2, 3):
            lo, hi, pmin, pmax, fsig = jackknife(x, y, d)
            row.update({f"jk{d}_min": lo, f"jk{d}_max": hi, f"jk{d}_pmin": pmin, f"jk{d}_pmax": pmax, f"jk{d}_fsig": fsig})
        row["ad_stat"], row["p_AD"] = ad_halves(x, y, rng)
        rows.append(row)
        print(f"{lab[:45]:45s} N={len(x):2d} rS={rs:+.2f} p={ps:.3f} | "
              + " | ".join(f"d={d}: [{row[f'jk{d}_min']:+.2f},{row[f'jk{d}_max']:+.2f}] p in [{row[f'jk{d}_pmin']:.3f},{row[f'jk{d}_pmax']:.3f}] sig {row[f'jk{d}_fsig']:.2f}"
                           for d in (1, 2, 3))
              + f" | AD p={row['p_AD']:.3f}", flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(RC.OUT / "robust_extras.csv", index=False)
    md = ["# Delete-d jackknife and Anderson-Darling companion tests", "",
          "Jackknife over all subsets with d hosts removed: range of $r_S$, range of its $p$, and fraction of subsets with $p < 0.05$. "
          f"AD: face-on half ($q$ above median) vs edge-on half, permutation $p$ ({N_PERM:,} resamples).", "",
          "| sample | N | $r_S$ ($p$) | d = 1 | d = 2 | d = 3 | frac. $p<0.05$ (d = 1/2/3) | AD $p$ |", "|---|---|---|---|---|---|---|---|"]
    for r in d.itertuples():
        md.append(f"| {r.label} | {r.N} | {r.r_S:+.2f} ({r.p_S:.3f}) | "
                  + " | ".join(f"{getattr(r, f'jk{k}_min'):+.2f} to {getattr(r, f'jk{k}_max'):+.2f}" for k in (1, 2, 3))
                  + f" | {r.jk1_fsig:.2f} / {r.jk2_fsig:.2f} / {r.jk3_fsig:.2f} | {r.p_AD:.3f} |")
    (RC.OUT / "ROBUST_EXTRAS.md").write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
