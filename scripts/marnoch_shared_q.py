"""Axis ratios of the hosts shared with Marnoch et al. (2026): agreement and dependence on host size.

Output: Publication/materials/marnoch_recheck/shared_q_comparison.csv (+ printed summary)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_paper_inclination_figures as P  # noqa: E402

ROOT = P.REPO_ROOT
MR = ROOT / "Publication" / "materials" / "marnoch_recheck"


def main() -> None:
    m = pd.read_csv(MR / "marnoch2026_table1.csv", dtype={"frb": str})[["frb", "q", "n", "m_galfit"]]
    h = P.load_hosts()[["frb", "q", "re", "fwhm", "re_over_fwhm", "mag"]]
    s = pd.read_csv(ROOT / "Publication" / "materials" / "host_properties" / "seed_56.csv")[["frb", "re_arcsec_galfit", "q_astrophot", "imaging"]]
    d = m.merge(h, on="frb", suffixes=("_marnoch", "_ours")).merge(s, on="frb")
    d["dq"] = d.q_ours - d.q_marnoch
    d = d.sort_values("re_arcsec_galfit")
    d.to_csv(MR / "shared_q_comparison.csv", index=False)
    print(d[["frb", "q_marnoch", "q_ours", "q_astrophot", "dq", "re_arcsec_galfit", "re_over_fwhm", "imaging", "mag"]].to_string(index=False))
    rs = stats.spearmanr(d.q_marnoch, d.q_ours)
    print(f"N={len(d)} r_S={rs[0]:.3f} p={rs[1]:.2g} median|dq|={d.dq.abs().median():.3f} mean dq={d.dq.mean():+.3f}")
    out = d.dq.abs() > 3 * 1.4826 * np.median(np.abs(d.dq - d.dq.median()))
    print("outliers (>3 MAD):", d.loc[out, "frb"].tolist())
    rest = d[~out]
    print(f"without outliers: N={len(rest)} r_S={stats.spearmanr(rest.q_marnoch, rest.q_ours)[0]:.3f} median|dq|={rest.dq.abs().median():.3f}")


if __name__ == "__main__":
    main()
