"""Sensitivity of the AD and MW tests to the edge-on deficit of Bhardwaj et al. (2024).

For each sample size, draw N_BOOT samples with replacement from the 23 published
inclinations and test each against the full HSC control pool. The detection rate is
the fraction of samples with p < 0.05.

Output: Publication/materials/statistics/power_bhardwaj_deficit.json
"""
from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_paper_inclination_figures as P  # noqa: E402

warnings.filterwarnings("ignore")
N_BOOT = 2_000
OUT = P.REPO_ROOT / "Publication" / "materials" / "statistics" / "power_bhardwaj_deficit.json"


def published_cosi() -> np.ndarray:
    inc = []
    for line in P.NATURE.read_text(encoding="utf-8").splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) >= 3:
            inc.append(int(parts[2].replace('"', "").split(",")[0]))
    assert len(inc) == 23, len(inc)
    return np.cos(np.radians(np.asarray(inc, dtype=float)))


def main() -> None:
    hsc, _ = P.load_hsc_pool(mag_limit=P.MAG_LIMIT, q0=P.Q0)
    pub = published_cosi()
    n_hosts = len(P.load_hosts())
    rng = np.random.default_rng(P.SEED + 7)
    out = {"n_boot": N_BOOT, "n_pool": int(len(hsc)), "seed": P.SEED + 7}
    for n in (n_hosts, 23):
        ad = mw = 0
        for _ in range(N_BOOT):
            s = rng.choice(pub, n, replace=True)
            ad += stats.anderson_ksamp([hsc, s]).pvalue < 0.05
            mw += stats.mannwhitneyu(s, hsc, alternative="two-sided").pvalue < 0.05
        out[f"N{n}"] = {"ad_detect": ad / N_BOOT, "mw_detect": mw / N_BOOT}
        print(n, out[f"N{n}"])
    OUT.write_text(json.dumps(out, indent=1), encoding="utf-8")


if __name__ == "__main__":
    main()
