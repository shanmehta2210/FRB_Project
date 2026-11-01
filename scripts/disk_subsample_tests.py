"""Main-text inclination tests restricted by the disk evidence of disk_evidence.py.

Output: Publication/materials/appendix/morphology/disk_subsample_tests.{csv,md}
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import appendix_robustness as A  # noqa: E402
import plot_paper_inclination_figures as P  # noqa: E402

OUT = P.REPO_ROOT / "Publication" / "materials" / "appendix" / "morphology"


def main() -> None:
    h = P.load_hosts().set_index("frb")
    ev = pd.read_csv(OUT / "disk_evidence.csv", dtype={"frb": str}).set_index("frb")
    cls = pd.Series("reliable", index=h.index)
    cls.loc[ev.index] = ev.cls
    pool = A.cosi_from_q(A.select(A.load_parent())["ba"].to_numpy(), A.Q0)
    variants = [("All hosts", cls.index),
                ("Without the early-type host", cls.index[cls != "early"]),
                ("Without early-type and ambiguous hosts", cls.index[~cls.isin(["early", "ambig"])]),
                ("Reliable Sb--Irr or literature disks only", cls.index[cls.isin(["reliable", "lit"])])]
    rows = []
    for lab, idx in variants:
        r = A.matched_tests(h.loc[idx, "cosi_point"].to_numpy(float), pool, None, 10_000)
        r.update(label=lab)
        rows.append(r)
        print(lab, r["N_frb"], r["quartile_counts"], round(r["ad_pm"], 3), round(r["mw2_pm"], 3), round(r["mw1_pm"], 3), flush=True)
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "disk_subsample_tests.csv", index=False)
    md = ["| sample | N | quartile counts | AD | MW 2-sided | MW 1-sided |", "|---|---|---|---|---|---|"]
    md += [f"| {r.label} | {r.N_frb} | {r.quartile_counts} | {r.ad_pm:.2f} | {r.mw2_pm:.2f} | {r.mw1_pm:.2f} |" for r in d.itertuples()]
    (OUT / "disk_subsample_tests.md").write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
