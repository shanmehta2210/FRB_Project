"""How many FRB sightlines intersect a foreground galaxy halo (output of foreground_halo_dm.py).

A galaxy contributes DM only when the sightline passes within r200 of its halo, so a galaxy intersects the
sightline when its Monte Carlo DM is non-zero. Levels:
  any     : DM > 0 in some draws (dm_mean > 0)
  likely  : DM > 0 in at least 84% of draws (dm_p16 > 0)
  spec-z  : likely, with a spectroscopic redshift below the host's
  > 10    : mean contribution above 10 pc cm^-3 (observer frame), weighted by the foreground probability

Output: Publication/materials/foreground_halos/SIGHTLINES.md, sightline_counts.csv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

FG = Path(__file__).resolve().parents[1] / "Publication" / "materials" / "foreground_halos"


def main() -> None:
    hosts = pd.read_csv(FG / "foreground_dm_per_host.csv", dtype={"frb": str})
    hosts = hosts[hosts.in_ls.astype(str).isin(["dr9", "dr10"])]
    g = pd.read_csv(FG / "foreground_galaxies.csv", dtype={"frb": str})
    g = g[g.frb.isin(hosts.frb)]
    g["any"] = g.dm_mean > 0
    g["likely"] = g.dm_p16 > 0
    g["specz"] = g.likely & np.isfinite(g.z_spec)
    g["big"] = g.dm_mean * g.p_fg > 10
    per = g.groupby("frb")[["any", "likely", "specz", "big"]].sum().reindex(hosts.frb, fill_value=0)
    per["dm_fg"] = hosts.set_index("frb").dm_fg
    per.to_csv(FG / "sightline_counts.csv")
    n = len(per)
    lines = ["# Foreground halos intersected by the FRB sightlines", "",
             f"{n} hosts in the LS DR9/DR10 footprint; {len(g):,} candidate foreground galaxies within 0.5 deg.", "",
             "| level | hosts with >= 1 | galaxies | median per host (hosts with >= 1) |", "|---|---|---|---|"]
    for col, lab in (("any", "intersected in some draws"), ("likely", "intersected in >= 84% of draws"),
                     ("specz", "likely, spectroscopic z"), ("big", "contributes > 10 pc cm^-3")):
        k = per[col] > 0
        lines.append(f"| {lab} | {k.sum()} / {n} | {int(per[col].sum())} | {per.loc[k, col].median():.0f} |")
    lines += ["", f"Hosts with total foreground DM > 50 pc cm^-3 (observer frame): {(per.dm_fg > 50).sum()} / {n}; "
              f"> 100: {(per.dm_fg > 100).sum()} / {n}."]
    (FG / "SIGHTLINES.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
