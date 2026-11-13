"""How many FRB sightlines intersect a foreground galaxy halo (output of foreground_halo_dm.py).

A galaxy contributes DM only when the sightline passes within r200 of its halo, so a galaxy intersects the
sightline when its Monte Carlo DM is non-zero. Photo-z draws are conditional on the galaxy lying in front of
the host, so a galaxy's intersection probability is p_fg times the fraction of its draws with DM > 0.
Levels:
  any     : DM > 0 in some draws (dm_mean > 0), ignoring p_fg
  likely  : DM > 0 in at least 84% of draws (dm_p16 > 0), ignoring p_fg
  in front: likely, with p_fg >= 0.84
  spec-z  : likely, with a spectroscopic redshift below the host's
  > 10    : mean contribution above 10 pc cm^-3 (observer frame), weighted by p_fg, counted in the total
            (not a satellite of a more massive candidate, not host environment)
The expected number of hosts with at least one intersection uses 1 - prod(1 - 0.84 p_fg) over likely
galaxies (a lower bound), and is compared with random sightlines and with the halo mass function.

Output: Publication/materials/foreground_halos/SIGHTLINES.md, sightline_counts.csv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from colossus.cosmology import cosmology
from colossus.halo import mass_so
from colossus.lss import mass_function

FG = Path(__file__).resolve().parents[1] / "Publication" / "materials" / "foreground_halos"
LOGMH_MIN = 11.0


def halo_crossings(z_host, logmh_min=LOGMH_MIN):
    """Mean number of halos above logmh_min whose r200c covers a random sightline to z_host (Tinker+08)."""
    cosmo = cosmology.setCosmology("planck18")
    h = cosmo.h
    m = 10 ** np.linspace(logmh_min, 15, 100) * h
    zg = np.linspace(0.002, max(z_host.max(), 0.01) + 0.01, 120)
    dn = []
    for z in zg:
        n = mass_function.massFunction(m, z, mdef="200c", model="tinker08", q_out="dndlnM") * h ** 3
        r = mass_so.M_to_R(m, z, "200c") / 1e3 / h
        dn.append(np.trapezoid(n * np.pi * r ** 2 * (1 + z) ** 2, np.log(m)) * 2.998e5 / cosmo.Hz(z))
    dn = np.array(dn)
    cum = np.concatenate([[0], np.cumsum(np.diff(zg) * (dn[1:] + dn[:-1]) / 2)])
    return np.interp(z_host, zg, cum)


def main() -> None:
    in_ls = lambda d: d[d.in_ls.astype(str).isin(["dr9", "dr10"])]
    hosts = in_ls(pd.read_csv(FG / "foreground_dm_per_host.csv", dtype={"frb": str}))
    controls = in_ls(pd.read_csv(FG / "foreground_dm_controls.csv", dtype={"frb": str}))
    g = pd.read_csv(FG / "foreground_galaxies.csv", dtype={"frb": str})
    g = g[g.frb.isin(hosts.frb)]
    g["any"] = g.dm_mean > 0
    g["likely"] = g.dm_p16 > 0
    g["front"] = g.likely & (g.p_fg >= 0.84)
    g["specz"] = g.likely & np.isfinite(g.z_spec)
    g["big"] = g.counted & (g.dm_mean * g.p_fg > 10)
    g["p_int"] = np.where(g.likely, 0.84 * g.p_fg, 0.0)
    per = g.groupby("frb")[["any", "likely", "front", "specz", "big"]].sum().reindex(hosts.frb, fill_value=0)
    per["p_any"] = (1 - (1 - g.p_int).groupby(g.frb).prod()).reindex(hosts.frb, fill_value=0)
    per["dm_fg"] = hosts.set_index("frb").dm_fg
    per.to_csv(FG / "sightline_counts.csv")

    n = len(per)
    theory = 1 - np.exp(-halo_crossings(hosts.z.to_numpy()))
    lines = ["# Foreground halos intersected by the FRB sightlines", "",
             f"{n} hosts in the LS DR9/DR10 footprint; {len(g):,} candidate foreground galaxies within 0.5 deg. "
             "A galaxy intersects when the sightline passes within its r200. Photo-z draws assume the galaxy is in "
             "front of the host, so the first two rows ignore the probability p_fg that it is; half of the 'likely' "
             f"galaxies have p_fg < {g.loc[g.likely, 'p_fg'].median():.2f}.", "",
             "| level | hosts with >= 1 | galaxies | median per host (hosts with >= 1) |", "|---|---|---|---|"]
    for col, lab in (("any", "intersected in some draws (ignores p_fg)"),
                     ("likely", "intersected in >= 84% of draws (ignores p_fg)"),
                     ("front", "likely, and p_fg >= 0.84"),
                     ("specz", "likely, spectroscopic z"),
                     ("big", "contributes > 10 pc cm^-3")):
        k = per[col] > 0
        lines.append(f"| {lab} | {k.sum()} / {n} | {int(per[col].sum())} | {per.loc[k, col].median():.0f} |")
    lines += ["", f"Expected number of hosts with at least one intersection, including p_fg (lower bound): "
              f"{per.p_any.sum():.1f} / {n}. Halo mass function (Tinker+08, halos above 10^{LOGMH_MIN:g} M_sun, "
              f"random sightlines to the same redshifts): {theory.sum():.1f} / {n}.", "",
              f"Hosts with total foreground DM > 50 pc cm^-3 (observer frame): {(per.dm_fg > 50).sum()} / {n}; "
              f"> 100: {(per.dm_fg > 100).sum()} / {n}.", "",
              "| total foreground DM (observer frame) | FRB sightlines | random sightlines |", "|---|---|---|",
              *[f"| > {t} pc cm^-3 | {(per.dm_fg > t).mean():.2f} | {(controls.dm_fg > t).mean():.2f} |"
                for t in (0, 10, 50, 100)],
              f"| median | {per.dm_fg.median():.0f} | {controls.dm_fg.median():.0f} |", "",
              f"N = {n} FRB sightlines, {len(controls)} random sightlines (0.75 deg from each host)."]
    (FG / "SIGHTLINES.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
