"""Standalone Gini-M20 morphology pass over production hosts.

Reads existing ``statmorph_results.json`` from each ``Output/<FRB>_all/`` —
does **not** re-run the verification suite or Phase Statmorph.

Classification follows Lotz et al. (2008), ApJ 672, 177
(https://doi.org/10.1086/523659):

    Mergers:   G >  -0.14 M20 + 0.33
    E/S0/Sa:   G <= -0.14 M20 + 0.33  and  G >  0.14 M20 + 0.80
    Sb-Irr:    G <= -0.14 M20 + 0.33  and  G <= 0.14 M20 + 0.80

Also reports ``gini_m20_bulge`` / ``gini_m20_merger`` as written by statmorph
(Rodriguez-Gomez et al. 2019 S/F statistics; positive bulge ⇒ early-type side).
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
OUT_ROOT = os.path.join(REPO, "pipeline_scripts", "Output")
TABLES = os.path.join(HERE, "outputs", "tables")


def lotz2008(gini: float, m20: float) -> str:
    if not (np.isfinite(gini) and np.isfinite(m20)):
        return "unknown"
    merger_line = -0.14 * m20 + 0.33
    early_line = 0.14 * m20 + 0.80
    if gini > merger_line:
        return "merger"
    if gini > early_line:
        return "early"  # E/S0/Sa
    return "late"  # Sb-Irr


def collect(cohort: str = "53") -> pd.DataFrame:
    metrics = pd.read_csv(
        os.path.join(TABLES, "fit_verification_metrics.csv"), dtype={"frb": str}
    )
    if cohort == "53":
        frbs = metrics.loc[metrics["in_53"].astype(str).isin(["True", "true", "1"]), "frb"]
    else:
        frbs = metrics["frb"]
    rows = []
    for frb in frbs.astype(str):
        path = os.path.join(OUT_ROOT, f"{frb}_all", "statmorph_results.json")
        row = {"frb": frb, "statmorph_path": path, "statmorph_ok": False}
        if os.path.isfile(path):
            with open(path, encoding="utf-8") as f:
                sm = json.load(f)
            row.update({k: sm.get(k) for k in (
                "gini", "m20", "concentration", "asymmetry", "smoothness",
                "gini_m20_merger", "gini_m20_bulge", "sn_per_pixel",
                "flag", "flag_sersic",
            )})
            row["statmorph_ok"] = True
            row["lotz2008_class"] = lotz2008(row.get("gini"), row.get("m20"))
            bulge = row.get("gini_m20_bulge")
            if bulge is not None and np.isfinite(float(bulge)):
                row["bulge_side"] = "early" if float(bulge) > 0 else "late"
            else:
                row["bulge_side"] = "unknown"
        else:
            row["lotz2008_class"] = "missing"
            row["bulge_side"] = "missing"
        rows.append(row)
    return pd.DataFrame(rows)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cohort", choices=("53", "all64"), default="53")
    args = ap.parse_args(argv)

    os.makedirs(TABLES, exist_ok=True)
    df = collect(args.cohort)
    out = os.path.join(TABLES, f"gini_m20_{args.cohort}.csv")
    df.to_csv(out, index=False)

    ok = df[df.statmorph_ok]
    print(f"Hosts requested ({args.cohort}): {len(df)}")
    print(f"With statmorph_results.json: {len(ok)}")
    print(f"Missing: {len(df) - len(ok)}")
    if len(ok):
        print("\nLotz+2008 class counts:")
        print(ok["lotz2008_class"].value_counts().to_string())
        print("\nstatmorph gini_m20_bulge side (S>0 early, S<=0 late):")
        print(ok["bulge_side"].value_counts().to_string())
        print("\nGini / M20 summary:")
        for col in ("gini", "m20", "gini_m20_bulge", "gini_m20_merger"):
            v = pd.to_numeric(ok[col], errors="coerce").dropna()
            if len(v):
                print(f"  {col:18s} median={v.median():+.4f}  "
                      f"[p16 {v.quantile(0.16):+.4f}, p84 {v.quantile(0.84):+.4f}]")
        flagged = ok[pd.to_numeric(ok["flag"], errors="coerce").fillna(0) > 0]
        print(f"\nstatmorph flag > 0: {len(flagged)}/{len(ok)}")
        print("\nPer-host (frb, G, M20, Lotz, bulge_side):")
        show = ok[["frb", "gini", "m20", "lotz2008_class", "bulge_side",
                   "gini_m20_bulge", "flag"]].copy()
        show = show.sort_values("frb")
        pd.set_option("display.width", 200)
        pd.set_option("display.max_rows", 80)
        print(show.to_string(index=False, float_format=lambda x: f"{x:+.4f}"))
    print(f"\nWrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
