"""Audit every production `Output/<FRB>_all/` against the confirmed cohort.

Reports, per host, whether it is in the confirmed-literature cohort, whether it
appears in `pipeline_galfit_results.csv`, and which Phase 3a neighbor policy its
`cutout_meta.json` records. Anything that is not both *in the cohort* and
*re_separation* is listed as unconfirmed.
"""
import glob
import json
import os

import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXPECTED_POLICY = "re_separation"


def main() -> int:
    cohort = set(
        pd.read_csv(os.path.join(REPO, "production_confirmed_lit_hosts.csv"),
                    dtype={"frb": str}).frb
    )
    results = set(
        pd.read_csv(os.path.join(REPO, "pipeline_galfit_results.csv"),
                    dtype={"frb": str}).frb
    )
    dirs = sorted(glob.glob(os.path.join(REPO, "pipeline_scripts", "Output", "*_all")))

    print(f"folders={len(dirs)}  results_rows={len(results)}  cohort={len(cohort)}\n")
    header = f"{'FRB':<12} {'cohort':<7} {'results':<8} {'policy':<15} {'fit.log':<8}"
    print(header)
    print("-" * len(header))

    unconfirmed = []
    for d in dirs:
        frb = os.path.basename(d)[: -len("_all")]
        meta_path = os.path.join(d, "cutout_meta.json")
        if os.path.isfile(meta_path):
            with open(meta_path, encoding="utf-8") as fh:
                policy = json.load(fh).get("neighbor_policy") or "null"
        else:
            policy = "NO META"
        has_fit = os.path.isfile(os.path.join(d, "fit.log"))
        ok = frb in cohort and policy == EXPECTED_POLICY
        if ok:
            continue
        unconfirmed.append(frb)
        print(f"{frb:<12} {str(frb in cohort):<7} {str(frb in results):<8} "
              f"{policy:<15} {str(has_fit):<8}")

    print(f"\nconfirmed (cohort + {EXPECTED_POLICY}): {len(dirs) - len(unconfirmed)}")
    print(f"unconfirmed                        : {len(unconfirmed)}")
    if unconfirmed:
        print("  " + ", ".join(unconfirmed))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
