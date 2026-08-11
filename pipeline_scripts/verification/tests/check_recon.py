"""Ad-hoc: how well the analytic rebuild matches GALFIT's model, cohort-wide."""

import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import vercommon as vc  # noqa: E402

old = pd.read_csv(os.path.join(vc.TABLES_ROOT, "fit_verification_metrics.csv"),
                  dtype={"frb": str})
old["_r"] = pd.to_numeric(old["model_recon_max_frac"], errors="coerce")
rows = []
for _, r in old.sort_values("_r", ascending=False).iterrows():
    host = vc.load_host(r.frb)
    res = vc.model_reconstruction_error(host)
    rows.append((r.frb, host.re, host.n, r._r, res["model_recon_max_frac"],
                 res["model_recon_flux_frac"], res["model_recon_oversample"]))
df = pd.DataFrame(rows, columns=["frb", "re", "n", "before", "after",
                                 "flux_err", "k"])
print(df.head(12).to_string(index=False, float_format=lambda v: f"{v:9.5f}"))
print()
print(f"max frac error   before {df.before.max():.4f}  ->  after {df.after.max():.5f}")
print(f"median           before {df.before.median():.4f}  ->  after {df.after.median():.5f}")
print(f"n above 2%       before {(df.before > 0.02).sum()}  ->  after {(df.after > 0.02).sum()}")
print(f"n above 5%       before {(df.before > 0.05).sum()}  ->  after {(df.after > 0.05).sum()}")
print("\nchosen oversampling:", dict(df.k.value_counts().sort_index()))
