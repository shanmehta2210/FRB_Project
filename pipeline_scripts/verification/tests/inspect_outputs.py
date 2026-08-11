"""Ad-hoc: describe the shape of everything the suite writes."""

import glob
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import vercommon as vc  # noqa: E402

d = os.path.join(vc.PER_HOST_ROOT, "20240210A")

print("--- per_host/<FRB>/ ---")
for p in sorted(glob.glob(d + "/*")):
    tag = "dir" if os.path.isdir(p) else f"{os.path.getsize(p) / 1024:.0f} KB"
    print(f"  {os.path.basename(p):28s} {tag}")

for name in ("fourier_profiles.npz", "isophote_profiles.npz"):
    print(f"\n--- {name} ---")
    z = np.load(os.path.join(d, name))
    for k in sorted(z.files):
        print(f"  {k:30s} {str(z[k].shape):10s} {z[k].dtype}")

print("\n--- staged sky refit dir ---")
for p in sorted(glob.glob(os.path.join(d, "galfit_sky_plus", "*"))):
    print(f"  {os.path.basename(p):28s} {os.path.getsize(p) / 1024:.0f} KB")

print("\n--- per-check JSON ---")
for c in ["chi2", "rff", "fourier", "psf", "mag", "isophote", "sky", "astrophot",
          "visual"]:
    j = json.load(open(os.path.join(d, f"{c}.json")))
    print(f"  {c:10s} {len(j):3d} keys  status={j.get('status')}  "
          f"runtime={j.get('_runtime_s')}s")

print("\n--- metrics CSV columns by check prefix ---")
m = pd.read_csv(os.path.join(vc.TABLES_ROOT, "fit_verification_metrics.csv"),
                dtype={"frb": str})
print(f"  {len(m)} rows x {len(m.columns)} columns")
print("  all columns:")
for i in range(0, len(m.columns), 4):
    print("    " + "  ".join(f"{c:34s}" for c in m.columns[i:i + 4]))

print("\n--- flags CSV columns ---")
f = pd.read_csv(os.path.join(vc.TABLES_ROOT, "fit_verification_flags.csv"))
print("   " + ", ".join(f.columns))

print("\n--- population_summary.json key families ---")
pop = json.load(open(os.path.join(vc.TABLES_ROOT, "population_summary.json")))
print(f"  {len(pop)} keys")
fams = sorted({k.rsplit("_", 1)[0] for k in pop})
for k in fams:
    print("   ", k)
