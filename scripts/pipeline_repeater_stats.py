"""Repeater stats for pipeline GALFIT hosts."""
import csv
from collections import Counter

with open("master_frb_localization.csv", newline="", encoding="utf-8") as f:
    loc = {r["frb"]: r for r in csv.DictReader(f)}
with open("pipeline_galfit_results.csv", newline="", encoding="utf-8") as f:
    pipe = list(csv.DictReader(f))

hosts_pipe = [
    r["frb"]
    for r in pipe
    if loc.get(r["frb"], {}).get("coord_semantics") == "host"
]
print(f"pipeline hosts: {len(hosts_pipe)}")

rep = Counter()
for frb in hosts_pipe:
    v = (loc[frb].get("repeater") or "").strip().lower() or "unknown"
    rep[v] += 1
print("repeater counts:", dict(rep))
known = rep["yes"] + rep["no"]
if known:
    print(f"repeater fraction (known only): {rep['yes']/known*100:.1f}% yes of {known}")

print("\nby survey (pipeline hosts):")
surv = {}
for frb in hosts_pipe:
    s = loc[frb].get("survey") or "unknown"
    r = (loc[frb].get("repeater") or "").strip().lower() or "unknown"
    surv.setdefault(s, Counter())[r] += 1
for s, c in sorted(surv.items()):
    print(f"  {s}: {dict(c)}")

print("\nrepeaters in pipeline sample:")
for frb in sorted(hosts_pipe):
    if loc[frb].get("repeater", "").lower() == "yes":
        print(
            f"  {frb}  z={loc[frb].get('z', '')}  survey={loc[frb].get('survey', '')}"
        )

print("\nunknown repeater in pipeline:")
for frb in sorted(hosts_pipe):
    if not (loc[frb].get("repeater") or "").strip():
        print(f"  {frb}")

# morphology from pipeline if present
num_cols = [c for c in ("n", "b_a", "inc", "mag", "re") if c in pipe[0]]
if num_cols:
    pipe_by_frb = {r["frb"]: r for r in pipe}
    print(f"\nmorphology columns: {num_cols}")
    for label in ("yes", "no"):
        frbs = [f for f in hosts_pipe if loc[f].get("repeater", "").lower() == label]
        print(f"\n  repeater={label} (n={len(frbs)}):")
        for col in num_cols:
            vals = []
            for f in frbs:
                v = pipe_by_frb.get(f, {}).get(col, "")
                if v not in ("", None):
                    try:
                        vals.append(float(v))
                    except ValueError:
                        pass
            if vals:
                import statistics as stats

                print(
                    f"    {col}: median={stats.median(vals):.3f}  "
                    f"mean={stats.mean(vals):.3f}  n={len(vals)}"
                )
