import csv

master = list(csv.DictReader(open("master_frb_localization.csv", encoding="utf-8")))
sig11 = [
    "20181112A",
    "20200627A",
    "20210214G",
    "20210407E",
    "20210809C",
    "20210912A",
    "20220531A",
    "20230521A",
    "20230718A",
    "20230731A",
    "20231006A",
]
sig_pub = ["20191228A", "20220501C", "20220918A"]
print("=== 11 no-pub signal: ellipse errors ===")
for frb in sig11:
    r = next(x for x in master if x["frb"] == frb)
    maj = (r.get("major_sigma_as") or "").strip()
    minor = (r.get("minor_sigma_as") or "").strip()
    rae = (r.get("ra_err_as") or "").strip()
    dece = (r.get("dec_err_as") or "").strip()
    has_ell = bool(maj and minor)
    has_rae = bool(rae and dece)
    print(
        f"{frb}: maj={maj!r} min={minor!r} ra_err={rae!r} dec_err={dece!r} "
        f"pa={r.get('pa_deg')!r} has_ellipse={has_ell} has_ra_dec_err={has_rae} "
        f"status={r.get('status')!r}"
    )
print()
print("=== signal-but-published-host candidates ===")
for frb in sig_pub:
    r = next(x for x in master if x["frb"] == frb)
    print(
        f"{frb}: sem={r['coord_semantics']} maj={r.get('major_sigma_as')!r} "
        f"src={r.get('repeater_source')!r}"
    )
