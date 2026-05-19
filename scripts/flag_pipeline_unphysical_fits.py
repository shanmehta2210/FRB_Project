"""
Flag pipeline GALFIT results that look unphysical for manual review.

Reads pipeline_galfit_results.csv, enriches with fit.log / feedme checks, and writes
pipeline_unphysical_fits_review.csv at the repo root.

Run from repo root:
    python scripts/flag_pipeline_unphysical_fits.py
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from scripts.galfit_fitlog_parse import parse_fitlog_sky_level  # noqa: E402
OUTPUT_ROOT = REPO_ROOT / "pipeline_scripts" / "Output"
RESULTS_CSV = REPO_ROOT / "pipeline_galfit_results.csv"
OUT_CSV = REPO_ROOT / "pipeline_unphysical_fits_review.csv"

MAG_BRIGHT = 15.0  # AB mag at J)=22.5 — brighter than this is suspicious for hosts
MAG_FAINT = 25.0
CHI2_HIGH = 5.0
RE_CEILING = 99.0
N_CEILING = 5.95
BA_FLOOR = 0.12


def _load_sky_audit(odir: Path) -> dict:
    path = odir / "sky_fit_audit.json"
    if not path.is_file():
        return {}
    try:
        import json

        data = json.loads(path.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except (json.JSONDecodeError, OSError):
        return {}


def _parse_log_summary(log_path: Path) -> dict:
    txt = log_path.read_text(encoding="utf-8", errors="replace")
    blocks = [b for b in txt.split("-------------") if "sersic" in b and "Chi^2/nu" in b]
    last = blocks[-1] if blocks else ""
    sersic_lines = [ln for ln in last.splitlines() if ln.strip().startswith("sersic")]
    n_sersic = len(sersic_lines) // 2 if sersic_lines else 0
    chi_m = re.search(r"Chi\^2/nu\s*=\s*([\d.]+)", last)
    chi2nu = float(chi_m.group(1)) if chi_m else np.nan
    sky_parsed = parse_fitlog_sky_level(log_path)
    sky = float(sky_parsed) if sky_parsed is not None else np.nan
    host_mag = host_re = np.nan
    if sersic_lines:
        clean = sersic_lines[0].replace("(", " ").replace(")", " ").replace(",", " ")
        parts = clean.split()
        if len(parts) >= 9 and parts[0] == "sersic":
            host_mag = float(parts[4])
            host_re = float(parts[5])
    feedme = log_path.parent / "galfit.feedme"
    init_mag = np.nan
    if feedme.is_file():
        m = re.search(r"^\s*3\)\s*([-\d.]+)", feedme.read_text(errors="replace"), re.M)
        if m:
            init_mag = float(m.group(1))
    return {
        "n_sersic": n_sersic,
        "chi2nu_log": chi2nu,
        "sky": sky,
        "host_mag_log": host_mag,
        "host_re_log": host_re,
        "feedme_init_mag": init_mag,
    }


def _classify_row(
    r: pd.Series, log_extra: dict, sky_audit: dict | None = None
) -> tuple[str, list[str], str]:
    flags: list[str] = []
    notes: list[str] = []

    if pd.notna(r.get("mag")) and r["mag"] < MAG_BRIGHT:
        flags.append("mag_too_bright")
        notes.append(f"mag={r['mag']:.2f} (typical hosts ~17–24 at J)=22.5)")
    if pd.notna(r.get("mag")) and r["mag"] > MAG_FAINT:
        flags.append("mag_too_faint")
    if pd.notna(r.get("re")) and r["re"] >= RE_CEILING:
        flags.append("re_at_ceiling")
        notes.append(f"Re={r['re']:.1f} px (constraint max 100)")
    if pd.notna(r.get("n")) and r["n"] >= N_CEILING:
        flags.append("n_at_ceiling")
        notes.append(f"n={r['n']:.2f} at constraint max 6")
    if pd.notna(r.get("b_a")) and r["b_a"] <= BA_FLOOR:
        flags.append("b_a_floor")
        notes.append(f"b/a={r['b_a']:.3f} → inclination pinned at 90°")
    if pd.notna(r.get("chi2nu")) and r["chi2nu"] > CHI2_HIGH:
        flags.append("high_chi2nu")
        notes.append(f"chi2/nu={r['chi2nu']:.2f}")
    if pd.notna(r.get("inc")) and r["inc"] >= 89.5:
        flags.append("inc_face_on")
    if str(r.get("parse_strategy", "")) != "last_sane_single_sersic":
        flags.append("non_sane_parser_block")
        notes.append(f"parser used {r.get('parse_strategy')}")
    if pd.notna(r.get("mag_err")) and r["mag_err"] > 0.3:
        flags.append("large_mag_uncert")
    if pd.notna(r.get("n_err")) and r["n_err"] > 5:
        flags.append("huge_n_err")
    if pd.notna(r.get("pa_err")) and r["pa_err"] > 40:
        flags.append("huge_pa_err")

    sky = log_extra.get("sky", np.nan)
    if pd.notna(sky) and (sky < -10 or sky > 100):
        flags.append("bad_sky")
        notes.append(f"sky parameter {sky:.4g} in fit.log")

    audit = sky_audit or {}
    if audit and audit.get("sky_check_enabled", True) is not False:
        if audit.get("passed") is False:
            flags.append("sky_qa_failed")
            ref = audit.get("sky_ref_adu")
            final = audit.get("sky_final_adu")
            tol = audit.get("sky_tolerance_adu", 3.0)
            notes.append(
                f"sky QA failed: ref={ref}, final={final}, tol=±{tol} ADU"
            )
        elif audit.get("retried"):
            notes.append("sky QA retried with constrained sky")

    if log_extra.get("n_sersic", 1) > 1:
        flags.append("multi_sersic_final")
        notes.append(f"{log_extra['n_sersic']} Sérsic components in last block (host may not be component 1)")

    init_mag = log_extra.get("feedme_init_mag", np.nan)
    if pd.notna(init_mag) and init_mag < 12:
        flags.append("bright_feedme_seed")
        notes.append(f"feedme init mag={init_mag:.2f}")

    tier = "ok"
    if any(
        f in flags
        for f in (
            "mag_too_bright",
            "re_at_ceiling",
            "high_chi2nu",
            "bad_sky",
            "sky_qa_failed",
            "multi_sersic_final",
        )
    ):
        tier = "A_degenerate"
    elif flags:
        tier = "B_suspect"

    return tier, flags, "; ".join(notes)


def main() -> None:
    df = pd.read_csv(RESULTS_CSV)
    num_cols = [
        "chi2nu", "mag", "mag_err", "re", "re_err", "n", "n_err",
        "b_a", "b_a_err", "pa", "pa_err", "inc", "inc_err",
    ]
    for c in num_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    rows: list[dict] = []
    in_csv: set[str] = set()

    for _, r in df.iterrows():
        frb = str(r["frb"])
        in_csv.add(frb)
        odir = OUTPUT_ROOT / f"{frb}_all"
        log_path = odir / "fit.log"
        log_extra = _parse_log_summary(log_path) if log_path.is_file() else {}
        sky_audit = _load_sky_audit(odir)
        tier, flags, notes = _classify_row(r, log_extra, sky_audit)
        if tier == "ok":
            continue
        rows.append({
            "frb": frb,
            "review_tier": tier,
            "flags": "; ".join(flags),
            "notes": notes,
            "chi2nu": r.get("chi2nu"),
            "mag": r.get("mag"),
            "mag_err": r.get("mag_err"),
            "re": r.get("re"),
            "n": r.get("n"),
            "b_a": r.get("b_a"),
            "inc": r.get("inc"),
            "pa": r.get("pa"),
            "parse_strategy": r.get("parse_strategy"),
            "feedme_init_mag": log_extra.get("feedme_init_mag"),
            "log_n_sersic": log_extra.get("n_sersic"),
            "log_sky": log_extra.get("sky"),
            "sky_qa_passed": sky_audit.get("passed") if sky_audit else np.nan,
            "output_dir": f"pipeline_scripts/Output/{frb}_all",
        })

    # Pipeline outputs not in pipeline_galfit_results.csv (e.g. benchmark exclusions)
    for odir in sorted(OUTPUT_ROOT.glob("*_all")):
        frb = odir.name.replace("_all", "")
        if frb in in_csv:
            continue
        log_path = odir / "fit.log"
        if not log_path.is_file():
            rows.append({
                "frb": frb,
                "review_tier": "C_missing_log",
                "flags": "no_fit_log",
                "notes": "Output folder exists but no fit.log",
                "output_dir": f"pipeline_scripts/Output/{frb}_all",
            })
            continue
        ex = _parse_log_summary(log_path)
        sky_audit = _load_sky_audit(odir)
        notes = ["not in pipeline_galfit_results.csv — parse fit.log directly"]
        flags = ["excluded_from_results_csv"]
        tier = "B_suspect"
        if ex.get("chi2nu_log", 0) > CHI2_HIGH:
            flags.append("high_chi2nu")
            notes.append(f"chi2/nu={ex['chi2nu_log']:.2f}")
        if pd.notna(ex.get("host_mag_log")) and ex["host_mag_log"] < MAG_BRIGHT:
            flags.append("mag_too_bright")
            notes.append(f"first sersic mag={ex['host_mag_log']:.2f}")
        if ex.get("host_re_log", 0) >= RE_CEILING:
            flags.append("re_at_ceiling")
        if pd.notna(ex.get("sky")) and ex["sky"] < -10:
            flags.append("bad_sky")
            notes.append(f"sky={ex['sky']:.4g}")
        if ex.get("n_sersic", 1) > 1:
            flags.append("multi_sersic_final")
            notes.append(f"{ex['n_sersic']} Sérsic in last block")
        if ex.get("sky", 0) < -1000 or ex.get("chi2nu_log", 0) > 20:
            tier = "A_degenerate"
        rows.append({
            "frb": frb,
            "review_tier": tier,
            "flags": "; ".join(flags),
            "notes": "; ".join(notes),
            "chi2nu": ex.get("chi2nu_log"),
            "mag": ex.get("host_mag_log"),
            "re": ex.get("host_re_log"),
            "feedme_init_mag": ex.get("feedme_init_mag"),
            "log_n_sersic": ex.get("n_sersic"),
            "log_sky": ex.get("sky"),
            "output_dir": f"pipeline_scripts/Output/{frb}_all",
        })

    out = pd.DataFrame(rows)
    if not out.empty:
        order = {"A_degenerate": 0, "B_suspect": 1, "C_missing_log": 2}
        out["_sort"] = out["review_tier"].map(order)
        out = out.sort_values(["_sort", "frb"]).drop(columns="_sort")

    out.to_csv(OUT_CSV, index=False)
    n_total = len(df)
    print(f"Wrote {OUT_CSV} ({len(out)} flagged of {n_total} in {RESULTS_CSV.name})")
    if not out.empty:
        print(out["review_tier"].value_counts().to_string())
        print("\nTier A — review first:")
        for frb in out.loc[out["review_tier"] == "A_degenerate", "frb"]:
            note = out.loc[out["frb"] == frb, "notes"].iloc[0]
            print(f"  {frb}: {note}")


if __name__ == "__main__":
    main()
