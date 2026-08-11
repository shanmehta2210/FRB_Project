#!/usr/bin/env python3
"""Build / refresh GTC pipeline-trial cohort records under GTC data/.

Trial cohort: 13 FRBs that need imaging, are not in the production 62 fitted
hosts, pass >=1 GTC night (2026-06-24 .. 2026-07-24), and have paired cutouts.

Run after pipeline trials or QA updates:
    python "GTC data/pipeline_trial/consolidate_trial_cohort.py"
"""

from __future__ import annotations

import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

try:
    from astropy.coordinates import SkyCoord
except ImportError:
    SkyCoord = None  # type: ignore

REPO = Path(__file__).resolve().parents[2]
TRIAL_DIR = Path(__file__).resolve().parent
LOC_CSV = REPO / "master_frb_localization.csv"
NEEDS_CSV = REPO / "GTC data/visibility/master_frb_localization_needs_imaging.csv"
EXCLUDE_CSV = REPO / "GTC data/visibility/exclude_pipeline_fitted_frbs.csv"
GTC_SUMMARY = (
    REPO
    / "GTC data/visibility/summaries/gtc_availability_by_frb_2026-06-24_to_2026-07-24.csv"
)
REGISTRY = REPO / "large_cutouts/cutout_registry.csv"
OUT_PIPE = REPO / "pipeline_scripts/Output"
PREVIEW_DIR = REPO / "GTC data/visibility/preview_cutouts_1arcmin"

COHORT_ID = "gtc_jun2026_pipeline_trial"
CREATED_UTC = "2026-06-24 00:00:00 UTC"

MANIFEST_CSV = TRIAL_DIR / "cohort_manifest.csv"
MANIFEST_MD = TRIAL_DIR / "cohort_manifest.md"
FRB_LIST = TRIAL_DIR / "frb_list.txt"
EXCLUDED_CSV = TRIAL_DIR / "excluded_bad_fits.csv"

FRBS_TRIAL = [
    "20210117A",
    "20210214G",
    "20210809C",
    "20220204A",
    "20220506D",
    "20221116A",
    "20230501A",
    "20230521A",
    "20230521B",
    "20230814B",
    "20230913",
    "20230930A",
    "20240203",
]

MANIFEST_COLS = [
    "cohort_id",
    "frb",
    "ra_deg",
    "dec_deg",
    "coord_semantics",
    "selection_needs_imaging",
    "selection_not_production_fitted",
    "selection_gtc_visible",
    "cutout_source",
    "cutout_status",
    "paired_cutout_on_disk",
    "gtc_nights_pass",
    "gtc_pass_fraction",
    "preview_png",
    "pipeline_trial_status",
    "has_output_dir",
    "has_fit_log",
    "has_galfit_png",
    "P_O_max",
    "fit_disposition",
    "fit_disposition_reason",
    "exclude_from_production",
    "pipeline_output_rel",
    "notes",
    "created_utc",
    "updated_utc",
]


def _ts() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")


def paired(frb: str) -> bool:
    cut = REPO / "large_cutouts"
    return (cut / f"{frb}_flux.fits").is_file() and (
        cut / f"{frb}_invvar.fits"
    ).is_file()


def pipeline_trial_status(frb: str) -> dict:
    out_dir = OUT_PIPE / f"{frb}_all"
    row = {
        "pipeline_trial_status": "not_run",
        "has_output_dir": out_dir.is_dir(),
        "has_fit_log": False,
        "has_galfit_png": False,
        "P_O_max": None,
    }
    if not row["has_output_dir"]:
        return row

    row["has_fit_log"] = (out_dir / "fit.log").is_file()
    row["has_galfit_png"] = (out_dir / "galfit_results.png").is_file()
    row["pipeline_trial_status"] = "partial"

    log_path = out_dir / "master_run.log"
    if log_path.is_file():
        log = log_path.read_text(encoding="utf-8", errors="replace")
        if "exceeds --max-host-sep" in log or "No galaxy (SPREAD" in log:
            row["pipeline_trial_status"] = "host_missing"

    if (out_dir / "astropath_posteriors.csv").is_file() and SkyCoord is not None:
        try:
            post = pd.read_csv(out_dir / "astropath_posteriors.csv")
            if len(post) and "posterior_O" in post.columns:
                row["P_O_max"] = float(post["posterior_O"].max())
        except Exception:
            pass

    if row["has_fit_log"] and row["has_galfit_png"]:
        row["pipeline_trial_status"] = "complete"
    elif row["has_output_dir"] and row["pipeline_trial_status"] == "not_run":
        row["pipeline_trial_status"] = "partial"

    return row


def _clean_note(val) -> str:
    if val is None or (isinstance(val, float) and val != val):
        return ""
    s = str(val).strip()
    return "" if s.lower() == "nan" else s


def load_prior_dispositions() -> dict[str, dict]:
    if not MANIFEST_CSV.is_file():
        return {}
    old = pd.read_csv(MANIFEST_CSV)
    out: dict[str, dict] = {}
    for _, r in old.iterrows():
        frb = str(r["frb"])
        out[frb] = {
            "fit_disposition": str(r.get("fit_disposition", "pending") or "pending"),
            "fit_disposition_reason": _clean_note(r.get("fit_disposition_reason", "")),
            "notes": _clean_note(r.get("notes", "")),
            "created_utc": str(r.get("created_utc", CREATED_UTC) or CREATED_UTC),
        }
    return out


def build_rows() -> list[dict]:
    loc = pd.read_csv(LOC_CSV).set_index("frb")
    needs = set(pd.read_csv(NEEDS_CSV)["frb"])
    excluded_prod = set(pd.read_csv(EXCLUDE_CSV)["frb"])
    gtc = pd.read_csv(GTC_SUMMARY).set_index("frb")
    reg = pd.read_csv(REGISTRY) if REGISTRY.is_file() else pd.DataFrame()
    reg_by = {str(r["frb"]): r for _, r in reg.iterrows()} if len(reg) else {}
    prior = load_prior_dispositions()
    now = _ts()

    rows: list[dict] = []
    for frb in FRBS_TRIAL:
        if frb not in needs:
            raise SystemExit(f"{frb} not in needs_imaging CSV")
        if frb in excluded_prod:
            raise SystemExit(f"{frb} is in production-fitted exclude list")

        lr = loc.loc[frb]
        gr = gtc.loc[frb]
        rr = reg_by.get(frb, {})
        pipe = pipeline_trial_status(frb)
        prev = prior.get(frb, {})

        preview = PREVIEW_DIR / f"{frb}_1arcmin.png"
        out_rel = f"pipeline_scripts/Output/{frb}_all"

        rows.append(
            {
                "cohort_id": COHORT_ID,
                "frb": frb,
                "ra_deg": float(lr["ra_deg"]),
                "dec_deg": float(lr["dec_deg"]),
                "coord_semantics": str(lr.get("coord_semantics", "")),
                "selection_needs_imaging": True,
                "selection_not_production_fitted": True,
                "selection_gtc_visible": int(gr["n_pass_nights"]) >= 1,
                "cutout_source": rr.get("source", ""),
                "cutout_status": rr.get("status", "ok" if paired(frb) else "missing"),
                "paired_cutout_on_disk": paired(frb),
                "gtc_nights_pass": int(gr["n_pass_nights"]),
                "gtc_pass_fraction": float(gr["pass_fraction"]),
                "preview_png": str(preview.relative_to(REPO)).replace("\\", "/")
                if preview.is_file()
                else "",
                **pipe,
                "fit_disposition": prev.get("fit_disposition", "pending"),
                "fit_disposition_reason": prev.get("fit_disposition_reason", ""),
                "exclude_from_production": True,
                "pipeline_output_rel": out_rel if pipe["has_output_dir"] else "",
                "notes": prev.get("notes", ""),
                "created_utc": prev.get("created_utc", CREATED_UTC),
                "updated_utc": now,
            }
        )
    return rows


def write_excluded(df: pd.DataFrame) -> None:
    bad = df[df["fit_disposition"].isin(["bad_fit", "exclude"])].copy()
    cols = [
        "frb",
        "fit_disposition",
        "fit_disposition_reason",
        "pipeline_trial_status",
        "has_fit_log",
        "pipeline_output_rel",
        "notes",
        "updated_utc",
    ]
    if len(bad):
        bad[cols].to_csv(EXCLUDED_CSV, index=False)
    else:
        EXCLUDED_CSV.write_text(
            "frb,fit_disposition,fit_disposition_reason,pipeline_trial_status,"
            "has_fit_log,pipeline_output_rel,notes,updated_utc\n",
            encoding="utf-8",
        )


def write_md(df: pd.DataFrame) -> None:
  n_complete = (df["pipeline_trial_status"] == "complete").sum()
  n_pending = (df["fit_disposition"] == "pending").sum()
  n_bad = df["fit_disposition"].isin(["bad_fit", "exclude"]).sum()

  lines = [
    "# GTC pipeline trial cohort (13 FRBs)",
    "",
    f"Last consolidated: **{df['updated_utc'].iloc[0]}** "
    f"(`python \"GTC data/pipeline_trial/consolidate_trial_cohort.py\"`)",
    "",
    "Machine-readable table: `GTC data/pipeline_trial/cohort_manifest.csv`",
    "",
    "## Purpose",
    "",
    "Archival imaging trial runs for GTC-target FRBs that still need host imaging.",
    "These are **not** part of the production 62-host `pipeline_galfit_results.csv` set.",
    "Mark bad fits in `cohort_manifest.csv` (`fit_disposition=bad_fit`) then re-run",
    "consolidate; removable paths are listed in `excluded_bad_fits.csv`.",
    "",
    "## Selection (frozen)",
    "",
    "1. In `master_frb_localization_needs_imaging.csv` (no prior pipeline GALFIT fit)",
    "2. Not in `visibility/exclude_pipeline_fitted_frbs.csv` (production 62)",
    "3. Passes >=1 rigorous GTC night, 2026-06-24 .. 2026-07-24",
    "4. Paired `large_cutouts/{FRB}_flux.fits` + `_invvar.fits` on disk",
    "",
    "## Run pipeline batch",
    "",
    "```bash",
    "python pipeline_scripts/run_all_frbs.py \\",
    "  --list-file \"GTC data/pipeline_trial/frb_list.txt\" \\",
    "  --include-signal \\",
    "  --skip-existing",
    "```",
    "",
    "(`--include-signal` required: 4/13 have `coord_semantics=signal`.)",
    "",
    "## Summary",
    "",
    f"| Metric | Count |",
    f"|--------|------:|",
    f"| Cohort FRBs | {len(df)} |",
    f"| Paired cutouts | {df['paired_cutout_on_disk'].sum()} |",
    f"| Pipeline complete (fit.log + galfit_results.png) | {n_complete} |",
    f"| Fit disposition pending review | {n_pending} |",
    f"| Marked bad_fit / exclude | {n_bad} |",
    "",
    "## Cohort table",
    "",
    "| FRB | Dec | Semantics | Cutout | GTC nights | Pipeline | Disposition | Preview |",
    "|-----|-----|-----------|--------|------------|----------|-------------|---------|",
  ]
  for _, r in df.iterrows():
    prev = r["preview_png"].split("/")[-1] if r["preview_png"] else "—"
    lines.append(
      f"| {r['frb']} | {r['dec_deg']:+.2f} | {r['coord_semantics']} | "
      f"{r['cutout_source'] or r['cutout_status']} | {int(r['gtc_nights_pass'])} | "
      f"{r['pipeline_trial_status']} | {r['fit_disposition']} | {prev} |"
    )
  lines.append("")
  MANIFEST_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    rows = build_rows()
    df = pd.DataFrame(rows, columns=MANIFEST_COLS)
    df.to_csv(MANIFEST_CSV, index=False)
    FRB_LIST.write_text("\n".join(FRBS_TRIAL) + "\n", encoding="utf-8")
    write_excluded(df)
    write_md(df)
    print(f"Wrote {MANIFEST_CSV} ({len(df)} rows)")
    print(f"Wrote {FRB_LIST}")
    print(f"Wrote {EXCLUDED_CSV}")
    print(f"Wrote {MANIFEST_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
