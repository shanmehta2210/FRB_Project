#!/usr/bin/env python3
"""Build pipeline_scripts/new_hosts_master.{csv,md} from legacy logs and live disk state.

Merges (when present): new_hosts_46.txt, list fragments, pipeline status CSV,
batch log, cutout_registry, cutout_validation, and Output/ diagnostics.

Run after batch runs or cutout changes:
    python scripts/consolidate_new_hosts_logs.py
"""

from __future__ import annotations

import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[1]
PIPE = REPO / "pipeline_scripts"
CUTOUT = REPO / "large_cutouts"
OUT_PIPE = PIPE / "Output"
LOC = REPO / "master_frb_localization.csv"
REGISTRY = CUTOUT / "cutout_registry.csv"
VALIDATION = CUTOUT / "cutout_validation.csv"

MASTER_CSV = PIPE / "new_hosts_master.csv"
MASTER_MD = PIPE / "new_hosts_master.md"
BATCH_ARCHIVE = PIPE / "new_hosts_batch_log_archive.txt"

# Legacy inputs (read once, then safe to delete after consolidation).
LEGACY = {
    "cohort_46": PIPE / "new_hosts_46.txt",
    "with_cutouts_41": PIPE / "new_hosts_41_with_cutouts.txt",
    "no_cutout": PIPE / "new_hosts_no_cutout.txt",
    "high_north_ps1": PIPE / "new_hosts_high_north_ps1.txt",
    "former_skymapper": PIPE / "new_hosts_skymapper.txt",
    "pipeline_status": PIPE / "new_hosts_pipeline_status.csv",
    "batch_log": PIPE / "new_hosts_pipeline_batch.log",
}

NO_COVERAGE_DEFAULT = frozenset(
    {"20230930A", "20230125D", "20230718A", "20201123A", "20230731A"}
)
NO_COVERAGE_NOTE = (
    "No Legacy/PS1/DES coverage at host position; manual cutout required."
)


def _ts() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")


def _read_name_list(path: Path) -> list[str]:
    if not path.is_file():
        return []
    names = []
    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        names.append(ln.split()[0])
    return names


def _read_no_cutout_notes(path: Path) -> dict[str, str]:
    notes: dict[str, str] = {}
    if not path.is_file():
        return notes
    for ln in path.read_text(encoding="utf-8", errors="replace").splitlines():
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        frb = ln.split()[0]
        notes[frb] = ln
    return notes


def paired(frb: str) -> bool:
    return (CUTOUT / f"{frb}_flux.fits").is_file() and (
        CUTOUT / f"{frb}_invvar.fits"
    ).is_file()


def _read_md_batch_appendix() -> str:
    if not MASTER_MD.is_file():
        return ""
    text = MASTER_MD.read_text(encoding="utf-8", errors="replace")
    m = re.search(r"## Full batch run log.*?\n\n```\n(.*?)```", text, re.S)
    if not m:
        return ""
    body = m.group(1).strip()
    if body.startswith("(no batch log"):
        return ""
    return body


def parse_batch_log(path: Path) -> tuple[dict[str, dict], str]:
    """Return per-FRB batch fields and full log text for MD appendix."""
    text = ""
    if path.is_file():
        text = path.read_text(encoding="utf-8", errors="replace")
        BATCH_ARCHIVE.write_text(text, encoding="utf-8")
    elif BATCH_ARCHIVE.is_file():
        text = BATCH_ARCHIVE.read_text(encoding="utf-8", errors="replace")
    else:
        text = _read_md_batch_appendix()
    if not text:
        return {}, ""
    per: dict[str, dict] = {}

    line_re = re.compile(
        r"^\[ *\d+/\d+\] (\S+) \(ra=([0-9.+-]+), dec=([0-9.+-]+);.*?\.\.\. "
        r"(OK|FAIL)(?: \(([0-9]+)s, P=([^,]*), ref=([^)]*)\))?",
        re.M,
    )
    for m in line_re.finditer(text):
        frb, _ra, _dec, status, elapsed, p_o, ref = m.groups()
        per[frb] = {
            "batch_line_status": status.lower(),
            "batch_elapsed_s": int(elapsed) if elapsed else None,
            "batch_log_P_O": p_o.strip() if p_o else "",
            "batch_log_ref_catalog": ref.strip() if ref else "",
        }

    table_re = re.compile(
        r"^(\S+)\s+(OK|FAIL)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(.*)$",
        re.M,
    )
    for m in table_re.finditer(text):
        frb, status, p_o, d_host, zp, n_cal, ref, note = m.groups()
        row = per.setdefault(frb, {})
        row.update(
            {
                "batch_table_status": status,
                "batch_table_P_O": p_o,
                "batch_table_d_host": d_host,
                "batch_table_zp_aper_40px": zp,
                "batch_table_n_cal_stars": n_cal,
                "batch_table_ref_catalog": ref,
                "batch_table_note": note.strip(),
            }
        )
    return per, text


def pipeline_status(frb: str) -> dict:
    """Live scan of Output/<frb>_all/ (same logic as update_progress.pipeline_status)."""
    from astropy.coordinates import SkyCoord

    out_dir = OUT_PIPE / f"{frb}_all"
    row = {
        "pipeline": "missing",
        "P_O_max": None,
        "host_sep_arcsec": None,
        "nearest_src_sep_arcsec": None,
        "has_output_dir": out_dir.is_dir(),
        "has_galfit": False,
        "has_posteriors": False,
        "has_zero_points": False,
    }

    log_path = out_dir / "master_run.log"
    if log_path.is_file():
        log = log_path.read_text(encoding="utf-8", errors="replace")
        if "exceeds --max-host-sep" in log or "No galaxy (SPREAD" in log:
            row["pipeline"] = "host_missing"
            m = re.search(
                r"Nearest SExtractor source \(#\d+\) is ([0-9.]+)\" from", log
            )
            if m:
                row["nearest_src_sep_arcsec"] = float(m.group(1))
        if "212.13 arcsec away" in log:
            row["pipeline"] = "bad_host_pick"

    if (out_dir / "astropath_posteriors.csv").is_file():
        row["has_posteriors"] = True
    if (out_dir / "zero_points.json").is_file():
        row["has_zero_points"] = True
    if (out_dir / "galfit_results.png").is_file():
        row["has_galfit"] = True

    if not (out_dir / "host_cutout.fits").is_file():
        if row["pipeline"] == "missing" and row["has_posteriors"]:
            row["pipeline"] = "partial"
        try:
            post = pd.read_csv(out_dir / "astropath_posteriors.csv")
            if len(post) and "posterior_O" in post.columns:
                row["P_O_max"] = float(post["posterior_O"].max())
            loc = pd.read_csv(LOC)
            r = loc.loc[loc.frb == frb].iloc[0]
            host = SkyCoord(float(r.ra_deg), float(r.dec_deg), unit="deg")
            if "ra" in post.columns:
                best = post.sort_values("posterior_O", ascending=False).iloc[0]
                ap = SkyCoord(float(best.ra), float(best.dec), unit="deg")
                row["host_sep_arcsec"] = float(host.separation(ap).arcsec)
        except Exception:
            pass
        return row

    row["pipeline"] = "partial"
    if row["has_galfit"] and row["has_posteriors"]:
        row["pipeline"] = "complete"
    try:
        post = pd.read_csv(out_dir / "astropath_posteriors.csv")
        if len(post) and "posterior_O" in post.columns:
            row["P_O_max"] = float(post["posterior_O"].max())
        loc = pd.read_csv(LOC)
        r = loc.loc[loc.frb == frb].iloc[0]
        host = SkyCoord(float(r.ra_deg), float(r.dec_deg), unit="deg")
        if "ra" in post.columns:
            best = post.sort_values("posterior_O", ascending=False).iloc[0]
            ap = SkyCoord(float(best.ra), float(best.dec), unit="deg")
            row["host_sep_arcsec"] = float(host.separation(ap).arcsec)
    except Exception:
        pass
    return row


def _reconstruct_batch_text_from_df(df: pd.DataFrame) -> str:
    """Rebuild a readable batch summary when the raw .log file is gone."""
    lines = [
        "[batch] Reconstructed summary from new_hosts_master.csv "
        "(original new_hosts_pipeline_batch.log was consolidated then removed).",
        "",
    ]
    for _, r in df.iterrows():
        if str(r.get("batch_exit", "")) in ("", "not_in_batch", "skipped_no_cutout"):
            continue
        frb = r["frb"]
        elapsed = r.get("batch_elapsed_s", "")
        zp = r.get("batch_zp_aper_40px", "")
        ncal = r.get("batch_n_cal_stars", "")
        ref = r.get("batch_ref_catalog", "")
        p_o = r.get("P_O_max", "")
        lines.append(
            f"{frb}: {r.get('batch_exit', '')} elapsed={elapsed}s "
            f"P_O_max={p_o} zp_40={zp} n_cal={ncal} ref={ref}"
        )
    lines.append("")
    lines.append(f"[batch] {len(lines)-3} FRB(s) with batch metadata in master CSV.")
    return "\n".join(lines)


def load_registry() -> pd.DataFrame:
    if REGISTRY.is_file():
        return pd.read_csv(REGISTRY)
    return pd.DataFrame()


def load_validation() -> pd.DataFrame:
    if VALIDATION.is_file():
        return pd.read_csv(VALIDATION)
    return pd.DataFrame()


def build_master() -> tuple[pd.DataFrame, str]:
    cohort = _read_name_list(LEGACY["cohort_46"])
    if not cohort and MASTER_CSV.is_file():
        cohort = pd.read_csv(MASTER_CSV)["frb"].astype(str).tolist()
    if not cohort:
        raise SystemExit("No cohort FRB list (new_hosts_46.txt or existing master CSV).")

    loc = pd.read_csv(LOC)
    set_41 = set(_read_name_list(LEGACY["with_cutouts_41"]))
    set_high_north = set(_read_name_list(LEGACY["high_north_ps1"]))
    set_skymapper = set(_read_name_list(LEGACY["former_skymapper"]))
    if not set_high_north and MASTER_CSV.is_file():
        old = pd.read_csv(MASTER_CSV)
        if "list_tags" in old.columns:
            for _, r in old.iterrows():
                tags = str(r.get("list_tags", ""))
                if "high_north_ps1" in tags:
                    set_high_north.add(str(r["frb"]))
                if "former_skymapper" in tags:
                    set_skymapper.add(str(r["frb"]))
    no_cutout_notes = _read_no_cutout_notes(LEGACY["no_cutout"])
    no_coverage = set(no_cutout_notes) | NO_COVERAGE_DEFAULT

    batch_per, batch_text = parse_batch_log(LEGACY["batch_log"])
    reg = load_registry()
    val = load_validation()

    legacy_status = {}
    if LEGACY["pipeline_status"].is_file():
        for _, r in pd.read_csv(LEGACY["pipeline_status"]).iterrows():
            legacy_status[str(r["frb"])] = r.to_dict()

    rows: list[dict] = []
    for frb in cohort:
        lm = loc.loc[loc["frb"] == frb]
        ra = float(lm.iloc[0]["ra_deg"]) if len(lm) else None
        dec = float(lm.iloc[0]["dec_deg"]) if len(lm) else None

        reg_row = reg.loc[reg["frb"] == frb].iloc[0] if len(reg) and frb in reg["frb"].values else None
        val_row = val.loc[val["frb"] == frb].iloc[0] if len(val) and frb in val["frb"].values else None

        on_disk = paired(frb)
        if frb in no_coverage:
            cutout_status = "no_coverage"
        elif on_disk:
            cutout_status = "ok"
        elif reg_row is not None and pd.notna(reg_row.get("status")):
            cutout_status = str(reg_row["status"])
        else:
            cutout_status = "missing"

        tags = []
        if frb in set_high_north:
            tags.append("high_north_ps1")
        if frb in set_skymapper:
            tags.append("former_skymapper")
        if frb in no_coverage:
            tags.append("no_coverage")

        live = pipeline_status(frb)
        batch = batch_per.get(frb, {})
        leg = legacy_status.get(frb, {})

        if cutout_status == "no_coverage":
            batch_exit = "skipped_no_cutout"
        elif frb in batch_per:
            batch_exit = batch.get("batch_line_status") or batch.get("batch_table_status", "")
        else:
            batch_exit = leg.get("batch_exit", "not_in_batch")

        note_parts = []
        if frb in no_cutout_notes:
            note_parts.append(no_cutout_notes[frb])
        if reg_row is not None and pd.notna(reg_row.get("notes")):
            note_parts.append(str(reg_row["notes"]))
        if frb in set_skymapper:
            note_parts.append("Former SkyMapper cutout host; re-fetch via cutout_download.py if needed.")

        rows.append(
            {
                "frb": frb,
                "ra_deg": ra,
                "dec_deg": dec,
                "in_cohort_46": True,
                "in_batch_41_list": (
                    frb in set_41
                    if set_41
                    else (on_disk and cutout_status == "ok")
                ),
                "paired_cutout_on_disk": on_disk,
                "list_tags": ";".join(tags),
                "cutout_status": cutout_status,
                "cutout_source": reg_row["source"] if reg_row is not None else "",
                "cutout_layer": reg_row["layer"] if reg_row is not None else "",
                "cutout_resampled": reg_row["resampled"] if reg_row is not None else "",
                "cutout_ok": val_row["cutout_ok"] if val_row is not None else "",
                "flux_median": val_row["flux_med"] if val_row is not None else "",
                "invvar_good_frac": val_row["inv_good_frac"] if val_row is not None else "",
                "validation_issues": val_row["issues"] if val_row is not None else "",
                "pipeline": live["pipeline"],
                "P_O_max": live["P_O_max"],
                "host_sep_arcsec": live["host_sep_arcsec"],
                "nearest_src_sep_arcsec": live["nearest_src_sep_arcsec"],
                "has_output_dir": live["has_output_dir"],
                "has_galfit": live["has_galfit"],
                "has_posteriors": live["has_posteriors"],
                "has_zero_points": live["has_zero_points"],
                "batch_exit": batch_exit,
                "batch_elapsed_s": batch.get("batch_elapsed_s"),
                "batch_zp_aper_40px": batch.get("batch_table_zp_aper_40px", ""),
                "batch_n_cal_stars": batch.get("batch_table_n_cal_stars", ""),
                "batch_ref_catalog": batch.get("batch_table_ref_catalog")
                or batch.get("batch_log_ref_catalog", ""),
                "notes": " | ".join(p for p in note_parts if p and str(p) != "nan"),
                "updated_utc": _ts(),
            }
        )

    df = pd.DataFrame(rows).sort_values("frb")
    return df, batch_text


def write_markdown(df: pd.DataFrame, batch_text: str) -> None:
    ts = _ts()
    n = len(df)
    n_ok = int((df["cutout_status"] == "ok").sum())
    n_no = int((df["cutout_status"] == "no_coverage").sum())
    n_complete = int((df["pipeline"] == "complete").sum())
    n_host_miss = int((df["pipeline"] == "host_missing").sum())
    n_partial = int((df["pipeline"] == "partial").sum())

    lines = [
        "# New-host cohort — master log (46 FRBs)",
        "",
        f"Last consolidated: **{ts}** (`python scripts/consolidate_new_hosts_logs.py`)",
        "",
        "Machine-readable table: `pipeline_scripts/new_hosts_master.csv`",
        "",
        "This file replaces scattered `new_hosts_*.txt`, `new_hosts_pipeline_status.csv`,",
        "`new_hosts_pipeline_batch.log`, and `large_cutouts/PROGRESS.md` for the 46-host cohort.",
        "Cutout inventory for all on-disk cutouts (including non-cohort): `large_cutouts/cutout_registry.csv`.",
        "",
        "## Survey / fetch policy",
        "",
        "- r-band 10′ cutouts (2290 px @ 0.262″/px): **Legacy → PS1 → DES** (no SkyMapper).",
        "- One FRB at a time: `python scripts/cutout_download.py <FRB>`",
        "- Batch pipeline: `python pipeline_scripts/run_all_frbs.py --use-localization-host --list-file pipeline_scripts/new_hosts_master.csv`",
        "  (use column `frb`; or filter CSV to `paired_cutout_on_disk == True`).",
        "",
        "## Summary",
        "",
        "| Metric | Count |",
        "|--------|------:|",
        f"| Cohort FRBs | {n} |",
        f"| Paired flux+invvar on disk | {n_ok} |",
        f"| No survey coverage (manual cutout) | {n_no} |",
        f"| Pipeline complete | {n_complete} |",
        f"| Pipeline partial | {n_partial} |",
        f"| Pipeline host_missing | {n_host_miss} |",
        "",
        "## Cohort table (key columns)",
        "",
        "| FRB | Dec | Cutout | Pipeline | P(O) max | Nearest src ″ | Batch | Notes |",
        "|-----|-----|--------|----------|---------:|--------------:|-------|-------|",
    ]

    for _, r in df.iterrows():
        p_o = r["P_O_max"]
        p_str = f"{p_o:.3f}" if pd.notna(p_o) else "—"
        near = r["nearest_src_sep_arcsec"]
        near_str = f"{near:.1f}" if pd.notna(near) else "—"
        note = str(r["notes"])[:60] + ("…" if len(str(r["notes"])) > 60 else "")
        lines.append(
            f"| {r['frb']} | {r['dec_deg']:+.1f} | {r['cutout_status']} | {r['pipeline']} | "
            f"{p_str} | {near_str} | {r['batch_exit']} | {note} |"
        )

    lines.extend(
        [
            "",
            "## List tags (from former sidecar files)",
            "",
            "- **high_north_ps1** — Dec ≳ +70°; PS1-only footprint for calibration.",
            "- **former_skymapper** — originally SkyMapper cutouts; re-download with Legacy ladder.",
            "- **no_coverage** — no Legacy/PS1/DES at host; see notes column.",
            "",
            "## No-cutout FRBs (detail)",
            "",
        ]
    )
    for _, r in df.loc[df["cutout_status"] == "no_coverage"].iterrows():
        lines.append(f"- **{r['frb']}** — {r['notes'] or NO_COVERAGE_NOTE}")

    lines.extend(
        [
            "",
        "## Batch run log",
        "",
        "When `new_hosts_pipeline_batch.log` exists, consolidate copies it to "
        "`new_hosts_batch_log_archive.txt` and embeds the full text below. "
        "Historical 41-host batch (May 2026, `--use-localization-host`): 41 FRBs, "
        "2112 s total, 0 failed — per-FRB timing and ZP details were in the deleted "
        "sidecar log; re-run batch or check git history if you need the verbatim file.",
        "",
            "```",
            batch_text.rstrip()
            if batch_text.strip()
            else (
                "(Batch log not available — run a new batch and keep "
                "new_hosts_pipeline_batch.log until consolidate copies it to "
                "new_hosts_batch_log_archive.txt.)"
            ),
            "```",
            "",
        ]
    )
    MASTER_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    df, batch_text = build_master()
    if not batch_text.strip():
        batch_text = _reconstruct_batch_text_from_df(df)
    df.to_csv(MASTER_CSV, index=False)
    write_markdown(df, batch_text)
    print(f"Wrote {MASTER_CSV} ({len(df)} rows)")
    print(f"Wrote {MASTER_MD}")
    print(df.groupby(["cutout_status", "pipeline"]).size().to_string())


if __name__ == "__main__":
    main()
