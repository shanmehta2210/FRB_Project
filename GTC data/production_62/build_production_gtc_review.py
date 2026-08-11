#!/usr/bin/env python3
"""Tiered GTC science review for the 62 production pipeline hosts.

Reads pipeline_galfit_results.csv, astropath_posteriors.csv, pipeline_unphysical_fits_review.csv,
and the Jun–Jul 2026 GTC visibility rollup. Writes reports under GTC data/production_62/.

No pipeline code changes. Refresh after new fits or visibility scans:
    python "GTC data/production_62/build_production_gtc_review.py"
"""

from __future__ import annotations

import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from astropy.coordinates import SkyCoord

REPO = Path(__file__).resolve().parents[2]
OUT_DIR = Path(__file__).resolve().parent

PROD_CSV = REPO / "pipeline_galfit_results.csv"
LOC_CSV = REPO / "master_frb_localization.csv"
UNPHYS_CSV = REPO / "pipeline_unphysical_fits_review.csv"
CONFIDENT_HOSTS_TXT = REPO / "Archive/notes/new_confident_hosts.txt"
NEW_HOSTS_MASTER = REPO / "pipeline_scripts/new_hosts_master.csv"
GTC_SUMMARY = (
    REPO
    / "GTC data/visibility/summaries/gtc_availability_by_frb_2026-06-24_to_2026-07-24.csv"
)
OUT_PIPE = REPO / "pipeline_scripts/Output"

CANDIDATES_CSV = OUT_DIR / "gtc_science_candidates.csv"
CANDIDATES_MD = OUT_DIR / "gtc_science_candidates.md"
VISIBLE_CSV = OUT_DIR / "gtc_visible_intersection.csv"
LIT_HOSTS_CSV = OUT_DIR / "literature_confident_hosts.csv"

_CONFIDENT_ROW_RE = re.compile(
    r"^\s*(\d{8}[A-Z]?)\s*&\s*"
    r"([-\d.]+)\s*&\s*"
    r"([-\d.]+)\s*&\s*"
    r"([\d.]+)\s*&\s*"
    r"([\d.]+)\s*&\s*"
    r"([\d.]+)\s*&\s*"
    r"([\d.]+)\s*&\s*"
    r"([\d.]+)\s*&\s*"
    r"([^&\\]+)",
    re.IGNORECASE,
)

P_O_THRESH = 0.85
HOST_OFFSET_ARCSEC = 2.0
SNR_DEPTH_THRESH = 5.0
MAG_FAINT = 25.0
LIT_P_HOST_THRESH = 0.85
ASSOC_REASONS = frozenset({"low_P_O", "host_offset_gt2as"})

MORPH_FLAGS = frozenset(
    {"b_a_floor", "inc_face_on", "n_at_ceiling", "re_at_ceiling"}
)
DEPTH_FLAGS = frozenset(
    {"sky_qa_failed", "bad_sky", "mag_too_faint", "high_chi2nu"}
)


def _ts() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")


def _flag_set(flags_str: str) -> set[str]:
    if not flags_str or str(flags_str).lower() == "nan":
        return set()
    return {f.strip() for f in str(flags_str).split(";") if f.strip()}


def _load_literature_confident_hosts() -> pd.DataFrame:
    """Parse Archive/notes/new_confident_hosts.txt (Sharma/Verdi host table)."""
    if not CONFIDENT_HOSTS_TXT.is_file():
        return pd.DataFrame(
            columns=["frb", "lit_P_host", "survey", "sharma_nature_ref", "in_cohort_46"]
        )
    text = CONFIDENT_HOSTS_TXT.read_text(encoding="utf-8", errors="replace")
    rows: list[dict] = []
    for line in text.splitlines():
        m = _CONFIDENT_ROW_RE.match(line)
        if not m:
            continue
        frb = m.group(1)
        tail = line.strip()
        rows.append(
            {
                "frb": frb,
                "lit_P_host": float(m.group(4)),
                "survey": m.group(9).strip(),
                "sharma_nature_ref": bool(
                    re.search(r"635|964\.\.131S|Sharma", tail, re.I)
                ),
            }
        )
    df = pd.DataFrame(rows)
    if not len(df):
        return df
    cohort_46: set[str] = set()
    if NEW_HOSTS_MASTER.is_file():
        nh = pd.read_csv(NEW_HOSTS_MASTER)
        if "in_cohort_46" in nh.columns:
            cohort_46 = set(nh.loc[nh["in_cohort_46"] == True, "frb"].astype(str))
    df["in_cohort_46"] = df["frb"].isin(cohort_46)
    return df


def _galfit_ok_for_assoc(review_tier: str, flags: set[str]) -> bool:
    """Archival fit acceptable for trusting literature host over AstroPath mismatch."""
    if review_tier == "ok":
        return True
    if review_tier == "A_degenerate":
        return False
    if review_tier == "B_suspect":
        return True
    return False


def _apply_lit_host_downgrade(
    tiers: list[str],
    reasons: list[str],
    lit_p: float | None,
    galfit_ok: bool,
) -> tuple[list[str], list[str], bool]:
    """Drop Tier A when literature P_host is secure and archival GALFIT is fine."""
    if (
        "A" not in tiers
        or lit_p is None
        or lit_p < LIT_P_HOST_THRESH
        or not galfit_ok
    ):
        return tiers, reasons, False
    new_tiers = [t for t in tiers if t != "A"]
    new_reasons = [r for r in reasons if r not in ASSOC_REASONS]
    return new_tiers, new_reasons, True


def _astropath_metrics(frb: str, loc: pd.DataFrame) -> tuple[float | None, float | None]:
    post = OUT_PIPE / f"{frb}_all" / "astropath_posteriors.csv"
    if not post.is_file():
        return None, None
    p = pd.read_csv(post)
    if not len(p) or "posterior_O" not in p.columns:
        return None, None
    p_o = float(p["posterior_O"].max())
    best = p.loc[p["posterior_O"].idxmax()]
    if frb not in loc.index:
        return p_o, None
    host = SkyCoord(
        float(loc.loc[frb, "ra_deg"]),
        float(loc.loc[frb, "dec_deg"]),
        unit="deg",
    )
    ap = SkyCoord(float(best["ra"]), float(best["dec"]), unit="deg")
    sep = float(host.separation(ap).arcsec)
    return p_o, sep


def _tier_a(p_o: float | None, sep: float | None) -> tuple[bool, list[str]]:
    reasons: list[str] = []
    if p_o is not None and p_o < P_O_THRESH:
        reasons.append("low_P_O")
    if sep is not None and sep > HOST_OFFSET_ARCSEC:
        reasons.append("host_offset_gt2as")
    return bool(reasons), reasons


def _tier_b(
    tier: str,
    flags: set[str],
    n_sersic: int,
    p_o: float | None,
    sep: float | None,
    snr_win: float | None,
    mag: float | None,
) -> tuple[bool, list[str]]:
    secure_assoc = (
        p_o is not None
        and p_o >= 0.9
        and (sep is None or sep < 1.0)
    )
    if flags & MORPH_FLAGS:
        return False, []
    if "multi_sersic_final" in flags and n_sersic > 1 and secure_assoc:
        return False, []

    reasons: list[str] = []
    if tier == "A_degenerate":
        if flags & DEPTH_FLAGS:
            reasons.append("degenerate_depth_snr")
        elif "multi_sersic_final" in flags and not secure_assoc:
            reasons.append("degenerate_multi_sersic")
        elif not flags:
            reasons.append("degenerate_unspecified")
        else:
            return False, []

    elif tier == "B_suspect":
        if "sky_qa_failed" in flags or "bad_sky" in flags:
            reasons.append("sky_calibration")
        if mag is not None and mag > MAG_FAINT:
            reasons.append("mag_too_faint")
        if (
            "high_chi2nu" in flags
            and n_sersic <= 1
            and snr_win is not None
            and snr_win < SNR_DEPTH_THRESH
        ):
            reasons.append("high_chi2nu_low_snr")
        if not reasons:
            return False, []

    else:
        return False, []

    return bool(reasons), reasons


def _tier_c_note(flags: set[str], n_sersic: int, p_o: float | None, sep: float | None) -> str:
    notes: list[str] = []
    if flags & {"b_a_floor", "inc_face_on"}:
        notes.append("face_on_pin")
    secure = p_o is not None and p_o >= 0.9 and (sep is None or sep < 1.0)
    if n_sersic > 1 and secure and "multi_sersic_final" in flags:
        notes.append("multi_sersic_secure_assoc")
    return ";".join(notes)


def build_rows(lit_hosts: pd.DataFrame) -> pd.DataFrame:
    prod = pd.read_csv(PROD_CSV)
    loc = pd.read_csv(LOC_CSV).set_index("frb")
    vis = pd.read_csv(GTC_SUMMARY).set_index("frb")
    unphys = pd.read_csv(UNPHYS_CSV) if UNPHYS_CSV.is_file() else pd.DataFrame()
    unphys_by = {str(r["frb"]): r for _, r in unphys.iterrows()} if len(unphys) else {}
    lit_by = {str(r["frb"]): r for _, r in lit_hosts.iterrows()} if len(lit_hosts) else {}

    rows: list[dict] = []
    for _, r in prod.iterrows():
        frb = str(r["frb"])
        p_o, sep = _astropath_metrics(frb, loc)
        u = unphys_by.get(frb, {})
        review_tier = str(u.get("review_tier", "ok") or "ok")
        flags = _flag_set(str(u.get("flags", "")))

        n_sersic = int(r.get("n_sersic_components", 1))
        snr = float(r["snr_win"]) if pd.notna(r.get("snr_win")) else None
        mag = float(r["mag"]) if pd.notna(r.get("mag")) else None

        a_ok, a_reasons = _tier_a(p_o, sep)
        b_ok, b_reasons = _tier_b(review_tier, flags, n_sersic, p_o, sep, snr, mag)

        tiers: list[str] = []
        all_reasons: list[str] = []
        if a_ok:
            tiers.append("A")
            all_reasons.extend(a_reasons)
        if b_ok:
            tiers.append("B")
            all_reasons.extend(b_reasons)
        if not tiers:
            continue

        lit = lit_by.get(frb)
        lit_p = float(lit["lit_P_host"]) if lit is not None else None
        in_lit = lit is not None
        in_c46 = bool(lit["in_cohort_46"]) if lit is not None else False
        galfit_ok = _galfit_ok_for_assoc(review_tier, flags)
        tiers, all_reasons, downgraded = _apply_lit_host_downgrade(
            tiers, all_reasons, lit_p, galfit_ok
        )
        if not tiers:
            continue

        gr = vis.loc[frb] if frb in vis.index else None
        nights = int(gr["n_pass_nights"]) if gr is not None else 0
        pass_frac = float(gr["pass_fraction"]) if gr is not None else 0.0

        rows.append(
            {
                "frb": frb,
                "gtc_tiers": "+".join(tiers),
                "gtc_reasons": "|".join(all_reasons),
                "P_O_max": p_o,
                "host_offset_arcsec": sep,
                "lit_P_host": lit_p,
                "in_confident_hosts_table": in_lit,
                "in_cohort_46": in_c46,
                "galfit_ok_for_assoc": galfit_ok,
                "gtc_assoc_downgraded": downgraded,
                "mag_ab": mag,
                "chi2nu": float(r["chi2nu"]) if pd.notna(r.get("chi2nu")) else None,
                "inc_deg": float(r["inc"]) if pd.notna(r.get("inc")) else None,
                "b_a": float(r["b_a"]) if pd.notna(r.get("b_a")) else None,
                "n_sersic": n_sersic,
                "snr_win": snr,
                "unphys_tier": review_tier,
                "unphys_flags": ";".join(sorted(flags)) if flags else "",
                "tier_c_note": _tier_c_note(flags, n_sersic, p_o, sep),
                "gtc_n_pass_nights": nights,
                "gtc_pass_fraction": pass_frac,
                "gtc_visible": nights >= 1,
                "coord_semantics": str(loc.loc[frb, "coord_semantics"])
                if frb in loc.index
                else "",
                "dec_deg": float(loc.loc[frb, "dec_deg"]) if frb in loc.index else None,
                "updated_utc": _ts(),
            }
        )

    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values(
            ["gtc_visible", "gtc_n_pass_nights", "gtc_tiers", "frb"],
            ascending=[False, False, True, True],
        )
    return df


def write_md(df: pd.DataFrame, n_prod: int, lit_hosts: pd.DataFrame) -> None:
    n_a = df["gtc_tiers"].str.contains("A").sum() if len(df) else 0
    n_b = df["gtc_tiers"].str.contains("B").sum() if len(df) else 0
    n_vis = int(df["gtc_visible"].sum()) if len(df) else 0
    n_down = int(df["gtc_assoc_downgraded"].sum()) if len(df) and "gtc_assoc_downgraded" in df.columns else 0
    n_c46_lit = int(lit_hosts["in_cohort_46"].sum()) if len(lit_hosts) and "in_cohort_46" in lit_hosts.columns else 0

    vis_all = pd.read_csv(GTC_SUMMARY)
    prod_frbs = set(pd.read_csv(PROD_CSV)["frb"].astype(str))
    prod_vis = vis_all[vis_all["frb"].isin(prod_frbs)]
    n_prod_visible = int((prod_vis["n_pass_nights"] >= 1).sum())

    flagged_frbs = set(df["frb"]) if len(df) else set()
    clean = prod_frbs - flagged_frbs
    clean_vis = prod_vis[prod_vis["frb"].isin(clean) & (prod_vis["n_pass_nights"] >= 1)]

    lines = [
        "# Production-62 GTC science review",
        "",
        f"Last built: **{df['updated_utc'].iloc[0] if len(df) else _ts()}**",
        f"(`python \"GTC data/production_62/build_production_gtc_review.py\"`)",
        "",
        "Machine-readable: `GTC data/production_62/gtc_science_candidates.csv`",
        "",
        "## Summary",
        "",
        f"| Metric | Count |",
        f"|--------|------:|",
        f"| Production hosts (`pipeline_galfit_results.csv`) | {n_prod} |",
        f"| GTC-visible (Jun 24–Jul 24 2026) | {n_prod_visible} |",
        f"| Tier A or B flagged | {len(df)} |",
        f"| — Tier A (association) | {n_a} |",
        f"| — Tier B (depth/SNR degenerate) | {n_b} |",
        f"| Flagged and GTC-visible | {n_vis} |",
        f"| Clean (not flagged) and GTC-visible | {len(clean_vis)} |",
        f"| Literature confident hosts (`new_confident_hosts.txt`) | {len(lit_hosts)} |",
        f"| — overlap 46-host cohort | {n_c46_lit} |",
        f"| Tier A downgraded (lit host + good GALFIT) | {n_down} |",
        "",
        "## Tier definitions",
        "",
        "- **Tier A**: `P(O) < 0.85` and/or AstroPath host > 2″ from published localization.",
        "  Downgraded when FRB is in `Archive/notes/new_confident_hosts.txt` with",
        f"  `P_host ≥ {LIT_P_HOST_THRESH}` and archival GALFIT is acceptable (AP mismatch alone",
        "  is not a GTC driver when Sharma/Verdi association is already secure).",
        "- **Tier B**: Degenerate/suspect archival fit driven by depth/SNR/calibration,",
        "  excluding intrinsic morphology (face-on pins, n/Re ceilings, secure multi-Sérsic blends).",
        "- **Tier C** (notes only): face-on pins, secure multi-Sérsic — not GTC drivers.",
        "",
        "## Tier A + B candidates",
        "",
        "| FRB | Tiers | Reasons | P(O) | Lit P_host | Offset″ | mag | χ²/ν | GTC nights | Visible |",
        "|-----|-------|---------|------|------------|---------|-----|------|------------|---------|",
    ]
    for _, r in df.iterrows():
        p = f"{r['P_O_max']:.3g}" if pd.notna(r["P_O_max"]) else "—"
        lp = f"{r['lit_P_host']:.3g}" if pd.notna(r.get("lit_P_host")) else "—"
        s = f"{r['host_offset_arcsec']:.2f}" if pd.notna(r["host_offset_arcsec"]) else "—"
        m = f"{r['mag_ab']:.2f}" if pd.notna(r["mag_ab"]) else "—"
        c = f"{r['chi2nu']:.2f}" if pd.notna(r["chi2nu"]) else "—"
        vis = "yes" if r["gtc_visible"] else "no"
        lines.append(
            f"| {r['frb']} | {r['gtc_tiers']} | {r['gtc_reasons']} | {p} | {lp} | {s} | "
            f"{m} | {c} | {int(r['gtc_n_pass_nights'])} | {vis} |"
        )

    lines.extend(
        [
            "",
            "## GTC-visible intersection (scheduling list)",
            "",
            "See `gtc_visible_intersection.csv`.",
            "",
        ]
    )
    if n_vis:
        vis_df = df[df["gtc_visible"]]
        for _, r in vis_df.iterrows():
            lines.append(
                f"- **{r['frb']}** [{r['gtc_tiers']}] {r['gtc_reasons']} "
                f"({int(r['gtc_n_pass_nights'])} nights)"
            )
    else:
        lines.append("(none)")

    lines.append("")
    CANDIDATES_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    n_prod = len(pd.read_csv(PROD_CSV))
    lit_hosts = _load_literature_confident_hosts()
    lit_hosts.to_csv(LIT_HOSTS_CSV, index=False)

    df = build_rows(lit_hosts)
    df.to_csv(CANDIDATES_CSV, index=False)

    vis_df = df[df["gtc_visible"]].copy() if len(df) else df
    vis_df.to_csv(VISIBLE_CSV, index=False)

    write_md(df, n_prod, lit_hosts)

    print(f"Wrote {LIT_HOSTS_CSV} ({len(lit_hosts)} literature hosts)")
    print(f"Wrote {CANDIDATES_CSV} ({len(df)} flagged / {n_prod} production)")
    print(f"Wrote {VISIBLE_CSV} ({len(vis_df)} GTC-visible)")
    print(f"Wrote {CANDIDATES_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
