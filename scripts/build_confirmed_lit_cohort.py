"""Build the confirmed-literature-host production cohort.

An FRB qualifies when it already has a production GALFIT fit
(``pipeline_galfit_results.csv``) *and* its host association is backed by a
**published** localization paper. Hosts whose only citation is an in-prep
manuscript (Verdi+2025) or that rest on pipeline AstroPath alone are dropped —
see ``pipeline_scripts/docs/WEAK_ASSOCIATIONS_PRODUCTION67.md``.

Writes ``production_confirmed_lit_hosts.csv`` (``--list-file`` compatible).
"""
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
RESULTS = REPO / "pipeline_galfit_results.csv"
AUDIT = REPO / "production67_lit_astropath_audit.csv"
LOC = REPO / "master_frb_localization.csv"
OUT = REPO / "production_confirmed_lit_hosts.csv"

IN_PREP_MARKERS = ("in prep", "in prep.", "in preparation")


def has_published_cite(text) -> bool:
    if not isinstance(text, str) or not text.strip():
        return False
    return not any(m in text.lower() for m in IN_PREP_MARKERS)


def main() -> int:
    results = pd.read_csv(RESULTS, dtype={"frb": str})
    audit = pd.read_csv(AUDIT, dtype={"frb": str}).set_index("frb")
    loc = pd.read_csv(LOC, dtype={"frb": str}).set_index("frb")

    rows = []
    for frb in results["frb"]:
        if frb in audit.index:
            published = bool(audit.at[frb, "has_published_lit_cite"])
            citation = audit.at[frb, "citation"]
            source = "production67_lit_astropath_audit.csv"
        else:
            citation = loc.at[frb, "repeater_source"] if frb in loc.index else None
            published = has_published_cite(citation)
            source = "master_frb_localization.csv"
        rows.append(
            {
                "frb": frb,
                "confirmed_lit_host": published,
                "citation": citation,
                "cite_source": source,
                "coord_semantics": loc.at[frb, "coord_semantics"] if frb in loc.index else None,
            }
        )

    df = pd.DataFrame(rows)
    keep = df[df["confirmed_lit_host"]].drop(columns=["confirmed_lit_host"])
    dropped = df[~df["confirmed_lit_host"]]

    keep.to_csv(OUT, index=False)
    print(f"[cohort] production results : {len(df)}")
    print(f"[cohort] confirmed-lit hosts: {len(keep)}  -> {OUT.name}")
    print(f"[cohort] dropped            : {len(dropped)}")
    for _, r in dropped.iterrows():
        print(f"          - {r['frb']:<11} {r['coord_semantics'] or '?':<7} {r['citation'] or 'no citation'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
