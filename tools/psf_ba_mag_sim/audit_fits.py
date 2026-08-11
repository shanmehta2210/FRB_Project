#!/usr/bin/env python3
"""Diagnose GALFIT fit quality across the simulation grid."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from sim_utils import TOOL_DIR, ensure_output_layout, load_config

_DCHI_RE = re.compile(r"dChi2/Chi2:\s*([-+eE0-9.]+)")


def parse_fitlog_dchi(log_path: Path) -> float | None:
    if not log_path.is_file():
        return None
    text = log_path.read_text(encoding="utf-8", errors="replace")
    for line in text.splitlines():
        if "Iteration : 1" in line and "dChi2/Chi2" in line:
            m = _DCHI_RE.search(line)
            if m:
                try:
                    return float(m.group(1))
                except ValueError:
                    return None
    return None


def log_size(log_path: Path) -> int:
    try:
        return log_path.stat().st_size if log_path.is_file() else 0
    except OSError:
        return 0


def build_audit_table(cfg: dict, layout: dict[str, Path]) -> pd.DataFrame:
    fit = pd.read_csv(layout["catalogs"] / "fit_results.csv")
    truth = pd.read_csv(layout["catalogs"] / "truth_catalog.csv")
    m = fit.merge(truth, on=["galaxy_id", "realization"], how="left")

    for col in ("ba_fit", "mag_fit", "re_fit_pix", "chi2nu", "mag_true", "ba_true", "re_arcsec_true"):
        m[col] = pd.to_numeric(m[col], errors="coerce")
    m["converged"] = m["converged"].astype(str).str.lower().eq("true")

    rows = []
    for _, r in m.iterrows():
        log_path = layout["fits"] / r["galaxy_id"] / r["mode"] / "fit.log"
        dchi1 = parse_fitlog_dchi(log_path)
        mag_err = r["mag_fit"] - r["mag_true"] if pd.notna(r["mag_fit"]) else np.nan
        ba_err = r["ba_fit"] - r["ba_true"] if pd.notna(r["ba_fit"]) else np.nan

        numeric_blowup = dchi1 is not None and abs(dchi1) > 1e6
        mag_stuck_25 = pd.notna(r["mag_fit"]) and abs(r["mag_fit"] - 25.0) < 0.05
        mag_off = pd.notna(mag_err) and abs(mag_err) > 1.5
        ba_collapse = pd.notna(r["ba_fit"]) and r["ba_fit"] < 0.15
        empty_log = log_size(log_path) < 200

        if empty_log:
            status = "crash_empty_log"
        elif ba_collapse and not bool(r.get("converged")):
            status = "ba_collapse"
        elif mag_off and not bool(r.get("converged")):
            status = "mag_degenerate"
        elif bool(r.get("converged")):
            status = "ok"
        elif numeric_blowup:
            status = "numeric_blowup_iter1"
        elif r.get("converged"):
            status = "ok"
        else:
            status = "parse_failed"

        rows.append(
            {
                **r.to_dict(),
                "dchi2_iter1": dchi1,
                "mag_err": mag_err,
                "ba_err": ba_err,
                "log_bytes": log_size(log_path),
                "numeric_blowup": numeric_blowup,
                "mag_stuck_near_25": mag_stuck_25,
                "quality_status": status,
            }
        )
    return pd.DataFrame(rows)


def write_report(audit: pd.DataFrame, out_path: Path) -> None:
    lines = [
        "# GALFIT simulation fit audit",
        "",
        "## Summary",
        "",
        f"- Total fit rows: **{len(audit)}** (expected 270)",
        f"- Quality `ok`: **{(audit['quality_status'] == 'ok').sum()}**",
        f"- Numeric blow-up (iter 1 |dChi2/Chi2| > 1e6): **{audit['numeric_blowup'].sum()}**",
        f"- Magnitude stuck near 25.0 (sky-dominated degeneracy): **{audit['mag_stuck_near_25'].sum()}**",
        f"- |mag_fit - mag_true| > 1.5 mag: **{(audit['mag_err'].abs() > 1.5).sum()}**",
        f"- Empty / crashed fit.log: **{(audit['quality_status'] == 'crash_empty_log').sum()}**",
        f"- b/a collapse (fit < 0.15): **{(audit['quality_status'] == 'ba_collapse').sum()}**",
        "",
        "## Root cause (v1 — fixed in v2)",
        "",
        "v1 used ZP-tuned sky ~26.5 e-/pix, custom sigma, and mag constraints.",
        "v2: sky-subtracted stamps, C) none, ZP 22.5 + coadd, sky ±~1 ADU constraint.",
        "",
        "## Failures by mode",
        "",
        audit.groupby("mode")["quality_status"].value_counts().to_string(),
        "",
        "## Failures vs magnitude (bright = low mag)",
        "",
    ]

    for mode in ("psf", "nopsf"):
        sub = audit[audit["mode"] == mode]
        lines.append(f"### {mode}")
        tbl = (
            sub.groupby(pd.cut(sub["mag_true"], [16.5, 18, 20, 22, 24.5]))["quality_status"]
            .value_counts()
            .unstack(fill_value=0)
        )
        lines.append(tbl.to_string())
        lines.append("")

    lines.extend(
        [
            "## Recommended simulation fixes (applied in v2)",
            "",
            "1. **Fixed zeropoint** 22.5 + **coadd_exptime** scaling (not ZP tuning).",
            "2. **Sky-subtracted stamps** for GALFIT; sky component seeded at 0, still free.",
            "3. **C) none** — do not feed ideal Poisson sigma (breaks mag on sky-sub images).",
            "4. **No `1 mag` constraint** in constraints.txt (caused mag snap to ~25).",
            "5. **Quality gates**: |dmag| < 1.5, sane b/a, non-empty fit.log.",
            "",
            "Re-run: `generate_mocks.py` → `run_galfit_grid.py` → `analyze.py` → `audit_fits.py`.",
            "",
        ]
    )

    out_path.write_text("\n".join(lines), encoding="utf-8")


def plot_audit(audit: pd.DataFrame, plots_dir: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    for ax, mode, title in zip(
        axes.flat[:2],
        ["psf", "nopsf"],
        ["PSF fit", "No PSF"],
    ):
        sub = audit[(audit["mode"] == mode) & audit["mag_err"].notna()]
        sc = ax.scatter(
            sub["mag_true"],
            sub["mag_err"],
            c=sub["ba_err"],
            cmap="coolwarm",
            vmin=-0.6,
            vmax=0.6,
            s=12,
            alpha=0.7,
        )
        ax.axhline(0, color="k", lw=0.8)
        ax.axhline(8, color="r", ls="--", lw=0.8, label="±8 mag (stuck-at-sky)")
        ax.axhline(-8, color="r", ls="--", lw=0.8)
        ax.set_xlabel(r"True $m_r$")
        ax.set_ylabel(r"$\Delta$mag (fit $-$ true)")
        ax.set_title(title)
        ax.invert_xaxis()
        ax.grid(True, alpha=0.3)
        plt.colorbar(sc, ax=ax, label=r"$\Delta(b/a)$")

    ax = axes[1, 0]
    status_counts = audit.groupby(["mode", "quality_status"]).size().unstack(fill_value=0)
    status_counts.plot(kind="bar", ax=ax, rot=45)
    ax.set_title("Quality status counts")
    ax.set_ylabel("N fits")

    ax = axes[1, 1]
    ok = audit[audit["quality_status"] == "ok"]
    bad = audit[audit["quality_status"] != "ok"]
    ax.scatter(ok["mag_true"], ok["ba_fit"], s=8, alpha=0.4, label="ok", c="C0")
    ax.scatter(bad["mag_true"], bad["ba_fit"], s=8, alpha=0.6, label="bad", c="C3")
    ax.set_xlabel(r"True $m_r$")
    ax.set_ylabel(r"Reported $b/a$")
    ax.invert_xaxis()
    ax.legend()
    ax.set_title("All fits: reported b/a vs true mag")
    ax.grid(True, alpha=0.3)

    fig.suptitle("GALFIT simulation quality audit", fontsize=13)
    fig.tight_layout()
    fig.savefig(plots_dir / "fit_quality_audit.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=TOOL_DIR / "config.yaml")
    args = parser.parse_args()

    cfg = load_config(args.config)
    layout = ensure_output_layout(cfg)
    reports = layout["root"] / "reports"
    reports.mkdir(exist_ok=True)

    audit = build_audit_table(cfg, layout)
    audit.to_csv(reports / "fit_audit.csv", index=False)
    write_report(audit, reports / "FIT_AUDIT.md")
    plot_audit(audit, layout["plots"])

    print(f"Wrote {reports / 'fit_audit.csv'}")
    print(f"Wrote {reports / 'FIT_AUDIT.md'}")
    print(f"Wrote {layout['plots'] / 'fit_quality_audit.png'}")
    print("\nStatus counts:")
    print(audit["quality_status"].value_counts().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
