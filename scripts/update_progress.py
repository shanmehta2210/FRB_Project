#!/usr/bin/env python3
"""Refresh cutout validation and merge into new_hosts_master.csv via consolidate script."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table

REPO = Path(__file__).resolve().parents[1]
MASTER_CSV = REPO / "pipeline_scripts" / "new_hosts_master.csv"
CUTOUT = REPO / "large_cutouts"
OUT_PIPE = REPO / "pipeline_scripts" / "Output"
VALIDATION = CUTOUT / "cutout_validation.csv"

FLUX_MED_OK = (1e-6, 5.0)
INV_FRAC_OK = 0.3
INV_MED_OK = 1.0


def paired_frbs() -> set[str]:
    out = set()
    for p in CUTOUT.glob("*_flux.fits"):
        frb = p.stem.replace("_flux", "")
        if (CUTOUT / f"{frb}_invvar.fits").is_file():
            out.add(frb)
    return out


def validate_cutout(frb: str) -> dict:
    row = {"frb": frb, "cutout_ok": False, "flux_med": None, "inv_good_frac": None, "inv_med": None, "issues": ""}
    try:
        with fits.open(CUTOUT / f"{frb}_flux.fits") as hf:
            flux = np.squeeze(hf[0].data).astype(float)
        with fits.open(CUTOUT / f"{frb}_invvar.fits") as hi:
            inv = np.squeeze(hi[0].data).astype(float)
        row["flux_med"] = float(np.nanmedian(flux))
        good = inv > 0
        row["inv_good_frac"] = float(np.mean(good))
        row["inv_med"] = float(np.nanmedian(inv[good])) if np.any(good) else 0.0
        issues = []
        if flux.shape != (2290, 2290):
            issues.append("shape")
        if not (FLUX_MED_OK[0] <= row["flux_med"] <= FLUX_MED_OK[1]):
            issues.append("flux_scale")
        if row["inv_good_frac"] < INV_FRAC_OK:
            issues.append("invvar_sparse")
        if row["inv_med"] < INV_MED_OK:
            issues.append("invvar_low")
        row["issues"] = ";".join(issues)
        row["cutout_ok"] = len(issues) == 0
    except Exception as exc:
        row["issues"] = str(exc)
    return row


def _read_posterior_summary(out_dir: Path, frb: str, row: dict) -> None:
    try:
        post = pd.read_csv(out_dir / "astropath_posteriors.csv")
        if len(post) and "posterior_O" in post.columns:
            row["P_O_max"] = float(post["posterior_O"].max())
        loc = pd.read_csv(REPO / "master_frb_localization.csv")
        r = loc.loc[loc.frb == frb].iloc[0]
        host = SkyCoord(float(r.ra_deg), float(r.dec_deg), unit="deg")
        if "ra_deg" in post.columns:
            best = post.sort_values("posterior_O", ascending=False).iloc[0]
            ap = SkyCoord(float(best.ra_deg), float(best.dec_deg), unit="deg")
            row["host_sep_arcsec"] = float(host.separation(ap).arcsec)
    except Exception:
        pass


def pipeline_status(frb: str) -> dict:
    import re

    out_dir = OUT_PIPE / f"{frb}_all"
    row = {
        "frb": frb,
        "pipeline": "missing",
        "P_O_max": None,
        "host_sep_arcsec": None,
        "nearest_src_sep_arcsec": None,
    }

    log_path = out_dir / "master_run.log"
    if log_path.is_file():
        log = log_path.read_text(encoding="utf-8", errors="replace")
        if "exceeds --max-host-sep" in log:
            row["pipeline"] = "host_missing"
            m = re.search(
                r"Nearest SExtractor source \(#\d+\) is ([0-9.]+)\" from",
                log,
            )
            if m:
                row["nearest_src_sep_arcsec"] = float(m.group(1))
            _read_posterior_summary(out_dir, frb, row)
            return row
        if "212.13 arcsec away" in log:
            row["pipeline"] = "bad_host_pick"
            _read_posterior_summary(out_dir, frb, row)
            return row

    if not (out_dir / "host_cutout.fits").exists():
        if (out_dir / "astropath_posteriors.csv").exists():
            _read_posterior_summary(out_dir, frb, row)
        return row

    row["pipeline"] = "partial"
    if (out_dir / "galfit_results.png").exists() and (
        out_dir / "astropath_posteriors.csv"
    ).exists():
        row["pipeline"] = "complete"
    _read_posterior_summary(out_dir, frb, row)
    return row


def main():
    if not MASTER_CSV.is_file():
        raise SystemExit(f"Missing {MASTER_CSV}; run consolidate_new_hosts_logs.py first.")
    names = pd.read_csv(MASTER_CSV)["frb"].astype(str).tolist()
    have = paired_frbs()
    val_rows = [validate_cutout(f) for f in names]
    val_df = pd.DataFrame(val_rows)
    val_df.to_csv(VALIDATION, index=False)

    n_have = len([f for f in names if f in have])
    n_ok = int(val_df["cutout_ok"].sum())
    print(f"Wrote {VALIDATION} ({n_ok}/{len(names)} pass validation)")
    print(f"cutouts on disk {n_have}/{len(names)}")
    print("Run: python scripts/consolidate_new_hosts_logs.py  (refresh master CSV/MD)")


if __name__ == "__main__":
    main()
