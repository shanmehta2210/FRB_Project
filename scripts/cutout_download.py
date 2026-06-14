#!/usr/bin/env python3
"""Download one r-band 10' cutout at a time (Legacy -> PS1 -> DES).

Updates a single registry: large_cutouts/cutout_registry.csv

Examples:
  python cutout_download.py 20180301A
  python cutout_download.py 20180301A --force
  python cutout_download.py --scan          # rebuild registry from files on disk
  python cutout_download.py --list          # print registry
"""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from astropy.io import fits

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from cutout_fetch_common import (
    DES_LAYERS,
    LS_LAYERS_GLOBAL,
    PS1_DEC_MIN,
    PS1_PIXSCALE,
    fetch_des_pair,
    fetch_legacy_pair,
    fetch_ps1_pair,
    legacy_layers_for_dec,
    preflight_coverage,
    probe_des,
)
from cutout_resample import write_standardized

REPO = Path(__file__).resolve().parents[1]
CUTOUT_DIR = REPO / "large_cutouts"
LOC_CSV = REPO / "master_frb_localization.csv"
REGISTRY = CUTOUT_DIR / "cutout_registry.csv"
AUDIT_CSV = CUTOUT_DIR / "coverage_audit.csv"

REGISTRY_COLS = [
    "frb",
    "ra_deg",
    "dec_deg",
    "source",
    "layer",
    "resampled",
    "downloaded_utc",
    "status",
]


def _ts() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")


def paired_on_disk() -> set[str]:
    """Only count standard {FRB}_flux.fits + {FRB}_invvar.fits pairs."""
    import re

    pat = re.compile(r"^(\d{4}[0-9A-Za-z]+)_(flux|invvar)\.fits$")
    flux_frbs = set()
    for p in CUTOUT_DIR.glob("*_flux.fits"):
        m = pat.match(p.name)
        if m and (CUTOUT_DIR / f"{m.group(1)}_invvar.fits").is_file():
            flux_frbs.add(m.group(1))
    return flux_frbs


def _print_preflight(pf: dict[str, bool], dec: float) -> None:
    print("  Preflight:")
    print(f"    Legacy global     : {'yes' if pf['legacy_global'] else 'no'}")
    print(f"    Legacy hemisphere : {'yes' if pf['legacy_hemisphere'] else 'no'} "
          f"({pf.get('hemisphere_layer', '')})")
    if dec <= PS1_DEC_MIN:
        print(f"    PS1               : no (dec <= {PS1_DEC_MIN})")
    else:
        print(f"    PS1               : {'yes' if pf['ps1'] else 'no'}")
    print(f"    DES               : {'yes' if pf['des'] else 'no'}")


def _audit_fast_path(frb: str, ra: float, dec: float, allow_ps1: bool) -> dict[str, bool] | None:
    """When audit already rules out Legacy/PS1, only probe DES (~30s not ~15min)."""
    if not AUDIT_CSV.is_file():
        return None
    m = pd.read_csv(AUDIT_CSV)
    m = m.loc[m["frb"] == frb]
    if m.empty:
        return None
    row = m.iloc[0]
    leg_g = bool(row.get("legacy_global_ok", False))
    leg_h = bool(row.get("legacy_hemisphere_ok", False))
    ps1_ok = bool(row.get("ps1_ok", False)) and allow_ps1 and dec > PS1_DEC_MIN
    if leg_g or leg_h or ps1_ok:
        return None
    print("  (audit: Legacy + PS1 unavailable — probing DES only)")
    des_ok = probe_des(ra, dec)
    hem = legacy_layers_for_dec(dec)
    return {
        "legacy_global": False,
        "legacy_hemisphere": False,
        "ps1": False,
        "des": des_ok,
        "any": des_ok,
        "hemisphere_layer": hem[0],
    }


def load_registry() -> pd.DataFrame:
    if REGISTRY.is_file():
        df = pd.read_csv(REGISTRY)
        for c in REGISTRY_COLS:
            if c not in df.columns:
                df[c] = ""
        return df[REGISTRY_COLS]
    return pd.DataFrame(columns=REGISTRY_COLS)


def save_registry(df: pd.DataFrame) -> None:
    CUTOUT_DIR.mkdir(parents=True, exist_ok=True)
    df.sort_values("frb").to_csv(REGISTRY, index=False)


def upsert_row(
    df: pd.DataFrame,
    frb: str,
    ra: float,
    dec: float,
    *,
    source: str,
    layer: str,
    resampled: bool,
    status: str,
) -> pd.DataFrame:
    row = {
        "frb": frb,
        "ra_deg": ra,
        "dec_deg": dec,
        "source": source,
        "layer": layer,
        "resampled": resampled,
        "downloaded_utc": _ts(),
        "status": status,
    }
    df = df[df["frb"] != frb]
    return pd.concat([df, pd.DataFrame([row])], ignore_index=True)


def remove_cutout(frb: str) -> None:
    for suffix in ("_flux.fits", "_invvar.fits"):
        p = CUTOUT_DIR / f"{frb}{suffix}"
        if p.is_file():
            p.unlink()


def _save_native(frb, flux, fh, inv, ih, source, layer) -> None:
    fits.writeto(CUTOUT_DIR / f"{frb}_flux.fits", flux, fh, overwrite=True)
    fits.writeto(CUTOUT_DIR / f"{frb}_invvar.fits", inv, ih, overwrite=True)


def _save_ps1(frb, ra, dec, flux, fh, inv, ih) -> None:
    write_standardized(
        str(CUTOUT_DIR / f"{frb}_flux.fits"),
        str(CUTOUT_DIR / f"{frb}_invvar.fits"),
        flux, fh, inv, ih, ra, dec, PS1_PIXSCALE,
        center_on_array=False,
    )


def fetch_frb(
    frb: str,
    ra: float,
    dec: float,
    *,
    force: bool = False,
    allow_ps1: bool = True,
    skip_preflight: bool = False,
) -> tuple[bool, str, str, bool]:
    """Return (ok, source, layer, resampled)."""
    if force:
        remove_cutout(frb)
    elif frb in paired_on_disk():
        print(f"{frb}: cutouts already on disk (use --force to replace)")
        return True, "", "", False

    print(f"{frb}  ra={ra:.4f} dec={dec:.4f}")

    if skip_preflight:
        pf = {
            "legacy_global": True,
            "legacy_hemisphere": True,
            "ps1": allow_ps1 and dec > PS1_DEC_MIN,
            "des": True,
            "any": True,
            "hemisphere_layer": legacy_layers_for_dec(dec)[0],
        }
    else:
        pf = _audit_fast_path(frb, ra, dec, allow_ps1) or preflight_coverage(
            ra, dec, allow_ps1=allow_ps1
        )
    _print_preflight(pf, dec)

    if not pf["any"]:
        print("  Abort: no survey has coverage at this position.")
        return False, "", "", False

    if pf["legacy_global"]:
        print("  Legacy global (full) ...", end=" ", flush=True)
        flux, fh, inv, ih, layer = fetch_legacy_pair(ra, dec, LS_LAYERS_GLOBAL)
        if flux is not None:
            _save_native(frb, flux, fh, inv, ih, "legacy", layer)
            print("OK")
            return True, "legacy", layer, False
        print("fail (probe passed but full cutout failed)")

    elif pf["legacy_hemisphere"]:
        hem = legacy_layers_for_dec(dec)
        print(f"  Legacy {hem[0]} (full) ...", end=" ", flush=True)
        flux, fh, inv, ih, layer = fetch_legacy_pair(ra, dec, hem)
        if flux is not None:
            _save_native(frb, flux, fh, inv, ih, "legacy", layer)
            print("OK")
            return True, "legacy", layer, False
        print("fail (probe passed but full cutout failed)")
    else:
        print("  Legacy: skip (preflight)")

    if pf["ps1"]:
        print("  PS1 r (full, slow) ...", end=" ", flush=True)
        flux, fh, inv, ih = fetch_ps1_pair(ra, dec)
        if flux is not None and inv is not None:
            _save_ps1(frb, ra, dec, flux, fh, inv, ih)
            print("OK")
            return True, "ps1", "ps1", True
        print("fail (probe passed but full cutout failed)")
    elif dec <= PS1_DEC_MIN:
        print(f"  PS1: skip (dec <= {PS1_DEC_MIN})")
    else:
        print("  PS1: skip (preflight)")

    if pf["des"]:
        print(f"  DES {DES_LAYERS[0]} (full) ...", end=" ", flush=True)
        flux, fh, inv, ih, layer = fetch_des_pair(ra, dec)
        if flux is not None:
            _save_native(frb, flux, fh, inv, ih, "des", layer)
            print("OK")
            return True, "des", layer, False
        print("fail (probe passed but full cutout failed)")
    else:
        print("  DES: skip (preflight)")

    print("  Abort: all probed tiers failed on full download.")
    return False, "", "", False


def scan_registry() -> pd.DataFrame:
    """Rebuild registry rows for every paired cutout on disk."""
    loc = pd.read_csv(LOC_CSV)
    have = paired_on_disk()
    old = load_registry()
    known = {r["frb"]: r for r in old.to_dict("records")} if len(old) else {}

    rows = []
    for frb in sorted(have):
        m = loc.loc[loc["frb"] == frb]
        ra = float(m.iloc[0]["ra_deg"]) if len(m) else ""
        dec = float(m.iloc[0]["dec_deg"]) if len(m) else ""
        prev = known.get(frb, {})
        rows.append(
            {
                "frb": frb,
                "ra_deg": ra,
                "dec_deg": dec,
                "source": prev.get("source", ""),
                "layer": prev.get("layer", ""),
                "resampled": prev.get("resampled", ""),
                "downloaded_utc": prev.get("downloaded_utc", ""),
                "status": "on_disk",
            }
        )
    df = pd.DataFrame(rows, columns=REGISTRY_COLS)
    save_registry(df)
    return df


def cmd_list(df: pd.DataFrame) -> None:
    loc = pd.read_csv(LOC_CSV)
    have = paired_on_disk()
    print(f"Registry: {REGISTRY}")
    print(f"Paired on disk: {len(have)}")
    if len(df):
        print(df.to_string(index=False))
    missing = sorted(set(loc["frb"]) - have)
    if missing:
        print(f"\nIn localization CSV but no cutout ({len(missing)}):")
        for frb in missing[:20]:
            print(f"  {frb}")
        if len(missing) > 20:
            print(f"  ... +{len(missing) - 20} more")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("frb", nargs="?", help="FRB name to download")
    parser.add_argument("--force", action="store_true", help="Replace existing files")
    parser.add_argument("--no-ps1", action="store_true", help="Skip PS1 tier")
    parser.add_argument("--scan", action="store_true", help="Rebuild registry from disk")
    parser.add_argument("--list", action="store_true", help="Show registry")
    parser.add_argument(
        "--preflight-only",
        action="store_true",
        help="Run probes only; do not download",
    )
    parser.add_argument(
        "--no-preflight",
        action="store_true",
        help="Skip probes and try every tier (slow)",
    )
    args = parser.parse_args()

    df = load_registry()

    if args.scan:
        df = scan_registry()
        print(f"Scanned {len(df)} paired cutout(s) -> {REGISTRY}")
        return

    if args.list:
        cmd_list(df)
        return

    if not args.frb:
        parser.error("provide FRB name, or use --scan / --list")

    loc = pd.read_csv(LOC_CSV)
    match = loc.loc[loc["frb"] == args.frb]
    if match.empty:
        raise SystemExit(f"{args.frb!r} not in {LOC_CSV}")

    row = match.iloc[0]
    ra, dec = float(row["ra_deg"]), float(row["dec_deg"])

    if args.preflight_only:
        pf = _audit_fast_path(args.frb, ra, dec, not args.no_ps1) or preflight_coverage(
            ra, dec, allow_ps1=not args.no_ps1
        )
        _print_preflight(pf, dec)
        raise SystemExit(0 if pf["any"] else 1)

    ok, source, layer, resampled = fetch_frb(
        args.frb,
        ra,
        dec,
        force=args.force,
        allow_ps1=not args.no_ps1,
        skip_preflight=args.no_preflight,
    )

    if ok and args.frb in paired_on_disk():
        df = upsert_row(
            df, args.frb, ra, dec,
            source=source or "existing",
            layer=layer,
            resampled=resampled,
            status="ok",
        )
    else:
        df = upsert_row(
            df, args.frb, ra, dec,
            source="none",
            layer="",
            resampled=False,
            status="failed",
        )

    save_registry(df)
    if not ok:
        raise SystemExit(1)
    print(f"Updated {REGISTRY}")


if __name__ == "__main__":
    main()
