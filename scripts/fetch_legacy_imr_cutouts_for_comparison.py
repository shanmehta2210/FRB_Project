import argparse
from pathlib import Path

import pandas as pd
import requests


BASE_URL = "https://www.legacysurvey.org/viewer/cutout.fits"
LAYER_MAP = {
    "image": "ls-dr10",
    "model": "ls-dr10-model",
    "resid": "ls-dr10-resid",
}


def fetch_layer(ra: float, dec: float, layer: str, pixscale: float, size: int) -> bytes:
    params = {
        "ra": f"{ra:.8f}",
        "dec": f"{dec:.8f}",
        "layer": layer,
        "pixscale": f"{pixscale}",
        "size": str(size),
        "invvar": "True",
    }
    r = requests.get(BASE_URL, params=params, timeout=60)
    r.raise_for_status()
    return r.content


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fetch Legacy Survey image/model/residual FITS cutouts for FRBs in comparison table."
    )
    parser.add_argument("--comparison-csv", default="legacy_vs_galfit_two_inclinations.csv")
    parser.add_argument("--master-csv", default="master_frb_summary.csv")
    parser.add_argument("--exclude-types", default="", help="Comma-separated LS types to exclude (e.g., REX)")
    parser.add_argument("--size", type=int, default=128, help="Cutout size in pixels")
    parser.add_argument("--pixscale", type=float, default=0.262, help="Pixel scale (arcsec/pixel)")
    parser.add_argument("--out-dir", default="tools/legacy/imr_fits")
    args = parser.parse_args()

    comp = pd.read_csv(args.comparison_csv)
    master = pd.read_csv(args.master_csv)

    exclude_types = [t.strip().upper() for t in str(args.exclude_types).split(",") if t.strip()]
    if exclude_types and "type_ls" in comp.columns:
        comp = comp[~comp["type_ls"].astype(str).str.upper().isin(exclude_types)].copy()

    keep_frbs = set(comp["FRB"].astype(str))
    m = master[master["FRB"].astype(str).isin(keep_frbs)].copy()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    ok = 0
    fail = 0
    for row in m.itertuples(index=False):
        frb = str(row.FRB)
        ra = float(row.RA_deg)
        dec = float(row.DEC_deg)
        frb_dir = out_dir / frb
        frb_dir.mkdir(parents=True, exist_ok=True)

        print(f"Fetching {frb}...")
        for tag, layer in LAYER_MAP.items():
            out_path = frb_dir / f"{frb}_{tag}.fits"
            try:
                payload = fetch_layer(ra, dec, layer=layer, pixscale=args.pixscale, size=args.size)
                out_path.write_bytes(payload)
                ok += 1
            except Exception as exc:
                print(f"  Failed {tag}: {exc}")
                fail += 1

    print("Done")
    print(f"FRBs requested: {len(m)}")
    print(f"Files downloaded: {ok}")
    print(f"Files failed: {fail}")
    print(f"Output dir: {out_dir}")


if __name__ == "__main__":
    main()
