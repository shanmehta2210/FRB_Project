import os
import time
from io import BytesIO

import numpy as np
import pandas as pd
import requests
from astropy.io import fits

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
MASTER_CSV = os.path.join(ROOT, "master_frb_summary.csv")
BASE_DIR = os.path.join(ROOT, "psfs", "PSFEx + SExtractor")
RUNS_DIR = os.path.join(BASE_DIR, "runs")

# 10 arcmin = 600 arcsec. Legacy Survey pixel scale = 0.262 arcsec / pixel.
PIXSCALE = 0.262
CUTOUT_SIZE_PX = int(round(600.0 / PIXSCALE))  # ~2290 px
BAND = "r"
LAYERS = ["ls-dr10", "ls-dr9"]


def ensure_structure(frb: str) -> str:
    run_dir = os.path.join(RUNS_DIR, frb)
    subdirs = [
        os.path.join(run_dir, "input", "sextractor"),
        os.path.join(run_dir, "input", "psfex"),
        os.path.join(run_dir, "output", "sextractor"),
        os.path.join(run_dir, "output", "psfex"),
        os.path.join(run_dir, "cutouts"),
        os.path.join(run_dir, "logs"),
    ]
    for d in subdirs:
        os.makedirs(d, exist_ok=True)
    return run_dir


def fetch_legacy(ra: float, dec: float):
    for layer in LAYERS:
        url = (
            "https://www.legacysurvey.org/viewer/cutout.fits?"
            f"ra={ra}&dec={dec}&size={CUTOUT_SIZE_PX}&layer={layer}"
            f"&pixscale={PIXSCALE}&bands={BAND}&invvar=True"
        )
        try:
            resp = requests.get(url, timeout=180)
            if resp.status_code != 200:
                continue
            with fits.open(BytesIO(resp.content)) as hdul:
                flux = hdul[0].data if len(hdul) >= 1 else None
                flux_header = hdul[0].header.copy() if len(hdul) >= 1 else None
                invvar = hdul[1].data if len(hdul) >= 2 else None
                invvar_header = hdul[1].header.copy() if len(hdul) >= 2 else None

                if flux is None or np.all(flux == 0):
                    continue
                return {
                    "layer": layer,
                    "flux": flux,
                    "flux_header": flux_header,
                    "invvar": invvar,
                    "invvar_header": invvar_header,
                }
        except requests.RequestException:
            continue
        except Exception:
            continue
    return None


def fetch_ps1_flux(ra: float, dec: float):
    if dec < -30:
        return None

    filenames_url = (
        "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py?"
        f"ra={ra}&dec={dec}&filters=r&sep=comma"
    )
    try:
        f_resp = requests.get(filenames_url, timeout=60)
        if f_resp.status_code != 200:
            return None
        lines = f_resp.text.strip().split("\n")
        if len(lines) < 2:
            return None
        keys = lines[0].split(",")
        values = lines[1].split(",")
        row = dict(zip(keys, values))
        filename = row.get("filename")
        if not filename:
            return None

        fits_url = (
            "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
            f"ra={ra}&dec={dec}&size={CUTOUT_SIZE_PX}&format=fits&red={filename}"
        )
        resp = requests.get(fits_url, timeout=240)
        if resp.status_code != 200:
            return None
        if b"<html>" in resp.content[:256].lower():
            return None

        with fits.open(BytesIO(resp.content)) as hdul:
            flux = hdul[0].data if len(hdul) >= 1 else None
            flux_header = hdul[0].header.copy() if len(hdul) >= 1 else None
            if flux is None:
                return None
            return {"flux": flux, "flux_header": flux_header}
    except Exception:
        return None


def write_cutouts(frb: str, run_dir: str, product: dict, source: str):
    cutout_dir = os.path.join(run_dir, "cutouts")
    flux_path = os.path.join(cutout_dir, f"{frb}_10arcmin_flux.fits")

    flux_header = product["flux_header"] if product.get("flux_header") is not None else fits.Header()
    flux_header["OBJECT"] = frb
    flux_header["SURVEY"] = source
    flux_header["BAND"] = BAND
    flux_header["CUTTYPE"] = "10arcmin"
    fits.writeto(flux_path, product["flux"], flux_header, overwrite=True)

    invvar_saved = False
    if product.get("invvar") is not None:
        invvar_path = os.path.join(cutout_dir, f"{frb}_10arcmin_invvar.fits")
        invvar_header = product["invvar_header"] if product.get("invvar_header") is not None else fits.Header()
        invvar_header["OBJECT"] = frb
        invvar_header["SURVEY"] = source
        invvar_header["BAND"] = BAND
        invvar_header["CUTTYPE"] = "10arcmin"
        fits.writeto(invvar_path, product["invvar"], invvar_header, overwrite=True)
        invvar_saved = True

    return flux_path, invvar_saved


def main():
    if not os.path.exists(MASTER_CSV):
        raise FileNotFoundError(f"Missing required file: {MASTER_CSV}")

    os.makedirs(RUNS_DIR, exist_ok=True)
    logs_dir = os.path.join(BASE_DIR, "logs")
    os.makedirs(logs_dir, exist_ok=True)

    df = pd.read_csv(MASTER_CSV)
    records = []

    print("=" * 72)
    print("Bootstrap PSFEx + SExtractor structure and start 10-arcmin cutout fetch")
    print(f"Total FRBs: {len(df)}")
    print(f"Output root: {BASE_DIR}")
    print(f"Cutout size: {CUTOUT_SIZE_PX} px (~10 arcmin at {PIXSCALE} arcsec/px)")
    print("=" * 72)

    start = time.time()
    for idx, row in df.iterrows():
        frb = str(row["FRB"])
        ra = float(row["RA_deg"])
        dec = float(row["DEC_deg"])

        run_dir = ensure_structure(frb)
        print(f"[{idx + 1}/{len(df)}] {frb}  RA={ra:.6f} DEC={dec:.6f}")

        status = "failed"
        source = "none"
        flux_file = ""
        invvar_saved = False

        legacy = fetch_legacy(ra, dec)
        if legacy is not None:
            source = f"Legacy Survey ({legacy['layer']})"
            flux_file, invvar_saved = write_cutouts(frb, run_dir, legacy, source)
            status = "ok"
            print(f"    Saved Legacy cutout -> {flux_file}")
        else:
            ps1 = fetch_ps1_flux(ra, dec)
            if ps1 is not None:
                source = "Pan-STARRS DR1"
                flux_file, invvar_saved = write_cutouts(frb, run_dir, ps1, source)
                status = "ok"
                print(f"    Saved PS1 cutout -> {flux_file}")
            else:
                print("    No cutout available from Legacy or PS1")

        records.append(
            {
                "FRB": frb,
                "RA_deg": ra,
                "DEC_deg": dec,
                "status": status,
                "source": source,
                "flux_file": flux_file,
                "invvar_saved": invvar_saved,
                "run_dir": run_dir,
            }
        )

    manifest = os.path.join(BASE_DIR, "download_manifest_10arcmin.csv")
    pd.DataFrame.from_records(records).to_csv(manifest, index=False)

    elapsed = time.time() - start
    ok_count = sum(1 for r in records if r["status"] == "ok")
    print("=" * 72)
    print(f"Finished in {elapsed:.1f}s")
    print(f"Success: {ok_count}/{len(records)}")
    print(f"Manifest: {manifest}")
    print("=" * 72)


if __name__ == "__main__":
    main()
