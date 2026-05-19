import os
from io import BytesIO

import pandas as pd
import requests
from astropy.io import fits

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BASE = os.path.join(ROOT, "psfs", "PSFEx + SExtractor")
MANIFEST = os.path.join(BASE, "download_manifest_10arcmin.csv")


def get_ps1_filename(ra: float, dec: float) -> str | None:
    url = (
        "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py?"
        f"ra={ra}&dec={dec}&filters=r&sep=comma"
    )
    try:
        resp = requests.get(url, timeout=45)
        if resp.status_code != 200:
            return None
        lines = resp.text.strip().split("\n")
        if len(lines) < 2:
            return None
        keys = lines[0].split(",")
        vals = lines[1].split(",")
        row = dict(zip(keys, vals))
        return row.get("filename")
    except Exception:
        return None


def fetch_ps1_weight_cutout(ra: float, dec: float, size: int):
    filename = get_ps1_filename(ra, dec)
    if not filename:
        return None, None

    wt_filename = filename.replace(".fits", ".wt.fits") if filename.endswith(".fits") else filename + ".wt.fits"

    url = (
        "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
        f"ra={ra}&dec={dec}&size={size}&format=fits&red={wt_filename}"
    )
    try:
        resp = requests.get(url, timeout=180)
        if resp.status_code != 200:
            return None, None
        if b"<html" in resp.content[:200].lower():
            return None, None
        with fits.open(BytesIO(resp.content)) as hdul:
            if len(hdul) < 1 or hdul[0].data is None:
                return None, None
            return hdul[0].data, hdul[0].header.copy()
    except Exception:
        return None, None


def ensure_cutout_invvar(row: pd.Series) -> tuple[bool, str]:
    frb = row["FRB"]
    ra = float(row["RA_deg"])
    dec = float(row["DEC_deg"])
    source = str(row["source"])
    run_dir = row["run_dir"]

    cutouts_dir = os.path.join(run_dir, "cutouts")
    flux_path = os.path.join(cutouts_dir, f"{frb}_10arcmin_flux.fits")
    inv_path = os.path.join(cutouts_dir, f"{frb}_10arcmin_invvar.fits")

    if os.path.exists(inv_path):
        return True, "present"
    if not os.path.exists(flux_path):
        return False, "missing_flux"

    # Use flux dimensions to enforce same 10-arcmin pixel grid in fetched weight map.
    with fits.open(flux_path) as hdul:
        flux_header = hdul[0].header.copy()
        ny, nx = hdul[0].data.shape
        size = int(nx)

    if "Pan-STARRS" in source:
        data, hdr = fetch_ps1_weight_cutout(ra, dec, size=size)
        if data is None:
            return False, "ps1_fetch_failed"

        hdr["OBJECT"] = frb
        hdr["SURVEY"] = "Pan-STARRS DR1"
        hdr["BAND"] = "r"
        hdr["CUTTYPE"] = "10arcmin"
        hdr["HISTORY"] = "Pan-STARRS weight map used as inverse-variance map"
        fits.writeto(inv_path, data, hdr, overwrite=True)
        return True, "downloaded_ps1"

    if "Legacy Survey" in source:
        # Safety path if any LS invvar is missing: refetch with invvar=True.
        layers = ["ls-dr10", "ls-dr9"]
        for layer in layers:
            url = (
                "https://www.legacysurvey.org/viewer/cutout.fits?"
                f"ra={ra}&dec={dec}&size={size}&layer={layer}&pixscale=0.262&bands=r&invvar=True"
            )
            try:
                resp = requests.get(url, timeout=180)
                if resp.status_code != 200:
                    continue
                with fits.open(BytesIO(resp.content)) as hdul:
                    if len(hdul) < 2 or hdul[1].data is None:
                        continue
                    hdr = hdul[1].header.copy()
                    hdr["OBJECT"] = frb
                    hdr["SURVEY"] = f"Legacy Survey ({layer})"
                    hdr["BAND"] = "r"
                    hdr["CUTTYPE"] = "10arcmin"
                    fits.writeto(inv_path, hdul[1].data, hdr, overwrite=True)
                    return True, "downloaded_ls"
            except Exception:
                continue
        return False, "ls_fetch_failed"

    return False, "unknown_source"


def sync_input_links(row: pd.Series) -> tuple[bool, str]:
    frb = row["FRB"]
    run_dir = row["run_dir"]

    cutouts_dir = os.path.join(run_dir, "cutouts")
    in_dir = os.path.join(run_dir, "input", "sextractor")

    src_flux = os.path.join(cutouts_dir, f"{frb}_10arcmin_flux.fits")
    src_inv = os.path.join(cutouts_dir, f"{frb}_10arcmin_invvar.fits")
    dst_flux = os.path.join(in_dir, "image.fits")
    dst_inv = os.path.join(in_dir, "invvar.fits")

    if not os.path.exists(src_flux) or not os.path.exists(src_inv):
        return False, "missing_cutout_inputs"

    os.makedirs(in_dir, exist_ok=True)

    # Copy (not symlink) for cross-platform simplicity.
    with fits.open(src_flux) as hdul:
        fits.writeto(dst_flux, hdul[0].data, hdul[0].header, overwrite=True)
    with fits.open(src_inv) as hdul:
        fits.writeto(dst_inv, hdul[0].data, hdul[0].header, overwrite=True)

    return True, "synced"


def main() -> None:
    if not os.path.exists(MANIFEST):
        raise FileNotFoundError(f"Missing manifest: {MANIFEST}")

    df = pd.read_csv(MANIFEST)

    fetched_ok = 0
    fetched_fail = []
    sync_ok = 0
    sync_fail = []

    for _, row in df.iterrows():
        frb = row["FRB"]
        ok, msg = ensure_cutout_invvar(row)
        if ok:
            fetched_ok += 1
        else:
            fetched_fail.append(f"{frb}:{msg}")
            continue

        ok2, msg2 = sync_input_links(row)
        if ok2:
            sync_ok += 1
        else:
            sync_fail.append(f"{frb}:{msg2}")

    print(f"CUTOUT_INVVAR_READY={fetched_ok}/{len(df)}")
    print(f"INPUT_SYNC_READY={sync_ok}/{len(df)}")
    if fetched_fail:
        print("FETCH_FAIL_LIST=" + ",".join(fetched_fail))
    if sync_fail:
        print("SYNC_FAIL_LIST=" + ",".join(sync_fail))


if __name__ == "__main__":
    main()
