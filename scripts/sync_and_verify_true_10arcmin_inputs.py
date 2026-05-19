import os
from io import BytesIO

import pandas as pd
import requests
from astropy.io import fits

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BASE = os.path.join(ROOT, "psfs", "PSFEx + SExtractor")
MANIFEST = os.path.join(BASE, "download_manifest_10arcmin.csv")

LS_PIXSCALE = 0.262
PS1_PIXSCALE = 0.25
LS_SIZE = int(round(600.0 / LS_PIXSCALE))   # 2290 px
PS1_SIZE = int(round(600.0 / PS1_PIXSCALE)) # 2400 px


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


def fetch_ps1_cutout(ra: float, dec: float, size: int, filename: str) -> tuple[object, object]:
    url = (
        "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
        f"ra={ra}&dec={dec}&size={size}&format=fits&red={filename}"
    )
    try:
        resp = requests.get(url, timeout=240)
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


def fetch_ps1_weight(ra: float, dec: float, size: int, flux_filename: str) -> tuple[object, object]:
    wt_filename = flux_filename.replace(".fits", ".wt.fits") if flux_filename.endswith(".fits") else flux_filename + ".wt.fits"
    return fetch_ps1_cutout(ra, dec, size, wt_filename)


def fetch_ls_pair(ra: float, dec: float, size: int) -> tuple[object, object, object, object, str]:
    for layer in ["ls-dr10", "ls-dr9"]:
        url = (
            "https://www.legacysurvey.org/viewer/cutout.fits?"
            f"ra={ra}&dec={dec}&size={size}&layer={layer}&pixscale={LS_PIXSCALE}&bands=r&invvar=True"
        )
        try:
            resp = requests.get(url, timeout=240)
            if resp.status_code != 200:
                continue
            with fits.open(BytesIO(resp.content)) as hdul:
                if len(hdul) < 2:
                    continue
                flux = hdul[0].data
                invv = hdul[1].data
                if flux is None or invv is None:
                    continue
                return flux, hdul[0].header.copy(), invv, hdul[1].header.copy(), layer
        except Exception:
            continue
    return None, None, None, None, ""


def shape2d(path: str):
    with fits.open(path) as hdul:
        return hdul[0].data.shape


def write_cutouts_and_sync(row: pd.Series) -> dict:
    frb = str(row["FRB"])
    ra = float(row["RA_deg"])
    dec = float(row["DEC_deg"])
    source = str(row["source"])
    run_dir = str(row["run_dir"])

    cutouts_dir = os.path.join(run_dir, "cutouts")
    in_dir = os.path.join(run_dir, "input", "sextractor")
    os.makedirs(cutouts_dir, exist_ok=True)
    os.makedirs(in_dir, exist_ok=True)

    flux_path = os.path.join(cutouts_dir, f"{frb}_10arcmin_flux.fits")
    invv_path = os.path.join(cutouts_dir, f"{frb}_10arcmin_invvar.fits")

    status = "ok"
    notes = []

    if "Pan-STARRS" in source:
        if dec < -30:
            return {"FRB": frb, "status": "fail", "notes": "ps1_no_coverage_dec", "source": source}

        fn = get_ps1_filename(ra, dec)
        if not fn:
            return {"FRB": frb, "status": "fail", "notes": "ps1_filename_missing", "source": source}

        redownload_flux = True
        if os.path.exists(flux_path):
            try:
                redownload_flux = shape2d(flux_path) != (PS1_SIZE, PS1_SIZE)
            except Exception:
                redownload_flux = True

        if redownload_flux:
            flux, hdr = fetch_ps1_cutout(ra, dec, PS1_SIZE, fn)
            if flux is None:
                return {"FRB": frb, "status": "fail", "notes": "ps1_flux_fetch_failed", "source": source}
            hdr["OBJECT"] = frb
            hdr["SURVEY"] = "Pan-STARRS DR1"
            hdr["BAND"] = "r"
            hdr["CUTTYPE"] = "10arcmin"
            fits.writeto(flux_path, flux, hdr, overwrite=True)
            notes.append("flux_redownloaded_ps1")

        redownload_invv = True
        if os.path.exists(invv_path):
            try:
                redownload_invv = shape2d(invv_path) != (PS1_SIZE, PS1_SIZE)
            except Exception:
                redownload_invv = True

        if redownload_invv:
            wt, hdrw = fetch_ps1_weight(ra, dec, PS1_SIZE, fn)
            if wt is None:
                return {"FRB": frb, "status": "fail", "notes": "ps1_invvar_fetch_failed", "source": source}
            hdrw["OBJECT"] = frb
            hdrw["SURVEY"] = "Pan-STARRS DR1"
            hdrw["BAND"] = "r"
            hdrw["CUTTYPE"] = "10arcmin"
            hdrw["HISTORY"] = "Pan-STARRS weight map used as inverse-variance map"
            fits.writeto(invv_path, wt, hdrw, overwrite=True)
            notes.append("invvar_redownloaded_ps1")

        expected_size = PS1_SIZE
        pixscale = PS1_PIXSCALE

    else:
        # Legacy Survey target
        redownload = True
        if os.path.exists(flux_path) and os.path.exists(invv_path):
            try:
                redownload = (shape2d(flux_path) != (LS_SIZE, LS_SIZE)) or (shape2d(invv_path) != (LS_SIZE, LS_SIZE))
            except Exception:
                redownload = True

        if redownload:
            flux, hdrf, invv, hdri, layer = fetch_ls_pair(ra, dec, LS_SIZE)
            if flux is None:
                return {"FRB": frb, "status": "fail", "notes": "ls_fetch_failed", "source": source}
            hdrf["OBJECT"] = frb
            hdrf["SURVEY"] = f"Legacy Survey ({layer})"
            hdrf["BAND"] = "r"
            hdrf["CUTTYPE"] = "10arcmin"
            fits.writeto(flux_path, flux, hdrf, overwrite=True)

            hdri["OBJECT"] = frb
            hdri["SURVEY"] = f"Legacy Survey ({layer})"
            hdri["BAND"] = "r"
            hdri["CUTTYPE"] = "10arcmin"
            fits.writeto(invv_path, invv, hdri, overwrite=True)
            notes.append("ls_pair_redownloaded")

        expected_size = LS_SIZE
        pixscale = LS_PIXSCALE

    # Sync run input files from canonical cutouts.
    input_flux = os.path.join(in_dir, "image.fits")
    input_invv = os.path.join(in_dir, "invvar.fits")

    with fits.open(flux_path) as hdul:
        fits.writeto(input_flux, hdul[0].data, hdul[0].header, overwrite=True)
    with fits.open(invv_path) as hdul:
        fits.writeto(input_invv, hdul[0].data, hdul[0].header, overwrite=True)

    # Verify dimensions and FOV consistency.
    fshape = shape2d(flux_path)
    ishape = shape2d(invv_path)
    infshape = shape2d(input_flux)
    iishape = shape2d(input_invv)

    fov_arcmin = (fshape[1] * pixscale) / 60.0
    ok = (fshape == ishape == infshape == iishape == (expected_size, expected_size))

    return {
        "FRB": frb,
        "status": "ok" if ok else "fail",
        "notes": ";".join(notes) if notes else "existing_ok",
        "source": source,
        "expected_size_px": expected_size,
        "flux_shape": f"{fshape[0]}x{fshape[1]}",
        "invvar_shape": f"{ishape[0]}x{ishape[1]}",
        "input_image_shape": f"{infshape[0]}x{infshape[1]}",
        "input_invvar_shape": f"{iishape[0]}x{iishape[1]}",
        "effective_fov_arcmin": round(fov_arcmin, 4),
    }


def main() -> None:
    if not os.path.exists(MANIFEST):
        raise FileNotFoundError(f"Missing manifest: {MANIFEST}")

    df = pd.read_csv(MANIFEST)
    out = []

    for _, row in df.iterrows():
        out.append(write_cutouts_and_sync(row))

    rep = pd.DataFrame(out)
    report_path = os.path.join(BASE, "verification_10arcmin_inputs.csv")
    rep.to_csv(report_path, index=False)

    total = len(rep)
    ok = int((rep["status"] == "ok").sum())
    fail = total - ok

    print(f"TOTAL={total}")
    print(f"OK={ok}")
    print(f"FAIL={fail}")
    print(f"REPORT={report_path}")
    if fail:
        bad = rep.loc[rep["status"] != "ok", ["FRB", "notes"]]
        print("FAIL_ROWS=" + ";".join(f"{r.FRB}:{r.notes}" for _, r in bad.iterrows()))


if __name__ == "__main__":
    main()
