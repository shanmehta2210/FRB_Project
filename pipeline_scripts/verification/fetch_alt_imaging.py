"""Inspection stamps from deeper surveys (and g-band for Fourier spirals).

Uses CDS hips2fits. Default stamps are ~1' FITS + PNG under ``alt_imaging/``.
Those are inspection only. DES DR2 is **not** independent of LS DR10 (DECam
is already in the LS coadd).

HiPS + the pipeline 0.262"/px grid **cannot** recover native HSC (~0.75") or
Euclid VIS (~0.16") FWHM. hips2fits serves an already-resampled pyramid; then
``cutout_resample.write_standardized`` bilinear-zooms onto 2290 px / 10'.
Native-PSF science cutouts need HSC DAS / Euclid archive pixels and must
**not** be resampled to 0.262". ``--install-large-cutout`` is that HiPS path
(depth only). Production ``Output/`` is never touched by the inspection
cohorts.

Cohorts
-------
* mag 20–22 — GALFIT ``mag`` in [20, 22]. Probes **every** survey in
  ``SCOUR_SURVEYS`` (not first-hit). FITS hits get 1' stamps.
* mag 21–22 — GALFIT ``mag`` in (21, 22]. First-hit ladder (legacy).
* spirals — Fourier winding: reliable, |ψ₂'| ≥ 30 °/Re, |slope σ| ≥ 3.

Examples
--------
    python fetch_alt_imaging.py --cohort mag20_22
    python fetch_alt_imaging.py --cohort mag21_22
    python fetch_alt_imaging.py --cohort spirals --force
    python fetch_alt_imaging.py --install-large-cutout 20240208A
"""
from __future__ import annotations

import argparse
import sys
import time
from io import BytesIO
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from astropy.io import fits

import vercommon as vc

HIPS2FITS = "https://alasky.cds.unistra.fr/hips-image-services/hips2fits"
HIPS2FITS_ALT = "https://alaskybis.cds.unistra.fr/hips-image-services/hips2fits"

# 1 arcmin on the sky.
FOV_DEG = 1.0 / 60.0
STAMP_PX = 256
PROBE_PX = 32
PROBE_TIMEOUT_S = 12
STAMP_TIMEOUT_S = 60
RETRIES = 2

# First hit in each list that has real coverage is downloaded.
# UNIONS/CFIS IDs are probed even though CDS currently has none (logged miss).
DEEP_R_LADDER: list[tuple[str, str]] = [
    ("hsc_deep_r", "CDS/P/HSC/DR2/deep/r"),
    ("hsc_wide_r", "CDS/P/HSC/DR2/wide/r"),
    ("hsc_wide_i", "CDS/P/HSC/DR2/wide/i"),
    ("unions_r", "CDS/P/CFIS/r"),
    ("des_r", "CDS/P/DES-DR2/r"),
    ("ls_r", "CDS/P/DESI-Legacy-Surveys/DR10/r"),
]

# Mag [20, 22] scour: probe every layer. kind=
#   fits        — hips2fits FITS; write 1' stamp on ok
#   fits_status — hips2fits FITS; status only (LS baseline, not an upgrade)
#   color       — jpeg/png coverage only (no science FITS HiPS)
#   das         — unauthenticated HSC PDR3 DAS (expect 401)
#   unions      — CFIS HiPS (usually 400) + optional CADC preview GET
SCOUR_SURVEYS: list[tuple[str, str, str]] = [
    ("hsc_deep_r", "CDS/P/HSC/DR2/deep/r", "fits"),
    ("hsc_wide_r", "CDS/P/HSC/DR2/wide/r", "fits"),
    ("hsc_wide_i", "CDS/P/HSC/DR2/wide/i", "fits"),
    ("hscla_r", "CDS/P/HSCLA/2016/r", "fits"),
    ("cfhtls_d_r", "CDS/P/CFHTLS/D/r", "fits"),
    ("cfhtls_w_r", "CDS/P/CFHTLS/W/r", "fits"),
    ("euclid_q1_vis", "CDS/P/Euclid/Q1/VIS", "fits"),
    ("euclid_q2_vis", "CDS/P/Euclid/Q2/VIS", "fits"),
    ("hst_sdsr", "CDS/P/HST/SDSSr", "fits"),
    ("hst_i", "CDS/P/HST/I", "fits"),
    ("hst_wideV", "CDS/P/HST/wideV", "fits"),
    ("jwst_f150w", "CDS/P/JWST/F150W", "fits"),
    ("jwst_nircam", "ESAVO/P/JWST/NIRCam_Imaging", "fits"),
    ("cosmosweb_acs", "CDS/P/COSMOS-Web/DR1/HST/ACS/F814W", "fits"),
    ("ps1_r", "CDS/P/PanSTARRS/DR1/r", "fits"),
    ("ls_r", "CDS/P/DESI-Legacy-Surveys/DR10/r", "fits_status"),
    ("kids_gri", "CDS/P/KiDS/DR5/color-gri", "color"),
    ("unions_r", "CDS/P/CFIS/r", "unions"),
    ("hsc_das_pdr3", "", "das"),
]

HSC_DAS_CUTOUT = (
    "https://hsc-release.mtk.nao.ac.jp/das_cutout/pdr3/cgi-bin/cutout"
)
CADC_UNIONS_PREVIEW = (
    "https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/community/unions/preview.html"
)
DEEPER_IDS = frozenset({
    "hsc_deep_r", "hsc_wide_r", "hsc_wide_i", "hscla_r",
    "cfhtls_d_r", "cfhtls_w_r",
    "euclid_q1_vis", "euclid_q2_vis",
    "hst_sdsr", "hst_i", "hst_wideV",
    "jwst_f150w", "jwst_nircam", "cosmosweb_acs",
    "kids_gri",
})

G_LADDER: list[tuple[str, str]] = [
    ("hsc_deep_g", "CDS/P/HSC/DR2/deep/g"),
    ("hsc_wide_g", "CDS/P/HSC/DR2/wide/g"),
    ("unions_g", "CDS/P/CFIS/g"),
    ("des_g", "CDS/P/DES-DR2/g"),
    ("ls_g", "CDS/P/DESI-Legacy-Surveys/DR10/g"),
]

MAG_LO = 21.0
MAG_HI = 22.0
MAG20_LO = 20.0
MAG20_HI = 22.0
SPIRAL_SLOPE_MIN = 30.0
SPIRAL_SIG_MIN = 3.0

ALT_ROOT = Path(vc.VER_DIR) / "alt_imaging"
LOC_CSV = Path(vc.REPO) / "master_frb_localization.csv"
METRICS_CSV = Path(vc.TABLES_ROOT) / "fit_verification_metrics.csv"


def _flux_ok(arr: np.ndarray | None) -> bool:
    if arr is None:
        return False
    a = np.asarray(arr, dtype=float)
    finite = np.isfinite(a)
    if int(np.count_nonzero(finite)) < max(8, int(0.05 * a.size)):
        return False
    return bool(np.nanmax(np.abs(a[finite])) > 0)


def _hips2fits(ra: float, dec: float, hips: str, *, width: int, height: int,
               fov: float, timeout: float, fmt: str = "fits"
               ) -> tuple[np.ndarray | None, object | None, str]:
    params = {
        "hips": hips,
        "ra": f"{ra:.8f}",
        "dec": f"{dec:.8f}",
        "fov": f"{fov:.8f}",
        "width": str(int(width)),
        "height": str(int(height)),
        "projection": "TAN",
        "coordsys": "icrs",
        "format": fmt,
    }
    last = "no_response"
    for base in (HIPS2FITS, HIPS2FITS_ALT):
        for attempt in range(1, RETRIES + 1):
            try:
                resp = requests.get(base, params=params, timeout=timeout)
                if resp.status_code in (400, 404, 422):
                    return None, None, f"http_{resp.status_code}"
                if resp.status_code != 200:
                    last = f"http_{resp.status_code}"
                    time.sleep(1 * attempt)
                    continue
                if fmt != "fits":
                    if len(resp.content) <= 200:
                        return None, None, "empty"
                    try:
                        from PIL import Image
                        arr = np.asarray(Image.open(BytesIO(resp.content)),
                                         dtype=float)
                    except (OSError, ValueError):
                        return None, None, "empty"
                    if (not np.isfinite(arr).any()) or float(np.nanstd(arr)) < 1.0:
                        return None, None, "empty"
                    return None, None, "ok"
                with fits.open(BytesIO(resp.content)) as hdul:
                    data = np.squeeze(hdul[0].data)
                    hdr = hdul[0].header.copy()
                if not _flux_ok(data):
                    return None, None, "empty"
                return np.asarray(data, dtype=float), hdr, "ok"
            except (requests.Timeout, requests.ConnectTimeout):
                last = "timeout"
                continue
            except (requests.RequestException, OSError, ValueError) as exc:
                last = type(exc).__name__
                time.sleep(1 * attempt)
        # only fall through to the mirror on connection/timeout failures
        if last not in ("timeout", "no_response") and not last.startswith("HTTP"):
            break
    return None, None, last


def probe(ra: float, dec: float, hips: str, *, fmt: str = "fits") -> str:
    _, _, status = _hips2fits(
        ra, dec, hips, width=PROBE_PX, height=PROBE_PX, fov=FOV_DEG,
        timeout=PROBE_TIMEOUT_S, fmt=fmt,
    )
    return status


def probe_hsc_das(ra: float, dec: float) -> str:
    """Unauthenticated PDR3 DAS cutout. Native PSF needs a STARS account."""
    params = {
        "ra": f"{ra:.8f}",
        "dec": f"{dec:.8f}",
        "sw": "30asec",
        "sh": "30asec",
        "type": "coadd",
        "image": "on",
        "filter": "HSC-R2",
        "rerun": "pdr3_wide",
    }
    try:
        resp = requests.get(HSC_DAS_CUTOUT, params=params, timeout=PROBE_TIMEOUT_S)
        return f"http_{resp.status_code}"
    except (requests.Timeout, requests.ConnectTimeout):
        return "timeout"
    except requests.RequestException as exc:
        return type(exc).__name__


def probe_unions(ra: float, dec: float) -> str:
    hips_status = probe(ra, dec, "CDS/P/CFIS/r")
    if hips_status == "ok":
        return "ok"
    try:
        resp = requests.get(
            CADC_UNIONS_PREVIEW,
            params={"ra": f"{ra:.8f}", "dec": f"{dec:.8f}"},
            timeout=PROBE_TIMEOUT_S,
            allow_redirects=True,
        )
        if resp.status_code in (401, 403):
            return f"cadc_{resp.status_code}"
    except (requests.Timeout, requests.ConnectTimeout, requests.RequestException):
        pass
    return hips_status


def fetch_stamp(ra: float, dec: float, hips: str
                ) -> tuple[np.ndarray | None, object | None, str]:
    return _hips2fits(
        ra, dec, hips, width=STAMP_PX, height=STAMP_PX, fov=FOV_DEG,
        timeout=STAMP_TIMEOUT_S,
    )


def _asinh_png(data: np.ndarray, dest: Path, title: str) -> None:
    finite = np.isfinite(data)
    if not finite.any():
        return
    lo, hi = np.nanpercentile(data[finite], [1, 99])
    x = np.clip((data - lo) / (hi - lo + 1e-30), 0.0, 1.0)
    img = np.arcsinh(10.0 * x) / np.arcsinh(10.0)
    fig, ax = plt.subplots(figsize=(4.2, 4.2), constrained_layout=True)
    ax.imshow(img, origin="lower", cmap="gray", vmin=0, vmax=1,
              interpolation="nearest")
    ax.set_title(title, fontsize=9)
    ax.set_xticks([])
    ax.set_yticks([])
    fig.savefig(dest, dpi=140)
    plt.close(fig)


def _write_stamp(dest_fits: Path, dest_png: Path, data: np.ndarray, hdr,
                 title: str) -> None:
    dest_fits.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(data=np.asarray(data, dtype=np.float32),
                    header=hdr).writeto(dest_fits, overwrite=True)
    _asinh_png(data, dest_png, title)


def _localizations() -> pd.DataFrame:
    loc = pd.read_csv(LOC_CSV, dtype={"frb": str})
    return loc[["frb", "ra_deg", "dec_deg"]].drop_duplicates("frb")


def mag21_22_hosts() -> pd.DataFrame:
    df = vc.cohort("all64").merge(_localizations(), on="frb", how="left")
    mag = pd.to_numeric(df["mag"], errors="coerce")
    return df.loc[(mag > MAG_LO) & (mag <= MAG_HI)].reset_index(drop=True)


def mag20_22_hosts() -> pd.DataFrame:
    df = vc.cohort("all64").merge(_localizations(), on="frb", how="left")
    mag = pd.to_numeric(df["mag"], errors="coerce")
    return df.loc[(mag >= MAG20_LO) & (mag <= MAG20_HI)].reset_index(drop=True)


def spiral_hosts() -> pd.DataFrame:
    mets = pd.read_csv(METRICS_CSV, dtype={"frb": str})
    rel = mets["fourier_reliable"].astype(str).str.lower().eq("true")
    slope = pd.to_numeric(mets["fourier_m2_phase_slope_deg_per_re"],
                          errors="coerce")
    sig = pd.to_numeric(mets["fourier_m2_phase_slope_sig"], errors="coerce")
    flagged = rel & slope.abs().ge(SPIRAL_SLOPE_MIN) & sig.abs().ge(SPIRAL_SIG_MIN)
    cols = ["frb", "fourier_reliable", "fourier_m2_phase_slope_deg_per_re",
            "fourier_m2_phase_slope_sig", "mag"]
    keep = [c for c in cols if c in mets.columns]
    out = mets.loc[flagged, keep].copy()
    return out.merge(_localizations(), on="frb", how="left")


def _probe_all(ra: float, dec: float, ladder: list[tuple[str, str]]) -> dict:
    hits: dict[str, str] = {}
    for name, hips in ladder:
        hits[name] = probe(ra, dec, hips)
        print(f"    probe {name}: {hits[name]}", flush=True)
    return hits


def _first_ok(hits: dict[str, str], ladder: list[tuple[str, str]]
              ) -> tuple[str, str] | None:
    for name, hips in ladder:
        if hits.get(name) == "ok":
            return name, hips
    return None


def process_cohort(df: pd.DataFrame, outdir: Path, ladder: list[tuple[str, str]],
                   *, band: str, force: bool) -> pd.DataFrame:
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for i, row in df.iterrows():
        frb = str(row["frb"])
        ra, dec = float(row["ra_deg"]), float(row["dec_deg"])
        print(f"[{i + 1}/{len(df)}] {frb}  ra={ra:.4f} dec={dec:.4f}", flush=True)
        rec: dict = {
            "frb": frb,
            "ra_deg": ra,
            "dec_deg": dec,
            "mag": row.get("mag", float("nan")),
            "band": band,
        }
        if not (np.isfinite(ra) and np.isfinite(dec)):
            rec["fetch_status"] = "no_coords"
            rec["chosen_survey"] = ""
            rec["chosen_hips"] = ""
            rec["stamp_fits"] = ""
            rec["stamp_png"] = ""
            rows.append(rec)
            print("    skip: missing RA/Dec")
            continue
        if "fourier_m2_phase_slope_deg_per_re" in row.index:
            rec["psi2_prime"] = row["fourier_m2_phase_slope_deg_per_re"]
            rec["psi2_prime_sig"] = row.get("fourier_m2_phase_slope_sig")
        hits = _probe_all(ra, dec, ladder)
        for name, _hips in ladder:
            rec[f"{name}_status"] = hits.get(name, "")
        chosen = _first_ok(hits, ladder)
        rec["chosen_survey"] = chosen[0] if chosen else ""
        rec["chosen_hips"] = chosen[1] if chosen else ""
        rec["stamp_fits"] = ""
        rec["stamp_png"] = ""
        rec["fetch_status"] = "no_coverage"
        if chosen is None:
            rows.append(rec)
            continue
        name, hips = chosen
        stem = outdir / f"{frb}_{name}"
        fits_path, png_path = Path(f"{stem}.fits"), Path(f"{stem}.png")
        if fits_path.is_file() and png_path.is_file() and not force:
            rec["stamp_fits"] = str(fits_path.relative_to(Path(vc.VER_DIR)))
            rec["stamp_png"] = str(png_path.relative_to(Path(vc.VER_DIR)))
            rec["fetch_status"] = "cached"
            rows.append(rec)
            print(f"    cached {fits_path.name}", flush=True)
            continue
        data, hdr, status = fetch_stamp(ra, dec, hips)
        rec["fetch_status"] = status
        if data is None:
            rows.append(rec)
            continue
        title = f"{frb}  {name}  {band}"
        _write_stamp(fits_path, png_path, data, hdr, title)
        rec["stamp_fits"] = str(fits_path.relative_to(Path(vc.VER_DIR)))
        rec["stamp_png"] = str(png_path.relative_to(Path(vc.VER_DIR)))
        print(f"    wrote {fits_path.name}", flush=True)
        rows.append(rec)
    return pd.DataFrame(rows)


def _probe_survey(ra: float, dec: float, name: str, hips: str, kind: str) -> str:
    if kind == "das":
        return probe_hsc_das(ra, dec)
    if kind == "unions":
        return probe_unions(ra, dec)
    if kind == "color":
        status = probe(ra, dec, hips, fmt="jpg")
        if status.startswith("http_"):
            status = probe(ra, dec, hips, fmt="png")
        return status
    return probe(ra, dec, hips, fmt="fits")


def process_scour(df: pd.DataFrame, outdir: Path, *, force: bool) -> pd.DataFrame:
    """Probe every SCOUR_SURVEYS layer; stamp all FITS hits except LS baseline."""
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for i, row in df.iterrows():
        frb = str(row["frb"])
        ra, dec = float(row["ra_deg"]), float(row["dec_deg"])
        print(f"[{i + 1}/{len(df)}] {frb}  ra={ra:.4f} dec={dec:.4f}  mag={row.get('mag', '')}",
              flush=True)
        rec: dict = {
            "frb": frb,
            "ra_deg": ra,
            "dec_deg": dec,
            "mag": row.get("mag", float("nan")),
        }
        stamp_fits: list[str] = []
        stamp_png: list[str] = []
        hits: list[str] = []
        if not (np.isfinite(ra) and np.isfinite(dec)):
            for name, _hips, _kind in SCOUR_SURVEYS:
                rec[f"{name}_status"] = "no_coords"
            rec["hits"] = ""
            rec["deeper_hits"] = ""
            rec["stamp_fits"] = ""
            rec["stamp_png"] = ""
            rows.append(rec)
            print("    skip: missing RA/Dec")
            continue
        for name, hips, kind in SCOUR_SURVEYS:
            status = _probe_survey(ra, dec, name, hips, kind)
            rec[f"{name}_status"] = status
            print(f"    probe {name}: {status}", flush=True)
            if status != "ok":
                continue
            hits.append(name)
            if kind != "fits":
                continue
            stem = outdir / f"{frb}_{name}"
            fits_path, png_path = Path(f"{stem}.fits"), Path(f"{stem}.png")
            if fits_path.is_file() and png_path.is_file() and not force:
                stamp_fits.append(str(fits_path.relative_to(Path(vc.VER_DIR))))
                stamp_png.append(str(png_path.relative_to(Path(vc.VER_DIR))))
                print(f"    cached {fits_path.name}", flush=True)
                continue
            data, hdr, fetch_status = fetch_stamp(ra, dec, hips)
            if data is None:
                rec[f"{name}_status"] = fetch_status
                if name in hits:
                    hits.remove(name)
                print(f"    stamp {name}: {fetch_status}", flush=True)
                continue
            title = f"{frb}  {name}"
            _write_stamp(fits_path, png_path, data, hdr, title)
            stamp_fits.append(str(fits_path.relative_to(Path(vc.VER_DIR))))
            stamp_png.append(str(png_path.relative_to(Path(vc.VER_DIR))))
            print(f"    wrote {fits_path.name}", flush=True)
        rec["hits"] = ",".join(hits)
        rec["deeper_hits"] = ",".join(n for n in hits if n in DEEPER_IDS)
        rec["stamp_fits"] = ";".join(stamp_fits)
        rec["stamp_png"] = ";".join(stamp_png)
        rows.append(rec)
    return pd.DataFrame(rows)


def install_large_cutout(
    frb: str,
    hips: str,
    *,
    ra: float,
    dec: float,
    survey: str,
) -> Path:
    """Write a pipeline-standard 10′ flux+invvar pair into ``large_cutouts/``.

    Backs up any existing pair to ``large_cutouts/_pre_<survey>/``. Invvar is
    the sky-MAD map (hips2fits has no variance plane). Phase 2 still sets the
    photometric ZP from PS1/LS stars, so HiPS linear flux is enough.
    """
    scripts = Path(vc.REPO) / "scripts"
    if str(scripts) not in sys.path:
        sys.path.insert(0, str(scripts))
    from cutout_download import load_registry, save_registry, upsert_row
    from cutout_resample import (  # noqa: WPS433
        FOV_ARCMIN,
        TARGET_PIXSCALE,
        TARGET_SIZE,
        invvar_from_sky,
        write_standardized,
    )

    cut_dir = Path(vc.REPO) / "large_cutouts"
    flux_path = cut_dir / f"{frb}_flux.fits"
    inv_path = cut_dir / f"{frb}_invvar.fits"
    bak = cut_dir / f"_pre_{survey}"
    bak.mkdir(parents=True, exist_ok=True)
    for src in (flux_path, inv_path):
        if src.is_file():
            dest = bak / src.name
            if not dest.is_file():
                dest.write_bytes(src.read_bytes())
                print(f"backed up {src.name} -> {dest}", flush=True)

    fov = FOV_ARCMIN / 60.0
    print(f"hips2fits 10' {hips}  {frb}  {TARGET_SIZE}px ...", flush=True)
    data, hdr, status = _hips2fits(
        ra, dec, hips,
        width=TARGET_SIZE, height=TARGET_SIZE, fov=fov,
        timeout=max(STAMP_TIMEOUT_S * 5, 300),
    )
    if data is None:
        raise SystemExit(f"10′ fetch failed ({status}) for {hips}")
    flux = np.asarray(data, dtype=np.float64)
    good = np.isfinite(flux) & (np.abs(flux) > 0)
    if good.mean() < 0.2:
        good = np.isfinite(flux)
    inv = invvar_from_sky(flux, good)
    if hdr.get("CDELT2") is not None:
        native = abs(float(hdr["CDELT2"])) * 3600.0
    elif hdr.get("CD2_2") is not None:
        native = abs(float(hdr["CD2_2"])) * 3600.0
    else:
        native = TARGET_PIXSCALE
    write_standardized(
        str(flux_path), str(inv_path),
        flux, hdr, inv, hdr, ra, dec, native,
        center_on_array=False,
    )
    for path in (flux_path, inv_path):
        with fits.open(path, mode="update") as hdul:
            hdul[0].header["SURVEY"] = survey.upper()
            hdul[0].header["HIPS"] = hips
            hdul.flush()
    df = upsert_row(
        load_registry(), frb, ra, dec,
        source=survey, layer=hips, resampled=True, status="ok",
    )
    save_registry(df)
    print(f"wrote {flux_path} and {inv_path}", flush=True)
    return flux_path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cohort", choices=["all", "mag20_22", "mag21_22", "spirals"],
                    default="all")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--install-large-cutout", metavar="FRB", default=None,
                    help="replace large_cutouts/<FRB>_{flux,invvar}.fits with a 10′ HiPS stamp")
    ap.add_argument("--hips", default="CDS/P/HSC/DR2/wide/r",
                    help="HiPS id for --install-large-cutout")
    args = ap.parse_args(argv)

    if args.install_large_cutout:
        frb = str(args.install_large_cutout)
        loc = _localizations()
        row = loc.loc[loc["frb"] == frb]
        if row.empty:
            raise SystemExit(f"no RA/Dec for {frb} in master_frb_localization.csv")
        ra, dec = float(row.iloc[0]["ra_deg"]), float(row.iloc[0]["dec_deg"])
        hips_u = args.hips.upper()
        if "HSC" in hips_u:
            survey = "hsc"
        elif "DES" in hips_u and "LEGACY" not in hips_u:
            survey = "des"
        else:
            survey = "ls"
        install_large_cutout(frb, args.hips, ra=ra, dec=dec, survey=survey)
        return 0

    do_scour = args.cohort == "mag20_22"
    do_mag = args.cohort in ("all", "mag21_22")
    do_sp = args.cohort in ("all", "spirals")
    rc = 0

    if do_scour:
        hosts = mag20_22_hosts()
        print(f"mag 20-22 scour: {len(hosts)} hosts x {len(SCOUR_SURVEYS)} surveys",
              flush=True)
        cov = process_scour(hosts, ALT_ROOT / "mag20_22", force=args.force)
        dest = ALT_ROOT / "mag20_22" / "coverage.csv"
        dest.parent.mkdir(parents=True, exist_ok=True)
        cov.to_csv(dest, index=False)
        print(f"wrote {dest}")
        return 0

    if do_mag:
        hosts = mag21_22_hosts()
        print(f"mag 21–22: {len(hosts)} hosts", flush=True)
        cov = process_cohort(hosts, ALT_ROOT / "mag21_22", DEEP_R_LADDER,
                             band="r/i", force=args.force)
        dest = ALT_ROOT / "mag21_22" / "coverage.csv"
        dest.parent.mkdir(parents=True, exist_ok=True)
        cov.to_csv(dest, index=False)
        print(f"wrote {dest}")

    if do_sp:
        hosts = spiral_hosts()
        print(f"spirals: {len(hosts)} hosts", flush=True)
        if hosts.empty:
            print("  (no Fourier-winding hosts at the |ψ2'|≥30, |σ|≥3 gate)")
        cov = process_cohort(hosts, ALT_ROOT / "spirals", G_LADDER,
                             band="g", force=args.force)
        dest = ALT_ROOT / "spirals" / "coverage.csv"
        dest.parent.mkdir(parents=True, exist_ok=True)
        cov.to_csv(dest, index=False)
        print(f"wrote {dest}")

    return rc


if __name__ == "__main__":
    raise SystemExit(main())
