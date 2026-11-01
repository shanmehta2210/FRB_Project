#!/usr/bin/env python3
"""PS1 image-provenance refits for the 4 Bhardwaj+2024 hosts they fit on PS1.

Re-runs the production chain (master_run.py: SExtractor + PSFEx, photometry,
cutouts, GALFIT with PSF) on PS1 DR1 r-band imaging instead of LS DR10, plus a
GALFIT refit of the same PS1 stamp with no PSF (``D) none``). Everything is
written under Publication/materials/image_provenance/ps1_refits/<FRB>/;
production large_cutouts/, pipeline_scripts/Output/ and seed_56.csv are only
read.

Steps (each idempotent, run in order):
  python scripts/ps1_provenance_refits.py download  [--frb ...]
  python scripts/ps1_provenance_refits.py pipeline  [--frb ...]
  python scripts/ps1_provenance_refits.py variants  [--frb ...]
  python scripts/ps1_provenance_refits.py table
"""
from __future__ import annotations

import argparse
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))
sys.path.insert(0, str(REPO / "pipeline_scripts" / "galfit_fitting"))

from cutout_fetch_common import PS1_PIXSCALE, fetch_ps1_pair  # noqa: E402
from cutout_resample import write_standardized  # noqa: E402

FRBS = ["20190714A", "20210320C", "20220509G", "20220914A"]
BASE = REPO / "Publication" / "materials" / "image_provenance"
OUT = BASE / "ps1_refits"
PROD = REPO / "pipeline_scripts" / "Output"
SEED = REPO / "Publication" / "materials" / "host_properties" / "seed_56.csv"
LOC = REPO / "master_frb_localization.csv"
MASTER = REPO / "pipeline_scripts" / "master_run.py"
WINNING_PSF_LEG = REPO / "pipeline_scripts" / "verification" / "Re-fits"
Q0 = 0.2
GALFIT_INPUTS = ("host_cutout.fits", "host_sigma.fits", "host_mask.fits",
                 "proto_image.fits", "constraints.txt")


def hubble_i(q: float) -> float:
    if q is None or not math.isfinite(q):
        return float("nan")
    if q >= 1.0:
        return 0.0
    if q <= Q0:
        return 90.0
    return math.degrees(math.acos(math.sqrt((q * q - Q0 * Q0) / (1.0 - Q0 * Q0))))


def coords(frb: str) -> tuple[float, float]:
    loc = pd.read_csv(LOC, dtype={"frb": str})
    r = loc.loc[loc.frb == frb].iloc[0]
    return float(r.ra_deg), float(r.dec_deg)


def input_paths(frb: str) -> tuple[Path, Path]:
    d = OUT / frb / "input"
    return d / f"{frb}_flux.fits", d / f"{frb}_invvar.fits"


def run_dir(frb: str) -> Path:
    return OUT / frb / f"{frb}_all"


# ---------------------------------------------------------------- download
def cmd_download(frbs: list[str], force: bool) -> None:
    for frb in frbs:
        flux_p, inv_p = input_paths(frb)
        if flux_p.is_file() and inv_p.is_file() and not force:
            print(f"{frb}: PS1 cutout already present")
            continue
        flux_p.parent.mkdir(parents=True, exist_ok=True)
        ra, dec = coords(frb)
        print(f"{frb}: fetching PS1 r (ra={ra:.5f} dec={dec:.5f}) ...", flush=True)
        flux, fh, inv, ih = fetch_ps1_pair(ra, dec)
        if flux is None or inv is None:
            print(f"{frb}: PS1 download FAILED", flush=True)
            continue
        write_standardized(str(flux_p), str(inv_p), flux, fh, inv, ih, ra, dec,
                           PS1_PIXSCALE, center_on_array=False)
        print(f"{frb}: wrote {flux_p.relative_to(REPO)}", flush=True)


# ---------------------------------------------------------------- pipeline
def production_overrides(frb: str) -> list[str]:
    """Replay the CLI overrides the production run used (from its master_run.log)."""
    log = (PROD / f"{frb}_all" / "master_run.log").read_text(encoding="utf-8", errors="replace")
    m = re.search(r"YAML overrides: (\{.*\})", log)
    ov = eval(m.group(1)) if m else {}  # noqa: S307 - our own log, dict literal
    args: list[str] = []
    for key, flag in (("err_a_arcsec", "--err-a-arcsec"), ("err_b_arcsec", "--err-b-arcsec"),
                      ("err_theta_deg", "--err-theta-deg"), ("p_u", "--p-u")):
        if key in ov:
            args += [flag, str(ov[key])]
    if ov.get("use_localization_host"):
        args.append("--use-localization-host")
    return args


def cmd_pipeline(frbs: list[str]) -> None:
    for frb in frbs:
        flux_p, inv_p = input_paths(frb)
        if not flux_p.is_file():
            print(f"{frb}: no PS1 input, skipping")
            continue
        ra, dec = coords(frb)
        cmd = [sys.executable, str(MASTER), "--image", str(flux_p), "--invvar", str(inv_p),
               "--ra", str(ra), "--dec", str(dec), "--outputs", "all",
               "--frb-name", frb, "--output-dir", str(OUT / frb),
               *production_overrides(frb)]
        log_path = OUT / frb / "master_run.log"
        print(f"{frb}: {' '.join(cmd)}", flush=True)
        with open(log_path, "w", encoding="utf-8") as fh:
            rc = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT,
                                stdin=subprocess.DEVNULL, cwd=str(REPO)).returncode
        print(f"{frb}: master_run exit {rc} (log {log_path.relative_to(REPO)})", flush=True)


# ---------------------------------------------------------------- variants
def _run_galfit(wkdir: Path) -> bool:
    from run_galfit_fitting import run_galfit
    for stale in ("fit.log", "out.fits"):
        (wkdir / stale).unlink(missing_ok=True)
    for g in wkdir.glob("galfit.[0-9][0-9]"):
        g.unlink()
    return bool(run_galfit(str(wkdir)))


def _components(text: str) -> tuple[str, list[str]]:
    """Split a feedme into (header, [component blocks]) ignoring trailing ====."""
    parts = re.split(r"(?=^# Component number:)", text, flags=re.M)
    header, comps = parts[0], parts[1:]
    comps = [re.sub(r"\n=+\s*$", "\n", c.rstrip() + "\n") for c in comps]
    return header, comps


def _assemble(header: str, comps: list[str]) -> str:
    out = header.rstrip() + "\n\n"
    for i, c in enumerate(comps, 1):
        c = re.sub(r"^# Component number:\s*\d+", f"# Component number: {i}", c, flags=re.M)
        out += c.rstrip() + "\n\n"
    return out + "=" * 80 + "\n"


def _set_line(header: str, key: str, value: str) -> str:
    return re.sub(rf"^{key}\).*$", f"{key}) {value}", header, count=1, flags=re.M)


def _star_component_text(frb: str, src: Path, use_psf: bool) -> str | None:
    """Extra point-source component at the star the winning LS 'psf' leg added."""
    leg = WINNING_PSF_LEG / frb / "sandbox"
    g01 = leg / "galfit.01"
    if not g01.is_file():
        return None
    _, comps = _components(g01.read_text(encoding="utf-8", errors="replace"))
    star = next((c for c in comps if re.search(r"^\s*0\)\s*psf", c, re.M)), None)
    if star is None:
        return None
    xs, ys = map(float, re.search(r"^\s*1\)\s+(\S+)\s+(\S+)", star, re.M).groups())
    host_ls = next(c for c in comps if re.search(r"^\s*0\)\s*sersic", c, re.M))
    dmag = (float(re.search(r"^\s*3\)\s+(\S+)", star, re.M).group(1))
            - float(re.search(r"^\s*3\)\s+(\S+)", host_ls, re.M).group(1)))
    mag = host_fit(src)["mag"] + dmag
    # LS sandbox stamp pixel -> sky -> PS1 stamp pixel (both stamp headers carry WCS).
    w_ls = WCS(fits.getheader(leg / "host_cutout.fits")).celestial
    w_ps = WCS(fits.getheader(src / "host_cutout.fits")).celestial
    ra, dec = w_ls.wcs_pix2world(xs, ys, 1)
    xp, yp = w_ps.wcs_world2pix(ra, dec, 1)
    ny, nx = fits.getdata(src / "host_cutout.fits").shape
    if not (1 <= xp <= nx and 1 <= yp <= ny):
        return None
    if use_psf:
        return (f"# Component number: 0\n 0) psf  # star\n 1) {xp:.4f} {yp:.4f} 1 1\n"
                f" 3) {mag:.4f} 1\n Z) 0\n")
    # No PSF available: model the star as a free Gaussian (FWHM seed = PSF FWHM in px).
    fwhm_px = psf_fwhm_px(src / "proto_image.fits")
    return (f"# Component number: 0\n 0) gaussian  # star\n 1) {xp:.4f} {yp:.4f} 1 1\n"
            f" 3) {mag:.4f} 1\n 4) {fwhm_px:.3f} 1\n 9) 1.0000 1\n10) 0.0000 1\n Z) 0\n")


def _true_sigma(frb: str, src: Path) -> np.ndarray:
    """sigma = sqrt(PS1 variance), globally rescaled to the stamp's sky MAD.

    Production feeds the PS1 ``.wt`` (variance) map in as if it were invvar, so
    its host_sigma is k/sqrt(var); this rebuilds it with the correct spatial weighting.
    """
    import json
    meta = json.loads((src / "cutout_meta.json").read_text())
    x0, x1, y0, y1 = meta["cutout_bounds"]
    var = np.asarray(fits.getdata(input_paths(frb)[1]), float)[y0:y1, x0:x1]
    seg = np.asarray(fits.getdata(src / "segmentation_map.fits"))[y0:y1, x0:x1]
    data = np.asarray(fits.getdata(src / "host_cutout.fits"), float)
    good = np.isfinite(var) & (var > 0)
    sig = np.full(var.shape, 1e30)
    sig[good] = np.sqrt(var[good])
    sky = good & (seg == 0) & np.isfinite(data)
    pix = data[sky]
    mad = 1.4826 * np.median(np.abs(pix - np.median(pix)))
    sig[good] *= mad / np.median(sig[sky])
    return sig.astype(np.float32)


def _make_variant(frb: str, name: str, *, psf: bool, add_star: bool = False,
                  sky_free: bool = False, true_sigma: bool = False) -> Path | None:
    src = run_dir(frb)
    dest = OUT / frb / name
    dest.mkdir(parents=True, exist_ok=True)
    for f in GALFIT_INPUTS:
        if (src / f).is_file():
            shutil.copy2(src / f, dest / f)
    if true_sigma:
        fits.writeto(dest / "host_sigma.fits", _true_sigma(frb, src),
                     fits.getheader(src / "host_sigma.fits"), overwrite=True)
    header, comps = _components((src / "galfit.feedme").read_text(encoding="utf-8"))
    if not psf:
        header = _set_line(header, "D", "none  # no PSF convolution (Bhardwaj-like)")
    sky_i = next(i for i, c in enumerate(comps) if re.search(r"^\s*0\)\s*sky", c, re.M))
    cons = dest / "constraints.txt"
    lines = cons.read_text().splitlines() if cons.is_file() else []
    sky_line = re.compile(rf"^\s*{sky_i + 1}\s+1\s")
    if sky_free:
        lines = [ln for ln in lines if not sky_line.match(ln)]
    if add_star:
        star = _star_component_text(frb, src, use_psf=psf)
        if star is None:
            print(f"{frb}/{name}: star position not available, skipping")
            return None
        comps.insert(sky_i, star)
        # Phase 3a masks the star on PS1 (it deblends there); free its mask island.
        from scipy.ndimage import label
        xs, ys = map(float, re.search(r"^\s*1\)\s+(\S+)\s+(\S+)", star, re.M).groups())
        mask = fits.getdata(dest / "host_mask.fits")
        lab, _ = label(mask > 0)
        isl = lab[int(round(ys)) - 1, int(round(xs)) - 1]
        if isl:
            mask = mask.copy()
            mask[lab == isl] = 0
            fits.writeto(dest / "host_mask.fits", mask, fits.getheader(dest / "host_mask.fits"),
                         overwrite=True)
        lines = [re.sub(rf"^(\s*){sky_i + 1}(\s+1\s)", rf"\g<1>{sky_i + 2}\2", ln)
                 if sky_line.match(ln) else ln for ln in lines]
    # GALFIT cannot parse CRLF constraint files; keep LF like production.
    if cons.is_file():
        cons.write_text("\n".join(lines) + "\n", newline="\n")
    (dest / "galfit.feedme").write_text(_assemble(header, comps), encoding="utf-8", newline="\n")
    ok = _run_galfit(dest)
    print(f"{frb}/{name}: GALFIT {'OK' if ok else 'FAILED'}", flush=True)
    return dest


def cmd_variants(frbs: list[str]) -> None:
    for frb in frbs:
        if not (run_dir(frb) / "galfit.feedme").is_file():
            print(f"{frb}: no PS1 production run, skipping variants")
            continue
        _make_variant(frb, "nopsf", psf=False)
        _make_variant(frb, "psf_skyfree", psf=True, sky_free=True)
        _make_variant(frb, "nopsf_skyfree", psf=False, sky_free=True)
        _make_variant(frb, "psf_truesigma", psf=True, true_sigma=True)
        _make_variant(frb, "nopsf_truesigma", psf=False, true_sigma=True)
        if (WINNING_PSF_LEG / frb / "sandbox" / "galfit.01").is_file() and \
                _seed_row(frb).get("q_winning_leg") == "psf":
            _make_variant(frb, "psf_leg", psf=True, add_star=True)
            _make_variant(frb, "psf_leg_nopsf", psf=False, add_star=True)


# ---------------------------------------------------------------- measurements
def psf_fwhm_px(proto: Path) -> float:
    """FWHM of an elliptical 2D Gaussian fit to a PSFEx proto_image (geometric mean axis)."""
    from astropy.modeling import fitting, models
    img = np.asarray(fits.getdata(proto), float)
    img = np.squeeze(img)
    ny, nx = img.shape
    yy, xx = np.mgrid[:ny, :nx]
    yc, xc = np.unravel_index(np.nanargmax(img), img.shape)
    g0 = models.Gaussian2D(img.max(), xc, yc, 1.5, 1.5, 0.0)
    g = fitting.LevMarLSQFitter()(g0, xx, yy, img)
    s = math.sqrt(abs(g.x_stddev.value * g.y_stddev.value))
    return 2.0 * math.sqrt(2.0 * math.log(2.0)) * s


def psfex_fwhm_px(xml: Path) -> float:
    try:
        import xml.etree.ElementTree as ET
        root = ET.parse(xml).getroot()
        fields = root.findall(".//{*}FIELD") or root.findall(".//FIELD")
        tds = root.findall(".//{*}TR/{*}TD") or root.findall(".//TR/TD")
        names = [f.attrib.get("name", "") for f in fields]
        return float(tds[names.index("FWHM_FromFluxRadius_Mean")].text)
    except Exception:
        return float("nan")


def plate_scale(fits_path: Path) -> float:
    from astropy.wcs.utils import proj_plane_pixel_scales
    w = WCS(fits.getheader(fits_path)).celestial
    return float(np.mean(proj_plane_pixel_scales(w)) * 3600.0)


def host_fit(wkdir: Path) -> dict:
    """Component-1 Sérsic parameters from GALFIT's galfit.01 (final iteration)."""
    g01 = wkdir / "galfit.01"
    if not g01.is_file():
        return {}
    text = g01.read_text(encoding="utf-8", errors="replace")
    _, comps = _components(text)
    host = next((c for c in comps if re.search(r"^\s*0\)\s*sersic", c, re.M)), None)
    if host is None:
        return {}
    def val(k):
        m = re.search(rf"^\s*{k}\)\s+(\S+)", host, re.M)
        return float(m.group(1)) if m else float("nan")
    chi = re.search(r"Chi\^2/nu\s*=\s*([0-9.eE+-]+)", text)
    return {"re_px": val(4), "n": val(5), "q": val(9), "pa": val(10), "mag": val(3),
            "chi2nu": float(chi.group(1)) if chi else float("nan")}


def _seed_row(frb: str) -> dict:
    s = pd.read_csv(SEED, dtype={"frb": str})
    return s.loc[s.frb == frb].iloc[0].to_dict()


def cmd_table() -> None:
    import json
    nan = float("nan")
    rows = []
    for frb in FRBS:
        seed = _seed_row(frb)
        prod = PROD / f"{frb}_all"
        ps = run_dir(frb)
        fp, fn = host_fit(ps), host_fit(OUT / frb / "nopsf")
        var = {v: host_fit(OUT / frb / v).get("q", nan)
               for v in ("psf_skyfree", "nopsf_skyfree", "psf_truesigma", "nopsf_truesigma",
                         "psf_leg", "psf_leg_nopsf")}
        summ = json.loads((ps / "pipeline_summary.json").read_text())["galfit_host"]
        scale_ps1 = plate_scale(ps / "host_cutout.fits")
        scale_ls = plate_scale(prod / "host_cutout.fits")
        q_ls = float(seed["q_winning"])
        q_p, q_n = fp.get("q", nan), fn.get("q", nan)
        i_ls, i_p, i_n = hubble_i(q_ls), hubble_i(q_p), hubble_i(q_n)
        rows.append({
            "frb": frb,
            "q_ls10_psf": q_ls,
            "q_ps1_psf": q_p,
            "q_ps1_nopsf": q_n,
            "n_ps1_psf": fp.get("n", nan),
            "re_ps1_psf_arcsec": fp.get("re_px", nan) * scale_ps1,
            "fwhm_ls10_arcsec": psf_fwhm_px(prod / "proto_image.fits") * scale_ls,
            "fwhm_ps1_arcsec": psf_fwhm_px(ps / "proto_image.fits") * scale_ps1,
            "i_bhardwaj2024": seed["i_bhardwaj2024"],
            "i_ls10_psf_deg": i_ls,
            "i_ps1_psf_deg": i_p,
            "i_ps1_nopsf_deg": i_n,
            "di_ps1psf_minus_ls10_deg": i_p - i_ls,
            "di_ps1nopsf_minus_ps1psf_deg": i_n - i_p,
            "q_ls10_leg": seed["q_winning_leg"],
            "q_ls10_psf_err": seed["q_winning_err"],
            "q_ps1_psf_err": summ.get("b_a_err"),
            "chi2nu_ps1_psf": fp.get("chi2nu", nan),
            "chi2nu_ps1_nopsf": fn.get("chi2nu", nan),
            "n_ps1_nopsf": fn.get("n", nan),
            "re_ps1_nopsf_arcsec": fn.get("re_px", nan) * scale_ps1,
            "fwhm_ls10_psfex_arcsec": psfex_fwhm_px(prod / "psfex.xml") * scale_ls,
            "fwhm_ps1_psfex_arcsec": psfex_fwhm_px(ps / "psfex.xml") * scale_ps1,
            "q_ps1_psf_skyfree": var["psf_skyfree"],
            "q_ps1_nopsf_skyfree": var["nopsf_skyfree"],
            "q_ps1_psf_truesigma": var["psf_truesigma"],
            "q_ps1_nopsf_truesigma": var["nopsf_truesigma"],
            "q_ps1_psf_starcomp": var["psf_leg"],
            "q_ps1_nopsf_starcomp": var["psf_leg_nopsf"],
        })
    df = pd.DataFrame(rows)
    num = df.select_dtypes("number").columns
    df[num] = df[num].round(4)
    dest = BASE / "ps1_refits.csv"
    df.to_csv(dest, index=False)
    print(df.to_string(index=False))
    print(f"wrote {dest.relative_to(REPO)}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=["download", "pipeline", "variants", "table", "all"])
    ap.add_argument("--frb", nargs="+", default=FRBS)
    ap.add_argument("--force", action="store_true", help="re-download PS1 cutouts")
    a = ap.parse_args()
    if a.step in ("download", "all"):
        cmd_download(a.frb, a.force)
    if a.step in ("pipeline", "all"):
        cmd_pipeline(a.frb)
    if a.step in ("variants", "all"):
        cmd_variants(a.frb)
    if a.step in ("table", "all"):
        cmd_table()


if __name__ == "__main__":
    main()
