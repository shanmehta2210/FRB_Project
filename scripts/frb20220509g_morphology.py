"""FRB 20220509G host in deeper imaging: is it early-type (Sharma+23, Connor+23) or a barred spiral (Bhardwaj+24, App. D)?

Data: CFIS r band, 3 x 100 s MegaCam exposures (CADC programme 21AP99, public since 2023-08), reprojected onto the
first frame and averaged with outlier rejection. Gaia DR3 identifies point sources, which are masked.
Tests:
  unsharp mask        image minus a 1.5" Gaussian smoothing: bars and arms appear as elongated ridges
  isophote residual   image minus a photutils elliptical-isophote model: non-axisymmetric structure
  profile fits        PSF-convolved single Sersic vs Sersic bulge + exponential disk (Gaussian PSF from field stars)

Outputs: Publication/materials/appendix/morphology/frb20220509g/{cfis_morphology.pdf, .png, fits.csv, gaia.csv}
"""
from __future__ import annotations

import urllib.request
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.convolution import Gaussian2DKernel, convolve_fft
from astropy.io import fits
from astropy.modeling.models import Sersic2D
from astropy.stats import sigma_clipped_stats
from astropy.visualization import AsinhStretch, ImageNormalize, PercentileInterval
from astropy.wcs import WCS
from PIL import Image
from photutils.isophote import Ellipse, EllipseGeometry, build_ellipse_model
from photutils.segmentation import SourceCatalog, detect_sources
from reproject import reproject_interp
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Publication" / "materials" / "appendix" / "morphology" / "frb20220509g"
CACHE = ROOT / "Publication" / "materials" / "galaxy_cache"
RA, DEC = 282.67481, 70.24274  # LS DR9 Tractor centre of the host (Gaia DR3 2265362396989859456)
FRAMES = (2617289, 2617298, 2617326)
PIX = 0.185
CADC = "https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/caom2ops/sync?ID=cadc%3ACFHT%2F{}p.fits.fz&POS=CIRCLE+{}+{}+0.0083"
LS_GALFIT_N = 0.5


def frames() -> list:
    out = []
    for i in FRAMES:
        p = CACHE / f"cfis_20220509G_{i}.fits"
        if not p.exists():
            p.write_bytes(urllib.request.urlopen(CADC.format(i, RA, DEC), timeout=120).read())
        out.append(fits.open(p)[0])
    return out


def stack() -> tuple[np.ndarray, WCS, float]:
    fs = frames()
    h0 = fs[0].header
    ims = []
    for f in fs:
        im = reproject_interp((f.data.astype(float), f.header), h0)[0]
        ims.append(im - sigma_clipped_stats(im, sigma=3)[1])
    cube = np.array(ims)
    med = np.nanmedian(cube, 0)
    rms = sigma_clipped_stats(med, sigma=3)[2]
    cube[np.abs(cube - med) > 5 * rms + 0.3 * np.abs(med)] = np.nan
    img = np.nanmean(cube, 0)
    return img, WCS(h0), sigma_clipped_stats(img, sigma=3)[2]


def gaia() -> pd.DataFrame:
    p = OUT / "gaia.csv"
    if not p.exists():
        from astroquery.gaia import Gaia
        q = (f"SELECT source_id, ra, dec, phot_g_mean_mag, bp_rp, parallax, parallax_over_error, pmra, pmdec, ruwe "
             f"FROM gaiadr3.gaia_source WHERE CONTAINS(POINT(ra, dec), CIRCLE({RA}, {DEC}, {40 / 3600}))=1")
        Gaia.launch_job(q).get_results().to_pandas().to_csv(p, index=False)
    return pd.read_csv(p)


def psf_fwhm(img: np.ndarray, rms: float, x0: float, y0: float) -> float:
    seg = detect_sources(img, 10 * rms, npixels=10)
    sc = SourceCatalog(img, seg)
    cat = pd.DataFrame(dict(x=sc.xcentroid, y=sc.ycentroid, ecc=np.asarray(sc.eccentricity), area=np.asarray(sc.area),
                            fwhm=np.asarray(sc.fwhm)))
    stars = cat[(np.hypot(cat.x - x0, cat.y - y0) > 40) & (cat.ecc < 0.4) & (cat.area < 400)]
    return float(np.median(stars.fwhm.nsmallest(5)))


def masks(img, rms, wcs, x0, y0, g, fwhm) -> np.ndarray:
    seg = detect_sources(gaussian_filter(img, 1.5), 2 * rms, npixels=8)
    host = seg.data[int(round(y0)), int(round(x0))]
    m = (seg.data > 0) & (seg.data != host)
    yy, xx = np.indices(img.shape)
    for r in g.itertuples():
        x, y = wcs.world_to_pixel_values(r.ra, r.dec)
        if np.hypot(x - x0, y - y0) > 3 and (r.parallax_over_error > 3 or abs(r.pmra) + abs(r.pmdec) > 5):
            m |= np.hypot(xx - x, yy - y) < 2 * fwhm
    return m | ~np.isfinite(img)


def sersic_model(p, shape, kernel, n_fixed=None):
    yy, xx = np.indices(shape)
    amp, reff, n, x0, y0, ell, th = p[:7]
    m = Sersic2D(amp, reff, n_fixed if n_fixed is not None else n, x0, y0, ell, th)(xx, yy)
    if len(p) > 7:
        amp2, reff2, ell2, th2 = p[7:]
        m = m + Sersic2D(amp2, reff2, 1.0, x0, y0, ell2, th2)(xx, yy)
    return convolve_fft(m, kernel, normalize_kernel=True)


def fit(img, mask, rms, kernel, x0, y0, two):
    ok = ~mask
    w = 1 / np.sqrt(rms ** 2 + np.clip(img, 0, None) / 1.66)
    p0 = [img[int(y0), int(x0)] / 5, 15, 2.0, x0, y0, 0.3, 0.7]
    lo = [0, 1, 0.3, x0 - 5, y0 - 5, 0, -np.pi]
    hi = [np.inf, 150, 8, x0 + 5, y0 + 5, 0.9, np.pi]
    if two:
        p0 = [img[int(y0), int(x0)] / 2, 4, 2.0, x0, y0, 0.2, 0.7, img[int(y0), int(x0)] / 20, 25, 0.5, 0.7]
        lo += [0, 3, 0, -np.pi]
        hi += [np.inf, 200, 0.9, np.pi]
    r = least_squares(lambda p: ((img - sersic_model(p, img.shape, kernel)) * w)[ok], p0, bounds=(lo, hi), max_nfev=4000)
    model = sersic_model(r.x, img.shape, kernel)
    chi2 = float(np.sum(((img - model) * w)[ok] ** 2) / (ok.sum() - len(p0)))
    out = dict(model=model, p=r.x, chi2_nu=chi2, bic=float(np.sum(((img - model) * w)[ok] ** 2) + len(p0) * np.log(ok.sum())))
    if two:
        yy, xx = np.indices(img.shape)
        b = Sersic2D(*r.x[:7])(xx, yy).sum()
        d = Sersic2D(r.x[7], r.x[8], 1.0, r.x[3], r.x[4], r.x[9], r.x[10])(xx, yy).sum()
        out["bt"] = b / (b + d)
    return out


def along(model, iso) -> np.ndarray:
    t = np.linspace(0, 2 * np.pi, 72)
    out = []
    for x0, y0, a, eps, pa in zip(iso.x0, iso.y0, iso.sma, iso.eps, iso.pa):
        b = a * (1 - eps)
        x = x0 + a * np.cos(t) * np.cos(pa) - b * np.sin(t) * np.sin(pa)
        y = y0 + a * np.cos(t) * np.sin(pa) + b * np.sin(t) * np.cos(pa)
        out.append(map_coordinates(model, [y, x], order=1).mean())
    return np.array(out)


def main() -> None:
    warnings.filterwarnings("ignore")
    OUT.mkdir(parents=True, exist_ok=True)
    full, wcs, rms = stack()
    x0, y0 = (float(v) for v in wcs.world_to_pixel_values(RA, DEC))
    fwhm = psf_fwhm(full, rms, x0, y0)
    g = gaia()
    mask_full = masks(full, rms, wcs, x0, y0, g, fwhm)
    h = 70
    sl = (slice(int(y0) - h, int(y0) + h), slice(int(x0) - h, int(x0) + h))
    img, mask = full[sl], mask_full[sl]
    cx, cy = x0 - (int(x0) - h), y0 - (int(y0) - h)
    yy, xx = np.indices(img.shape)
    peak = np.unravel_index(np.nanargmax(np.where(mask | (np.hypot(xx - cx, yy - cy) > 15), -np.inf, img)), img.shape)
    cy, cx = float(peak[0]), float(peak[1])
    kernel = Gaussian2DKernel(fwhm / 2.355, x_size=31, y_size=31)

    unsharp = img - gaussian_filter(np.where(mask, np.nanmedian(img), img), 1.5 / PIX)
    filled = np.where(mask, np.nan, img)
    geo = EllipseGeometry(cx, cy, sma=6, eps=0.2, pa=0.5)
    iso = Ellipse(np.ma.masked_array(filled, mask | ~np.isfinite(filled)), geo).fit_image(maxsma=60, nclip=2)
    iso_model = build_ellipse_model(img.shape, iso)
    iso_resid = np.where(iso_model > 0, img - iso_model, np.nan)

    one = fit(img, mask, rms, kernel, cx, cy, two=False)
    bd = fit(img, mask, rms, kernel, cx, cy, two=True)
    res = pd.DataFrame([
        dict(model="single Sersic", n=one["p"][2], reff_arcsec=one["p"][1] * PIX, bt=np.nan, chi2_nu=one["chi2_nu"], bic=one["bic"]),
        dict(model="Sersic bulge + exponential disk", n=bd["p"][2], reff_arcsec=bd["p"][1] * PIX, bt=bd["bt"],
             chi2_nu=bd["chi2_nu"], bic=bd["bic"], disk_reff_arcsec=bd["p"][8] * PIX),
    ])
    res["psf_fwhm_arcsec"] = fwhm * PIX
    res.to_csv(OUT / "fits.csv", index=False)
    print(res.to_string())
    print(g[["phot_g_mean_mag", "bp_rp", "parallax", "parallax_over_error", "pmra", "pmdec"]].to_string())

    fig, ax = plt.subplots(2, 3, figsize=(11, 7.4))
    ax = ax.ravel()
    jpg = CACHE / "20220509G_ls-dr10.jpg"
    if jpg.exists():
        ax[0].imshow(np.asarray(Image.open(jpg).convert("RGB")))
    ax[0].set_title("(a) Legacy Survey $grz$ (our imaging)", fontsize=9)
    norm = ImageNormalize(img[~mask], interval=PercentileInterval(99.8), stretch=AsinhStretch(0.02))
    ax[1].imshow(img, origin="lower", cmap="gray_r", norm=norm)
    ax[1].contour(mask, levels=[0.5], colors="c", linewidths=0.5)
    ax[1].set_title(f"(b) CFIS $r$, 3$\\times$100 s (FWHM {fwhm * PIX:.2f}\")", fontsize=9)
    for a, im, t in ((ax[2], unsharp, "(c) CFIS unsharp mask (1.5\" Gaussian)"),
                     (ax[3], iso_resid, "(d) CFIS minus elliptical-isophote model"),
                     (ax[4], np.where(mask, np.nan, img - bd["model"]), "(e) CFIS minus bulge + disk model")):
        v = np.nanpercentile(np.abs(im[~mask & np.isfinite(im) & (np.hypot(xx - cx, yy - cy) < 40)]), 98)
        a.imshow(np.where(mask, np.nan, im), origin="lower", cmap="RdBu_r", vmin=-v, vmax=v)
        a.set_title(t, fontsize=9)
    for a in ax[:5]:
        a.set_xticks([]); a.set_yticks([])
    sma = iso.sma * PIX
    good = iso.intens > 0
    ax[5].semilogy(sma[good], iso.intens[good], "k.", ms=4, label="CFIS isophotes")
    for f, lab, c in ((one, f"single Sersic, n = {one['p'][2]:.1f}", "tab:red"),
                      (bd, f"bulge n = {bd['p'][2]:.1f} + disk, B/T = {bd['bt']:.2f}", "tab:blue")):
        ax[5].semilogy(sma[good], along(f["model"], iso)[good], color=c, lw=1, label=lab)
    ax[5].set(xlabel="semi-major axis (arcsec)", ylabel="intensity (counts/pix)", ylim=(rms / 2, None),
              title="(f) Surface-brightness profile")
    ax[5].text(0.97, 0.75, f"LS GALFIT (ours): n = {LS_GALFIT_N}\n$\\Delta$BIC (1 comp $-$ 2 comp) = {one['bic'] - bd['bic']:.0f}",
               transform=ax[5].transAxes, ha="right", fontsize=8)
    ax[5].legend(fontsize=7, loc="upper right")
    fig.suptitle("Host of FRB 20220509G in CFIS $r$ (the deeper imaging used by Bhardwaj et al. 2024)", fontsize=10)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"cfis_morphology.{ext}", dpi=200)


if __name__ == "__main__":
    main()
