"""Inject a known axis-ratio error and confirm the Fourier estimator recovers it.

This is the check that makes check 3 trustworthy. Everything else about the
m=2 analysis is interpretation; this is the only place the estimator is held
against a known answer, including its sign.

Run standalone (``python test_fourier_recovery.py``) or under pytest.
"""

from __future__ import annotations

import math
import os
import sys

import numpy as np
from astropy.io import fits
from scipy.signal import fftconvolve

VER_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for path in (VER_DIR, os.path.join(VER_DIR, "checks")):
    if path not in sys.path:
        sys.path.insert(0, path)

import vercommon as vc  # noqa: E402
from checks.fourier import decompose  # noqa: E402

REFERENCE_FRB = "20240210A"


def sersic_image(shape, xc, yc, q, pa, re, n, total=1.0):
    """Analytic Sersic on the same elliptical coordinates the estimator uses."""
    a, _ = vc.elliptical_coords(shape, xc, yc, q, pa)
    bn = 2.0 * n - 1.0 / 3.0 + 4.0 / (405.0 * n)
    img = np.exp(-bn * ((np.maximum(a, 1e-6) / re) ** (1.0 / n) - 1.0))
    return img * (total / img.sum())


def _psf(shape_hint: int = 25):
    path = os.path.join(vc.host_dir(REFERENCE_FRB), "proto_image.fits")
    if os.path.isfile(path):
        psf = np.squeeze(np.asarray(fits.getdata(path), dtype=float))
        psf = np.clip(np.nan_to_num(psf), 0.0, None)
        if psf.sum() > 0:
            return psf / psf.sum(), "proto_image"
    yy, xx = np.mgrid[0:shape_hint, 0:shape_hint]
    c = (shape_hint - 1) / 2.0
    sigma = 1.5
    psf = np.exp(-((xx - c) ** 2 + (yy - c) ** 2) / (2 * sigma**2))
    return psf / psf.sum(), "gaussian_fallback"


def recover(dq_true: float = 0.0, dpa_true: float = 0.0, q0: float = 0.60,
            pa: float = 30.0, re: float = 20.0, n: float = 1.0, size: int = 241,
            calibrated: bool = True) -> tuple[float, float]:
    """Inject a known (dq, dPA) and return what the estimator recovers."""
    psf, _ = _psf()
    fwhm = vc.psf_second_moments(psf).get("fwhm_px", 3.0)
    shape = (size, size)
    xc = yc = (size - 1) / 2.0

    def conv(q, pa_deg):
        return fftconvolve(sersic_image(shape, xc, yc, q, pa_deg, re, n), psf,
                           mode="same")

    model = conv(q0, pa)
    truth = conv(q0 + dq_true, pa + dpa_true)

    kernel_q = kernel_pa = None
    if calibrated:
        eps_q, eps_pa = 0.02, 2.0
        kernel_q = (conv(q0 + eps_q, pa) - conv(q0 - eps_q, pa)) / (2 * eps_q)
        kernel_pa = (conv(q0, pa + eps_pa) - conv(q0, pa - eps_pa)) / (2 * eps_pa)

    sigma = np.full(shape, 1e-6 * float(np.max(model)))
    prof = decompose(
        truth - model, model, sigma, np.ones(shape, dtype=bool),
        xc=xc, yc=yc, q=q0, pa=pa, re=re, fwhm=fwhm,
        kernel_q=kernel_q, kernel_pa=kernel_pa,
    )
    inner = np.isfinite(prof["a_mid_re"]) & (prof["a_mid_re"] <= 2.0)
    dq, _ = vc.weighted_mean(prof["dq"][inner], prof["dq_err"][inner])
    dpa, _ = vc.weighted_mean(prof["dpa_deg"][inner], prof["dpa_err_deg"][inner])
    return float(dq), float(dpa)


def test_recovers_injected_dq():
    dq, _ = recover(dq_true=0.05)
    assert math.isfinite(dq), "estimator returned no finite dq"
    assert abs(dq - 0.05) < 0.005, f"recovered dq={dq:.4f}, expected 0.05"


def test_recovers_injected_dpa():
    _, dpa = recover(dpa_true=3.0)
    assert abs(dpa - 3.0) < 0.4, f"recovered dPA={dpa:.3f} deg, expected 3.0"


def test_modes_do_not_mix():
    """dq must not leak into dPA or vice versa; that is the whole premise."""
    dq_only, dpa_leak = recover(dq_true=0.05)
    dq_leak, dpa_only = recover(dpa_true=3.0)
    assert abs(dpa_leak) < 0.3, f"dq=0.05 produced a spurious dPA={dpa_leak:.3f} deg"
    assert abs(dq_leak) < 0.005, f"dPA=3 deg produced a spurious dq={dq_leak:.4f}"


def test_sign_and_symmetry():
    plus, _ = recover(dq_true=+0.04)
    minus, _ = recover(dq_true=-0.04)
    assert plus > 0 > minus, f"sign wrong: +0.04 -> {plus:.4f}, -0.04 -> {minus:.4f}"
    assert abs(abs(plus) - abs(minus)) < 0.006, "estimator strongly asymmetric"


def test_null_case():
    dq, dpa = recover()
    assert abs(dq) < 1e-3, f"non-zero dq={dq:.2e} on an identical pair"
    assert abs(dpa) < 1e-2, f"non-zero dPA={dpa:.2e} on an identical pair"


def test_calibration_beats_analytic():
    """The PSF-blind analytic conversion is biased; the calibrated one is not."""
    cal, _ = recover(dq_true=0.05, calibrated=True)
    raw, _ = recover(dq_true=0.05, calibrated=False)
    assert abs(cal - 0.05) < abs(raw - 0.05), (
        f"calibration did not help: calibrated={cal:.4f}, analytic={raw:.4f}"
    )


if __name__ == "__main__":
    psf, source = _psf()
    print(f"PSF source: {source}, FWHM = {vc.psf_second_moments(psf)['fwhm_px']:.2f} px")
    print("  injected dq   calibrated   analytic")
    for truth in (0.0, 0.02, 0.05, -0.04, 0.10):
        cal, _ = recover(dq_true=truth, calibrated=True)
        raw, _ = recover(dq_true=truth, calibrated=False)
        print(f"   {truth:+.3f}        {cal:+.4f}      {raw:+.4f}")
    print("  injected dPA  recovered")
    for truth in (0.0, 3.0, -5.0):
        _, dpa = recover(dpa_true=truth)
        print(f"   {truth:+.1f}          {dpa:+.3f}")
    for fn in (test_recovers_injected_dq, test_recovers_injected_dpa,
               test_modes_do_not_mix, test_sign_and_symmetry, test_null_case,
               test_calibration_beats_analytic):
        fn()
        print(f"  pass: {fn.__name__}")
    print("all assertions passed")
