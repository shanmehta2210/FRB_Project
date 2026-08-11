"""Check 6 — PSF leakage.

The direct test of whether the measured shapes are partly tracking the
instrument rather than the galaxies. Per-host this only gathers the PSF's own
shape; the correlations against the cohort are built in ``aggregate.py``, since
a leakage test is meaningless one galaxy at a time.

Two independent PSF ellipticity estimates are recorded: PSFEx's Stokes
components from ``psfex.xml``, and the flux-weighted second moments of
``proto_image.fits`` — the stamp GALFIT actually convolved with. They probe
different things (PSFEx's model fit versus the delivered stamp), so a
disagreement is itself informative.
"""

from __future__ import annotations

import math
import os

import numpy as np

import vercommon as vc

NAME = "psf"


def run(host: vc.HostData, outdir: str) -> dict:
    xml = vc.parse_psfex_xml(os.path.join(host.dir, "psfex.xml"))
    e1 = xml.get("Ellipticity1_Mean", float("nan"))
    e2 = xml.get("Ellipticity2_Mean", float("nan"))
    if not (math.isfinite(e1) and math.isfinite(e2)):
        e1 = xml.get("Ellipticity1_PixelFree_Mean", float("nan"))
        e2 = xml.get("Ellipticity2_PixelFree_Mean", float("nan"))

    out: dict = {
        "psf_e1": e1,
        "psf_e2": e2,
        "psf_ellipticity_xml": xml.get("Ellipticity_Mean", float("nan")),
        "psf_ellipticity_stdev": xml.get("Ellipticity_StDev", float("nan")),
        "psf_fwhm_xml_px": xml.get("FWHM_FromFluxRadius_Mean", float("nan")),
        "psf_nstars": xml.get("NStars_Accepted_Total", float("nan")),
        "psfex_chi2": xml.get("Chi2_Mean", float("nan")),
    }

    if math.isfinite(e1) and math.isfinite(e2):
        out["psf_ellipticity"] = math.hypot(e1, e2)
        # Stokes angle is measured from +x; GALFIT PA is measured from +y.
        out["psf_pa_deg"] = vc.wrap_pa(0.5 * math.degrees(math.atan2(e2, e1)) - 90.0)
    else:
        out["psf_ellipticity"] = out["psf_ellipticity_xml"]
        out["psf_pa_deg"] = float("nan")

    mom = vc.psf_second_moments(host.psf) if host.psf.size > 1 else {}
    out["psf_e1_moments"] = mom.get("e1", float("nan"))
    out["psf_e2_moments"] = mom.get("e2", float("nan"))
    out["psf_ellipticity_moments"] = mom.get("ellipticity", float("nan"))
    out["psf_pa_moments_deg"] = mom.get("pa_deg", float("nan"))
    out["psf_fwhm_moments_px"] = mom.get("fwhm_px", float("nan"))

    # Fall back to the moment-based PA when PSFEx did not write the Stokes pair.
    if not math.isfinite(out["psf_pa_deg"]):
        out["psf_pa_deg"] = out["psf_pa_moments_deg"]
        out["psf_pa_source"] = "moments"
    else:
        out["psf_pa_source"] = "psfex_xml"

    out["psf_fwhm_px"] = host.psf_fwhm
    out["q_host"] = host.q
    out["pa_host_deg"] = vc.wrap_pa(host.pa)
    out["re_px"] = host.re
    out["re_over_fwhm"] = (host.re / host.psf_fwhm
                           if math.isfinite(host.psf_fwhm) and host.psf_fwhm > 0
                           else float("nan"))
    out["fwhm_over_re"] = (host.psf_fwhm / host.re
                           if math.isfinite(host.re) and host.re > 0
                           else float("nan"))

    # Wrapped to [0, 90): PA is defined modulo 180, and the alignment question
    # is direction-agnostic, so 89 deg and 91 deg are the same misalignment.
    dpa = out["pa_host_deg"] - out["psf_pa_deg"]
    out["dpa_host_psf_deg"] = (abs(vc.wrap_pa(dpa)) if math.isfinite(dpa)
                               else float("nan"))

    # Ellipticity the PSF alone would impose on a perfectly round source of this
    # size: the expected size of any leakage, for scale.
    ell = out["psf_ellipticity"]
    if math.isfinite(ell) and math.isfinite(out["fwhm_over_re"]):
        out["psf_leak_scale"] = ell * out["fwhm_over_re"] ** 2
    else:
        out["psf_leak_scale"] = float("nan")

    if np.count_nonzero(np.isfinite(host.psf)) < 4:
        out["status"] = "no_psf_stamp"
    return out
