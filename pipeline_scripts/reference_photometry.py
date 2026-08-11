"""
Trusted reference photometry for the FRB host, independent of the pipeline ZP.

Why
---
When Phase 2 zero-point calibration fails (too few PS1/Legacy calibration
stars) GALFIT still runs with a fallback J) zeropoint, so ``fit.log`` carries a
wrong host magnitude while all structural parameters (re, n, b/a, PA and hence
inclination) remain valid. Rather than discarding those hosts, we always query
an external survey for the host's r-band magnitude and store it alongside the
pipeline products. Downstream, ``mag_final`` is:

* the reference-survey magnitude when the pipeline ZP is untrusted **or**
  GALFIT ``mag_err > 1``;
* the GALFIT magnitude when the pipeline zero-point is trusted
  (``zp_aper`` finite and ``zp_aper_std`` <= ``MAX_ZP_APER_STD``) and
  ``mag_err <= 1``;
* empty (``unavailable``) when neither a survey match nor a trusted / rescaled
  pipeline ZP is available — those FRBs are listed under
  ``docs/WEAK_ASSOCIATIONS_PRODUCTION67.md`` §A.

Survey order
------------
1. **Legacy Surveys DR10 Tractor** (NOIRLab Data Lab TAP): model ``flux_r``
   (nanomaggies) -> AB mag; error from ``flux_ivar_r``. Nearest source within
   the search radius, preferring non-PSF (galaxy) Tractor types.
2. **PS1 DR1** (VizieR ``II/349``): Kron ``rKmag``/``e_rKmag`` (galaxy-suited),
   falling back to mean PSF ``rmag``/``e_rmag``. Used when the position is
   outside the LS DR10 footprint or the TAP query fails.

Both fail (e.g. dec < -30 with no LS coverage) -> status is recorded and
``mag_final`` falls back to the GALFIT value or stays empty.

Usage
-----
Called automatically by ``master_run.py`` after the fit phases. Manual backfill
of existing Output folders (writes ``reference_photometry.json`` and updates
``pipeline_summary.json`` in each folder):

    python pipeline_scripts/reference_photometry.py --backfill
    python pipeline_scripts/reference_photometry.py --backfill --frb 20230913 --force
"""

from __future__ import annotations

import argparse
import datetime
import json
import math
import time
from pathlib import Path

# Pipeline ZP is considered unreliable above this scatter (mag). Catches e.g.
# a 3-star calibration with zp_aper_std ~ 52 mag (20221116A).
MAX_ZP_APER_STD = 1.0

# Prefer the external survey mag whenever GALFIT's formal mag error exceeds
# this threshold, even if the pipeline aperture ZP itself looks fine.
MAX_GALFIT_MAG_ERR = 1.0

DEFAULT_SEARCH_RADIUS_ARCSEC = 3.0

_LS_TAP_URL = "https://datalab.noirlab.edu/tap"
_LS_SURVEY = "LS_DR10"
_PS1_SURVEY = "PS1_DR1"


def _utcnow() -> str:
    return datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")


def _retry(fn, attempts: int = 2, delay_s: float = 3.0, label: str = "query"):
    last = None
    for i in range(attempts):
        try:
            return fn()
        except Exception as exc:  # network / service errors
            last = exc
            if i + 1 < attempts:
                time.sleep(delay_s)
    raise RuntimeError(f"{label} failed after {attempts} attempts: {last}")


def _sep_arcsec(ra1, dec1, ra2, dec2) -> float:
    """Small-angle separation (arcsec); fine for arcsec-scale radii."""
    cosd = math.cos(math.radians((dec1 + dec2) / 2.0))
    dra = (ra1 - ra2) * cosd
    ddec = dec1 - dec2
    return math.hypot(dra, ddec) * 3600.0


def _query_ls_dr10(ra: float, dec: float, radius_arcsec: float) -> dict | None:
    """Nearest LS DR10 Tractor source; prefers non-PSF (galaxy) types."""
    import numpy as np
    import pyvo

    r_deg = radius_arcsec / 3600.0
    cos_dec = math.cos(math.radians(max(-85.0, min(85.0, dec))))
    ra_half = r_deg / max(cos_dec, 1e-6)
    ra_min, ra_max = ra - ra_half, ra + ra_half
    if ra_min < 0:
        ra_clause = f"(ra > {ra_min + 360:.8f} OR ra < {ra_max:.8f})"
    elif ra_max > 360:
        ra_clause = f"(ra > {ra_min:.8f} OR ra < {ra_max - 360:.8f})"
    else:
        ra_clause = f"ra > {ra_min:.8f} AND ra < {ra_max:.8f}"
    query = (
        "SELECT ra, dec, type, flux_r, flux_ivar_r FROM ls_dr10.tractor WHERE "
        + ra_clause
        + f" AND dec > {dec - r_deg:.8f} AND dec < {dec + r_deg:.8f}"
    )
    service = pyvo.dal.TAPService(_LS_TAP_URL)
    tab = _retry(lambda: service.search(query).to_table(), label="LS DR10 TAP")
    if tab is None or len(tab) == 0:
        return None

    rows = []
    for r in tab:
        flux = float(r["flux_r"]) if r["flux_r"] is not None else float("nan")
        if not (math.isfinite(flux) and flux > 0):
            continue
        sra, sdec = float(r["ra"]), float(r["dec"])
        sep = _sep_arcsec(ra, dec, sra, sdec)
        if sep > radius_arcsec:
            continue
        ivar = float(r["flux_ivar_r"]) if r["flux_ivar_r"] is not None else 0.0
        mag = 22.5 - 2.5 * math.log10(flux)
        mag_err = None
        if ivar > 0:
            mag_err = 2.5 / math.log(10) / (flux * math.sqrt(ivar))
        ls_type = str(r["type"]).strip()
        rows.append((sep, sra, sdec, mag, mag_err, ls_type))
    if not rows:
        return None
    rows.sort(key=lambda t: t[0])
    galaxies = [t for t in rows if t[5].upper() != "PSF"]
    sep, sra, sdec, mag, mag_err, ls_type = (galaxies or rows)[0]
    return {
        "survey": _LS_SURVEY,
        "band": "r",
        "mag": round(mag, 4),
        "mag_err": round(mag_err, 4) if mag_err is not None else None,
        "mag_type": "tractor_model_flux_r",
        "ra_deg": sra,
        "dec_deg": sdec,
        "sep_arcsec": round(sep, 3),
        "ls_type": ls_type,
    }


def _query_ps1(ra: float, dec: float, radius_arcsec: float) -> dict | None:
    """Nearest PS1 DR1 source (VizieR II/349); Kron mag preferred."""
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier

    v = Vizier(
        columns=["RAJ2000", "DEJ2000", "rmag", "e_rmag", "rKmag", "e_rKmag"],
        row_limit=50,
    )
    center = SkyCoord(ra=ra, dec=dec, unit="deg")
    result = _retry(
        lambda: v.query_region(center, radius=radius_arcsec * u.arcsec, catalog="II/349"),
        label="PS1 VizieR",
    )
    if not result or len(result[0]) == 0:
        return None

    best = None
    for r in result[0]:
        try:
            sra, sdec = float(r["RAJ2000"]), float(r["DEJ2000"])
        except (TypeError, ValueError):
            continue
        sep = _sep_arcsec(ra, dec, sra, sdec)

        def _f(col):
            try:
                import numpy as np

                if np.ma.is_masked(r[col]):
                    return None
                val = float(r[col])
                return val if math.isfinite(val) else None
            except (TypeError, ValueError, KeyError):
                return None

        mag, mag_err, mag_type = _f("rKmag"), _f("e_rKmag"), "kron_rKmag"
        if mag is None:
            mag, mag_err, mag_type = _f("rmag"), _f("e_rmag"), "mean_psf_rmag"
        if mag is None:
            continue
        if best is None or sep < best[0]:
            best = (sep, sra, sdec, mag, mag_err, mag_type)
    if best is None:
        return None
    sep, sra, sdec, mag, mag_err, mag_type = best
    return {
        "survey": _PS1_SURVEY,
        "band": "r",
        "mag": round(mag, 4),
        "mag_err": round(mag_err, 4) if mag_err is not None else None,
        "mag_type": mag_type,
        "ra_deg": sra,
        "dec_deg": sdec,
        "sep_arcsec": round(sep, 3),
    }


def query_reference_mag(
    ra: float,
    dec: float,
    radius_arcsec: float = DEFAULT_SEARCH_RADIUS_ARCSEC,
) -> dict:
    """LS DR10 first, PS1 fallback. Never raises; ``status`` reports outcome."""
    out = {
        "query_ra_deg": ra,
        "query_dec_deg": dec,
        "search_radius_arcsec": radius_arcsec,
        "queried_utc": _utcnow(),
        "status": None,
    }
    errors = []
    for name, fn in (("LS DR10", _query_ls_dr10), ("PS1", _query_ps1)):
        try:
            hit = fn(ra, dec, radius_arcsec)
        except Exception as exc:
            errors.append(f"{name}: {exc}")
            continue
        if hit is not None:
            out.update(hit)
            out["status"] = "ok"
            return out
        errors.append(f"{name}: no source within {radius_arcsec:.1f} arcsec")
    out["status"] = "no_match: " + " | ".join(errors)
    return out


def resolve_host_coords(folder: Path) -> tuple[float, float, str] | None:
    """Best host position from a run folder (workdir or collected Output).

    Preference: Phase 3a host row (host_components.csv row 0, SExtractor
    windowed centroid) -> best AstroPath candidate -> None.
    """
    comp = folder / "host_components.csv"
    if comp.is_file():
        try:
            import pandas as pd

            df = pd.read_csv(comp, nrows=1)
            if not df.empty and {"ALPHAWIN_J2000", "DELTAWIN_J2000"} <= set(df.columns):
                ra = float(df.iloc[0]["ALPHAWIN_J2000"])
                dec = float(df.iloc[0]["DELTAWIN_J2000"])
                if math.isfinite(ra) and math.isfinite(dec):
                    return ra, dec, "host_components.csv"
        except Exception:
            pass
    post = folder / "astropath_posteriors.csv"
    if post.is_file():
        try:
            import pandas as pd

            df = pd.read_csv(post)
            if not df.empty and {"ra_deg", "dec_deg", "posterior_O"} <= set(df.columns):
                best = df.sort_values("posterior_O", ascending=False).iloc[0]
                ra, dec = float(best["ra_deg"]), float(best["dec_deg"])
                if math.isfinite(ra) and math.isfinite(dec):
                    return ra, dec, "astropath_posteriors.csv"
        except Exception:
            pass
    return None


def fetch_reference_photometry(
    folder: Path,
    fallback_ra: float | None = None,
    fallback_dec: float | None = None,
    radius_arcsec: float = DEFAULT_SEARCH_RADIUS_ARCSEC,
) -> dict:
    """Resolve host coords in ``folder``, query surveys, write JSON there."""
    resolved = resolve_host_coords(folder)
    if resolved is not None:
        ra, dec, coord_source = resolved
    elif fallback_ra is not None and fallback_dec is not None:
        ra, dec, coord_source = fallback_ra, fallback_dec, "frb_localization"
    else:
        ref = {"status": "error: no host coordinates available", "queried_utc": _utcnow()}
        (folder / "reference_photometry.json").write_text(
            json.dumps(ref, indent=2), encoding="utf-8"
        )
        return ref

    ref = query_reference_mag(ra, dec, radius_arcsec=radius_arcsec)
    ref["coord_source"] = coord_source
    (folder / "reference_photometry.json").write_text(
        json.dumps(ref, indent=2), encoding="utf-8"
    )
    return ref


def load_photometry(folder: Path) -> dict | None:
    """Pipeline photometry dict for a run folder.

    Prefers the ``photometry`` section of ``pipeline_summary.json``; falls back
    to ``zero_points.json`` (older Output folders predate the consolidated
    summary). Keys are normalised so ``zp_aper``/``zp_aper_std`` always exist.
    """
    for name, getter in (
        ("pipeline_summary.json", lambda d: d.get("photometry")),
        ("zero_points.json", lambda d: d),
    ):
        path = folder / name
        if not path.is_file():
            continue
        try:
            data = getter(json.loads(path.read_text(encoding="utf-8")))
        except (json.JSONDecodeError, OSError):
            continue
        if not isinstance(data, dict):
            continue
        phot = dict(data)
        if phot.get("zp_aper") is None:
            phot["zp_aper"] = phot.get("zp_aper_40px")
        if phot.get("zp_aper_std") is None:
            phot["zp_aper_std"] = phot.get("zp_aper_40px_std")
        if phot.get("zp_aper") is not None:
            return phot
    return None


def pipeline_zp_ok(photometry: dict | None) -> bool:
    """True when the pipeline aperture zero-point is trustworthy."""
    if not isinstance(photometry, dict):
        return False
    try:
        zp = float(photometry.get("zp_aper"))
    except (TypeError, ValueError):
        return False
    if not math.isfinite(zp):
        return False
    std = photometry.get("zp_aper_std")
    if std is not None:
        try:
            std = float(std)
            if math.isfinite(std) and std > MAX_ZP_APER_STD:
                return False
        except (TypeError, ValueError):
            pass
    return True


def galfit_zp_used(folder: Path) -> float | None:
    """The J) zeropoint GALFIT actually ran with, parsed from galfit.feedme."""
    feedme = folder / "galfit.feedme"
    if not feedme.is_file():
        return None
    try:
        for line in feedme.read_text(encoding="utf-8", errors="replace").splitlines():
            stripped = line.strip()
            if stripped.startswith("J)"):
                val = float(stripped.split()[1])
                return val if math.isfinite(val) else None
    except (OSError, ValueError, IndexError):
        pass
    return None


def _best_alternative_zp(photometry: dict | None) -> tuple[float, float, str] | None:
    """Most reliable non-aperture pipeline ZP: (zp, std, name) or None.

    zp_psf / zp_auto come from independent flux estimators; when the aperture
    ZP scatter blows up (e.g. one bad reference star dominating a 3-star
    calibration) one of them is often still tight.
    """
    if not isinstance(photometry, dict):
        return None
    best = None
    for name in ("zp_psf", "zp_auto"):
        try:
            zp = float(photometry.get(name))
            std = float(photometry.get(f"{name}_std"))
        except (TypeError, ValueError):
            continue
        if not (math.isfinite(zp) and math.isfinite(std)):
            continue
        if std > MAX_ZP_APER_STD:
            continue
        if best is None or std < best[1]:
            best = (zp, std, name)
    return best


def _galfit_mag_err_large(galfit_mag_err) -> bool:
    """True when GALFIT's formal mag error exceeds ``MAX_GALFIT_MAG_ERR``."""
    try:
        err = float(galfit_mag_err)
    except (TypeError, ValueError):
        return False
    return math.isfinite(err) and err > MAX_GALFIT_MAG_ERR


def final_host_mag(
    photometry: dict | None,
    reference: dict | None,
    galfit_mag,
    galfit_mag_err,
    zp_used: float | None = None,
) -> tuple[float | None, float | None, str]:
    """(mag_final, mag_final_err, source) applying the substitution policy.

    1. survey preferred when ZP untrusted **or** GALFIT ``mag_err > 1``
       → LS DR10 / PS1 r mag;
    2. pipeline aperture ZP trusted and ``mag_err ≤ 1`` → GALFIT mag as-is;
    3. a non-aperture pipeline ZP is tight → GALFIT mag rescaled by
       (zp_alt - zp_used), needs the J) value GALFIT ran with;
    4. no usable magnitude (no survey, no trusted/rescaled ZP)
       → ``(None, None, "unavailable")`` — do not ship an uncalibrated
       GALFIT mag for science.
    """
    ref_ok = (
        isinstance(reference, dict)
        and reference.get("status") == "ok"
        and reference.get("mag") is not None
    )
    zp_ok = pipeline_zp_ok(photometry)
    prefer_survey = (not zp_ok) or _galfit_mag_err_large(galfit_mag_err)

    if prefer_survey and ref_ok:
        return (
            reference["mag"],
            reference.get("mag_err"),
            f"reference_{reference.get('survey', 'unknown')}",
        )
    if zp_ok and galfit_mag is not None and not _galfit_mag_err_large(galfit_mag_err):
        return galfit_mag, galfit_mag_err, "galfit_calibrated_zp"
    # Survey preferred but missing: try rescaling onto a tight non-aperture ZP
    # only when the formal GALFIT mag_err is still usable (≤ 1 mag).
    if (
        galfit_mag is not None
        and zp_used is not None
        and not _galfit_mag_err_large(galfit_mag_err)
    ):
        alt = _best_alternative_zp(photometry)
        if alt is not None:
            zp_alt, std_alt, name = alt
            mag = round(float(galfit_mag) - zp_used + zp_alt, 4)
            err = None
            try:
                err = round(math.hypot(float(galfit_mag_err), std_alt), 4)
            except (TypeError, ValueError):
                err = round(std_alt, 4)
            return mag, err, f"galfit_rescaled_{name}"
    return None, None, "unavailable"


def attach_reference_photometry(summary: dict, folder: Path) -> None:
    """Embed reference_photometry.json into a pipeline summary dict and set
    ``galfit_host.mag_final`` / ``mag_final_err`` / ``mag_final_source``."""
    ref = None
    ref_path = folder / "reference_photometry.json"
    if ref_path.is_file():
        try:
            ref = json.loads(ref_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            ref = None
    summary["reference_photometry"] = ref

    galfit = summary.get("galfit_host")
    if isinstance(galfit, dict):
        mag_final, mag_final_err, source = final_host_mag(
            summary.get("photometry"), ref, galfit.get("mag"), galfit.get("mag_err"),
            zp_used=galfit_zp_used(folder),
        )
        galfit["mag_final"] = mag_final
        galfit["mag_final_err"] = mag_final_err
        galfit["mag_final_source"] = source


# --------------------------------------------------------------------------
# Backfill CLI for already-collected Output/<FRB>_<tag>/ folders
# --------------------------------------------------------------------------

def _backfill_folder(folder: Path, force: bool, radius_arcsec: float) -> str:
    ref_path = folder / "reference_photometry.json"
    if ref_path.is_file() and not force:
        try:
            existing = json.loads(ref_path.read_text(encoding="utf-8"))
            if existing.get("status") == "ok":
                return "skipped (already ok; use --force to requery)"
        except (json.JSONDecodeError, OSError):
            pass

    summary_path = folder / "pipeline_summary.json"
    summary = {}
    if summary_path.is_file():
        try:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, OSError):
            summary = {}

    ref = fetch_reference_photometry(
        folder,
        fallback_ra=summary.get("ra_deg"),
        fallback_dec=summary.get("dec_deg"),
        radius_arcsec=radius_arcsec,
    )
    if summary:
        attach_reference_photometry(summary, folder)
        summary_path.write_text(
            json.dumps(summary, indent=2, default=str), encoding="utf-8"
        )
    if ref.get("status") == "ok":
        err = ref.get("mag_err")
        err_s = f"+-{err:.3f}" if isinstance(err, (int, float)) else "+-?"
        return (
            f"{ref['survey']} r={ref['mag']:.3f}{err_s} "
            f"sep={ref['sep_arcsec']}\" ({ref.get('coord_source')})"
        )
    return ref.get("status", "unknown failure")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--backfill", action="store_true",
                        help="Query and store reference photometry for existing Output folders.")
    parser.add_argument("--output-root", type=Path,
                        default=Path(__file__).resolve().parent / "Output",
                        help="Root of <FRB>_<tag> folders (default: pipeline_scripts/Output).")
    parser.add_argument("--tag", default="all", help="Folder suffix (default: all).")
    parser.add_argument("--frb", nargs="+", default=None,
                        help="Only these FRB names (default: every folder).")
    parser.add_argument("--force", action="store_true",
                        help="Requery even when reference_photometry.json already exists.")
    parser.add_argument("--radius-arcsec", type=float, default=DEFAULT_SEARCH_RADIUS_ARCSEC,
                        help=f"Search radius (default: {DEFAULT_SEARCH_RADIUS_ARCSEC}).")
    args = parser.parse_args()

    if not args.backfill:
        parser.error("nothing to do: pass --backfill (module is otherwise used as a library)")

    suffix = f"_{args.tag}"
    folders = sorted(
        d for d in args.output_root.iterdir()
        if d.is_dir() and d.name.endswith(suffix)
    )
    if args.frb:
        wanted = set(args.frb)
        folders = [d for d in folders if d.name[: -len(suffix)] in wanted]
    if not folders:
        print(f"[!] No matching folders under {args.output_root}")
        return

    n_ok = 0
    for folder in folders:
        frb = folder.name[: -len(suffix)]
        outcome = _backfill_folder(folder, args.force, args.radius_arcsec)
        if not outcome.startswith(("skipped", "no_match", "error")):
            n_ok += 1
        print(f"  {frb}: {outcome}")
    print(f"[*] Backfill done: {n_ok}/{len(folders)} with reference mag.")


if __name__ == "__main__":
    main()
