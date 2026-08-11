"""
Build Legacy Survey DR10 Tractor null catalog (v2) — full LS footprint, EXP-only.

Unlike v1 (``build_legacy_catalog_csv.py``), v2:
  - Covers the full Legacy footprint (no joint-SDSS Dec clip).
  - Keeps only Tractor ``type = 'EXP'`` (true exponential disks; not REX/DEV/SER).
  - Targets ~2M unique galaxies with SDSS-v2-parallel morphology columns.

Random sampling (the correct, simple way):
  ``ls_dr10.tractor`` ships a precomputed, indexed ``random_id`` column that is
  uniform on [0, 100] and *independent of sky position and galaxy shape*. Selecting
  ``random_id`` in a window therefore gives a uniform random sample spread across the
  whole footprint with no positional or shape clustering — no sky tiling required.

  We pull the sample in successive ``random_id`` windows via **async UWS jobs**
  (each window ~50-60k rows; the ``/tap/sync`` endpoint has a hard ~60s limit that
  a full 2M — or even 236k — pull would blow past). Each window is cached, so the
  build is fully resumable.

Density (measured 2026-07): ~23.6k EXP per 0.01 of ``random_id`` (~2.36M per unit),
so ``random_id < ~0.85`` yields ~2M galaxies.

Run from repo root::

    python scripts/build_legacy_catalog_v2_exp.py
"""

from __future__ import annotations

import argparse
import concurrent.futures
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyvo

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (
    LS_CATALOG_V2_EXP_DEFAULT,
    LS_V2_EXP_TARGET_ROWS,
    Q0,
    footprint_summary,
)

TAP_URL = "https://datalab.noirlab.edu/tap"
TABLE = "ls_dr10.tractor"
DEFAULT_OUT = LS_CATALOG_V2_EXP_DEFAULT
DEFAULT_CACHE_DIR = "LS_v2_exp_fetch_cache"

# random_id sampling window. ~23.6k EXP per 0.01; 0.015 -> ~35k rows/window.
# Kept lean because the Data Lab async service is flaky under large results.
DEFAULT_RID_WINDOW = 0.015
# Per-sync-call slice: ~0.004 -> ~9.5k rows, sized to (usually) clear the ~60s cap.
DEFAULT_SYNC_SUBWINDOW = 0.004
DEFAULT_RID_START = 0.0
DEFAULT_RID_MAX = 100.0
DEFAULT_CELL_CAP = 200_000  # per-window safety TOP; keep >> rows/window
DEFAULT_RETRIES = 6
DEFAULT_POLL_SEC = 8.0
DEFAULT_JOB_TIMEOUT_SEC = 900.0
# Hard client-side cap on a single sync TAP call. Without this, a hung socket can
# stall the whole build for hours with no log output (seen 2026-07-11).
DEFAULT_SYNC_CALL_TIMEOUT_SEC = 120.0

# Lean column set: every SDSS-v2 morphology parallel we actually use, nothing more
# (extra WISE/Gaia/per-band columns bloat the result and destabilise the service).
#   position: ra, dec              magnitudes: flux_{griz} (-> mags+colors), dered_mag_r
#   shape/incl: shape_e1/e2(+ivar) size: shape_r   profile: sersic, type
#   morphology evidence: dchisq_1..5 (PSF,REX,DEV,EXP,SER ~ SDSS lnL)
#   QA: ebv, snr_r, nobs_r, maskbits   ids/sampling: ls_id, brickid, release, random_id
_SELECT_COLS = """
    ls_id, objid, brickid, release, ra, dec, type, brick_primary,
    flux_g, flux_r, flux_i, flux_z, flux_ivar_r,
    dered_mag_r,
    sersic, shape_r,
    shape_e1, shape_e2, shape_e1_ivar, shape_e2_ivar,
    dchisq_1, dchisq_2, dchisq_3, dchisq_4, dchisq_5,
    ebv, snr_r, nobs_r, maskbits, random_id
"""


def flux_to_mag(flux_nmgy: np.ndarray) -> np.ndarray:
    flux = np.asarray(flux_nmgy, dtype=float)
    out = np.full(flux.shape, np.nan, dtype=float)
    good = np.isfinite(flux) & (flux > 0)
    out[good] = 22.5 - 2.5 * np.log10(flux[good])
    return out


def sigma_from_ivar(ivar: np.ndarray) -> np.ndarray:
    iv = np.asarray(ivar, dtype=float)
    out = np.full(iv.shape, np.nan, dtype=float)
    good = np.isfinite(iv) & (iv > 0)
    out[good] = 1.0 / np.sqrt(iv[good])
    return out


def q_from_e1e2(e1: np.ndarray, e2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    e1a = np.asarray(e1, dtype=float)
    e2a = np.asarray(e2, dtype=float)
    eabs = np.hypot(e1a, e2a)
    q = np.full(eabs.shape, np.nan, dtype=float)
    good = np.isfinite(eabs) & (eabs < 1.0)
    q[good] = (1.0 - eabs[good]) / (1.0 + eabs[good])
    return q, eabs


def q_err_from_e1e2_ivar(
    e1: np.ndarray,
    e2: np.ndarray,
    e1_ivar: np.ndarray,
    e2_ivar: np.ndarray,
) -> np.ndarray:
    e1a = np.asarray(e1, dtype=float)
    e2a = np.asarray(e2, dtype=float)
    s1 = sigma_from_ivar(e1_ivar)
    s2 = sigma_from_ivar(e2_ivar)
    eabs = np.hypot(e1a, e2a)
    out = np.full(eabs.shape, np.nan, dtype=float)
    valid = np.isfinite(eabs) & (eabs < 1.0) & np.isfinite(s1) & np.isfinite(s2)
    if not np.any(valid):
        return out
    idx = np.where(valid)[0]
    ev = eabs[idx]
    sigma_e = np.zeros_like(ev)
    nonzero = ev > 1e-12
    if np.any(nonzero):
        j = idx[nonzero]
        ev_nz = eabs[j]
        term1 = (e1a[j] / ev_nz) ** 2 * (s1[j] ** 2)
        term2 = (e2a[j] / ev_nz) ** 2 * (s2[j] ** 2)
        sigma_e[nonzero] = np.sqrt(term1 + term2)
    if np.any(~nonzero):
        j0 = idx[~nonzero]
        sigma_e[~nonzero] = np.sqrt(0.5 * (s1[j0] ** 2 + s2[j0] ** 2))
    dq_de = 2.0 / (1.0 + ev) ** 2
    out[idx] = dq_de * sigma_e
    return out


def incl_from_q(q: np.ndarray, q0: float = Q0) -> np.ndarray:
    qa = np.asarray(q, dtype=float)
    out = np.full(qa.shape, np.nan, dtype=float)
    good = np.isfinite(qa)
    if not np.any(good):
        return out
    out[good & (qa <= q0)] = 90.0
    hi = good & (qa > q0)
    val = (qa[hi] ** 2 - q0**2) / (1.0 - q0**2)
    val = np.clip(val, 0.0, 1.0)
    out[hi] = np.degrees(np.arccos(np.sqrt(val)))
    return out


def _sql_rid_window(rid_lo: float, rid_hi: float, top_n: int) -> str:
    """EXP query for one random_id window (uniform random sample, index-backed)."""
    return f"""
    SELECT TOP {int(top_n)}
    {_SELECT_COLS.strip()}
    FROM {TABLE}
    WHERE brick_primary = 1
      AND type = 'EXP'
      AND flux_r > 0
      AND shape_e1 IS NOT NULL
      AND shape_e2 IS NOT NULL
      AND shape_e1_ivar IS NOT NULL
      AND shape_e2_ivar IS NOT NULL
      AND shape_e1 > -0.999999 AND shape_e1 < 0.999999
      AND shape_e2 > -0.999999 AND shape_e2 < 0.999999
      AND random_id >= {rid_lo:.10g} AND random_id < {rid_hi:.10g}
    """


def _poll_phase(job: pyvo.dal.AsyncTAPJob) -> str | None:
    """Return the job phase, or ``None`` if the phase check failed transiently."""
    try:
        return str(job.phase).upper()
    except Exception:  # noqa: BLE001 - Data Lab 502/503/504 while polling is transient
        return None


def fetch_rid_window_async(
    service: pyvo.dal.TAPService,
    rid_lo: float,
    rid_hi: float,
    top_n: int,
    *,
    retries: int,
    poll_sec: float,
    job_timeout_sec: float,
) -> pd.DataFrame:
    """
    Fetch one random_id window via an async UWS job (avoids the ~60s sync cap).

    Data Lab is intermittently flaky (502/503/504). This is hardened so that:
      - transient errors while *polling* the job phase do not abandon the (possibly
        still-running) job — we keep polling until ``job_timeout_sec``;
      - the result fetch is itself retried;
      - the whole submit->run->fetch cycle is retried ``retries`` times.
    """
    sql = _sql_rid_window(rid_lo, rid_hi, top_n)
    last_err: Exception | None = None

    for attempt in range(1, retries + 1):
        job = None
        try:
            job = pyvo.dal.AsyncTAPJob.create(service.baseurl, sql, session=service._session)
            job.run()
            t0 = time.time()
            phase = None
            consecutive_poll_errors = 0
            while True:
                phase = _poll_phase(job)
                if phase in {"COMPLETED", "ERROR", "ABORTED"}:
                    break
                if phase is None:
                    consecutive_poll_errors += 1
                    # Bail this attempt only if polling has been failing for a long time
                    if consecutive_poll_errors * poll_sec > job_timeout_sec:
                        raise RuntimeError("phase polling failed repeatedly")
                else:
                    consecutive_poll_errors = 0
                if time.time() - t0 > job_timeout_sec:
                    raise TimeoutError(f"job exceeded {job_timeout_sec:.0f}s (phase={phase})")
                time.sleep(poll_sec)

            if phase != "COMPLETED":
                raise RuntimeError(f"job phase={phase}")

            # Result fetch can also 5xx transiently; retry a few times.
            fetch_err: Exception | None = None
            for fetch_try in range(1, 5):
                try:
                    return job.fetch_result().to_table().to_pandas()
                except Exception as fexc:  # noqa: BLE001
                    fetch_err = fexc
                    time.sleep(min(20.0, 5.0 * fetch_try))
            raise RuntimeError(f"result fetch failed: {fetch_err}")
        except Exception as exc:  # noqa: BLE001 - transient TAP failures are retried
            last_err = exc
            print(
                f"    [retry] window [{rid_lo:.4g},{rid_hi:.4g}) attempt {attempt}/{retries}: "
                f"{type(exc).__name__} {str(exc)[:120]}",
                flush=True,
            )
            if attempt < retries:
                time.sleep(min(60.0, 8.0 * attempt))
        finally:
            if job is not None:
                try:
                    job.delete()
                except Exception:  # noqa: BLE001
                    pass

    raise RuntimeError(
        f"random_id window [{rid_lo:.5g}, {rid_hi:.5g}) failed after {retries} attempts: {last_err}"
    )


def fetch_rid_window_sync(
    service: pyvo.dal.TAPService,
    rid_lo: float,
    rid_hi: float,
    top_n: int,
    *,
    retries: int,
    call_timeout_sec: float = DEFAULT_SYNC_CALL_TIMEOUT_SEC,
) -> pd.DataFrame:
    """
    Fetch one (small) random_id window via the sync endpoint.

    Used when the async/UWS subsystem is unavailable. The sync gateway has a soft
    ~60s cap (it sometimes 504s, sometimes lets a query run longer), so callers must
    keep the window small (a few thousand rows). Transient 5xx/timeouts are retried.
    A hard client-side ``call_timeout_sec`` aborts hung sockets so the build cannot
    stall silently for hours.
    """
    sql = _sql_rid_window(rid_lo, rid_hi, top_n)
    last_err: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
                fut = pool.submit(lambda: service.search(sql).to_table().to_pandas())
                return fut.result(timeout=call_timeout_sec)
        except Exception as exc:  # noqa: BLE001 - transient TAP failures are retried
            last_err = exc
            print(
                f"    [retry-sync] window [{rid_lo:.5g},{rid_hi:.5g}) attempt {attempt}/{retries}: "
                f"{type(exc).__name__} {str(exc)[:100]}",
                flush=True,
            )
            if attempt < retries:
                time.sleep(min(45.0, 6.0 * attempt))
    raise RuntimeError(
        f"random_id window [{rid_lo:.5g}, {rid_hi:.5g}) sync failed after {retries} attempts: {last_err}"
    )


def fetch_window(
    service: pyvo.dal.TAPService,
    rid_lo: float,
    rid_hi: float,
    top_n: int,
    *,
    mode: str,
    sync_subwindow: float,
    retries: int,
    poll_sec: float,
    job_timeout_sec: float,
) -> pd.DataFrame:
    """
    Fetch a random_id window using the requested transport.

    mode='async' : one async UWS job (efficient; needs a healthy async subsystem).
    mode='sync'  : split the window into <=``sync_subwindow`` slices, each a sync call
                   sized to fit under the ~60s gateway cap.
    mode='auto'  : try async once; on failure fall back to sync slices for the window.
    """
    if mode == "async":
        return fetch_rid_window_async(
            service, rid_lo, rid_hi, top_n,
            retries=retries, poll_sec=poll_sec, job_timeout_sec=job_timeout_sec,
        )

    if mode == "auto":
        try:
            return fetch_rid_window_async(
                service, rid_lo, rid_hi, top_n,
                retries=2, poll_sec=poll_sec, job_timeout_sec=job_timeout_sec,
            )
        except RuntimeError as exc:
            print(f"    [auto] async failed, falling back to sync: {str(exc)[:80]}", flush=True)

    # sync (or auto-fallback): split into small slices under the 60s cap
    parts: list[pd.DataFrame] = []
    lo = rid_lo
    while lo < rid_hi:
        hi = min(lo + sync_subwindow, rid_hi)
        parts.append(
            fetch_rid_window_sync(
                service, lo, hi, top_n,
                retries=retries,
                call_timeout_sec=DEFAULT_SYNC_CALL_TIMEOUT_SEC,
            )
        )
        lo = hi
    if not parts:
        return pd.DataFrame()
    return pd.concat(parts, ignore_index=True)


def build_catalog(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()

    def col(name: str) -> pd.Series:
        """Numeric column, or an all-NaN series if it was not selected."""
        if name in df.columns:
            return pd.to_numeric(df[name], errors="coerce")
        return pd.Series(np.nan, index=df.index, dtype=float)

    e1 = pd.to_numeric(df["shape_e1"], errors="coerce").to_numpy(dtype=float)
    e2 = pd.to_numeric(df["shape_e2"], errors="coerce").to_numpy(dtype=float)
    e1_ivar = pd.to_numeric(df["shape_e1_ivar"], errors="coerce").to_numpy(dtype=float)
    e2_ivar = pd.to_numeric(df["shape_e2_ivar"], errors="coerce").to_numpy(dtype=float)

    q, eabs = q_from_e1e2(e1, e2)
    q_err = q_err_from_e1e2_ivar(e1, e2, e1_ivar, e2_ivar)
    inc = incl_from_q(q)

    # Flux-derived mags (matches v1 convention: 22.5 - 2.5 log10 flux, no extinction).
    gmag = flux_to_mag(pd.to_numeric(df["flux_g"], errors="coerce").to_numpy(dtype=float))
    rmag = flux_to_mag(pd.to_numeric(df["flux_r"], errors="coerce").to_numpy(dtype=float))
    imag = flux_to_mag(pd.to_numeric(df["flux_i"], errors="coerce").to_numpy(dtype=float))
    zmag = flux_to_mag(pd.to_numeric(df["flux_z"], errors="coerce").to_numpy(dtype=float))

    shape_r = pd.to_numeric(df["shape_r"], errors="coerce")
    sersic = pd.to_numeric(df["sersic"], errors="coerce")

    ls_id = pd.to_numeric(df["ls_id"], errors="coerce")

    out = pd.DataFrame(
        {
            "ls_id": ls_id.astype("Int64"),
            "release": pd.to_numeric(df["release"], errors="coerce"),
            "brickid": pd.to_numeric(df["brickid"], errors="coerce"),
            "objid": pd.to_numeric(df["objid"], errors="coerce"),
            "RA_ICRS": pd.to_numeric(df["ra"], errors="coerce"),
            "DE_ICRS": pd.to_numeric(df["dec"], errors="coerce"),
            # Magnitude parallels to SDSS modelMag_* (flux-derived, no extinction)
            "tractor_mag_r": rmag,
            "rmag": rmag,
            "modelMag_r": rmag,
            "gmag": gmag,
            "modelMag_g": gmag,
            "imag": imag,
            "zmag": zmag,
            "g_r": gmag - rmag,
            "r_i": rmag - imag,
            "i_z": imag - zmag,
            # Table-provided extinction-corrected r (SDSS modelMag_r is not dered;
            # kept for a dust-sensitivity check).
            "dered_mag_r": col("dered_mag_r"),
            # Shape / inclination (EXP profile by construction)
            "expAB_r": q,
            "b_a": q,
            "b_a_err": q_err,
            "best_model_ba_r": q,
            "q_lt_q0": q <= Q0,
            "inclination_deg_q0_0p2": inc,
            "ellipticity_abs": eabs,
            "shape_e1": e1,
            "shape_e2": e2,
            "shape_e1_ivar": e1_ivar,
            "shape_e2_ivar": e2_ivar,
            # Size / Sérsic parallels to expRad_r / n_eff_r
            "shape_r": shape_r,
            "expRad_r": shape_r,
            "best_model_re_r": shape_r,
            "rPrad": shape_r,
            "sersic": sersic,
            "n_eff_r": sersic,
            "rdVrad": sersic,
            # Morphology: type='EXP' means the exponential model won over
            # PSF/REX/DEV/SER. dchisq_{1..5} are the per-model chi2 improvements
            # (PSF, REX, DEV, EXP, SER) — Tractor's analogue of SDSS lnL evidence.
            "tractor_type": df["type"].astype(str),
            "model_winner_is_exp": True,
            "dchisq_psf": col("dchisq_1"),
            "dchisq_rex": col("dchisq_2"),
            "dchisq_dev": col("dchisq_3"),
            "dchisq_exp": col("dchisq_4"),
            "dchisq_ser": col("dchisq_5"),
            "brick_primary": col("brick_primary"),
            "ebv": col("ebv"),
            "snr_r": col("snr_r"),
            "nobs_r": col("nobs_r"),
            "maskbits": col("maskbits"),
            "flux_g_nmgy": col("flux_g"),
            "flux_r_nmgy": col("flux_r"),
            "flux_i_nmgy": col("flux_i"),
            "flux_z_nmgy": col("flux_z"),
            "flux_ivar_r": col("flux_ivar_r"),
            "random_id": col("random_id"),
        }
    )

    good = (
        out["ls_id"].notna()
        & np.isfinite(out["RA_ICRS"])
        & np.isfinite(out["DE_ICRS"])
        & np.isfinite(out["tractor_mag_r"])
        & np.isfinite(out["expAB_r"])
        & (out["expAB_r"] >= 0.0)
        & (out["expAB_r"] <= 1.0)
    )
    out = out.loc[good].copy().reset_index(drop=True)
    return out


def query_random(
    *,
    target_rows: int,
    rid_window: float,
    rid_start: float,
    rid_max: float,
    cell_cap: int,
    mode: str,
    sync_subwindow: float,
    retries: int,
    poll_sec: float,
    job_timeout_sec: float,
    cache_dir: Path | None,
    resume: bool,
) -> tuple[pd.DataFrame, list[dict]]:
    """
    Fill ``target_rows`` EXP galaxies from successive ``random_id`` windows.

    Because ``random_id`` is uniform and independent of position/shape, the sample
    is automatically spread across the whole sky with no clustering.
    """
    service = pyvo.dal.TAPService(TAP_URL)
    batch_log: list[dict] = []
    frames: list[pd.DataFrame] = []
    seen: set[int] = set()
    n_unique = 0

    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)

    n_windows = int(np.ceil((rid_max - rid_start) / rid_window))
    print(
        f"[*] Full-sky EXP sample via random_id: target={target_rows:,}, "
        f"mode={mode}, window={rid_window}"
        + (f" (sync sub={sync_subwindow})" if mode != "async" else "")
        + f", start={rid_start}, up to {n_windows} windows",
        flush=True,
    )

    widx = 0
    rid_lo = rid_start
    while n_unique < target_rows and rid_lo < rid_max:
        rid_hi = min(rid_lo + rid_window, rid_max)
        cache_path = cache_dir / f"rid_{widx:04d}.csv" if cache_dir else None
        raw: pd.DataFrame | None = None
        source = "query"

        if resume and cache_path is not None and cache_path.is_file():
            raw = pd.read_csv(cache_path)
            source = "cache"
        else:
            try:
                raw = fetch_window(
                    service,
                    rid_lo,
                    rid_hi,
                    cell_cap,
                    mode=mode,
                    sync_subwindow=sync_subwindow,
                    retries=retries,
                    poll_sec=poll_sec,
                    job_timeout_sec=job_timeout_sec,
                )
            except RuntimeError as exc:
                print(f"  [!] window {widx} [{rid_lo:.4g},{rid_hi:.4g}): {exc}", flush=True)
                batch_log.append(
                    {
                        "window": widx,
                        "rid_lo": rid_lo,
                        "rid_hi": rid_hi,
                        "n_raw": 0,
                        "n_kept": 0,
                        "status": "failed",
                    }
                )
                widx += 1
                rid_lo = rid_hi
                continue
            if cache_path is not None:
                raw.to_csv(cache_path, index=False)

        n_raw = 0 if raw is None else len(raw)
        if n_raw >= cell_cap:
            print(
                f"  [warn] window {widx} hit cap {cell_cap:,}; shrink --rid-window",
                flush=True,
            )

        kept = 0
        if n_raw:
            built = build_catalog(raw)
            if not built.empty:
                ids = built["ls_id"].astype("int64")
                mask = ~ids.isin(seen)
                new = built.loc[mask.to_numpy()]
                if len(new):
                    frames.append(new)
                    seen.update(new["ls_id"].astype("int64").tolist())
                    n_unique += len(new)
                    kept = len(new)

        batch_log.append(
            {
                "window": widx,
                "rid_lo": rid_lo,
                "rid_hi": rid_hi,
                "n_raw": n_raw,
                "n_kept": kept,
                "status": source,
            }
        )
        print(
            f"  window {widx} rid[{rid_lo:.4g},{rid_hi:.4g}) [{source}] "
            f"raw={n_raw} kept={kept} -> unique={n_unique:,}/{target_rows:,}",
            flush=True,
        )

        widx += 1
        rid_lo = rid_hi

    if not frames:
        raise RuntimeError("No rows returned from any random_id window.")

    catalog = pd.concat(frames, ignore_index=True)
    catalog = catalog.drop_duplicates(subset=["ls_id"], keep="first").reset_index(drop=True)
    if len(catalog) > target_rows:
        catalog = catalog.sample(n=target_rows, random_state=42).reset_index(drop=True)

    return catalog, batch_log


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Build Legacy DR10 full-sky EXP-only null catalog "
            "(v2, ~2M, uniform random_id sampling)."
        )
    )
    parser.add_argument(
        "--target-rows",
        type=int,
        default=LS_V2_EXP_TARGET_ROWS,
        help="Target unique EXP galaxies (default 2_000_000).",
    )
    parser.add_argument(
        "--rid-window",
        type=float,
        default=DEFAULT_RID_WINDOW,
        help="random_id window width per async job (~23.6k EXP per 0.01).",
    )
    parser.add_argument("--rid-start", type=float, default=DEFAULT_RID_START)
    parser.add_argument("--rid-max", type=float, default=DEFAULT_RID_MAX)
    parser.add_argument(
        "--mode",
        choices=["async", "sync", "auto"],
        default="auto",
        help=(
            "Transport: 'async' (efficient, needs healthy UWS), "
            "'sync' (small slices under the 60s cap), "
            "'auto' (try async, fall back to sync). Default: auto."
        ),
    )
    parser.add_argument(
        "--sync-subwindow",
        type=float,
        default=DEFAULT_SYNC_SUBWINDOW,
        help="random_id slice width per sync call (~23.6k EXP per 0.01; keep small).",
    )
    parser.add_argument(
        "--cell-cap",
        type=int,
        default=DEFAULT_CELL_CAP,
        help="Safety TOP per window (keep well above rows/window).",
    )
    parser.add_argument("--retries", type=int, default=DEFAULT_RETRIES)
    parser.add_argument("--poll-sec", type=float, default=DEFAULT_POLL_SEC)
    parser.add_argument("--job-timeout-sec", type=float, default=DEFAULT_JOB_TIMEOUT_SEC)
    parser.add_argument("--out", default=DEFAULT_OUT, help="Output CSV path.")
    parser.add_argument(
        "--cache-dir",
        default=DEFAULT_CACHE_DIR,
        help="Per-window raw TAP cache directory.",
    )
    parser.add_argument("--no-cache", action="store_true", help="Do not read/write cache.")
    parser.add_argument(
        "--no-resume", action="store_true", help="Ignore existing cache; re-query."
    )
    args = parser.parse_args()

    if args.target_rows <= 0:
        raise ValueError("--target-rows must be > 0")
    if args.rid_window <= 0:
        raise ValueError("--rid-window must be > 0")

    cache_dir = None if args.no_cache else Path(args.cache_dir)
    catalog, batch_log = query_random(
        target_rows=args.target_rows,
        rid_window=args.rid_window,
        rid_start=args.rid_start,
        rid_max=args.rid_max,
        cell_cap=args.cell_cap,
        mode=args.mode,
        sync_subwindow=args.sync_subwindow,
        retries=args.retries,
        poll_sec=args.poll_sec,
        job_timeout_sec=args.job_timeout_sec,
        cache_dir=cache_dir,
        resume=not args.no_resume,
    )

    footprint_summary(catalog, "Legacy v2 EXP full-sky")
    print(
        f"  expAB_r median={catalog['expAB_r'].median():.4f} "
        f"mean={catalog['expAB_r'].mean():.4f}",
        flush=True,
    )
    print(
        f"  RA span={catalog['RA_ICRS'].min():.2f}..{catalog['RA_ICRS'].max():.2f}  "
        f"Dec={catalog['DE_ICRS'].min():.2f}..{catalog['DE_ICRS'].max():.2f}",
        flush=True,
    )
    print(f"  types: {catalog['tractor_type'].value_counts().to_dict()}", flush=True)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    catalog.to_csv(out_path, index=False)
    log_path = out_path.with_suffix(".batch_log.csv")
    pd.DataFrame(batch_log).to_csv(log_path, index=False)
    print(f"Wrote: {out_path} ({len(catalog):,} rows)", flush=True)
    print(f"Wrote: {log_path}", flush=True)


if __name__ == "__main__":
    main()
