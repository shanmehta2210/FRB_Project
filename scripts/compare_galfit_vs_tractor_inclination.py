import math
import time
from typing import Dict, Optional

import numpy as np
import pandas as pd
import pyvo

# Per project note: these failed Legacy Survey coverage and used Pan-STARRS fallback.
FAILED_LEGACY_FRBS = {"20171020A", "20210807D", "20211127I"}
Q0 = 0.2


def incl_from_q(q: float, q0: float = Q0) -> float:
    if not np.isfinite(q):
        return np.nan
    if q <= q0:
        return 90.0
    val = (q * q - q0 * q0) / (1.0 - q0 * q0)
    val = min(1.0, max(0.0, val))
    return math.degrees(math.acos(math.sqrt(val)))


def q_from_e1e2(e1: float, e2: float) -> tuple[float, float]:
    # Tractor ellipticity convention: |e| = (1 - q) / (1 + q)
    eabs = float(np.hypot(e1, e2))
    if eabs >= 1.0:
        return np.nan, eabs
    q = (1.0 - eabs) / (1.0 + eabs)
    return q, eabs


def sigma_from_ivar(ivar: float) -> float:
    if not np.isfinite(ivar) or ivar <= 0:
        return np.nan
    return float(1.0 / np.sqrt(ivar))


def incl_from_q_array(q: np.ndarray, q0: float = Q0) -> np.ndarray:
    q = np.asarray(q, dtype=float)
    val = (q * q - q0 * q0) / (1.0 - q0 * q0)
    val = np.clip(val, 0.0, 1.0)
    return np.degrees(np.arccos(np.sqrt(val)))


def incl_err_from_shape_ivar_mc(
    e1: float,
    e2: float,
    e1_ivar: float,
    e2_ivar: float,
    n_draws: int = 4000,
) -> float:
    s1 = sigma_from_ivar(e1_ivar)
    s2 = sigma_from_ivar(e2_ivar)
    if not np.isfinite(s1) or not np.isfinite(s2):
        return np.nan

    rng = np.random.default_rng(12345)
    e1_draw = rng.normal(float(e1), s1, int(n_draws))
    e2_draw = rng.normal(float(e2), s2, int(n_draws))
    eabs = np.hypot(e1_draw, e2_draw)

    valid = eabs < 1.0
    if np.count_nonzero(valid) < 50:
        return np.nan

    q = (1.0 - eabs[valid]) / (1.0 + eabs[valid])
    i_draw = incl_from_q_array(q)
    return float(np.nanstd(i_draw, ddof=1))


def _ra_clause(ra: float, dra: float) -> str:
    ra_min = ra - dra
    ra_max = ra + dra
    if ra_min < 0:
        return f"(ra > {ra_min + 360:.8f} OR ra < {ra_max:.8f})"
    if ra_max > 360:
        return f"(ra > {ra_min:.8f} OR ra < {ra_max - 360:.8f})"
    return f"ra > {ra_min:.8f} AND ra < {ra_max:.8f}"


def query_nearest_tractor(
    svc: pyvo.dal.TAPService,
    ra: float,
    dec: float,
    table: str = "ls_dr10.tractor",
    radius_arcsec: float = 10.0,
) -> Optional[Dict]:
    dec_clip = max(-85.0, min(85.0, dec))
    dra = (radius_arcsec / 3600.0) / math.cos(math.radians(dec_clip))
    ddec = radius_arcsec / 3600.0

    query = f"""
    SELECT TOP 300 objid, ra, dec, type, flux_r, sersic, shape_e1, shape_e2, shape_e1_ivar, shape_e2_ivar
    FROM {table}
    WHERE {_ra_clause(ra, dra)}
      AND dec > {dec - ddec:.8f} AND dec < {dec + ddec:.8f}
      AND flux_r > 0
    """

    tab = svc.search(query).to_table()
    if len(tab) == 0:
        return None

    ra_arr = np.array(tab["ra"], dtype=float)
    dec_arr = np.array(tab["dec"], dtype=float)
    dra_as = (ra_arr - ra) * np.cos(np.radians(dec)) * 3600.0
    ddec_as = (dec_arr - dec) * 3600.0
    sep = np.hypot(dra_as, ddec_as)

    types = np.array(tab["type"]).astype(str)
    non_psf_idx = np.where(types != "PSF")[0]
    if len(non_psf_idx) > 0:
        idx = int(non_psf_idx[np.argmin(sep[non_psf_idx])])
    else:
        idx = int(np.argmin(sep))

    r = tab[idx]
    return {
        "objid": int(r["objid"]),
        "ra_ls": float(r["ra"]),
        "dec_ls": float(r["dec"]),
        "type_ls": str(r["type"]),
        "flux_r_ls": float(r["flux_r"]),
        "sersic_ls_fit": float(r["sersic"]),
        "shape_e1": float(r["shape_e1"]),
        "shape_e2": float(r["shape_e2"]),
        "shape_e1_ivar": float(r["shape_e1_ivar"]),
        "shape_e2_ivar": float(r["shape_e2_ivar"]),
        "sep_arcsec": float(sep[idx]),
        "n_candidates": int(len(tab)),
    }


def main() -> None:
    frb = pd.read_csv("master_frb_summary.csv")
    frb = frb[~frb["FRB"].isin(FAILED_LEGACY_FRBS)].copy()

    svc = pyvo.dal.TAPService("https://datalab.noirlab.edu/tap")

    rows = []
    for i, row in enumerate(frb.itertuples(index=False), start=1):
        rec = {
            "FRB": row.FRB,
            "RA_deg": float(row.RA_deg),
            "DEC_deg": float(row.DEC_deg),
            "galfit_inc_psf_deg": float(row.inc_psf),
        }

        print(f"[{i}/{len(frb)}] Querying {row.FRB}...")
        match = None
        last_err = None
        for _ in range(3):
            try:
                match = query_nearest_tractor(svc, float(row.RA_deg), float(row.DEC_deg))
                break
            except Exception as exc:
                last_err = str(exc)
                time.sleep(1.5)

        if match is None:
            rec["match_found"] = False
            if last_err:
                rec["query_error"] = last_err
            rows.append(rec)
            continue

        rec["match_found"] = True
        rec.update(match)

        q_ls, eabs = q_from_e1e2(rec["shape_e1"], rec["shape_e2"])
        rec["ellipticity_abs"] = eabs
        rec["q_ls_from_e"] = q_ls
        rec["shape_e1_sigma"] = sigma_from_ivar(rec["shape_e1_ivar"])
        rec["shape_e2_sigma"] = sigma_from_ivar(rec["shape_e2_ivar"])
        rec["ls_inc_deg"] = incl_from_q(q_ls)
        rec["ls_inc_err_deg"] = incl_err_from_shape_ivar_mc(
            rec["shape_e1"],
            rec["shape_e2"],
            rec["shape_e1_ivar"],
            rec["shape_e2_ivar"],
        )
        if np.isfinite(rec["ls_inc_deg"]):
            rec["delta_deg_ls_minus_galfit"] = rec["ls_inc_deg"] - rec["galfit_inc_psf_deg"]
        else:
            rec["delta_deg_ls_minus_galfit"] = np.nan

        rows.append(rec)

    out = pd.DataFrame(rows)
    out_path = "legacy_vs_galfit_inclination_comparison.csv"
    out.to_csv(out_path, index=False)

    matched = out[out["match_found"] == True].copy()
    valid = matched[np.isfinite(matched["ls_inc_deg"])].copy()

    print("\nSaved:", out_path)
    print("Total FRBs in master_frb_summary:", len(pd.read_csv("master_frb_summary.csv")))
    print("Legacy-eligible FRBs:", len(frb))
    print("Matched rows:", len(matched))
    print("Valid LS inclinations:", len(valid))

    if len(valid) > 0:
        d = valid["delta_deg_ls_minus_galfit"].astype(float)
        print("Mean delta (LS-GALFIT) deg:", round(float(d.mean()), 3))
        print("Median delta deg:", round(float(d.median()), 3))
        print("RMS delta deg:", round(float(np.sqrt(np.mean(d**2))), 3))
        print("MAE delta deg:", round(float(np.mean(np.abs(d))), 3))


if __name__ == "__main__":
    main()
