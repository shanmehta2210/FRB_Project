import math
import time
from typing import Dict, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyvo
from scipy.stats import pearsonr, spearmanr


FAILED_LEGACY_FRBS = {"20171020A", "20210807D", "20211127I"}
PIX_SCALE_ARCSEC = 0.262


def parse_float(value):
    if value is None:
        return np.nan
    try:
        text = str(value).replace("*", "").strip()
        return float(text)
    except Exception:
        return np.nan


def _ra_clause(ra: float, dra: float) -> str:
    ra_min = ra - dra
    ra_max = ra + dra
    if ra_min < 0:
        return f"(ra > {ra_min + 360:.8f} OR ra < {ra_max:.8f})"
    if ra_max > 360:
        return f"(ra > {ra_min:.8f} OR ra < {ra_max - 360:.8f})"
    return f"ra > {ra_min:.8f} AND ra < {ra_max:.8f}"


def query_nearest_tractor_reff(
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
    SELECT TOP 300 objid, ra, dec, type, flux_r, sersic, shape_r, shape_e1, shape_e2
    FROM {table}
    WHERE {_ra_clause(ra, dra)}
      AND dec > {dec - ddec:.8f} AND dec < {dec + ddec:.8f}
      AND flux_r > 0
      AND shape_r IS NOT NULL
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
        "tractor_objid": int(r["objid"]),
        "type_ls": str(r["type"]),
        "ra_ls": float(r["ra"]),
        "dec_ls": float(r["dec"]),
        "sep_arcsec": float(sep[idx]),
        "shape_r_ls_arcsec": float(r["shape_r"]),
        "sersic_n_ls_fit": float(r["sersic"]),
        "shape_e1": float(r["shape_e1"]),
        "shape_e2": float(r["shape_e2"]),
        "n_candidates": int(len(tab)),
    }


def make_plot(df: pd.DataFrame, out_png: str, out_pdf: str) -> None:
    fit_df = df[
        np.isfinite(df["re_psf_arcsec"]) & np.isfinite(df["shape_r_ls_arcsec"])
    ].copy()

    x = fit_df["re_psf_arcsec"].to_numpy(dtype=float)
    y = fit_df["shape_r_ls_arcsec"].to_numpy(dtype=float)

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(x, y, s=36, alpha=0.9)

    if len(x) >= 2:
        lo = min(np.min(x), np.min(y))
        hi = max(np.max(x), np.max(y))
        ax.plot([lo, hi], [lo, hi], "k--", linewidth=1.2, label="1:1")

        m, b = np.polyfit(x, y, 1)
        xx = np.linspace(lo, hi, 100)
        ax.plot(xx, m * xx + b, color="tab:red", linewidth=1.4, label=f"fit: y={m:.2f}x+{b:.2f}")

        r_p, p_p = pearsonr(x, y)
        r_s, p_s = spearmanr(x, y)
        txt = (
            f"N={len(x)}\n"
            f"Pearson r={r_p:.3f} (p={p_p:.2e})\n"
            f"Spearman rho={r_s:.3f} (p={p_s:.2e})"
        )
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, va="top", fontsize=9,
                bbox=dict(facecolor="white", alpha=0.8, edgecolor="0.8"))

    ax.set_xlabel("GALFIT R_e (PSF, arcsec)")
    ax.set_ylabel("Legacy Tractor shape_r (arcsec)")
    ax.set_title("Legacy vs GALFIT Size Comparison")
    ax.grid(alpha=0.25)
    ax.legend(loc="lower right", fontsize=8)

    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    master = pd.read_csv("master_frb_summary.csv")
    galfit = pd.read_csv("galfit_sigma_metrics_summary.csv")

    frb = master[~master["FRB"].isin(FAILED_LEGACY_FRBS)].copy()
    galfit["re_psf_pix"] = galfit["re_psf"].apply(parse_float)
    galfit["re_psf_arcsec"] = galfit["re_psf_pix"] * PIX_SCALE_ARCSEC
    galfit_map = galfit.set_index("FRB")["re_psf_arcsec"].to_dict()

    svc = pyvo.dal.TAPService("https://datalab.noirlab.edu/tap")

    rows = []
    for i, row in enumerate(frb.itertuples(index=False), start=1):
        rec = {
            "FRB": row.FRB,
            "RA_deg": float(row.RA_deg),
            "DEC_deg": float(row.DEC_deg),
            "re_psf_arcsec": parse_float(galfit_map.get(row.FRB)),
        }

        print(f"[{i}/{len(frb)}] Querying Legacy Tractor for {row.FRB}...")
        match = None
        last_err = None
        for _ in range(3):
            try:
                match = query_nearest_tractor_reff(svc, float(row.RA_deg), float(row.DEC_deg))
                break
            except Exception as exc:
                last_err = str(exc)
                time.sleep(1.5)

        if match is None:
            rec["match_found"] = False
            rec["query_error"] = last_err
            rows.append(rec)
            continue

        rec["match_found"] = True
        rec.update(match)
        rows.append(rec)

    out = pd.DataFrame(rows)
    out_csv = "legacy_vs_galfit_reff_comparison.csv"
    out.to_csv(out_csv, index=False)

    plot_dir = "plots/plots_legacy_cdf"
    os.makedirs(plot_dir, exist_ok=True)
    out_png = f"{plot_dir}/legacy_vs_galfit_reff_scatter.png"
    out_pdf = f"{plot_dir}/legacy_vs_galfit_reff_scatter.pdf"
    make_plot(out, out_png=out_png, out_pdf=out_pdf)

    valid = out[np.isfinite(out["re_psf_arcsec"]) & np.isfinite(out["shape_r_ls_arcsec"])].copy()
    print("\nSaved:", out_csv)
    print("Saved:", out_png)
    print("Saved:", out_pdf)
    print("Matched FRBs:", int(out["match_found"].fillna(False).sum()))
    print("Valid size pairs:", len(valid))
    if len(valid) > 0:
        d = valid["shape_r_ls_arcsec"].to_numpy(dtype=float) - valid["re_psf_arcsec"].to_numpy(dtype=float)
        print("Mean delta (Legacy-GALFIT) arcsec:", round(float(np.mean(d)), 4))
        print("Median delta arcsec:", round(float(np.median(d)), 4))
        print("MAE delta arcsec:", round(float(np.mean(np.abs(d))), 4))


if __name__ == "__main__":
    import os
    main()
