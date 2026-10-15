"""Full statmorph pass over the 51-host working sample.

Re-runs statmorph on the adopted cutout of each host (``wkdir`` in
``plots/paper_figures/fig_host_inputs.csv``), using the same segmentation
recipe as ``pipeline_scripts/galfit_fitting/run_statmorph_pipeline.py`` plus
the GALFIT star mask, and stores the complete output: Sersic fit (n, ellip,
r_half), moment ellipticities, Gini, M20, C, A, S, F(G,M20), S(G,M20) and flags.

Then compares the statmorph-derived inclinations with the GALFIT ones and
summarises the disk/bulge classification.

Outputs: pipeline_scripts/verification/outputs/tables/statmorph_51_full.csv
         pipeline_scripts/verification/outputs/tables/statmorph_51_summary.md
"""
from __future__ import annotations

import json
import math
import os
import sys
import warnings

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "pipeline_scripts", "galfit_fitting"))
sys.path.insert(0, os.path.join(ROOT, "pipeline_scripts"))
from run_statmorph_pipeline import _clamp_weightmap, _load_fits_data, _make_segmap  # noqa: E402

OUT = os.path.join(ROOT, "pipeline_scripts", "verification", "outputs", "tables")
Q0 = 0.2


def inc_deg(q: float, q0: float = Q0) -> float:
    if not np.isfinite(q) or q <= q0:
        return 90.0
    if q >= 1:
        return 0.0
    return math.degrees(math.acos(math.sqrt((q * q - q0 * q0) / (1 - q0 * q0))))


def lotz2008(g: float, m20: float) -> str:
    if not (np.isfinite(g) and np.isfinite(m20)):
        return "unknown"
    if g > -0.14 * m20 + 0.33:
        return "merger"
    if g > 0.14 * m20 + 0.80:
        return "early"
    return "late"


def run_one(wkdir: str) -> dict:
    import statmorph
    from astropy.stats import sigma_clipped_stats

    image = _load_fits_data(os.path.join(wkdir, "host_cutout.fits"))
    _, med, std = sigma_clipped_stats(image, sigma=3.0)
    image_sub = image - med
    segmap = _make_segmap(image_sub, std)
    kw: dict = {}
    sp = os.path.join(wkdir, "host_sigma.fits")
    if os.path.isfile(sp):
        kw["weightmap"] = _clamp_weightmap(_load_fits_data(sp))
    pp = os.path.join(wkdir, "proto_image.fits")
    if os.path.isfile(pp):
        kw["psf"] = _load_fits_data(pp)
    mp = os.path.join(wkdir, "host_mask.fits")
    if os.path.isfile(mp):
        kw["mask"] = _load_fits_data(mp) > 0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        m = statmorph.source_morphology(image_sub, segmap, **kw)[0]

    def f(v):
        try:
            v = float(v)
            return v if math.isfinite(v) else np.nan
        except (TypeError, ValueError):
            return np.nan

    return {
        "sm_sersic_n": f(m.sersic_n), "sm_sersic_ellip": f(m.sersic_ellip),
        "sm_sersic_rhalf_px": f(m.sersic_rhalf), "sm_sersic_theta": f(m.sersic_theta),
        "sm_sersic_chi2_dof": f(m.sersic_chi2_dof),
        "sm_ellip_centroid": f(m.ellipticity_centroid), "sm_ellip_asym": f(m.ellipticity_asymmetry),
        "sm_elong_asym": f(m.elongation_asymmetry),
        "sm_gini": f(m.gini), "sm_m20": f(m.m20), "sm_C": f(m.concentration), "sm_A": f(m.asymmetry),
        "sm_S": f(m.smoothness), "sm_F_bulge": f(m.gini_m20_bulge), "sm_S_merger": f(m.gini_m20_merger),
        "sm_rpetro_circ_px": f(m.rpetro_circ), "sm_rhalf_circ_px": f(m.rhalf_circ),
        "sm_sn_per_pixel": f(m.sn_per_pixel), "sm_flag": int(m.flag), "sm_flag_sersic": int(m.flag_sersic),
    }


def main() -> None:
    hosts = pd.read_csv(os.path.join(ROOT, "plots", "paper_figures", "fig_host_inputs.csv"), dtype={"frb": str})
    prod = pd.read_csv(os.path.join(ROOT, "pipeline_galfit_results.csv"), dtype={"frb": str}).set_index("frb")
    rows = []
    for _, h in hosts.iterrows():
        wk = os.path.join(ROOT, str(h.wkdir).replace("/", os.sep))
        row = {"frb": h.frb, "wkdir": h.wkdir, "q_galfit": float(h.q), "q_galfit_err": float(h.q_err),
               "i_galfit": float(h.i_deg), "re_over_fwhm": float(h.re_over_fwhm), "m_r": float(h.mag)}
        if h.frb in prod.index:
            for c in ("n", "sersic_n", "n_fixed", "n_sersic"):
                if c in prod.columns:
                    row["n_galfit"] = prod.loc[h.frb, c]
                    break
        try:
            row.update(run_one(wk))
        except Exception as exc:  # noqa: BLE001
            row["error"] = str(exc)
        rows.append(row)
        print(h.frb, "ok" if "error" not in row else row["error"])
    d = pd.DataFrame(rows)
    d["q_sm_sersic"] = 1 - d.sm_sersic_ellip
    d["q_sm_moment"] = 1 - d.sm_ellip_asym
    d["i_sm_sersic"] = d.q_sm_sersic.map(inc_deg)
    d["i_sm_moment"] = d.q_sm_moment.map(inc_deg)
    d["di_sersic"] = d.i_sm_sersic - d.i_galfit
    d["di_moment"] = d.i_sm_moment - d.i_galfit
    d["lotz2008"] = [lotz2008(g, m) for g, m in zip(d.sm_gini, d.sm_m20)]
    os.makedirs(OUT, exist_ok=True)
    d.to_csv(os.path.join(OUT, "statmorph_51_full.csv"), index=False)

    def stats(x: pd.Series) -> str:
        x = x.dropna()
        rms = math.sqrt(np.mean(x ** 2))
        mad = 1.4826 * np.median(np.abs(x - np.median(x)))
        return (f"N={len(x)}, mean={x.mean():+.2f}, median={x.median():+.2f}, RMS={rms:.2f}, "
                f"std={x.std(ddof=1):.2f}, 1.4826*MAD={mad:.2f}, |d|>10 deg: {(x.abs() > 10).sum()}")

    good = d[(d.sm_flag == 0)]
    good_s = d[(d.sm_flag == 0) & (d.sm_flag_sersic == 0)]
    lines = ["# statmorph vs GALFIT, 51-host working sample", "",
             "Inclinations from cos^2 i = (q^2 - q0^2)/(1 - q0^2), q0 = 0.2, q = 1 - ellipticity.", "",
             f"- Sersic-fit ellipticity, all hosts: {stats(d.di_sersic)}",
             f"- Sersic-fit ellipticity, flag=0 & flag_sersic=0 ({len(good_s)}): {stats(good_s.di_sersic)}",
             f"- Moment (asymmetry-centre) ellipticity, all: {stats(d.di_moment)}",
             f"- Moment ellipticity, flag=0 ({len(good)}): {stats(good.di_moment)}", "",
             "## Classification counts", "",
             f"- Lotz+2008 G-M20: {d.lotz2008.value_counts().to_dict()}",
             f"- F(G,M20) bulge statistic < 0 (disk side): {(d.sm_F_bulge < 0).sum()} / {d.sm_F_bulge.notna().sum()}",
             f"- statmorph Sersic n < 2.5: {(d.sm_sersic_n < 2.5).sum()}; GALFIT n < 2.5: {(pd.to_numeric(d.get('n_galfit'), errors='coerce') < 2.5).sum()}",
             f"- Concentration C < 3.5 (Conselice 2003 late-type side): {(d.sm_C < 3.5).sum()}",
             ""]
    (open(os.path.join(OUT, "statmorph_51_summary.md"), "w", encoding="utf-8")).write("\n".join(lines))
    print("\n".join(lines))


if __name__ == "__main__":
    main()
