"""WISE colour-colour classification of the 56 hosts (Bhardwaj et al. 2024, Sec. III.1.1; boundaries of Jarrett et al. 2017).

W2-W3 (Vega) < 2.0 spheroid, 2.0-3.5 intermediate disk, > 3.5 star-forming disk; W1-W2 > 0.8 AGN-like (Stern et al. 2012).
Photometry at the GALFIT host centroid:
  LS    Legacy Survey DR10/DR9 forced unWISE photometry (deblended with the optical model): the nearest extended
        source within 2.5", else the nearest point source
  ALLW  AllWISE profile-fit photometry (IRSA), within 3"; used where LS is missing
AB -> Vega: W1 2.67, W2 3.32, W3 5.24 mag (Jarrett et al. 2017). W3 below S/N 3 gives an upper limit on W2-W3.
Colours are observed-frame; Jarrett's divisions are calibrated on nearby galaxies.

Outputs: Publication/materials/appendix/morphology/{wise_colours.csv, wise_colours.pdf, wise_colours.png}
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from foreground_halo_dm import dl_query  # noqa: E402
from host_gallery import host_radec  # noqa: E402

OUT = ROOT / "Publication" / "materials" / "appendix" / "morphology"
CACHE = OUT / "wise_cache"
VEGA_AB = {"w1": 2.67, "w2": 3.32, "w3": 5.24}
SPHEROID, SF_DISK, AGN = 2.0, 3.5, 0.8
CLASS_STYLE = {
    "Sb-Irr (Gini-M20)": ("tab:blue", "o"),
    "Disk (literature)": ("tab:green", "s"),
    "Disk-like profile": ("tab:cyan", "D"),
    "Ambiguous": ("tab:orange", "^"),
    "Early-type": ("tab:red", "*"),
}


def our_class() -> pd.Series:
    g = pd.read_csv(ROOT / "Publication/materials/appendix/robustness/gini_m20_56.csv", dtype={"frb": str}).set_index("frb")
    ev = pd.read_csv(OUT / "disk_evidence.csv", dtype={"frb": str}).set_index("frb")
    cls = pd.Series("Sb-Irr (Gini-M20)", index=g.index)
    cls.loc[ev.index] = ev.label
    return cls


def positions() -> pd.DataFrame:
    seed = pd.read_csv(ROOT / "Publication/materials/host_properties/seed_56.csv", dtype={"frb": str})
    rows = []
    for r in seed.itertuples():
        try:
            ra, dec, _ = host_radec(r.frb)
            src = "galfit"
        except Exception:  # noqa: BLE001
            ra, dec, src = r.ra_deg, r.dec_deg, "frb"
        rows.append(dict(frb=r.frb, z=r.z, ra=ra, dec=dec, pos=src))
    return pd.DataFrame(rows)


def ls_phot(frb, ra, dec) -> dict:
    f = CACHE / f"{frb}_ls.csv"
    if not f.exists():
        d = pd.DataFrame()
        for dr in ("dr10", "dr9"):
            d = dl_query(f"""SELECT ls_id, ra, dec, type, flux_w1, flux_w2, flux_w3, flux_ivar_w1, flux_ivar_w2, flux_ivar_w3,
                             mw_transmission_w1, mw_transmission_w2, mw_transmission_w3,
                             q3c_dist(ra, dec, {ra}, {dec}) * 3600 AS sep FROM ls_{dr}.tractor
                             WHERE q3c_radial_query(ra, dec, {ra}, {dec}, {2.5 / 3600})
                             ORDER BY (type = 'PSF'), sep LIMIT 1""")
            if len(d):
                d["dr"] = dr
                break
        d.to_csv(f, index=False)
    d = pd.read_csv(f)
    if d.empty:
        return {}
    r = d.iloc[0]
    out = {"ls_dr": r.dr, "ls_sep": r.sep}
    for b in ("w1", "w2", "w3"):
        flux = r[f"flux_{b}"] / r[f"mw_transmission_{b}"]
        err = 1 / np.sqrt(r[f"flux_ivar_{b}"]) / r[f"mw_transmission_{b}"] if r[f"flux_ivar_{b}"] > 0 else np.nan
        out[f"ls_{b}_snr"] = flux / err
        out[f"ls_{b}"] = 22.5 - 2.5 * np.log10(max(flux, 3 * err)) - VEGA_AB[b]
    return out


def allwise_phot(frb, ra, dec) -> dict:
    f = CACHE / f"{frb}_allwise.csv"
    if not f.exists():
        from astroquery.ipac.irsa import Irsa
        import astropy.units as u
        from astropy.coordinates import SkyCoord
        t = Irsa.query_region(SkyCoord(ra, dec, unit="deg"), catalog="allwise_p3as_psd", spatial="Cone",
                              radius=3 * u.arcsec, columns="ra,dec,w1mpro,w1snr,w2mpro,w2snr,w3mpro,w3snr,ph_qual,ext_flg")
        d = t.to_pandas()
        if len(d):
            d["sep"] = np.hypot((d.ra - ra) * np.cos(np.radians(dec)), d.dec - dec) * 3600
            d = d.sort_values("sep").head(1)
        d.to_csv(f, index=False)
    d = pd.read_csv(f)
    if d.empty:
        return {}
    r = d.iloc[0]
    out = {"aw_sep": r.sep, "aw_ph_qual": r.ph_qual}
    for b in ("w1", "w2", "w3"):
        out[f"aw_{b}"] = r[f"{b}mpro"]
        out[f"aw_{b}_snr"] = r[f"{b}snr"]
    return out


def classify(w12, w23, upper) -> str:
    if not np.isfinite(w23):
        return "no data"
    if w12 > AGN:
        return "AGN-like"
    if upper:
        return "spheroid" if w23 < SPHEROID else "unconstrained (W3 limit)"
    return "spheroid" if w23 < SPHEROID else "intermediate disk" if w23 < SF_DISK else "star-forming disk"


def build() -> pd.DataFrame:
    CACHE.mkdir(parents=True, exist_ok=True)
    cls = our_class()
    rows = []
    for p in positions().itertuples():
        r = dict(frb=p.frb, z=p.z, ra=p.ra, dec=p.dec, pos=p.pos, our_class=cls[p.frb])
        r.update(ls_phot(p.frb, p.ra, p.dec))
        r.update(allwise_phot(p.frb, p.ra, p.dec))
        for src in ("ls", "aw"):
            if f"{src}_w3" in r and np.isfinite(r[f"{src}_w2"]) and min(r[f"{src}_w1_snr"], r[f"{src}_w2_snr"]) >= 3:
                r[f"{src}_w12"] = r[f"{src}_w1"] - r[f"{src}_w2"]
                r[f"{src}_w23"] = r[f"{src}_w2"] - r[f"{src}_w3"]
                r[f"{src}_upper"] = not r[f"{src}_w3_snr"] >= 3
        src = "ls" if "ls_w23" in r else "aw" if "aw_w23" in r else None
        r["phot"] = {"ls": "LS", "aw": "AllWISE", None: ""}[src]
        r["w12"] = r.get(f"{src}_w12", np.nan)
        r["w23"] = r.get(f"{src}_w23", np.nan)
        r["w23_upper"] = bool(r.get(f"{src}_upper", False))
        r["wise_class"] = classify(r["w12"], r["w23"], r["w23_upper"])
        rows.append(r)
        print(p.frb, r["phot"], f"{r['w23']:.2f}", r["wise_class"], flush=True)
    return pd.DataFrame(rows)


def figure(d: pd.DataFrame) -> None:
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4.2))
    for ax in (a, b):
        ax.axvline(SPHEROID, color="b", ls="--", lw=1) if ax is a else ax.axhline(SPHEROID, color="b", ls="--", lw=1)
        ax.axvline(SF_DISK, color="g", ls="--", lw=1) if ax is a else ax.axhline(SF_DISK, color="g", ls="--", lw=1)
    a.axhline(AGN, color="orangered", ls="--", lw=1)
    for lab, (c, m) in CLASS_STYLE.items():
        s = d[(d.our_class == lab) & d.w23.notna()]
        if s.empty:
            continue
        for up, face in ((False, c), (True, "none")):
            t = s[s.w23_upper == up]
            a.scatter(t.w23, t.w12, c=face, edgecolors=c, marker=m, s=70 if m == "*" else 30,
                      label=lab if not up else None, zorder=3)
            b.scatter(t.z, t.w23, c=face, edgecolors=c, marker=m, s=70 if m == "*" else 30, zorder=3)
            if up:
                a.quiver(t.w23, t.w12, -0.35, 0, angles="xy", scale_units="xy", scale=1, color=c, width=0.003)
                b.quiver(t.z, t.w23, 0, -0.35, angles="xy", scale_units="xy", scale=1, color=c, width=0.003)
    for r in d[(d.our_class.isin(["Early-type", "Ambiguous"]) | (d.frb == "20220509G")) & d.w23.notna()].itertuples():
        a.annotate(r.frb, (r.w23, r.w12), fontsize=6, xytext=(3, 3), textcoords="offset points")
        b.annotate(r.frb, (r.z, r.w23), fontsize=6, xytext=(3, 3), textcoords="offset points")
    a.text(0.6, -0.42, "spheroids", color="b", fontsize=8)
    a.text(2.2, -0.42, "intermediate disks", color="g", fontsize=8)
    a.text(3.9, -0.42, "SF disks", color="firebrick", fontsize=8)
    a.text(0.1, AGN + 0.03, "AGN (Stern+12)", color="orangered", fontsize=7)
    a.set(xlim=(0, 5.5), ylim=(-0.5, 2.0), xlabel="W2 $-$ W3 (Vega)", ylabel="W1 $-$ W2 (Vega)",
          title="(a) WISE colours, Jarrett et al. (2017) divisions")
    a.legend(fontsize=7, loc="upper right", title="our morphology", title_fontsize=7)
    b.set(xlabel="redshift", ylabel="W2 $-$ W3 (Vega)", ylim=(0, 5.5), title="(b) Observed W2 $-$ W3 against redshift")
    b.text(0.98, 0.02, "open: W3 S/N < 3 (upper limit)", transform=b.transAxes, ha="right", fontsize=7)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(OUT / f"wise_colours.{ext}", dpi=200)


def main() -> None:
    warnings.filterwarnings("ignore")
    d = build()
    d.to_csv(OUT / "wise_colours.csv", index=False)
    figure(d)
    print(pd.crosstab(d.our_class, d.wise_class))


if __name__ == "__main__":
    main()
