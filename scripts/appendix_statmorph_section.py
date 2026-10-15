"""Write the statmorph appendix material from statmorph_51_full.csv.

Outputs (Publication/materials/appendix/):
  STATMORPH_MORPHOLOGY.md   narrative + summary tables + per-host table
  statmorph_51.csv          per-host numbers used in the tables
  statmorph_51_table.tex    AASTeX deluxetable of the per-host values
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "pipeline_scripts" / "verification" / "outputs" / "tables" / "statmorph_51_full.csv"
SEED = ROOT / "Publication" / "materials" / "host_properties" / "seed_51.csv"
OUT = ROOT / "Publication" / "materials" / "appendix"

LIT = {
    "20190608B": "barred spiral, FRB on a spiral arm (Chittidi+2021; Mannings+2021 HST)",
    "20231120A": "2MFGC 07440: 2MASS Flat Galaxy Catalogue, edge-on by selection ($a/b \\ge 3$; Mitronova+2004); "
                 "spiral arms (Sharma+2024); edge-on disk (Leung+2025)",
    "20220319D": "IRAS 02044+7048, barred spiral SBa, FRB on the eastern arm (Ravi+2025; Law+2024)",
    "20220509G": "quiescent early-type in Abell 2311 (Sharma+2023; Connor+2023); literature contradicts a disk",
}


def stat_line(x: pd.Series) -> dict:
    x = x.dropna()
    return {"N": len(x), "median": float(x.median()), "RMS": float(np.sqrt(np.mean(x ** 2))),
            "MADn": float(1.4826 * np.median(np.abs(x - x.median()))), "n_gt10": int((x.abs() > 10).sum())}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(SRC, dtype={"frb": str})
    s = pd.read_csv(SEED, dtype={"frb": str}).set_index("frb")
    d["i_astrophot"] = d.frb.map(s.i_astrophot_deg)
    d["di_ap"] = d.i_astrophot - d.i_galfit
    d["di_ap_sm"] = d.i_astrophot - d.i_sm_sersic
    d["n_galfit"] = pd.to_numeric(d.n_galfit, errors="coerce")
    d["disk_F"] = d.sm_F_bulge < 0
    d["lit_note"] = d.frb.map(LIT).fillna("")
    d = d.sort_values("frb")

    cols = ["frb", "m_r", "re_over_fwhm", "q_galfit", "i_galfit", "q_sm_sersic", "i_sm_sersic", "di_sersic",
            "i_astrophot", "di_ap", "n_galfit", "sm_sersic_n", "sm_gini", "sm_m20", "sm_F_bulge", "sm_S_merger",
            "sm_C", "sm_A", "sm_S", "lotz2008", "sm_flag", "sm_flag_sersic", "disk_F", "lit_note"]
    d[cols].to_csv(OUT / "statmorph_51.csv", index=False)

    pairs = {"statmorph Sérsic − GALFIT": d.di_sersic, "AstroPhot − GALFIT": d.di_ap,
             "AstroPhot − statmorph": d.di_ap_sm}
    good = d[(d.sm_flag == 0) & (d.sm_flag_sersic == 0)]
    pairs["statmorph − GALFIT, flag = flag_sersic = 0"] = good.di_sersic

    md = ["# Appendix: non-parametric morphology (statmorph) of the 51 hosts", "",
          "Source: `scripts/statmorph_51_full.py` (statmorph v0.5 run on the adopted GALFIT cutout of each host, "
          "`host_cutout.fits`, with `host_sigma.fits` as weight map, the PSFEx prototype as PSF and the GALFIT star "
          "mask; segmentation = 1.5σ detection on a Gaussian-smoothed image, central segment kept). "
          "Table generator: `scripts/appendix_statmorph_section.py`. Per-host numbers: `statmorph_51.csv`; "
          "LaTeX: `statmorph_51_table.tex`.", "",
          "## A. Independent inclination check", "",
          "statmorph fits its own single Sérsic model (different optimiser, PSF convolution and masking from GALFIT); "
          "AstroPhot is the third, fully independent fitter already used in the paper. All three are converted to "
          "inclination with $\\cos^2 i = (q^2 - q_0^2)/(1 - q_0^2)$, $q_0 = 0.2$.", "",
          "| pair | $N$ | median $\\Delta i$ | RMS | $1.4826\\,$MAD | $|\\Delta i| > 10^\\circ$ |",
          "|---|---|---|---|---|---|"]
    for k, v in pairs.items():
        st = stat_line(v)
        md.append(f"| {k} | {st['N']} | {st['median']:+.2f}° | {st['RMS']:.2f}° | {st['MADn']:.2f}° | {st['n_gt10']} |")
    tri = d[["i_galfit", "i_astrophot", "i_sm_sersic"]].std(axis=1, ddof=1)
    md += ["", f"Per-host standard deviation of the three Sérsic inclinations: median {tri.median():.2f}°, "
               f"mean {tri.mean():.2f}°.", "",
           "The RMS of 4–5° between independent Sérsic codes on the same images is what motivates the $5^\\circ$ "
           "inclination floor; the core scatter is below 1° and the RMS is set entirely by the marginally resolved "
           "tail ($R_e/{\\rm FWHM} < 1$):", ""]
    out = d[(d.di_sersic.abs() > 10) | (d.di_ap.abs() > 10)]
    md += ["| FRB | $R_e$/FWHM | $i$ GALFIT | $i$ statmorph | $i$ AstroPhot | note |", "|---|---|---|---|---|---|"]
    for _, r in out.iterrows():
        note = {"20210410D": "GALFIT $n$ at the 6 cap; statmorph $q = 0.18$",
                "20220509G": "statmorph `flag_sersic = 2`; see B",
                "20220501C": "statmorph `flag = 1`; faint ($m_r = 21.8$)",
                "20190523A": "faint ($m_r = 21.7$), $R_e$/FWHM = 0.55"}.get(r.frb, "")
        ap = f"{r.i_astrophot:.0f}" if np.isfinite(r.i_astrophot) else "—"
        md.append(f"| {r.frb} | {r.re_over_fwhm:.2f} | {r.i_galfit:.0f} | {r.i_sm_sersic:.0f} | {ap} | {note} |")
    md += ["", "Moment-based ellipticities (statmorph `ellipticity_asymmetry`, isophote fits at $1R_e$) are not "
               "PSF-corrected and sit $\\sim 8$–$11^\\circ$ low systematically; they are not used.", "",
           "## B. Disk classification", "",
           "Adopted metric: the Gini–$M_{20}$ bulge statistic $F(G, M_{20})$ of Rodriguez-Gomez et al. (2019), "
           "$F < 0$ = disk-dominated side. Supporting metrics: the Lotz et al. (2008) $G$–$M_{20}$ classes, "
           "concentration $C$ and the Sérsic index.", "",
           "| metric | disk side | of |", "|---|---|---|",
           f"| $F(G,M_{{20}}) < 0$ | {int(d.disk_F.sum())} | {int(d.sm_F_bulge.notna().sum())} |",
           f"| Lotz+2008 late (Sb–Irr) | {int((d.lotz2008 == 'late').sum())} | 51 "
           f"(early {int((d.lotz2008 == 'early').sum())}, merger {int((d.lotz2008 == 'merger').sum())}) |",
           f"| $C < 3.5$ | {int((d.sm_C < 3.5).sum())} | 51 |",
           f"| statmorph $n < 2.5$ | {int((d.sm_sersic_n < 2.5).sum())} | 51 |",
           f"| GALFIT $n < 2.5$ | {int((d.n_galfit < 2.5).sum())} | {int(d.n_galfit.notna().sum())} |", "",
           "Hosts not on the disk side of $F$, with the literature morphology:", "",
           "| FRB | $F$ | Lotz | $C$ | $n$ GALFIT / statmorph | literature |", "|---|---|---|---|---|---|"]
    for _, r in d[~d.disk_F].iterrows():
        md.append(f"| {r.frb} | {r.sm_F_bulge:+.2f} | {r.lotz2008} | {r.sm_C:.2f} | {r.n_galfit:.1f} / {r.sm_sersic_n:.1f} | {LIT.get(r.frb, '')} |")
    m = d[d.lotz2008 == "merger"]
    for _, r in m.iterrows():
        md.append(f"| {r.frb} | {r.sm_F_bulge:+.2f} | {r.lotz2008} ($A = {r.sm_A:.2f}$) | {r.sm_C:.2f} | "
                  f"{r.n_galfit:.1f} / {r.sm_sersic_n:.1f} | {LIT.get(r.frb, '')} |")
    md += ["", "All three $F \\ge 0$ hosts are spirals in the literature (two barred spirals, one catalogued flat "
               "edge-on galaxy whose concentration is raised by projection). The one host whose disk nature is "
               "contradicted rather than supported is FRB 20220509G: statmorph puts it on the disk side "
               "($F = -0.58$) but with an unconstrained Sérsic shape (GALFIT $n$ at the 0.5 floor, statmorph "
               "$n = 5.2$, `flag_sersic = 2`, $A = 0.60$), and Sharma+2023 / Connor+2023 classify it as a "
               "quiescent early-type galaxy. It was part of the 2024 sample ($i = 60^\\circ$).", "",
           "## C. Per-host table", "",
           "| FRB | $m_r$ | $R_e$/FWHM | $q$ GALFIT | $i$ GALFIT | $q$ sm | $i$ sm | $i$ AP | $n$ GALFIT | $n$ sm | $G$ | $M_{20}$ | $F$ | $C$ | $A$ | Lotz | flags |",
           "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for _, r in d.iterrows():
        ap = f"{r.i_astrophot:.0f}" if np.isfinite(r.i_astrophot) else "—"
        md.append(f"| {r.frb} | {r.m_r:.2f} | {r.re_over_fwhm:.2f} | {r.q_galfit:.2f} | {r.i_galfit:.0f} | {r.q_sm_sersic:.2f} | "
                  f"{r.i_sm_sersic:.0f} | {ap} | {r.n_galfit:.2f} | {r.sm_sersic_n:.2f} | {r.sm_gini:.2f} | {r.sm_m20:.2f} | "
                  f"{r.sm_F_bulge:+.2f} | {r.sm_C:.2f} | {r.sm_A:.2f} | {r.lotz2008} | {int(r.sm_flag)}/{int(r.sm_flag_sersic)} |")
    (OUT / "STATMORPH_MORPHOLOGY.md").write_text("\n".join(md) + "\n", encoding="utf-8")

    # LaTeX deluxetable
    tex = [r"\begin{longrotatetable}", r"\begin{deluxetable*}{lcccccccccccccc}", r"\tabletypesize{\scriptsize}",
           r"\tablecaption{Non-parametric morphology of the 51 hosts (statmorph)\label{tab:statmorph}}",
           r"\tablehead{\colhead{FRB} & \colhead{$m_r$} & \colhead{$R_e$/FWHM} & \colhead{$i_{\rm GALFIT}$} & "
           r"\colhead{$i_{\rm statmorph}$} & \colhead{$i_{\rm AstroPhot}$} & \colhead{$n_{\rm GALFIT}$} & "
           r"\colhead{$n_{\rm statmorph}$} & \colhead{$G$} & \colhead{$M_{20}$} & \colhead{$F(G,M_{20})$} & "
           r"\colhead{$C$} & \colhead{$A$} & \colhead{Lotz} & \colhead{flags} \\ "
           r"\colhead{} & \colhead{AB} & \colhead{} & \colhead{deg} & \colhead{deg} & \colhead{deg} & \colhead{} & "
           r"\colhead{} & \colhead{} & \colhead{} & \colhead{} & \colhead{} & \colhead{} & \colhead{} & \colhead{}}",
           r"\startdata"]
    for _, r in d.iterrows():
        ap = f"{r.i_astrophot:.0f}" if np.isfinite(r.i_astrophot) else r"\nodata"
        tex.append(f"{r.frb} & {r.m_r:.2f} & {r.re_over_fwhm:.2f} & {r.i_galfit:.0f} & {r.i_sm_sersic:.0f} & {ap} & "
                   f"{r.n_galfit:.2f} & {r.sm_sersic_n:.2f} & {r.sm_gini:.2f} & {r.sm_m20:.2f} & {r.sm_F_bulge:+.2f} & "
                   f"{r.sm_C:.2f} & {r.sm_A:.2f} & {r.lotz2008} & {int(r.sm_flag)}/{int(r.sm_flag_sersic)} \\\\")
    tex += [r"\enddata",
            r"\tablecomments{statmorph \citep{2019MNRAS.483.4140R} run on the adopted $r$-band cutout of each host with the "
            r"GALFIT star mask and PSFEx PSF. Inclinations use $\cos^2 i = (q^2-q_0^2)/(1-q_0^2)$ with $q_0 = 0.2$ and "
            r"$q = 1 - $ ellipticity of the statmorph S\'ersic fit. $F(G,M_{20}) < 0$ is the disk-dominated side; "
            r"Lotz classes follow \citet{2008ApJ...672..177L}. Flags are statmorph \texttt{flag}/\texttt{flag\_sersic} "
            r"(0 = clean).}",
            r"\end{deluxetable*}", r"\end{longrotatetable}"]
    (OUT / "statmorph_51_table.tex").write_text("\n".join(tex) + "\n", encoding="utf-8")
    print("written", OUT)


if __name__ == "__main__":
    main()
