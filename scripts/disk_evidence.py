"""Case-by-case disk evidence for the 32 hosts whose Gini-M20 class is unreliable or not Sb-Irr.

Classes
  lit      literature morphology or a profile fit to deeper imaging shows a disk / spiral
  profile  our free-n PSF-convolved Sersic index n < 2.5 (no contrary literature)
  ambig    n >= 2.5 in an unresolved host with no morphological information, or quiescent with n ~ 2
  early    early-type in the literature

n_free is the production fit with n free; the adopted fit for some hosts fixes n = 1.

Outputs: Publication/materials/appendix/morphology/{disk_evidence.csv, disk_evidence.tex}
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "Publication" / "materials" / "appendix" / "morphology"

# frb: (class, literature evidence (LaTeX), short note for gallery)
EVIDENCE = {
    "20180924B": ("lit", r"HST: spiral arms and bar in residuals \citep{mannings2021}; $b/a = 0.62$ disk-like profile \citep{bannister2019}", "spiral (HST)"),
    "20181112A": ("lit", r"$n = 1.9$ \citep{prochaska2019}", "disk-like (lit. n=1.9)"),
    "20190102C": ("profile", r"star-forming \citep{Gordon2023}", "profile"),
    "20190523A": ("profile", r"massive, low SFR \citep{20190523AHost}", "profile"),
    "20190608B": ("lit", r"barred spiral; rotating gas disk, $i_{\rm gas} = 37^\circ$ \citep{chittidi2021}", "barred spiral"),
    "20210410D": ("ambig", r"transitioning \citep{Gordon2023}", "ambiguous"),
    "20211203C": ("lit", r"$n = 0.6$ in deeper $R$-band imaging \citep{gordon2025}", "disk (lit. n=0.6)"),
    "20220207C": ("lit", r"disk-dominated, nearly edge-on \citep{law2024}", "edge-on disk"),
    "20220310F": ("profile", r"star-forming \citep{sharma2024}", "profile"),
    "20220319D": ("lit", r"barred spiral (SBa), FRB on a spiral arm \citep{ravi2025}", "barred spiral"),
    "20220501C": ("profile", r"star-forming \citep{glowacki2025}", "profile"),
    "20220509G": ("early", r"early-type, quiescent \citep{sharma2023, connor2023}", "early-type"),
    "20220825A": ("profile", r"star-forming \citep{sharma2024}", "profile"),
    "20220914A": ("lit", r"late-type spiral \citep{connor2023}", "spiral"),
    "20221012A": ("ambig", r"massive, quiescent \citep{law2024, sharma2024}", "ambiguous"),
    "20221219A": ("lit", r"Milky Way-like late type \citep{faber2024}", "late type"),
    "20230124A": ("profile", r"star-forming \citep{sharma2024}", "profile"),
    "20230526A": ("lit", r"disk, $n = 0.3$ ($R$) and 1.4 ($K$) \citep{gordon2025}", "disk (lit. n=0.3)"),
    "20230613A": ("profile", r"\nodata", "profile"),
    "20230626A": ("ambig", r"transition \citep{sharma2024}", "ambiguous"),
    "20230628A": ("profile", r"star-forming \citep{sharma2024}", "profile"),
    "20230902A": ("lit", r"disk, $n = 1.0$ \citep{gordon2025}", "disk (lit. n=1.0)"),
    "20231020B": ("profile", r"\nodata", "profile"),
    "20231120A": ("lit", r"spiral arms \citep{sharma2024}; edge-on disk \citep{leung2025}; flat-galaxy catalog \citep{mitronova2004}", "edge-on spiral"),
    "20231220A": ("profile", r"\nodata", "profile"),
    "20240119A": ("profile", r"\nodata", "profile"),
    "20240203D": ("ambig", r"\nodata", "ambiguous"),
    "20240208A": ("profile", r"star-forming \citep{glowacki2025}", "profile"),
    "20240215A": ("profile", r"\nodata", "profile"),
    "20240224A": ("profile", r"\nodata", "profile"),
    "20240229A": ("ambig", r"\nodata", "ambiguous"),
    "20240304A": ("lit", r"disk, $n = 1.2$ \citep{gordon2025}", "disk (lit. n=1.2)"),
}
LABEL = {"lit": "Disk (literature)", "profile": "Disk-like profile", "ambig": "Ambiguous", "early": "Early-type"}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    seed = pd.read_csv(ROOT / "Publication/materials/host_properties/seed_56.csv", dtype={"frb": str}).set_index("frb")
    sm = pd.read_csv(ROOT / "pipeline_scripts/verification/outputs/tables/statmorph_56_full.csv", dtype={"frb": str}).set_index("frb")
    g = pd.read_csv(ROOT / "Publication/materials/appendix/robustness/gini_m20_56.csv", dtype={"frb": str}).set_index("frb")
    rows = []
    for frb, (cls, lit, note) in EVIDENCE.items():
        rows.append(dict(frb=frb, n_free=seed.loc[frb, "n_galfit"], n_statmorph=sm.loc[frb, "sm_sersic_n"],
                         re_over_fwhm=g.loc[frb, "re_over_fwhm"], gm20=g.loc[frb, "lotz2008"],
                         why=g.loc[frb, "why_unreliable"], cls=cls, label=LABEL[cls], literature=lit, note=note))
    d = pd.DataFrame(rows)
    assert set(d.frb) == set(g.index[(~g.reliable) | (g.lotz2008 != "late")])
    bad = d[(d.cls == "profile") & (d.n_free >= 2.5)]
    assert bad.empty, bad
    d.rename(columns={"note": "class"}).to_csv(OUT / "disk_evidence.csv", index=False)
    print(d.cls.value_counts().to_dict())
    gm = {"late": "Sb--Irr", "early": "E/S0/Sa", "merger": "Merger"}
    tex = []
    for r in d.sort_values(["cls", "frb"], key=lambda s: s.map({"lit": 0, "profile": 1, "ambig": 2, "early": 3}) if s.name == "cls" else s).itertuples():
        tex.append(f"{r.frb} & {r.n_free:.1f} & {r.n_statmorph:.1f} & {r.re_over_fwhm:.2f} & {gm[r.gm20]} & "
                   f"{r.label} & {r.literature} \\\\")
    (OUT / "disk_evidence.tex").write_text("\n".join(tex) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
