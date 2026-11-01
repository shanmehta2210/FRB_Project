"""Main-text inclination tests repeated with alternative axis-ratio measurements of the same hosts.

Variants: adopted GALFIT (PSF-convolved), AstroPhot (PSF-convolved), GALFIT and AstroPhot without PSF
convolution (pipeline_scripts/verification/outputs/nopsf/nopsf_refits.csv), and, for the hosts shared
with Bhardwaj et al. (2024), their published AutoProf inclinations against ours. Each sample is tested
against the HSC disk control with the matched-draw p_m of the main text (same seeds).

Output: Publication/materials/appendix/fit_variants/{fit_variants.csv, FIT_VARIANTS.md}
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import appendix_robustness as A  # noqa: E402
import plot_paper_inclination_figures as P  # noqa: E402

ROOT = P.REPO_ROOT
SEED56 = ROOT / "Publication" / "materials" / "host_properties" / "seed_56.csv"
NOPSF = ROOT / "pipeline_scripts" / "verification" / "outputs" / "nopsf" / "nopsf_refits.csv"
OUT = ROOT / "Publication" / "materials" / "appendix" / "fit_variants"
N_DRAWS = 10_000


def cosi_from_i(i_deg):
    return np.cos(np.radians(np.asarray(i_deg, float)))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    hosts = P.load_hosts()
    seed = pd.read_csv(SEED56).set_index("frb")
    h = hosts.set_index("frb")
    pool = A.cosi_from_q(A.select(A.load_parent())["ba"].to_numpy(), A.Q0)
    variants = [("GALFIT, PSF-convolved (adopted)", h.cosi_point.to_numpy(float), h.index.to_numpy())]
    q_ap = seed.loc[h.index, "q_astrophot"].to_numpy(float)
    variants.append(("AstroPhot, PSF-convolved", A.cosi_from_q(q_ap, A.Q0), h.index.to_numpy()))
    if NOPSF.exists():
        n = pd.read_csv(NOPSF).set_index("frb")
        for col, lab in (("q_galfit_nopsf", "GALFIT, no PSF"), ("q_astrophot_nopsf", "AstroPhot, no PSF")):
            q = n.reindex(h.index)[col].to_numpy(float)
            ok = np.isfinite(q)
            variants.append((lab, A.cosi_from_q(q[ok], A.Q0), h.index.to_numpy()[ok]))
    shared = seed.loc[h.index, "i_bhardwaj2024"].notna().to_numpy()
    i24 = seed.loc[h.index[shared], "i_bhardwaj2024"].astype(str).str.split(",").str[0].astype(float)
    variants.append((f"Shared with Bhardwaj et al. (2024): this work", h.cosi_point.to_numpy(float)[shared], h.index[shared].to_numpy()))
    variants.append((f"Shared with Bhardwaj et al. (2024): their AutoProf values", cosi_from_i(i24), h.index[shared].to_numpy()))
    rows = []
    for lab, c, frbs in variants:
        r = A.matched_tests(c, pool, None, N_DRAWS)
        r.update(label=lab)
        rows.append(r)
        print(lab, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()
                    if k in ("N_frb", "median_i_frb", "quartile_counts", "ad_pm", "mw2_pm", "mw1_pm", "chi2_pm")})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "fit_variants.csv", index=False)
    md = ["# Inclination tests with alternative axis-ratio measurements", "",
          f"Control: HSC disk pool ($N = {len(pool):,}$). Matched-draw $p_{{\\rm m}}$ from {N_DRAWS:,} draws (main-text seeds).", "",
          "| measurement | $N$ | median $i$ | quartile counts (face-on to edge-on) | AD $p_{\\rm m}$ | MW 2-sided $p_{\\rm m}$ | MW 1-sided $p_{\\rm m}$ | $\\chi^2$ $p_{\\rm m}$ |",
          "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        md.append(f"| {r['label']} | {r['N_frb']} | {r['median_i_frb']:.1f}° | {r['quartile_counts']} | {r['ad_pm']:.2g} | "
                  f"{r['mw2_pm']:.2g} | {r['mw1_pm']:.2g} | {r['chi2_pm']:.2g} |")
    if NOPSF.exists():
        nn = pd.read_csv(NOPSF).set_index("frb").reindex(i24.index)
        cmp = [("GALFIT, PSF-convolved", seed.loc[i24.index, "i_galfit_deg"]),
               ("GALFIT, no PSF", nn.i_galfit_nopsf_deg), ("AstroPhot, no PSF", nn.i_astrophot_nopsf_deg)]
        md += ["", f"## Agreement with the published inclinations of Bhardwaj et al. (2024) ({len(i24)} shared hosts)", "",
               "| our measurement | median $i_{\\rm ours} - i_{\\rm B24}$ | RMS difference |", "|---|---|---|"]
        out = []
        for lab, v in cmp:
            d = v.to_numpy(float) - i24.to_numpy(float)
            out.append(dict(measurement=lab, median_diff=np.nanmedian(d), rms=np.sqrt(np.nanmean(d ** 2))))
            md.append(f"| {lab} | {out[-1]['median_diff']:+.1f}° | {out[-1]['rms']:.1f}° |")
        pd.DataFrame(out).to_csv(OUT / "bhardwaj2024_agreement.csv", index=False)
    (OUT / "FIT_VARIANTS.md").write_text("\n".join(md) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
