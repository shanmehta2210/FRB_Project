from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    root = Path(__file__).resolve().parents[2]
    outdir = root / "Bias_Report"
    outdir.mkdir(exist_ok=True)
    nl = r"\\"

    ap = pd.read_csv(root / "AstroPhot_Analysis/results/astrophot_psf_sigma_n1_inclination_angles.csv")
    gf = pd.read_csv(root / "galfit_sigma_metrics_summary.csv").rename(columns={"FRB": "frb_name"})
    gf["b_a_psf_num"] = pd.to_numeric(gf["b_a_psf"].astype(str).str.replace("*", "", regex=False), errors="coerce")

    q0 = 0.2
    val = (gf["b_a_psf_num"] ** 2 - q0**2) / (1 - q0**2)
    val = np.clip(val, 0, 1)
    gf["inc_psf_deg"] = np.degrees(np.arccos(np.sqrt(val)))
    gf.loc[gf["b_a_psf_num"] <= q0, "inc_psf_deg"] = 90.0

    merged = ap.merge(gf[["frb_name", "re_psf", "b_a_psf_num", "inc_psf_deg"]], on="frb_name", how="inner")
    merged = merged.sort_values("frb_name").copy()
    for col in ["Re", "q", "inclination_angle", "re_psf", "b_a_psf_num", "inc_psf_deg"]:
        merged[col] = pd.to_numeric(merged[col], errors="coerce")

    rows = []
    for _, row in merged.iterrows():
        rows.append(
            f"{row['frb_name']} & {row['Re']:.2f} & {row['re_psf']:.2f} & {row['q']:.3f} & "
            f"{row['b_a_psf_num']:.3f} & {row['inclination_angle']:.2f} & {row['inc_psf_deg']:.2f} {nl}"
        )

    table1_lines = [
        "\\begin{table}[H]",
        "\\centering",
        "\\scriptsize",
        "\\caption{AstroPhot (PSF+sigma, fixed $n=1$) vs GALFIT (PSF+sigma): side-by-side structural parameters.}",
        "\\label{tab:astrophot_galfit_side_by_side}",
        "\\begin{tabular}{lrrrrrr}",
        "\\toprule",
        f"FRB & $R_{{e,\\mathrm{{AP}}}}$ & $R_{{e,\\mathrm{{GF}}}}$ & $(b/a)_{{\\mathrm{{AP}}}}$ & $(b/a)_{{\\mathrm{{GF}}}}$ & $i_{{\\mathrm{{AP}}}}$ (deg) & $i_{{\\mathrm{{GF}}}}$ (deg) {nl}",
        "\\midrule",
        *rows,
        "\\bottomrule",
        "\\end{tabular}",
        "\\end{table}",
    ]
    (outdir / "table_astrophot_galfit_side_by_side.tex").write_text("\n".join(table1_lines), encoding="ascii")

    sm = pd.read_csv(root / "statmorph_nonparam_results.csv").copy().sort_values("frb_name")
    for col in ["concentration", "asymmetry", "smoothness", "gini", "m20"]:
        sm[col] = pd.to_numeric(sm[col], errors="coerce")

    sm_rows = []
    for _, row in sm.iterrows():
        sm_rows.append(
            f"{row['frb_name']} & {row['concentration']:.3f} & {row['asymmetry']:.3f} & "
            f"{row['smoothness']:.3f} & {row['gini']:.3f} & {row['m20']:.3f} {nl}"
        )

    table2_lines = [
        "\\begin{longtable}{lrrrrr}",
        f"\\caption{{Statmorph non-parametric morphology for all FRB hosts (CAS, Gini, M20).}}\\label{{tab:statmorph_casginim20}}{nl}",
        "\\toprule",
        f"FRB & C & A & S & Gini & M20 {nl}",
        "\\midrule",
        "\\endfirsthead",
        "\\toprule",
        f"FRB & C & A & S & Gini & M20 {nl}",
        "\\midrule",
        "\\endhead",
        "\\midrule",
        f"\\multicolumn{{6}}{{r}}{{Continued on next page}} {nl}",
        "\\midrule",
        "\\endfoot",
        "\\bottomrule",
        "\\endlastfoot",
        *sm_rows,
        "\\end{longtable}",
    ]
    (outdir / "table_statmorph_casginim20.tex").write_text("\n".join(table2_lines), encoding="ascii")

    print(outdir / "table_astrophot_galfit_side_by_side.tex")
    print(outdir / "table_statmorph_casginim20.tex")


if __name__ == "__main__":
    main()
