from pathlib import Path

import numpy as np
import pandas as pd


Q0 = 0.2


def inclination_from_q(q: float, q0: float = Q0) -> float:
    if not np.isfinite(q):
        return np.nan
    val = (q * q - q0 * q0) / (1.0 - q0 * q0)
    val = min(1.0, max(0.0, val))
    return float(np.degrees(np.arccos(np.sqrt(val))))


def inclination_error_from_q(q: float, q_err: float, q0: float = Q0) -> float:
    if not np.isfinite(q) or not np.isfinite(q_err) or q_err < 0:
        return np.nan
    q_lo = min(1.0, max(0.0, q - q_err))
    q_hi = min(1.0, max(0.0, q + q_err))
    i0 = inclination_from_q(q, q0=q0)
    ilo = inclination_from_q(q_hi, q0=q0)
    ihi = inclination_from_q(q_lo, q0=q0)
    return float(max(abs(i0 - ilo), abs(ihi - i0)))


def to_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def fmt(val, nd=2):
    if pd.isna(val):
        return "--"
    return f"{float(val):.{nd}f}"


def latex_escape_frb(v: str) -> str:
    return str(v).replace("_", "\\_")


def write_simple_table(rows, headers, out_path: Path, align=None):
    if align is None:
        align = "l" * len(headers)
    lines = []
    lines.append("\\begin{tabular}{%s}" % align)
    lines.append("\\hline")
    lines.append(" {} \\\\".format(" & ".join(headers)))
    lines.append("\\hline")
    for r in rows:
        lines.append(" {} \\\\".format(" & ".join(r)))
    lines.append("\\hline")
    lines.append("\\end{tabular}")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def write_longtable(rows, headers, out_path: Path, align=None, caption="", label=""):
    if align is None:
        align = "l" * len(headers)
    lines = []
    lines.append("\\begin{longtable}{%s}" % align)
    lines.append("\\caption{%s}\\label{%s} \\\\" % (caption, label))
    lines.append("\\hline")
    lines.append(" {} \\\\".format(" & ".join(headers)))
    lines.append("\\hline")
    lines.append("\\endfirsthead")
    lines.append("\\hline")
    lines.append(" {} \\\\".format(" & ".join(headers)))
    lines.append("\\hline")
    lines.append("\\endhead")
    for r in rows:
        lines.append(" {} \\\\".format(" & ".join(r)))
    lines.append("\\hline")
    lines.append("\\end{longtable}")
    out_path.write_text("\n".join(lines), encoding="utf-8")


def main():
    repo = Path(__file__).resolve().parents[1]
    report_dir = repo / "Reports" / "01 PSFEx_Comparison_Report"
    tables_dir = report_dir / "tables"
    report_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)

    cmp_csv = repo / "Archive" / "csv" / "psf" / "psf_constant_comparison.csv"
    galfit_csv = repo / "galfit_sigma_metrics_summary.csv"
    legacy_csv = repo / "legacy_vs_galfit_two_inclinations.csv"

    cmp_df = pd.read_csv(cmp_csv)
    galfit_df = pd.read_csv(galfit_csv)
    legacy_df = pd.read_csv(legacy_csv)

    # Condensed PSFEx-vs-ePSF metrics table
    metric_rows = []
    metric_rows.append(["N(FRB)", fmt(len(cmp_df), 0)])
    metric_rows.append(["Mean |PSFEx Avg FWHM - ePSF Avg FWHM| [pix]", fmt(cmp_df["Delta_Avg_FWHM_PSFEx_minus_ePSF"].abs().mean(), 3)])
    metric_rows.append(["Median |PSFEx Avg FWHM - ePSF Avg FWHM| [pix]", fmt(cmp_df["Delta_Avg_FWHM_PSFEx_minus_ePSF"].abs().median(), 3)])
    metric_rows.append(["Mean |PSFEx Moffat FWHM - ePSF Moffat FWHM| [pix]", fmt(cmp_df["Delta_moffat_fwhm_PSFEx_minus_ePSF"].abs().mean(), 3)])
    metric_rows.append(["Mean |ePSF(exp - Moffat)| [pix]", fmt(cmp_df["Abs_ePSF_exp_minus_moffat_fwhm"].mean(), 3)])
    metric_rows.append(["Mean |PSFEx(exp - Moffat)| [pix]", fmt(cmp_df["Abs_PSFEx_exp_minus_moffat_fwhm"].mean(), 3)])
    metric_rows.append(["Better fit by max frac residual: ePSF count", fmt((cmp_df["Better_Moffat_Fit_by_MaxFracResid"] == "ePSF").sum(), 0)])
    metric_rows.append(["Better fit by max frac residual: PSFEx count", fmt((cmp_df["Better_Moffat_Fit_by_MaxFracResid"] == "PSFEx").sum(), 0)])
    metric_rows.append(["Better exp-vs-Moffat agreement: ePSF count", fmt((cmp_df["Better_ExpVsMoffat_Agreement"] == "ePSF").sum(), 0)])
    metric_rows.append(["Better exp-vs-Moffat agreement: PSFEx count", fmt((cmp_df["Better_ExpVsMoffat_Agreement"] == "PSFEx").sum(), 0)])

    write_simple_table(
        metric_rows,
        headers=["Metric", "Value"],
        out_path=tables_dir / "psfex_summary_metrics.tex",
        align="lp{2.2cm}",
    )

    # Full per-FRB PSFEx assessment table (all FRBs, compact columns to fit page).
    all_headers = [
        "FRB",
        "eAvg",
        "pAvg",
        "dAvg",
        "eMof",
        "pMof",
        "dMof",
        "eResid\\%",
        "pResid\\%",
        "BestFit",
        "BestAgr",
    ]
    all_rows = []
    for _, r in cmp_df.sort_values("FRB").iterrows():
        all_rows.append([
            latex_escape_frb(r["FRB"]),
            fmt(r["Avg_FWHM"], 2),
            fmt(r["PSFEx_Avg_FWHM"], 2),
            fmt(r["Delta_Avg_FWHM_PSFEx_minus_ePSF"], 2),
            fmt(r["ePSF_moffat_fwhm"], 2),
            fmt(r["PSFEx_moffat_fwhm"], 2),
            fmt(r["Delta_moffat_fwhm_PSFEx_minus_ePSF"], 2),
            fmt(r["ePSF_max_frac_resid_pct"], 2),
            fmt(r["PSFEx_max_frac_resid_pct"], 2),
            str(r["Better_Moffat_Fit_by_MaxFracResid"]),
            str(r["Better_ExpVsMoffat_Agreement"]),
        ])

    write_longtable(
        all_rows,
        headers=all_headers,
        out_path=tables_dir / "psfex_all_frbs_compact.tex",
        align="l" + "r" * 8 + "ll",
        caption="All-FRB compact PSFEx versus ePSF comparison table.",
        label="tab:psfex_all_frbs_compact",
    )

    # GALFIT inclination table with propagated errors
    galfit_df["b_a_psf_num"] = to_num(galfit_df["b_a_psf"])
    galfit_df["b_a_err_psf_num"] = to_num(galfit_df["b_a_err_psf"])
    galfit_df["i_psf_deg"] = galfit_df["b_a_psf_num"].apply(inclination_from_q)
    galfit_df["i_psf_err_deg"] = galfit_df.apply(
        lambda r: inclination_error_from_q(r["b_a_psf_num"], r["b_a_err_psf_num"]), axis=1
    )

    galfit_rows = []
    for _, r in galfit_df.sort_values("FRB").iterrows():
        galfit_rows.append([
            latex_escape_frb(r["FRB"]),
            fmt(r["b_a_psf_num"], 3),
            fmt(r["b_a_err_psf_num"], 3),
            fmt(r["i_psf_deg"], 2),
            fmt(r["i_psf_err_deg"], 2),
        ])

    write_longtable(
        galfit_rows,
        headers=["FRB", "$b/a$", "$\\sigma_{b/a}$", "$i_{\\mathrm{GALFIT}}$ [deg]", "$\\sigma_i$ [deg]"],
        out_path=tables_dir / "galfit_inclinations.tex",
        align="lrrrr",
        caption="GALFIT PSF+sigma inclination angles with propagated uncertainties.",
        label="tab:galfit_inclinations",
    )

    # Legacy vs GALFIT comparison table
    legacy_df["galfit_inc_psf_deg"] = to_num(legacy_df["galfit_inc_psf_deg"])
    legacy_df["ls_inc_deg"] = to_num(legacy_df["ls_inc_deg"])
    legacy_df = legacy_df[legacy_df["ls_inc_deg"].notna()].copy()

    # attach galfit propagated uncertainty
    err_map = galfit_df.set_index("FRB")["i_psf_err_deg"].to_dict()
    legacy_df["galfit_i_err_deg"] = legacy_df["FRB"].map(err_map)
    legacy_df["delta_ls_minus_galfit_deg"] = legacy_df["ls_inc_deg"] - legacy_df["galfit_inc_psf_deg"]

    legacy_rows = []
    for _, r in legacy_df.sort_values("FRB").iterrows():
        n_val = r.get("sersic_n_fit")
        n_txt = "--" if pd.isna(n_val) else fmt(n_val, 2)
        legacy_rows.append([
            latex_escape_frb(r["FRB"]),
            str(r.get("type_ls", "--")),
            n_txt,
            fmt(r["ls_inc_deg"], 2),
            fmt(r.get("ls_inc_err_deg"), 2),
            fmt(r["galfit_inc_psf_deg"], 2),
            fmt(r["galfit_i_err_deg"], 2),
            fmt(r["delta_ls_minus_galfit_deg"], 2),
        ])

    write_longtable(
        legacy_rows,
        headers=[
            "FRB",
            "Type",
            "$n_{LS}$",
            "$i_{LS}$",
            "$\\sigma_{i,LS}$",
            "$i_{GF}$",
            "$\\sigma_{i,GF}$",
            "$\\Delta i$",
        ],
        out_path=tables_dir / "legacy_vs_galfit_inclinations.tex",
        align="llrrrrrr",
        caption="Legacy Survey inclination angles versus GALFIT inclinations, with Legacy inclination uncertainty propagated from Tractor shape ivars.",
        label="tab:legacy_vs_galfit",
    )

    tex_lines = [
        "\\documentclass[11pt]{article}",
        "\\usepackage[a4paper,margin=1in]{geometry}",
        "\\usepackage{graphicx}",
        "\\usepackage{booktabs}",
        "\\usepackage{longtable}",
        "\\usepackage{caption}",
        "\\usepackage{float}",
        "\\usepackage{hyperref}",
        "\\setlength{\\textfloatsep}{8pt plus 2pt minus 2pt}",
        "\\setlength{\\intextsep}{8pt plus 2pt minus 2pt}",
        "\\setlength{\\abovecaptionskip}{4pt}",
        "\\setlength{\\belowcaptionskip}{2pt}",
        "\\setlength{\\LTpre}{4pt}",
        "\\setlength{\\LTpost}{4pt}",
        "\\hypersetup{colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue}",
        "",
        "\\title{PSFEx-Centered Inclination Analysis Report (Draft)}",
        "\\author{}",
        "\\date{\\today}",
        "",
        "\\begin{document}",
        "\\maketitle",
        "",
        "\\section{PSFEx Assessment}",
        "\\input{tables/psfex_summary_metrics.tex}",
        "",
        "\\footnotesize",
        "\\setlength{\\tabcolsep}{4pt}",
        "\\input{tables/psfex_all_frbs_compact.tex}",
        "\\normalsize",
        "",
        "\\section{Plots}",
        "The $B$-band proxy used in the multiband comparison follows:",
        "\\[",
        "    B = g + 0.3130\\,(g-r) + 0.2271.",
        "\\]",
        "",
        "\\begin{figure}[htbp]",
        "    \\centering",
        "    \\includegraphics[width=0.92\\linewidth]{../../plots/plots_multiband_cdf/CDF_bias_multiband_rgb_bands_mc_policy.png}",
        "    \\caption{Multi-band CDF comparison (r/g/b nulls) with FRB host Monte Carlo envelope.}",
        "\\end{figure}",
        "",
        "\\begin{figure}[htbp]",
        "    \\centering",
        "    \\includegraphics[width=0.92\\linewidth]{../../plots/plots_astrophot_psf_sigma/CDF_bias_comparison_psf_sigma_updated_mc_both.png}",
        "    \\caption{GALFIT versus AstroPhot comparison where both methods are shown with Monte Carlo uncertainty draws.}",
        "\\end{figure}",
        "",
        "\\section{GALFIT Inclination Angles}",
        "\\input{tables/galfit_inclinations.tex}",
        "",
        "\\section{Legacy Survey Inclinations vs GALFIT}",
        "Legacy/Tractor fields queried for host matching and shape-based inclination were \\texttt{objid}, \\texttt{ra}, \\texttt{dec}, \\texttt{type}, \\texttt{flux\\_r}, \\texttt{sersic}, \\texttt{shape\\_e1}, and \\texttt{shape\\_e2}.",
        "Legacy inclination was derived from shape terms as $|e|=\\sqrt{e_1^2+e_2^2}$, $q=(1-|e|)/(1+|e|)$, and $i=\\cos^{-1}\\!\\left(\\sqrt{(q^2-q_0^2)/(1-q_0^2)}\\right)$ with $q_0=0.2$.",
        "A direct TAP schema query found \\texttt{shape\\_e1\\_ivar} and \\texttt{shape\\_e2\\_ivar}; we convert ivar to standard deviation with $\\sigma=1/\\sqrt{\\mathrm{ivar}}$ and propagate to $\\sigma_{i,\\mathrm{LS}}$ using Monte Carlo draws in $(e_1,e_2)$ space.",
        "\\small",
        "\\setlength{\\tabcolsep}{4pt}",
        "\\renewcommand{\\arraystretch}{1.08}",
        "\\input{tables/legacy_vs_galfit_inclinations.tex}",
        "\\renewcommand{\\arraystretch}{1.0}",
        "\\normalsize",
        "",
        "\\end{document}",
    ]
    tex = "\n".join(tex_lines) + "\n"

    (report_dir / "psfex_galfit_report.tex").write_text(tex, encoding="utf-8")

    print(f"Wrote {report_dir / 'psfex_galfit_report.tex'}")
    print(f"Wrote tables in {tables_dir}")


if __name__ == "__main__":
    main()
