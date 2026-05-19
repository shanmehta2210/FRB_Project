import math
import numpy as np
import pandas as pd


def hubble_inc_from_ba(b_a, q0=0.2):
    if b_a <= q0:
        return 90.0
    val = (b_a**2 - q0**2) / (1 - q0**2)
    val = max(0.0, min(1.0, val))
    return math.degrees(math.acos(math.sqrt(val)))


ap = pd.read_csv("AstroPhot_Analysis/results/astrophot_psf_sigma_inclination_angles.csv")
gf = pd.read_csv("galfit_sigma_metrics_summary.csv")

# Keep only PSF-fit columns from GALFIT and compute inclination consistently from b/a.
gf = gf.rename(columns={"FRB": "frb_name"})
gf["inc_psf"] = gf["b_a_psf"].apply(lambda x: hubble_inc_from_ba(float(str(x).replace("*", ""))))

merged = ap.merge(
    gf[["frb_name", "b_a_psf", "inc_psf", "n_psf", "re_psf", "chi2nu_psf"]],
    on="frb_name",
    how="inner",
)

merged["delta_q"] = merged["q"] - merged["b_a_psf"]
merged["delta_inc_deg"] = merged["inclination_angle"] - merged["inc_psf"]
merged["delta_n"] = merged["n"] - merged["n_psf"]
merged["delta_Re"] = merged["Re"] - merged["re_psf"]
merged["delta_chi2nu"] = merged["chi2_nu"] - merged["chi2nu_psf"]

# Correlation checks for similarity (higher absolute agreement for direct quantities).
summary_rows = []
for ap_col, gf_col, label in [
    ("q", "b_a_psf", "q vs b_a_psf"),
    ("inclination_angle", "inc_psf", "inclination vs inc_psf"),
    ("n", "n_psf", "n vs n_psf"),
    ("Re", "re_psf", "Re vs re_psf"),
]:
    use = merged[[ap_col, gf_col]].replace([np.inf, -np.inf], np.nan).dropna()
    summary_rows.append(
        {
            "comparison": label,
            "N": len(use),
            "pearson_r": use[ap_col].corr(use[gf_col], method="pearson"),
            "spearman_rho": use[ap_col].corr(use[gf_col], method="spearman"),
            "mean_abs_diff": (use[ap_col] - use[gf_col]).abs().mean(),
            "median_abs_diff": (use[ap_col] - use[gf_col]).abs().median(),
        }
    )

summary = pd.DataFrame(summary_rows)
summary.to_csv("AstroPhot_Analysis/results/astrophot_vs_galfit_psf_similarity_summary.csv", index=False)

per_object = merged[
    [
        "frb_name",
        "q",
        "b_a_psf",
        "delta_q",
        "inclination_angle",
        "inc_psf",
        "delta_inc_deg",
        "n",
        "n_psf",
        "delta_n",
        "Re",
        "re_psf",
        "delta_Re",
        "chi2_nu",
        "chi2nu_psf",
        "delta_chi2nu",
    ]
].copy()
per_object = per_object.sort_values("delta_inc_deg", key=lambda s: s.abs(), ascending=False)
per_object.to_csv("AstroPhot_Analysis/results/astrophot_vs_galfit_psf_per_object.csv", index=False)

print("Saved AstroPhot_Analysis/results/astrophot_vs_galfit_psf_similarity_summary.csv")
print("Saved AstroPhot_Analysis/results/astrophot_vs_galfit_psf_per_object.csv")
print("\nSimilarity summary:")
print(summary.to_string(index=False))

print("\nLargest inclination disagreements (top 8):")
print(
    per_object[["frb_name", "inclination_angle", "inc_psf", "delta_inc_deg", "q", "b_a_psf", "delta_q"]]
    .head(8)
    .to_string(index=False)
)
