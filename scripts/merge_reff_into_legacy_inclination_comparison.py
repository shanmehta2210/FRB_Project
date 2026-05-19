import pandas as pd


def main() -> None:
    incl_path = "legacy_vs_galfit_inclination_comparison.csv"
    reff_path = "legacy_vs_galfit_reff_comparison.csv"

    incl = pd.read_csv(incl_path)
    reff = pd.read_csv(reff_path)

    keep_cols = [
        "FRB",
        "re_psf_arcsec",
        "shape_r_ls_arcsec",
        "sep_arcsec",
        "tractor_objid",
    ]
    keep_cols = [c for c in keep_cols if c in reff.columns]

    reff_small = reff[keep_cols].copy()
    rename_map = {
        "sep_arcsec": "sep_arcsec_reff",
        "tractor_objid": "tractor_objid_reff",
    }
    reff_small = reff_small.rename(columns=rename_map)

    merged = incl.drop(
        columns=[c for c in ["re_psf_arcsec", "shape_r_ls_arcsec", "delta_reff_arcsec", "sep_arcsec_reff", "tractor_objid_reff"] if c in incl.columns],
        errors="ignore",
    ).merge(reff_small, on="FRB", how="left")

    if {"shape_r_ls_arcsec", "re_psf_arcsec"}.issubset(merged.columns):
        merged["delta_reff_arcsec"] = merged["shape_r_ls_arcsec"] - merged["re_psf_arcsec"]

    merged.to_csv(incl_path, index=False)

    n_with_reff = int(merged["shape_r_ls_arcsec"].notna().sum()) if "shape_r_ls_arcsec" in merged.columns else 0
    print(f"Updated: {incl_path}")
    print(f"Rows: {len(merged)}")
    print(f"Rows with Reff pair: {n_with_reff}")


if __name__ == "__main__":
    main()
