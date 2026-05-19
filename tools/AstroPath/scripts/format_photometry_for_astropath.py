import argparse
import os
import numpy as np
import pandas as pd


def build_candidates(df: pd.DataFrame, min_ang_size: float) -> pd.DataFrame:
    required = ["RA", "Dec", "MAG_CALIB_APER_40PX", "FLUX_RADIUS"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    out = pd.DataFrame(
        {
            "objid": df.get("objid_ls", pd.Series([-1] * len(df))).fillna(-1).astype("int64"),
            "ra": pd.to_numeric(df["RA"], errors="coerce"),
            "dec": pd.to_numeric(df["Dec"], errors="coerce"),
            "type": df.get("type_ls", pd.Series(["UNKNOWN"] * len(df))).fillna("UNKNOWN"),
            "mag": pd.to_numeric(df["MAG_CALIB_APER_40PX"], errors="coerce"),
            "ang_size": pd.to_numeric(df["FLUX_RADIUS"], errors="coerce"),
        }
    )

    # Keep likely galaxies (non-stellar) and valid coordinates/magnitudes.
    out = out[out["type"].astype(str).str.upper() != "PSF"].copy()
    out = out[np.isfinite(out["ra"]) & np.isfinite(out["dec"]) & np.isfinite(out["mag"])].copy()

    # PATH needs strictly positive sizes; fallback for invalid entries.
    out.loc[~np.isfinite(out["ang_size"]) | (out["ang_size"] <= 0), "ang_size"] = min_ang_size
    out["ang_size"] = out["ang_size"].astype(float)
    out["source_catalog"] = "LS_DR10"

    out = out.sort_values("mag", ascending=True).reset_index(drop=True)
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Format photometry CSV into AstroPath candidate schema.")
    parser.add_argument(
        "--input",
        default="tools/Photometry/r70_target_comparison_photometry.csv",
        help="Input photometry CSV path.",
    )
    parser.add_argument(
        "--field-id",
        default="r70",
        help="Field identifier used in output filename.",
    )
    parser.add_argument(
        "--output-dir",
        default="tools/astropath/data",
        help="Directory for AstroPath candidate CSV output.",
    )
    parser.add_argument(
        "--min-ang-size",
        type=float,
        default=0.5,
        help="Minimum fallback angular size used when FLUX_RADIUS is invalid.",
    )
    args = parser.parse_args()

    input_path = os.path.abspath(args.input)
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input CSV not found: {input_path}")

    df = pd.read_csv(input_path)
    out = build_candidates(df, min_ang_size=args.min_ang_size)

    os.makedirs(args.output_dir, exist_ok=True)
    output_path = os.path.join(args.output_dir, f"{args.field_id}_candidates.csv")
    out.to_csv(output_path, index=False)

    print(f"Wrote {len(out)} candidates to {output_path}")
    print(f"Columns: {', '.join(out.columns)}")


if __name__ == "__main__":
    main()
