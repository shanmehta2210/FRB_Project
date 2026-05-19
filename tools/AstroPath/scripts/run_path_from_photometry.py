import argparse
import os
import sys
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u


sys.path.append("tools/AstroPath/astropath_pkg")
from astropath import path  # type: ignore


def run_path(
    candidates: pd.DataFrame,
    ra_deg: float,
    dec_deg: float,
    a_arcsec: float,
    b_arcsec: float,
    pa_deg: float,
    theta_max: float,
    theta_scale: float,
    step_size: float,
    max_radius: float,
) -> tuple[pd.DataFrame, dict]:
    frb_coord = SkyCoord(ra=ra_deg, dec=dec_deg, unit="deg", frame="icrs")

    mypath = path.PATH()
    mypath.init_localization(
        "eellipse",
        center_coord=frb_coord,
        eellipse={"a": a_arcsec, "b": b_arcsec, "theta": pa_deg},
    )

    mypath.init_candidates(
        ra=candidates["ra"].to_numpy(dtype=float),
        dec=candidates["dec"].to_numpy(dtype=float),
        ang_size=candidates["ang_size"].to_numpy(dtype=float),
        mag=candidates["mag"].to_numpy(dtype=float),
    )
    mypath.init_cand_prior(P_O_method="inverse", P_U=0.1)
    mypath.init_theta_prior(PDF="exp", max=theta_max, scale=theta_scale)
    mypath.calc_priors()
    p_oix, p_ux = mypath.calc_posteriors(method="local", step_size=step_size, max_radius=max_radius)

    out = candidates.copy()
    out["prior_O"] = mypath.prior_Oi
    out["posterior_O"] = p_oix
    out = out.sort_values("posterior_O", ascending=False).reset_index(drop=True)

    best = out.iloc[0]
    best_coord = SkyCoord(ra=float(best["ra"]), dec=float(best["dec"]), unit="deg")
    sep_arcsec = frb_coord.separation(best_coord).to(u.arcsec).value

    summary = {
        "field_id": "",
        "frb_ra_deg": ra_deg,
        "frb_dec_deg": dec_deg,
        "loc_a_arcsec": a_arcsec,
        "loc_b_arcsec": b_arcsec,
        "loc_pa_deg": pa_deg,
        "P_U": float(p_ux),
        "best_objid": int(best["objid"]) if np.isfinite(best["objid"]) else -1,
        "best_mag": float(best["mag"]),
        "best_ang_size": float(best["ang_size"]),
        "best_posterior": float(best["posterior_O"]),
        "best_sep_arcsec": float(sep_arcsec),
    }
    return out, summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Run AstroPath PATH on formatted photometry candidates.")
    parser.add_argument("--field-id", default="r70", help="Field identifier for I/O names.")
    parser.add_argument(
        "--input-csv",
        default="tools/astropath/data/r70_candidates.csv",
        help="Formatted AstroPath candidate CSV path.",
    )
    parser.add_argument("--ra", type=float, default=64.3996075, help="FRB/field center RA in deg.")
    parser.add_argument("--dec", type=float, default=7.9311060, help="FRB/field center Dec in deg.")
    parser.add_argument("--loc-a", type=float, default=0.5, help="Localization semi-major axis in arcsec.")
    parser.add_argument("--loc-b", type=float, default=0.5, help="Localization semi-minor axis in arcsec.")
    parser.add_argument("--loc-pa", type=float, default=0.0, help="Localization PA in degrees.")
    parser.add_argument("--theta-max", type=float, default=60.0, help="Theta prior max arcsec.")
    parser.add_argument("--theta-scale", type=float, default=2.0, help="Theta prior scale.")
    parser.add_argument("--step-size", type=float, default=0.1, help="PATH local posterior step size.")
    parser.add_argument("--max-radius", type=float, default=60.0, help="PATH max radius arcsec.")
    parser.add_argument(
        "--results-dir",
        default="tools/astropath/results",
        help="Output directory for posterior and summary CSVs.",
    )
    args = parser.parse_args()

    input_csv = os.path.abspath(args.input_csv)
    if not os.path.exists(input_csv):
        raise FileNotFoundError(f"Candidate CSV not found: {input_csv}")

    candidates = pd.read_csv(input_csv)
    needed = {"objid", "ra", "dec", "mag", "ang_size", "source_catalog"}
    missing = sorted(needed - set(candidates.columns))
    if missing:
        raise ValueError(f"Input CSV missing required columns: {missing}")
    if candidates.empty:
        raise ValueError("Input CSV has no candidates after filtering.")

    posterior_df, summary = run_path(
        candidates=candidates,
        ra_deg=args.ra,
        dec_deg=args.dec,
        a_arcsec=args.loc_a,
        b_arcsec=args.loc_b,
        pa_deg=args.loc_pa,
        theta_max=args.theta_max,
        theta_scale=args.theta_scale,
        step_size=args.step_size,
        max_radius=args.max_radius,
    )
    summary["field_id"] = args.field_id

    os.makedirs(args.results_dir, exist_ok=True)
    posterior_path = os.path.join(args.results_dir, f"{args.field_id}_posterior.csv")
    summary_path = os.path.join(args.results_dir, f"{args.field_id}_summary.csv")
    posterior_df.to_csv(posterior_path, index=False)
    pd.DataFrame([summary]).to_csv(summary_path, index=False)

    print(f"Wrote posterior table: {posterior_path}")
    print(f"Wrote summary table: {summary_path}")


if __name__ == "__main__":
    main()
