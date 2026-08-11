"""Ad-hoc: headline numbers and outliers from the aggregated tables."""

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import vercommon as vc  # noqa: E402

pd.set_option("display.width", 220)
m = pd.read_csv(os.path.join(vc.TABLES_ROOT, "fit_verification_metrics.csv"),
                dtype={"frb": str})
f = pd.read_csv(os.path.join(vc.TABLES_ROOT, "fit_verification_flags.csv"),
                dtype={"frb": str})
pop = json.load(open(os.path.join(vc.TABLES_ROOT, "population_summary.json")))
d = m[m.in_53]


def num(frame, col):
    return pd.to_numeric(frame[col], errors="coerce")


print("=" * 100)
print("DISTRIBUTIONS (53-host science cut)")
print("=" * 100)
for c in ["chi2nu_global", "chi2nu_local_2re", "sigma_calibration_ratio",
          "rff_2re", "rff_outer_minus_inner", "fourier_dq", "fourier_dpa_deg",
          "dq_sky", "dq_astrophot", "iso_dq_2re", "dmag_ref",
          "model_recon_max_frac", "re_over_fwhm", "sky_sigma_over_pixel_rms"]:
    v = num(d, c).dropna()
    if not len(v):
        continue
    print(f"{c:28s} n={len(v):2d}  median={v.median():+.4f}  "
          f"[p16 {v.quantile(.16):+.4f}, p84 {v.quantile(.84):+.4f}]  "
          f"range [{v.min():+.4f}, {v.max():+.4f}]")

rel = f[f.in_53 & f.fourier_reliable]
print(f"\nFourier estimator usable on {len(rel)}/{int(f.in_53.sum())} science hosts")
if "fourier_unreliable_reasons" in m.columns:
    reasons = m[m.in_53]["fourier_unreliable_reasons"].dropna()
    reasons = reasons[reasons.astype(str) != "nan"]
    counts = {}
    for r in reasons:
        for part in str(r).split(","):
            if part:
                counts[part] = counts.get(part, 0) + 1
    print("  reasons:", counts)
sub = d[[x in set(rel.frb) for x in d.frb]]
v = num(sub, "fourier_dq").dropna()
if len(v):
    print(f"  dq on those: median={v.median():+.4f} "
          f"[p16 {v.quantile(.16):+.4f}, p84 {v.quantile(.84):+.4f}] "
          f"max|dq|={v.abs().max():.4f}")

print("\n" + "=" * 100)
print("POPULATION TESTS (53-host cut)")
print("=" * 100)
for label, keys in [
    ("q vs PSF ellipticity (Spearman)", ("psf_q_vs_epsf_spearman_in53", "psf_q_vs_epsf_p_in53")),
    ("q vs FWHM/Re      (Spearman)", ("psf_q_vs_fwhm_over_re_spearman_in53", "psf_q_vs_fwhm_over_re_p_in53")),
    ("PA_host - PA_psf uniform (KS)", ("psf_dpa_ks_stat_in53", "psf_dpa_ks_p_in53")),
    ("chi2nu_local vs SNR (Spearman)", ("chi2nu_local_vs_snr_spearman_in53", "chi2nu_local_vs_snr_p_in53")),
    ("chi2nu_local vs Re/FWHM", ("chi2nu_local_vs_re_over_fwhm_spearman_in53", "chi2nu_local_vs_re_over_fwhm_p_in53")),
]:
    a, b = pop.get(keys[0]), pop.get(keys[1])
    if a is not None:
        print(f"{label:34s} r={a:+.3f}  p={b:.3g}")
print()
for name in ("re", "n", "sky", "chi2nu"):
    s = pop.get(f"dmag_vs_{name}_slope_in53")
    sig = pop.get(f"dmag_vs_{name}_slope_sig_in53")
    if s is not None:
        print(f"dmag vs {name:7s} slope={s:+.4f}  significance={sig:+.2f} sigma")
print()
for name in ("isophote", "astrophot"):
    o = pop.get(f"q_vs_{name}_median_offset_in53")
    s = pop.get(f"q_vs_{name}_scatter_in53")
    if o is not None:
        print(f"q vs {name:10s} median offset={o:+.4f}  scatter={s:.4f}")

print("\n" + "=" * 100)
print("TRUST TIERS")
print("=" * 100)
print(f.groupby(["in_53", "trust_tier"]).size().unstack(fill_value=0).to_string())

print("\nFlag counts (53-host cut):")
for c in sorted(x for x in f.columns if x.startswith("flag_")):
    n = int(f.loc[f.in_53, c].sum())
    if n:
        print(f"  {c:34s} {n}")

print("\n" + "=" * 100)
print("HOSTS NEEDING A LOOK (tier C or '?', science cut)")
print("=" * 100)
bad = f[f.in_53 & f.trust_tier.isin(["C", "?"])].frb.tolist()
cols = ["frb", "b_a", "re_over_fwhm", "chi2nu_local_2re", "rff_2re",
        "fourier_dq", "fourier_reliable", "iso_dq_2re", "dq_sky",
        "dq_astrophot", "dmag_ref"]
cols = [c for c in cols if c in m.columns]
print(m[m.frb.isin(bad)][cols].to_string(index=False,
                                         float_format=lambda v: f"{v:9.4f}"))

print("\nTop |dq_astrophot| disagreements:")
t = m[m.in_53].reindex(num(m[m.in_53], "dq_astrophot").abs().sort_values(
    ascending=False).index).head(5)
print(t[["frb", "b_a", "ap_q", "dq_astrophot", "re_over_fwhm",
         "chi2nu_local_2re"]].to_string(index=False,
                                        float_format=lambda v: f"{v:9.4f}"))

print("\nWorst dmag:")
t = m[m.in_53].reindex(num(m[m.in_53], "dmag_ref").abs().sort_values(
    ascending=False).index).head(5)
print(t[["frb", "mag", "ref_mag", "dmag_ref", "ref_survey", "n_sersic",
         "re_arcsec"]].to_string(index=False,
                                 float_format=lambda v: f"{v:9.3f}"))

print("\nMost sky-sensitive:")
t = m[m.in_53].reindex(num(m[m.in_53], "dq_sky").sort_values(
    ascending=False).index).head(5)
print(t[["frb", "b_a", "dq_sky", "dq_sky_over_q_err", "re_over_fwhm",
         "sky_sigma_over_pixel_rms"]].to_string(index=False,
                                                float_format=lambda v: f"{v:9.4f}"))
