"""
Extract EXP-analogue (exponential-disk) subsamples from DES Y1 morph and
Kawinwanichakij HSC catalogs.

Neither product has Tractor ``type=EXP`` or SDSS ``lnLExp``. Both ship free-n
single-Sérsic fits, so the EXP analogue is a Sérsic-index window around n=1
(the exponential profile), plus quality flags where available.

Primary cut (written under ``catalog/``):
  DES:  FIT_AVAILABLE_R=1  AND  0.4 < n_r < 1.5
  HSC:  goodfits_flag=1    AND  0.4 < fitted_sersic < 1.5
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

_SCRIPTS = Path(__file__).resolve().parent
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))

from null_catalog_utils import (
    DES_Y1_MORPH_EXP_DEFAULT,
    DES_Y1_MORPH_SAMPLE_DEFAULT,
    EXP_ANALOGUE_N_MAX,
    EXP_ANALOGUE_N_MIN,
    HSC_KAWIN_EXP_DEFAULT,
    HSC_KAWIN_SAMPLE_DEFAULT,
    LS_CATALOG_V2_EXP_DEFAULT,
    REPO_ROOT,
)


def pcts(x, qs=(5, 25, 50, 75, 95)):
    x = pd.to_numeric(x, errors="coerce").to_numpy(dtype=float)
    x = x[np.isfinite(x)]
    return np.round(np.percentile(x, qs), 4) if len(x) else None


def block(title: str, df: pd.DataFrame, ba: pd.Series, n: pd.Series, mag: pd.Series, n0: int):
    ba = pd.to_numeric(ba, errors="coerce")
    n = pd.to_numeric(n, errors="coerce")
    mag = pd.to_numeric(mag, errors="coerce")
    print(f"\n{title}")
    print(f"  N={len(df):,}  ({100.0 * len(df) / n0:.1f}% of parent 500k)")
    print(f"  ba   p5/25/50/75/95 = {pcts(ba)}")
    print(f"  n    p5/25/50/75/95 = {pcts(n)}")
    print(
        f"  mag  p25/50/75      = "
        f"{mag.quantile(0.25):.3f} / {mag.median():.3f} / {mag.quantile(0.75):.3f}"
    )


def main() -> None:
    n0 = 500_000
    des_in = REPO_ROOT / DES_Y1_MORPH_SAMPLE_DEFAULT
    hsc_in = REPO_ROOT / HSC_KAWIN_SAMPLE_DEFAULT
    des_out = REPO_ROOT / DES_Y1_MORPH_EXP_DEFAULT
    hsc_out = REPO_ROOT / HSC_KAWIN_EXP_DEFAULT
    des_out.parent.mkdir(parents=True, exist_ok=True)

    nmin, nmax = EXP_ANALOGUE_N_MIN, EXP_ANALOGUE_N_MAX
    primary_label = f"{nmin}<n<{nmax}  [PRIMARY]"

    # ---- DES ----
    d = pd.read_csv(des_in)
    n = pd.to_numeric(d["n_r"], errors="coerce")
    ba = pd.to_numeric(d["ba_r"], errors="coerce")
    mag = pd.to_numeric(d["mag_r"], errors="coerce")
    gr = pd.to_numeric(d["MAG_SERSIC_G"], errors="coerce") - pd.to_numeric(
        d["MAG_SERSIC_R"], errors="coerce"
    )
    base = (d["FIT_AVAILABLE_R"] == 1) & ba.notna() & n.notna() & (n > 0) & (n < 20)

    print("=" * 60)
    print(f"DES Y1 morph — EXP analogue via Sérsic n in ({nmin}, {nmax})")
    print("=" * 60)
    block("parent (valid fit)", d.loc[base], ba[base], n[base], mag[base], n0)

    cuts = {
        f"n<{nmax}": base & (n < nmax),
        primary_label: base & (n > nmin) & (n < nmax),
        f"{nmin}<n<{nmax} + (g-r)<0.8": base
        & (n > nmin)
        & (n < nmax)
        & gr.notna()
        & (gr < 0.8),
        f"{nmin}<n<{nmax} + (g-r)<0.6": base
        & (n > nmin)
        & (n < nmax)
        & gr.notna()
        & (gr < 0.6),
    }
    for label, m in cuts.items():
        block(label, d.loc[m], ba[m], n[m], mag[m], n0)

    primary = cuts[primary_label]
    d.loc[primary].to_csv(des_out, index=False)
    print(f"\nWrote {des_out}  ({primary.sum():,} rows)")

    # ---- HSC ----
    h = pd.read_csv(hsc_in)
    ns = pd.to_numeric(h["fitted_sersic"], errors="coerce")
    q = pd.to_numeric(h["fitted_q"], errors="coerce")
    magh = pd.to_numeric(h["fitted_mag"], errors="coerce")
    gr_h = pd.to_numeric(h["gmag"], errors="coerce") - pd.to_numeric(h["rmag"], errors="coerce")
    good = pd.to_numeric(h["goodfits_flag"], errors="coerce") == 1
    sf = pd.to_numeric(h["quiescent_flag"], errors="coerce") == 0
    base_h = q.notna() & ns.notna() & (ns > 0) & (ns < 10)

    print("\n" + "=" * 60)
    print(f"HSC Kawinwanichakij — EXP analogue via fitted_sersic in ({nmin}, {nmax})")
    print("=" * 60)
    block("parent (valid n,q)", h.loc[base_h], q[base_h], ns[base_h], magh[base_h], n0)
    block(
        "goodfits=1",
        h.loc[base_h & good],
        q[base_h & good],
        ns[base_h & good],
        magh[base_h & good],
        n0,
    )

    cuts_h = {
        f"goodfits + n<{nmax}": base_h & good & (ns < nmax),
        f"goodfits + {primary_label}": base_h & good & (ns > nmin) & (ns < nmax),
        f"goodfits + {nmin}<n<{nmax} + SF": base_h
        & good
        & (ns > nmin)
        & (ns < nmax)
        & sf,
        f"goodfits + {nmin}<n<{nmax} + SF + (g-r)<0.8": base_h
        & good
        & (ns > nmin)
        & (ns < nmax)
        & sf
        & gr_h.notna()
        & (gr_h < 0.8),
        f"goodfits + {nmin}<n<{nmax} + (g-r)<0.6": base_h
        & good
        & (ns > nmin)
        & (ns < nmax)
        & gr_h.notna()
        & (gr_h < 0.6),
    }
    for label, m in cuts_h.items():
        block(label, h.loc[m], q[m], ns[m], magh[m], n0)

    primary_h = cuts_h[f"goodfits + {primary_label}"]
    h.loc[primary_h].to_csv(hsc_out, index=False)
    print(f"\nWrote {hsc_out}  ({primary_h.sum():,} rows)")

    # ---- LS reference ----
    ls_path = REPO_ROOT / LS_CATALOG_V2_EXP_DEFAULT
    ba_ls_med = None
    if ls_path.is_file():
        ls = pd.read_csv(ls_path, nrows=500_000)
        ba_ls = pd.to_numeric(ls["expAB_r"], errors="coerce")
        ba_ls_med = ba_ls.median()
        print("\n" + "=" * 60)
        print("LS Tractor type=EXP reference (first 500k of 2M catalog)")
        print("=" * 60)
        print("  N=500,000  (true EXP by construction)")
        print(f"  ba (expAB_r) p5/25/50/75/95 = {pcts(ba_ls)}")

    print("\n" + "=" * 60)
    print("SUMMARY (median b/a)")
    print("=" * 60)
    if ba_ls_med is not None:
        print(f"  LS  type=EXP                              : {ba_ls_med:.4f}")
    print(f"  HSC goodfits + {nmin}<n<{nmax}  [PRIMARY]  : {q[primary_h].median():.4f}")
    print(f"  DES {nmin}<n<{nmax}             [PRIMARY]  : {ba[primary].median():.4f}")


if __name__ == "__main__":
    main()
