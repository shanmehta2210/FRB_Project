# LS REX bias and elliptical-disk scaling (Ryden / Padilla)

Merged from former `REX_INCLINATION_RESEARCH.md` + `ELLIPTICAL_DISK_PLAN.md` (those paths are stubs).

---

## Part A — REX inclination bias (why EXP lacks face-ons)

### A.1 What REX actually does (root cause, confirmed)

Legacy Survey Tractor assigns a source the type **REX** ("round exponential galaxy")
— an exponential profile **forced to zero ellipticity (b/a = 1)** — *unless* an
elliptical model (EXP/DEV/SER) improves the penalised chi2 by >~9 (about 3σ)
[legacysurvey.org/dr10/description, /catalogs; legacypipe `RexGalaxy`]. Consequences:

- A round / near-face-on disk is well fit by the round REX model, so it **never gets
  promoted to EXP**. Selecting `type=EXP` therefore **systematically discards face-on
  disks** — they are sitting in the REX bucket with no shape measurement.
- REX is also the default for faint / marginally resolved objects, so the excluded
  face-on population is huge: DR10 has **1.12e9 REX vs 2.25e8 EXP** (~5× more REX).

**Empirical proof in our catalog** (`LS_catalog_v2_fullsky_exp.csv`, N=2e6):

| b/a percentile | 50 | 90 | 99 | 99.9 | max |
|---------------:|---:|---:|---:|-----:|----:|
| value | 0.445 | 0.595 | 0.722 | **0.803** | 0.988 |

Only **0.11%** of EXP galaxies have b/a > 0.8 and **0.01%** have b/a > 0.9. The EXP
sample is effectively truncated at b/a ≈ 0.8 — the face-on edge is missing.

### A.2 How the literature (incl. this paper) avoids it

**Bhardwaj et al. 2024 (this repo's paper; Nature 634, 1065 / arXiv:2408.01876):**

- Inclination via the **Hubble formula** with **q0 = 0.2**; **discard b/a < 0.2**.
- **They never use Tractor `type=EXP`.** FRB-host b/a comes from their own
  **Sérsic + elliptical-isophote fitting** on the images — no REX truncation.
- **The null is NOT a theoretical uniform.** They compare against a **control sample of
  SDSS DR16 field disk galaxies** with the identical Hubble formula. Shared biases cancel
  in the two-sample test. SDSS exponential-fit b/a is measured for all galaxies (REX-free).

Takeaway: solve REX by (a) measuring b/a with a method that fits all galaxies, or
(b) using a **matched control** so residual shape bias cancels — not by trusting a
Tractor-EXP inclination in absolute terms.

### A.3 Options for a good inclination distribution from LS

1. **Match the paper (recommended).** Control whose b/a is measured like the FRB hosts
   (SDSS expAB, or re-fit Sérsic). Compare FRB-vs-control with AD/KS; not vs uniform.
2. **Un-truncate Tractor by adding REX back.** Assign REX b/a ≈ 1 and merge with EXP.
3. **Band-aid: `scaled`** — renormalise to the empirical face-on edge (below).

### A.4 Is `scaled` (b/a → (b/a)/0.8) physically justified? — Yes, first-order

`scaled` declares **b/a = 0.8 to be face-on** and stretches cos(i) back to [0,1]. This is
justified *because REX truncates the EXP sample exactly at b/a ≈ 0.8* (99.9th pct = 0.803).

Why `scaled` beats Ryden *here*:

- **Refit Ryden** absorbs the truncation → cos(i) collapses to ~uniform
  (`scaled_ryden/CIRCULARITY_CHECK.md`).
- **Fixed-lit Ryden** reads truncation as real inclination → too edge-on (~0.41).
- **`scaled`** encodes the verified face-on edge at 0.8.

**Caveats:** linear stretch, does not recover real face-ons, assumes sharp edge at 0.8.
Adequate as a **matched transform** on FRB + LS control; not a substitute for measuring
face-ons. See also `scaled_ryden_fixed/SCALED_IS_DEGENERATE_RYDEN.md`,
`scaled_ryden_capped/CAPPED_RYDEN_EXPERIMENT.md`.

### A.5 Recommendation

Apply one b/a→cos(i) map **identically** to FRB hosts and LS (or SDSS) control, then
two-sample test. Use **`scaled`** for REX-truncated LS EXP. Treat Ryden fits as shape
sanity checks, not the inclination null.

---

## Part B — Elliptical-disk scaling: Ryden vs Padilla

Physically motivated replacements for the ad-hoc [`scaled/`](scaled/) stretch
(`b/a ≤ 0.8`, `(b/a)/0.8`). Two parallel tracks:

| Directory | Shape | Dust |
|-----------|-------|------|
| [`scaled_ryden/`](scaled_ryden/) | Ryden (2004) 4-param elliptical disk | Unterborn \(\Delta m=1.27(\log_{10} q)^2\) re-cut |
| [`scaled_padilla/`](scaled_padilla/) | Same triaxial disk + Padilla \(E_0\) | Joint via inclination weight \(\psi(\theta)\) |

Ad-hoc `scaled/` kept for comparison.

### B.1 Literature

**Ryden 2004** (ApJ 601, 214): disks intrinsically elliptical. Thickness \(\gamma=C/A\) ~
Gaussian; face-on ellipticity \(\varepsilon=1-B/A\) ~ lognormal. Project with Binney (1985).
SDSS \(i\)-band seed: \(\mu_\gamma\approx 0.22\), \(\sigma_\gamma\approx 0.06\),
\(\mu=\ln\varepsilon\approx -1.85\), \(\sigma\approx 0.89\).

**Padilla & Strauss 2008** (MNRAS 388, 1321): same geometry + planar dust.
\(E(\theta)=E_0(1+y-\cos\theta)\) (capped at \(E_0\)); viewing angles with completeness
\(\psi(\theta)\). Spirals: \(\ln e\approx -2.33\), \(E_0\simeq 0.45\) mag. Full paper uses
\(1/V_{\max}\) + LF; **LS has no z**, so we use \(\psi(\theta)\propto 10^{-0.4 E(\theta)}\).

**Unterborn & Ryden 2008**: selection re-cut \(m^f=m-1.27(\log_{10} q)^2\le m_{\lim}\) —
used on the Ryden track after the shape fit (not inside Padilla).

### B.2 Scripts

```bash
python scripts/fit_ls_scaled_elliptical.py --mode ryden
python scripts/fit_ls_scaled_elliptical.py --mode padilla
python scripts/fit_ls_scaled_elliptical.py --mode both
```

Core: `scripts/elliptical_disk_model.py`.

### B.3 FRB reuse

Load `scaled_ryden/fit_params.json` or `scaled_padilla/fit_params.json` for
\(P(\cos i\mid q)\). Apply Unterborn only on the Ryden track; Padilla folds dust into \(E_0\).

### B.4 Knobs (keep separate)

1. Elliptical-disk geometry (why \(q\sim 1\) is rare)
2. Dust selection (Unterborn A1 or Padilla \(E_0\)) — see [`../DUST_AND_MEDIAN_BA.md`](../DUST_AND_MEDIAN_BA.md)
3. Ad-hoc REX stretch (`scaled/`) — legacy only
