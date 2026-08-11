# Dust selection and the ~0.6 median axis-ratio / cos(i) trend

## 1. Observed result

Across LS EXP (scaled), HSC, and DES Y1 morph, mag-limited disk pools sit **systematically face-on of isotropic**.

Plot: [`comparisons/median_ba_vs_mag_quad.png`](comparisons/median_ba_vs_mag_quad.png) — median \(b/a\) vs \(m_r\) plateaus near ~0.55–0.7 for LS scaled / HSC / DES over a wide mag range (LS unscaled is lower: REX, not dust).

Mag-cut Hubble cos(i) medians ([`comparisons/summary.csv`](comparisons/summary.csv); \(q_0=0.2\); LS uses ba≤0.8 scaling):

| \(m_r\le\) | LS scaled | HSC | DES |
|------------|-----------|-----|-----|
| 20 | 0.616 | 0.587 | 0.647 |
| 21 | 0.573 | 0.585 | 0.624 |
| 22 | 0.537 | 0.576 | 0.628 |

Isotropic thin disks → median cos(i) ≈ 0.5. Observed ~0.54–0.65 is the shared optical signature to correct.

## 2. Why this happens

Internal dust dims edge-ons: \(m_{\rm obs}=m^f+A(q)\) with \(A\to 0\) as \(q\to 1\). A cut on \(m_{\rm obs}\) under-represents edge-ons and leaves median \(b/a\) / cos(i) high. That is **selection**, not a broken Hubble formula.

Keep separate:

| Effect | Direction | Knob |
|--------|-----------|------|
| Dust + mag limit | too face-on | face-on mag re-cut (below) |
| LS Tractor REX | too edge-on | ba≤0.8 / cos(i) scaling |

Do not treat REX scaling as a dust correction.

## 3. Methods (what works vs what wastes sample)

| Method | Idea | Verdict |
|--------|------|---------|
| **Survival / inclination-complete trim** | Keep only objects that would still pass \(m_{\lim}\) if they suffered max edge-on extinction | Reject as primary — kills a huge fraction near the limit |
| **Unterborn & Ryden A1** | \(\Delta m_r=1.27(\log_{10} q)^2\); select on \(m^f=m-\Delta m\le m_{\lim}\) | **Baseline** — designed for mag-limited exponential disks; adds/retains edge-ons |
| TF linear \(A=\gamma\log(a/b)\) | Same re-cut, linear law | Stage-2 sensitivity only (worse \(q\)-shape; UR08) |
| Volume-limited (\(M\), \(z\)) | Absolute-mag + redshift box | Best where reliable \(z\) exists (HSC subset later) |
| NIR/MIR selection (WISE/2MASS) | Select on dust-weak band; optical for \(b/a\) | Strong; needs cross-match |
| Shao LF / Tuffs RT | Full opacity + LF or RT | Deferred — needs \(z\) / B/T |

Do **not** refit \(\gamma\) from our optical mag-limited pools (Devour & Bell 2016: estimator is itself selection-biased).

**Clarify “face-on trim”:** the costly survival cut discards near-limit face-ons. Unterborn A1 is different: \(m^f=m-\Delta m\) with \(\Delta m\ge 0\), so \(m^f\le m\). A cut \(m^f\le m_{\lim}\) is a **superset** of \(m\le m_{\lim}\): edge-ons that failed the observed cut can enter; face-ons are unchanged.

### What A1 does *not* fix

A1 is still a **single magnitude limit** — choosing \(m_{\lim}=21\) just builds a different (inclination-fairer) sample at that depth. It does **not** automatically debias every differential observed-mag bin:

- Galaxies still “move” in \(m_{\rm obs}\) when dusty; A1 reassigns membership by \(m^f\), not by removing bin migration in \(m_{\rm obs}\).
- Median \(b/a\) vs **observed** mag can stay tilted; the fairer quantity is membership / CDFs in **face-on** mag (or a volume limit).
- Parent catalogs remain flux-limited: A1 can only recover edge-ons that are still in the CSV but outside our science cut, not objects never detected.

So run A1 **separately at each science cut** (20 / 21 / 22). That equalizes inclination mix *up to that cut*; it is not a global fix of the full mag–ba plane.

## 4. Implementation — strict pools, A1 + survival

Script: `scripts/apply_inclination_extinction.py` → [`extinction/`](extinction/).

**Strict:** every CDF object has \(b/a > q_0=0.2\) (asserted; `min_ba` in summary). LS also \(b/a\le 0.8\) + cos(i) rescale (REX). Scaling does **not** change \(\Delta m(q)\), which uses native \(b/a\) — LS remains more edge-on in dust space than DES/HSC.

**Survival:** keep only if \(m_{\rm edge}=m^f+\Delta m(q_0)\le m_{\lim}\) (would still pass if viewed at \(q=q_0\)).

Median cos(i) ([`extinction/funnel.csv`](extinction/funnel.csv)):

| limit | survey | raw | A1 | survival | surv/raw N |
|-------|--------|-----|-----|----------|------------|
| 20 | LS | 0.616 | 0.527 | 0.540 | 0.44 |
| 20 | HSC | 0.587 | 0.521 | 0.480 | 0.40 |
| 20 | DES | 0.647 | 0.590 | 0.596 | 0.53 |
| 21 | LS | 0.573 | 0.509 | 0.520 | 0.49 |
| 21 | HSC | 0.585 | 0.538 | 0.525 | 0.49 |
| 21 | DES | 0.624 | 0.578 | 0.582 | 0.54 |
| 22 | LS | 0.537 | 0.506 | 0.505 | 0.60 |
| 22 | HSC | 0.576 | 0.545 | 0.544 | 0.59 |
| 22 | DES | 0.628 | 0.605 | 0.579 | 0.60 |

Survival cuts **~40–60%** of raw N; A1 **adds** galaxies. Prefer A1 for science CDFs; survival is the costly sanity bound. Plots: `cdfs/*/mag*_before_after.png` (raw/A1/surv), `mag*_survival.png`, `compare/`, `overlay_all_{a1,survival}.png`.

Full competitor list and WP table (archived plan): [`Archive/notes/EXTINCTION_CORRECTION_PLAN.md`](../../../Archive/notes/EXTINCTION_CORRECTION_PLAN.md) — live stub at [`ls_audit/EXTINCTION_CORRECTION_PLAN.md`](ls_audit/EXTINCTION_CORRECTION_PLAN.md). DES morph/calibration context: [`ls_audit/DES_calibration_and_dust_research.md`](ls_audit/DES_calibration_and_dust_research.md).

**LS elliptical-disk geometry (separate from mag-cut dust):** Ryden shape + Unterborn A1 vs Padilla joint \(E_0\) — see [`ls_audit/REX_AND_ELLIPTICAL_DISK.md`](ls_audit/REX_AND_ELLIPTICAL_DISK.md), outputs [`ls_audit/scaled_ryden/`](ls_audit/scaled_ryden/), [`ls_audit/scaled_padilla/`](ls_audit/scaled_padilla/).

## 5. Key references

- Unterborn & Ryden 2008, ApJ 687, 976 — \(\Delta M_r=1.27(\log q)^2\)
- Ryden 2004, ApJ 601, 214 — elliptical disks; lognormal face-on ellipticity
- Padilla & Strauss 2008, MNRAS 388, 1321 — shape + dust self-consistent
- Devour & Bell 2016 — do not refit \(\gamma\) on optically selected samples
- Tully et al. 1998; Shao et al. 2007; Driver et al. 2007 — TF / LF / RT alternatives
