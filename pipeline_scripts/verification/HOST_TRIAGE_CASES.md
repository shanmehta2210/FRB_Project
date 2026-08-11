# Host triage cases

Detailed per-FRB arguments for (or against) including GALFIT geometry in the
paper. Short user notes live in `host_confirmation.csv` (`notes` column);
this file is the long-form case record.

Workflow / CSV conventions:
[`HOST_CONFIRMATION_WORKFLOW.md`](HOST_CONFIRMATION_WORKFLOW.md).
Doc index: [`VERIFICATION_README.md`](VERIFICATION_README.md).

`confirmed = True` means residuals / discrepancies are judged **physically
reasonable** and the host’s \(q\) (etc.) may go into the paper directly.

Triage walks **all fitted panels** (production 64). The mag / \(b/a\) science
cut is applied later — do not skip faint hosts during this pass
([`SCIENCE_CUT_AND_COHORT.md`](SCIENCE_CUT_AND_COHORT.md)).

---

## 20171020A — CONFIRMED

**Panel:** `outputs/panels/20171020A.png`

### Host character
Bright, huge, fully resolved disk: \(m=15.10\), \(R_e=45.9\) px (\(11.5''\)),
\(R_e/\mathrm{FWHM}\approx 12\), SNR \(\sim 800\).
\(q=0.399\pm 0.001\), PA \(=-50.8^\circ\), \(n=0.54\). Trust tier A.

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.000\) | Independent code agrees |
| Sky \(\pm1\sigma\) | \(q_+=0.398\), \(q_-=0.399\) | Sky-insensitive |
| Isophotes | \(q(1R_e)\approx 0.41\) data / \(0.40\) model; `iso_dq_2re`$=+0.003$ | Observed shape matches convolved model |
| Fourier | reliable; \(\delta q=-0.006\pm 0.001\) | No global \(q\) bias |

### Residual dirtiness (expected, not anti-\(q\))
- \(\chi^2/\nu_{\rm corr}\sim 2.2\) (global), \(\sim 3.4\) (\(2R_e\)): normal for a mag-15 structured galaxy.
- \(\mathrm{RFF}_{2R_e}=+0.080\): real residual light; below the 0.10 flag.
- Residual map: spiral / ring pattern; central positive core under-fit by \(n=0.54\).
- \(\psi_2'\approx +76^\circ/R_e\); \(m_{1..4}\approx 0.06,0.08,0.19,0.10\): winding \(m=2\) = spiral, so collapsed \(\delta q\) stays tiny while modes are non-zero.

### Flags that fire but do not kill \(q\)
- `flag_dmag`: \(\Delta m_{\rm ref-m}=+0.68\) vs PS1 — aperture mismatch for an \(11''\) galaxy, not shape failure.
- `flag_fourier_dq_significant`: \(|\delta q|/\sigma\sim 4\) but \(|\delta q|=0.006\) — negligible absolute bias.

### Decision
**CONFIRMED.** Include \(q\) (and PA, \(R_e\), \(n\)) in the paper. Caveat only on quoting GALFIT \(m\) vs PS1 without an aperture note.

---

## 20180924B — CONFIRMED (resolution caveat)

**Panel:** `outputs/panels/20180924B.png`

### Host character
Faint, **barely resolved**: \(m=20.24\), \(R_e=3.5\) px (\(0.92''\)),
\(R_e/\mathrm{FWHM}\approx 0.79\), SNR \(\sim 2.2\).
\(q=0.679\pm 0.014\), PA \(=-24.8^\circ\), \(n=3.91\). Trust tier A.

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q_{\rm AP-G}=-0.019\) (\(q_{\rm AP}=0.660\)) | Independent fitter agrees within \(\sim 0.02\) |
| Sky \(\pm1\sigma\) | \(q_+=0.674\), \(q_-=0.684\) | \(\lvert\Delta q_{\rm sky}\rvert\approx 0.005\) — stable |
| Isophotes | \(q_{\rm data}(1R_e)\approx 0.82\), \(q_{\rm model}\approx 0.89\); `iso_dq_2re`\(\approx -0.046\) | Data vs convolved model agree; both rounder than GALFIT \(q\) via PSF |
| Fourier | unreliable; ignore \(\delta q\) | no leverage |

### Decision
**CONFIRMED.** Caveat: barely resolved; do not quote Fourier \(\delta q\); formal
\(\sigma_q\) is a lower bound. Residual clump is a neighbour/knot, not a wrong
global ellipse.

---

## 20181112A — CONFIRMED (n=1 fixed re-fit)

**Production panel:** `outputs/panels/20181112A.png` (unchanged; free-\(n\) / free-sky).
**Confirmed panel:** `outputs/panels/20181112A_n1.png` ← `Re-fits/20181112A/panel_n1.png`.

### Host character
Very faint, **barely resolved**: production \(m=21.98\), \(R_e=2.7\) px,
\(R_e/\mathrm{FWHM}\approx 0.56\). Bright neighbour masked. Lotz+2008: **merger**.

### Re-fit (`n1`: \(n=1\) fixed, sky free)
| test | result | reading |
|---|---|---|
| Geometry | \(q=0.862\pm 0.049\), PA \(=-11.9^\circ\), \(R_e=2.7\) px, \(n=1\) | Matches production \(q\) with \(n\) locked |
| AstroPhot | \(\Delta q_{\rm AP-G}=-0.007\) (\(n\) also locked) | Agrees |
| Sky \(\pm1\sigma\) | \(q_{\rm sky+}\!=\!0.858\), \(q_{\rm sky-}\!=\!0.865\) | \(\lvert\Delta q_{\rm sky}\rvert\approx 0.004\) — stable with \(n\) fixed |
| RFF\(_{2R_e}\) | \(-0.011\) | Clean |

### Decision
**CONFIRMED** (`confirmed=True`) on the **n=1 fixed** re-fit. Production free-\(n\)
run was sky-sensitive (\(\Delta q_{\rm sky}\sim 0.10\)); locking \(n=1\) kills that
trade. Do not overwrite the production panel; use `20181112A_n1.png` for the
confirmed geometry.

---

## 20190102C — CONFIRMED (protocol-sky re-fit)

**Panel:** `outputs/panels/20190102C.png`

### Host character
Faint, barely resolved: \(m=21.67\), \(R_e=3.9\) px (\(1.02''\)),
\(R_e/\mathrm{FWHM}\approx 0.66\), SNR \(\sim 1.3\).
\(q=0.417\pm 0.034\), PA \(=+4.0^\circ\), \(n=0.50\) (**at GALFIT floor**).
Lotz+2008: **late**. Mag vs LS DR10: \(\Delta m_{\rm ref-m}=-0.04\) (excellent).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q_{\rm AP-G}=-0.022\) (\(q_{\rm AP}=0.395\)) | Agrees on elongated \(q\) |
| Sky \(\pm1\sigma\) | \(q_+=0.327\), \(q_-=0.493\) | \(\lvert\Delta q_{\rm sky}\rvert\approx 0.09\); **sky-minus runaway**: \(n_-\!=\!6\), \(R_{e,-}\!=\!20\) px |
| Isophotes | data \(q(1R_e)\approx 0.63\), model \(\approx 0.69\); `iso_dq_2re`\(\approx -0.038\) | Data≈model; PSF trap vs intrinsic \(q\) |
| Fourier | unreliable; ignore | no leverage |

### Decision
**CONFIRMED** on protocol-sky re-fit (`panel`: `outputs/panels/20190102C_sky.png`). Free sky always degenerates; B/E/F agree at same OOM (span \(\sim 5\times 10^{-4}\); consensus \(6.7\times 10^{-4}\)). Production sky was \(\sim 20\times\) lower.

---

## 20190523A — CONFIRMED (n=1 fixed re-fit)

**Panel:** `outputs/panels/20190523A.png`

### Host character
Very faint / unresolved: \(m=21.72\), \(R_e=2.2\) px (\(0.59''\)),
\(R_e/\mathrm{FWHM}\approx 0.33\), SNR \(\sim 1.5\).
\(q=0.732\pm 0.079\), PA \(=-2.5^\circ\), \(n=0.62\).
Lotz+2008: **merger**. Mag vs LS: \(\Delta m_{\rm ref-m}=-0.18\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(q_{\rm AP}=0.468\), \(\Delta q=-0.264\); \(n_{\rm AP}\approx 7.7\) | Hard disagree — different solution |
| Sky \(\pm1\sigma\) | \(q_+=0.706\), \(q_-=0.577\); \(n\) jumps \(0.5\leftrightarrow 6\) | \(\lvert\Delta q_{\rm sky}\rvert\approx 0.15\) — **sky issue** (typical of these unresolved hosts) |
| Isophotes | data \(q(1R_e)\approx 0.78\), model \(\approx 0.89\); `iso_dq`\(\approx -0.113\) | Disagree |
| Fourier | unreliable (`too_few_annuli,unresolved`); \(n_{\rm annuli}=2\) | **Cannot trust** — too few points |

### Decision
**CONFIRMED** on **n=1** (`outputs/panels/20190523A_n1.png`).

---

## 20190608B — CONFIRMED (spiral residual; Fourier not a \(q\) fix)

**Panel:** `outputs/panels/20190608B.png`

### Host character
Bright, resolved enough: \(m=17.46\), \(R_e=7.4\) px (\(1.95''\)),
\(R_e/\mathrm{FWHM}\approx 1.45\), SNR \(\sim 4.6\).
\(q=0.699\), PA \(=-29.3^\circ\), \(n=3.13\). Lotz+2008: **early**.
Mag vs LS: \(\Delta m_{\rm ref-m}=+0.09\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.006\) | Agrees |
| Sky \(\pm1\sigma\) | \(q_\pm\) within \(0.001\) | Rock stable |
| Isophotes | `iso_dq_2re`\(\approx -0.046\) | Mild gap |
| Fourier | reliable; \(\delta q=-0.097\); \(\psi_2'\approx +51^\circ/R_e\) | Spiral — do not apply \(\delta q\) as \(q\) fix |

### Decision
**CONFIRMED.** Sky+AP lock \(q\); winding \(\psi_2'\) explains dirty residuals / large Fourier \(\delta q\).

---

## 20190611B — REJECTED / DOOMED (do not revisit)

**Panel:** `outputs/panels/20190611B.png`

### Host character
Outside the eventual mag cut (\(m=23.01\)), unresolved:
\(R_e=1.5\) px, \(R_e/\mathrm{FWHM}\approx 0.26\), \(q=0.378\), \(n=0.5\) at floor.

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(q_{\rm AP}=0.088\), \(\Delta q=-0.29\) | Collapse |
| Sky | \(q_+=0.10\), \(q_-=0.78\) (\(\Delta q\sim 0.40\)) | Totally unconstrained |
| Fourier | 2 annuli; \(\delta q\sim 20\pm 30\) | Physically failed |

### Decision
**REJECTED — DOOMED.** Do **not** revisit; no Re-fits. Cannot be fit.

---

## 20190711A — REJECTED (looks like a star)

**Panel:** `outputs/panels/20190711A.png`

### Host character
\(m=18.09\), \(R_e=1.5\) px, \(R_e/\mathrm{FWHM}\approx 0.26\), \(n=6\) at ceiling,
\(R_e\) pinned. Morphologically a **PSF-like point source**, not an extended galaxy.
AstroPhot collapses (\(q_{\rm AP}=0.055\)); local \(\chi^2/\nu\sim 380\), RFF\(+0.22\);
Fourier unusable.

### Decision
**REJECTED.** Note: **looks like a star**, not a galaxy. Not a Re-fits candidate for host geometry.

---

## 20190714A — recommend CONFIRMED (\(n\) floor; barely resolved)

**Panel:** `outputs/panels/20190714A.png`

### Host character
Faint disk, just resolved: \(m=20.24\), \(R_e=4.0\) px (\(1.04''\)),
\(R_e/\mathrm{FWHM}\approx 1.03\), SNR \(\sim 2.2\).
\(q=0.400\pm 0.014\), PA \(=-46.5^\circ\), \(n=0.50\) (**at floor**).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=-0.12\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.025\) (\(q_{\rm AP}=0.376\)) | Agrees on elongated \(q\) |
| Sky \(\pm1\sigma\) | \(q_+=0.397\), \(q_-=0.402\) | \(\lvert\Delta q_{\rm sky}\rvert\approx 0.003\) — stable |
| Isophotes | data \(q(1R_e)\approx 0.586\), model \(\approx 0.590\); `iso_dq_2re`\(\approx -0.010\) | **Data≈model**; gap to intrinsic \(0.40\) is the expected PSF trap |
| Fourier | unreliable (`too_few_annuli`); \(\delta q=+0.033\pm 0.046\) | Consistent with zero; ignore as a correction |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 0.97\); \(\mathrm{RFF}_{2R_e}\approx -0.017\) — clean / noise-like.
- \(n\) at floor: profile shape pinned, but \(q\) is locked by sky + AstroPhot + isophote data–model match.

### Decision
**CONFIRMED.** Caveats: \(n=0.5\) at bound; barely resolved; ignore Fourier \(\delta q\).

---

## 20191001A — recommend CONFIRMED (spiral / merger residual; geometry locked)

**Panel:** `outputs/panels/20191001A.png`

### Host character
Bright, resolved: \(m=18.41\), \(R_e=6.0\) px (\(1.56''\)),
\(R_e/\mathrm{FWHM}\approx 1.41\), SNR \(\sim 4.7\).
\(q=0.542\), PA \(=-41.8^\circ\), \(n=0.51\).
Lotz+2008: **merger**. Mag vs LS: \(\Delta m_{\rm ref-m}=+0.03\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.006\) | Agrees |
| Sky \(\pm1\sigma\) | \(q_\pm\) within \(0.001\) | Rock stable |
| Isophotes | data \(q(1R_e)\approx 0.52\), model \(\approx 0.62\); `iso_dq_2re`\(\approx +0.007\) | Mild data–model gap from structure; IVW \(\Delta q\) tiny |
| Fourier | **reliable**; \(\delta q=-0.036\pm 0.002\); \(\psi_2'\approx -99^\circ/R_e\) | **Winding = spiral/arms** — do not apply \(\delta q\) as \(q\) fix |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\sim 30\); \(\mathrm{RFF}_{2R_e}\approx +0.083\) — dirty, expected for a bright structured / Lotz-merger host.
- Sky + AstroPhot lock \(q\approx 0.54\); dirt is morphology, not wrong ellipse.

### Decision
**CONFIRMED.** Trust \(q\); ignore Fourier \(\delta q\) (\(\psi_2'\) winds).

---

## 20200430A — recommend REJECT (needle \(q\) unconstrained; looks clean but isn’t)

**Panel:** `outputs/panels/20200430A.png`

### Host character
Tiny / unresolved: \(m=21.31\), \(R_e=1.53\) px (\(0.40''\)),
\(R_e/\mathrm{FWHM}\approx 0.30\), SNR \(\sim 0.8\).
\(q=0.043\pm 0.18\), PA \(=+18^\circ\), \(n=2.30\).
Fails eventual \(b/a>0.2\) cut. No Lotz entry (outside 53).

### Why it *looks* perfect (and why that’s a trap)
- Residuals flat, \(\chi^2/\nu\sim 1\), profiles/isophotes can track because there are almost **no independent pixels** inside \(2R_e\) (RFF \(n_{\rm pix}=2\)).
- A clean residual on a PSF-sized blob means “model ≈ data blob,” not “\(q\) is measured.”
- Sky \(\Delta q\sim 0.0001\) and AstroPhot \(\Delta q=+0.010\) (\(q_{\rm AP}=0.052\)) agree on a needle — they can share the same under-constrained solution.
- AstroPhot \(n\approx 7.9\) vs GALFIT \(n=2.3\): profile shape already disagrees.
- Formal \(\sigma_q=0.18\gg q\): the error bar spans nearly the full physical range of \(q\).
- Fourier: **0 annuli** — no vote.

### Decision
**REJECTED.** Way too unresolved / way too small \(q\) (\(0.04\pm 0.18\)).

---

## 20200906A — recommend CONFIRMED

**Panel:** `outputs/panels/20200906A.png`

### Host character
\(m=20.13\), \(R_e=6.0\) px (\(1.56''\)), \(R_e/\mathrm{FWHM}\approx 1.50\),
SNR \(\sim 2.4\). \(q=0.338\), PA \(=+32.9^\circ\), \(n=0.65\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=-0.12\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.016\) (\(q_{\rm AP}=0.322\)) | Agrees |
| Sky \(\pm1\sigma\) | \(q_\pm\) within \(0.001\) | Rock stable |
| Isophotes | data \(q(1R_e)\approx 0.515\), model \(\approx 0.497\); `iso_dq_2re`\(\approx +0.005\) | **Data≈model**; PSF trap vs intrinsic \(0.34\) |
| Fourier | unreliable (`too_few_annuli`); \(\delta q=-0.037\pm 0.015\); \(\psi_2'\approx -8^\circ/R_e\) (flat) | Mild; not a \(q\) veto |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 5.5\); \(\mathrm{RFF}_{2R_e}\approx +0.056\) — some real residual light, modest for a late-type disk.
- Ellipse locked by sky + AP + isophote match.

### Decision
**CONFIRMED.**

---

## 20210320C — recommend CONFIRMED (barely resolved; clean)

**Panel:** `outputs/panels/20210320C.png`

### Host character
\(m=19.80\), \(R_e=3.7\) px (\(0.96''\)), \(R_e/\mathrm{FWHM}\approx 0.94\),
SNR \(\sim 2.6\). \(q=0.702\pm 0.018\), PA \(=-25.2^\circ\), \(n=1.13\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=-0.02\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.013\) | Agrees |
| Sky \(\pm1\sigma\) | \(q_\pm\) identical at \(0.702\) | Rock stable |
| Isophotes | data \(q(1R_e)\approx 0.860\), model \(\approx 0.860\); `iso_dq_2re`\(\approx -0.013\) | **Perfect data–model match**; gap to intrinsic \(0.70\) = PSF trap |
| Fourier | unreliable (`too_few_annuli`); \(\delta q=+0.019\pm 0.047\) | Consistent with zero |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 1.12\); \(\mathrm{RFF}_{2R_e}\approx -0.014\) — clean.

### Decision
**CONFIRMED.** Barely resolved caveat only.

---

## 20210410D — CONFIRMED (n1+sky)

**Panel:** `outputs/panels/20210410D.png`

### Host character
\(m=21.12\), \(R_e=5.4\) px (\(1.43''\)), \(R_e/\mathrm{FWHM}\approx 0.96\),
SNR \(\sim 1.5\). \(q=0.336\pm 0.038\), PA \(=+56.1^\circ\), \(n=6.0\) (**ceiling**).
Lotz+2008: **merger**. Mag vs LS: \(\Delta m_{\rm ref-m}=+0.47\) (poor).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.016\) (\(q_{\rm AP}=0.320\)) | Agrees on \(q\) — but \(n_{\rm AP}\approx 7.7\) also at high-\(n\) |
| Sky \(\pm1\sigma\) | \(q_+=0.211\), \(q_-=0.477\); \(R_{e,-}=14\) px | \(\lvert\Delta q_{\rm sky}\rvert\approx 0.14\) — **sky runaway** |
| Isophotes | data \(q(1R_e)\approx 0.74\), model \(\approx 0.81\); `iso_dq`\(\approx -0.061\) | Mild gap |
| Fourier | unreliable (`too_few_annuli`); ignore | no leverage |

### Residual morphology
- \(\chi^2/\nu\sim 0.9\), RFF\(\approx -0.013\) — looks clean, but parameters are sky-unstable and \(n\) is pinned.

### Decision
**CONFIRMED** on **n1+sky** (`outputs/panels/20210410D_n1_sky.png`). B/E/F span \(\sim 4\times 10^{-4}\) (`agree=False`) but same OOM; consensus \(\sim 4.5\times 10^{-4}\).

---

## 20210807D — CONFIRMED (n=1 re-fit; tames \(R_e\))

**Panel:** `outputs/panels/20210807D.png`

### Host character
Bright extended galaxy: \(m=15.93\), reported \(R_e=92\) px (\(23''\)),
\(R_e/\mathrm{FWHM}\approx 19\), \(q=0.720\), \(n=3.58\). Lotz: merger.
\(1R_e\)/\(2R_e\) overlays look wrong / off-scale because **\(R_e\) itself is nonsense**
(sky-minus also hits the \(R_e=100\) px bound) — not a missing-ellipse plot bug.

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.004\) | \(q\) agrees |
| Sky | \(q\) stable (\(\Delta q\sim 0.01\)); \(R_e\) unbound | Size/sky degeneracy |
| Fourier | winding \(\psi_2'\) | Structure OK; not the issue |
| Mag | \(\Delta m\sim +1.9\) vs PS1 | Consistent with wrong/overlarge \(R_e\) + aperture |

### Decision
**CONFIRMED** on **n=1** (`outputs/panels/20210807D_n1.png`): \(R_e\sim 15\) px (was \(\sim 92\)).

---

## 20211127I — CONFIRMED (n=1 re-fit; tames \(R_e=100\))

**Panel:** `outputs/panels/20211127I.png`

### Host character
Bright giant: \(m=14.25\), \(R_e=100\) px (**exact ceiling**), \(R_e/\mathrm{FWHM}\approx 18\),
\(q=0.968\), PA \(=-30.9^\circ\), \(n=3.14\). Lotz: merger.
Mag vs PS1: \(\Delta m_{\rm ref-m}=+1.24\) — aperture / size issue.

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.002\) | \(q\) agrees (near-round) |
| Sky | \(q\) stable; \(R_e\) pinned at 100 on nominal & sky-minus | **Size at bound** — same failure mode as 20210807D |
| Isophotes | `iso_dq_2re`\(\approx -0.17\) | Large data–model shape tension (partly from wrong size scale) |
| Fourier | reliable; \(\delta q=+0.094\); \(\psi_2'\approx +81^\circ/R_e\) | Structure / winding; \(1R_e\)/\(2R_e\) overlays meaningless if \(R_e\) is fake |

### Decision
**CONFIRMED** on **n=1** (`outputs/panels/20211127I_n1.png`): \(R_e\sim 23\) px (was 100).

---

## 20211203C — CONFIRMED (faint; maybe sky problems)

**Panel:** `outputs/panels/20211203C.png`

### Host character
Outside mag cut (\(m=22.17\)), barely resolved: \(R_e=2.8\) px (\(0.72''\)),
\(R_e/\mathrm{FWHM}\approx 0.69\), SNR \(\lesssim 1\).
\(q=0.437\pm 0.09\), \(n=1.63\). Mag vs LS OK (\(\Delta m\sim -0.06\)).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.041\) | Agrees within large \(\sigma_q\) |
| Sky | \(q_+=0.381\), \(q_-=0.480\); \(n\) jumps | \(\Delta q_{\rm sky}\approx 0.056\) — **maybe sky problems**, not crazy for this faint |
| Isophotes / residuals | modest gap; clean \(\chi^2\)/RFF | Other metrics look OK |
| Fourier | unreliable | ignore |

### Decision
**CONFIRMED for now.** AP + overall metrics carry it; note possible sky issues.

---

## 20211212A — recommend CONFIRMED (spiral residual; geometry locked)

**Panel:** `outputs/panels/20211212A.png`

### Host character
Bright resolved: \(m=16.46\), \(R_e=12.9\) px (\(3.37''\)), \(R_e/\mathrm{FWHM}\approx 2.54\),
SNR \(\sim 9\). \(q=0.829\), PA \(=-55.1^\circ\), \(n=1.09\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=+0.09\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=+0.000\) | Perfect agree |
| Sky \(\pm1\sigma\) | \(q_\pm\) within \(0.0002\) | Rock stable |
| Isophotes | data \(q(1R_e)\approx 0.78\), model \(\approx 0.85\); `iso_dq`\(\approx -0.055\) | Structure-driven gap |
| Fourier | **reliable**; \(\delta q=-0.024\); \(\psi_2'\approx +99^\circ/R_e\) | **Winding = spiral** — not a \(q\) fix |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\sim 34\); \(\mathrm{RFF}_{2R_e}\approx +0.14\) — dirty late-type disk; expected.
- Ellipse locked by sky + AP.

### Decision
**CONFIRMED.** Trust \(q\); ignore Fourier \(\delta q\) / high RFF as morphology.

---

## 20220105A — recommend REJECT (too unresolved / too small \(q\))

**Panel:** `outputs/panels/20220105A.png`

### Host character
\(m=21.47\), \(R_e=2.5\) px (\(0.66''\)), \(R_e/\mathrm{FWHM}\approx 0.43\),
\(q=0.105\pm 0.14\), \(n=1.87\). Fails \(b/a>0.2\) cut. No Lotz.

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(q_{\rm AP}=0.051\), \(\Delta q=-0.054\) | Also needle; \(\sigma_q\) still huge |
| Sky | \(q\sim 0.09\) both sides | Stable needle — shared under-constraint |
| Fourier / RFF | 0 annuli; 7 pix in \(2R_e\) | No leverage |
| Isophotes | observed \(q\sim 0.7\) vs intrinsic \(0.11\) | Classic PSF trap on unresolved needle |

### Decision
**REJECTED.** Way too unresolved / way too small \(q\) (same note as 20200430A).

---

## 20220207C — recommend CONFIRMED (edge-on disk; \(n\) floor)

**Panel:** `outputs/panels/20220207C.png`

### Host character
Bright elongated: \(m=18.12\), \(R_e=23.4\) px (\(5.86''\)), \(R_e/\mathrm{FWHM}\approx 6.2\),
SNR \(\sim 61\). \(q=0.229\), PA \(=+43.7^\circ\), \(n=0.50\) (**floor**).
Lotz+2008: **late**. Mag vs PS1: \(\Delta m_{\rm ref-m}=+0.60\) — aperture/catalog mismatch possible.

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.001\) | Perfect agree |
| Sky \(\pm1\sigma\) | \(q_+=0.227\), \(q_-=0.231\) | Stable |
| Isophotes | data \(q(1R_e)\approx 0.266\), model \(\approx 0.265\); `iso_dq`\(\approx +0.030\) | **Data≈model** |
| Fourier | unreliable (`model_rebuild_poor`); \(\delta q=+0.005\pm 0.010\); flat \(\psi_2'\) | Consistent with zero — ignore |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 1.18\); RFF negative — clean / slightly over-bright model.
- \(n\) pinned at floor: profile shape not free, but \(q\) locked hard.

### Decision
**CONFIRMED.** Strong edge-on geometry; caveat \(n=0.5\) at bound.

---

## 20220307B — recommend CONFIRMED

**Panel:** `outputs/panels/20220307B.png`

### Host character
\(m=19.87\), \(R_e=5.6\) px (\(1.41''\)), \(R_e/\mathrm{FWHM}\approx 1.15\),
SNR \(\sim 49\). \(q=0.632\pm 0.029\), PA \(=+80.9^\circ\), \(n=0.58\).
Lotz+2008: **late**. Mag vs PS1: \(\Delta m_{\rm ref-m}=+0.03\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.005\) | Agrees |
| Sky \(\pm1\sigma\) | \(q_+=0.617\), \(q_-=0.647\) | \(\lvert\Delta q_{\rm sky}\rvert\approx 0.015\) — mild, OK |
| Isophotes | data \(q(1R_e)\approx 0.83\), model \(\approx 0.73\); `iso_dq`\(\approx +0.058\) | Modest data–model gap |
| Fourier | reliable; \(\delta q=-0.054\pm 0.061\); \(\psi_2'\approx -47^\circ/R_e\) (\(\sim 2\sigma\)) | Consistent with zero \(\delta q\); mild winding |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 1.14\); \(\mathrm{RFF}_{2R_e}\approx +0.006\) — clean.

### Decision
**CONFIRMED.**

---

## 20220310F — recommend CONFIRMED (barely resolved; clean)

**Panel:** `outputs/panels/20220310F.png`

### Host character
\(m=20.83\), \(R_e=4.4\) px (\(1.14''\)), \(R_e/\mathrm{FWHM}\approx 0.84\),
SNR \(\sim 2.1\). \(q=0.687\pm 0.027\), PA \(=-61.1^\circ\), \(n=1.50\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=+0.08\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.010\) | Agrees |
| Sky \(\pm1\sigma\) | \(q_+=0.684\), \(q_-=0.692\) | Stable |
| Isophotes | data \(q(1R_e)\approx 0.828\), model \(\approx 0.835\); `iso_dq`\(\approx -0.018\) | **Data≈model**; PSF trap vs \(0.69\) |
| Fourier | unreliable (`too_few_annuli`); \(\delta q\approx 0\); \(\psi_2'\approx -95^\circ/R_e\) | Ignore (few annuli) |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 0.82\); \(\mathrm{RFF}_{2R_e}\approx -0.035\) — clean.

### Decision
**CONFIRMED.** Barely-resolved caveat; geometry locked.

---

## 20220319D — recommend CONFIRMED (\(n\) ceiling; structured)

**Panel:** `outputs/panels/20220319D.png`

### Host character
Bright resolved: \(m=16.08\), \(R_e=26.6\) px (\(6.66''\)), \(R_e/\mathrm{FWHM}\approx 6.2\),
SNR \(\sim 500\). \(q=0.801\pm 0.007\), PA \(=+60.6^\circ\), \(n=6.0\) (**ceiling**).
Lotz+2008: **early**. Mag vs PS1: \(\Delta m_{\rm ref-m}=+0.12\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=+0.002\); \(n_{\rm AP}\approx 8\) | \(q\) agrees; both want high \(n\) |
| Sky \(\pm1\sigma\) | \(q_+=0.791\), \(q_-=0.814\); \(n\) pinned at 6 | \(\lvert\Delta q_{\rm sky}\rvert\approx 0.013\) — \(q\) OK |
| Isophotes | data \(q(1R_e)\approx 0.77\), model \(\approx 0.81\); `iso_dq`\(\approx -0.040\) | Mild gap |
| Fourier | reliable; \(\delta q=-0.047\); \(\psi_2'\approx -43^\circ/R_e\) | Winding structure — do not treat \(\delta q\) as \(q\) fix |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\sim 3.7\); \(\mathrm{RFF}_{2R_e}\approx +0.10\) — real residual light (flag-level); expected for bright structured host with \(n\) pinned.

### Decision
**CONFIRMED.** \(q\) locked by AP+sky; caveat \(n\) at ceiling.

---

## 20220418A — recommend CONFIRMED (\(n\) floor; clean)

**Panel:** `outputs/panels/20220418A.png`

### Host character
\(m=20.88\), \(R_e=5.1\) px (\(1.33''\)), \(R_e/\mathrm{FWHM}\approx 1.00\),
SNR \(\sim 2.0\). \(q=0.757\pm 0.021\), PA \(=+35.6^\circ\), \(n=0.50\) (**floor**).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=-0.13\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.002\) | Agrees |
| Sky \(\pm1\sigma\) | \(q_+=0.752\), \(q_-=0.758\) | Stable |
| Isophotes | data \(q(1R_e)\approx 0.74\), model \(\approx 0.82\); `iso_dq`\(\approx -0.004\) | IVW OK; mild local gap |
| Fourier | reliable; \(\delta q=+0.049\pm 0.043\); \(\psi_2'\approx -40^\circ/R_e\) | Consistent with zero \(\delta q\) |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 1.11\); \(\mathrm{RFF}_{2R_e}\approx -0.031\) — clean.

### Decision
**CONFIRMED.** \(n\)-floor caveat; \(q\) locked.

---

## 20220501C — CONFIRMED (n=1)

**Panel:** `outputs/panels/20220501C.png`

### Host character
\(m=21.80\), \(R_e=2.4\) px (\(0.63''\)), \(R_e/\mathrm{FWHM}\approx 0.47\),
SNR \(\lesssim 1\). \(q=0.523\pm 0.21\), \(n=1.31\). Lotz: late.

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(q_{\rm AP}=0.225\), \(\Delta q=-0.30\); \(n_{\rm AP}\approx 7\) | **Hard disagree** |
| Sky | \(q_+=0.54\), \(q_-=0.35\); \(n\) jumps \(1\leftrightarrow 5\) | \(\Delta q_{\rm sky}\approx 0.17\) — sky |
| Isophotes | `status: unresolved` | No leverage |
| Fourier | few annuli / unresolved | ignore |

### Decision
**CONFIRMED** on **n=1** (`outputs/panels/20220501C_n1.png`).

---

## 20220509G — recommend REJECT (catastrophic residuals; photutils collapse)

**Panel:** `outputs/panels/20220509G.png`

### Host character
Bright: \(m=16.66\), \(R_e=12.6\) px (\(3.31''\)), \(R_e/\mathrm{FWHM}\approx 1.98\),
\(q=0.380\), \(n=0.50\) (**floor**). Lotz: **merger**. Mag vs LS OK (\(\Delta m\sim -0.08\)).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.003\) | \(q\) agrees — but see residuals |
| Sky | \(q_\pm\) within \(0.001\) | Rock stable on a bad model |
| Isophotes | data \(q(1R_e)\approx 0.063\), model \(\approx 0.46\); `iso_dq`\(\approx -0.14\) | **Photutils collapse** — data≪model |
| Fourier | reliable; \(\delta q=-0.047\); \(\psi_2'\approx +15^\circ/R_e\) | Structure present; not the main issue |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\sim 1000\); \(\mathrm{RFF}_{2R_e}\approx +0.48\) — **fit is destroyed**, not “structured but OK.”
- Single-Sérsic cannot describe this host; \(q\) agreement alone is not enough.

### Decision
**REJECTED — WEIRD, investigate in detail.** Centroid/residual catastrophe; photutils collapse. Park `Re-fits/20220509G/`.

---

## 20220717A — recommend REJECT (needle \(q\); AP disagree; mag mismatch)

**Panel:** `outputs/panels/20220717A.png`

### Host character
\(m=18.29\), \(R_e=3.6\) px (\(0.93''\)), \(R_e/\mathrm{FWHM}\approx 0.83\),
\(q=0.071\pm 0.09\), \(n=0.66\). Fails \(b/a>0.2\). No Lotz.
Mag vs PS1: \(\Delta m_{\rm ref-m}=+3.12\), sep \(1.3''\) — wrong/confused counterpart or aperture.

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(q_{\rm AP}=0.200\), \(\Delta q=+0.13\) | Disagrees |
| Sky | \(q\sim 0.07\) both sides | Stable needle |
| Fourier | 0 annuli | no leverage |
| Isophotes | data \(q(1R_e)\approx 0.72\) vs intrinsic \(0.07\) | PSF trap + unconstrained thinness |

### Decision
**REJECTED.** Too unresolved / way too small \(q\).

---

## 20220725A — WEIRD (circle by eye, ellipse in fit)

**Panel:** `outputs/panels/20220725A.png`

### Host character
\(m=17.76\), \(R_e=6.7\) px (\(1.75''\)), \(R_e/\mathrm{FWHM}\approx 1.08\),
SNR \(\sim 3.8\). \(q=0.596\), PA \(=-60.2^\circ\), \(n=0.96\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=+0.05\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.002\) | Agrees |
| Sky \(\pm1\sigma\) | \(q_\pm\) within \(0.001\) | Rock stable |
| Isophotes | data \(q(1R_e)\approx 0.60\), model \(\approx 0.74\); `iso_dq`\(\approx -0.059\) | Structure-driven gap |
| Fourier | reliable; \(\delta q=-0.114\); \(\psi_2'\approx +24^\circ/R_e\) | **Winding = spiral** — not a \(q\) fix |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\sim 47\); \(\mathrm{RFF}_{2R_e}\approx +0.061\) — dirty late-type; expected.
- Ellipse locked by sky + AP.

### Decision
**REJECTED — WEIRD, investigate.** Data look circular; model forces \(q\approx 0.60\)
with a strong quadrupole residual. PSF is round (\(q_{\rm PSF}\approx 0.95\)) — not a
PSF-ellipticity leak. Major/minor data profiles nearly coincide.

---

## 20220825A — recommend CONFIRMED (mild structure; \(q\) locked)

**Panel:** `outputs/panels/20220825A.png`

### Host character
\(m=19.98\), \(R_e=5.4\) px (\(1.36''\)), \(R_e/\mathrm{FWHM}\approx 1.42\),
SNR \(\sim 40\). \(q=0.617\pm 0.050\), PA \(=+44.3^\circ\), \(n=0.68\).
Lotz+2008: **late**. Mag vs PS1: \(\Delta m_{\rm ref-m}=+0.07\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.006\) | Agrees |
| Sky \(\pm1\sigma\) | \(q_+=0.607\), \(q_-=0.621\) | Mild / OK |
| Isophotes | data \(q(1R_e)\approx 0.87\), model \(\approx 0.72\); `iso_dq`\(\approx +0.16\) | Data rounder than model — watch visually |
| Fourier | reliable; \(\delta q=-0.12\pm 0.10\); \(\psi_2'\approx +142^\circ/R_e\) | Winding structure; \(\delta q\) consistent with 0 within err |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 1.09\); \(\mathrm{RFF}_{2R_e}\approx -0.008\) — clean.

### Decision
**CONFIRMED for now.** Looks OK; **too low SNR** (visual / panel). Revisit if SNR cut tightens.

---

## 20220912A — recommend CONFIRMED (mild sky; iso gap)

**Panel:** `outputs/panels/20220912A.png`

### Host character
\(m=19.52\), \(R_e=5.2\) px (\(1.31''\)), \(R_e/\mathrm{FWHM}\approx 1.26\),
SNR \(\sim 53\). \(q=0.634\pm 0.027\), PA \(=+36.3^\circ\), \(n=0.86\).
Lotz+2008: **late**. Mag vs PS1: \(\Delta m_{\rm ref-m}=+0.04\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.004\) | Agrees |
| Sky \(\pm1\sigma\) | \(q_+=0.599\), \(q_-=0.659\); \(n\)/\(R_e\) move | \(\lvert\Delta q_{\rm sky}\rvert\approx 0.035\) — mild |
| Isophotes | data \(q(1R_e)\approx 0.58\), model \(\approx 0.76\); `iso_dq`\(\approx -0.13\) | Data flatter — watch |
| Fourier | reliable; \(\delta q=-0.054\pm 0.052\); \(\psi_2'\approx +103^\circ/R_e\) | Winding; \(\delta q\sim 0\) within err |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 1.16\); \(\mathrm{RFF}_{2R_e}\approx +0.001\) — clean.

### Decision
**CONFIRMED.**

---

## 20220914A — recommend CONFIRMED (barely resolved; clean)

**Panel:** `outputs/panels/20220914A.png`

### Host character
\(m=19.65\), \(R_e=4.8\) px (\(1.25''\)), \(R_e/\mathrm{FWHM}\approx 0.67\),
SNR \(\sim 3.5\). \(q=0.509\pm 0.011\), PA \(=-63.6^\circ\), \(n=2.55\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=-0.04\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.012\) | Agrees |
| Sky \(\pm1\sigma\) | \(q_+=0.505\), \(q_-=0.513\) | Stable |
| Isophotes | data \(q(1R_e)\approx 0.82\), model \(\approx 0.88\); `iso_dq`\(\approx -0.041\) | Data≈model; PSF trap vs \(0.51\) |
| Fourier | unreliable (`too_few_annuli,unresolved`); ignore | no leverage |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 0.95\); \(\mathrm{RFF}_{2R_e}\approx -0.005\) — clean.

### Decision
**CONFIRMED.** Barely-resolved caveat; \(q\) locked.

---

## 20220918A — recommend REJECT (too unresolved / too small \(q\))

**Panel:** `outputs/panels/20220918A.png`

### Host character
\(m=23.08\) (outside mag cut), \(R_e=1.5\) px pinned, \(R_e/\mathrm{FWHM}\approx 0.35\),
\(q=0.095\pm 1.52\), \(n=0.5\) at floor. No Lotz. No ref mag.

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(q_{\rm AP}=0.010\), \(\Delta q=-0.085\) | Collapsed needle |
| Sky | \(q\sim 0.10\); \(R_e,n\) frozen at bounds | Stable nonsense |
| Fourier / RFF | 0 annuli; 3 pix in \(2R_e\) | No leverage |

### Decision
**REJECTED.** Barely resolved; too small \(q\) / unconstrained.

---

## 20220920A — recommend CONFIRMED

**Panel:** `outputs/panels/20220920A.png`

### Host character
\(m=18.13\), \(R_e=7.2\) px (\(1.90''\)), \(R_e/\mathrm{FWHM}\approx 1.24\),
SNR \(\sim 6.5\). \(q=0.860\), PA \(=+47.9^\circ\), \(n=1.12\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=-0.10\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.002\) | Agrees |
| Sky \(\pm1\sigma\) | \(q_\pm\) within \(0.001\) | Rock stable |
| Isophotes | data \(q(1R_e)\approx 0.903\), model \(\approx 0.909\); `iso_dq`\(\approx -0.016\) | **Data≈model** |
| Fourier | reliable; \(\delta q=-0.014\pm 0.005\); mild \(\psi_2'\) | Tiny absolute bias |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 1.47\); \(\mathrm{RFF}_{2R_e}\approx +0.005\) — clean.

### Decision
**CONFIRMED.**

---

## 20221012A — recommend CONFIRMED (barely resolved; clean)

**Panel:** `outputs/panels/20221012A.png`

### Host character
\(m=18.96\), \(R_e=5.6\) px (\(1.47''\)), \(R_e/\mathrm{FWHM}\approx 0.99\),
SNR \(\sim 3.9\). \(q=0.533\), PA \(=+10.0^\circ\), \(n=2.01\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=-0.04\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.010\) | Agrees |
| Sky \(\pm1\sigma\) | \(q_\pm\) within \(0.001\) | Rock stable |
| Isophotes | data \(q(1R_e)\approx 0.747\), model \(\approx 0.760\); `iso_dq`\(\approx +0.007\) | **Data≈model** |
| Fourier | unreliable (`too_few_annuli`); \(\delta q=-0.051\pm 0.017\); \(\psi_2'\approx -76^\circ/R_e\) | Mild; ignore as \(q\) fix |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 1.36\); \(\mathrm{RFF}_{2R_e}\approx -0.001\) — clean.

### Decision
**CONFIRMED.** Barely-resolved caveat; \(q\) locked.

---

## 20221101B — recommend CONFIRMED (\(q\) locked; \(R_e\) sky-sensitive)

**Panel:** `outputs/panels/20221101B.png`

### Host character
\(m=18.32\), \(R_e=17.5\) px (\(4.58''\)), \(R_e/\mathrm{FWHM}\approx 4.3\),
SNR \(\sim 88\). \(q=0.752\pm 0.032\), PA \(=-34.1^\circ\), \(n=4.18\).
Lotz+2008: **merger**. Mag vs PS1: \(\Delta m_{\rm ref-m}=+1.10\), sep \(1.3''\) — aperture / counterpart mismatch.

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.008\) | Agrees |
| Sky \(\pm1\sigma\) | \(q_+=0.760\), \(q_-=0.741\); \(R_{e,-}=47\) px | \(q\) stable; **\(R_e\) sky-runaway** (envelope) |
| Isophotes | data \(q(1R_e)\approx 0.88\), model \(\approx 0.79\); `iso_dq`\(\approx +0.049\) | Mild |
| Fourier | reliable; \(\delta q=-0.071\pm 0.045\); \(\psi_2'\approx -52^\circ/R_e\) | Winding — not a \(q\) fix |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 1.16\); \(\mathrm{RFF}_{2R_e}\approx -0.035\) — clean enough.

### Decision
**CONFIRMED** for \(q\). Note: \(R_e\) may be too big; possible \(R_e\)–\(n\) degeneracy (sky-minus blows \(R_e\), \(n\to 6\)).

---

## 20221106A — recommend CONFIRMED

**Panel:** `outputs/panels/20221106A.png`

### Host character
\(m=18.09\), \(R_e=9.4\) px (\(2.46''\)), \(R_e/\mathrm{FWHM}\approx 1.68\),
SNR \(\sim 3.2\). \(q=0.672\), PA \(=+26.6^\circ\), \(n=1.79\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=+0.19\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.004\) | Agrees |
| Sky \(\pm1\sigma\) | \(q_\pm\) within \(0.001\) | Rock stable |
| Isophotes | data \(q(1R_e)\approx 0.780\), model \(\approx 0.775\); `iso_dq`\(\approx -0.029\) | **Data≈model** |
| Fourier | reliable; \(\delta q=-0.010\pm 0.006\); mild \(\psi_2'\) | Tiny |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 1.55\); \(\mathrm{RFF}_{2R_e}\approx +0.003\) — clean.

### Decision
**CONFIRMED.**

---

## 20221116A — recommend REJECT (bullshit \(R_e\); \(n\) ceiling; sky; no isophotes)

**Panel:** `outputs/panels/20221116A.png`

### Host character
\(m=17.01\), \(R_e=74\) px (\(19''\)), \(R_e/\mathrm{FWHM}\approx 18\),
\(q=0.165\pm 0.036\), \(n=6.0\) (**ceiling**). Fails \(b/a>0.2\). No Lotz / no ref mag.

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.011\) | Thin \(q\) agrees |
| Sky | \(q_+=0.097\), \(q_-=0.251\); \(R_{e,-}=100\) | \(\Delta q_{\rm sky}\approx 0.09\); **\(R_e\) at bound** |
| Isophotes | `isophote_fit_failed` | No photutils |
| Fourier | unreliable (`model_rebuild_poor`); huge \(\psi_2'\) | Don't trust |
| RFF | \(\mathrm{RFF}_{2R_e}\approx -0.68\) | Model badly wrong in flux |

### Decision
**REJECTED.** Bullshit \(R_e\); sky-dominated / almost nothing real in the stamp
(\(\sim 62\%\) model flux outside stamp; RFF \(\sim -0.68\)).
**No ref mag:** production `zp_ok=false`, `mag_final_source=unavailable` — catalog ZP/match
never landed (photometry broken by the oversize/\(R_e\) mess), so `ref_mag` is null.

---

## 20221219A — CONFIRMED (n1+sky re-fit)

**Panel:** `outputs/panels/20221219A.png`

### Host character
\(m=21.97\), \(R_e=4.0\) px (\(1.04''\)), \(R_e/\mathrm{FWHM}\approx 0.53\),
\(q=0.427\pm 0.069\), \(n=0.5\) at floor. Mag vs LS: \(\Delta m\sim -0.30\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.027\) | Mild agree |
| Sky | \(q_+=0.380\), \(q_-=0.452\) | \(\Delta q_{\rm sky}\approx 0.047\) |
| Isophotes | data \(0.53\) vs model \(0.72\) at \(1R_e\); `iso_dq`\(\approx -0.17\) | Large gap |
| Fourier | unreliable / unresolved | ignore |

### Decision
**CONFIRMED** on **n1+sky** (outputs/panels/20221219A_n1_sky.png). B/E/F poorly agree (F near 0); consensus used as most physical sky.

---

## 20230124A — recommend CONFIRMED (barely resolved; clean)

**Panel:** `outputs/panels/20230124A.png`

### Host character
\(m=18.58\), \(R_e=3.8\) px (\(1.00''\)), \(R_e/\mathrm{FWHM}\approx 0.58\),
SNR \(\sim 4.4\). \(q=0.903\pm 0.006\), PA \(=+87.6^\circ\), \(n=1.93\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=-0.04\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.004\) | Agrees |
| Sky \(\pm1\sigma\) | \(q_\pm\) within \(0.002\) | Rock stable |
| Isophotes | data \(q(1R_e)\approx 0.982\), model \(\approx 0.967\); `iso_dq`\(\approx +0.011\) | Data≈model; PSF trap vs \(0.90\) |
| Fourier | unreliable (`too_few_annuli,unresolved`); \(\delta q\sim 0\) | ignore |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 1.05\); \(\mathrm{RFF}_{2R_e}\approx -0.001\) — clean.

### Decision
**CONFIRMED.** Barely-resolved caveat; \(q\) locked (near-round).

---

## 20230307A — recommend CONFIRMED (\(n\) ceiling; \(q\) locked)

**Panel:** `outputs/panels/20230307A.png`

### Host character
\(m=18.18\), \(R_e=25.1\) px (\(6.58''\)), \(R_e/\mathrm{FWHM}\approx 3.9\),
SNR \(\sim 4.1\). \(q=0.804\pm 0.011\), PA \(=+13.7^\circ\), \(n=6.0\) (**ceiling**).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=+0.32\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=+0.001\); \(n_{\rm AP}\approx 7.2\) | \(q\) agrees; both high-\(n\) |
| Sky \(\pm1\sigma\) | \(q_+=0.799\), \(q_-=0.807\); \(n\) pinned at 6 | \(q\) stable |
| Isophotes | data \(q(1R_e)\approx 0.73\), model \(\approx 0.81\); `iso_dq`\(\approx -0.034\) | Mild gap |
| Fourier | reliable; \(\delta q\approx 0\); mild \(\psi_2'\) | No geometry bias |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 1.28\); \(\mathrm{RFF}_{2R_e}\approx -0.035\) — clean enough.

### Decision
**CONFIRMED.** Caveat \(n\) at ceiling; \(q\) locked.

---

## 20230526A — REJECTED (star-like; \(R_e\) at floor)

**Panel:** `outputs/panels/20230526A.png`

### Host character
\(m=20.09\), \(R_e=1.5\) px (**floor**), \(R_e/\mathrm{FWHM}\approx 0.29\),
SNR \(\sim 1.3\). \(q=0.709\pm 0.029\), \(n=1.85\). Lotz: late.
Mag vs LS: \(\Delta m_{\rm ref-m}=+0.37\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=+0.011\) | Agrees on \(q\) |
| Sky | \(q\) stable; \(R_e\) frozen at 1.5 | Size not measured |
| Fourier | 2 annuli / unresolved | ignore |
| Isophotes | no \(q(1R_e)\) | Almost no leverage |

### Decision
**REJECTED.** Looks like a star; \(R_e\) at floor / too unresolved for a constrained host geometry.

---

## 20230613A — recommend CONFIRMED (compact disk; \(q\) locked)

**Panel:** `outputs/panels/20230613A.png`

### Host character
\(m=20.31\), \(R_e=4.2\) px (\(1.11''\)), \(R_e/\mathrm{FWHM}\approx 0.83\),
SNR \(\sim 3.0\). \(q=0.838\pm 0.006\), PA \(=-78.2^\circ\), \(n=0.62\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=-0.18\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.001\) | Rock agree |
| Sky \(\pm1\sigma\) | \(q_\pm\) within \(0.001\) | Rock stable |
| Isophotes | data \(q(1R_e)\approx 0.90\), model \(\approx 0.91\); `iso_dq`\(\approx -0.021\) | Mild rounder isophotes vs \(q=0.84\) |
| Fourier | unreliable (`too_few_annuli`); \(\delta q\sim -0.026\) | ignore for geometry |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 1.42\); \(\mathrm{RFF}_{2R_e}\approx +0.004\) — clean overall; mild core dip/peak.
- Profiles track Sérsic well; no spiral/bar signature.

### Decision
**CONFIRMED.** Compact but resolved; \(q\) locked by sky+AP. Caveat: low \(n\), mild central residual, barely above PSF scale.

---

## 20230626A — REJECTED (star-like; too unresolved)

**Panel:** `outputs/panels/20230626A.png`

### Host character
\(m=20.00\), \(R_e=1.8\) px (\(0.48''\)), \(R_e/\mathrm{FWHM}\approx 0.28\),
SNR \(\sim 2.9\). \(q=0.738\pm 0.025\), PA \(=-22.7^\circ\), \(n=6.0\) (**ceiling**).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=+0.10\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.020\); \(n_{\rm AP}\approx 8\) | Mild \(q\) offset; both high-\(n\) |
| Sky \(\pm1\sigma\) | \(q_+=0.722\), \(q_-=0.753\); \(R_e\) swings \(1.6\leftrightarrow 1.9\) | \(q\) soft; size soft |
| Isophotes | no \(q(1R_e)\); data \(q(2R_e)\approx 0.92\) vs GALFIT \(0.74\) | Classic PSF trap |
| Fourier | unreliable (`too_few_annuli,unresolved`) | ignore |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 1.36\); \(\mathrm{RFF}_{2R_e}\approx -0.002\) — clean but unconstrained.

### Decision
**REJECTED.** Looks like a star; too unresolved; \(n\) ceiling; PSF trap on \(q\).

---

## 20230628A — recommend CONFIRMED (disk; \(q\) locked; neighbor)

**Panel:** `outputs/panels/20230628A.png`

### Host character
\(m=19.56\), \(R_e=4.5\) px (\(1.18''\)), \(R_e/\mathrm{FWHM}\approx 0.72\),
SNR \(\sim 3.8\). \(q=0.727\pm 0.008\), PA \(=+82.9^\circ\), \(n=1.04\).
Lotz+2008: **merger** (likely neighbor). Mag vs LS: \(\Delta m_{\rm ref-m}=-0.09\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.004\) | Rock agree |
| Sky \(\pm1\sigma\) | \(q_\pm\) within \(0.005\) | Rock stable |
| Isophotes | data \(q(1R_e)\approx 0.84\), model \(\approx 0.89\); `iso_dq`\(\approx -0.032\) | Mild rounder data vs \(q=0.73\) |
| Fourier | unreliable (`too_few_annuli,unresolved`) | ignore |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 1.02\); \(\mathrm{RFF}_{2R_e}\approx -0.006\) — clean on host.
- Unmodeled neighbor NE of host (drives Lotz merger).

### Decision
**CONFIRMED.** Compact disk; \(q\) locked by sky+AP. Caveats: neighbor in stamp; Lotz=merger; barely resolved.

---

## 20230708A — recommend REJECT (too unresolved / too small \(q\))

**Panel:** `outputs/panels/20230708A.png`

### Host character
\(m=21.89\) (outside mag cut), \(R_e=1.5\) px (**floor**), \(R_e/\mathrm{FWHM}\approx 0.25\),
SNR \(\sim 0.6\). \(q=0.051\pm 0.16\), \(n=6.0\) (**ceiling**). Not in 53-cut.
Mag vs LS: \(\Delta m_{\rm ref-m}=-0.76\) (**dmag flag**).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(q_{\rm AP}\approx 0.01\); \(\Delta q=-0.040\) | Both nonsense needles |
| Sky | \(R_e\), \(n\) frozen; \(q\sim 0.05\) | No leverage |
| Isophotes | `isophote_fit_failed` | No leverage |
| Fourier | 0 annuli / unresolved | ignore |

### Decision
**REJECTED.** Way too unresolved / way too small \(q\) (same class as 20200430A etc.).

---

## 20230712A — recommend REJECT (too unresolved / too small \(q\))

**Panel:** `outputs/panels/20230712A.png`

### Host character
\(m=22.17\) (outside mag cut), \(R_e=3.2\) px (\(0.83''\)), \(R_e/\mathrm{FWHM}\approx 0.45\),
SNR \(\sim 0.8\). \(q=0.117\pm 0.47\), \(n=1.07\). Not in 53-cut.
Mag vs LS: \(\Delta m_{\rm ref-m}=-0.28\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(q_{\rm AP}\approx 0.02\); \(\Delta q=-0.094\) | Both needles |
| Sky \(\pm1\sigma\) | \(q_+=0.089\), \(q_-=0.117\) | Soft; still tiny |
| Isophotes | data \(q(1R_e)\approx 0.43\) vs model \(\approx 0.79\); GALFIT \(q=0.12\) | PSF / unconstrained mess |
| Fourier | 1 annulus / unresolved | ignore |

### Decision
**REJECTED.** Too unresolved / too small \(q\) with huge \(q\) error — same class as prior needles.

---

## 20230902A — CONFIRMED (barely resolved; PSF-trap caveat)

**Panel:** `outputs/panels/20230902A.png`

### Host character
\(m=21.54\), \(R_e=1.6\) px (\(0.42''\)), \(R_e/\mathrm{FWHM}\approx 0.34\),
SNR \(\sim 1.0\). \(q=0.628\pm 0.069\), PA \(=+24.5^\circ\), \(n=0.88\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=-0.05\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.015\) | Soft agree |
| Sky \(\pm1\sigma\) | \(q_+\approx 0.614\), \(q_-\approx 0.642\); \(n\) swings | \(q\) soft; \(n\) soft |
| Isophotes | no \(q(1R_e)\); data \(q(2R_e)\approx 0.91\) vs GALFIT \(0.63\) | Classic PSF trap |
| Fourier | unreliable (`too_few_annuli,unresolved`) | ignore |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 1.71\); \(\mathrm{RFF}_{2R_e}\approx -0.007\) — clean but unconstrained.

### Decision
**CONFIRMED** (user override). Barely resolved; \(q\) soft-locked by sky+AP. Keep PSF-trap / low-resolution caveat.

---

## 20230907D — recommend CONFIRMED (resolved disk; \(q\) locked)

**Panel:** `outputs/panels/20230907D.png`

### Host character
\(m=19.41\), \(R_e=4.6\) px (\(1.22''\)), \(R_e/\mathrm{FWHM}\approx 1.18\),
SNR \(\sim 2.7\). \(q=0.855\pm 0.015\), PA \(=+82.0^\circ\), \(n=1.58\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=+0.09\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=+0.002\) | Rock agree |
| Sky \(\pm1\sigma\) | \(q_\pm\) within \(0.003\) | Rock stable |
| Isophotes | data \(q(1R_e)\approx 0.916\), model \(\approx 0.914\); `iso_dq`\(\approx +0.006\) | Excellent match |
| Fourier | unreliable (`too_few_annuli`); \(\delta q\sim 0\) | ignore; no geometry bias |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 1.17\); \(\mathrm{RFF}_{2R_e}\approx -0.006\) — clean.

### Decision
**CONFIRMED.** Well resolved for this cohort; \(q\) locked by sky+AP+isophotes.

---

## 20230930A — WEIRD (investigate)

**Panel:** `outputs/panels/20230930A.png`

### Host character
\(m=8.90\) (**not physical** — Phase-2 ZP failure; \(\Delta m_{\rm ref-m}=-8.57\) vs PS1),
\(R_e=80.8\) px (\(21.2''\)), \(R_e/\mathrm{FWHM}\approx 18.7\), SNR \(\sim 184\).
\(q=0.841\pm 0.013\), PA \(=+70.9^\circ\), \(n=4.70\). Lotz+2008: **merger**.

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.003\) | \(q\) agrees |
| Sky \(\pm1\sigma\) | \(q_+=0.734\), \(q_-=0.892\); \(R_e\) \(10\leftrightarrow 100\); \(n\) \(1.6\leftrightarrow 3.5\); \(\Delta m_{\rm sky}\sim 1.7\) | **Sky runaway** |
| Isophotes | `iso_dq_2re`\(\approx -0.21\); no \(q(1R_e)\) | Strong geometry mismatch |
| Fourier | reliable; \(\delta q\sim -0.03\); wild outer \(\delta q(a)\) | Outer structure / stamp issues |

### Residual morphology
- Core bipolar residual; profiles under-predict central light.
- Documented in `FIT_VERIFICATION_CHECKS.md`: ZP failure — mag should not be trusted.

### Decision
**REJECTED — WEIRD, investigate.** Visual/model mismatch; bullshit \(R_e\); ZP failure; sky runaway.

---

## 20231020B — recommend CONFIRMED (barely resolved; PSF-trap caveat)

**Panel:** `outputs/panels/20231020B.png`

### Host character
\(m=21.75\), \(R_e=2.2\) px (\(0.59''\)), \(R_e/\mathrm{FWHM}\approx 0.56\),
SNR \(\sim 1.3\). \(q=0.591\pm 0.028\), PA \(=+1.6^\circ\), \(n=0.67\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=+0.04\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.017\) | Soft agree |
| Sky \(\pm1\sigma\) | \(q_\pm\) within \(0.013\) | Stable |
| Isophotes | data \(q(1R_e)\approx 0.73\), model \(\approx 0.82\); GALFIT \(0.59\) | Mild PSF trap |
| Fourier | unreliable (`too_few_annuli,unresolved`) | ignore |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 1.35\); \(\mathrm{RFF}_{2R_e}\approx -0.007\) — clean; mild core ring.

### Decision
**CONFIRMED.** Barely resolved; \(q\) soft-locked by sky+AP. Keep PSF-trap / low-\(R_e\) caveat (same class as 20230124A / 20230902A).

---

## 20231120A — recommend CONFIRMED (edge-on disk; spiral residual; \(q\) locked)

**Panel:** `outputs/panels/20231120A.png`

### Host character
\(m=15.47\), \(R_e=30.0\) px (\(7.87''\)), \(R_e/\mathrm{FWHM}\approx 5.4\),
SNR \(\sim 15.5\). \(q=0.251\pm 0.000\), PA \(=+28.4^\circ\), \(n=1.91\).
Lotz+2008: **early**. Mag vs LS: \(\Delta m_{\rm ref-m}=-0.03\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.002\) | Rock agree |
| Sky \(\pm1\sigma\) | \(q\) frozen at \(0.251\) | Rock stable |
| Isophotes | data \(q(1R_e)\approx 0.30\), model \(\approx 0.31\); outer \(q\to 0.25\) | Matches thin disk |
| Fourier | **reliable**; \(\delta q=-0.034\); \(\psi_2'\approx -89^\circ/R_e\) | **Winding = spiral** — not a \(q\) fix |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\sim 91\); \(\mathrm{RFF}_{2R_e}\approx +0.13\) — dirty; S/spiral residual expected for bright structured host.
- Same class as 20190608B / 20211212A.

### Decision
**CONFIRMED.** Trust \(q\) (sky+AP lock); ignore high RFF / Fourier \(\delta q\) as morphology. Caveat: \(q\) near science cut (\(0.2\)).

---

## 20231123B — recommend CONFIRMED (inclined disk; mild structure; \(q\) locked)

**Panel:** `outputs/panels/20231123B.png`

### Host character
\(m=18.59\), \(R_e=12.0\) px (\(3.15''\)), \(R_e/\mathrm{FWHM}\approx 2.0\),
SNR \(\sim 5.3\). \(q=0.378\pm 0.002\), PA \(=+21.5^\circ\), \(n=0.70\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=-0.19\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.004\) | Rock agree |
| Sky \(\pm1\sigma\) | \(q_\pm\) within \(0.002\) | Rock stable |
| Isophotes | data \(q(1R_e)\approx 0.45\), model \(\approx 0.47\); `iso_dq`\(\approx 0\) | Good match (PSF-rounder core) |
| Fourier | **reliable**; \(\delta q=-0.041\); \(\psi_2'\approx -89^\circ/R_e\) | Winding / central \(m=2\) — not a \(q\) fix |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 3.17\); \(\mathrm{RFF}_{2R_e}\approx +0.034\) — mild structure; acceptable.

### Decision
**CONFIRMED.** Inclined late-type; \(q\) locked by sky+AP.

---

## 20231220A — recommend REJECT (sky runaway; small \(q\); PSF trap)

**Panel:** `outputs/panels/20231220A.png`

### Host character
\(m=20.22\), \(R_e=4.1\) px (\(1.08''\)), \(R_e/\mathrm{FWHM}\approx 0.55\),
SNR \(\sim 2.4\). \(q=0.214\pm 0.039\) (near cut), PA \(=-27.2^\circ\), \(n=0.57\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=-0.15\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(q_{\rm AP}\approx 0.17\); \(\Delta q=-0.044\) | Disagrees / tinier |
| Sky \(\pm1\sigma\) | \(q_+=0.135\), \(q_-=0.275\); \(n\) \(0.5\leftrightarrow 0.79\) | **Sky runaway** (\(\lvert\Delta q\rvert\sim 0.08\)) |
| Isophotes | data \(q(1R_e)\approx 0.72\) vs GALFIT \(0.21\) | Classic PSF trap |
| Fourier | unreliable (`too_few_annuli,unresolved`) | ignore |

### Residual morphology
- Clean \(\chi^2\)/RFF — unconstrained geometry, not a pass.

### Decision
**REJECTED** for now. Sky runaway plus \(q\) too small (near cut); PSF trap.

---

## 20231226A — recommend CONFIRMED (disk; \(q\) locked; isophote dip)

**Panel:** `outputs/panels/20231226A.png`

### Host character
\(m=18.54\), \(R_e=8.0\) px (\(2.10''\)), \(R_e/\mathrm{FWHM}\approx 1.54\),
SNR \(\sim 2.7\). \(q=0.758\pm 0.007\), PA \(=-17.4^\circ\), \(n=0.65\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=+0.31\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=+0.000\) | Rock agree |
| Sky \(\pm1\sigma\) | \(q_\pm\) within \(0.004\) | Rock stable |
| Isophotes | data \(q(1R_e)\approx 0.53\), model \(\approx 0.81\); `iso_dq`\(\approx -0.19\) | Strong radial \(q\) dip |
| Fourier | **reliable**; \(\delta q=-0.035\); \(\psi_2'\approx +81^\circ/R_e\) | Winding / structure — not a global \(q\) fix |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 2.68\); \(\mathrm{RFF}_{2R_e}\approx +0.063\) — mild structure; central dipole.

### Decision
**CONFIRMED.** Sky+AP lock \(q\); treat isophote dip / Fourier winding as morphology. Caveat: large `iso_dq`.

---

## 20240119A — recommend CONFIRMED (barely resolved; \(n\) floor; \(q\) soft-locked)

**Panel:** `outputs/panels/20240119A.png`

### Host character
\(m=21.44\), \(R_e=5.2\) px (\(1.36''\)), \(R_e/\mathrm{FWHM}\approx 0.69\),
SNR \(\sim 1.2\). \(q=0.750\pm 0.052\), PA \(=-75.9^\circ\), \(n=0.50\) (**floor**).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=-0.18\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=+0.010\) | Agrees |
| Sky \(\pm1\sigma\) | \(q_+\approx 0.727\), \(q_-\approx 0.745\); \(n\) pinned | Soft; \(n\) frozen |
| Isophotes | data \(q(1R_e)\approx 0.38\) (noisy), model \(\approx 0.86\) | Unstable inner isophotes |
| Fourier | unreliable (`too_few_annuli,unresolved`) | ignore |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 0.93\); \(\mathrm{RFF}_{2R_e}\approx -0.094\) — slightly over-bright model; clumpy core.

### Decision
**CONFIRMED.** Good fit visually; sky response solid. Caveats: **SNR very low** (\(\sim 1.2\)); \(n\) at floor; noisy isophotes.

---

## 20240201A — recommend CONFIRMED (edge-on; \(q\) locked; below \(q\) cut)

**Panel:** `outputs/panels/20240201A.png`

### Host character
\(m=16.94\), \(R_e=18.5\) px (\(4.85''\)), \(R_e/\mathrm{FWHM}\approx 4.0\),
SNR \(\sim 6.2\). \(q=0.179\pm 0.001\) (**below** \(b/a>0.2\) cut → not in 53),
PA \(=+62.2^\circ\), \(n=0.61\). Mag vs LS: \(\Delta m_{\rm ref-m}=+0.00\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.004\) | Rock agree |
| Sky \(\pm1\sigma\) | \(q\) frozen \(\approx 0.179\) | Rock stable |
| Isophotes | data \(q(1R_e)\approx 0.22\), model \(\approx 0.24\) | Good match for thin disk |
| Fourier | **reliable**; \(\delta q=-0.019\); mild \(\psi_2'\) | Small; not a \(q\) fix |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 6.9\); \(\mathrm{RFF}_{2R_e}\approx +0.059\) — structured residual (dust/arms) expected for bright edge-on.

### Decision
**CONFIRMED** for geometry. Caveat: outside science cut (\(q<0.2\)).

---

## 20240208A — CONFIRMED (user override; barely resolved)

**Panel:** `outputs/panels/20240208A.png`

### Host character
\(m=21.62\), \(R_e=2.8\) px (\(0.73''\)), \(R_e/\mathrm{FWHM}\approx 0.60\),
SNR \(\sim 1.0\). \(q=0.277\pm 0.037\), PA \(=+66.2^\circ\), \(n=1.45\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=+0.10\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(q_{\rm AP}\approx 0.23\); \(\Delta q=-0.044\) | Mild offset / tinier |
| Sky \(\pm1\sigma\) | \(q_+\approx 0.265\), \(q_-\approx 0.288\); \(n\) swings | Soft but bounded |
| Isophotes | data \(q(1R_e)\approx 0.69\) vs GALFIT \(0.28\) | PSF trap caveat |
| Fourier | unreliable (`too_few_annuli,unresolved`) | ignore |

### Decision
**CONFIRMED** (user override). Barely resolved; keep PSF-trap + AP-offset + low-SNR caveats.

---

## 20240210A — recommend CONFIRMED (bright disk; spiral residual; \(q\) locked)

**Panel:** `outputs/panels/20240210A.png`

### Host character
\(m=14.50\), \(R_e=27.7\) px (\(7.26''\)), \(R_e/\mathrm{FWHM}\approx 5.8\),
SNR \(\sim 16\). \(q=0.594\pm 0.001\), PA \(=+57.3^\circ\), \(n=1.00\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=+0.40\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.000\) | Rock agree |
| Sky \(\pm1\sigma\) | \(q\) frozen at \(0.594\) | Rock stable |
| Isophotes | data \(q(1R_e)\approx 0.65\), model \(\approx 0.61\); `iso_dq`\(\approx 0\) | Good |
| Fourier | **reliable**; \(\delta q=-0.013\); \(\psi_2'\approx -80^\circ/R_e\) | **Winding = spiral** — not a \(q\) fix |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\sim 180\); \(\mathrm{RFF}_{2R_e}\approx +0.19\) — dirty bright late-type (same class as 20211212A / 20231120A).

### Decision
**CONFIRMED.** Trust \(q\); ignore high RFF / spiral residual for geometry.

---

## 20240213A — recommend CONFIRMED (inclined disk; \(q\) locked; near cut)

**Panel:** `outputs/panels/20240213A.png`

### Host character
\(m=18.35\), \(R_e=9.6\) px (\(2.52''\)), \(R_e/\mathrm{FWHM}\approx 1.67\),
SNR \(\sim 6.1\). \(q=0.238\pm 0.001\) (near cut), PA \(=-40.9^\circ\), \(n=0.51\) (near floor).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=-0.07\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.008\) | Soft agree |
| Sky \(\pm1\sigma\) | \(q_\pm\) within \(0.002\) | Rock stable |
| Isophotes | data \(q(1R_e)\approx 0.35\), model \(\approx 0.38\) | Reasonable for thin disk |
| Fourier | **reliable**; \(\delta q=-0.010\); \(\psi_2'\approx +99^\circ/R_e\) | Winding / structure — not a \(q\) fix |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 27.8\); \(\mathrm{RFF}_{2R_e}\approx +0.073\) — structured residual expected.

### Decision
**CONFIRMED.** \(q\) locked; caveats near \(q\) cut and low \(n\).

---

## 20240215A — recommend CONFIRMED (barely resolved; \(q\) soft-locked; PSF-trap caveat)

**Panel:** `outputs/panels/20240215A.png`

### Host character
\(m=20.06\), \(R_e=3.1\) px (\(0.81''\)), \(R_e/\mathrm{FWHM}\approx 0.45\),
SNR \(\sim 2.7\). \(q=0.482\pm 0.016\), PA \(=+20.7^\circ\), \(n=2.03\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=-0.06\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.007\) | Agrees |
| Sky \(\pm1\sigma\) | \(q_\pm\) within \(0.018\); \(n\) \(1.6\leftrightarrow 2.7\) | \(q\) soft-stable; \(n\) soft |
| Isophotes | data \(q(1R_e)\approx 0.70\), model \(\approx 0.80\); GALFIT \(0.48\) | PSF trap |
| Fourier | unreliable (`too_few_annuli,unresolved`) | ignore |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 2.18\); \(\mathrm{RFF}_{2R_e}\approx +0.007\) — clean.

### Decision
**CONFIRMED.** Barely resolved; \(q\) soft-locked by sky+AP. Keep PSF-trap / low-\(R_e\) caveat.

---

## 20240229A — CONFIRMED (user override; unresolved + \(n\) ceiling)

**Panel:** `outputs/panels/20240229A.png`

### Host character
\(m=20.84\), \(R_e=2.1\) px (\(0.56''\)), \(R_e/\mathrm{FWHM}\approx 0.32\),
SNR \(\sim 1.9\). \(q=0.626\pm 0.069\), PA \(=-78.3^\circ\), \(n=6.0\) (**ceiling**).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=+0.10\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.011\); \(n_{\rm AP}\approx 8\) | Soft agree; both high-\(n\) |
| Sky \(\pm1\sigma\) | \(q_+=0.567\), \(q_-=0.683\); \(R_e\) \(1.7\leftrightarrow 2.7\) | Soft / near-runaway |
| Isophotes | data \(q(1R_e)\approx 0.92\) vs GALFIT \(0.63\) | Classic PSF trap |
| Fourier | unreliable (`too_few_annuli,unresolved`) | ignore |

### Decision
**CONFIRMED** (user override). Keep unresolved / \(n\)-ceiling / sky-soft / PSF-trap caveats.

---

## 20240304A — recommend CONFIRMED (compact disk; \(q\) locked)

**Panel:** `outputs/panels/20240304A.png`

### Host character
\(m=20.81\), \(R_e=3.6\) px (\(0.93''\)), \(R_e/\mathrm{FWHM}\approx 0.80\),
SNR \(\sim 1.3\). \(q=0.665\pm 0.038\), PA \(=+20.0^\circ\), \(n=1.24\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=-0.00\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.014\) | Soft agree |
| Sky \(\pm1\sigma\) | \(q\) frozen \(\approx 0.664\); \(n\) soft | \(q\) rock stable |
| Isophotes | data \(q(1R_e)\approx 0.69\), model \(\approx 0.82\); mild gap | Acceptable |
| Fourier | unreliable (`too_few_annuli,unresolved`) | ignore |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 1.06\); \(\mathrm{RFF}_{2R_e}\approx -0.029\) — clean.

### Decision
**CONFIRMED.** Compact but resolved; \(q\) locked by sky. Caveat: low SNR.

---

## 20240310A — recommend REJECT (too small \(q\); PSF trap; AP disagree)

**Panel:** `outputs/panels/20240310A.png`

### Host character
\(m=19.67\), \(R_e=3.5\) px (\(0.92''\)), \(R_e/\mathrm{FWHM}\approx 0.70\),
SNR \(\sim 1.6\). \(q=0.103\pm 0.028\) (**far below** cut → not in 53),
PA \(=-13.6^\circ\), \(n=0.76\). Mag vs LS: \(\Delta m_{\rm ref-m}=+0.38\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(q_{\rm AP}\approx 0.042\); \(\Delta q=-0.061\) | Needle; disagrees |
| Sky \(\pm1\sigma\) | \(q\sim 0.10\)–\(0.11\) | Stable but nonsense \(q\) |
| Isophotes | data \(q(1R_e)\approx 0.58\) vs GALFIT \(0.10\) | Classic PSF trap |
| Fourier | 1 annulus / unresolved | ignore |

### Decision
**REJECTED.** Very needle-like; PSF trap; AP tinier. Same reason as other needles.

---

## 20240318A — CONFIRMED (inclined disk; \(q\) locked) — **last host**

**Panel:** `outputs/panels/20240318A.png`

### Host character
\(m=18.47\), \(R_e=9.6\) px (\(2.53''\)), \(R_e/\mathrm{FWHM}\approx 1.80\),
SNR \(\sim 3.4\). \(q=0.427\pm 0.003\), PA \(=-6.9^\circ\), \(n=0.70\).
Lotz+2008: **late**. Mag vs LS: \(\Delta m_{\rm ref-m}=+0.33\).

### Geometry
| test | result | reading |
|---|---|---|
| AstroPhot | \(\Delta q=-0.004\) | Rock agree |
| Sky \(\pm1\sigma\) | \(q_\pm\) within \(0.002\) | Rock stable |
| Isophotes | data \(q(1R_e)\approx 0.50\), model \(\approx 0.53\) | Good match |
| Fourier | **reliable**; \(\delta q=-0.026\); mild \(\psi_2'\) | Mild structure — not a \(q\) fix |

### Residual morphology
- \(\chi^2/\nu|_{2R_e,{\rm corr}}\approx 1.27\); \(\mathrm{RFF}_{2R_e}\approx +0.000\) — clean overall; mild central residual.

### Decision
**CONFIRMED** (user re-confirm after brief reject). Inclined late-type; \(q\) locked. Completes the 64-host triage.
