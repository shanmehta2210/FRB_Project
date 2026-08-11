# PSF normalization (PSFEx → GALFIT → verification)

## What PSFEx writes

Pipeline Phase 1 (`SExtractor + PSFEx/run_psf_pipeline.py`):

- `PHOTFLUX_KEY = FLUX_APER(1)` with first aperture diameter **4 px**
  (`phot_apertures_px` in config).
- `PSF_SAMPLING = 1` → native pixel scale (`PSF_SAMP=1` in headers).
- Checkimage `PROTOTYPES` → `proto_image.fits` (25×25). With
  `PSFVAR_DEGREES=0` this stamp is **bit-identical** to the reconstructed
  constant term of `image.psf`.

**Normalization:** PSFEx scales the model so flux inside the photometric
aperture ≈ 1 — **not** so \(\sum_{\rm stamp}=1\), and **not**
\(\sum \mathrm{PSF}^2=1\).

Cohort check (64 hosts): median stamp sum ~3.4 (range ~2.2–5.7), strongly
correlated with FWHM; median sum inside diam=4 px ≈ **1.035**. Sum of squares
~0.15 — not \(L^2\)-normalized.

## What GALFIT does

Peng et al. (2002/2010) / FAQ: for **convolution**, GALFIT **renormalizes the
PSF image to unit total flux** before the FFT. A stamp that sums to 4 does
**not** multiply the Sérsic model by 4; amplitude is absorbed into mag.

Caveats:

- Object type `psf` (pure star) assumes the stamp contains ~100% of the light
  for magnitude; truncated wings bias that mag.
- Item `E)` fine-sampling factor must match `PSF_SAMP` (here both **1**).
- Wrong **shape** / FWHM still biases \(n\), \(R_e\), \(q\) (mismatched PSF).

## What verification does

| routine | behavior |
|---|---|
| `vercommon.convolve_psf` | divides kernel by `sum` before FFT convolve |
| Analytic Sérsic rebuild | normalized to total flux analytically, then convolved |
| `psf_second_moments` | flux-weighted moments of the stamp (shape check) |
| AstroPhot cross-fit | uses the same `proto_image.fits` |

So structural checks and the visual model are consistent with unit-sum
convolution kernels even when the on-disk stamp is aperture-normalized.

## Spotcheck artifacts

- `Re-fits/psf_spotcheck.png` — visual stamps for compact / weird hosts.
- `Re-fits/psf_spotcheck_n1sky_sky.{csv,json}` — sky + n1_sky fit summary for
  the unresolved spotcheck set.

## Trusting tiny \(R_e\)

Literature (e.g. Davari et al. 2014; Peng FAQ): forward-modeling can recover
sizes \(\lesssim\mathrm{FWHM}\) when PSF + sky are good.

For **inclination / \(q\)**: when \(R_e/\mathrm{FWHM}\sim 0.2\) (Re floor,
\(n\to\) ceiling), the source is star-like — \(q\) is PSF-trapped. Confirm only
with unresolved caveats, or reject needles below the \(b/a\) cut.

Tiny ellipses on `data/σ` or asinh panels are **intrinsic** \(1\text{–}2R_e\),
not the observed light size ([`VISUAL_PANELS.md`](VISUAL_PANELS.md)).

## Practical checklist

1. `PSF_SAMP` / GALFIT `E)` both 1.  
2. Stamp looks round, FWHM ~5–8 px, no crazy ellipticity.  
3. Do **not** panic if \(\sum\neq 1\); check aperture sum (~4 px diam) ≈ 1.  
4. If shapes look wrong, rebuild from `image.psf` and compare to `proto_image`.  
5. Sub-FWHM \(q\): treat as soft even if the residual looks “clean”.
