# Visual panels

Implementation: `checks/visual.py`. Interpretation of bottom plots and header
metrics: [`FIT_VERIFICATION_CHECKS.md` §4.9](FIT_VERIFICATION_CHECKS.md)
(physics unchanged; **top-row stretch updated** — see below).

## Canonical stretch (verification panel)

| panel | colormap | scaling |
|---|---|---|
| **data** | grayscale | sky-subtracted flux → asinh over data’s **[1, 99]%** window; soft factor 10 |
| **model** | grayscale | **same** lo/hi window as data (fair comparison) |
| **residual** | `RdBu_r` | residual / σ, clipped **±5** |

```text
x = clip((img - lo) / (hi - lo), 0, 1)
display = asinh(10·x) / asinh(10)
```

Black ellipses on images = **intrinsic** GALFIT \(1R_e\) and \(2R_e\)
(semi-major), not free isophotes. For \(R_e\ll\mathrm{FWHM}\) they look tiny
vs the PSF-smeared blob — expected ([`PSF_NORMALIZATION.md`](PSF_NORMALIZATION.md)).

Sandbox quick panels use asinh for data/model but **linear** residual clip
([`SANDBOX.md`](SANDBOX.md)) — different on purpose.

## Layout (2×3)

| | |
|---|---|
| top | data · model · residual (/σ) |
| bottom | major/minor \(\mu(a)\) · isophotal \(q(a)\) · Fourier \(\delta q(a)\) |

Empty bottom plots mean that check’s `.npz` was missing when visual ran —
re-run `isophote` / `fourier` (or the full suite) then `visual --force`.

## Paths

| path | meaning |
|---|---|
| `outputs/per_host/<FRB>/panel.png` | suite output for production host |
| `outputs/panels/<FRB>.png` | published copy (sync after regen) |
| `outputs/panels/<FRB>_<leg>.png` | confirmed alternate leg |
| `Re-fits/<FRB>/<leg>/panel.png` | re-fit suite panel |
| `Re-fits/<FRB>/panel_*.png` | reject-grid published copies |
| `outputs/plots/contact_sheet.png` | residual contact sheet (confirmed-leg / PSF-only where gated) |
| `outputs/confirmed_fit_panels.pptx` | one confirmed panel per slide |

## Regenerating panels

```bash
# All 64 production visuals
python run_verification.py --checks visual --force --jobs 4 --no-aggregate

# Sync to outputs/panels/
python -c "import glob,os,shutil,vercommon as vc
for p in glob.glob(os.path.join(vc.PER_HOST_ROOT,'*','panel.png')):
    frb=os.path.basename(os.path.dirname(p))
    shutil.copy2(p, os.path.join(vc.OUT_ROOT,'panels',f'{frb}.png'))"
```

Re-fit legs: re-run visual via `run_refit` or `vis.run(load_host_from_dir(...))`
into that leg, then re-copy confirmed `outputs/panels/<FRB>_<leg>.png`.

## PowerPoint deck

```bash
python build_confirmed_pptx.py
```

Builds `outputs/confirmed_fit_panels.pptx`:

1. Title slide  
2. Population overview (key numbers from `population_summary.json`)  
3. One slide each: `population_diagnostics.png`, `mag_leakage.png`,
   `dq_comparison.png`, `contact_sheet.png`  
4. Section divider → one slide per **in_53 ∩ confirmed** host (alternate panel
   if notes cite it)  
5. Final slide(s) for **in_53 ∩ rejected** (red banner + note + production panel)

Panel / plot PNGs on disk are **not** recomputed. The builder embeds
full-resolution PNGs at **200 DPI** (letterboxed when aspect differs from the
host-panel slide), flattening RGBA→RGB only in memory so PowerPoint does not
soft-render transparency.

## Reading traps (short)

- GALFIT \(q\) (blue dotted) is **intrinsic**; isophotal \(q(a)\) is
  PSF-rounded — gap at low \(R_e/\mathrm{FWHM}\) is the PSF trap, not
  automatic reject.
- Localized \(\chi^2/\nu|_{2R_e}\) ignores bright junk outside \(2R_e\).
- Lotz class in the header comes from production `statmorph_results.json`
  when present.
