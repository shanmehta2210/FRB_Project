# median expAB_r mag-bin plots — SDSS provenance and cut variants

## Is `expAB_r` exactly what SDSS reports?

**Yes, in this pipeline.** The v2 catalog stores the PhotoObj field **unchanged**:

1. **SQL query** (`scripts/build_sdss_null_catalog_v2.py`): `p.expAB_r AS expAB_r`
2. **Catalog build**: `expab = pd.to_numeric(df["expAB_r"], errors="coerce")` → written to column `expAB_r` with **no rescaling, inversion, or Hubble mapping**
3. **These plots**: `median(expAB_r)` per 0.5 mag bin — **no transform** on the axis ratio

The only filters on `expAB_r` before it enters the CSV are:

| Stage | Filter | Applied by |
|-------|--------|------------|
| SDSS SQL | `expAB_r > 0 AND expAB_r <= 1` | Catalog build query |
| CSV row drop | finite, `0 ≤ expAB_r ≤ 1` | `build_catalog()` |
| Plot pools (below) | optional u-r, lnL, strict b/a | Audit script only |

We do **not** re-fit profiles or recompute b/a in this repo.

### What SDSS means by `expAB_r`

From the SDSS imaging pipeline ([model magnitude algorithm](https://cas.sdss.org/dr7/en/help/docs/algorithm.asp?key=mag_model)):

- For each object, SDSS fits a **2D exponential profile**  
  \(I(r) = I_0 \exp(-1.68\, r/r_e)\) (plus a separate de Vaucouleurs fit).
- Each fit has free **axis ratio** and position angle; the exponential axis ratio is stored as **`expAB_r`** in `PhotoObj` (“Exponential fit a/b”, [DR7 schema](https://cas.sdss.org/dr7/en/help/browser/description.asp?n=PhotoObj&t=V)).
- **`lnLExp_r`** and **`lnLDeV_r`** are the log-likelihoods of the exponential vs de Vaucouleurs fits.

SDSS reports these for **all** galaxies meeting basic photo flags; it does **not** apply our u-r, lnL-winner, or `b/a > 0.2` cuts. Those are **null-pool selections** applied after download.

**Note:** PhotoObj stores axis ratios with a practical floor near **0.05** (~22% of v1 sample at exactly 0.05; see `SDSS_BA_FLOOR_MIN` in `null_catalog_utils.py`). That is an SDSS storage/quantization limit, not something we impose in the median plots.

---

## Three plot variants

All use **median raw `expAB_r`** vs `modelMag_r` (0.5 mag bins, min N = 30). No inclination math.

| Pool | Cuts | CSV | PNG |
|------|------|-----|-----|
| **Full catalog** | Finite `expAB_r ∈ [0,1]` only | `ba_per_mag_bin.csv` | `ba_mag_joint_panel.png` |
| **Morphology** | `u-r < 2.3`, `lnLExp_r > lnLDeV_r`, all b/a | `ba_per_mag_bin_ur_lnl.csv` | `ba_mag_joint_panel_ur_lnl.png` |
| **Strict production** | Above + `expAB_r > 0.2` | `ba_per_mag_bin_strict.csv` | `ba_mag_joint_panel_strict.png` |

Code: `run_ba_mag_bin_plots()` in `scripts/audit_sdss_v2_cosi_formal.py`.

---

## How to read the three curves together

| Mag range | Full catalog | ur + lnL | ur + lnL + strict |
|-----------|--------------|----------|-------------------|
| Bright (mag ≲ 20) | Median ~0.6 | Similar | Similar (most pass strict) |
| mag 21–22 | Falls sharply | Falls, slightly higher | Higher still (low-b/a removed) |
| mag 22–23 | Median ~0.10 (floor pileup) | Median ~0.13 | Median ~0.40–0.45 |

- **Full → morphology:** lnL exp-winner removes bulge-dominated systems; median faint `expAB_r` changes modestly.
- **Morphology → strict:** Removing `expAB_r ≤ 0.2` **drops the floor-dominated majority** at faint mags; the strict curve shows the axis ratio of **survivors only**, not the raw SDSS column for all exp-winners.

For CDF work, the **strict** b/a plot is closest to the sample that enters production `cos(i)` mapping; the **full** plot shows what SDSS reports before any science cuts.

---

## Why strict median `expAB_r` looks like median cos(i)

**There is no bug.** The Hubble map (q₀ = 0.2) is **strictly increasing** in `q` for `q > q₀`:

\[
\cos i = \sqrt{\frac{q^2 - q_0^2}{1 - q_0^2}}
\]

For any **monotonic increasing** function \(f\), **median(f(q)) = f(median(q))**. So per mag bin:

\[
\text{median}(\cos i) = \text{Hubble}(\text{median}(\texttt{expAB\_r}))
\]

The curves have the **same shape**; cos(i) sits **below** b/a by ~0.02–0.06 (larger offset when b/a is smaller). Example mag 22.0–22.5:

| Quantity | Value |
|----------|-------|
| median `expAB_r` | 0.449 |
| median cos(i) | 0.410 |
| Difference | 0.039 |

They are **not** numerically identical — but on a 0–1 axis the tracks look parallel.

**Diagnostic:** `ba_cosi_strict_overlay.png` overlays both medians; lower panel plots `median(b/a) − median(cos i)`. Data: `ba_cosi_strict_comparison.csv`.

**Note:** `mean(cos i) ≠ Hubble(mean(b/a))` in general; only the **median** commutes with this transform.

---

## Reproduction

```bash
python scripts/audit_sdss_v2_cosi_formal.py
# or only the b/a block:
python -c "
import sys; sys.path.insert(0,'scripts')
from pathlib import Path
from audit_sdss_v2_cosi_formal import read_sdss_null_catalog, run_ba_mag_bin_plots, OUT_DIR
from pipeline_null_plot_utils import DEFAULT_SDSS_V2
from null_catalog_utils import read_sdss_null_catalog
df = read_sdss_null_catalog(DEFAULT_SDSS_V2)
run_ba_mag_bin_plots(df, OUT_DIR)
"
```

See also [`THREE_PLOT_DEEP_ANALYSIS.md`](THREE_PLOT_DEEP_ANALYSIS.md) for faint-bin floor statistics.
