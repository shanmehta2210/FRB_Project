# SDSS b/a quantization audit (catalog vs binning)

## Conclusion

**The horizontal stacking at b/a ≈ n×0.05 is from SDSS PhotoObj storage/fitting, not from hexbin or scatter binning.**

- **Scatter panel:** raw `best_model_ba_r` values — no y-binning.
- **Hexbin y spacing:** ≈0.0175 in b/a over [0,1] at `gridsize=50` — **not** 0.05.
- **Exactly** b/a = n×0.05: 22.5% of galaxies.
- **Within 0.01 of some n×0.05:** 58.9% — explains **multiple** visible bands.
- Dominant floor: **22.3%** at exactly 0.05.


The **scatter panel uses raw (mag, b/a) points — no binning.** If bands appear there, they are catalog values.

## Full v1 SDSS sample

| Metric | Value |
|--------|-------|
| Rows scanned | 457,190 |
| Finite b/a (0–1) | 457,190 |
| On n×0.05 grid | 102,723 (22.47%) |
| Exactly b/a = 0.05 | 101,979 (22.31%) |

### Counts at each 0.05 rung (full sample)

| b/a | N | % of finite |
|-----|---|-------------|
| 0.05 | 101,979 | 22.31% |
| 0.10 | 75 | 0.02% |
| 0.15 | 114 | 0.02% |
| 0.20 | 88 | 0.02% |
| 0.25 | 70 | 0.02% |
| 0.30 | 49 | 0.01% |
| 0.35 | 52 | 0.01% |
| 0.40 | 49 | 0.01% |
| 0.45 | 41 | 0.01% |
| 0.50 | 34 | 0.01% |
| 0.55 | 39 | 0.01% |
| 0.60 | 25 | 0.01% |
| 0.65 | 19 | 0.00% |
| 0.70 | 23 | 0.01% |
| 0.75 | 16 | 0.00% |
| 0.80 | 15 | 0.00% |
| 0.85 | 10 | 0.00% |
| 0.90 | 8 | 0.00% |
| 0.95 | 4 | 0.00% |
| 1.00 | 8 | 0.00% |

### Off-grid values (sample)

- n_off_grid_sampled: 50000
- off_grid_min: 0.05000202
- off_grid_max: 0.999997
- residual_abs_max: 0.02499979999999999
- residual_abs_median: 0.010818949999999938

## 15k scatter subsample (same draw as plots)

- u−r < 3.5: N=15,000, on 0.05 grid 23.1%, at 0.05 floor 22.9%
- u−r < 2.2: N=15,000, on 0.05 grid 24.1%, at 0.05 floor 23.9%
- u−r < 1.5: N=15,000, on 0.05 grid 26.9%, at 0.05 floor 26.8%

## Per u−r color cut (full population passing cut)

- u−r < 3.5: N=326,568, on grid 23.0%, at 0.05 22.8%
- u−r < 2.2: N=198,447, on grid 23.8%, at 0.05 23.6%
- u−r < 1.5: N=104,921, on grid 26.9%, at 0.05 26.7%

### Galaxies within 0.01 of each rung (explains multi-band look)

| b/a rung | N (near) | % of sample |
|----------|----------|-------------|
| 0.05 | 110,433 | 24.15% |
| 0.10 | 17,680 | 3.87% |
| 0.15 | 11,812 | 2.58% |
| 0.20 | 10,570 | 2.31% |
| 0.25 | 10,640 | 2.33% |
| 0.30 | 10,139 | 2.22% |
| 0.35 | 10,599 | 2.32% |
| 0.40 | 9,788 | 2.14% |
| 0.45 | 9,593 | 2.10% |
| 0.50 | 9,359 | 2.05% |
| 0.55 | 8,814 | 1.93% |
| 0.60 | 8,309 | 1.82% |
| 0.65 | 8,106 | 1.77% |
| 0.70 | 7,174 | 1.57% |
| 0.75 | 6,367 | 1.39% |
| 0.80 | 5,560 | 1.22% |
| 0.85 | 4,195 | 0.92% |
| 0.90 | 3,138 | 0.69% |
| 0.95 | 2,081 | 0.46% |
| 1.00 | 4,805 | 1.05% |

## Implication for SDSS email

1. `best_model_ba_r` shows a strong **floor at 0.05** (~22% exactly 0.05).
2. Most other galaxies sit **within ~0.01–0.025** of some n×0.05 value — not continuous face-on→edge-on.
3. Ask SDSS whether exp/deV axis ratios are quantized, rounded, or bounded in PhotoObj.
4. State that **our plots do not impose 0.05 spacing**; scatter is unbinned; hex cell height ≈0.0175.

## v2 formal follow-up

Mag-binned median `expAB_r` and floor fractions on the 1.9M v2 catalog: [`plots/plots_null/v2/sdss_audit/formal/THREE_PLOT_DEEP_ANALYSIS.md`](../../../v2/sdss_audit/formal/THREE_PLOT_DEEP_ANALYSIS.md) §2.
