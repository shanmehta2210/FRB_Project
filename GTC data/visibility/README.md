# GTC visibility

Deterministic five-gate filter for inclination-safe GTC imaging. Reads
`master_frb_localization.csv` at repo root.

## Usage

```bash
conda activate gtc_visibility

# Single night
python "GTC data/visibility/gtc_visibility_batch.py" --date 2026-06-24

# Month scan (nightly CSVs + rollups)
python "GTC data/visibility/gtc_visibility_batch.py" --date 2026-06-24 --end-date 2026-07-24 --quiet
```

## Output

| Directory | Contents |
|-----------|----------|
| `nightly/` | `gtc_visibility_YYYY-MM-DD.csv` per night |
| `summaries/` | `gtc_availability_by_frb_*`, `gtc_availability_by_night_*`, `gtc_visibility_long_*` |

## Gates

1. δ ≥ −36.25° (GTC dome lower limit 25°)
2. X ≤ 1.5 (Bernstein & Jarvis 2002 PSF smearing)
3. Moon separation ≥ 30° (Lotz et al. 2004 outer isophotes)
4. Astronomical dark (Sun ≤ −18°)
5. h ≤ 72° (GTC upper-shutter vignetting)

Default: ≥30 min contiguous window where gates 2–5 hold inside dark time.
