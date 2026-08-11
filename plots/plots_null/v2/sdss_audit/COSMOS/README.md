# COSMOS HST vs SDSS b/a audit

See [CATALOG_DECISIONS.md](CATALOG_DECISIONS.md) and [ZURICH_CATALOG_NOTES.md](ZURICH_CATALOG_NOTES.md).

## Run order

From repo root:

```bash
python scripts/build_cosmos_hst_zurich_catalog.py
python scripts/build_cosmos_sdss_catalog.py --no-color-cut
python scripts/plot_cosmos_null_cdfs.py --clean
```

CDF plots use **disk-only** pools by default:
- **HST:** GIM2D `TYPE = 2` (disk/spiral)
- **SDSS:** `lnLExp_r > lnLDeV_r` (exponential profile wins)

Optional pass-1 b/a vs mag panels:

```bash
python scripts/plot_cosmos_ba_mag_audit.py
```

## Catalogs

| File | Description |
|------|-------------|
| `cosmos_hst_zurich_disk_entire.csv` | HST GIM2D TYPE=2, no b/a floor |
| `cosmos_hst_zurich_disk_strict.csv` | HST disk, b/a > 0.2 |
| `cosmos_sdss_dr17_nocolor_disk_entire.csv` | SDSS exp-winner disks in ACS box |
| `cosmos_sdss_dr17_nocolor_disk_strict.csv` | SDSS disk, expAB_r > 0.2 |

SDSS entire pool includes all footprint sources from the SQL query (type=3, clean=1, mode=1) with finite `modelMag_r` and `expAB_r`; no lnL exp-winner or colour cuts.

## CDF plots (`plots/cdfs/`)

| Folder | Mag cut | Files |
|--------|---------|-------|
| `mag20/` | modelMag_r / ACS_MAG_AUTO <= 20 | `entire.png`, `strict.png` |
| `mag21/` | <= 21 | `entire.png`, `strict.png` |
| `mag22/` | <= 22 | `entire.png`, `strict.png` |

Summary: `plots/cdfs/ks_summary.csv`

**entire** — all sources with valid b/a in footprint (no b/a > 0.2 cut).  
**strict** — same pools with b/a > 0.2.
