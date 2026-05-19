# Null catalog v0 (archived)

**Archived:** 2026-05-19. Do not use for new science; retained for reproduction of old figures only.

## Files

| File | Notes |
|------|--------|
| `LS_catalog.csv` | Legacy DR10 TAP `TOP 200000` without spatial `ORDER BY` — landed in RA ~284–288°, Dec ~44–46° (not COSMOS). Column `petroMag_r` is Tractor model r mag, not Petrosian. |
| `SDSS_catalog.csv` | COSMOS-only export (`new_SDSS_DR16_cosmos`); **no RA column**; `petroMag_r` is Petrosian; `rmag` is model r. |
| `SDSS_catalogue.txt` | Parent export for `convert_catalog.py`. |

## Known issues (see repo `tasks.md`)

- Magnitude system mismatch at plot time (`petroMag_r` cut on both surveys).
- Zero spatial overlap between SDSS (COSMOS) and Legacy v0 patch.

## Replacement

Use repo-root v1 catalogs:

- `LS_catalog_v1_allsky_modelmr.csv`
- `SDSS_catalog_v1_allsky_modelmr.csv`

Scripts not yet migrated to v1 defaults: `generate_sigma_plots.py`, `generate_all_plots.py`, `generate_multiband_cdf_null_plot.py` — pass `--sdss-csv` explicitly or point to this archive.
