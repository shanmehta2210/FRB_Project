# Host-property audit: protocol (round 2)

Purpose: every value that will enter the paper's host table must be verified against the primary paper, and every gap must be filled from the most recent properly published source. Accuracy over speed. Read `SWEEP_PROTOCOL.md` first for file conventions, fetching rules and the controlled quantity list.

## Inputs

- `audit_seed_picks.csv`: the value and source currently selected for each Table-1 column of the 51-host working sample. Your job is to confirm or correct every cell for your hosts.
- `host_properties_long.csv`: every value collected in round 1, with `source_short`, `method`, `location`.
- `papers/`: PDFs and text extracts already downloaded. Reuse them; download anything else you need into the same folder.
- `../../../master_frb_localization.csv`: repo coordinates, $z$, DM, NE2001 for all hosts (including the six extra hosts listed below).

## Output (one set per group, do not overwrite round-1 files)

- `audit_group_<N>.csv`, UTF-8, columns exactly:
  `frb,quantity,current_value,current_source,verdict,verified_value,err_plus,err_minus,unit,method,source_short,arxiv_id,doi,location,notes`
  - `verdict` is one of `confirmed`, `corrected`, `new`, `unavailable`, `conflict`.
  - `confirmed`: you opened the cited paper and saw the value at `location` (table/section/figure). Quote the location precisely (e.g. `Table 3, row 20210320C`).
  - `corrected`: the cited paper says something different; put the paper's value in `verified_value`.
  - `new`: the cell was empty and you found a published value.
  - `unavailable`: you searched (see search rules) and nothing published exists. Say in `notes` which queries/papers you checked.
  - `conflict`: two or more refereed papers disagree beyond their quoted errors; put all values in separate rows with `verdict=conflict` and explain in notes which one you recommend and why (newest refereed paper wins unless it explicitly adopts an older value or is a compilation).
- `audit_group_<N>_notes.md`: per-host prose summary and anything odd.
- `audit_group_<N>.bib`: BibTeX for any paper not already in `host_properties.bib` (key `firstauthorYYYY`, with doi and eprint).

## Quantities to verify for every host

`z` (and `z_type`: spec/photo), `dm_obs`, `dm_mw` (NE2001, state YMW16 separately), `logmstar`, `sfr` (state timescale and indicator in `method`), `m_r` (literature host $r$ magnitude, say which survey/filter), `offset_kpc`, `inclination_deg` or `axis_ratio` if any paper quotes one, `re_kpc`, `host_name`, `galaxy_type`. Also record `tm_gyr`, `metallicity_logz`, `av`, `u_r_color`, `rm_rad_m2`, `tau_scatt_ms` if the paper has them.

Rules for mass/SFR: record the SED code, IMF and timescale. Prospector 0-100 Myr SFR is the target column; Bagpipes/CIGALE SFRs acceptable when nothing else exists; H$\alpha$ SFR recorded with `method=halpha`. Do not accept values copied by a compilation paper without checking the original (cite the original, note the compilation).

## Search rules for missing values

For any host with a missing `logmstar`, `sfr`, `m_r` or `z` do an independent search, in this order, and log each query in the notes:
1. ADS full-text / abstract search via `https://ui.adsabs.harvard.edu/search/q=%22FRB%2020240119A%22` style queries (also try the short form, e.g. `FRB 20240119`, and the host name if any).
2. arXiv listing search: `https://arxiv.org/a/?searchtype=all&query=%22FRB+20240119A%22` (use `curl.exe -A "Mozilla/5.0" -L`).
3. Google Scholar is not scriptable; use `https://www.semanticscholar.org/search?q=` or `https://api.semanticscholar.org/graph/v1/paper/search?query=FRB%2020240119A&fields=title,year,externalIds` instead.
4. Check the TNS entry and the CHIME/DSA/CRAFT catalogue pages if relevant.
Prefer refereed papers; arXiv-only papers may be used but mark `notes=arXiv only` and record the arXiv id.

## Hosts

Working sample (51): see `audit_seed_picks.csv`.
Six additional hosts that are in the 57-host pre-cut sample but were not swept in round 1: 20200430A, 20220105A, 20220717A, 20230708A, 20240310A, 20240201A. These need a full round-1-style sweep (every published property), written in the same `audit_group_<N>.csv` format with `verdict=new`.

## Known open questions to settle

- 20231120A redshift: Sharma et al. 2024 (LRIS spectrum) and Verdi et al. 2026 give 0.0368; Connor et al. 2025 and Hussaini et al. 2025 give 0.0700; NED lists a photometric 0.0664. Find which value the most recent refereed DSA-110 papers adopt, whether any paper comments on the discrepancy, and whether an independent spectrum exists (e.g. SDSS/DESI/LAMOST/6dF/2MRS for the host 2MFGC 07440 at RA 143.98, Dec +73.28). Report cz, the instrument, and the paper.
- 20231220A, 20240119A, 20240213A, 20240215A, 20240229A: no published mass/SFR in round 1. Search for newer DSA-110 host papers (2025-2026), the DSA-110 host catalogue, and any other follow-up.
- 20220501C, 20230902A, 20240208A, 20240304A, 20240318A, 20240210A: CRAFT hosts with missing SED values; search for the CRAFT host-galaxy catalogue papers (Marnoch et al., Gordon et al. 2025, Shannon et al. 2025 companion papers, Glowacki et al. 2025) and any 2026 follow-ups.
- 20220912A: only an H$\alpha$ lower limit and a radio SFR exist; check Hewitt et al. 2024, Ravi et al. 2023, Sheikh et al. 2024 and any 2025-2026 host papers for an SED SFR.
- 20190523A: SED mass exists (Ravi et al. 2019); SFR is an upper limit. Confirm there is no later SED SFR.

## Windows notes

PowerShell: use `;` not `&&`. Avoid printing non-ASCII in `python -c`; write to files with `encoding="utf-8"`. Do not edit any file outside `Publication/materials/host_properties/`. Do not run git commands that change state.
