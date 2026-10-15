# Host-property literature sweep: protocol

Target: the 51 FRB hosts in `seed_51.csv` (the working sample of the inclination paper). Goal is a complete per-host property table mirroring Table 1 of Bhardwaj, Lee & Ji (2024, Nature 634, 1065) and an extended table with every published quantity we can find.

## Output files (one set per group)

- `group_<N>.csv` long format, one row per (host, quantity), UTF-8, columns exactly:
  `frb,quantity,value,err_plus,err_minus,unit,method,source_short,arxiv_id,doi,location,notes`
  - `quantity` from the controlled list below (add new ones if needed, lower_snake_case).
  - `value` numeric where possible; strings for names/classes. Upper limits as `<value` in `value` with `method=upper_limit`.
  - `location` = table number / section / figure where the value appears in the paper.
  - `method` = how the value was derived in the paper (e.g. `prospector_sed`, `halpha`, `cigale`, `NE2001`, `YMW16`, `spec_z`, `photo_z`, `PATH`, `galfit`, `sextractor`).
- `group_<N>_notes.md` per-host prose: host name, which papers were read, conflicts between papers, anything odd.
- `group_<N>.bib` BibTeX for every paper used (key = `firstauthorYYYY`, e.g. `gordon2023`). Include doi and eprint.

## Controlled quantity list

host_name, host_ra_deg, host_dec_deg, z, z_type, DM_obs, DM_MW (give model in method), DM_halo, DM_host (as published; method = the IGM model used), DM_exgal, tau_scatt_ms (give frequency in notes), RM_rad_m2, fluence_Jy_ms, snr, width_ms, m_r, m_g, m_i, m_z, M_r_abs, logMstar (log10 Msun), SFR (Msun/yr; say timescale and indicator in method), log_sSFR, age_mass_weighted_Gyr, tm_Gyr, metallicity_logZ, AV, u_r_color, g_r_color, Re_arcsec, Re_kpc, offset_arcsec, offset_kpc, offset_normalized, PATH_posterior, morphology, galaxy_type (e.g. SF/transition/quiescent), inclination_deg (if any paper quotes one), axis_ratio, position_angle_deg, environment (cluster/group/field), telescope_localization, repeater, luminosity_distance_Mpc.

## Rules

- Read what the repo already holds first: `seed_51.csv`, `../notes/HOST_COORDINATE_CHECK.md`, `../../../production67_lit_astropath_audit.csv`, `../../../Archive/notes/repeater_audit_notes.md`, `../../../repeater_localizations_README.md`, `../../../CHIME/catalog/host_association_papers.md` if present.
- Primary sources: the host-association paper cited in `seed_51.csv` plus any host-property compilations (Gordon et al. 2023 ApJ 952 122; Gordon et al. 2025 ApJ 993 119; Sharma et al. 2024 Nature 635 61; Law et al. 2024 ApJ 967 29; Hussaini et al. 2025 ApJL 993 L27; Heintz et al. 2020 ApJ 903 152; Bhandari et al. 2022 AJ 163 69; Shannon et al. 2025 PASA 42 e036; Connor et al. 2023 ApJL 949 L26; Pastor-Marazuela et al. 2026 MNRAS 545; Caleb et al. 2023 MNRAS 524 2064).
- Record every value you see for a host, even when two papers disagree; put both rows in with their own `source_short`.
- Fetching: arXiv via `curl.exe -A "Mozilla/5.0" -L https://arxiv.org/abs/<id>` or the HTML rendering `https://arxiv.org/html/<id>`; plain `urllib` returns HTTP 406. PDFs: `curl.exe -A "Mozilla/5.0" -L -o <file> https://arxiv.org/pdf/<id>` then extract text with Python (`pypdf` or `pdfplumber`); install into the current environment if missing. Keep downloaded PDFs under `Publication/materials/host_properties/papers/`.
- Windows PowerShell: `&&` does not work, use `;`. Avoid printing non-ASCII in `python -c`; write to files with `encoding="utf-8"`.
- Do not edit any file outside `Publication/materials/host_properties/`. Do not run git commands that change state.
