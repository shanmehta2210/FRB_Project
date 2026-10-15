# ApJL publication workspace

Local drafting and decision records for the expanded FRB host-inclination study.
This directory is excluded by the repository `.gitignore` and must not be
uploaded to GitHub.

## Working manuscript structure

1. Abstract — deferred
2. Introduction — drafted
3. Methodology — drafted
4. Results and Analysis — pending
5. Discussion — pending
6. Conclusions — pending
7. Appendices — fit and data-health checks

ApJL guidance: keep the abstract within 250 words, main text near 3500 words,
and figures plus tables near five combined.

## Files

| File | Purpose |
|---|---|
| `manuscript/manuscript.tex` | Standalone Introduction and Methodology source |
| `manuscript/manuscript.pdf` | Compiled review copy |
| `paper_decisions.md` | Locked choices, unresolved issues, and evidence map |
| `Bhardwaj et al. - 2024 - Selection bias obfuscates the discovery of fast radio burst sources.pdf` | Persistent scientific and stylistic reference for the new manuscript |

Compile from `Publication/manuscript/` with:

```powershell
& "$env:LOCALAPPDATA\Programs\MiKTeX\miktex\bin\x64\xelatex.exe" manuscript.tex
```

Run twice after changes that affect references. The current review PDF uses
the self-contained two-column `article` class because the MiKTeX package
repository timed out while retrieving `aastex`; switch the class to
`aastex631` before journal submission.

## Current sample chain

The repository currently encodes:

$$
64\ {\rm production\ fits}
\rightarrow 60\ (m_r\leq22)
\rightarrow 53\ (m_r\leq22,\ q>0.2)
\rightarrow 50\ {\rm confirmed}.
$$

The manuscript should use the final 50-host cohort unless the sample policy is
changed. The intermediate magnitude count is 60 with the current GALFIT `mag`
column and 59 with `mag_final`; the combined magnitude-plus-axis-ratio cohort is
53 either way. Resolve the reported magnitude definition before submission.

## Drafting principles

- Cite Bhardwaj et al. (2024) in the opening paragraph.
- Use the Bhardwaj et al. (2024) PDF as the continuing reference for argument
  flow and scientific tone, without carrying forward unsupported sample claims.
- Present the new work as a larger, independently reduced test of the
  inclination-selection result.
- Do not repeat the old paper's CRAFT/DSA-only or non-repeater selection unless
  the working sample is changed to satisfy it.
- Do not claim a uniform PATH $P(O)>0.95$ cut. The present cohort uses published
  host localizations; pipeline AstroPath is an independent cross-check.
- Use the sky-aware uncertainty protocol. Do not use the obsolete
  $10\times\sigma_q$ prescription.
- Describe HSC correctly: its structural catalog uses Lenstronomy single-Sérsic
  fits, not GALFIT. The parity is in the model class, PSF convolution, and fitted
  axis ratio.
