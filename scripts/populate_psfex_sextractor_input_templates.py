import os
import shutil
import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BASE = os.path.join(ROOT, "psfs", "PSFEx + SExtractor")
TEMPLATES = os.path.join(BASE, "templates")
MANIFEST = os.path.join(BASE, "download_manifest_10arcmin.csv")

SEX_TEMPLATE = os.path.join(TEMPLATES, "sextractor_psfex.sex")
PARAM_TEMPLATE = os.path.join(TEMPLATES, "sextractor_psfex.param")
CONV_TEMPLATE = os.path.join(TEMPLATES, "default.conv")
NNW_TEMPLATE = os.path.join(TEMPLATES, "default.nnw")
PSFEX_TEMPLATE = os.path.join(TEMPLATES, "psfex_default.psfex")


def write_text(path: str, text: str) -> None:
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)


def main() -> None:
    if not os.path.exists(MANIFEST):
        raise FileNotFoundError(f"Missing manifest: {MANIFEST}")

    for p in [SEX_TEMPLATE, PARAM_TEMPLATE, CONV_TEMPLATE, NNW_TEMPLATE, PSFEX_TEMPLATE]:
        if not os.path.exists(p):
            raise FileNotFoundError(f"Missing template: {p}")

    df = pd.read_csv(MANIFEST)
    updated = 0

    for _, row in df.iterrows():
        run_dir = row["run_dir"]
        frb = row["FRB"]

        sex_in = os.path.join(run_dir, "input", "sextractor")
        psfex_in = os.path.join(run_dir, "input", "psfex")
        cutouts = os.path.join(run_dir, "cutouts")

        os.makedirs(sex_in, exist_ok=True)
        os.makedirs(psfex_in, exist_ok=True)

        # Copy static templates.
        shutil.copy2(SEX_TEMPLATE, os.path.join(sex_in, "default.sex"))
        shutil.copy2(PARAM_TEMPLATE, os.path.join(sex_in, "default.param"))
        shutil.copy2(CONV_TEMPLATE, os.path.join(sex_in, "default.conv"))
        shutil.copy2(NNW_TEMPLATE, os.path.join(sex_in, "default.nnw"))
        shutil.copy2(PSFEX_TEMPLATE, os.path.join(psfex_in, "default.psfex"))

        # Create expected local symlink targets via copy for simple command lines.
        flux_src = os.path.join(cutouts, f"{frb}_10arcmin_flux.fits")
        invvar_src = os.path.join(cutouts, f"{frb}_10arcmin_invvar.fits")

        flux_dst = os.path.join(sex_in, "image.fits")
        invvar_dst = os.path.join(sex_in, "invvar.fits")

        if os.path.exists(flux_src):
            shutil.copy2(flux_src, flux_dst)
        if os.path.exists(invvar_src):
            shutil.copy2(invvar_src, invvar_dst)

        # WSL-friendly run snippets for reference.
        run_sex = """#!/usr/bin/env bash
set -euo pipefail

# Run from this directory: input/sextractor
source-extractor image.fits -c default.sex
"""
        write_text(os.path.join(sex_in, "run_source_extractor.sh"), run_sex)

        run_psfex = """#!/usr/bin/env bash
set -euo pipefail

# Run from this directory: input/psfex
# Assumes sextractor_catalog.fits has been copied here from input/sextractor output.
psfex sextractor_catalog.fits -c default.psfex
"""
        write_text(os.path.join(psfex_in, "run_psfex.sh"), run_psfex)

        updated += 1

    print(f"Populated input templates for {updated} FRB runs.")


if __name__ == "__main__":
    main()
