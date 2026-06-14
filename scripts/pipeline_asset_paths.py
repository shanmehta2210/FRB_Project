"""Paths to pipeline host cutouts (replaces removed cropped_host_galaxies/)."""

from __future__ import annotations

import glob
import os
from typing import Iterator


def repo_root(start: str | None = None) -> str:
    if start:
        return start
    return os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


def pipeline_output(root: str | None = None) -> str:
    return os.path.join(repo_root(root), "pipeline_scripts", "Output")


def host_run_dir(frb: str, root: str | None = None) -> str:
    return os.path.join(pipeline_output(root), f"{frb}_all")


def host_cutout_path(frb: str, root: str | None = None) -> str:
    return os.path.join(host_run_dir(frb, root), "host_cutout.fits")


def host_sigma_path(frb: str, root: str | None = None) -> str:
    return os.path.join(host_run_dir(frb, root), "host_sigma.fits")


def frb_from_guess_filename(filename: str) -> str:
    return filename.replace("_flux.fits", "").replace("host_cutout.fits", "")


def resolve_host_files(frb: str, root: str | None = None) -> tuple[str | None, str | None]:
    cutout = host_cutout_path(frb, root)
    sigma = host_sigma_path(frb, root)
    cutout_ok = os.path.isfile(cutout)
    sigma_ok = os.path.isfile(sigma)
    return (cutout if cutout_ok else None, sigma if sigma_ok else None)


def iter_host_cutouts(root: str | None = None) -> Iterator[tuple[str, str]]:
    out_root = pipeline_output(root)
    if not os.path.isdir(out_root):
        return
    for run_dir in sorted(glob.glob(os.path.join(out_root, "*_all"))):
        cutout = os.path.join(run_dir, "host_cutout.fits")
        if os.path.isfile(cutout):
            frb = os.path.basename(run_dir).replace("_all", "")
            yield frb, cutout
