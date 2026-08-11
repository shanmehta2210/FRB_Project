"""Regression tests for Phase 1 / Phase 2 catalog NUMBER cross-matching."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
from astropy.table import Table

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "pipeline_scripts" / "galfit_fitting"))

import generate_galfit_cutouts as ggc  # noqa: E402


def _psf_table():
    """Phase 2 catalog: NUMBER 1176 is a star at different sky position than cat #1176."""
    return Table(
        names=[
            "NUMBER",
            "ALPHAWIN_J2000",
            "DELTAWIN_J2000",
            "SPREAD_MODEL",
            "SPREADERR_MODEL",
        ],
        data=[
            [1176, 1206],
            [76.953, 77.015],
            [26.066, 26.061],
            [0.0005, 0.0158],
            [0.0001, 0.0003],
        ],
    )


def _phase1_cat():
    return Table(
        names=["NUMBER", "ALPHAWIN_J2000", "DELTAWIN_J2000"],
        data=[[1176], [77.015], [26.061]],
    )


@pytest.fixture()
def psf_cat():
    return _psf_table()


@pytest.fixture()
def spread_by_number(psf_cat):
    return {
        int(row["NUMBER"]): (float(row["SPREAD_MODEL"]), float(row["SPREADERR_MODEL"]))
        for row in psf_cat
    }


def test_spread_for_catalog_index_uses_sky_match(psf_cat, spread_by_number, monkeypatch):
    cat = _phase1_cat()
    monkeypatch.setattr(ggc, "_load_psf_catalog", lambda _path: psf_cat)

    spread = ggc._spread_for_catalog_index(0, cat, spread_by_number, "image.cat")
    assert spread == pytest.approx((0.0158, 0.0003))
    assert not ggc._spread_is_point_source(*spread)


def test_spread_for_seg_number_uses_sky_match(psf_cat, spread_by_number, monkeypatch):
    cat = _phase1_cat()
    monkeypatch.setattr(ggc, "_load_psf_catalog", lambda _path: psf_cat)

    spread = ggc._spread_for_seg_number(1176, cat, spread_by_number, "image.cat")
    assert spread == pytest.approx((0.0158, 0.0003))


def test_direct_number_key_would_misclassify_host(spread_by_number):
    """Document the bug: naive NUMBER key returns the wrong (stellar) SPREAD."""
    wrong = spread_by_number[1176]
    assert ggc._spread_is_point_source(*wrong)
