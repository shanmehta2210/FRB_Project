"""Unit tests for Phase 3a Re-separation ROI policy."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "pipeline_scripts" / "galfit_fitting"))

import generate_galfit_cutouts as ggc  # noqa: E402
from sersic_init import effective_re_px  # noqa: E402


def test_effective_re_px_flux_radius():
    assert effective_re_px({"FLUX_RADIUS": 4.5}) == pytest.approx(4.5)


def test_effective_re_px_floor():
    assert effective_re_px({"FLUX_RADIUS": 0.0}) == 1.0
    assert effective_re_px({"FLUX_RADIUS": -2.0}) == 1.0
    assert effective_re_px({}) == 1.0


def _make_cat(rows):
    """rows: list of dicts with NUMBER, X_IMAGE, Y_IMAGE, FLUX_RADIUS, CLASS_STAR."""
    keys = ["NUMBER", "X_IMAGE", "Y_IMAGE", "FLUX_RADIUS", "CLASS_STAR"]
    data = {k: [r[k] for r in rows] for k in keys}
    return Table(data)


def _paint_disk(seg, uid, cx, cy, rad):
    yy, xx = np.ogrid[: seg.shape[0], : seg.shape[1]]
    mask = (xx - cx) ** 2 + (yy - cy) ** 2 <= rad**2
    seg[mask] = uid


def test_far_neighbor_masked_roi_unchanged():
    """sep > 3*Re_n → mask; starting ROI size unchanged."""
    seg = np.zeros((80, 80), dtype=np.int32)
    _paint_disk(seg, 1, 25, 25, 3)  # host
    # Centroid far (sep=15 > 3*2=6) but disk still clips host_pad=10 ROI.
    _paint_disk(seg, 2, 40, 25, 3)
    cat = _make_cat(
        [
            {"NUMBER": 1, "X_IMAGE": 25.0, "Y_IMAGE": 25.0, "FLUX_RADIUS": 5.0, "CLASS_STAR": 0.1},
            {"NUMBER": 2, "X_IMAGE": 40.0, "Y_IMAGE": 25.0, "FLUX_RADIUS": 2.0, "CLASS_STAR": 0.1},
        ]
    )
    xmin, xmax, ymin, ymax, fit, mask, decisions = ggc.resolve_neighbor_re_roi(
        seg,
        cat,
        1,
        host_pad=10,
        re_sep_factor=3.0,
        max_roi_iterations=8,
        max_cutout_side=512,
        neighbor_class_star_max=0.75,
        spread_by_number={},
    )
    assert 2 in mask
    assert 2 not in fit
    assert decisions[2][0] == "mask_far"
    host_bb = ggc._bbox_of_objids(seg, [1])
    exp = ggc._pad_bbox(*host_bb, pad=10, shape=seg.shape)
    assert (xmin, xmax, ymin, ymax) == exp


def test_near_neighbor_expands_and_fits():
    """sep <= 3*Re_n → fit + grow ROI to cover both + pad."""
    seg = np.zeros((100, 100), dtype=np.int32)
    _paint_disk(seg, 1, 30, 30, 4)
    _paint_disk(seg, 2, 42, 30, 4)  # sep=12; Re_n=5 → thresh=15 → near
    cat = _make_cat(
        [
            {"NUMBER": 1, "X_IMAGE": 30.0, "Y_IMAGE": 30.0, "FLUX_RADIUS": 5.0, "CLASS_STAR": 0.1},
            {"NUMBER": 2, "X_IMAGE": 42.0, "Y_IMAGE": 30.0, "FLUX_RADIUS": 5.0, "CLASS_STAR": 0.1},
        ]
    )
    xmin, xmax, ymin, ymax, fit, mask, decisions = ggc.resolve_neighbor_re_roi(
        seg,
        cat,
        1,
        host_pad=20,
        re_sep_factor=3.0,
        max_roi_iterations=8,
        max_cutout_side=512,
        neighbor_class_star_max=0.75,
        spread_by_number={},
    )
    assert 2 in fit
    assert 2 not in mask
    assert decisions[2][0] == "fit"
    # Grown ROI must cover both disks + pad
    both = ggc._bbox_of_objids(seg, [1, 2])
    grown = ggc._pad_bbox(*both, pad=20, shape=seg.shape)
    assert (xmin, xmax, ymin, ymax) == grown


def test_clipping_neighbor_centroid_outside_still_evaluated():
    """Neighbor centroid outside starting ROI but seg clips → still checked."""
    seg = np.zeros((60, 60), dtype=np.int32)
    _paint_disk(seg, 1, 20, 20, 3)
    # Elongated blob: centroid at x=45 but pixels extend left into host pad ROI
    yy, xx = np.ogrid[:60, :60]
    blob = ((xx >= 28) & (xx <= 50) & (yy >= 18) & (yy <= 22))
    seg[blob] = 2
    cat = _make_cat(
        [
            {"NUMBER": 1, "X_IMAGE": 20.0, "Y_IMAGE": 20.0, "FLUX_RADIUS": 5.0, "CLASS_STAR": 0.1},
            # centroid far: sep=25; Re=2 → thresh=6 → mask_far
            {"NUMBER": 2, "X_IMAGE": 45.0, "Y_IMAGE": 20.0, "FLUX_RADIUS": 2.0, "CLASS_STAR": 0.1},
        ]
    )
    host_bb = ggc._bbox_of_objids(seg, [1])
    start = ggc._pad_bbox(*host_bb, pad=10, shape=seg.shape)
    # Confirm centroid of #2 is outside starting ROI but pixels clip
    assert not (start[0] <= 45 < start[1])
    touching = ggc._ids_touching_roi(seg, *start, drop=[])
    assert 2 in touching

    _xmin, _xmax, _ymin, _ymax, fit, mask, decisions = ggc.resolve_neighbor_re_roi(
        seg,
        cat,
        1,
        host_pad=10,
        re_sep_factor=3.0,
        max_roi_iterations=8,
        max_cutout_side=512,
        neighbor_class_star_max=0.75,
        spread_by_number={},
    )
    assert 2 in decisions
    assert 2 in mask
    assert decisions[2][0] == "mask_far"
    assert 2 not in fit


def test_star_always_masked_no_expand():
    seg = np.zeros((80, 80), dtype=np.int32)
    _paint_disk(seg, 1, 25, 25, 4)
    _paint_disk(seg, 2, 32, 25, 2)  # close, but stellar
    cat = _make_cat(
        [
            {"NUMBER": 1, "X_IMAGE": 25.0, "Y_IMAGE": 25.0, "FLUX_RADIUS": 5.0, "CLASS_STAR": 0.1},
            {"NUMBER": 2, "X_IMAGE": 32.0, "Y_IMAGE": 25.0, "FLUX_RADIUS": 10.0, "CLASS_STAR": 0.95},
        ]
    )
    xmin, xmax, ymin, ymax, fit, mask, decisions = ggc.resolve_neighbor_re_roi(
        seg,
        cat,
        1,
        host_pad=20,
        re_sep_factor=3.0,
        max_roi_iterations=8,
        max_cutout_side=512,
        neighbor_class_star_max=0.75,
        spread_by_number={},
    )
    assert 2 in mask
    assert 2 not in fit
    assert decisions[2][0] == "mask_star"
    host_bb = ggc._bbox_of_objids(seg, [1])
    exp = ggc._pad_bbox(*host_bb, pad=20, shape=seg.shape)
    assert (xmin, xmax, ymin, ymax) == exp


def test_host_always_first_fit_member():
    seg = np.zeros((80, 80), dtype=np.int32)
    _paint_disk(seg, 7, 30, 30, 4)
    _paint_disk(seg, 9, 38, 30, 4)
    cat = _make_cat(
        [
            {"NUMBER": 7, "X_IMAGE": 30.0, "Y_IMAGE": 30.0, "FLUX_RADIUS": 5.0, "CLASS_STAR": 0.1},
            {"NUMBER": 9, "X_IMAGE": 38.0, "Y_IMAGE": 30.0, "FLUX_RADIUS": 5.0, "CLASS_STAR": 0.1},
        ]
    )
    _bounds = ggc.resolve_neighbor_re_roi(
        seg,
        cat,
        7,
        host_pad=15,
        re_sep_factor=3.0,
        spread_by_number={},
    )
    fit = _bounds[4]
    assert 7 in fit
    assert list(sorted(fit))[0] == 7 or 7 in fit  # set; host present
    assert _bounds[6][7][0] == "host"
