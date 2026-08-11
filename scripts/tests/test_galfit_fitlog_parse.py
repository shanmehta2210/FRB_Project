"""Unit tests for scripts/galfit_fitlog_parse.py."""

from __future__ import annotations

import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, "scripts"))

from galfit_fitlog_parse import (
    inclination_err_from_b_a_err,
    inclination_from_b_a,
    parse_fitlog_block,
)


class TestInclinationFromBA:
    def test_face_on(self):
        assert inclination_from_b_a(1.0) == 0.0

    def test_edge_on_at_q0(self):
        assert inclination_from_b_a(0.2) == 90.0

    def test_below_q0_clamped(self):
        assert inclination_from_b_a(0.1) == 90.0

    def test_none_returns_none(self):
        assert inclination_from_b_a(None) is None


class TestInclinationErrFromBAErr:
    def test_finite_positive(self):
        result = inclination_err_from_b_a_err(0.5, 0.05)
        assert isinstance(result, float)
        assert result > 0.0

    def test_none_err_returns_zero(self):
        assert inclination_err_from_b_a_err(0.5, None) == 0.0


class TestParseFitlogBlock:
    def test_empty_string(self):
        data, strategy = parse_fitlog_block("")
        assert data == {}
        assert strategy == "empty"
