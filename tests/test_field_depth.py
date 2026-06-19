"""Unit tests for pipeline_scripts/photometry + astropath/field_depth.py."""

from __future__ import annotations

import importlib
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir, "pipeline_scripts"))
sys.path.insert(
    0,
    os.path.join(os.path.dirname(__file__), os.pardir, "pipeline_scripts", "photometry + astropath"),
)

field_depth = importlib.import_module("field_depth")
m_lim_5sigma_aperture = field_depth.m_lim_5sigma_aperture
production_aperture_diameter_px = field_depth.production_aperture_diameter_px


class TestMLim5SigmaAperture:
    def test_returns_tuple_of_three_floats(self):
        result = m_lim_5sigma_aperture(25.0, 0.01, 40.0)
        assert isinstance(result, tuple)
        assert len(result) == 3
        assert all(isinstance(v, float) for v in result)


class TestProductionApertureDiameterPx:
    def test_none_none_returns_default_max(self):
        assert production_aperture_diameter_px(None, None) == 40.0

    def test_explicit_config_diam(self):
        assert production_aperture_diameter_px(30.0, None) == 30.0
