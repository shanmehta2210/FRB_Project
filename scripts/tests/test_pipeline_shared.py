"""Unit tests for pipeline_scripts/pipeline_shared.py."""

from __future__ import annotations

import logging
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), os.pardir, os.pardir, "pipeline_scripts"))

from pipeline_shared import (
    DEFAULT_APERTURE_DIAMS_PX,
    format_phot_apertures,
    get_logger,
    header_mag_zeropoint_from_fits,
    render_param_template,
    resolve_apertures,
)


class TestResolveApertures:
    def test_empty_config_returns_defaults(self):
        apertures, prod_idx, prod_diam = resolve_apertures({})
        assert len(apertures) == 15
        assert prod_idx == 14
        assert prod_diam == 40.0

    def test_explicit_list_and_production(self):
        cfg = {"phot_apertures_px": [5, 10, 20], "production_aperture_px": 10}
        apertures, prod_idx, prod_diam = resolve_apertures(cfg)
        assert apertures == [5.0, 10.0, 20.0]
        assert prod_idx == 1
        assert prod_diam == 10.0

    def test_legacy_scalar_appends_to_base(self):
        cfg = {"phot_apertures_px": 40.0}
        apertures, prod_idx, prod_diam = resolve_apertures(cfg)
        assert len(apertures) == 15
        assert prod_idx == 14
        assert prod_diam == 40.0

    def test_list_without_production_defaults_to_largest(self):
        cfg = {"phot_apertures_px": [5, 10, 20]}
        apertures, prod_idx, prod_diam = resolve_apertures(cfg)
        assert prod_idx == 2
        assert prod_diam == 20.0


class TestFormatPhotApertures:
    def test_basic_formatting(self):
        result = format_phot_apertures([4.0, 5.0, 40.0])
        assert result == "4, 5, 40"


class TestRenderParamTemplate:
    def test_substitution(self):
        template = "FLUX_APER({NAPER})\nMAG_APER({NAPER})"
        result = render_param_template(template, 15)
        assert result == "FLUX_APER(15)\nMAG_APER(15)"


class TestGetLogger:
    def test_returns_child_of_frb_pipeline(self):
        logger = get_logger("test")
        assert isinstance(logger, logging.Logger)
        assert logger.name.startswith("frb_pipeline.")


class TestHeaderMagZeropointFromFits:
    def test_missing_file_returns_none(self):
        assert header_mag_zeropoint_from_fits("/nonexistent/path.fits") is None

    def test_ps1_stack_header(self):
        repo = os.path.join(os.path.dirname(__file__), os.pardir)
        ps1 = os.path.join(
            repo,
            "CHIME",
            "large_cutouts",
            "20180814A_flux.fits",
        )
        if not os.path.isfile(ps1):
            return  # optional asset — skip in minimal checkouts
        zp = header_mag_zeropoint_from_fits(ps1)
        assert zp is not None
        assert 24.0 <= zp <= 26.0
