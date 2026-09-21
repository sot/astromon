from unittest.mock import patch

import numpy as np
import pytest
import scipy.optimize
from astropy.table import Table

from astromon.source_detection import (
    _fit,
    concentration_ratio,
    fit_gaussian_2d,
)


def _point_source_events(yag0=0.0, zag0=0.0, sigma=0.5, n=500, seed=42):
    rng = np.random.default_rng(seed)
    yag = rng.normal(yag0, sigma, n)
    zag = rng.normal(zag0, sigma, n)
    return yag, zag


def _extended_events(yag0=0.0, zag0=0.0, sigma=8.0, n=5000, seed=42):
    rng = np.random.default_rng(seed)
    yag = rng.normal(yag0, sigma, n)
    zag = rng.normal(zag0, sigma, n)
    return yag, zag


def _non_converging_result(x0):
    """A scipy.optimize.OptimizeResult as returned by a failed minimize() call."""
    return scipy.optimize.OptimizeResult(
        x=np.asarray(x0, dtype=float),
        success=False,
        message="Non-convergent (forced for test)",
        hess_inv=1e10 * np.eye(len(x0)),
        fun=np.inf,
        nit=0,
    )


class TestConcentrationRatio:
    def test_point_source_near_one(self):
        yag, zag = _point_source_events(n=2000)
        cr = concentration_ratio(yag, zag, 0.0, 0.0, r_core_as=2.0, r_extract_as=10.0)
        assert cr > 0.85

    def test_extended_source_near_zero(self):
        yag, zag = _extended_events(n=5000)
        cr = concentration_ratio(yag, zag, 0.0, 0.0, r_core_as=2.0, r_extract_as=10.0)
        assert cr < 0.1

    def test_returns_nan_when_no_events_in_aperture(self):
        yag = np.array([100.0, 101.0])
        zag = np.array([100.0, 101.0])
        cr = concentration_ratio(yag, zag, 0.0, 0.0)
        assert np.isnan(cr)

    def test_ratio_in_unit_interval(self):
        rng = np.random.default_rng(1)
        yag = rng.normal(0, 2.0, 1000)
        zag = rng.normal(0, 2.0, 1000)
        cr = concentration_ratio(yag, zag, 0.0, 0.0)
        assert 0.0 <= cr <= 1.0


class TestFitGaussian2dPreservesComponent:
    """gaussian_detect and peak_gaussian_detect both call fit_gaussian_2d exactly once
    per celldetect source, and both rely on it returning that source's COMPONENT value
    unchanged (see astromon/observation.py: both pass "COMPONENT": source["COMPONENT"]
    straight through).  That invariant is what keeps a catalog match's x_id -- assigned
    from celldetect's source IDs in get_cat_obs_data.py -- valid when joined against
    gaussian_detect/peak_gaussian_detect's xray sources without any remapping.  If this
    ever stops holding, astromon_21/icrf3/rfc/tycho2 cross-matches for those two detect
    methods would join against the wrong (or no) x-ray source.
    """

    def test_component_preserved_on_successful_fit(self):
        yag, zag = _point_source_events(yag0=0.0, zag0=0.0, sigma=0.5, n=500)
        events = Table({"y_angle": yag, "z_angle": zag})
        source = {"y_angle": 0.0, "z_angle": 0.0, "COMPONENT": 7}

        result = fit_gaussian_2d(events, source)

        assert result["fit_ok"] is True
        assert result["COMPONENT"] == 7

    def test_component_preserved_when_too_few_events(self):
        """The <10-events early return (fail_value) must also preserve COMPONENT."""
        events = Table({"y_angle": [0.1, 0.2], "z_angle": [0.1, 0.2]})
        source = {"y_angle": 0.0, "z_angle": 0.0, "COMPONENT": 42}

        result = fit_gaussian_2d(events, source)

        assert result["fit_ok"] is False
        assert result["COMPONENT"] == 42

    @pytest.mark.parametrize("component", [0, 1, 5, 999])
    def test_component_value_round_trips_exactly(self, component):
        yag, zag = _point_source_events(n=500)
        events = Table({"y_angle": yag, "z_angle": zag})
        source = {"y_angle": 0.0, "z_angle": 0.0, "COMPONENT": component}

        result = fit_gaussian_2d(events, source)

        assert result["COMPONENT"] == component


class TestFitDoesNotRaiseOnNonConvergence:
    """_fit used to raise RuntimeError when scipy.optimize.minimize did not converge,
    which meant fit_gaussian_2d's own `if not result.success: return fail_value` check
    was dead code -- the exception had already propagated past it. That in turn crashed
    _fit_gaussian_sources's per-source loop (astromon/observation.py), aborting fit
    detection for every other source in the same obsid instead of recording just the
    one failed fit. _fit must return the failed OptimizeResult instead, so the caller's
    existing fail_value branch can do its job.
    """

    def test_fit_returns_failed_result_instead_of_raising(self):
        yag, zag = _point_source_events(n=500)
        events = Table({"y_angle": yag, "z_angle": zag})
        source = {"y_angle": 0.0, "z_angle": 0.0}

        with patch.object(
            scipy.optimize,
            "minimize",
            side_effect=lambda _fn, x0, **kwargs: _non_converging_result(x0),
        ):
            result = _fit(events, source)

        assert result.success is False

    def test_fit_gaussian_2d_returns_fail_value_on_non_convergence(self):
        yag, zag = _point_source_events(n=500)
        events = Table({"y_angle": yag, "z_angle": zag})
        source = {"y_angle": 0.0, "z_angle": 0.0, "COMPONENT": 3}

        with patch.object(
            scipy.optimize,
            "minimize",
            side_effect=lambda _fn, x0, **kwargs: _non_converging_result(x0),
        ):
            result = fit_gaussian_2d(events, source)

        assert result["fit_ok"] is False
        assert result["COMPONENT"] == 3
