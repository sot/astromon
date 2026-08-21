import numpy as np
import pytest
from astropy.table import Table

from astromon.source_detection import (
    concentration_ratio,
    find_local_peak,
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


class TestFindLocalPeak:
    def test_recovers_point_source_peak(self):
        yag, zag = _point_source_events(yag0=1.0, zag0=-0.5)
        peak_yag, peak_zag = find_local_peak(yag, zag, seed_yag=0.0, seed_zag=0.0)
        assert abs(peak_yag - 1.0) < 0.5
        assert abs(peak_zag - (-0.5)) < 0.5

    def test_falls_back_on_too_few_events(self):
        yag = np.array([0.1, 0.2])
        zag = np.array([0.1, 0.2])
        peak_yag, peak_zag = find_local_peak(yag, zag, seed_yag=5.0, seed_zag=5.0)
        assert peak_yag == 5.0
        assert peak_zag == 5.0

    def test_stays_within_box(self):
        rng = np.random.default_rng(0)
        yag = rng.uniform(-20, 20, 2000)
        zag = rng.uniform(-20, 20, 2000)
        box = 4.0
        peak_yag, peak_zag = find_local_peak(
            yag, zag, seed_yag=0.0, seed_zag=0.0, box_size=box
        )
        assert abs(peak_yag) <= box
        assert abs(peak_zag) <= box


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
