import numpy as np
import pytest
from astromon.source_detection import find_local_peak, concentration_ratio


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
        peak_yag, peak_zag = find_local_peak(yag, zag, seed_yag=0.0, seed_zag=0.0, box_size=box)
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
