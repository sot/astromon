import numpy as np
import pytest
from astropy.table import Table

from astromon.source_detection import (
    _fit,
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


class TestFitBoundsMatchBoxSize:
    """`_fit`'s centroid bounds must match the box the caller actually drew events
    from (`box_size`), not some other fixed value -- see astromon/observation.py's
    `_seed_and_select_events`, which selects events within `box_size` of the seed
    before handing them to the fit. A looser bound lets the optimizer converge on a
    centroid outside the data it was given: at obsid 15669 COMPONENT 1 (a faint,
    SNR~3.2 celldetect source), the old hardcoded +/-10" bound let a
    background-dominated fit (snr=0.09) wander 10.5" from its seed and outside its
    own 4" fit box, while still reporting fit_ok=True.
    """

    def test_bounds_scale_with_box_size(self, monkeypatch):
        captured = {}
        real_minimize = __import__("scipy").optimize.minimize

        def spy_minimize(fun, x0, bounds, **kwargs):
            captured["bounds"] = bounds
            return real_minimize(fun, x0, bounds=bounds, **kwargs)

        monkeypatch.setattr("scipy.optimize.minimize", spy_minimize)

        yag, zag = _point_source_events(n=200)
        events = Table({"y_angle": yag, "z_angle": zag})
        source = {"y_angle": 0.0, "z_angle": 0.0, "COMPONENT": 1}

        for box_size in (2, 4, 9):
            captured.clear()
            _fit(events, source, box_size=box_size)
            assert captured["bounds"][0] == (-box_size, box_size)
            assert captured["bounds"][1] == (-box_size, box_size)

    @pytest.mark.parametrize("box_size", [2, 4, 8])
    def test_fitted_centroid_never_exceeds_box_size(self, box_size):
        """Even when the only "signal" in the box sits right at its edge, the fitted
        centroid must not be pushed past box_size from the seed -- it has nowhere
        else to go for data it wasn't given.
        """
        rng = np.random.default_rng(3)
        n_edge = 150
        n_bkg = 150
        yag = np.concatenate(
            [
                rng.normal(box_size - 0.2, 0.2, n_edge),
                rng.uniform(-box_size, box_size, n_bkg),
            ]
        )
        zag = np.concatenate(
            [
                rng.normal(0, 0.2, n_edge),
                rng.uniform(-box_size, box_size, n_bkg),
            ]
        )
        events = Table({"y_angle": yag, "z_angle": zag})
        source = {"y_angle": 0.0, "z_angle": 0.0, "COMPONENT": 1}

        result = fit_gaussian_2d(events, source, box_size=box_size)

        if result["fit_ok"]:
            assert abs(result["y_angle"]) <= box_size
            assert abs(result["z_angle"]) <= box_size

    def test_rejects_a_fit_pinned_at_the_box_boundary(self, monkeypatch):
        """A centroid that lands (to numerical precision) exactly on a fit bound
        means the optimizer wanted to move further but couldn't -- there is no
        localized signal in the box, so `fit_gaussian_2d` must not report this as a
        usable fit (fit_ok=True) even though scipy itself calls it "converged"
        (result.success=True). This is precisely how obsid 15669 COMPONENT 1 (a
        faint, SNR~3.2 celldetect source) got a fit_ok=True fit 10.5" from its own
        seed and outside its own 4" fit box before this bounds/guard fix.

        `_fit`'s own optimizer dynamics on a low-SNR dataset are not reliably
        reproducible in a fast, deterministic unit test, so this drives the guard
        directly: it forges the `_fit` return value that a pinned-at-bound
        optimizer would produce and checks `fit_gaussian_2d` catches it.
        """
        box_size = 4
        source = {"y_angle": 0.0, "z_angle": 0.0, "COMPONENT": 1}
        # Events sit right at the edge of the box (not at the origin), so the
        # forged fitted centroid below has real, nearby events to "explain" -- this
        # keeps the pre-existing `n == 0` fallback from masking whether the new
        # boundary guard is the thing actually rejecting the fit.
        yag, zag = _point_source_events(yag0=box_size - 0.1, zag0=0.0, n=50)
        events = Table({"y_angle": yag, "z_angle": zag})

        class _FakeResult:
            success = True
            # centroid pinned exactly at the y_angle bound; everything else benign.
            x = np.array([box_size, 0.0, 1.0, 1.0, 0.0, 3.0])
            hess_inv = type("_H", (), {"todense": lambda self: np.eye(6)})()

        monkeypatch.setattr(
            "astromon.source_detection._fit", lambda *args, **kwargs: _FakeResult()
        )

        result = fit_gaussian_2d(events, source, box_size=box_size)

        assert result["fit_ok"] is False
