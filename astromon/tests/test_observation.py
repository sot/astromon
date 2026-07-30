import re
from unittest.mock import Mock, patch

import numpy as np
import pytest
from astropy import table

from astromon import observation

# The real HRC plate scale, per the Chandra Proposers' Observatory Guide. filter_events and
# filter_sources both filter events/sources in a circle around the optical axis, and need to
# convert the filter radius from arcsec to detector pixels. Using the wrong scale here does not
# raise an exception: it silently changes the filtered region size (see astromon/observation.py).
HRC_ARCSEC_PER_PIXEL = 0.13180
ACIS_ARCSEC_PER_PIXEL = 0.5

CIRCLE_RE = re.compile(r"circle\([^,]+,[^,]+,([^)]+)\)")


def _sources_table(**extra_cols):
    """Minimal sources table as returned by _get_sources() after renaming.

    Includes all columns _get_sources() always produces. Pass extra_cols to add
    columns only present for certain detect versions (e.g. psfratio, concentration_ratio).
    """
    n = 2
    cols = {
        "id": np.array([0, 1], dtype=np.int32),
        "ra": np.array([83.8, 83.9]),
        "dec": np.array([-5.4, -5.5]),
        "net_counts": np.array([120.0, 55.0], dtype=np.float32),
        "snr": np.array([9.0, 5.5], dtype=np.float32),
        "y_angle": np.array([4.0, -2.0]),
        "z_angle": np.array([1.5, 0.8]),
        "r_angle": np.array([4.3, 2.2]),
        "near_neighbor_dist": np.array([80.0, 80.0]),
        "pileup": np.array([0.0, 0.0], dtype=np.float32),
        "acis_streak": np.array([False, False]),
        "caldb_version": np.array(["4.9.4", "4.9.4"]),
        "obsid": np.array([1234, 1234], dtype=np.int32),
    }
    cols.update(extra_cols)
    return table.Table(cols)


_GET_SOURCES_FUNC = observation.Observation.__dict__["get_sources"].func


def test_get_sources_fills_gaussian_only_columns_for_celldetect():
    """Columns only produced by gaussian-family detectors (psfratio, concentration_ratio)
    are filled with NaN when get_sources() is called with version='celldetect'."""
    celldetect_sources = _sources_table()  # no psfratio or concentration_ratio

    obs = Mock()
    obs._get_sources.return_value = celldetect_sources

    with patch.object(observation.TASKS.tasks["celldetect"], "get_result", return_value=None):
        result = _GET_SOURCES_FUNC(obs, version="celldetect")

    assert "concentration_ratio" in result.colnames
    assert "psfratio" in result.colnames
    assert np.all(np.isnan(result["concentration_ratio"]))
    assert np.all(np.isnan(result["psfratio"]))
    assert result["detect_method"].tolist() == ["celldetect", "celldetect"]


def test_get_sources_passes_through_gaussian_columns():
    """concentration_ratio and psfratio are preserved as-is for gaussian_detect sources."""
    gaussian_sources = _sources_table(
        psfratio=np.array([0.85, 0.72], dtype=np.float32),
        concentration_ratio=np.array([0.6, 0.5], dtype=np.float32),
    )

    obs = Mock()
    obs._get_sources.return_value = gaussian_sources

    with patch.object(observation.TASKS.tasks["gaussian_detect"], "get_result", return_value=None):
        result = _GET_SOURCES_FUNC(obs, version="gaussian_detect")

    np.testing.assert_array_almost_equal(result["psfratio"], [0.85, 0.72])
    np.testing.assert_array_almost_equal(result["concentration_ratio"], [0.6, 0.5])
    assert result["detect_method"].tolist() == ["gaussian_detect", "gaussian_detect"]


def _fake_obs(is_hrc):
    obs = Mock()
    obs.is_hrc = is_hrc
    obs.ciao = Mock()
    obs.ciao.pget = Mock(side_effect=["10.0", "-20.0", "4096.5", "4096.5"])
    return obs


def _dmcopy_circle_radius(ciao_mock):
    (dmcopy_call,) = [c for c in ciao_mock.call_args_list if c.args[0] == "dmcopy"]
    match = CIRCLE_RE.search(dmcopy_call.args[1])
    return float(match.group(1))


@pytest.mark.parametrize(
    "is_hrc,arcsec_per_pixel",
    [(True, HRC_ARCSEC_PER_PIXEL), (False, ACIS_ARCSEC_PER_PIXEL)],
)
def test_filter_events_pixel_scale(is_hrc, arcsec_per_pixel):
    obs = _fake_obs(is_hrc)
    inputs = {"events": ["evt2.fits"], "fov": ["fov1.fits"]}
    outputs = {"events": "evt2_filtered.fits.gz"}

    observation.filter_events.func(obs, inputs, outputs)

    radius = 180  # matches the fixed radius in filter_events
    assert _dmcopy_circle_radius(obs.ciao) == pytest.approx(radius / arcsec_per_pixel)


@pytest.mark.parametrize(
    "is_hrc,arcsec_per_pixel",
    [(True, HRC_ARCSEC_PER_PIXEL), (False, ACIS_ARCSEC_PER_PIXEL)],
)
def test_filter_sources_pixel_scale(is_hrc, arcsec_per_pixel):
    obs = _fake_obs(is_hrc)
    inputs = {"events": "evt2_filtered.fits.gz", "src": "baseline.src"}
    outputs = {"src": "filtered.src"}

    observation.filter_sources.func(obs, inputs, outputs)

    radius = 180  # matches the fixed radius in filter_sources
    assert _dmcopy_circle_radius(obs.ciao) == pytest.approx(radius / arcsec_per_pixel)
