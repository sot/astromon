import re
from unittest.mock import Mock, patch

import numpy as np
import pytest
from astropy import table

from astromon import observation

HRC_ARCSEC_PER_PIXEL = 0.13180
ACIS_ARCSEC_PER_PIXEL = 0.5

CIRCLE_RE = re.compile(r"circle\([^,]+,[^,]+,([^)]+)\)")


def _sources_table(**extra_cols):
    """Minimal sources table as returned by _get_sources() after renaming.

    Includes all columns _get_sources() always produces. Pass extra_cols to add
    columns only present for certain detect versions (e.g. psfratio, concentration_ratio).
    """
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

    with patch.object(
        observation.TASKS.tasks["celldetect"], "get_result", return_value=None
    ):
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

    with patch.object(
        observation.TASKS.tasks["gaussian_detect"], "get_result", return_value=None
    ):
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


# ---------------------------------------------------------------------------
# _flag_brightest_source tests
# ---------------------------------------------------------------------------


def test_flag_brightest_source_marks_max_snr():
    """The source with the highest SNR gets brightest=True, all others False."""
    snr = np.array([3.0, 10.0, 5.0, 7.0], dtype=np.float32)
    result = observation._flag_brightest_source(snr)
    assert result.tolist() == [False, True, False, False]


def test_flag_brightest_source_single_source():
    """A single-source observation always gets brightest=True."""
    snr = np.array([6.5], dtype=np.float32)
    result = observation._flag_brightest_source(snr)
    assert result.tolist() == [True]


def test_flag_brightest_source_all_nan_returns_all_false():
    """When all SNR values are NaN, no source is marked brightest."""
    snr = np.array([float("nan"), float("nan"), float("nan")])
    result = observation._flag_brightest_source(snr)
    assert not any(result)
    assert len(result) == 3


def test_flag_brightest_source_ignores_nan_picks_max_finite():
    """NaN values are skipped; the max among finite values is chosen."""
    snr = np.array([float("nan"), 4.0, float("nan"), 9.0, 2.0])
    result = observation._flag_brightest_source(snr)
    assert result.tolist() == [False, False, False, True, False]


def test_flag_brightest_source_uses_uppercase_snr_column():
    """_flag_brightest_source works on a raw celldetect-style column (uppercase SNR)
    by accepting the column values directly — the caller resolves the column name."""
    # In _get_sources the caller does: sources["SNR"] or sources["snr"]
    # _flag_brightest_source just gets the numeric values, so this is a passthrough test.
    snr = np.array([15.0, 3.0, 8.0])
    result = observation._flag_brightest_source(snr)
    assert result[0]  # index 0 has the highest SNR


# --- peak_offset ------------------------------------------------------------
#
# peak_offset records how far the peak-seeded Gaussian fit landed from this
# row's own position. It replaces storing peak_gaussian_detect as a separate
# detect_method: the peak position is a stability diagnostic, not a competing
# answer, so a scalar per row carries the signal without a duplicate table.


def _peak_src(path, components, ras, decs):
    """Write a minimal peak_gaussian_detect .src file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    table.Table(
        {
            "COMPONENT": np.array(components, dtype=np.int32),
            "RA": np.array(ras, dtype=float),
            "DEC": np.array(decs, dtype=float),
        }
    ).write(path, format="fits", overwrite=True)


def _obs_for_peak(tmp_path):
    return observation.Observation(
        7001, workdir=tmp_path / "work", archive_dir=tmp_path / "arch", use_ciao=False
    )


def _sources_with_components(components, ras, decs):
    return table.Table(
        {
            "COMPONENT": np.array(components, dtype=np.int32),
            "RA": np.array(ras, dtype=float),
            "DEC": np.array(decs, dtype=float),
        }
    )


def test_peak_offset_is_zero_when_the_peak_fit_coincides(tmp_path):
    obs = _obs_for_peak(tmp_path)
    _peak_src(
        obs.workdir / "sources" / "7001_peak_gaussian_detect.src",
        [1, 2],
        [150.0, 151.0],
        [20.0, 21.0],
    )
    src = _sources_with_components([1, 2], [150.0, 151.0], [20.0, 21.0])
    offset = obs._peak_offset(src)
    assert np.allclose(offset, 0.0, atol=1e-6)


def test_peak_offset_matches_a_known_separation(tmp_path):
    """A 2 arcsec offset in declination must read back as 2 arcsec."""
    obs = _obs_for_peak(tmp_path)
    _peak_src(
        obs.workdir / "sources" / "7001_peak_gaussian_detect.src",
        [1],
        [150.0],
        [20.0 + 2.0 / 3600.0],
    )
    src = _sources_with_components([1], [150.0], [20.0])
    offset = obs._peak_offset(src)
    assert offset[0] == pytest.approx(2.0, abs=1e-3)


def test_peak_offset_is_nan_when_the_peak_file_is_absent(tmp_path):
    """NaN, not zero: an unmeasured offset must not look like perfect agreement."""
    obs = _obs_for_peak(tmp_path)
    src = _sources_with_components([1, 2], [150.0, 151.0], [20.0, 21.0])
    offset = obs._peak_offset(src)
    assert np.all(np.isnan(offset))


def test_peak_offset_is_nan_for_a_component_absent_from_the_peak_fit(tmp_path):
    """The gaussian fit can fail on a source celldetect found; that row gets NaN."""
    obs = _obs_for_peak(tmp_path)
    _peak_src(
        obs.workdir / "sources" / "7001_peak_gaussian_detect.src", [1], [150.0], [20.0]
    )
    src = _sources_with_components([1, 2], [150.0, 151.0], [20.0, 21.0])
    offset = obs._peak_offset(src)
    assert offset[0] == pytest.approx(0.0, abs=1e-6)
    assert np.isnan(offset[1])


def test_peak_offset_matches_on_component_not_row_order(tmp_path):
    """COMPONENT is the shared key across all three .src files; order must not matter.

    The offset is applied in declination so the separation is exactly 1 arcsec --
    an RA offset would subtend 1" * cos(dec).
    """
    obs = _obs_for_peak(tmp_path)
    _peak_src(
        obs.workdir / "sources" / "7001_peak_gaussian_detect.src",
        [2, 1],
        [151.0, 150.0],
        [21.0, 20.0 + 1.0 / 3600.0],
    )
    src = _sources_with_components([1, 2], [150.0, 151.0], [20.0, 21.0])
    offset = obs._peak_offset(src)
    assert offset[0] == pytest.approx(1.0, abs=1e-2), (
        "component 1 should pair with component 1"
    )
    assert offset[1] == pytest.approx(0.0, abs=1e-6)


def test_peak_offset_in_xray_src_schema():
    from astromon import db

    assert "peak_offset" in db.ASTROMON_XRAY_SRC_DTYPE.names
    assert db.ASTROMON_XRAY_SRC_DTYPE["peak_offset"].kind == "f", (
        "must be a float so the missing-value fill is NaN rather than a real-looking 0"
    )
