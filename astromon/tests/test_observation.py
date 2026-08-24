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


def _offset_peak_events(seed_yag=0.0, seed_zag=0.0, peak_offset=3.0, n_peak=400):
    """Events with a tight cluster offset from the seed, plus a sparse background.

    The cluster sits `peak_offset` arcsec in +y_angle from the seed, i.e. near the
    edge of a box_size=4 window centred on the seed but at the centre of one
    centred on the peak.
    """
    rng = np.random.default_rng(20260821)
    peak_y = rng.normal(seed_yag + peak_offset, 0.3, n_peak)
    peak_z = rng.normal(seed_zag, 0.3, n_peak)
    bkg_y = rng.uniform(seed_yag - 8, seed_yag + 8, 200)
    bkg_z = rng.uniform(seed_zag - 8, seed_zag + 8, 200)
    return table.Table(
        {
            "y_angle": np.concatenate([peak_y, bkg_y]),
            "z_angle": np.concatenate([peak_z, bkg_z]),
        }
    )


def test_seed_and_select_events_without_peak_uses_source_position():
    """gaussian_detect (seed_from_peak=False) keeps the celldetect seed and its box."""
    events = _offset_peak_events()
    source = {"y_angle": 0.0, "z_angle": 0.0}

    fit_source, sel = observation._seed_and_select_events(
        np.asarray(events["y_angle"]),
        np.asarray(events["z_angle"]),
        source,
        box_size=4,
        seed_from_peak=False,
    )

    assert fit_source is source
    assert np.all(np.abs(events["y_angle"][sel] - 0.0) < 4)
    assert np.all(np.abs(events["z_angle"][sel] - 0.0) < 4)


def test_seed_and_select_events_centres_box_on_peak():
    """peak_gaussian_detect centres the event box on the peak it seeds from.

    Selecting events around the celldetect position while seeding the fit from an
    offset peak leaves the fit window asymmetric about its own seed: events on the
    far side of the peak are missing, which biases the centroid back toward the
    celldetect position the re-seed was meant to escape.
    """
    events = _offset_peak_events(peak_offset=3.0)
    source = {"y_angle": 0.0, "z_angle": 0.0}
    box_size = 4

    fit_source, sel = observation._seed_and_select_events(
        np.asarray(events["y_angle"]),
        np.asarray(events["z_angle"]),
        source,
        box_size=box_size,
        seed_from_peak=True,
    )

    # The seed moved to the peak.
    assert fit_source["y_angle"] == pytest.approx(3.0, abs=0.5)
    assert fit_source["z_angle"] == pytest.approx(0.0, abs=0.5)

    # The event window is centred on that same seed, symmetric about it.
    selected_y = np.asarray(events["y_angle"])[sel]
    selected_z = np.asarray(events["z_angle"])[sel]
    assert np.all(np.abs(selected_y - fit_source["y_angle"]) < box_size)
    assert np.all(np.abs(selected_z - fit_source["z_angle"]) < box_size)

    # And it reaches past the peak on the far side, which a box anchored to the
    # celldetect position could not do.
    assert selected_y.max() > 4.0


def test_seed_and_select_events_preserves_other_source_columns():
    """Re-seeding copies the source rather than mutating it, keeping other columns."""
    events = _offset_peak_events()
    source = {"y_angle": 0.0, "z_angle": 0.0, "COMPONENT": 7, "RA": 83.8}

    fit_source, _ = observation._seed_and_select_events(
        np.asarray(events["y_angle"]),
        np.asarray(events["z_angle"]),
        source,
        box_size=4,
        seed_from_peak=True,
    )

    assert fit_source["COMPONENT"] == 7
    assert fit_source["RA"] == 83.8
    assert source["y_angle"] == 0.0, "input source must not be mutated"


# --- task download declarations ---------------------------------------------


def test_gaussian_tasks_do_not_declare_a_raw_evt2_download():
    """The Gaussian fits consume derived products, so they must not demand a download.

    Their inputs are the filtered event file and the celldetect source list, both
    produced upstream. Declaring download=["evt2"] made every run call
    obs.download() for the raw event file they never open -- which fails outright
    on a working directory that cleanup_downloads has pruned, and needs arc5gl or
    CIAO to satisfy a dependency that does not exist. filter_events is the task
    that actually reads the raw file, and it declares the download itself.
    """
    from astromon.task import TASKS

    for name in ("gaussian_detect", "peak_gaussian_detect"):
        task = TASKS.tasks[name]
        declared = set(task._inputs.values())
        assert not any("_evt2.fits" in v and "filtered" not in v for v in declared), (
            f"{name} unexpectedly reads the raw evt2"
        )
        assert not task._download, (
            f"{name} declares download={task._download} but reads only derived"
            " products; filter_events declares the raw evt2 download"
        )


def test_filter_events_still_declares_the_raw_evt2_download():
    """The task that does read the raw file must keep asking for it."""
    from astromon.task import TASKS

    assert "evt2" in TASKS.tasks["filter_events"]._download


# ---- dropping crowded seeds before a gaussian fit is attempted ----


def _celldetect_pair(y1, z1, y2, z2, *, component1=1, component2=2):
    """Two celldetect sources, as the .src file gaussian_detect reads."""
    return table.Table(
        {
            "COMPONENT": np.array([component1, component2], dtype=np.int32),
            "y_angle": np.array([y1, y2]),
            "z_angle": np.array([z1, z2]),
        }
    )


def test_drop_crowded_seeds_removes_a_close_pair():
    """obsid 7263 in miniature: two seeds 3.71" apart, both must go."""
    sources = _celldetect_pair(-11.52, 5.77, -15.22, 5.55)

    result = observation._drop_crowded_seeds(sources)

    assert len(result) == 0


def test_drop_crowded_seeds_keeps_isolated_sources():
    """Two seeds well past the crowding radius survive untouched."""
    sources = _celldetect_pair(0.0, 0.0, 100.0, 0.0)

    result = observation._drop_crowded_seeds(sources)

    assert sorted(result["COMPONENT"].tolist()) == [1, 2]


def test_drop_crowded_seeds_boundary_matches_simple_cross_match():
    """Exactly at the threshold is dropped, matching simple_cross_match's own
    `near_neighbor_dist > near_neighbor_dist` cut: a value equal to the
    threshold does not satisfy a strict ">" there either, so a source sitting
    exactly on the boundary is excluded downstream regardless -- the pre-filter
    must agree, or it would keep something selection then throws away anyway."""
    sources = _celldetect_pair(0.0, 0.0, 6.0, 0.0)

    result = observation._drop_crowded_seeds(sources)

    assert len(result) == 0

    just_outside = _celldetect_pair(0.0, 0.0, 6.001, 0.0)
    assert len(observation._drop_crowded_seeds(just_outside)) == 2


def test_drop_crowded_seeds_handles_a_single_source():
    """One source has no neighbour at all, so it is never crowded."""
    sources = table.Table(
        {
            "COMPONENT": np.array([1], dtype=np.int32),
            "y_angle": np.array([0.0]),
            "z_angle": np.array([0.0]),
        }
    )

    result = observation._drop_crowded_seeds(sources)

    assert len(result) == 1


def test_drop_crowded_seeds_handles_empty_input():
    sources = table.Table(
        {
            "COMPONENT": np.array([], dtype=np.int32),
            "y_angle": np.array([]),
            "z_angle": np.array([]),
        }
    )

    result = observation._drop_crowded_seeds(sources)

    assert len(result) == 0


def test_drop_crowded_seeds_three_sources_only_the_close_pair_goes():
    """A third, isolated source in the same field must not be affected."""
    sources = table.Table(
        {
            "COMPONENT": np.array([1, 2, 3], dtype=np.int32),
            "y_angle": np.array([0.0, 3.0, 100.0]),
            "z_angle": np.array([0.0, 0.0, 0.0]),
        }
    )

    result = observation._drop_crowded_seeds(sources)

    assert sorted(result["COMPONENT"].tolist()) == [3]
