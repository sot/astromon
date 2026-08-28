import contextlib
import os
import time
from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np
import pytest
import requests
from astropy import coordinates as coords
from astropy import units as u
from astropy.table import Table, vstack
from cxotime import CxoTime
from Quaternion import Quat
from ska_helpers.retry import retry
from testr.test_helper import has_internet

from astromon import cross_match, db

DATA_DIR = Path(__file__).parent / "data"

HAS_INTERNET = has_internet()

# Vizier queries occasionally hit a transient connection reset with no indication of an actual
# problem with the query itself (observed as requests.exceptions.ConnectionError). Retry those
# rather than failing the test outright.
retry_on_connection_error = retry(
    exceptions=(requests.exceptions.ConnectionError, requests.exceptions.Timeout),
    tries=3,
    delay=1,
    backoff=2,
)


@contextlib.contextmanager
def skip_on_gaia_outage():
    """Turn a Gaia archive outage inside the block into a skip.

    get_gaia_var_stars and get_milliquas_gaia deliberately degrade gracefully: they log a
    warning and return an empty table when a TAP query fails. That means an unreachable
    Gaia archive reaches the test as "expected at least one source, got 0" — a science
    failure rather than an outage. Record any TAP exception and skip on it instead.
    """
    tap_errors = []
    real_tap_query = cross_match._execute_gaia_tap_query

    def recording_tap_query(query):
        try:
            return real_tap_query(query)
        except Exception as exc:
            tap_errors.append(exc)
            raise

    with patch.object(
        cross_match, "_execute_gaia_tap_query", side_effect=recording_tap_query
    ):
        yield

    if tap_errors:
        pytest.skip(f"Gaia archive unreachable: {tap_errors[0]!r}")


def _get_table(name, *args, **kawargs):
    filename = DATA_DIR / f"{name}.ecsv"
    assert filename.exists()
    return Table.read(filename)


def test_db_agreement(monkeypatch):
    monkeypatch.setattr(db, "get_table", _get_table)

    matches = cross_match.compute_cross_matches("astromon_21")
    matches2 = db.get_cross_matches()
    # the cross_match method can potentially return more columns, that's ok.
    assert not [m for m in matches2.colnames if m not in matches.colnames]
    assert len(matches) > 0
    assert len(matches) == len(matches2)
    assert np.all(matches2["category"] == matches["category"])


def test_default_cross_match(monkeypatch):
    monkeypatch.setattr(db, "get_table", _get_table)

    # cross_match without arguments uses the default standard cross_match
    matches2 = cross_match.compute_cross_matches()
    matches = cross_match.compute_cross_matches("astromon_21")

    assert matches2.colnames == matches.colnames
    assert len(matches) > 0
    assert len(matches) == len(matches2)
    assert np.all(matches2["category"] == matches["category"])


def test_custom_args(monkeypatch):
    monkeypatch.setattr(db, "get_table", _get_table)

    # cross_match with custom arguments uses the "simple" algorithm
    matches = cross_match.compute_cross_matches(
        catalogs=["Tycho2"],
        snr=5,
        r_angle=120.0,
        r_angle_grating=0.0,
        near_neighbor_dist=6.0,
        dr=3.0,
    )
    assert len(matches) == 4


# Real, known cross-matches to check against live Vizier queries, one per catalog in
# cross_match.VIZIER_CATALOGS. Each entry in cross_match.VIZIER_CATALOGS maps 'ra'/'dec'/'mag'
# to specific column names in the raw astroquery/Vizier response, and those names are neither
# obvious nor guaranteed to stay the same: astroquery > 0.4.8 changed the separator used in
# proper-motion-corrected column names (e.g. '_RAJ2000_2020.500' -> '_RAJ2000/2020.500') and
# renamed the 2MASS name column ('_2MASS' -> '2MASS'). A wrong column name for 'ra'/'dec'/'mag'
# does not raise an exception: cross_match.get_vizier silently fills the column with masked
# values instead (see astromon/cross_match.py), so a regression here is easy to miss without a
# test that checks actual values.
#
# The Tycho2/2MASS/USNO-B1.0/UCAC4 entries below are the real cross-matches for OBSID 13276
# recorded in astromon/tests/data/astromon_cat_src.ecsv. HIP/HIP2 use Sirius, and ICRS uses
# 3C 273 (a defining ICRF source), since those catalogs have no match in that field.
#
# 'search_ra'/'search_dec' are passed to the query. 'ra'/'dec' are the expected result: these
# differ for HIP/HIP2 because Sirius has a large enough proper motion (~1.3 arcsec/year) that its
# J2000 catalog position (the search center) and its position propagated to VIZIER_QUERY_TIME
# (the expected 'ra'/'dec' column values) are about 17 arcsec apart.
KNOWN_VIZIER_SOURCES = {
    "Tycho2": {
        "search_ra": 258.0868754233,
        "search_dec": -38.4918276755,
        "radius": 2 * u.arcsec,
        "name": "7870-965-1",
        "ra": 258.0868754233,
        "dec": -38.4918276755,
        "mag": 10.982,
    },
    "2MASS": {
        "search_ra": 258.044673,
        "search_dec": -38.507812,
        "radius": 2 * u.arcsec,
        "name": "17121072-3830281",
        "ra": 258.044673,
        "dec": -38.507812,
        "mag": 13.087,
    },
    "USNO-B1.0": {
        "search_ra": 258.120025,
        "search_dec": -38.471181,
        "radius": 2 * u.arcsec,
        "name": "0515-0502275",
        "ra": 258.120025,
        "dec": -38.471181,
        "mag": 15.13,
    },
    "UCAC4": {
        "search_ra": 258.084256925,
        "search_dec": -38.498867073,
        "radius": 2 * u.arcsec,
        "name": "258-115041",
        "ra": 258.084256925,
        "dec": -38.498867073,
        "mag": 14.871,
    },
    "HIP": {
        "search_ra": 101.28715533,
        "search_dec": -16.71611586,
        "radius": 25 * u.arcsec,
        "name": "32349",
        "ra": 101.2850787107,
        "dec": -16.7205708639,
        "mag": -1.44,
    },
    "HIP2": {
        "search_ra": 101.28715533,
        "search_dec": -16.71611586,
        "radius": 25 * u.arcsec,
        "name": "32349",
        "ra": 101.2850787128,
        "dec": -16.7205708058,
        "mag": -1.0876,
    },
    "RFC": {
        "search_ra": 187.2779155448750,
        "search_dec": 2.0523884103056,
        "radius": 2 * u.arcsec,
        "name": "RFC J1229+0203",
        "ra": 187.2779155448750,
        "dec": 2.0523884103056,
        "mag": None,  # RFC has no magnitude column, mag is expected to stay masked
    },
    "SDSS": {
        "search_ra": 184.99294600,
        "search_dec": 14.99606200,
        "radius": 2 * u.arcsec,
        "name": "J121958.30+145945.8",
        "ra": 184.99294600,
        "dec": 14.99606200,
        "mag": 25.920,
    },
    "Gaia2": {
        "search_ra": 258.0868754233,
        "search_dec": -38.4918276755,
        "radius": 3 * u.arcsec,
        "name": "5973148908477305600",
        "ra": 258.0868421690529,
        "dec": -38.4917874972768,
        "mag": 10.3580,
    },
}

VIZIER_QUERY_TIME = CxoTime("2013-02-11T07:39:09")
QUAT = Quat(equatorial=(150.0, 2.0, 30.0))


@pytest.mark.skipif(not HAS_INTERNET, reason="Requires network access")
@pytest.mark.parametrize("catalog", sorted(KNOWN_VIZIER_SOURCES))
@retry_on_connection_error
def test_vizier_catalog_columns(catalog):
    expected = KNOWN_VIZIER_SOURCES[catalog]
    result = cross_match._get(
        catalog,
        VIZIER_QUERY_TIME,
        ra=expected["search_ra"],
        dec=expected["search_dec"],
        radius=expected["radius"],
    )

    assert len(result) == 1
    assert result["name"][0] == expected["name"]
    assert not np.ma.is_masked(result["ra"][0])
    assert not np.ma.is_masked(result["dec"][0])
    assert np.isclose(result["ra"][0], expected["ra"], atol=1e-4)
    assert np.isclose(result["dec"][0], expected["dec"], atol=1e-4)

    if expected["mag"] is None:
        assert np.ma.is_masked(result["mag"][0])
    else:
        assert not np.ma.is_masked(result["mag"][0])
        assert np.isclose(result["mag"][0], expected["mag"], atol=0.01)


# test_vizier_catalog_columns above checks that cross_match._get() maps columns correctly, but
# nothing exercises rough_match() itself -- the function the production pipeline actually calls
# (see astromon/scripts/get_cat_obs_data.py) -- for any catalog. That matters most for HIP, HIP2
# and Gaia2: the standard cross-matches in CROSS_MATCHES_ARGS use only ICRS/Tycho2, and the
# default call to rough_match() in get_cat_obs_data.py uses rough_match's own default catalog
# list, which also excludes HIP/HIP2/Gaia2. So a bug in how rough_match() itself handles a
# catalog (as opposed to a wrong column name) would surface in production for RFC/Tycho2/
# USNO-B1.0/2MASS/SDSS, but never for HIP/HIP2/Gaia2. This exercises rough_match() end-to-end for
# every catalog in VIZIER_CATALOGS.
@pytest.mark.skipif(not HAS_INTERNET, reason="Requires network access")
@pytest.mark.parametrize("catalog", sorted(KNOWN_VIZIER_SOURCES))
@retry_on_connection_error
def test_rough_match_known_source(catalog):
    expected = KNOWN_VIZIER_SOURCES[catalog]
    sources = Table(
        {"id": [1], "ra": [expected["search_ra"]], "dec": [expected["search_dec"]]}
    )

    result = cross_match.rough_match(
        sources, VIZIER_QUERY_TIME, radius=expected["radius"], catalogs=[catalog]
    )

    assert len(result) == 1
    expected_catalog = (
        catalog
        if catalog not in cross_match.VIZIER_CATALOGS
        else cross_match.VIZIER_CATALOGS[catalog]["catalog"]
    )
    assert result["catalog"][0] == expected_catalog
    assert result["name"][0] == expected["name"]
    assert result["x_id"][0] == 1
    assert result["separation"][0] < expected["radius"].to_value(u.arcsec)


# ---- GaiaAGN tests ----


def _fake_tap_result(source_id, ra, dec, mag=18.5):
    """Minimal astropy Table mimicking a Gaia TAP response."""
    return Table(
        {
            "source_id": np.array([source_id], dtype=np.int64),
            "ra": np.array([ra]),
            "dec": np.array([dec]),
            "phot_g_mean_mag": np.array([mag], dtype=np.float32),
        }
    )


def _fake_cached_gaia_catalog(catalog_name, source_id, ra, dec, mag=18.5):
    """Minimal cached Gaia catalog table with the post-download schema."""
    return Table(
        {
            "name": [f"{catalog_name}-{source_id}"],
            "ra": np.array([ra], dtype=float),
            "dec": np.array([dec], dtype=float),
            "mag": np.array([mag], dtype=np.float32),
        }
    )


def test_get_gaia_agn_empty_pos():
    """Empty position array returns an empty table with the correct schema."""
    empty_pos = coords.SkyCoord([], [], unit="deg")
    result = cross_match.get_gaia_agn(empty_pos)
    assert len(result) == 0
    assert set(cross_match.CROSS_MATCH_DTYPE.names).issubset(set(result.colnames))


def test_get_gaia_agn_mocked_catalog():
    """get_gaia_agn structures output correctly from a mocked cached catalog."""
    pos = coords.SkyCoord([187.2779], [2.0524], unit="deg")
    fake = _fake_cached_gaia_catalog(
        "GaiaAGN", source_id=3726862715688776192, ra=187.2779, dec=2.0524, mag=12.9
    )

    with patch.object(cross_match, "get_gaia_agn_catalog", return_value=fake):
        result = cross_match.get_gaia_agn(pos, radius=3 * u.arcsec)

    assert len(result) == 1
    assert result["catalog"][0] == "GaiaAGN"
    assert result["name"][0] == "GaiaAGN-3726862715688776192"
    assert np.isclose(result["ra"][0], 187.2779)
    assert np.isclose(result["dec"][0], 2.0524)
    assert np.isclose(result["mag"][0], 12.9, atol=0.01)


def test_get_gaia_agn_returns_catalog_rows_once_for_multiple_positions():
    """One cached Gaia source is returned once even if multiple positions overlap it."""
    pos = coords.SkyCoord([187.2779, 187.2780], [2.0524, 2.0525], unit="deg")
    fake = _fake_cached_gaia_catalog(
        "GaiaAGN", source_id=3726862715688776192, ra=187.2779, dec=2.0524
    )

    with patch.object(cross_match, "get_gaia_agn_catalog", return_value=fake):
        result = cross_match.get_gaia_agn(pos, radius=3 * u.arcsec)

    assert len(result) == 1


def test_get_gaia_agn_catalog_failure_returns_empty():
    """A cached-catalog load/download failure returns an empty table without raising."""
    pos = coords.SkyCoord([187.2779], [2.0524], unit="deg")

    with patch.object(
        cross_match,
        "get_gaia_agn_catalog",
        side_effect=RuntimeError("catalog unavailable"),
    ):
        result = cross_match.get_gaia_agn(pos, radius=3 * u.arcsec)

    assert len(result) == 0
    assert set(cross_match.CROSS_MATCH_DTYPE.names).issubset(set(result.colnames))


def _minimal_obs(
    obsid: int = 99900,
    ra: float = 83.82,
    dec: float = -5.39,
    target: str = "test",
    category_id: int = 50,
) -> Table:
    """One-row astromon_obs table for use in cross_match tests."""
    return Table(
        {
            "obsid": [obsid],
            "detector": ["ACIS-S"],
            "target": [target],
            "grating": ["NONE"],
            "sim_z": [-190.143],
            "date_obs": ["2020-01-01T00:00:00"],
            "tstart": [0.0],
            "ascdsver": ["10.0"],
            "ra": [ra],
            "dec": [dec],
            "roll": [0.0],
            "category_id": [category_id],
            "version": ["10.0"],
        }
    )


def _minimal_xray_src(
    obsid: int = 99900,
    ra: float = 83.82,
    dec: float = -5.39,
    snr: float = 10.0,
    net_counts: float = 500.0,
) -> Table:
    """One celldetect X-ray source for use in cross_match tests."""
    return Table(
        {
            "obsid": [obsid],
            "id": [1],
            "ra": [ra],
            "dec": [dec],
            "net_counts": [net_counts],
            "y_angle": [0.0],
            "z_angle": [0.0],
            "r_angle": [10.0],
            "snr": [snr],
            "near_neighbor_dist": [60.0],
            "pileup": [0.0],
            "acis_streak": [False],
            "caldb_version": ["4.10.0"],
            "detect_method": ["celldetect"],
        }
    )


def _minimal_cat_src(
    obsid: int = 99900,
    celldetect_x_id: int = 1,
    catalog: str = "Tycho2",
    name: str = "test-star",
    ra: float = 83.82028,
    dec: float = -5.39,
    mag: float = 8.0,
) -> Table:
    """One catalog source positioned ~1 arcsec from the xray source."""
    return Table(
        {
            "obsid": [obsid],
            "id": [0],
            "celldetect_x_id": [celldetect_x_id],
            "catalog": [catalog],
            "name": [name],
            "ra": [ra],
            "dec": [dec],
            "mag": [mag],
            "y_angle": [0.0],
            "z_angle": [0.0],
            "separation": [1.0],
        }
    )


def test_compute_cross_matches_gaia_agn():
    """compute_cross_matches('gaia_agn') finds a GaiaAGN match from a synthetic cat_src."""
    matches = cross_match.compute_cross_matches(
        "gaia_agn",
        astromon_obs=_minimal_obs(99999, ra=187.2779, dec=2.0524),
        astromon_xray_src=_minimal_xray_src(99999, ra=187.2779, dec=2.0524),
        # GaiaAGN catalog source ~1 arcsec away; x_id=1 links it to xray_src
        # id=1, as simple_cross_match joins on (obsid, x_id).
        astromon_cat_src=_minimal_cat_src(
            99999,
            catalog="GaiaAGN",
            name="GaiaAGN-3726862715688776192",
            ra=187.2782,
            dec=2.0524,
            mag=12.9,
        ),
    )

    assert len(matches) == 1
    assert matches["select_name"][0] == "gaia_agn"
    assert matches["obsid"][0] == 99999
    assert matches["x_id"][0] == 1


# Gaia DR3 source confirmed present in gaiadr3.agn_cross_id (verified via TAP query).
# source_id=1535808699355081728, phot_g_mean_mag=13.38
_AGNCROSS_RA = 182.6357
_AGNCROSS_DEC = 39.4059
_AGNCROSS_SOURCE_ID = 1535808699355081728


# ---- GaiaVarStar tests ----


def _fake_var_tap_result(source_id, ra, dec, pmra=0.0, pmdec=0.0, mag=14.0):
    """Minimal astropy Table mimicking a Gaia TAP response with proper-motion columns."""
    return Table(
        {
            "source_id": np.array([source_id], dtype=np.int64),
            "ra": np.array([ra]),
            "dec": np.array([dec]),
            "pmra": np.array([pmra], dtype=np.float64),
            "pmdec": np.array([pmdec], dtype=np.float64),
            "phot_g_mean_mag": np.array([mag], dtype=np.float32),
        }
    )


def test_get_gaia_var_stars_ruwe_filter_in_query():
    """The TAP query includes the RUWE cut at the requested threshold."""
    pos = coords.SkyCoord(
        [237.252], [61.710], unit="deg", obstime=CxoTime("2020-01-01")
    )
    captured = {}

    def capture(query):
        captured["query"] = query

    with patch.object(cross_match, "_execute_gaia_tap_query", side_effect=capture):
        cross_match.get_gaia_var_stars(pos, radius=3 * u.arcsec, max_ruwe=1.2)

    assert "ruwe" in captured["query"].lower()
    assert "1.2" in captured["query"]


def test_get_gaia_var_stars_empty_pos():
    """Empty position array returns an empty table with the correct schema."""
    empty_pos = coords.SkyCoord([], [], unit="deg")
    result = cross_match.get_gaia_var_stars(empty_pos)
    assert len(result) == 0
    assert set(cross_match.CROSS_MATCH_DTYPE.names).issubset(set(result.colnames))


def test_get_gaia_var_stars_mocked_tap():
    """get_gaia_var_stars structures output correctly from a mocked TAP response."""
    pos = coords.SkyCoord(
        [237.252], [61.710], unit="deg", obstime=CxoTime("2020-01-01")
    )
    fake = _fake_var_tap_result(
        source_id=1627844733203662976,
        ra=237.252158,
        dec=61.710091,
        pmra=40.430,
        pmdec=-31.185,
        mag=13.14,
    )

    with patch.object(cross_match, "_execute_gaia_tap_query", return_value=fake):
        result = cross_match.get_gaia_var_stars(pos, radius=3 * u.arcsec)

    assert len(result) == 1
    assert result["catalog"][0] == "GaiaVarStar"
    assert result["name"][0] == "GaiaVarStar-1627844733203662976"
    assert np.isclose(result["mag"][0], 13.14, atol=0.01)


def test_get_gaia_var_stars_applies_proper_motion():
    """Returned positions are propagated from J2016.0 to the observation epoch."""
    # Use a source with known pmra=3600000 mas/yr (1 deg/yr) for easy arithmetic.
    # At dec=0, over 2 years from J2016 to J2018: delta_ra = +2 deg exactly.
    obs_epoch = CxoTime("2018-01-01T00:00:00")
    pos = coords.SkyCoord([10.0], [0.0], unit="deg", obstime=obs_epoch)
    fake = _fake_var_tap_result(
        source_id=99,
        ra=10.0,
        dec=0.0,
        pmra=3_600_000.0,
        pmdec=0.0,  # 1 deg/yr in RA at dec=0
    )

    delta_yr = obs_epoch.jyear - 2016.0
    expected_ra = 10.0 + 1.0 * delta_yr  # pmra/3600000 deg/yr * delta_yr yr

    with patch.object(cross_match, "_execute_gaia_tap_query", return_value=fake):
        result = cross_match.get_gaia_var_stars(pos, radius=3 * u.arcsec)

    assert len(result) == 1
    assert np.isclose(result["ra"][0], expected_ra, atol=1e-6)
    assert np.isclose(result["dec"][0], 0.0, atol=1e-9)


def test_get_gaia_var_stars_tap_failure_returns_empty():
    """A TAP failure returns an empty table without raising."""
    pos = coords.SkyCoord(
        [237.252], [61.710], unit="deg", obstime=CxoTime("2020-01-01")
    )

    with patch.object(
        cross_match,
        "_execute_gaia_tap_query",
        side_effect=RuntimeError("TAP unavailable"),
    ):
        result = cross_match.get_gaia_var_stars(pos, radius=3 * u.arcsec)

    assert len(result) == 0
    assert set(cross_match.CROSS_MATCH_DTYPE.names).issubset(set(result.colnames))


def test_compute_cross_matches_gaia_var_star():
    """compute_cross_matches('gaia_var_star') finds a GaiaVarStar match from a synthetic cat_src."""
    matches = cross_match.compute_cross_matches(
        "gaia_var_star",
        astromon_obs=_minimal_obs(
            99998, ra=237.252, dec=61.710, target="test active star", category_id=10
        ),
        astromon_xray_src=_minimal_xray_src(
            99998, ra=237.252, dec=61.710, snr=8.0, net_counts=300.0
        ),
        astromon_cat_src=_minimal_cat_src(
            99998,
            catalog="GaiaVarStar",
            name="GaiaVarStar-1627844733203662976",
            ra=237.2523,
            dec=61.7101,
            mag=13.14,
        ),
    )

    assert len(matches) == 1
    assert matches["select_name"][0] == "gaia_var_star"
    assert matches["obsid"][0] == 99998
    assert matches["x_id"][0] == 1


# Chandra obsid=26425 targeted '2RXS J163341.8-093307' on 2022-05-25 — a ROSAT-selected
# active star.  Gaia source_id=4338433743723131648 (G=10.52) is confirmed in
# gaiadr3.vari_rotation_modulation.  pmra=-64.9, pmdec=-177.9 mas/yr produces a ~1.2 arcsec
# total PM shift from J2016.0 to the observation epoch, large enough to verify the correction.
_RXSJ163_TARGET = "2RXS J163341.8-093307"
_RXSJ163_OBS_DATE = "2022-05-25T23:14:17"
_RXSJ163_SOURCE_ID = 4338433743723131648
_RXSJ163_GAIA_RA = 248.42308  # J2016.0
_RXSJ163_GAIA_DEC = -9.55411  # J2016.0
_RXSJ163_PMRA = -64.9  # mas/yr
_RXSJ163_PMDEC = -177.9  # mas/yr


@pytest.mark.skipif(not HAS_INTERNET, reason="Requires network access")
def test_get_gaia_var_stars_active_star_target():
    """get_gaia_var_stars finds the vari_rotation_modulation counterpart of a ROSAT/Chandra target.

    '2RXS J163341.8-093307' (Chandra obsid=26425) is an X-ray-selected active star.
    Its Gaia counterpart has pmra=-64.9, pmdec=-177.9 mas/yr, so the PM correction
    shifts the position by ~1.2 arcsec over the 6.4 years between J2016.0 and the
    2022 observation — large enough to clearly verify the correction is applied.
    """
    obs_time = CxoTime(_RXSJ163_OBS_DATE)
    pos = coords.SkyCoord(
        [_RXSJ163_GAIA_RA], [_RXSJ163_GAIA_DEC], unit="deg", obstime=obs_time
    )

    with skip_on_gaia_outage():
        result = cross_match.get_gaia_var_stars(pos, radius=5 * u.arcsec)

    assert len(result) >= 1, (
        f"Expected at least one vari_rotation_modulation source near "
        f"'{_RXSJ163_TARGET}' (ra={_RXSJ163_GAIA_RA}, dec={_RXSJ163_GAIA_DEC})"
    )
    source_ids = [int(n.split("-")[1]) for n in result["name"]]
    assert _RXSJ163_SOURCE_ID in source_ids, (
        f"source_id={_RXSJ163_SOURCE_ID} for target '{_RXSJ163_TARGET}' not found; "
        f"got source_ids={source_ids}"
    )
    idx = source_ids.index(_RXSJ163_SOURCE_ID)

    # Returned position should be PM-corrected away from the raw J2016.0 catalog position.
    # pmdec=-177.9 mas/yr over ~6.4 years gives ~1.1 arcsec southward shift — clearly detectable.
    dec_shift_arcsec = (result["dec"][idx] - _RXSJ163_GAIA_DEC) * 3600
    assert dec_shift_arcsec < -0.5, (
        f"Expected southward Dec shift > 0.5 arcsec from PM correction, got {dec_shift_arcsec:.3f} arcsec"
    )


# ---- GaiaQSO tests ----

# 3C 273 is the closest bright quasar (z=0.158, G=12.84), confirmed in
# gaiadr3.qso_candidates but NOT in gaiadr3.agn_cross_id — this specific pair
# demonstrates the extra coverage that get_gaia_qso_candidates provides over get_gaia_agn.
_3C273_RA = 187.27792
_3C273_DEC = 2.05239
_3C273_QSO_SOURCE_ID = 3700386905605055360


def test_get_gaia_qso_candidates_empty_pos():
    """Empty position array returns an empty table with the correct schema."""
    empty_pos = coords.SkyCoord([], [], unit="deg")
    result = cross_match.get_gaia_qso_candidates(empty_pos)
    assert len(result) == 0
    assert set(cross_match.CROSS_MATCH_DTYPE.names).issubset(set(result.colnames))


def test_get_gaia_qso_candidates_mocked_catalog():
    """get_gaia_qso_candidates structures output correctly from a mocked catalog."""
    pos = coords.SkyCoord([_3C273_RA], [_3C273_DEC], unit="deg")
    fake = _fake_cached_gaia_catalog(
        "GaiaQSO",
        source_id=_3C273_QSO_SOURCE_ID,
        ra=_3C273_RA,
        dec=_3C273_DEC,
        mag=12.84,
    )

    with patch.object(cross_match, "get_gaia_qso_catalog", return_value=fake):
        result = cross_match.get_gaia_qso_candidates(pos, radius=3 * u.arcsec)

    assert len(result) == 1
    assert result["catalog"][0] == "GaiaQSO"
    assert result["name"][0] == f"GaiaQSO-{_3C273_QSO_SOURCE_ID}"
    assert np.isclose(result["ra"][0], _3C273_RA)
    assert np.isclose(result["dec"][0], _3C273_DEC)
    assert np.isclose(result["mag"][0], 12.84, atol=0.01)


def test_get_gaia_qso_candidates_returns_catalog_rows_once_for_multiple_positions():
    """One cached QSO source is returned once even if multiple positions overlap it."""
    pos = coords.SkyCoord(
        [_3C273_RA, _3C273_RA + 0.001], [_3C273_DEC, _3C273_DEC], unit="deg"
    )
    fake = _fake_cached_gaia_catalog(
        "GaiaQSO", source_id=_3C273_QSO_SOURCE_ID, ra=_3C273_RA, dec=_3C273_DEC
    )

    with patch.object(cross_match, "get_gaia_qso_catalog", return_value=fake):
        result = cross_match.get_gaia_qso_candidates(pos, radius=3 * u.arcsec)

    assert len(result) == 1


def test_get_gaia_qso_candidates_catalog_failure_returns_empty():
    """A cached-catalog load/download failure returns an empty table without raising."""
    pos = coords.SkyCoord([_3C273_RA], [_3C273_DEC], unit="deg")

    with patch.object(
        cross_match,
        "get_gaia_qso_catalog",
        side_effect=RuntimeError("catalog unavailable"),
    ):
        result = cross_match.get_gaia_qso_candidates(pos, radius=3 * u.arcsec)

    assert len(result) == 0
    assert set(cross_match.CROSS_MATCH_DTYPE.names).issubset(set(result.colnames))


def test_compute_cross_matches_gaia_qso():
    """compute_cross_matches('gaia_qso') finds a GaiaQSO match from a synthetic cat_src."""
    matches = cross_match.compute_cross_matches(
        "gaia_qso",
        astromon_obs=_minimal_obs(99997, ra=_3C273_RA, dec=_3C273_DEC, target="3C 273"),
        astromon_xray_src=_minimal_xray_src(
            99997, ra=_3C273_RA, dec=_3C273_DEC, snr=20.0, net_counts=1000.0
        ),
        astromon_cat_src=_minimal_cat_src(
            99997,
            catalog="GaiaQSO",
            name=f"GaiaQSO-{_3C273_QSO_SOURCE_ID}",
            ra=_3C273_RA + 0.0003,  # ~1 arcsec offset
            dec=_3C273_DEC,
            mag=12.84,
        ),
    )

    assert len(matches) == 1
    assert matches["select_name"][0] == "gaia_qso"
    assert matches["obsid"][0] == 99997
    assert matches["x_id"][0] == 1


# ---- Quaia tests ----

# PKS 0537-441 is a well-known FSRQ blazar at z=0.894, confirmed present in the
# Quaia quasar catalog (Gaia DR3 source_id=4802492418750413568, sep~0.8 arcsec).
# G~17 puts it squarely within the Quaia G<20 selection (Gaia XP spectra).
_PKS0537_RA = 84.70984
_PKS0537_DEC = -44.08582
_PKS0537_QUAIA_NAME = "Gaia DR3 4802492418750413568"


def test_get_quaia_candidates_empty_pos():
    """Empty position array returns an empty table with the correct schema."""
    empty_pos = coords.SkyCoord([], [], unit="deg")
    result = cross_match.get_quaia_candidates(empty_pos)
    assert len(result) == 0
    assert set(cross_match.CROSS_MATCH_DTYPE.names).issubset(set(result.colnames))


def test_get_quaia_candidates_mocked_catalog():
    """get_quaia_candidates filters a mocked local catalog correctly."""
    pos = coords.SkyCoord([_PKS0537_RA], [_PKS0537_DEC], unit="deg")
    fake_catalog = Table(
        {
            "name": [_PKS0537_QUAIA_NAME, "Gaia DR3 9999999999999999999"],
            "ra": [_PKS0537_RA, _PKS0537_RA + 1.0],  # second source 1 deg away
            "dec": [_PKS0537_DEC, _PKS0537_DEC],
        }
    )

    with patch.object(cross_match, "get_quaia", return_value=fake_catalog):
        result = cross_match.get_quaia_candidates(pos, radius=3 * u.arcsec)

    assert len(result) == 1
    assert result["catalog"][0] == "Quaia"
    assert result["name"][0] == _PKS0537_QUAIA_NAME
    assert np.isclose(result["ra"][0], _PKS0537_RA, atol=1e-4)
    assert np.isclose(result["dec"][0], _PKS0537_DEC, atol=1e-4)


@pytest.mark.skipif(not HAS_INTERNET, reason="Requires network access")
def test_get_quaia_candidates_pks0537():
    """get_quaia_candidates finds PKS 0537-441 in the local Quaia catalog.

    PKS 0537-441 (FSRQ blazar, z=0.894, G~17) has a Quaia entry at
    source_id=4802492418750413568 within ~0.8 arcsec of the known position.
    Quaia uses Gaia XP spectra (G<20) and is stored as a local FITS cache.
    """
    pos = coords.SkyCoord([_PKS0537_RA], [_PKS0537_DEC], unit="deg")
    result = cross_match.get_quaia_candidates(pos, radius=3 * u.arcsec)

    assert len(result) >= 1, (
        f"Expected at least one Quaia entry near PKS 0537-441 "
        f"(ra={_PKS0537_RA}, dec={_PKS0537_DEC})"
    )
    assert _PKS0537_QUAIA_NAME in result["name"], (
        f"'{_PKS0537_QUAIA_NAME}' not found in Quaia results; got names={list(result['name'])}"
    )
    assert result["catalog"][0] == "Quaia"


# ---- MilliquasGaia tests ----

# 3C 273 is the canonical defining entry in the Million Quasar Catalog (Milliquas v8).
# Gaia DR3 source_id=3700386905605055360 is within 1.5 arcsec, satisfying the Gaia
# position-upgrade step in get_milliquas_gaia.  Milliquas name: "3C 273".
_3C273_MILLIQUAS_NAME = "3C 273"


def test_get_milliquas_gaia_empty_pos():
    """Empty position array returns an empty table with the correct schema."""
    empty_pos = coords.SkyCoord([], [], unit="deg")
    result = cross_match.get_milliquas_gaia(empty_pos)
    assert len(result) == 0
    assert set(cross_match.CROSS_MATCH_DTYPE.names).issubset(set(result.colnames))


def test_get_milliquas_gaia_mocked_response():
    """get_milliquas_gaia produces correct output from mocked VizieR + Gaia responses."""
    pos = coords.SkyCoord([_3C273_RA], [_3C273_DEC], unit="deg")

    # Mocked Milliquas VizieR row: spectroscopic type "Q" passes the type filter.
    mq_row = Table(
        {
            "RAJ2000": np.array([_3C273_RA], dtype=float),
            "DEJ2000": np.array([_3C273_DEC], dtype=float),
            "Name": [_3C273_MILLIQUAS_NAME],
            "Type": ["Q"],
            "Rmag": np.ma.MaskedArray([12.9], mask=[False]),
        }
    )
    # Mocked Gaia DR3 row confirming a counterpart within 1.5 arcsec.
    gaia_row = Table(
        {
            "source_id": np.array([_3C273_QSO_SOURCE_ID], dtype=np.int64),
            "ra": np.array([_3C273_RA]),
            "dec": np.array([_3C273_DEC]),
        }
    )

    with (
        patch.object(cross_match, "_query_vizier_region", return_value=[mq_row]),
        patch.object(cross_match, "_execute_gaia_tap_query", return_value=gaia_row),
    ):
        result = cross_match.get_milliquas_gaia(pos, radius=3 * u.arcsec)

    assert len(result) == 1
    assert result["catalog"][0] == "MilliquasGaia"
    assert _3C273_MILLIQUAS_NAME in result["name"][0]
    assert np.isclose(result["ra"][0], _3C273_RA, atol=1e-4)
    assert np.isclose(result["dec"][0], _3C273_DEC, atol=1e-4)


def test_get_milliquas_gaia_type_filter_excludes_photometric():
    """Photometric-only (K-type) Milliquas sources are excluded by the type filter."""
    pos = coords.SkyCoord([_3C273_RA], [_3C273_DEC], unit="deg")

    photometric_row = Table(
        {
            "RAJ2000": np.array([_3C273_RA], dtype=float),
            "DEJ2000": np.array([_3C273_DEC], dtype=float),
            "Name": ["PHOT-CANDIDATE"],
            "Type": ["K"],  # photometric candidate — should be excluded
            "Rmag": np.ma.MaskedArray([18.5], mask=[False]),
        }
    )

    with patch.object(
        cross_match, "_query_vizier_region", return_value=[photometric_row]
    ):
        result = cross_match.get_milliquas_gaia(pos, radius=3 * u.arcsec)

    assert len(result) == 0


# ---------------------------------------------------------------------------
# brightest / acis_streak filter tests for compute_cross_matches
# ---------------------------------------------------------------------------


# ---- DESIV161 tests ----

# Confirmed DESI EDR (V/161) QSO with ZWARN=0 in the COSMOS field.
# TargetID=39627670656389181, OType='QSO', ZWARN=0, in the DESI SV footprint.
# COSMOS area at RA~36.5, Dec~-4.7 is a well-covered DESI EDR tile.
_COSMOS_DESI_RA = 36.48102
_COSMOS_DESI_DEC = -4.65509
_COSMOS_DESI_TARGET_ID = 39627670656389181


def test_get_desi_v161_candidates_empty_pos():
    """Empty position array returns an empty table with the correct schema."""
    empty_pos = coords.SkyCoord([], [], unit="deg")
    result = cross_match.get_desi_v161_candidates(empty_pos)
    assert len(result) == 0
    assert set(cross_match.CROSS_MATCH_DTYPE.names).issubset(set(result.colnames))


def test_get_desi_v161_candidates_mocked_vizier():
    """get_desi_v161_candidates filters OType and ZWARN correctly from a mocked response."""
    pos = coords.SkyCoord([_COSMOS_DESI_RA], [_COSMOS_DESI_DEC], unit="deg")
    fake_rows = Table(
        {
            "RAICRS": np.array(
                [_COSMOS_DESI_RA, _COSMOS_DESI_RA + 0.0001], dtype=float
            ),
            "DEICRS": np.array([_COSMOS_DESI_DEC, _COSMOS_DESI_DEC], dtype=float),
            "TargetID": np.array(
                [_COSMOS_DESI_TARGET_ID, 99999999999999999], dtype=np.int64
            ),
            "OType": ["QSO", "STAR"],  # STAR should be filtered out
            "ZWARN": np.array([0, 0], dtype=int),
            "delChi2": np.array([500.0, 500.0], dtype=float),
        }
    )

    with patch.object(cross_match, "_query_vizier_region", return_value=[fake_rows]):
        result = cross_match.get_desi_v161_candidates(pos, radius=3 * u.arcsec)

    assert len(result) == 1, "STAR-type row should be excluded"
    assert result["catalog"][0] == "DESIV161"
    assert result["name"][0] == f"DESIV161-{_COSMOS_DESI_TARGET_ID}"
    assert np.isclose(result["ra"][0], _COSMOS_DESI_RA, atol=1e-4)
    assert np.isclose(result["dec"][0], _COSMOS_DESI_DEC, atol=1e-4)


def test_get_desi_v161_candidates_zwarn_filter():
    """ZWARN != 0 sources are excluded even if OType is QSO or GALAXY."""
    pos = coords.SkyCoord([_COSMOS_DESI_RA], [_COSMOS_DESI_DEC], unit="deg")
    fake_rows = Table(
        {
            "RAICRS": np.array([_COSMOS_DESI_RA, _COSMOS_DESI_RA], dtype=float),
            "DEICRS": np.array([_COSMOS_DESI_DEC, _COSMOS_DESI_DEC], dtype=float),
            "TargetID": np.array(
                [_COSMOS_DESI_TARGET_ID, 11111111111111111], dtype=np.int64
            ),
            "OType": ["QSO", "GALAXY"],
            "ZWARN": np.array([0, 4], dtype=int),  # second has SMALL_DELTA_CHI2 bit set
            "delChi2": np.array([500.0, 900.0], dtype=float),
        }
    )

    with patch.object(cross_match, "_query_vizier_region", return_value=[fake_rows]):
        result = cross_match.get_desi_v161_candidates(pos, radius=3 * u.arcsec)

    assert len(result) == 1, "ZWARN=4 row should be excluded"
    assert result["name"][0] == f"DESIV161-{_COSMOS_DESI_TARGET_ID}"


def test_dedupe_desi_by_delchi2_same_targetid_keeps_the_better_copy():
    """The same TargetID re-observed in two DESI programs: keep the better one.

    obsid 5722's DESIV161-39633335797940469 is DESI's own repeat targeting of
    one object across programs (e.g. "bright" then "dark"), and both copies pass
    ZWARN=0. Picking by array order (the old np.unique(TargetID) behaviour) can
    keep whichever happens to sort first, discarding a materially better fit --
    delChi2 (redrock's own redshift-confidence margin) reliably tells them apart:
    measured pairs in the live database differ by 3x-20x, never a close call.

    Tested against the helper directly: the two rows differ only in delChi2, so
    nothing in get_desi_v161_candidates' own output (name, ra, dec) could tell
    the winner apart from the loser.

    The worse copy is listed FIRST here, so a test that merely kept "the first
    row per TargetID" would pass by accident; this only passes if delChi2 is
    actually compared.
    """
    desi = Table(
        {
            "RAICRS": np.array([_COSMOS_DESI_RA, _COSMOS_DESI_RA], dtype=float),
            "DEICRS": np.array([_COSMOS_DESI_DEC, _COSMOS_DESI_DEC], dtype=float),
            "TargetID": np.array(
                [_COSMOS_DESI_TARGET_ID, _COSMOS_DESI_TARGET_ID], dtype=np.int64
            ),
            "delChi2": np.array([916.1, 6228.3], dtype=float),  # worse one first
        }
    )

    result = cross_match._dedupe_desi_by_delchi2(desi)

    assert len(result) == 1, "one physical target, one row"
    assert float(result["delChi2"][0]) == pytest.approx(6228.3), (
        "the higher-delChi2 copy must be the one kept"
    )


def test_get_desi_v161_candidates_position_dedup_handles_a_rounding_boundary():
    """Two positions straddling a rounding boundary must still be grouped.

    obsid 3000's real duplicate pair is exactly this shape: ra=190.544187 and
    ra=190.544196 are ~0.033" apart -- well inside the ~0.036" tolerance -- but
    190.544187 rounds to 190.54419 at 5 decimal places while 190.544196 rounds to
    190.54420, landing in different buckets under naive round-then-group. This
    is the failure that let 12 of 122 known duplicates survive the first version
    of the fix.
    """
    pos = coords.SkyCoord([_COSMOS_DESI_RA], [_COSMOS_DESI_DEC], unit="deg")
    fake_rows = Table(
        {
            "RAICRS": np.array([190.544187, 190.544196], dtype=float),
            "DEICRS": np.array([33.284073, 33.284073], dtype=float),
            "TargetID": np.array(
                [2305843020496374046, 39632945912218912], dtype=np.int64
            ),
            "OType": ["GALAXY", "GALAXY"],
            "ZWARN": np.array([0, 0], dtype=int),
            "delChi2": np.array([128268.3, 395323.8], dtype=float),
        }
    )

    with patch.object(cross_match, "_query_vizier_region", return_value=[fake_rows]):
        result = cross_match.get_desi_v161_candidates(pos, radius=3 * u.arcsec)

    assert len(result) == 1, '0.033" apart is the same object, not two'
    assert result["name"][0] == "DESIV161-39632945912218912"


def test_get_desi_v161_candidates_position_dedup_at_desis_own_accuracy_limit():
    """A genuine duplicate 0.0508" apart must merge -- and a real distinct
    object 1.6" away in the same anchor group must not.

    obsid 13994's DESIV161-2305843020093724629 and DESIV161-39628438616675080
    are 0.0508" apart, just over a naive 0.05" cutoff, while
    DESIV161-1092755086376973 in the same group is a genuinely different
    object 1.6" away. get_desi_v161_candidates' own docstring states DESI
    positions are accurate to <~0.1" -- narrower than that overstates the
    precision of a repeat observation and leaves real duplicates like this one
    unmerged.
    """
    pos = coords.SkyCoord([_COSMOS_DESI_RA], [_COSMOS_DESI_DEC], unit="deg")
    fake_rows = Table(
        {
            "RAICRS": np.array([194.915384, 194.915246, 194.915256], dtype=float),
            "DEICRS": np.array([27.954351, 27.953920, 27.953909], dtype=float),
            "TargetID": np.array(
                [1092755086376973, 2305843020093724629, 39628438616675080],
                dtype=np.int64,
            ),
            "OType": ["QSO", "QSO", "QSO"],
            "ZWARN": np.array([0, 0, 0], dtype=int),
            "delChi2": np.array([1000.0, 500.0, 5000.0], dtype=float),
        }
    )

    with patch.object(cross_match, "_query_vizier_region", return_value=[fake_rows]):
        result = cross_match.get_desi_v161_candidates(pos, radius=3 * u.arcsec)

    names = sorted(result["name"])
    assert names == ["DESIV161-1092755086376973", "DESIV161-39628438616675080"], (
        "the distant object survives on its own; the close pair merges to its "
        "higher-delChi2 member"
    )


def test_get_desi_v161_candidates_prefers_higher_delchi2_different_targetid():
    """Two different TargetIDs at the same position: DESI's own repeat targeting
    of one physical object across survey phases, not two real sources.

    This is the obsid 5722 shape: DESIV161-2305843020869671679 (Prog=backup) and
    DESIV161-39633335797940469 (Prog=dark) at an identical position, both
    ZWARN=0. The prior code deduplicated only by TargetID, so both survived as
    separate cat_src rows -- an exact dr tie downstream whose winner depended on
    table order. Position is not exposed as a query column here, so the
    dedup keys on (ra, dec) directly.
    """
    pos = coords.SkyCoord([_COSMOS_DESI_RA], [_COSMOS_DESI_DEC], unit="deg")
    fake_rows = Table(
        {
            "RAICRS": np.array([_COSMOS_DESI_RA, _COSMOS_DESI_RA], dtype=float),
            "DEICRS": np.array([_COSMOS_DESI_DEC, _COSMOS_DESI_DEC], dtype=float),
            "TargetID": np.array(
                [2305843020869671679, 39633335797940469], dtype=np.int64
            ),
            "OType": ["QSO", "QSO"],
            "ZWARN": np.array([0, 0], dtype=int),
            "delChi2": np.array([3517.9, 35610.5], dtype=float),
        }
    )

    with patch.object(cross_match, "_query_vizier_region", return_value=[fake_rows]):
        result = cross_match.get_desi_v161_candidates(pos, radius=3 * u.arcsec)

    assert len(result) == 1, "same physical position, one row"
    assert result["name"][0] == "DESIV161-39633335797940469", (
        "the higher-delChi2 TargetID must be the one kept"
    )


def test_get_desi_v161_candidates_position_dedup_leaves_real_pairs_alone():
    """Two genuinely different, well-separated sources must both survive."""
    pos = coords.SkyCoord([_COSMOS_DESI_RA], [_COSMOS_DESI_DEC], unit="deg")
    fake_rows = Table(
        {
            "RAICRS": np.array([_COSMOS_DESI_RA, _COSMOS_DESI_RA + 0.001], dtype=float),
            "DEICRS": np.array([_COSMOS_DESI_DEC, _COSMOS_DESI_DEC], dtype=float),
            "TargetID": np.array([111, 222], dtype=np.int64),
            "OType": ["QSO", "GALAXY"],
            "ZWARN": np.array([0, 0], dtype=int),
            "delChi2": np.array([500.0, 500.0], dtype=float),
        }
    )

    with patch.object(cross_match, "_query_vizier_region", return_value=[fake_rows]):
        result = cross_match.get_desi_v161_candidates(pos, radius=5 * u.arcsec)

    assert len(result) == 2, '0.001 deg (3.6") apart is not a duplicate position'


def test_get_desi_v161_candidates_filters_before_deduping_by_targetid():
    """A ZWARN failure on one repeat observation must not sink the whole TargetID.

    obsid 234.948429/43.228229 in the live database is exactly this: the same
    TargetID's "dark"-program row has ZWARN=2048 (excluded) while its
    "bright"-program row has ZWARN=0 and the far higher delChi2 (kept). Filtering
    after deduplication could pick the ZWARN=2048 row by table order, discard it,
    and lose the TargetID entirely even though a good observation exists.
    """
    pos = coords.SkyCoord([_COSMOS_DESI_RA], [_COSMOS_DESI_DEC], unit="deg")
    fake_rows = Table(
        {
            "RAICRS": np.array([_COSMOS_DESI_RA, _COSMOS_DESI_RA], dtype=float),
            "DEICRS": np.array([_COSMOS_DESI_DEC, _COSMOS_DESI_DEC], dtype=float),
            "TargetID": np.array(
                [_COSMOS_DESI_TARGET_ID, _COSMOS_DESI_TARGET_ID], dtype=np.int64
            ),
            "OType": ["QSO", "QSO"],
            "ZWARN": np.array([2048, 0], dtype=int),  # bad row listed first
            "delChi2": np.array([109260.6, 287295.6], dtype=float),
        }
    )

    with patch.object(cross_match, "_query_vizier_region", return_value=[fake_rows]):
        result = cross_match.get_desi_v161_candidates(pos, radius=3 * u.arcsec)

    assert len(result) == 1, "the ZWARN=0 observation must survive"
    assert result["name"][0] == f"DESIV161-{_COSMOS_DESI_TARGET_ID}"


def test_get_desi_v161_candidates_vizier_failure_returns_empty():
    """A VizieR failure returns an empty table without raising (graceful degradation)."""
    pos = coords.SkyCoord([_COSMOS_DESI_RA], [_COSMOS_DESI_DEC], unit="deg")

    with patch.object(
        cross_match,
        "_query_vizier_region",
        side_effect=RuntimeError("VizieR unavailable"),
    ):
        result = cross_match.get_desi_v161_candidates(pos, radius=3 * u.arcsec)

    assert len(result) == 0
    assert set(cross_match.CROSS_MATCH_DTYPE.names).issubset(set(result.colnames))


@pytest.mark.skipif(not HAS_INTERNET, reason="Requires network access")
@retry_on_connection_error
def test_get_desi_v161_candidates_cosmos_qso():
    """get_desi_v161_candidates finds a confirmed DESI EDR QSO in the COSMOS field.

    TargetID=39627670656389181 (OType='QSO', ZWARN=0) is confirmed present in
    VizieR V/161 within the COSMOS DESI survey-validation tile.
    """
    pos = coords.SkyCoord([_COSMOS_DESI_RA], [_COSMOS_DESI_DEC], unit="deg")
    result = cross_match.get_desi_v161_candidates(pos, radius=10 * u.arcsec)

    assert len(result) >= 1, (
        f"Expected at least one DESI EDR QSO near COSMOS field "
        f"(ra={_COSMOS_DESI_RA}, dec={_COSMOS_DESI_DEC})"
    )
    target_ids = [int(n.split("-")[1]) for n in result["name"]]
    assert _COSMOS_DESI_TARGET_ID in target_ids, (
        f"TargetID={_COSMOS_DESI_TARGET_ID} not found in DESIV161 results; "
        f"got target_ids={target_ids}"
    )
    assert result["catalog"][0] == "DESIV161"
    # All returned rows must pass the ZWARN=0 filter -- no STAR types.
    for name in result["name"]:
        assert name.startswith("DESIV161-")


def test_acis_streak_source_excluded_without_brightest():
    """A source with acis_streak=True and brightest=False is excluded from cross-matching."""
    obsid = 99901
    xray_src = Table(
        {
            "obsid": [obsid],
            "id": [1],
            "ra": [83.82],
            "dec": [-5.39],
            "net_counts": [500.0],
            "y_angle": [0.0],
            "z_angle": [0.0],
            "r_angle": [10.0],
            "snr": [15.0],
            "near_neighbor_dist": [60.0],
            "pileup": [0.0],
            "acis_streak": [1],  # on streak
            "brightest": [0],  # NOT the brightest — excluded
            "caldb_version": ["4.10.0"],
            "detect_method": ["celldetect"],
        }
    )

    matches = cross_match.compute_cross_matches(
        "astromon_21",
        astromon_obs=_minimal_obs(obsid),
        astromon_xray_src=xray_src,
        astromon_cat_src=_minimal_cat_src(obsid),
    )

    assert len(matches) == 0, "streak source without brightest should be excluded"


def test_acis_streak_source_included_when_brightest():
    """A source with acis_streak=True and brightest=True is included in cross-matching.

    This is the grating calibration target scenario: the dispersed spectrum
    triggers a streak detection, but the zeroth-order source (brightest) should
    still be cross-matched.
    """
    obsid = 99902
    xray_src = Table(
        {
            "obsid": [obsid],
            "id": [1],
            "ra": [83.82],
            "dec": [-5.39],
            "net_counts": [5000.0],
            "y_angle": [0.0],
            "z_angle": [0.0],
            "r_angle": [10.0],
            "snr": [120.0],
            "near_neighbor_dist": [60.0],
            "pileup": [0.0],
            "acis_streak": [1],  # flagged on streak (grating spectrum)
            "brightest": [1],  # BUT it's the brightest — include it
            "caldb_version": ["4.10.0"],
            "detect_method": ["celldetect"],
        }
    )

    matches = cross_match.compute_cross_matches(
        "astromon_21",
        astromon_obs=_minimal_obs(obsid),
        astromon_xray_src=xray_src,
        astromon_cat_src=_minimal_cat_src(obsid),
    )

    assert len(matches) == 1, (
        "streak source with brightest=True should be cross-matched"
    )
    assert matches["obsid"][0] == obsid
    assert matches["x_id"][0] == 1


def test_brightest_backward_compat_missing_column():
    """When 'brightest' column is absent, acis_streak=True sources are excluded (safe default).

    Ensures that old DB dumps without the brightest column still work: the
    cross-match treats every acis_streak source as non-brightest.
    """
    obsid = 99903
    xray_src = Table(
        {
            "obsid": [obsid],
            "id": [1],
            "ra": [83.82],
            "dec": [-5.39],
            "net_counts": [500.0],
            "y_angle": [0.0],
            "z_angle": [0.0],
            "r_angle": [10.0],
            "snr": [15.0],
            "near_neighbor_dist": [60.0],
            "pileup": [0.0],
            "acis_streak": [1],  # on streak
            # no 'brightest' column at all — backward compat case
            "caldb_version": ["4.10.0"],
            "detect_method": ["celldetect"],
        }
    )

    matches = cross_match.compute_cross_matches(
        "astromon_21",
        astromon_obs=_minimal_obs(obsid),
        astromon_xray_src=xray_src,
        astromon_cat_src=_minimal_cat_src(obsid),
    )

    assert len(matches) == 0, "acis_streak source excluded when brightest column absent"


# Real records from the RFC "Format version of 2024.10.14" fixed-width catalog.
# _parse_rfc slices by byte offset, so this fixture pins those offsets: if astrogeo
# changes the column layout the parse silently yields wrong coordinates rather than
# raising, and only a test against real records catches it.
RFC_SAMPLE = """\
# VLBI SOURCE POSITION CATALOGUE  Format version of 2024.10.14
# The Radio Fundamental Catalogue, data release: rfc_2026b
#
RFC J0000-0022  2357-006  00 00 01.659837 -00 22 10.02335    1.26   2.68    0.230
RFC J1229+0203  3C273B    12 29 06.699755 +02 03 08.59833    0.09   0.17   -0.013
"""


def test_parse_rfc_column_offsets():
    """_parse_rfc reads name/ra/dec at the documented byte offsets."""
    rfc = cross_match._parse_rfc(RFC_SAMPLE)

    assert len(rfc) == 2, "comment lines and blank lines are skipped"
    assert list(rfc["name"]) == ["RFC J0000-0022", "RFC J1229+0203"]

    # 3C 273: 12h29m06.699755s, +02d03m08.59833s
    assert rfc["ra"][1] == pytest.approx(187.27791559, abs=1e-7)
    assert rfc["dec"][1] == pytest.approx(2.05238842, abs=1e-7)


def test_parse_rfc_negative_dec_with_zero_degrees():
    """A '-00' declination keeps its sign, which the degrees field alone cannot carry."""
    rfc = cross_match._parse_rfc(RFC_SAMPLE)

    # -00d22m10.02335s -> the sign lives only in the separate sign byte
    assert rfc["dec"][0] == pytest.approx(-0.36945093, abs=1e-7)
    assert rfc["dec"][0] < 0, "sign byte must be applied when degrees parse to zero"


def test_parse_rfc_empty_input_has_correct_columns():
    """An all-comment catalog yields an empty table that still has the expected columns."""
    rfc = cross_match._parse_rfc("# only a header\n\n")

    assert len(rfc) == 0
    assert set(rfc.colnames) == {"name", "ra", "dec"}


# ---- _local_catalog_near tests (shared by RFC, ICRF3, Quaia) ----


def _synthetic_local_catalog():
    """A tiny (name, ra, dec) table shaped like get_rfc()/get_icrf3()/get_quaia()."""
    return Table(
        {
            "name": ["near-1", "near-2", "far-away"],
            "ra": [10.0, 10.0005, 200.0],
            "dec": [-5.0, -5.0003, 60.0],
        }
    )


def test_local_catalog_near_filters_by_cone_and_formats_output():
    """_local_catalog_near keeps only in-cone sources and fills CROSS_MATCH_DTYPE."""
    pos = coords.SkyCoord(ra=[10.0], dec=[-5.0], unit="deg")

    result = cross_match._local_catalog_near(
        _synthetic_local_catalog(), "RFC", pos, radius=3 * u.arcsec
    )

    assert set(result.colnames) == set(cross_match.CROSS_MATCH_DTYPE.names)
    assert sorted(result["name"]) == ["near-1", "near-2"]
    assert np.all(result["catalog"] == "RFC")
    assert np.all(np.asarray(result["mag"].mask))


def test_local_catalog_near_empty_pos_returns_empty():
    """An empty position array short-circuits without touching the catalog table."""
    empty_pos = coords.SkyCoord([], [], unit="deg")

    result = cross_match._local_catalog_near(
        _synthetic_local_catalog(), "ICRF3", empty_pos, radius=3 * u.arcsec
    )

    assert len(result) == 0
    assert set(result.colnames) == set(cross_match.CROSS_MATCH_DTYPE.names)


def test_local_catalog_near_no_nearby_sources_returns_empty():
    """A cone with no catalog sources inside it returns an empty, correctly-typed table."""
    pos = coords.SkyCoord(ra=[123.4], dec=[-77.0], unit="deg")

    result = cross_match._local_catalog_near(
        _synthetic_local_catalog(), "Quaia", pos, radius=3 * u.arcsec
    )

    assert len(result) == 0
    assert set(result.colnames) == set(cross_match.CROSS_MATCH_DTYPE.names)


def test_local_catalog_near_scalar_pos():
    """A scalar SkyCoord (no len()) is handled the same as a length-1 array."""
    pos = coords.SkyCoord(ra=10.0, dec=-5.0, unit="deg")

    result = cross_match._local_catalog_near(
        _synthetic_local_catalog(), "RFC", pos, radius=3 * u.arcsec
    )

    assert sorted(result["name"]) == ["near-1", "near-2"]


def test_cache_is_stale_absent_file(tmp_path):
    """A cache_path that does not exist is always stale."""
    assert cross_match._cache_is_stale(tmp_path / "missing.txt", max_age_days=30)


def test_cache_is_stale_respects_max_age(tmp_path):
    """An existing file is stale only once it exceeds max_age_days."""
    cache_path = tmp_path / "cached.txt"
    cache_path.write_text("data")

    assert not cross_match._cache_is_stale(cache_path, max_age_days=30)

    old_time = time.time() - 31 * 86400
    os.utime(cache_path, (old_time, old_time))
    assert cross_match._cache_is_stale(cache_path, max_age_days=30)


def test_cache_is_stale_none_max_age_never_refreshes(tmp_path):
    """max_age_days=None means only absence makes the cache stale."""
    cache_path = tmp_path / "cached.txt"
    cache_path.write_text("data")

    old_time = time.time() - 10_000 * 86400
    os.utime(cache_path, (old_time, old_time))
    assert not cross_match._cache_is_stale(cache_path, max_age_days=None)


# ---- the anchor invariant ----


def test_the_anchor_always_names_a_celldetect_source():
    """Every anchor resolves to a celldetect source in the same observation.

    This is what makes the anchor meaningful, and it held with no exceptions across
    the live database: 99,558 of 99,558. It would have failed the moment a second
    detection method started being written without giving cat_src a per-method
    dimension, which is the state that let _remap_x_id exist.
    """
    cat_src = Table.read(DATA_DIR / "astromon_cat_src.ecsv")
    xray_src = Table.read(DATA_DIR / "astromon_xray_src.ecsv")

    if "detect_method" in xray_src.colnames:
        xray_src = xray_src[xray_src["detect_method"] == "celldetect"]
    celldetect = set(
        zip(
            np.asarray(xray_src["obsid"]).tolist(),
            np.asarray(xray_src["id"]).tolist(),
            strict=True,
        )
    )

    unresolved = [
        (int(obsid), int(anchor))
        for obsid, anchor in zip(
            cat_src["obsid"], cat_src["celldetect_x_id"], strict=True
        )
        if (int(obsid), int(anchor)) not in celldetect
    ]

    assert not unresolved, f"anchors naming no celldetect source: {unresolved[:5]}"


# ---- the celldetect anchor is not a match key ----


def test_matches_are_unchanged_when_the_anchor_is_scrambled():
    """cat_src's anchor must not decide which X-ray source a catalogue row matches.

    This is the contract, not an incidental property: the anchor is celldetect's
    numbering, and joining on it made it a per-method key in a table that has no
    per-method dimension. Matching is decided by position and the dr cut, so
    scrambling the anchor must change nothing. If someone reinstates the join, this
    fails.
    """
    obs = _minimal_obs(99900)
    xray = _minimal_xray_src(99900)
    good = _minimal_cat_src(99900)

    scrambled = good.copy()
    scrambled["celldetect_x_id"] = 4242  # names no source that exists

    from_good = cross_match.compute_cross_matches(
        "tycho2", astromon_obs=obs, astromon_xray_src=xray, astromon_cat_src=good
    )
    from_scrambled = cross_match.compute_cross_matches(
        "tycho2", astromon_obs=obs, astromon_xray_src=xray, astromon_cat_src=scrambled
    )

    assert len(from_good) == 1
    for column in ("obsid", "x_id", "c_id", "dr", "detect_method"):
        assert list(from_good[column]) == list(from_scrambled[column]), column


def test_a_method_missing_the_source_yields_no_match_for_that_method():
    """gaussian dropping a source means no gaussian match, not a substituted one.

    Gaussian ids are a subset of celldetect ids in 8,445 of 8,445 obsids, so when
    the anchor names an id gaussian does not have, the source was dropped by the
    fit. _remap_x_id substituted the nearest gaussian source instead, which
    manufactured matches against different objects -- obsid 332 got dr=0.09" against
    a source 2.4" from the real counterpart.
    """
    xray = _two_method_xray_src(99900)
    # celldetect keeps id 1 on the catalogue source; gaussian has only id 2, well
    # outside the 3" dr cut. dr is measured in the y/z plane, so that is what has
    # to move -- ra/dec alone would leave dr at zero for both.
    xray["id"] = [1, 2]
    xray["y_angle"] = [0.0, 36.0]

    matches = cross_match.compute_cross_matches(
        "tycho2",
        astromon_obs=_minimal_obs(99900),
        astromon_xray_src=xray,
        astromon_cat_src=_minimal_cat_src(99900),
    )

    assert list(matches["detect_method"]) == ["celldetect"]
    assert list(matches["x_id"]) == [1]


def test_empty_inputs_still_return_the_joined_dtype():
    """The empty path must not reference a column the schema no longer has."""
    empty_cat = _minimal_cat_src(99900)[:0]

    matches = cross_match.compute_cross_matches(
        "tycho2",
        astromon_obs=_minimal_obs(99900),
        astromon_xray_src=_minimal_xray_src(99900),
        astromon_cat_src=empty_cat,
    )

    assert len(matches) == 0
    # dy/dz/dr are computed after the join, so the empty path never carried them.
    # What matters is that it builds its dtype without naming a column the schema
    # no longer has, and that the anchor and the match key stay distinct names.
    for column in ("obsid", "x_id", "c_id", "celldetect_x_id"):
        assert column in matches.colnames, column


# ---- per-method selection grouping ----


def _two_method_xray_src(obsid=99900, ra=83.82, dec=-5.39):
    """The same source id detected by both methods, as the pipeline stores it."""
    both = vstack([_minimal_xray_src(obsid, ra=ra, dec=dec)] * 2)
    both["detect_method"] = ["celldetect", "gaussian_detect"]
    return both


def test_selection_keeps_one_row_per_method():
    """A selection without detect_method_filter must not collapse the two methods.

    The grouping that picks the best catalogue match per X-ray source keyed on
    (obsid, x_id) alone. Handed a table holding both methods it kept whichever
    happened to fit better and silently dropped the other -- so the surviving
    detect_method varied per source, and dr was biased low by being the minimum of
    two measurements. Every caller today passes one method at a time, which is the
    only reason this was not visible; rebuild_xcorr did not, and hit it.
    """
    matches = cross_match.compute_cross_matches(
        "tycho2",
        astromon_obs=_minimal_obs(99900),
        astromon_xray_src=_two_method_xray_src(99900),
        astromon_cat_src=_minimal_cat_src(99900),
    )

    assert sorted(matches["detect_method"]) == ["celldetect", "gaussian_detect"]


def test_selection_still_picks_one_match_per_source_within_a_method():
    """The grouping itself is unchanged: still one best catalogue match per source."""
    cat = vstack(
        [
            _minimal_cat_src(99900, name="near", ra=83.82028),
            _minimal_cat_src(99900, name="far", ra=83.8206),
        ]
    )
    cat["id"] = [0, 1]

    matches = cross_match.compute_cross_matches(
        "tycho2",
        astromon_obs=_minimal_obs(99900),
        astromon_xray_src=_minimal_xray_src(99900),
        astromon_cat_src=cat,
    )

    assert len(matches) == 1, "one X-ray source, one method -> one row"
    assert matches["name"][0] == "near"


# ---- rough_match deduplication ----


def test_rough_match_deduplicates_repeated_catalog_rows():
    """One row per (catalog source, X-ray source) pair, not per returned copy.

    get_vizier is called with every X-ray position at once and returns a catalogue
    source once per input position that matched it. rough_match then cartesian-joins
    that against the source list, so each returned copy became its own cat_src row:
    17,595 rows in the live database sat in exact-duplicate groups, every one of them
    from a VizieR catalogue.
    """
    # Two X-ray sources ~2" apart, and one catalogue source between them -- returned
    # twice, as VizieR does when both positions match it.
    sources = Table({"id": [1, 2], "ra": [10.0, 10.0005], "dec": [20.0, 20.0]})
    returned_twice = Table(
        {
            "catalog": ["Tycho2", "Tycho2"],
            "name": ["dup-source", "dup-source"],
            "ra": [10.00025, 10.00025],
            "dec": [20.0, 20.0],
            "mag": [11.0, 11.0],
        }
    )

    with patch.object(cross_match, "_get", return_value=returned_twice):
        result = cross_match.rough_match(
            sources, VIZIER_QUERY_TIME, radius=5 * u.arcsec, catalogs=["Tycho2"]
        )

    pairs = sorted(zip(result["name"], result["x_id"], strict=True))
    assert pairs == [("dup-source", 1), ("dup-source", 2)], (
        f"expected one row per X-ray source, got {pairs}"
    )


def test_rough_match_keeps_same_named_sources_at_different_positions():
    """Deduplication keys on position too, because SDSS names are not unique.

    613 groups in the live database share (obsid, catalog, name) while holding
    genuinely different positions, all of them SDSS. Keying on name alone would
    merge distinct objects.
    """
    sources = Table({"id": [1], "ra": [10.0], "dec": [20.0]})
    two_objects_one_name = Table(
        {
            "catalog": ["SDSS", "SDSS"],
            "name": ["shared-name", "shared-name"],
            "ra": [10.0002, 10.0004],
            "dec": [20.0, 20.0],
            "mag": [19.0, 20.0],
        }
    )

    with patch.object(cross_match, "_get", return_value=two_objects_one_name):
        result = cross_match.rough_match(
            sources, VIZIER_QUERY_TIME, radius=5 * u.arcsec, catalogs=["SDSS"]
        )

    assert len(result) == 2, "distinct positions must survive as distinct rows"
    assert sorted(round(float(r), 4) for r in result["ra"]) == [10.0002, 10.0004]


def _xray_sources_for_cat_src_ids(
    obsid=7001, ids=(3, 7), ras=(150.0, 150.01), decs=(2.0, 2.01)
):
    return Table(
        {
            "obsid": np.full(len(ids), obsid, dtype=np.int32),
            "id": np.array(ids, dtype=np.int32),
            "ra": np.array(ras, dtype=np.float64),
            "dec": np.array(decs, dtype=np.float64),
        }
    )


def _candidates_for_cat_src_ids(ras, decs, names=None, catalog="RFC"):
    names = names if names is not None else [f"c{i}" for i in range(len(ras))]
    return Table(
        {
            "ra": np.array(ras, dtype=np.float64),
            "dec": np.array(decs, dtype=np.float64),
            "name": np.array(names),
            "catalog": np.array([catalog] * len(ras)),
            "separation": np.zeros(len(ras), dtype=np.float32),
            "mag": np.full(len(ras), np.nan, dtype=np.float32),
        }
    )


def test_assign_cat_src_ids_stamps_obsid_ids_and_angles():
    candidates = _candidates_for_cat_src_ids([150.0, 150.01], [2.0, 2.01])

    cross_match.assign_cat_src_ids(candidates, 7001, _xray_sources_for_cat_src_ids(), QUAT)

    assert list(candidates["obsid"]) == [7001, 7001]
    assert list(candidates["id"]) == [0, 1]
    assert np.all(np.isfinite(candidates["y_angle"]))
    assert np.all(np.isfinite(candidates["z_angle"]))


def test_assign_cat_src_ids_anchors_to_the_nearest_source():
    candidates = _candidates_for_cat_src_ids([150.0, 150.01], [2.0, 2.01])

    cross_match.assign_cat_src_ids(
        candidates, 7001, _xray_sources_for_cat_src_ids(ids=(3, 7)), QUAT
    )

    assert list(candidates["celldetect_x_id"]) == [3, 7]


def test_assign_cat_src_ids_offsets_ids():
    candidates = _candidates_for_cat_src_ids([150.0, 150.01], [2.0, 2.01])

    cross_match.assign_cat_src_ids(
        candidates, 7001, _xray_sources_for_cat_src_ids(), QUAT, id_offset=42
    )

    assert list(candidates["id"]) == [42, 43]


def test_assign_cat_src_ids_keeps_an_existing_anchor():
    candidates = _candidates_for_cat_src_ids([150.0, 150.0], [2.0, 2.0])
    candidates["celldetect_x_id"] = np.array([3, 7], dtype=np.int32)

    cross_match.assign_cat_src_ids(
        candidates, 7001, _xray_sources_for_cat_src_ids(ids=(3, 7)), QUAT
    )

    assert list(candidates["celldetect_x_id"]) == [3, 7]


def test_assign_cat_src_ids_renames_existing_x_id_anchor():
    candidates = _candidates_for_cat_src_ids([150.0, 150.0], [2.0, 2.0])
    candidates["x_id"] = np.array([3, 7], dtype=np.int32)

    cross_match.assign_cat_src_ids(
        candidates, 7001, _xray_sources_for_cat_src_ids(ids=(3, 7)), QUAT
    )

    assert "x_id" not in candidates.colnames
    assert list(candidates["celldetect_x_id"]) == [3, 7]


def test_assign_cat_src_ids_fills_separation_from_the_matched_source():
    candidates = _candidates_for_cat_src_ids([150.0], [2.01])
    candidates.remove_column("separation")

    cross_match.assign_cat_src_ids(
        candidates, 7001, _xray_sources_for_cat_src_ids(ids=(3, 7)), QUAT
    )

    assert candidates["separation"][0] == pytest.approx(36.0, abs=0.1)


def test_assign_cat_src_ids_keeps_a_separation_that_is_already_there():
    candidates = _candidates_for_cat_src_ids([150.0], [2.01])
    candidates["separation"] = np.array([1.25], dtype=np.float32)

    cross_match.assign_cat_src_ids(candidates, 7001, _xray_sources_for_cat_src_ids(), QUAT)

    assert candidates["separation"][0] == pytest.approx(1.25)


def test_assign_cat_src_ids_handles_empty_candidates():
    candidates = _candidates_for_cat_src_ids([], [])

    cross_match.assign_cat_src_ids(candidates, 7001, _xray_sources_for_cat_src_ids(), QUAT)

    assert len(candidates) == 0
    for name in ("obsid", "id", "celldetect_x_id", "y_angle", "z_angle"):
        assert name in candidates.colnames


# ---- catalog cache location ----


def _resolved_cache_paths(env_overrides):
    """Import cross_match in a subprocess and report where its caches resolve.

    A subprocess because the paths are module-level constants built at import;
    setting the environment variable inside this process is too late.
    """
    import json
    import os
    import subprocess
    import sys

    script = (
        "import json\n"
        "from astromon import cross_match\n"
        "from astromon import observation, utils\n"
        "print(json.dumps({\n"
        "    'dir': str(utils.ASTROMON_DATA_DIR),\n"
        "    'archive': str(observation.ARCHIVE_DIR),\n"
        "    'rfc': str(cross_match._RFC_CACHE_PATH),\n"
        "    'quaia': str(cross_match._QUAIA_CACHE_PATH),\n"
        # Enumerated by introspection rather than from CATALOG_CACHE_PATHS, which
        # does not exist on every branch, and so that a cache constant added later
        # is covered without anyone remembering to list it.
        "    'all': {n: str(getattr(cross_match, n))\n"
        "            for n in dir(cross_match) if n.endswith('_CACHE_PATH')},\n"
        "}))\n"
    )
    repo_root = Path(__file__).resolve().parents[2]
    # Start from a copy with the override cleared, so a value in the ambient
    # environment cannot decide the result -- running the suite itself under
    # ASTROMON_DATA_DIR would otherwise make the default case untestable.
    child_env = {k: v for k, v in os.environ.items() if k != "ASTROMON_DATA_DIR"}
    child_env["PYTHONPATH"] = str(repo_root)
    child_env.update(env_overrides)
    completed = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        check=True,
        env=child_env,
    )
    return json.loads(completed.stdout)


def test_cache_dir_defaults_under_ska():
    """Without an override the caches stay where they have always been."""
    resolved = _resolved_cache_paths({"SKA": "/tmp/fake-ska"})

    assert resolved["dir"] == "/tmp/fake-ska/data/astromon"
    assert resolved["rfc"] == "/tmp/fake-ska/data/astromon/rfc_catalog.txt"


def test_astromon_data_dir_overrides_the_cache_location():
    """ASTROMON_DATA_DIR isolates the catalog caches the way ASTROMON_FILE does the DB.

    Without this the only way to move the caches was to move SKA itself, which
    drags mica, CALDB and everything else along with it -- so a dev run could
    isolate its database but never its catalogs.
    """
    resolved = _resolved_cache_paths(
        {"SKA": "/tmp/fake-ska", "ASTROMON_DATA_DIR": "/tmp/isolated-catalogs"}
    )

    assert resolved["dir"] == "/tmp/isolated-catalogs"
    assert resolved["rfc"] == "/tmp/isolated-catalogs/rfc_catalog.txt"
    assert resolved["quaia"] == "/tmp/isolated-catalogs/quaia_catalog.fits"


def test_archive_dir_defaults_under_the_shared_data_dir():
    """The observation archive lives under the same root as the catalogs."""
    resolved = _resolved_cache_paths({"SKA": "/tmp/fake-ska"})

    assert resolved["archive"] == "/tmp/fake-ska/data/astromon/xray_observations"


def test_astromon_data_dir_moves_the_archive_too():
    """ARCHIVE_DIR follows ASTROMON_DATA_DIR, so one variable isolates a run.

    Before this it was pinned to $SKA regardless, so a run could redirect its
    database and its catalogs and still archive observations into the shared
    tree.
    """
    resolved = _resolved_cache_paths(
        {"SKA": "/tmp/fake-ska", "ASTROMON_DATA_DIR": "/tmp/isolated-catalogs"}
    )

    assert resolved["archive"] == "/tmp/isolated-catalogs/xray_observations"


def test_astromon_data_dir_moves_every_catalog_cache():
    """No cache may be left behind pointing at the default location."""
    resolved = _resolved_cache_paths(
        {"SKA": "/tmp/fake-ska", "ASTROMON_DATA_DIR": "/tmp/isolated-catalogs"}
    )

    assert resolved["all"], "no cache constants found to check"
    stragglers = {
        name: path
        for name, path in resolved["all"].items()
        if not path.startswith("/tmp/isolated-catalogs/")
    }
    assert not stragglers, f"still under the default location: {stragglers}"


# ---- RFC release discovery and fail-soft caching tests ----

# The announcement page as astrogeo.org actually serves it: the current release
# is stated in a link to /sol/rfc/<release>, and the same release is named again
# in the citation text. Other releases appear only as image/PDF paths under /rfc/.
_RFC_LANDING_HTML = """
<html><body>
<B><A HREF="/sol/rfc/rfc_2026b"> rfc_2026b catalogue of compact radio sources</A></B>.
Please cite data release rfc_2026b and include DOI 10.25966/dhrk-zh08.
<A HREF="/rfc/rfc_2026a_map.pdf"><IMG SRC="/rfc/rfc_2026a_map.png"></A>
</body></html>
"""


def test_discover_latest_rfc_release_reads_the_announced_release():
    """Discovery reports the release astrogeo announces, not the newest directory.

    The quarterly directory for a release appears under /sol/rfc/ before its data
    files do -- rfc_2026c had a landing page and no catalogue for weeks -- so
    taking the lexicographically-largest directory name adopts a release that has
    not been published and every download from it 404s.
    """
    fake_response = Mock(text=_RFC_LANDING_HTML)
    fake_response.raise_for_status = Mock()
    with patch.object(cross_match.requests, "get", return_value=fake_response):
        release = cross_match._discover_latest_rfc_release()
    assert release == "rfc_2026b"


def test_discover_latest_rfc_release_ignores_releases_named_only_as_assets():
    """A map PDF for another release must not be mistaken for the announcement."""
    html = (
        '<A HREF="/sol/rfc/rfc_2025d">rfc_2025d</A>'
        '<A HREF="/rfc/rfc_2026b_map.pdf">map</A>'
    )
    fake_response = Mock(text=html)
    fake_response.raise_for_status = Mock()
    with patch.object(cross_match.requests, "get", return_value=fake_response):
        assert cross_match._discover_latest_rfc_release() == "rfc_2025d"


def test_discover_latest_rfc_release_no_releases_found_raises():
    """A page announcing no release raises rather than returning junk."""
    fake_response = Mock(text="<html><body>nothing here</body></html>")
    fake_response.raise_for_status = Mock()
    with patch.object(cross_match.requests, "get", return_value=fake_response):
        with pytest.raises(RuntimeError, match="No RFC release"):
            cross_match._discover_latest_rfc_release()


def _write_fake_rfc_file(dest, *_args, **_kwargs):
    """Write a minimal valid RFC catalog to *dest* directly (not a urlretrieve stand-in)."""
    Path(dest).write_text(
        "RFC J1229+0203  3C273B    12 29 06.699755 +02 03 08.59833    0.09   0.17   -0.013\n"
    )


def _urlretrieve_writes_fake_rfc(_url, dest, *_args, **_kwargs):
    """Stand-in for urllib.request.urlretrieve(url, dest, ...): write to *dest*."""
    _write_fake_rfc_file(dest)


def test_get_rfc_no_cache_downloads_and_writes_release_marker(tmp_path):
    """First-ever call with no cache: discovers the release, downloads, records it."""
    cache = tmp_path / "rfc_catalog.txt"
    with (
        patch.object(
            cross_match, "_discover_latest_rfc_release", return_value="rfc_2026b"
        ),
        patch(
            "astromon.cross_match.urllib.request.urlretrieve",
            side_effect=_urlretrieve_writes_fake_rfc,
        ) as mock_urlretrieve,
    ):
        rfc = cross_match.get_rfc(cache_path=cache)

    assert len(rfc) == 1
    assert mock_urlretrieve.call_count == 1
    assert mock_urlretrieve.call_args[0][0] == (
        "http://astrogeo.org/sol/rfc/rfc_2026b/rfc_2026b_cat.txt"
    )
    marker = cache.with_name(cache.name + ".release")
    assert marker.read_text() == "rfc_2026b"


def test_get_rfc_fresh_cache_skips_network_entirely(tmp_path):
    """A cache newer than max_age_days never touches the network, not even discovery."""
    cache = tmp_path / "rfc_catalog.txt"
    _write_fake_rfc_file(cache)
    cache.with_name(cache.name + ".release").write_text("rfc_2026b")

    with (
        patch.object(cross_match, "_discover_latest_rfc_release") as mock_discover,
        patch("astromon.cross_match.urllib.request.urlretrieve") as mock_urlretrieve,
    ):
        rfc = cross_match.get_rfc(cache_path=cache, max_age_days=1)

    assert len(rfc) == 1
    mock_discover.assert_not_called()
    mock_urlretrieve.assert_not_called()


def test_get_rfc_stale_but_same_release_skips_download(tmp_path):
    """When the check interval elapses but the release hasn't changed, no download happens."""
    cache = tmp_path / "rfc_catalog.txt"
    _write_fake_rfc_file(cache)
    cache.with_name(cache.name + ".release").write_text("rfc_2026b")
    old_time = time.time() - 2 * 86400
    os.utime(cache, (old_time, old_time))

    with (
        patch.object(
            cross_match, "_discover_latest_rfc_release", return_value="rfc_2026b"
        ),
        patch("astromon.cross_match.urllib.request.urlretrieve") as mock_urlretrieve,
    ):
        rfc = cross_match.get_rfc(cache_path=cache, max_age_days=1)

    assert len(rfc) == 1
    mock_urlretrieve.assert_not_called()
    # the check-interval clock is reset even though nothing downloaded
    assert cache.stat().st_mtime > old_time


def test_get_rfc_new_release_downloads(tmp_path):
    """A genuinely newer release name triggers a real download and marker update."""
    cache = tmp_path / "rfc_catalog.txt"
    _write_fake_rfc_file(cache)
    cache.with_name(cache.name + ".release").write_text("rfc_2026a")
    old_time = time.time() - 2 * 86400
    os.utime(cache, (old_time, old_time))

    with (
        patch.object(
            cross_match, "_discover_latest_rfc_release", return_value="rfc_2026b"
        ),
        patch(
            "astromon.cross_match.urllib.request.urlretrieve",
            side_effect=_urlretrieve_writes_fake_rfc,
        ) as mock_urlretrieve,
    ):
        cross_match.get_rfc(cache_path=cache, max_age_days=1)

    assert mock_urlretrieve.call_count == 1
    marker = cache.with_name(cache.name + ".release")
    assert marker.read_text() == "rfc_2026b"


def test_get_rfc_download_failure_falls_back_to_stale_cache(tmp_path):
    """A download failure for a newly-discovered release (e.g. astrogeo publishing the
    release directory before the data file itself is live) falls back to the existing
    cache with a warning, rather than crashing the whole cross-match run."""
    cache = tmp_path / "rfc_catalog.txt"
    _write_fake_rfc_file(cache)
    cache.with_name(cache.name + ".release").write_text("rfc_2026b")
    old_time = time.time() - 2 * 86400
    os.utime(cache, (old_time, old_time))
    size_before = cache.stat().st_size

    with (
        patch.object(
            cross_match, "_discover_latest_rfc_release", return_value="rfc_2026c"
        ),
        patch(
            "astromon.cross_match.urllib.request.urlretrieve",
            side_effect=OSError("HTTP Error 404: Not Found"),
        ),
        patch.object(cross_match.logger, "warning") as mock_warning,
    ):
        rfc = cross_match.get_rfc(cache_path=cache, max_age_days=1)

    assert len(rfc) == 1  # stale data still served, not an empty/broken table
    assert cache.stat().st_size == size_before  # cache file untouched
    assert cache.with_name(cache.name + ".release").read_text() == "rfc_2026b"
    assert "RFC download failed" in mock_warning.call_args[0][0]


def test_get_rfc_discovery_failure_no_cache_raises(tmp_path):
    """With no cache to fall back to, a discovery failure has nothing to serve and raises."""
    cache = tmp_path / "rfc_catalog.txt"
    with patch.object(
        cross_match, "_discover_latest_rfc_release", side_effect=OSError("network down")
    ):
        with pytest.raises(OSError, match="network down"):
            cross_match.get_rfc(cache_path=cache)


def test_get_rfc_discovery_failure_with_cache_falls_back(tmp_path):
    """A transient failure checking for a new release falls back to the existing cache."""
    cache = tmp_path / "rfc_catalog.txt"
    _write_fake_rfc_file(cache)
    cache.with_name(cache.name + ".release").write_text("rfc_2026b")
    old_time = time.time() - 2 * 86400
    os.utime(cache, (old_time, old_time))

    with (
        patch.object(
            cross_match,
            "_discover_latest_rfc_release",
            side_effect=OSError("network down"),
        ),
        patch.object(cross_match.logger, "warning") as mock_warning,
    ):
        rfc = cross_match.get_rfc(cache_path=cache, max_age_days=1)

    assert len(rfc) == 1
    assert "RFC release check failed" in mock_warning.call_args[0][0]


def _local_catalog(name="SRC", ra=150.0, dec=20.0, mag=None):
    """Minimal locally-cached catalog table as consumed by _local_catalog_near."""
    cols = {
        "name": [name],
        "ra": np.array([ra], dtype=float),
        "dec": np.array([dec], dtype=float),
    }
    if mag is not None:
        cols["mag"] = np.array([mag], dtype=np.float32)
    return Table(cols)


def test_local_catalog_near_wraps_ra_across_zero_meridian():
    """A counterpart just across RA=0 is kept, not discarded by unwrapped RA arithmetic.

    The bounding-cone prefilter works in flat (dra, ddec) offsets. Without wrapping,
    a source at RA=359.99999 looks ~354 deg away from a position at RA=0.00001 rather
    than 1e-5 deg, so every counterpart across the meridian is silently dropped --
    losing all RFC/ICRF3/Quaia/GaiaAGN/GaiaQSO matches for fields straddling RA=0.
    """
    pos = coords.SkyCoord(ra=[0.00001], dec=[10.0], unit="deg")
    catalog = _local_catalog(name="ACROSS_MERIDIAN", ra=359.99999, dec=10.0)

    true_sep = pos.separation(
        coords.SkyCoord(catalog["ra"], catalog["dec"], unit="deg")
    ).arcsec
    assert true_sep[0] < 1.0, "test fixture must be a genuine sub-arcsec counterpart"

    result = cross_match._local_catalog_near(catalog, "TEST", pos, 3 * u.arcsec)
    assert len(result) == 1
    assert result["name"][0] == "ACROSS_MERIDIAN"


def test_local_catalog_near_wraps_ra_for_field_straddling_meridian():
    """Both halves of a field that straddles RA=0 keep their counterparts."""
    pos = coords.SkyCoord(ra=[359.98, 0.02], dec=[-30.0, -30.0], unit="deg")
    catalog = Table(
        {
            "name": ["west", "east"],
            "ra": np.array([359.980001, 0.020001]),
            "dec": np.array([-30.0, -30.0]),
        }
    )

    result = cross_match._local_catalog_near(catalog, "TEST", pos, 3 * u.arcsec)
    assert sorted(result["name"]) == ["east", "west"]


def test_local_catalog_near_applies_exact_radius_cut():
    """Sources beyond `radius` of every input position are excluded.

    The bounding cone is sized to the whole field, so on its own it admits sources
    arcminutes away from any X-ray source. get_gaia_agn/get_gaia_qso_candidates/
    get_quaia_candidates feed straight into astromon_cat_src with no later cut of
    their own, so the radius has to be enforced here.
    """
    pos = coords.SkyCoord(ra=[150.0, 150.1], dec=[20.0, 20.05], unit="deg")
    # Midway between the two X-ray sources: ~3.2 arcmin from the nearer one.
    catalog = _local_catalog(name="FAR", ra=150.05, dec=20.025)

    nearest = pos.separation(
        coords.SkyCoord(catalog["ra"], catalog["dec"], unit="deg")
    ).arcsec.min()
    assert nearest > 3.0, "test fixture must lie outside the search radius"

    result = cross_match._local_catalog_near(catalog, "TEST", pos, 3 * u.arcsec)
    assert len(result) == 0


def test_local_catalog_near_keeps_source_near_any_position():
    """The radius cut is per-position, not relative to the field centroid."""
    pos = coords.SkyCoord(ra=[150.0, 150.1], dec=[20.0, 20.05], unit="deg")
    catalog = Table(
        {
            "name": ["near_second", "far"],
            # 1 arcsec from the second position, and far from both.
            "ra": np.array([150.1 + 1.0 / 3600.0 / np.cos(np.radians(20.05)), 150.05]),
            "dec": np.array([20.05, 20.025]),
        }
    )

    result = cross_match._local_catalog_near(catalog, "TEST", pos, 3 * u.arcsec)
    assert list(result["name"]) == ["near_second"]


def test_local_catalog_near_scalar_pos_applies_radius_cut():
    """A scalar SkyCoord is handled and still gets the exact radius cut."""
    pos = coords.SkyCoord(ra=150.0, dec=20.0, unit="deg")
    catalog = Table(
        {
            "name": ["close", "far"],
            "ra": np.array([150.0 + 1.0 / 3600.0 / np.cos(np.radians(20.0)), 150.02]),
            "dec": np.array([20.0, 20.0]),
        }
    )

    result = cross_match._local_catalog_near(catalog, "TEST", pos, 3 * u.arcsec)
    assert list(result["name"]) == ["close"]


def test_get_gaia_agn_excludes_source_outside_radius():
    """The public getter does not emit candidates beyond `radius` of any position."""
    pos = coords.SkyCoord([187.2779, 187.3779], [2.0524, 2.0524], unit="deg")
    fake = _fake_cached_gaia_catalog(
        "GaiaAGN", source_id=3726862715688776192, ra=187.3279, dec=2.0524
    )

    with patch.object(cross_match, "get_gaia_agn_catalog", return_value=fake):
        result = cross_match.get_gaia_agn(pos, radius=3 * u.arcsec)

    assert len(result) == 0


# --- get_bad_source_mask ---------------------------------------------------


def _matches_for_mask(targets, names=None, obsids=None, x_ids=None):
    """Minimal matches table with the columns get_bad_source_mask inspects."""
    n = len(targets)
    return Table(
        {
            "target": np.array(targets, dtype="U40"),
            "name": np.array(names if names is not None else [""] * n, dtype="U40"),
            "obsid": np.array(
                obsids if obsids is not None else [0] * n, dtype=np.int32
            ),
            "x_id": np.array(x_ids if x_ids is not None else [0] * n, dtype=np.int32),
        }
    )


def test_get_bad_source_mask_empty_matches_returns_empty_mask():
    """An empty matches table yields an empty bool mask rather than raising.

    The list comprehensions produce a float64 empty array, so `exclude |= ...`
    raised TypeError. compute_cross_matches(exclude_bad_targets=True) hits this
    for any obsid whose matches are all removed by the dr/r_angle cuts.
    """
    empty = _matches_for_mask([])
    mask = cross_match.get_bad_source_mask(empty)
    assert mask.dtype == bool
    assert len(mask) == 0


def test_get_bad_target_mask_empty_matches_returns_empty_mask():
    """The deprecated alias is empty-safe too (it is what filter_matches calls)."""
    mask = cross_match.get_bad_target_mask(_matches_for_mask([]))
    assert mask.dtype == bool
    assert len(mask) == 0


def test_get_bad_source_mask_excludes_target_prefix():
    """Target-prefix rows match space-stripped, case-insensitive, on prefix."""
    matches = _matches_for_mask(["RW Aur A", "rwaur", "NGC 1234"])
    mask = cross_match.get_bad_source_mask(matches)
    assert list(mask) == [True, True, False]


def test_get_bad_source_mask_excludes_cat_source_name():
    """cat_source rows exclude by exact catalog source name, regardless of target."""
    matches = _matches_for_mask(
        ["NGC 1234", "NGC 1234"], names=["MCG 7-41-003", "something else"]
    )
    mask = cross_match.get_bad_source_mask(matches)
    assert list(mask) == [True, False]


def test_get_bad_source_mask_excludes_obsid_src(tmp_path):
    """obsid_src rows exclude one (obsid, x_id) pair only."""
    bad_file = tmp_path / "bad_sources.ecsv"
    bad_file.write_text(
        "# %ECSV 1.0\n"
        "# ---\n"
        "# datatype:\n"
        "# - {name: excl_type, datatype: string}\n"
        "# - {name: target_prefix, datatype: string}\n"
        "# - {name: obsid, datatype: int32}\n"
        "# - {name: x_id, datatype: int32}\n"
        "# - {name: cat_name, datatype: string}\n"
        "# - {name: reason, datatype: string}\n"
        "# schema: astropy-2.0\n"
        "excl_type target_prefix obsid x_id cat_name reason\n"
        'obsid_src "" 7001 3 "" testing\n'
    )
    matches = _matches_for_mask(
        ["A", "B", "C"], obsids=[7001, 7001, 7002], x_ids=[3, 4, 3]
    )
    mask = cross_match.get_bad_source_mask(matches, bad_sources_file=bad_file)
    assert list(mask) == [True, False, False]


def test_get_gaia_var_stars_accepts_a_scalar_position():
    """The docstring says a scalar SkyCoord works; the batching loop said otherwise.

    ``len()`` on a scalar SkyCoord raises TypeError, so ``range(0, len(pos), ...)``
    could not run at all for a single-position caller -- while the surrounding code
    had already normalised the coordinates with ``np.atleast_1d``, which is what
    made the mismatch easy to miss.
    """
    scalar = coords.SkyCoord(150.0, 2.0, unit="deg", obstime=CxoTime("2020:001"))
    assert scalar.isscalar

    with patch.object(
        cross_match,
        "_execute_gaia_tap_query",
        return_value=_fake_var_tap_result(4242, 150.0, 2.0),
    ):
        result = cross_match.get_gaia_var_stars(scalar, radius=3 * u.arcsec)

    assert len(result) == 1
    assert result["catalog"][0] == "GaiaVarStar"
