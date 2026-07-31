from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
import requests
from astropy import coordinates as coords
from astropy import units as u
from astropy.table import Table
from cxotime import CxoTime
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
    "ICRS": {
        "search_ra": 187.2779155448750,
        "search_dec": 2.0523884103056,
        "radius": 2 * u.arcsec,
        "name": "J122906.6+020308",
        "ra": 187.2779155448750,
        "dec": 2.0523884103056,
        "mag": None,  # ICRS has no magnitude column, mag is expected to stay masked
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
# catalog (as opposed to a wrong column name) would surface in production for ICRS/Tycho2/
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
    assert result["catalog"][0] == cross_match.VIZIER_CATALOGS[catalog]["catalog"]
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


def test_get_gaia_agn_empty_pos():
    """Empty position array returns an empty table with the correct schema."""
    empty_pos = coords.SkyCoord([], [], unit="deg")
    result = cross_match.get_gaia_agn(empty_pos)
    assert len(result) == 0
    assert set(cross_match.CROSS_MATCH_DTYPE.names).issubset(set(result.colnames))


def test_get_gaia_agn_mocked_tap():
    """get_gaia_agn structures output correctly from a mocked TAP response."""
    pos = coords.SkyCoord([187.2779], [2.0524], unit="deg")
    fake = _fake_tap_result(source_id=3726862715688776192, ra=187.2779, dec=2.0524, mag=12.9)

    with patch.object(cross_match, "_execute_gaia_tap_query", return_value=fake):
        result = cross_match.get_gaia_agn(pos, radius=3 * u.arcsec)

    assert len(result) == 1
    assert result["catalog"][0] == "GaiaAGN"
    assert result["name"][0] == "GaiaAGN-3726862715688776192"
    assert np.isclose(result["ra"][0], 187.2779)
    assert np.isclose(result["dec"][0], 2.0524)
    assert np.isclose(result["mag"][0], 12.9, atol=0.01)


def test_get_gaia_agn_deduplicates_across_batches():
    """The same Gaia source returned by two batches appears only once in the output."""
    pos = coords.SkyCoord([187.2779, 187.2780], [2.0524, 2.0525], unit="deg")
    # Both batch calls return the same source_id
    fake = _fake_tap_result(source_id=3726862715688776192, ra=187.2779, dec=2.0524)

    with patch.object(cross_match, "_execute_gaia_tap_query", return_value=fake):
        result = cross_match.get_gaia_agn(pos, radius=3 * u.arcsec)

    assert len(result) == 1


def test_get_gaia_agn_tap_failure_returns_empty():
    """A TAP query failure returns an empty table without raising (graceful degradation)."""
    pos = coords.SkyCoord([187.2779], [2.0524], unit="deg")

    with patch.object(
        cross_match, "_execute_gaia_tap_query", side_effect=RuntimeError("TAP unavailable")
    ):
        result = cross_match.get_gaia_agn(pos, radius=3 * u.arcsec)

    assert len(result) == 0
    assert set(cross_match.CROSS_MATCH_DTYPE.names).issubset(set(result.colnames))


def test_compute_cross_matches_gaia_agn():
    """compute_cross_matches('gaia_agn') finds a GaiaAGN match from a synthetic cat_src."""
    obs = Table(
        {
            "obsid": [99999],
            "detector": ["ACIS-S"],
            "target": ["test"],
            "grating": ["NONE"],
            "sim_z": [-190.143],
            "date_obs": ["2020-01-01T00:00:00"],
            "tstart": [0.0],
            "ascdsver": ["10.0"],
            "ra": [187.2779],
            "dec": [2.0524],
            "roll": [0.0],
            "category_id": [50],  # Active Galaxies and Quasars
            "version": ["10.0"],
        }
    )

    xray_src = Table(
        {
            "obsid": [99999],
            "id": [1],
            "ra": [187.2779],
            "dec": [2.0524],
            "net_counts": [500.0],
            "y_angle": [0.0],
            "z_angle": [0.0],
            "r_angle": [10.0],
            "snr": [10.0],
            "near_neighbor_dist": [60.0],
            "pileup": [0.0],
            "acis_streak": [False],
            "caldb_version": ["4.10.0"],
            "detect_method": ["celldetect"],
        }
    )

    # GaiaAGN catalog source ~1 arcsec away from xray source 1.
    # x_id=1 links this candidate to xray_src id=1, as simple_cross_match joins on (obsid, x_id).
    cat_src = Table(
        {
            "obsid": [99999],
            "id": [0],
            "x_id": [1],
            "catalog": ["GaiaAGN"],
            "name": ["GaiaAGN-3726862715688776192"],
            "ra": [187.2782],
            "dec": [2.0524],
            "mag": [12.9],
            "y_angle": [0.0],
            "z_angle": [0.0],
            "separation": [1.0],
        }
    )

    matches = cross_match.compute_cross_matches(
        "gaia_agn",
        astromon_obs=obs,
        astromon_xray_src=xray_src,
        astromon_cat_src=cat_src,
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


@pytest.mark.skipif(not HAS_INTERNET, reason="Requires network access")
def test_get_gaia_agn_live():
    """get_gaia_agn finds a known Gaia DR3 AGN cross-id source via live TAP query."""
    pos = coords.SkyCoord([_AGNCROSS_RA], [_AGNCROSS_DEC], unit="deg")
    result = cross_match.get_gaia_agn(pos, radius=3 * u.arcsec)

    assert len(result) >= 1
    assert result["catalog"][0] == "GaiaAGN"
    assert result["name"][0] == f"GaiaAGN-{_AGNCROSS_SOURCE_ID}"
    assert np.isclose(result["ra"][0], _AGNCROSS_RA, atol=1e-3)
    assert np.isclose(result["dec"][0], _AGNCROSS_DEC, atol=1e-3)
