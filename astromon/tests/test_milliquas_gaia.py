"""Tests for MilliquasGaia served from a precomputed all-sky crossmatch.

Queried live, this catalog cost a VizieR cone search plus a Gaia TAP position
upgrade per obsid -- ~17 s each, and 23.5 of the ~78 hours a full re-query of all
8,451 obsids would take. The upgrade does not depend on the observation at all,
so it can be done once for the whole catalogue: CDS X-Match crossmatches VII/294
against Gaia DR3 (I/355) server-side in one request, 664,428 pairs all-sky.

That reduction was checked against the live path before being adopted, not
derived from it: all 5,504 MilliquasGaia rows the live path had already written
came back with positions agreeing to below 1 mas and Rmag to 0.005 mag.
test_reproduces_what_the_live_path_produced keeps a sample of that pinned.
"""

from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
from astropy import coordinates as coords
from astropy import units as u
from astropy.table import Table

from astromon import cross_match

DATA_DIR = Path(__file__).parent / "data"
RADIUS = 3 * u.arcsec


def _xmatch_votable(rows):
    """An X-Match response VOTable: (angDist, Name, Type, Rmag, RAdeg, DEdeg)."""
    body = "".join(
        "<TR>"
        f"<TD>{ang}</TD><TD>{name}</TD><TD>{type_}</TD><TD>{rmag}</TD>"
        f"<TD>{ra}</TD><TD>{dec}</TD>"
        "</TR>"
        for ang, name, type_, rmag, ra, dec in rows
    )
    return (
        b'<?xml version="1.0" encoding="UTF-8"?>\n'
        b'<VOTABLE version="1.4" xmlns="http://www.ivoa.net/xml/VOTable/v1.3">\n'
        b'<INFO name="QUERY_STATUS" value="OK"/>\n'
        b"<RESOURCE><TABLE>"
        b'<FIELD name="angDist" datatype="double" unit="arcsec"/>'
        b'<FIELD name="Name" datatype="char" arraysize="*"/>'
        b'<FIELD name="Type" datatype="char" arraysize="*"/>'
        b'<FIELD name="Rmag" datatype="float"/>'
        b'<FIELD name="RAdeg" datatype="double"/>'
        b'<FIELD name="DEdeg" datatype="double"/>'
        b"<DATA><TABLEDATA>" + body.encode() + b"</TABLEDATA></DATA>"
        b"</TABLE></RESOURCE></VOTABLE>\n"
    )


class _StubResponse:
    def __init__(self, content):
        self.content = content
        self.status_code = 200

    def raise_for_status(self):
        return None


def _posting(content, calls=None):
    def post(url, data=None, timeout=None, **kwargs):
        if calls is not None:
            calls.append(data)
        return _StubResponse(content)

    return post


# ─── building the mapping ────────────────────────────────────────────────────


def test_catalog_keeps_only_spectroscopically_confirmed_types():
    """K-type photometric candidates and X/R associations are not counterparts."""
    votable = _xmatch_votable(
        [
            (0.10, "QSO A", "Q", 18.0, 10.0, 20.0),
            (0.10, "AGN B", "A", 19.0, 11.0, 21.0),
            (0.10, "BLLAC C", "B", 17.0, 12.0, 22.0),
            (0.10, "NLAGN D", "N", 20.0, 13.0, 23.0),
            (0.10, "PHOT E", "K", 21.0, 14.0, 24.0),
            (0.10, "XRAY F", "X", 21.0, 15.0, 25.0),
        ]
    )
    with patch.object(cross_match.requests, "post", _posting(votable)):
        catalog = cross_match._download_milliquas_gaia_xmatch()
    reduced = cross_match._reduce_milliquas_gaia_xmatch(catalog)

    assert sorted(reduced["name"]) == ["AGN B", "BLLAC C", "NLAGN D", "QSO A"]


def test_catalog_keeps_the_nearest_gaia_counterpart_per_source():
    """X-Match returns every pair within the radius; only the nearest is the match."""
    votable = _xmatch_votable(
        [
            (0.90, "QSO A", "Q", 18.0, 10.0009, 20.0),
            (0.05, "QSO A", "Q", 18.0, 10.0001, 20.0),
            (0.40, "QSO A", "Q", 18.0, 10.0004, 20.0),
        ]
    )
    with patch.object(cross_match.requests, "post", _posting(votable)):
        catalog = cross_match._download_milliquas_gaia_xmatch()
    reduced = cross_match._reduce_milliquas_gaia_xmatch(catalog)

    assert len(reduced) == 1
    assert reduced["ra"][0] == pytest.approx(10.0001)


def test_catalog_carries_the_gaia_position_not_the_milliquas_one():
    """The whole point of the upgrade is that the position comes from Gaia."""
    votable = _xmatch_votable([(0.10, "QSO A", "Q", 18.5, 10.00025, 20.00025)])
    with patch.object(cross_match.requests, "post", _posting(votable)):
        catalog = cross_match._download_milliquas_gaia_xmatch()
    reduced = cross_match._reduce_milliquas_gaia_xmatch(catalog)

    assert set(reduced.colnames) >= {"name", "ra", "dec", "mag"}
    assert reduced["ra"][0] == pytest.approx(10.00025)
    assert reduced["mag"][0] == pytest.approx(18.5)


def test_catalog_refuses_a_response_that_reports_an_error():
    """An X-Match error must not become an empty catalog written to the cache."""
    error = (
        b'<?xml version="1.0" encoding="UTF-8"?>\n'
        b'<VOTABLE version="1.4" xmlns="http://www.ivoa.net/xml/VOTable/v1.3">\n'
        b'<INFO name="QUERY_STATUS" value="ERROR">service overloaded</INFO>\n'
        b"</VOTABLE>\n"
    )
    with (
        patch.object(cross_match.requests, "post", _posting(error)),
        pytest.raises(cross_match.VizierServerError, match="overloaded"),
    ):
        cross_match._download_milliquas_gaia_xmatch()


def test_catalog_asks_for_only_the_columns_it_needs():
    """All 157 columns of both catalogues is 130 MB more than this needs."""
    calls = []
    votable = _xmatch_votable([(0.10, "QSO A", "Q", 18.0, 10.0, 20.0)])
    with patch.object(cross_match.requests, "post", _posting(votable, calls)):
        cross_match._download_milliquas_gaia_xmatch()

    (payload,) = calls
    assert payload["cols1"] == "Name,RAJ2000,DEJ2000,Rmag,Type"
    assert payload["cols2"] == "Source,RAdeg,DEdeg"
    assert payload["area"] == "allsky"
    assert float(payload["distMaxArcsec"]) == pytest.approx(1.5)


def test_catalog_is_cached_and_not_rebuilt_when_the_version_matches(tmp_path):
    cache = tmp_path / "milliquas_gaia_catalog.fits"
    Table(
        {
            "name": np.array(["QSO A"]),
            "ra": np.array([10.0]),
            "dec": np.array([20.0]),
            "mag": np.array([18.0], dtype=np.float32),
        }
    ).write(cache, format="fits")
    cross_match._record_catalog_version(cache, cross_match._MILLIQUAS_GAIA_VERSION)

    def explode(*args, **kwargs):
        raise AssertionError("a cached mapping must not be rebuilt")

    with patch.object(cross_match.requests, "post", explode):
        catalog = cross_match.get_milliquas_gaia_catalog(cache_path=cache)

    assert list(catalog["name"]) == ["QSO A"]


def test_catalog_is_rebuilt_when_the_version_pin_moves(tmp_path):
    """A new Milliquas release or Gaia data release invalidates the mapping."""
    cache = tmp_path / "milliquas_gaia_catalog.fits"
    Table(
        {
            "name": np.array(["STALE"]),
            "ra": np.array([1.0]),
            "dec": np.array([2.0]),
            "mag": np.array([9.0], dtype=np.float32),
        }
    ).write(cache, format="fits")
    cross_match._record_catalog_version(cache, "VII/290+I/350")

    votable = _xmatch_votable([(0.10, "QSO A", "Q", 18.0, 10.0, 20.0)])
    with patch.object(cross_match.requests, "post", _posting(votable)):
        catalog = cross_match.get_milliquas_gaia_catalog(cache_path=cache)

    assert list(catalog["name"]) == ["QSO A"]
    recorded = cross_match._release_marker(cache).read_text().strip()
    assert recorded == cross_match._MILLIQUAS_GAIA_VERSION


def test_catalog_is_registered_for_version_reporting():
    assert "MilliquasGaia" in cross_match.CATALOG_CACHE_PATHS


# ─── serving it per obsid ────────────────────────────────────────────────────


def _mini_catalog():
    return Table(
        {
            "name": np.array(["NEAR ONE", "NEAR TWO", "FAR AWAY"]),
            "ra": np.array([150.0, 150.0005, 220.0]),
            "dec": np.array([2.0, 2.0, -15.0]),
            "mag": np.array([18.0, 19.0, 20.0], dtype=np.float32),
        }
    )


def test_serves_candidates_from_the_precomputed_mapping():
    pos = coords.SkyCoord([150.0], [2.0], unit="deg")
    with patch.object(cross_match, "get_milliquas_gaia_catalog", _mini_catalog):
        result = cross_match.get_milliquas_gaia(pos, radius=RADIUS)

    assert sorted(result["name"]) == ["NEAR ONE", "NEAR TWO"]
    assert np.all(result["catalog"] == "MilliquasGaia")


def test_excludes_sources_beyond_the_radius():
    pos = coords.SkyCoord([220.0], [-15.0], unit="deg")
    with patch.object(cross_match, "get_milliquas_gaia_catalog", _mini_catalog):
        result = cross_match.get_milliquas_gaia(pos, radius=RADIUS)

    assert list(result["name"]) == ["FAR AWAY"]


def test_returns_an_empty_table_for_a_field_with_nothing_in_it():
    pos = coords.SkyCoord([300.0], [60.0], unit="deg")
    with patch.object(cross_match, "get_milliquas_gaia_catalog", _mini_catalog):
        result = cross_match.get_milliquas_gaia(pos, radius=RADIUS)

    assert len(result) == 0


def test_reproduces_what_the_live_path_produced():
    """Pinned against real rows the live VizieR + Gaia TAP path wrote.

    The fixture pairs 24 stored MilliquasGaia rows with the same sources as the
    all-sky X-Match returns them. All 5,504 stored rows agreed when this was
    adopted; these keep a spread across the sky honest.
    """
    reference = Table.read(DATA_DIR / "milliquas_gaia_reference.ecsv")
    catalog = Table(
        {
            "name": np.asarray(reference["name"]),
            "ra": np.asarray(reference["ra"]),
            "dec": np.asarray(reference["dec"]),
            "mag": np.asarray(reference["mag"]),
        }
    )

    with patch.object(cross_match, "get_milliquas_gaia_catalog", lambda: catalog):
        for row in reference:
            pos = coords.SkyCoord([row["stored_ra"]], [row["stored_dec"]], unit="deg")
            result = cross_match.get_milliquas_gaia(pos, radius=RADIUS)

            assert row["name"] in list(result["name"]), row["name"]
            found = result[result["name"] == row["name"]][0]
            separation = coords.SkyCoord(
                found["ra"], found["dec"], unit="deg"
            ).separation(
                coords.SkyCoord(row["stored_ra"], row["stored_dec"], unit="deg")
            )
            assert separation.to_value(u.mas) < 1.0, row["name"]
            assert found["mag"] == pytest.approx(row["stored_mag"], abs=0.005)
