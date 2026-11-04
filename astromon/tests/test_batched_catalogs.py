"""Tests for querying DESI once for many obsids at a time.

DESI is queried per obsid today, one network round trip each, and the round trip
is epoch-independent: DESI asks VizieR for plain RAICRS/DEICRS. So the positions of
many obsids can go in one query and be split up afterwards, which is what turns
thousands of round trips into a few hundred.
"""

from unittest.mock import patch

import numpy as np
import pytest
import requests
from astropy import coordinates as coords
from astropy import units as u
from astropy.table import Table
from cxotime import CxoTime

from astromon import cross_match

RADIUS = 3 * u.arcsec

# Two well-separated fields, so a source in one is nowhere near the other.
FIELD_A = coords.SkyCoord([150.0, 150.001], [2.0, 2.001], unit="deg")
FIELD_B = coords.SkyCoord([220.0], [-15.0], unit="deg")


def _gaia_rows(entries):
    """Raw Gaia-like rows (source_id, ra, dec, proper motion, G) for the radius filter."""
    return Table(
        {
            "source_id": np.array([e[0] for e in entries], dtype=np.int64),
            "ra": np.array([e[1] for e in entries], dtype=float),
            "dec": np.array([e[2] for e in entries], dtype=float),
            "pmra": np.array([e[3] for e in entries], dtype=float),
            "pmdec": np.array([e[4] for e in entries], dtype=float),
            "phot_g_mean_mag": np.array([e[5] for e in entries], dtype=np.float32),
        }
    )


# One source in field A, one in field B. Zero proper motion unless a test wants it.
GAIA_ROWS = _gaia_rows(
    [
        (1001, 150.0, 2.0, 0.0, 0.0, 15.0),
        (2002, 220.0, -15.0, 0.0, 0.0, 16.0),
    ]
)


def _positions_by_obsid():
    return {
        7001: coords.SkyCoord(FIELD_A, obstime=CxoTime("2020:001")),
        7002: coords.SkyCoord(FIELD_B, obstime=CxoTime("2010:001")),
    }


# ─── the shared radius filter ────────────────────────────────────────────────


def test_rows_near_positions_keeps_only_what_is_within_radius():
    catalog = Table(
        {
            "ra": np.array([150.0, 150.0, 220.0]),
            "dec": np.array([2.0, 2.5, -15.0]),  # second one is 0.5 deg away
        }
    )

    keep = cross_match.rows_near_positions(catalog, FIELD_A, RADIUS)

    assert list(keep) == [0]


def test_rows_near_positions_is_a_union_over_the_positions():
    """A source near any one of the positions is kept, not only the first."""
    catalog = Table({"ra": np.array([150.001]), "dec": np.array([2.001])})

    keep = cross_match.rows_near_positions(catalog, FIELD_A, RADIUS)

    assert list(keep) == [0]


def test_rows_near_positions_handles_an_empty_catalog():
    catalog = Table({"ra": np.array([]), "dec": np.array([])})

    assert len(cross_match.rows_near_positions(catalog, FIELD_A, RADIUS)) == 0


# ─── DESI ────────────────────────────────────────────────────────────────────


def _desi_tablelist():
    """A VizieR result shaped like V/161, one target per field."""
    return [
        Table(
            {
                "RAICRS": np.array([150.0, 220.0]),
                "DEICRS": np.array([2.0, -15.0]),
                "TargetID": np.array([111, 222], dtype=np.int64),
                "OType": np.array(["QSO", "GALAXY"]),
                "ZWARN": np.array([0, 0]),
                "delChi2": np.array([500.0, 500.0]),
            }
        )
    ]


def test_desi_by_obsid_uses_one_query_for_every_obsid():
    calls = []

    def stub(vizier, pos, radius, cat_identifier):
        calls.append(cat_identifier)
        return _desi_tablelist()

    with patch.object(cross_match, "_query_vizier_region", stub):
        result = cross_match.get_desi_v161_by_obsid(
            _positions_by_obsid(), radius=RADIUS
        )

    assert len(calls) == 1
    assert set(result) == {7001, 7002}


def test_desi_by_obsid_gives_each_obsid_only_its_own_field():
    with patch.object(
        cross_match, "_query_vizier_region", lambda *a: _desi_tablelist()
    ):
        result = cross_match.get_desi_v161_by_obsid(
            _positions_by_obsid(), radius=RADIUS
        )

    assert list(result[7001]["name"]) == ["DESIV161-111"]
    assert list(result[7002]["name"]) == ["DESIV161-222"]


def test_desi_by_obsid_agrees_with_the_per_obsid_function():
    """One obsid, so the batched and per-obsid results are directly comparable."""
    field_a_only = [
        Table(
            {
                "RAICRS": np.array([150.0]),
                "DEICRS": np.array([2.0]),
                "TargetID": np.array([111], dtype=np.int64),
                "OType": np.array(["QSO"]),
                "ZWARN": np.array([0]),
                "delChi2": np.array([500.0]),
            }
        )
    ]
    pos = coords.SkyCoord(FIELD_A, obstime=CxoTime("2020:001"))

    with patch.object(cross_match, "_query_vizier_region", lambda *a: field_a_only):
        batched = cross_match.get_desi_v161_by_obsid({7001: pos}, radius=RADIUS)
        one = cross_match.get_desi_v161_candidates(pos, radius=RADIUS)

    assert list(batched[7001]["name"]) == list(one["name"])
    assert np.allclose(batched[7001]["ra"], one["ra"])


# ─── a failing query must not look like an empty sky ─────────────────────────
#
# These exist because the first full run exposed what stubs cannot: under load,
# Gaia TAP retried and VizieR read-timed-out, and both getters turned that into
# "no candidates". Pooled across a hundred obsids, one timeout would have written
# an authoritative-looking absence for all of them -- the same defect this whole
# branch is repairing, one layer above astroquery and with a far larger reach.


def _boom(*args, **kwargs):
    raise requests.exceptions.ReadTimeout("Read timed out.")


def test_desi_candidates_raises_when_the_vizier_query_fails():
    pos = coords.SkyCoord(FIELD_A, obstime=CxoTime("2020:001"))

    with (
        patch.object(cross_match, "_query_vizier_region", _boom),
        pytest.raises(cross_match.CatalogQueryFailed, match="DESIV161"),
    ):
        cross_match.get_desi_v161_candidates(pos, radius=RADIUS)


def test_desi_by_obsid_raises_rather_than_returning_empty_tables():
    with (
        patch.object(cross_match, "_query_vizier_region", _boom),
        pytest.raises(cross_match.CatalogQueryFailed),
    ):
        cross_match.get_desi_v161_by_obsid(_positions_by_obsid(), radius=RADIUS)


# ─── sub-batching, so one timeout costs one sub-batch ────────────────────────


def _many_positions(count, obsids=4):
    """`count` positions spread over `obsids` keys, as a big chunk would produce."""
    per = count // obsids
    return {
        7000 + index: coords.SkyCoord(
            np.linspace(100.0 + index, 100.5 + index, per),
            np.linspace(10.0, 10.5, per),
            unit="deg",
            obstime=CxoTime("2020:001"),
        )
        for index in range(obsids)
    }


def test_desi_splits_positions_into_bounded_sub_batches():
    counts = []

    def counting(vizier, pos, radius, cat_identifier):
        counts.append(len(pos))
        return []

    with patch.object(cross_match, "_query_vizier_region", counting):
        cross_match.get_desi_v161_by_obsid(
            _many_positions(240), radius=RADIUS, batch_size=100
        )

    assert len(counts) == 3, f"240 positions at 100 per query should be 3, got {counts}"
    assert max(counts) <= 100


def test_batched_positions_per_query_is_well_under_what_the_services_refused():
    """900 positions per query timed out on both services; 500 was the old default."""
    assert cross_match.BATCHED_POSITIONS_PER_QUERY <= 100
