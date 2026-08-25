"""Tests for querying GaiaVarStar and DESI once for many obsids at a time.

Both are queried per obsid today, one network round trip each, and both round
trips are epoch-independent: the GaiaVarStar ADQL draws its cones around the
input positions against Gaia's J2016 catalog positions with no epoch in the
query, and DESI asks VizieR for plain RAICRS/DEICRS. So the positions of many
obsids can go in one query and be split up afterwards, which is what turns
8,451 round trips into a few hundred.

The proper-motion correction is what makes GaiaVarStar the harder of the two: it
is per-obsid, applied after the query, so the split has to happen before the
correction and the correction has to use each obsid's own epoch.
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
    """Raw Gaia rows shaped like the var-star TAP query's output."""
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


# ─── GaiaVarStar ─────────────────────────────────────────────────────────────


def test_gaia_var_stars_by_obsid_uses_one_query_for_every_obsid():
    """The whole point: round trips scale with batches, not with obsids."""
    calls = []

    def stub(query):
        calls.append(query)
        return GAIA_ROWS

    with patch.object(cross_match, "_execute_gaia_tap_query", stub):
        result = cross_match.get_gaia_var_stars_by_obsid(
            _positions_by_obsid(), radius=RADIUS, batch_size=500
        )

    assert len(calls) == 1, "two obsids must not cost two round trips"
    assert set(result) == {7001, 7002}


def test_gaia_var_stars_by_obsid_gives_each_obsid_only_its_own_field():
    with patch.object(cross_match, "_execute_gaia_tap_query", lambda q: GAIA_ROWS):
        result = cross_match.get_gaia_var_stars_by_obsid(
            _positions_by_obsid(), radius=RADIUS, batch_size=500
        )

    assert list(result[7001]["name"]) == ["GaiaVarStar-1001"]
    assert list(result[7002]["name"]) == ["GaiaVarStar-2002"]


def test_gaia_var_stars_by_obsid_corrects_to_each_obsids_own_epoch():
    """A shared source seen at two epochs must come back at two positions.

    This is the trap in batching: the query is epoch-free, but the answer is not.
    Correcting the pooled result once would stamp one obsid's epoch on all of them.
    """
    moving = _gaia_rows([(3003, 150.0, 2.0, 500.0, 500.0, 14.0)])  # 0.5 arcsec/yr
    positions = {
        7001: coords.SkyCoord(FIELD_A, obstime=CxoTime("2020:001")),
        7002: coords.SkyCoord(FIELD_A, obstime=CxoTime("2010:001")),
    }

    with patch.object(cross_match, "_execute_gaia_tap_query", lambda q: moving):
        result = cross_match.get_gaia_var_stars_by_obsid(
            positions, radius=RADIUS, batch_size=500
        )

    assert len(result[7001]) == 1
    assert len(result[7002]) == 1
    # 10 years of separation at 0.5 arcsec/yr in dec is ~5 arcsec of difference.
    separation = abs(result[7001]["dec"][0] - result[7002]["dec"][0]) * 3600
    assert separation == pytest.approx(5.0, abs=0.2)


def test_gaia_var_stars_by_obsid_shares_a_source_between_overlapping_fields():
    """Chandra revisits fields; the same var star is a candidate for each visit."""
    positions = {
        7001: coords.SkyCoord(FIELD_A, obstime=CxoTime("2020:001")),
        7002: coords.SkyCoord(FIELD_A, obstime=CxoTime("2020:001")),
    }

    with patch.object(cross_match, "_execute_gaia_tap_query", lambda q: GAIA_ROWS):
        result = cross_match.get_gaia_var_stars_by_obsid(
            positions, radius=RADIUS, batch_size=500
        )

    assert list(result[7001]["name"]) == ["GaiaVarStar-1001"]
    assert list(result[7002]["name"]) == ["GaiaVarStar-1001"]


def test_gaia_var_stars_by_obsid_agrees_with_the_per_obsid_function():
    """The regression guard: batching must not change what an obsid gets.

    One obsid, so the pooled positions are exactly that obsid's positions and the
    stub -- which cannot filter by position the way the real ADQL does -- is fed
    the same thing either way.
    """
    field_a_only = _gaia_rows([(1001, 150.0, 2.0, 400.0, -300.0, 15.0)])
    pos = coords.SkyCoord(FIELD_A, obstime=CxoTime("2020:001"))

    with patch.object(cross_match, "_execute_gaia_tap_query", lambda q: field_a_only):
        batched = cross_match.get_gaia_var_stars_by_obsid(
            {7001: pos}, radius=RADIUS, batch_size=500
        )
        one = cross_match.get_gaia_var_stars(pos, radius=RADIUS)

    assert list(batched[7001]["name"]) == list(one["name"])
    assert np.allclose(batched[7001]["ra"], one["ra"])
    assert np.allclose(batched[7001]["dec"], one["dec"])
    assert np.allclose(batched[7001]["mag"], one["mag"])


def test_gaia_var_stars_by_obsid_returns_an_empty_table_per_obsid_when_nothing_matches():
    with patch.object(cross_match, "_execute_gaia_tap_query", lambda q: _gaia_rows([])):
        result = cross_match.get_gaia_var_stars_by_obsid(
            _positions_by_obsid(), radius=RADIUS, batch_size=500
        )

    assert set(result) == {7001, 7002}
    assert all(len(candidates) == 0 for candidates in result.values())


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
    """One obsid, for the same reason as the GaiaVarStar equivalence test."""
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


def test_gaia_var_stars_raises_when_a_tap_batch_fails():
    pos = coords.SkyCoord(FIELD_A, obstime=CxoTime("2020:001"))

    with (
        patch.object(cross_match, "_execute_gaia_tap_query", _boom),
        pytest.raises(cross_match.CatalogQueryFailed, match="GaiaVarStar"),
    ):
        cross_match.get_gaia_var_stars(pos, radius=RADIUS)


def test_gaia_var_stars_by_obsid_raises_rather_than_returning_empty_tables():
    with (
        patch.object(cross_match, "_execute_gaia_tap_query", _boom),
        pytest.raises(cross_match.CatalogQueryFailed),
    ):
        cross_match.get_gaia_var_stars_by_obsid(
            _positions_by_obsid(), radius=RADIUS, batch_size=100
        )


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


def test_gaia_var_stars_splits_positions_into_bounded_sub_batches():
    """A query with ~900 OR'd CONTAINS clauses is what broke the first full run."""
    sizes = []

    def counting(query):
        sizes.append(query.count("CONTAINS("))
        return _gaia_rows([])

    with patch.object(cross_match, "_execute_gaia_tap_query", counting):
        cross_match.get_gaia_var_stars_by_obsid(
            _many_positions(240), radius=RADIUS, batch_size=100
        )

    assert len(sizes) == 3, f"240 positions at 100 per query should be 3, got {sizes}"
    assert max(sizes) <= 100


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


def test_gaia_var_stars_accepts_a_scalar_position():
    """The docstring has always claimed scalar support; the batching loop broke it.

    `len()` on a scalar SkyCoord raises TypeError, so the old loop over
    `range(0, len(pos), batch_size)` could not run. Splitting the query out fixed
    it incidentally by looping over the position array instead -- this pins that,
    since nothing else would notice it regressing.
    """
    scalar = coords.SkyCoord(150.0, 2.0, unit="deg", obstime=CxoTime("2020:001"))
    assert scalar.isscalar

    with patch.object(cross_match, "_execute_gaia_tap_query", lambda query: GAIA_ROWS):
        result = cross_match.get_gaia_var_stars(scalar, radius=RADIUS)

    # Both stub rows come back: the radius cut lives in the ADQL, which the stub
    # replaces, and get_gaia_var_stars does not re-filter what the query returned.
    # What is under test is that a scalar reaches the query at all.
    assert sorted(result["name"]) == ["GaiaVarStar-1001", "GaiaVarStar-2002"]
