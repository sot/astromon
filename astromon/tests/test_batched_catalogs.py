"""Tests for querying DESI once for many obsids at a time.

DESIV161 is now served from a local cache (astromon.cross_match's
``get_desi_v161_catalog``, built from DESI's own DR1 data release -- see
``_build_desi_v161_catalog``) rather than live VizieR, so the batching and
sub-batching ``get_desi_v161_by_obsid`` still does is no longer needed for
timeout avoidance; it is kept only because splitting pooled positions back out
per obsid is still useful bookkeeping. The tests below that existed purely to
pin VizieR-timeout-avoidance behavior (one pooled query, bounded sub-batches, a
failed query raising rather than looking like an empty sky) were removed along
with that VizieR call; what is left checks the still-relevant per-obsid
splitting.
"""

from unittest.mock import patch

import numpy as np
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


def _desi_catalog():
    """A local DESIV161 catalog shaped like the real cache, one source per field."""
    return Table(
        {
            "name": ["DESIV161-111", "DESIV161-222"],
            "ra": np.array([150.0, 220.0]),
            "dec": np.array([2.0, -15.0]),
        }
    )


def test_desi_by_obsid_gives_each_obsid_only_its_own_field():
    with patch.object(cross_match, "get_desi_v161_catalog", _desi_catalog):
        result = cross_match.get_desi_v161_by_obsid(
            _positions_by_obsid(), radius=RADIUS
        )

    assert list(result[7001]["name"]) == ["DESIV161-111"]
    assert list(result[7002]["name"]) == ["DESIV161-222"]


def test_desi_by_obsid_agrees_with_the_per_obsid_function():
    """One obsid, so the batched and per-obsid results are directly comparable."""
    field_a_only = Table(
        {"name": ["DESIV161-111"], "ra": np.array([150.0]), "dec": np.array([2.0])}
    )
    pos = coords.SkyCoord(FIELD_A, obstime=CxoTime("2020:001"))

    with patch.object(cross_match, "get_desi_v161_catalog", lambda: field_a_only):
        batched = cross_match.get_desi_v161_by_obsid({7001: pos}, radius=RADIUS)
        one = cross_match.get_desi_v161_candidates(pos, radius=RADIUS)

    assert list(batched[7001]["name"]) == list(one["name"])
    assert np.allclose(batched[7001]["ra"], one["ra"])


def test_batched_positions_per_query_is_well_under_what_the_services_refused():
    """900 positions per query timed out on both services; 500 was the old default."""
    assert cross_match.BATCHED_POSITIONS_PER_QUERY <= 100
