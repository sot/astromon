"""Tests for re-querying catalog candidates from positions already in the database.

The re-query exists because astromon_cat_src is missing catalogs for an unknown
set of obsids: a candidate set that was never queried leaves no trace, so the
only way to find out what is missing is to ask the catalogs again.  These tests
cover the two things that make that safe -- allocating ids that do not collide
with the rows already stored, and writing one catalog at a time so a re-query of
RFC cannot disturb the obsid's 2MASS rows.
"""

import numpy as np
import pytest
from astropy.table import Table, vstack
from Quaternion import Quat

from astromon import db
from astromon.scripts.maintenance import requery_cat_src

# A pointing quaternion for an arbitrary field; only used to turn ra/dec into
# y/z angles, so its exact value does not matter to any assertion below.
QUAT = Quat(equatorial=(150.0, 2.0, 30.0))


def _xray_sources(obsid=7001, ids=(3, 7), ras=(150.0, 150.01), decs=(2.0, 2.01)):
    """A minimal astromon_xray_src-like table: the columns the re-query reads."""
    return Table(
        {
            "obsid": np.full(len(ids), obsid, dtype=np.int32),
            "id": np.array(ids, dtype=np.int32),
            "ra": np.array(ras, dtype=np.float64),
            "dec": np.array(decs, dtype=np.float64),
            "detect_method": np.array(["celldetect"] * len(ids), dtype="U24"),
        }
    )


def _candidates(ras, decs, names=None, catalog="RFC"):
    """Catalog candidates as the cross_match getters return them."""
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


def _cat_src_rows(rows):
    """Build an astromon_cat_src table from (obsid, id, catalog, name) tuples."""
    out = Table(np.zeros(len(rows), dtype=db.ASTROMON_CAT_SRC_DTYPE))
    for i, (obsid, id_, catalog, name) in enumerate(rows):
        out["obsid"][i] = obsid
        out["id"][i] = id_
        out["catalog"][i] = catalog
        out["name"][i] = name
    return out

# ─── reading positions back out of the database ──────────────────────────────


def test_celldetect_sources_selects_one_obsid_and_method():
    xray = vstack(
        [
            _xray_sources(obsid=7001, ids=(1, 2)),
            _xray_sources(obsid=7002, ids=(1,), ras=(151.0,), decs=(3.0,)),
        ]
    )
    xray["detect_method"][1] = "gaussian_detect"

    sources = requery_cat_src.celldetect_sources(xray, 7001)

    assert list(sources["id"]) == [1]
    assert set(sources["obsid"]) == {7001}


def test_celldetect_sources_raises_when_obsid_absent():
    with pytest.raises(ValueError, match="no celldetect sources"):
        requery_cat_src.celldetect_sources(_xray_sources(obsid=7001), 9999)


def test_next_cat_src_id_is_above_every_stored_id():
    """Ids only have to be unique within an obsid, and existing c_id references
    must keep pointing at the same rows, so new rows go above the current max."""
    cat = _cat_src_rows(
        [
            (7001, 0, "2MASS", "a"),
            (7001, 9, "USNO-B1.0", "b"),
            (7002, 40, "2MASS", "c"),
        ]
    )

    assert requery_cat_src.next_cat_src_id(cat, 7001) == 10
    assert requery_cat_src.next_cat_src_id(cat, 7002) == 41
    assert requery_cat_src.next_cat_src_id(cat, 7003) == 0


# ─── assembling one obsid's re-query ─────────────────────────────────────────


def test_requery_obsid_stacks_catalogs_with_unique_ids():
    """Each catalog is queried separately but shares one id sequence per obsid."""
    getters = {
        "RFC": lambda sources, obs_time, logging_tag="": _candidates(
            [150.0], [2.0], catalog="RFC"
        ),
        "Quaia": lambda sources, obs_time, logging_tag="": _candidates(
            [150.01, 150.01], [2.01, 2.01], catalog="Quaia"
        ),
    }

    candidates = requery_cat_src.requery_obsid(
        obsid=7001,
        sources=_xray_sources(ids=(3, 7)),
        quat=QUAT,
        obs_time="2020:001",
        catalogs=("RFC", "Quaia"),
        id_offset=10,
        getters=getters,
    )

    assert list(candidates["catalog"]) == ["RFC", "Quaia", "Quaia"]
    assert list(candidates["id"]) == [10, 11, 12]
    assert list(candidates["celldetect_x_id"]) == [3, 7, 7]
    assert set(candidates.colnames) >= set(db.ASTROMON_CAT_SRC_DTYPE.names)


def test_requery_obsid_skips_catalogs_with_no_candidates():
    getters = {
        "RFC": lambda sources, obs_time, logging_tag="": _candidates([], []),
        "Quaia": lambda sources, obs_time, logging_tag="": _candidates(
            [150.0], [2.0], catalog="Quaia"
        ),
    }

    candidates = requery_cat_src.requery_obsid(
        obsid=7001,
        sources=_xray_sources(),
        quat=QUAT,
        obs_time="2020:001",
        catalogs=("RFC", "Quaia"),
        id_offset=0,
        getters=getters,
    )

    assert list(candidates["catalog"]) == ["Quaia"]
    assert list(candidates["id"]) == [0]


def test_requery_obsid_returns_empty_table_when_nothing_found():
    getters = {"RFC": lambda sources, obs_time, logging_tag="": _candidates([], [])}

    candidates = requery_cat_src.requery_obsid(
        obsid=7001,
        sources=_xray_sources(),
        quat=QUAT,
        obs_time="2020:001",
        catalogs=("RFC",),
        id_offset=0,
        getters=getters,
    )

    assert len(candidates) == 0
    assert set(candidates.colnames) >= set(db.ASTROMON_CAT_SRC_DTYPE.names)


def test_requery_obsid_reports_a_failing_catalog_without_losing_the_others():
    """One catalog server being down must not throw away the catalogs that worked."""

    def boom(sources, obs_time, logging_tag=""):
        raise RuntimeError("vizier is down")

    getters = {
        "RFC": boom,
        "Quaia": lambda sources, obs_time, logging_tag="": _candidates(
            [150.0], [2.0], catalog="Quaia"
        ),
    }

    with pytest.raises(requery_cat_src.CatalogQueryError, match="RFC"):
        requery_cat_src.requery_obsid(
            obsid=7001,
            sources=_xray_sources(),
            quat=QUAT,
            obs_time="2020:001",
            catalogs=("RFC", "Quaia"),
            id_offset=0,
            getters=getters,
        )


# ─── writing one catalog at a time ───────────────────────────────────────────


def _seeded_db(tmp_path):
    dbfile = tmp_path / "astromon.h5"
    db.create_empty_tables(dbfile)
    db.save(
        "astromon_cat_src",
        _cat_src_rows(
            [
                (7001, 0, "2MASS", "keep-2mass"),
                (7001, 1, "RFC", "stale-rfc"),
                (7002, 0, "RFC", "other-obsid-rfc"),
            ]
        ),
        dbfile,
        ignore_obsid=True,
    )
    return dbfile


def test_write_candidates_replaces_only_the_requeried_catalog(tmp_path):
    """This is the defect the whole script exists to avoid.

    The earlier backfills rewrote astromon_cat_src wholesale, so a run that was
    only rebuilding one catalog decided the fate of all of them.  Keying the
    write on (obsid, catalog) confines it to what was actually re-queried.
    """
    dbfile = _seeded_db(tmp_path)
    new = _cat_src_rows([(7001, 10, "RFC", "fresh-rfc")])

    requery_cat_src.write_candidates(dbfile, new)

    stored = db.get_table("astromon_cat_src", dbfile)
    by_name = {row["name"]: row for row in stored}
    assert set(by_name) == {"keep-2mass", "fresh-rfc", "other-obsid-rfc"}
    assert by_name["fresh-rfc"]["id"] == 10


def test_write_candidates_leaves_untouched_catalogs_alone(tmp_path):
    """A catalog that returned nothing is not evidence that its rows are wrong."""
    dbfile = _seeded_db(tmp_path)

    requery_cat_src.write_candidates(
        dbfile, _cat_src_rows([(7001, 10, "Quaia", "new-quaia")])
    )

    stored = db.get_table("astromon_cat_src", dbfile)
    assert "stale-rfc" in list(stored["name"])
    assert "keep-2mass" in list(stored["name"])


def test_write_candidates_dry_run_writes_nothing(tmp_path):
    dbfile = _seeded_db(tmp_path)
    before = db.get_table("astromon_cat_src", dbfile)

    requery_cat_src.write_candidates(
        dbfile, _cat_src_rows([(7001, 10, "RFC", "fresh-rfc")]), dry_run=True
    )

    after = db.get_table("astromon_cat_src", dbfile)
    assert list(after["name"]) == list(before["name"])


def test_write_candidates_requires_an_existing_table(tmp_path):
    """A missing table means a previous write was killed mid-flight, not a new db."""
    dbfile = tmp_path / "empty.h5"

    with pytest.raises(db.MissingTableException):
        requery_cat_src.write_candidates(dbfile, _cat_src_rows([(7001, 0, "RFC", "a")]))


def test_write_candidates_reports_obsids_needing_an_xcorr_rebuild(tmp_path):
    """Replacing a catalog renumbers its rows, so xcorr for that obsid is stale.

    Adding a catalog the obsid never had cannot invalidate anything, because the
    new ids are allocated above every id an xcorr row could already reference.
    """
    dbfile = _seeded_db(tmp_path)
    new = vstack(
        [
            _cat_src_rows([(7001, 10, "RFC", "fresh-rfc")]),  # replaces existing RFC
            _cat_src_rows([(7002, 10, "Quaia", "added")]),  # obsid 7002 had no Quaia
        ]
    )

    result = requery_cat_src.write_candidates(dbfile, new)

    assert result["replaced_obsids"] == [7001]
    assert result["added_rows"] == 2
