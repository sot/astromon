"""Tests for re-querying catalog candidates from positions already in the database.

The re-query exists because astromon_cat_src is missing catalogs for an unknown
set of obsids: a candidate set that was never queried leaves no trace, so the
only way to find out what is missing is to ask the catalogs again.  These tests
cover the two things that make that safe -- allocating ids that do not collide
with the rows already stored, and writing one catalog at a time so a re-query of
RFC cannot disturb the obsid's 2MASS rows.
"""

from unittest.mock import patch

import numpy as np
import pytest
from astropy.table import Table, vstack
from Quaternion import Quat

from astromon import cross_match, db
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


def test_requery_records_the_catalog_releases_it_queried_against(tmp_path):
    """A run's output should be traceable to the catalog releases behind it.

    astromon_cat_src has no per-row provenance, so the run summary is where this
    gets written down. It would not have detected the rows missing from the
    current database -- absence cannot be annotated -- but it makes the next
    question ("which release produced these?") answerable at all.
    """
    dbfile = _seeded_db(tmp_path)

    summary = requery_cat_src.requery(dbfile, [], catalogs=("RFC",))

    assert set(summary["catalog_versions"]) == set(cross_match.CATALOG_CACHE_PATHS)
# ─── catalog classification ──────────────────────────────────────────────────


def test_milliquas_gaia_is_served_locally_now():
    """It is a local lookup since the Gaia upgrade became a precomputed mapping."""
    assert "MilliquasGaia" in requery_cat_src.LOCAL_CATALOGS
    assert "MilliquasGaia" not in requery_cat_src.REMOTE_CATALOGS


def test_batched_catalogs_are_declared_apart_from_the_per_obsid_ones():
    """DESI is queried once for a whole chunk of obsids rather than once each.

    GaiaVarStar is queried the same way but is not in the default set at all -- see
    test_gaia_var_star_is_not_matched_by_default.
    """
    assert set(requery_cat_src.BATCHED_CATALOGS) == {"DESIV161"}
    assert set(requery_cat_src.BATCHED_CATALOGS) <= set(requery_cat_src.ALL_CATALOGS)
    assert set(requery_cat_src.BATCHED_GETTERS) == {"DESIV161", "GaiaVarStar"}
    # The four rough_match catalogs stay per obsid: their VizieR query asks for
    # positions at the observation epoch, so it cannot be pooled.
    assert set(requery_cat_src.REMOTE_CATALOGS) == {
        "Tycho2",
        "USNO-B1.0",
        "2MASS",
        "SDSS",
    }


# ─── folding a batched result into one obsid's candidates ────────────────────


def test_requery_obsid_uses_precomputed_candidates_instead_of_querying():
    def explode(*args, **kwargs):
        raise AssertionError("a precomputed catalog must not be queried again")

    candidates = requery_cat_src.requery_obsid(
        obsid=7001,
        sources=_xray_sources(),
        quat=QUAT,
        obs_time="2020:001",
        catalogs=("RFC",),
        id_offset=0,
        getters={"RFC": explode},
        precomputed={"RFC": _candidates([150.0], [2.0], catalog="RFC")},
    )

    assert list(candidates["catalog"]) == ["RFC"]


def test_requery_obsid_numbers_precomputed_and_queried_in_one_sequence():
    """Ids run consecutively per obsid however the candidates were obtained."""
    getters = {
        "Quaia": lambda sources, obs_time, logging_tag="": _candidates(
            [150.01], [2.01], catalog="Quaia"
        )
    }

    candidates = requery_cat_src.requery_obsid(
        obsid=7001,
        sources=_xray_sources(ids=(3, 7)),
        quat=QUAT,
        obs_time="2020:001",
        catalogs=("GaiaVarStar", "Quaia"),
        id_offset=5,
        getters=getters,
        precomputed={"GaiaVarStar": _candidates([150.0], [2.0], catalog="GaiaVarStar")},
    )

    assert list(candidates["catalog"]) == ["GaiaVarStar", "Quaia"]
    assert list(candidates["id"]) == [5, 6]


def test_requery_obsid_does_not_demand_a_getter_for_a_precomputed_catalog():
    candidates = requery_cat_src.requery_obsid(
        obsid=7001,
        sources=_xray_sources(),
        quat=QUAT,
        obs_time="2020:001",
        catalogs=("GaiaVarStar",),
        id_offset=0,
        getters={},
        precomputed={"GaiaVarStar": _candidates([150.0], [2.0], catalog="GaiaVarStar")},
    )

    assert len(candidates) == 1
# ─── resuming a long run ─────────────────────────────────────────────────────


def test_progress_records_which_catalogs_each_obsid_finished(tmp_path):
    progress = tmp_path / "done.txt"

    requery_cat_src.record_completed(progress, {7001: ["RFC", "Quaia"]})
    requery_cat_src.record_completed(progress, {7002: ["RFC"]})

    assert requery_cat_src.load_completed(progress) == {
        7001: {"RFC", "Quaia"},
        7002: {"RFC"},
    }


def test_progress_accumulates_catalogs_across_runs(tmp_path):
    """A second pass adds catalogs to an obsid rather than replacing them."""
    progress = tmp_path / "done.txt"

    requery_cat_src.record_completed(progress, {7001: ["RFC"]})
    requery_cat_src.record_completed(progress, {7001: ["DESIV161"]})

    assert requery_cat_src.load_completed(progress) == {7001: {"RFC", "DESIV161"}}


def test_load_completed_tolerates_an_absent_file(tmp_path):
    assert requery_cat_src.load_completed(tmp_path / "nope.txt") == {}


def test_load_completed_ignores_a_torn_final_line(tmp_path):
    """A run killed mid-append leaves a partial record; only whole ones count."""
    progress = tmp_path / "done.txt"
    progress.write_text("7001 RFC\n7002 RFC\n7003 QU")

    assert requery_cat_src.load_completed(progress) == {7001: {"RFC"}, 7002: {"RFC"}}


def test_an_obsid_is_skipped_only_when_every_requested_catalog_is_done():
    completed = {7001: {"RFC", "Quaia"}, 7002: {"RFC"}}

    assert requery_cat_src.is_complete(completed, 7001, ("RFC",))
    assert requery_cat_src.is_complete(completed, 7001, ("RFC", "Quaia"))
    # 7002 never got Quaia, so a run that wants Quaia must not skip it.
    assert not requery_cat_src.is_complete(completed, 7002, ("RFC", "Quaia"))
    assert not requery_cat_src.is_complete(completed, 7003, ("RFC",))


def test_a_legacy_obsid_only_record_is_not_treated_as_complete(tmp_path):
    """Bare-obsid files predate this and say nothing about which catalogs ran.

    Assuming they covered everything is the exact mistake this change removes, so
    they count as unknown and get re-queried.
    """
    progress = tmp_path / "legacy.txt"
    progress.write_text("7001\n7002\n")

    completed = requery_cat_src.load_completed(progress)

    assert completed == {}
    assert not requery_cat_src.is_complete(completed, 7001, ("RFC",))


def test_requery_skips_only_obsids_whose_requested_catalogs_are_all_done(tmp_path):
    dbfile = _seeded_db(tmp_path)
    progress = tmp_path / "done.txt"
    requery_cat_src.record_completed(progress, {7001: ["RFC"], 7002: ["DESIV161"]})

    summary = requery_cat_src.requery(
        dbfile,
        [7001, 7002],
        catalogs=("RFC",),
        progress_file=progress,
        resume=True,
    )

    # 7001 already has RFC; 7002 has only DESIV161, so it still needs RFC.
    assert summary["resumed"] == 1


# ─── GaiaVarStar is out of the default set ───────────────────────────────────


def test_gaia_var_star_is_not_matched_by_default():
    """It costs more than everything it finds.

    99 rows across the whole database, from 96 obsids, and its Gaia TAP query owns
    every stall in a full run: a 120 s hard timeout inside tries=3/delay=10/backoff=2
    is up to ~6.5 minutes per bad batch, and a measured chunk of 100 obsids hit four
    of them.
    """
    assert "GaiaVarStar" not in requery_cat_src.ALL_CATALOGS
    assert "GaiaVarStar" not in requery_cat_src.BATCHED_CATALOGS
    assert "GaiaVarStar" in requery_cat_src.OPTIONAL_CATALOGS


def test_catalogs_all_excludes_gaia_var_star():
    resolved = requery_cat_src._resolve_catalogs(["all"])

    assert "GaiaVarStar" not in resolved
    assert "DESIV161" in resolved
    assert "MilliquasGaia" in resolved


def test_gaia_var_star_can_still_be_asked_for_by_name():
    """Excluded from the default set, not removed -- a later pass can still run it."""
    assert requery_cat_src._resolve_catalogs(["GaiaVarStar"]) == ["GaiaVarStar"]
    assert "GaiaVarStar" in requery_cat_src.CATALOG_GETTERS
    assert "GaiaVarStar" in requery_cat_src.BATCHED_GETTERS


def test_named_gaia_var_star_still_goes_through_the_batched_path():
    """Asking for it explicitly must not fall back to a query per obsid."""
    calls = []

    def batched(positions, radius=None, logging_tag=""):
        calls.append(sorted(positions))
        return {obsid: _candidates([], []) for obsid in positions}

    with patch.dict(
        requery_cat_src.BATCHED_GETTERS, {"GaiaVarStar": batched}, clear=False
    ):
        result = requery_cat_src._batched_candidates(
            ["GaiaVarStar"], {7001: None, 7002: None}
        )

    assert calls == [[7001, 7002]]
    assert set(result) == {"GaiaVarStar"}
