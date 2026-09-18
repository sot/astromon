"""Tests for rebuild_xcorr, which recomputes astromon_xcorr from stored data.

Restored from backup/stack-09-before-cleanup-merge alongside the script itself
(astromon/scripts/maintenance/rebuild_xcorr.py) -- both were dropped in a
cleanup rebase that was meant for genuinely one-off migration scripts but swept
this general-purpose repair tool up with them. It is still needed: a bulk
detection-only backfill (--skip-catalog-match) followed by
requery_cat_src.py only ever populates astromon_cat_src candidates, never
astromon_xcorr -- this is the only script that computes the latter from data
already in the database, without re-running detection.
"""

import sys
from unittest.mock import patch

import numpy as np
import pytest
from astropy.table import Table, vstack

from astromon import db
from astromon.tests.test_db import _cat_src_row, _xcorr_row

# --- rebuilding xcorr from stored data --------------------------------------


def test_main_passes_force_to_rebuild(tmp_path):
    from astromon.scripts.maintenance import rebuild_xcorr

    dbfile = tmp_path / "astromon.h5"
    dbfile.touch()
    argv = ["rebuild_xcorr", "--db", str(dbfile), "--force"]

    with (
        patch.object(sys, "argv", argv),
        patch.object(rebuild_xcorr, "rebuild") as rebuild,
    ):
        rebuild_xcorr.main()

    rebuild.assert_called_once_with(
        dbfile, obsids=None, select_names=None, dry_run=False, force=True
    )


def test_find_obsids_with_dangling_c_id_flags_only_unresolvable_rows():
    """A c_id naming no cat_src row is broken; a differing anchor is not."""
    from astromon.scripts.maintenance.rebuild_xcorr import (
        find_obsids_with_dangling_c_id,
    )

    cat = vstack(
        [
            _cat_src_row(catalog="RFC", obsid=7001, celldetect_x_id=5, name="a"),
            _cat_src_row(catalog="RFC", obsid=7002, celldetect_x_id=9, name="b"),
        ],
        metadata_conflicts="silent",
    )
    cat["id"] = [1, 1]
    xcorr = vstack(
        [
            _xcorr_row(select_name="astromon_21", obsid=7001, c_id=1, x_id=5),
            # x_id 4 differs from this row's anchor of 9, which is legitimate:
            # the anchor is the nearest celldetect source, not the matched one.
            _xcorr_row(select_name="astromon_21", obsid=7002, c_id=1, x_id=4),
            _xcorr_row(select_name="astromon_21", obsid=7003, c_id=1, x_id=1),
        ],
        metadata_conflicts="silent",
    )

    bad = find_obsids_with_dangling_c_id(xcorr, cat)

    assert 7001 not in bad, "resolvable c_id must not be flagged"
    assert 7002 not in bad, "a differing anchor is not an inconsistency"
    assert bad == {7003: 1}


def test_select_names_in_reads_only_what_is_present():
    """A rebuild must not invent select_names the database never had."""
    from astromon.scripts.maintenance.rebuild_xcorr import select_names_in

    xcorr = vstack(
        [
            _xcorr_row(select_name="astromon_21"),
            _xcorr_row(select_name="rfc"),
            _xcorr_row(select_name="astromon_21"),
        ],
        metadata_conflicts="silent",
    )
    assert select_names_in(xcorr) == ["astromon_21", "rfc"]


def _obs_table(obsid):
    obs = db.create_table("astromon_obs")
    row = Table(np.zeros(1, dtype=db.ASTROMON_OBS_DTYPE))
    row["obsid"] = obsid
    return vstack([obs, row], metadata_conflicts="silent")


def test_rebuild_is_a_noop_when_everything_is_consistent(tmp_path):
    from astromon.scripts.maintenance.rebuild_xcorr import rebuild

    dbfile = tmp_path / "consistent.h5"
    cat = _cat_src_row(catalog="RFC", obsid=7001, celldetect_x_id=5)
    cat["id"] = 1
    db.save("astromon_cat_src", cat, dbfile, ignore_obsid=True)
    db.save(
        "astromon_xcorr",
        _xcorr_row(select_name="astromon_21", obsid=7001, c_id=1, x_id=5),
        dbfile,
        ignore_obsid=True,
    )
    db.save(
        "astromon_xray_src",
        db.create_table("astromon_xray_src"),
        dbfile,
        ignore_obsid=True,
    )
    db.save("astromon_obs", db.create_table("astromon_obs"), dbfile, ignore_obsid=True)

    result = rebuild(dbfile)
    assert result["obsids"] == []
    assert len(db.get_table("astromon_xcorr", dbfile)) == 1, "nothing touched"


def test_find_selections_with_fewer_matches_ignores_equal_or_more_rows():
    """No deficit -- same or more rows after -- is not a loss."""
    from astromon.scripts.maintenance.rebuild_xcorr import (
        find_selections_with_fewer_matches,
    )

    before = vstack(
        [
            _xcorr_row(select_name="astromon_25", obsid=7001, x_id=1),
            _xcorr_row(select_name="astromon_25", obsid=7001, x_id=2),
        ],
        metadata_conflicts="silent",
    )
    same = vstack(
        [
            _xcorr_row(select_name="astromon_25", obsid=7001, x_id=1),
            _xcorr_row(select_name="astromon_25", obsid=7001, x_id=2),
        ],
        metadata_conflicts="silent",
    )
    more = vstack([same, _xcorr_row(select_name="astromon_25", obsid=7001, x_id=3)])

    assert find_selections_with_fewer_matches(before, same) == {}
    assert find_selections_with_fewer_matches(before, more) == {}
    assert find_selections_with_fewer_matches(before, None) == {
        (7001, "astromon_25"): 2
    }


def test_find_selections_with_fewer_matches_flags_a_deficit():
    """The actual harm: a rebuild that would leave fewer rows than exist today.

    This is what a partial cat_src loss looks like in practice -- losing eight of
    astromon_25's nine hierarchy catalogs, say, and keeping only RFC. A
    catalog-presence check can't reliably tell that apart from every catalog but
    RFC simply never having had a candidate here, which is the common case, not
    loss. Comparing row counts sidesteps the question: whatever the cause, two
    matches recomputing down to one is a real deficit.
    """
    from astromon.scripts.maintenance.rebuild_xcorr import (
        find_selections_with_fewer_matches,
    )

    before = vstack(
        [
            _xcorr_row(select_name="astromon_25", obsid=7001, x_id=1),
            _xcorr_row(select_name="astromon_25", obsid=7001, x_id=2),
        ],
        metadata_conflicts="silent",
    )
    after = _xcorr_row(select_name="astromon_25", obsid=7001, x_id=1)

    assert find_selections_with_fewer_matches(before, after) == {
        (7001, "astromon_25"): 1
    }


def test_find_selections_with_fewer_matches_ignores_a_brand_new_selection():
    """A select_name recomputed for the first time is a backfill, not a loss."""
    from astromon.scripts.maintenance.rebuild_xcorr import (
        find_selections_with_fewer_matches,
    )

    before = _xcorr_row(select_name="astromon_21", obsid=7001)
    after = vstack(
        [
            _xcorr_row(select_name="astromon_21", obsid=7001),
            _xcorr_row(select_name="gaia_agn", obsid=7001),
        ],
        metadata_conflicts="silent",
    )
    assert find_selections_with_fewer_matches(before, after) == {}


def test_rebuild_refuses_when_recompute_finds_fewer_matches(tmp_path):
    """rebuild() itself refuses, not just the helper function in isolation.

    Recompute is mocked here to stand in for whatever caused a real deficit (a
    dropped catalog, a stricter cut, anything) -- the check does not need to know
    the cause, only that rows are about to disappear.
    """
    from astromon.scripts.maintenance import rebuild_xcorr

    dbfile = tmp_path / "fewer_matches.h5"
    cat = _cat_src_row(catalog="RFC", obsid=7001, celldetect_x_id=1)
    cat["id"] = 1
    db.save("astromon_cat_src", cat, dbfile, ignore_obsid=True)
    original = vstack(
        [
            _xcorr_row(select_name="astromon_25", obsid=7001, c_id=1, x_id=1),
            _xcorr_row(select_name="astromon_25", obsid=7001, c_id=1, x_id=2),
        ],
        metadata_conflicts="silent",
    )
    db.save("astromon_xcorr", original, dbfile, ignore_obsid=True)
    xray = Table(np.zeros(2, dtype=db.ASTROMON_XRAY_SRC_DTYPE))
    xray["obsid"] = 7001
    xray["id"] = [1, 2]
    xray["detect_method"] = "gaussian_detect"
    db.save("astromon_xray_src", xray, dbfile, ignore_obsid=True)
    db.save("astromon_obs", _obs_table(7001), dbfile, ignore_obsid=True)

    def one_row_only(select_names, obs, xray, cat):
        return _xcorr_row(select_name="astromon_25", obsid=7001, c_id=1, x_id=1)[
            list(db.ASTROMON_XCORR_DTYPE.names)
        ]

    with patch.object(rebuild_xcorr, "recompute", one_row_only):
        with pytest.raises(RuntimeError, match="would be lost, not repaired"):
            rebuild_xcorr.rebuild(dbfile, obsids=[7001])

    stored = db.get_table("astromon_xcorr", dbfile)
    assert len(stored) == 2, "the database must be untouched"


def test_rebuild_with_force_accepts_fewer_matches(tmp_path):
    """force=True is the documented escape hatch -- it must still write."""
    from astromon.scripts.maintenance import rebuild_xcorr

    dbfile = tmp_path / "fewer_matches_forced.h5"
    cat = _cat_src_row(catalog="RFC", obsid=7001, celldetect_x_id=1)
    cat["id"] = 1
    db.save("astromon_cat_src", cat, dbfile, ignore_obsid=True)
    original = vstack(
        [
            _xcorr_row(select_name="astromon_25", obsid=7001, c_id=1, x_id=1),
            _xcorr_row(select_name="astromon_25", obsid=7001, c_id=1, x_id=2),
        ],
        metadata_conflicts="silent",
    )
    db.save("astromon_xcorr", original, dbfile, ignore_obsid=True)
    xray = Table(np.zeros(2, dtype=db.ASTROMON_XRAY_SRC_DTYPE))
    xray["obsid"] = 7001
    xray["id"] = [1, 2]
    xray["detect_method"] = "gaussian_detect"
    db.save("astromon_xray_src", xray, dbfile, ignore_obsid=True)
    db.save("astromon_obs", _obs_table(7001), dbfile, ignore_obsid=True)

    def one_row_only(select_names, obs, xray, cat):
        return _xcorr_row(select_name="astromon_25", obsid=7001, c_id=1, x_id=1)[
            list(db.ASTROMON_XCORR_DTYPE.names)
        ]

    with patch.object(rebuild_xcorr, "recompute", one_row_only):
        rebuild_xcorr.rebuild(dbfile, obsids=[7001], force=True)

    stored = db.get_table("astromon_xcorr", dbfile)
    assert len(stored) == 1


# ─── rebuild_xcorr and the per-method x_id remap ─────────────────────────────


def test_recompute_makes_one_pass_over_every_detect_method():
    """One call per select_name, handed every method at once.

    This replaces three tests for remapped_cat_src, which re-pointed cat_src's
    anchor at each method's source ids and drove a pass per method. That existed
    only because the join keyed on (obsid, x_id), so a single stored anchor could
    join to just one method. The join is on obsid now and the grouping is per
    method, so one pass covers both -- and nothing rewrites the anchor.
    """
    from astromon.scripts.maintenance import rebuild_xcorr

    xray = db.create_table("astromon_xray_src")
    for method in ("celldetect", "gaussian_detect"):
        row = Table(np.zeros(1, dtype=db.ASTROMON_XRAY_SRC_DTYPE))
        row["obsid"] = 8001
        row["id"] = 1
        row["detect_method"] = method
        xray = vstack([xray, row], metadata_conflicts="silent")

    cat = _cat_src_row(catalog="RFC", obsid=8001, celldetect_x_id=1)
    cat["id"] = 1
    seen = []

    def record(name, astromon_obs, astromon_xray_src, astromon_cat_src):
        seen.append(
            (name, sorted(np.asarray(astromon_xray_src["detect_method"]).astype(str)))
        )
        return astromon_cat_src[:0]

    with patch.object(rebuild_xcorr, "compute_cross_matches", record):
        rebuild_xcorr.recompute(["rfc"], db.create_table("astromon_obs"), xray, cat)

    assert seen == [("rfc", ["celldetect", "gaussian_detect"])]


def test_rebuild_verifies_before_it_writes(tmp_path):
    """A rebuild that would leave dangling references must not land at all.

    The check used to run after db.save, re-reading the file, so a failure raised
    with the bad rows already written -- the caller saw an exception and the
    database was worse than before it started.
    """
    from astromon.scripts.maintenance import rebuild_xcorr

    dbfile = tmp_path / "verify.h5"
    cat = _cat_src_row(catalog="RFC", obsid=7001, celldetect_x_id=1)
    cat["id"] = 1
    db.save("astromon_cat_src", cat, dbfile, ignore_obsid=True)
    original = _xcorr_row(select_name="astromon_21", obsid=7001, c_id=1, x_id=1)
    db.save("astromon_xcorr", original, dbfile, ignore_obsid=True)
    xray = Table(np.zeros(1, dtype=db.ASTROMON_XRAY_SRC_DTYPE))
    xray["obsid"] = 7001
    xray["id"] = 1
    xray["detect_method"] = "celldetect"
    db.save("astromon_xray_src", xray, dbfile, ignore_obsid=True)
    db.save("astromon_obs", _obs_table(7001), dbfile, ignore_obsid=True)

    def rows_that_dangle(select_names, obs, xray, cat):
        bogus = _xcorr_row(select_name="astromon_21", obsid=7001, c_id=999, x_id=1)
        return bogus[list(db.ASTROMON_XCORR_DTYPE.names)]

    with patch.object(rebuild_xcorr, "recompute", rows_that_dangle):
        with pytest.raises(RuntimeError, match="dangling"):
            rebuild_xcorr.rebuild(dbfile, obsids=[7001], force=True)

    stored = db.get_table("astromon_xcorr", dbfile)
    assert list(stored["c_id"]) == [1], "the database must be untouched"
