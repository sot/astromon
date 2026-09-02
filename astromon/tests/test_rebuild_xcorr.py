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

from unittest.mock import patch

import numpy as np
import pytest
from astropy.table import Table, vstack

from astromon import db
from astromon.tests.test_db import _cat_src_row, _xcorr_row

# --- rebuilding xcorr from stored data --------------------------------------


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


def test_rebuild_refuses_when_cat_src_is_missing_catalogs(tmp_path):
    """A rebuild reads the stored cat_src, so an incomplete one turns the repair
    into data loss: matches whose catalog is no longer stored just vanish."""
    from astromon.scripts.maintenance.rebuild_xcorr import rebuild

    dbfile = tmp_path / "incomplete.h5"
    # cat_src holds only RFC, but an existing match came from tycho2
    cat = _cat_src_row(catalog="RFC", obsid=7001, celldetect_x_id=5)
    cat["id"] = 1
    db.save("astromon_cat_src", cat, dbfile, ignore_obsid=True)
    db.save(
        "astromon_xcorr",
        _xcorr_row(select_name="tycho2", obsid=7001, c_id=99, x_id=5),
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

    with pytest.raises(RuntimeError, match="absent from the stored cat_src"):
        rebuild(dbfile, obsids=[7001])


def test_find_orphaned_selections_ignores_selections_it_can_satisfy():
    from astromon.scripts.maintenance.rebuild_xcorr import find_orphaned_selections

    cat = _cat_src_row(catalog="Tycho2", obsid=7001, celldetect_x_id=5)
    xcorr = _xcorr_row(select_name="tycho2", obsid=7001, c_id=1, x_id=5)
    assert find_orphaned_selections(xcorr, cat) == {}


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
