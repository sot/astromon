"""Tests for astromon.scripts.backfill_icrf2.

Regression coverage for the bug caught in review: db.save(..., replace_keys=
("catalog",)) deletes every existing ICRF2 row globally, but the xcorr rebuild
only covers obsids backfill() actually produced a new row for. An obsid that
had an ICRF2 candidate before but gets none this time loses its cat_src row
with no matching xcorr cleanup, leaving a dangling icrf2 xcorr row.
"""

import sys
from unittest.mock import patch

import numpy as np
from astropy.table import Table

from astromon import db
from astromon.scripts import backfill_icrf2 as bf

OBSID_KEEP = 3001
OBSID_DROP = 3002


def _cat_src_row(obsid, source_id, catalog="ICRF2"):
    row = np.zeros(1, dtype=db.ASTROMON_CAT_SRC_DTYPE)
    row["obsid"] = obsid
    row["id"] = source_id
    row["catalog"] = catalog
    return Table(row)


def _xcorr_row(obsid, c_id, x_id, select_name="icrf2"):
    row = np.zeros(1, dtype=db.ASTROMON_XCORR_DTYPE)
    row["obsid"] = obsid
    row["c_id"] = c_id
    row["x_id"] = x_id
    row["select_name"] = select_name
    return Table(row)


def test_main_drops_stale_xcorr_for_an_obsid_that_lost_its_icrf2_candidate(
    tmp_path, monkeypatch
):
    dbfile = tmp_path / "astromon.h5"
    db.create_empty_tables(dbfile)

    db.save("astromon_cat_src", _cat_src_row(OBSID_KEEP, 0), dbfile)
    db.save("astromon_cat_src", _cat_src_row(OBSID_DROP, 0), dbfile)
    db.save("astromon_xcorr", _xcorr_row(OBSID_KEEP, 0, 1), dbfile)
    db.save("astromon_xcorr", _xcorr_row(OBSID_DROP, 0, 2), dbfile)

    # backfill() only produces a fresh row for OBSID_KEEP this run.
    new_cat_src = _cat_src_row(OBSID_KEEP, 0)
    new_xcorr = _xcorr_row(OBSID_KEEP, 0, 1)

    monkeypatch.setattr(sys, "argv", ["backfill_icrf2", "--dbfile", str(dbfile)])
    with patch.object(bf, "backfill", return_value=(new_cat_src, new_xcorr)):
        bf.main()

    xcorr = db.get_table("astromon_xcorr", dbfile)
    cat_src = db.get_table("astromon_cat_src", dbfile)

    assert OBSID_KEEP in np.asarray(xcorr["obsid"]), (
        "the obsid that got a real ICRF2 match must keep its xcorr row"
    )
    assert OBSID_DROP not in np.asarray(xcorr["obsid"]), (
        "the obsid that lost its ICRF2 candidate must not keep a dangling "
        "icrf2 xcorr row"
    )
    assert OBSID_DROP not in np.asarray(cat_src["obsid"])
