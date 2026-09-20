"""Tests for astromon.scripts.maintenance.backfill_icrf3.

Regression coverage for the bug caught in review: db.save(..., replace_keys=
("catalog",)) deletes every existing ICRS row globally, but the xcorr rebuild
only covers obsids that got a NEW ICRS match this run. An obsid that had an
ICRS candidate before but gets none this time (its celldetect sources are
gone, or ICRF3 no longer has an in-FOV source) loses its cat_src row with no
matching xcorr cleanup, leaving a dangling astromon_21 xcorr row.
"""

import sys
from unittest.mock import patch

import numpy as np
from astropy.table import Table

from astromon import cross_match, db
from astromon.scripts.maintenance import backfill_icrf3 as bf

OBSID_KEEP = 1001
OBSID_DROP = 1002


def _obs_row(obsid, ra=150.0, dec=2.0):
    row = np.zeros(1, dtype=db.ASTROMON_OBS_DTYPE)
    row["obsid"] = obsid
    row["ra"] = ra
    row["dec"] = dec
    return Table(row)


def _xray_row(obsid, source_id, ra, dec):
    row = np.zeros(1, dtype=db.ASTROMON_XRAY_SRC_DTYPE)
    row["obsid"] = obsid
    row["id"] = source_id
    row["ra"] = ra
    row["dec"] = dec
    row["detect_method"] = "celldetect"
    return Table(row)


def _icrs_cat_row(obsid, source_id):
    row = np.zeros(1, dtype=db.ASTROMON_CAT_SRC_DTYPE)
    row["obsid"] = obsid
    row["id"] = source_id
    row["catalog"] = "ICRS"
    return Table(row)


def _xcorr_row(obsid, c_id, x_id, select_name="astromon_21"):
    row = np.zeros(1, dtype=db.ASTROMON_XCORR_DTYPE)
    row["obsid"] = obsid
    row["c_id"] = c_id
    row["x_id"] = x_id
    row["select_name"] = select_name
    row["detect_method"] = "celldetect"
    return Table(row)


def test_main_drops_stale_xcorr_for_an_obsid_that_lost_its_icrs_candidate(
    tmp_path, monkeypatch
):
    dbfile = tmp_path / "astromon.h5"
    db.create_empty_tables(dbfile)

    ra, dec = 150.0, 2.0
    db.save("astromon_obs", _obs_row(OBSID_KEEP, ra, dec), dbfile)
    db.save("astromon_xray_src", _xray_row(OBSID_KEEP, 3, ra, dec), dbfile)
    db.save("astromon_cat_src", _icrs_cat_row(OBSID_KEEP, 0), dbfile)
    db.save("astromon_cat_src", _icrs_cat_row(OBSID_DROP, 0), dbfile)
    db.save("astromon_xcorr", _xcorr_row(OBSID_KEEP, 0, 3), dbfile)
    db.save("astromon_xcorr", _xcorr_row(OBSID_DROP, 0, 9), dbfile)
    # OBSID_DROP has no astromon_xray_src row at all (e.g. its celldetect
    # sources were dropped by an earlier backfill), so build_icrs_cat_src can
    # never produce a new ICRS row for it.

    icrf3 = Table({"name": ["ICRF3-test"], "ra": [ra + 1.0 / 3600], "dec": [dec]})

    monkeypatch.setattr(sys, "argv", ["backfill_icrf3", "--db", str(dbfile)])
    with patch.object(cross_match, "get_icrf3", return_value=icrf3):
        bf.main()

    xcorr = db.get_table("astromon_xcorr", dbfile)
    cat_src = db.get_table("astromon_cat_src", dbfile)

    assert OBSID_KEEP in np.asarray(xcorr["obsid"]), (
        "the obsid that got a real ICRS match must keep its xcorr row"
    )
    assert OBSID_DROP not in np.asarray(xcorr["obsid"]), (
        "the obsid that lost its ICRS candidate must not keep a dangling "
        "astromon_21 xcorr row"
    )
    assert OBSID_DROP not in np.asarray(cat_src["obsid"])
