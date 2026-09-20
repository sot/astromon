"""Tests for astromon.scripts.maintenance.backfill_rfc.

Regression coverage for the bug caught in review: merged_cat drops every
existing ICRS/RFC row globally (non_rfc_cat excludes them all), but the xcorr
rebuild only covers obsids that got a NEW RFC match this run. An obsid that
had an ICRS/RFC candidate before but gets none this time loses its cat_src
row with no matching xcorr cleanup, leaving dangling xcorr rows for every
select_name that includes RFC.
"""

import sys
from unittest.mock import patch

import numpy as np
from astropy.table import Table

from astromon import db
from astromon.scripts.maintenance import backfill_rfc as bf

OBSID_KEEP = 2001
OBSID_DROP = 2002


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


def _rfc_cat_row(obsid, source_id, catalog="RFC"):
    row = np.zeros(1, dtype=db.ASTROMON_CAT_SRC_DTYPE)
    row["obsid"] = obsid
    row["id"] = source_id
    row["catalog"] = catalog
    return Table(row)


def _xcorr_row(obsid, c_id, x_id, select_name):
    row = np.zeros(1, dtype=db.ASTROMON_XCORR_DTYPE)
    row["obsid"] = obsid
    row["c_id"] = c_id
    row["x_id"] = x_id
    row["select_name"] = select_name
    row["detect_method"] = "celldetect"
    return Table(row)


def test_main_drops_stale_xcorr_for_an_obsid_that_lost_its_rfc_candidate(
    tmp_path, monkeypatch
):
    dbfile = tmp_path / "astromon.h5"
    db.create_empty_tables(dbfile)

    ra, dec = 150.0, 2.0
    db.save("astromon_obs", _obs_row(OBSID_KEEP, ra, dec), dbfile)
    db.save("astromon_xray_src", _xray_row(OBSID_KEEP, 3, ra, dec), dbfile)
    db.save("astromon_cat_src", _rfc_cat_row(OBSID_KEEP, 0), dbfile)
    db.save("astromon_cat_src", _rfc_cat_row(OBSID_DROP, 0), dbfile)
    db.save("astromon_xcorr", _xcorr_row(OBSID_KEEP, 0, 3, "astromon_21"), dbfile)
    db.save("astromon_xcorr", _xcorr_row(OBSID_DROP, 0, 9, "astromon_21"), dbfile)
    # OBSID_DROP has no astromon_xray_src row at all, so build_rfc_cat_src can
    # never produce a new RFC row for it.

    rfc = Table({"name": ["RFC-test"], "ra": [ra + 1.0 / 3600], "dec": [dec]})

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "backfill_rfc",
            "--db",
            str(dbfile),
            "--catalog-cache",
            str(tmp_path / "rfc_cache.txt"),
        ],
    )
    with patch.object(bf, "get_rfc", return_value=rfc):
        bf.main()

    xcorr = db.get_table("astromon_xcorr", dbfile)
    cat_src = db.get_table("astromon_cat_src", dbfile)

    assert OBSID_KEEP in np.asarray(xcorr["obsid"]), (
        "the obsid that got a real RFC match must keep its xcorr row"
    )
    assert OBSID_DROP not in np.asarray(xcorr["obsid"]), (
        "the obsid that lost its RFC candidate must not keep a dangling xcorr row"
    )
    assert OBSID_DROP not in np.asarray(cat_src["obsid"])
