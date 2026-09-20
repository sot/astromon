"""Tests for astromon.scripts.maintenance.backfill_status_from_obs.

Regression coverage for the bug caught in review: build_backfill_rows()
constructed astromon_status rows with only 5 of the 7 required columns
(missing versions_done and catalog_matched), so db.save() raised on every
real (non-dry-run) invocation -- the script could never actually complete.
"""

import numpy as np
from astropy.table import Table

from astromon import db
from astromon.scripts.maintenance.backfill_status_from_obs import build_backfill_rows

OBSID = 3001


def _obs_row(obsid, ascdsver="10.8.3"):
    row = np.zeros(1, dtype=db.ASTROMON_OBS_DTYPE)
    row["obsid"] = obsid
    row["ascdsver"] = ascdsver
    return Table(row)


def _xray_row(obsid, source_id, detect_method):
    row = np.zeros(1, dtype=db.ASTROMON_XRAY_SRC_DTYPE)
    row["obsid"] = obsid
    row["id"] = source_id
    row["detect_method"] = detect_method
    return Table(row)


def test_build_backfill_rows_has_every_required_column(tmp_path):
    """The built table must satisfy db.save()'s missing-column check.

    Missing versions_done/catalog_matched made db.save("astromon_status", ...)
    raise "Saving table astromon_status with missing columns: versions_done,
    catalog_matched" on every real invocation.
    """
    dbfile = tmp_path / "astromon.h5"
    db.save("astromon_obs", _obs_row(OBSID), dbfile)
    db.save("astromon_xray_src", _xray_row(OBSID, 1, "celldetect"), dbfile)
    db.save("astromon_xcorr", db.create_table("astromon_xcorr"), dbfile)

    rows = build_backfill_rows(dbfile)

    assert set(rows.colnames) == set(db.ASTROMON_STATUS_DTYPE.names)

    db.create_empty_tables(dbfile, table_names=["astromon_status"])
    db.save("astromon_status", rows, dbfile, expect_existing=True)  # must not raise

    result = db.get_table("astromon_status", dbfile)
    assert len(result) == 1
    assert result["catalog_matched"][0] == 1


def test_build_backfill_rows_joins_completed_detect_methods(tmp_path):
    dbfile = tmp_path / "astromon.h5"
    db.save("astromon_obs", _obs_row(OBSID), dbfile)
    db.save("astromon_xray_src", _xray_row(OBSID, 1, "celldetect"), dbfile)
    db.save("astromon_xray_src", _xray_row(OBSID, 2, "gaussian_detect"), dbfile)
    db.save("astromon_xcorr", db.create_table("astromon_xcorr"), dbfile)

    rows = build_backfill_rows(dbfile)

    versions_done = str(rows["versions_done"][0])
    assert sorted(versions_done.split(",")) == ["celldetect", "gaussian_detect"]
