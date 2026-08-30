"""Tests for the per-obsid pipeline entry point in astromon.scripts.get_cat_obs_data."""

import tempfile
from pathlib import Path

import numpy as np
from astropy.table import Table

from astromon import db


def _obs_row(obsid):
    row = db.create_table("astromon_obs")
    row.add_row(np.zeros(1, dtype=db.DTYPES["astromon_obs"])[0])
    row["obsid"] = obsid
    return row


def _empty(name):
    return Table(dtype=db.DTYPES[name])


def test_split_versions_separates_peak_gaussian():
    """peak_gaussian_detect is run for its .src output, not saved as a method.

    Its fitted position never won on match count; what it carries is the
    peak_offset diagnostic, which get_sources reads off the .src file. So it is
    run but not persisted.
    """
    from astromon.scripts.get_cat_obs_data import _split_versions

    saved, diagnostic = _split_versions(
        ("celldetect", "gaussian_detect", "peak_gaussian_detect")
    )
    assert saved == ["celldetect", "gaussian_detect"]
    assert diagnostic == ["peak_gaussian_detect"]


def test_split_versions_preserves_order_and_handles_absence():
    from astromon.scripts.get_cat_obs_data import _split_versions

    assert _split_versions(("gaussian_detect", "celldetect")) == (
        ["gaussian_detect", "celldetect"],
        [],
    )
    assert _split_versions(("peak_gaussian_detect",)) == ([], ["peak_gaussian_detect"])
    assert _split_versions(()) == ([], [])


def test_status_for_result_success():
    from astromon.scripts.get_cat_obs_data import _status_for_result

    assert _status_for_result({"ok": True, "msg": ""}) == ("success", "")


def test_status_for_result_skipped_strips_the_prefix():
    from astromon.scripts.get_cat_obs_data import _status_for_result

    result = _status_for_result({"ok": True, "msg": "skipped: No x-ray sources found"})
    assert result == ("skipped", "No x-ray sources found")


def test_status_for_result_failure_strips_the_prefix():
    from astromon.scripts.get_cat_obs_data import _status_for_result

    result = _status_for_result({"ok": False, "msg": "error: CIAO crashed"})
    assert result == ("failure", "CIAO crashed")


def test_save_records_status_for_skipped_and_failed_obsids():
    """Obsids with nothing to write to the other tables still get a status row.

    This is the actual gap this whole feature closes: process()'s "skipped" and
    "failure" results carry astromon_obs=[] and were previously dropped entirely
    before save() ever looked at them, so a query against the DB alone could not
    tell "never attempted" apart from "attempted, skipped".
    """
    from astromon.scripts.get_cat_obs_data import save

    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "pipeline.h5"
        db.create_empty_tables(dbfile)

        results = [
            {
                "ok": True,
                "msg": "skipped: No x-ray sources found",
                "obsid": 5001,
                "astromon_obs": _empty("astromon_obs"),
                "astromon_xray_src": _empty("astromon_xray_src"),
                "astromon_cat_src": _empty("astromon_cat_src"),
                "astromon_xcorr": _empty("astromon_xcorr"),
            },
            {
                "ok": False,
                "msg": "error: CIAO crashed",
                "obsid": 5002,
                "astromon_obs": _empty("astromon_obs"),
                "astromon_xray_src": _empty("astromon_xray_src"),
                "astromon_cat_src": _empty("astromon_cat_src"),
                "astromon_xcorr": _empty("astromon_xcorr"),
            },
        ]

        save(results, dbfile)

        status = db.get_table("astromon_status", dbfile)
        by_obsid = {row["obsid"]: row for row in status}
        assert len(status) == 2
        assert by_obsid[5001]["status"] == "skipped"
        assert by_obsid[5001]["note"] == "No x-ray sources found"
        assert by_obsid[5002]["status"] == "failure"
        assert by_obsid[5002]["note"] == "CIAO crashed"
        assert len(db.get_table("astromon_obs", dbfile)) == 0


def test_save_records_status_alongside_a_successful_write():
    from astromon.scripts.get_cat_obs_data import save

    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "pipeline.h5"
        db.create_empty_tables(dbfile)

        results = [
            {
                "ok": True,
                "msg": "",
                "obsid": 5003,
                "astromon_obs": _obs_row(5003),
                "astromon_xray_src": _empty("astromon_xray_src"),
                "astromon_cat_src": _empty("astromon_cat_src"),
                "astromon_xcorr": _empty("astromon_xcorr"),
            }
        ]

        save(results, dbfile)

        assert len(db.get_table("astromon_obs", dbfile)) == 1
        status = db.get_table("astromon_status", dbfile)
        assert len(status) == 1
        assert status["obsid"][0] == 5003
        assert status["status"][0] == "success"
