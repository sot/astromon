"""Tests for the per-obsid pipeline entry point in astromon.scripts.get_cat_obs_data."""

import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
from astropy.table import Table, vstack

from astromon import db


def _obs_row(obsid):
    row = db.create_table("astromon_obs")
    row.add_row(np.zeros(1, dtype=db.DTYPES["astromon_obs"])[0])
    row["obsid"] = obsid
    return row


def _obs_row_with_ascdsver(ascdsver):
    row = _obs_row(5003)
    row["ascdsver"] = ascdsver
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


def _xray_src_row(*detect_methods, obsid=5003):
    """A minimal astromon_xray_src table with one row per detect_method given."""
    row = db.create_table("astromon_xray_src")
    for method in detect_methods:
        r = Table(np.zeros(1, dtype=db.DTYPES["astromon_xray_src"]))
        r["detect_method"] = method
        r["obsid"] = obsid
        row = vstack([row, r]) if len(row) else r
    return row


def test_status_for_result_success_with_no_astromon_obs():
    """A success result missing astromon_obs/astromon_xray_src reports unknowns."""
    from astromon.scripts.get_cat_obs_data import _status_for_result

    assert _status_for_result({"ok": True, "msg": ""}) == ("success", "", "", "")


def test_status_for_result_success_reads_ascdsver_from_astromon_obs():
    from astromon.scripts.get_cat_obs_data import _status_for_result

    result = _status_for_result(
        {"ok": True, "msg": "", "astromon_obs": _obs_row_with_ascdsver("10.8.3")}
    )
    assert result == ("success", "", "10.8.3", "")


def test_status_for_result_success_reads_versions_done_from_astromon_xray_src():
    from astromon.scripts.get_cat_obs_data import _status_for_result

    result = _status_for_result(
        {
            "ok": True,
            "msg": "",
            "astromon_xray_src": _xray_src_row("celldetect", "gaussian_detect"),
        }
    )
    assert result == ("success", "", "", "celldetect,gaussian_detect")


def test_status_for_result_skipped_strips_the_prefix_and_reads_ascdsver():
    from astromon.scripts.get_cat_obs_data import _status_for_result

    result = _status_for_result(
        {
            "ok": True,
            "msg": "skipped: No x-ray sources found",
            "ascdsver": "10.8.3",
        }
    )
    assert result == ("skipped", "No x-ray sources found", "10.8.3", "")


def test_status_for_result_skipped_defaults_to_unknown_ascdsver():
    """A skip with no ascdsver key (e.g. it happened before obspar existed)."""
    from astromon.scripts.get_cat_obs_data import _status_for_result

    result = _status_for_result({"ok": True, "msg": "skipped: not public"})
    assert result == ("skipped", "not public", "", "")


def test_status_for_result_failure_strips_the_prefix_and_has_no_ascdsver():
    from astromon.scripts.get_cat_obs_data import _status_for_result

    result = _status_for_result({"ok": False, "msg": "error: CIAO crashed"})
    assert result == ("failure", "CIAO crashed", "", "")


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
                "astromon_obs": _obs_row_with_ascdsver("10.8.3"),
                "astromon_xray_src": _xray_src_row("celldetect"),
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
        assert status["ascdsver"][0] == "10.8.3"
        assert status["versions_done"][0] == "celldetect"
        # This entry point has no --skip-catalog-match: success always means
        # catalog matching was attempted.
        assert status["catalog_matched"][0] == 1


# ---- cleanup=True runs Observation.cleanup_downloads() on every outcome ----
#
# cleanup_downloads() is itself selective (keeps the filtered evt2/asol/bpix,
# *.src, PSF-size tables, and cache/ state so a later rerun of e.g.
# gaussian_detect does not need to re-download from CDA) -- process_obsid()
# is what must guarantee it actually runs regardless of how the obsid turned
# out, since its own cleanup=True previously only fired on the normal-return
# (success) path: any exception raised earlier (Skipped, SkippedWithWarning,
# ObsidNotPubliclyAvailable, or a hard failure) skipped it entirely, leaving
# the full download in place for every non-success obsid.


def _mock_observation_raising(exc):
    from astromon.scripts import get_cat_obs_data

    mock_observation = MagicMock()
    mock_observation.process.side_effect = exc
    return patch.object(
        get_cat_obs_data, "Observation", return_value=mock_observation
    ), mock_observation


def test_process_obsid_cleanup_runs_when_process_obsid_raises_skipped():
    from astromon.scripts.get_cat_obs_data import Skipped, process_obsid

    patcher, mock_observation = _mock_observation_raising(
        Skipped("no x-ray sources found")
    )
    with patcher, pytest.raises(Skipped):
        process_obsid(7001, "/tmp/whatever", cleanup=True)

    mock_observation.cleanup_downloads.assert_called_once()


def test_process_obsid_cleanup_runs_on_a_hard_failure():
    from astromon.scripts.get_cat_obs_data import process_obsid

    patcher, mock_observation = _mock_observation_raising(
        RuntimeError("fluximage failed")
    )
    with patcher, pytest.raises(RuntimeError):
        process_obsid(7001, "/tmp/whatever", cleanup=True)

    mock_observation.cleanup_downloads.assert_called_once()


def test_process_obsid_without_cleanup_flag_never_calls_cleanup_downloads():
    """--cleanup is opt-in: a failure without it must not touch the workdir."""
    from astromon.scripts.get_cat_obs_data import process_obsid

    patcher, mock_observation = _mock_observation_raising(
        RuntimeError("fluximage failed")
    )
    with patcher, pytest.raises(RuntimeError):
        process_obsid(7001, "/tmp/whatever", cleanup=False)

    mock_observation.cleanup_downloads.assert_not_called()
