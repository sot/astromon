"""Tests for the bulk rerun orchestration scripts and the CIAO environment cache."""

import sys
import tempfile
from pathlib import Path
from unittest.mock import patch

import pytest
from astropy.table import vstack

from astromon import db, utils
from astromon.scripts.maintenance import run_all
from astromon.tests.test_db import _cat_src_row, _xcorr_row

# A worker that outlives the 1s timeout the test passes, so run_one takes its timeout
# branch, but exits on its own soon after so the follow-up communicate() does not hang.
_SLEEPING_WORKER = "import time; time.sleep(3)"


def test_run_one_timeout_note_uses_configured_timeout(tmp_path, monkeypatch):
    """The timeout note and log report the enforced timeout, not the module default.

    run_one enforces the `worker_timeout` argument (wired to --timeout), but the
    message interpolated _WORKER_TIMEOUT_SEC, so a run with --timeout 7200 recorded
    "timed out after 1800s" -- understating the real limit and making a genuinely
    hung obsid look merely slow.
    """
    log_dir = tmp_path / "logs"
    log_dir.mkdir()

    real_popen = run_all.subprocess.Popen

    def fake_popen(cmd, **kwargs):
        # Ignore the real worker module; run a process that just sleeps.
        return real_popen([sys.executable, "-c", _SLEEPING_WORKER], **kwargs)

    monkeypatch.setattr(run_all.subprocess, "Popen", fake_popen)
    # Stub the process-group teardown: this test is about which timeout value ends up
    # in the note and the log, not about signal delivery.
    monkeypatch.setattr(run_all, "_kill_process_group", lambda pgid: None)

    result = run_all.run_one(
        obsid=7001,
        db_file=tmp_path / "astromon.h5",
        workdir=tmp_path / "work",
        log_dir=log_dir,
        worker_timeout=1,
    )

    assert result["status"] == "failure"
    assert "after 1s" in result["note"], result["note"]
    assert str(run_all._WORKER_TIMEOUT_SEC) not in result["note"]

    log_text = (log_dir / "7001.log").read_text()
    assert "TIMED OUT after 1s" in log_text


def _make_tree(root: Path, marker: str) -> Path:
    """Create an obs14/14321 subtree containing a file with `marker` as its content."""
    obsid_dir = root / "obs14" / "14321"
    (obsid_dir / "primary").mkdir(parents=True)
    (obsid_dir / "primary" / "evt2.fits").write_text(marker)
    return obsid_dir


def test_preserve_workdir_replaces_stale_destination(tmp_path, monkeypatch):
    """A freshly processed workdir supersedes an older copy at the destination.

    Previously the fresh tree was rmtree'd and the stale copy kept, so reprocessing
    an obsid (e.g. to pick up ssigma=3 streak detection) was silently discarded and
    a later rerun from preserve-workdir reused the old cache and images.
    """
    workdir = tmp_path / "work"
    preserve = tmp_path / "preserve"
    _make_tree(workdir, "fresh")
    _make_tree(preserve, "stale")

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "run_all",
            "--db-file",
            str(tmp_path / "astromon.h5"),
            "--workdir",
            str(workdir),
            "--obsid-list",
            str(tmp_path / "obsids.txt"),
            "--log-dir",
            str(tmp_path / "logs"),
            "--tracking-csv",
            str(tmp_path / "tracking.csv"),
            "--preserve-workdir",
            str(preserve),
        ],
    )
    (tmp_path / "obsids.txt").write_text("14321\n")

    with patch.object(
        run_all,
        "run_one",
        return_value={
            "obsid": 14321,
            "status": "success",
            "note": "",
            "returncode": 0,
            "elapsed_sec": 1.0,
            "timestamp": "2026-08-21T00:00:00",
            "log_file": "x.log",
        },
    ):
        run_all.main()

    moved = preserve / "obs14" / "14321" / "primary" / "evt2.fits"
    assert moved.read_text() == "fresh"
    assert not (workdir / "obs14" / "14321").exists()


# --- DB write path of the maintenance scripts ------------------------------
#
# These exercise save_with_lock and the backfill save helper, so they live here
# rather than in test_db.py: they import astromon.scripts.maintenance, which
# test_db.py otherwise has no dependency on. The row builders stay in test_db
# where the rest of its tests use them.


def test_save_with_lock_bootstraps_fresh_db_then_guards():
    """save_with_lock creates all tables on a fresh DB, then treats a missing one as damage.

    End-to-end check of the pipeline write path: the first obsid into a
    nonexistent file must succeed (bootstrap), but once the file exists a
    destroyed table must raise rather than being silently recreated.
    """
    from astromon.scripts.maintenance.process_one_obsid import save_with_lock

    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "pipeline.h5"

        # First obsid into a brand-new file: bootstrap, then save.
        save_with_lock(dbfile, {"astromon_cat_src": _cat_src_row(catalog="RFC")})

        # Every table exists, even the ones this obsid had no rows for
        # (get_table raises MissingTableException if one is absent).
        with db.connect(dbfile) as h5:
            present = {node.name for node in h5.root}
        assert present == set(db.DTYPES), f"missing: {set(db.DTYPES) - present}"
        assert len(db.get_table("astromon_cat_src", dbfile)) == 1
        assert len(db.get_table("astromon_xcorr", dbfile)) == 0

        # A later obsid with rows for a table that started empty still works.
        save_with_lock(
            dbfile, {"astromon_xcorr": _xcorr_row(select_name="astromon_21")}
        )
        assert len(db.get_table("astromon_xcorr", dbfile)) == 1

        # Now simulate the crash window and confirm the next write refuses.
        with db.connect(dbfile, mode="r+") as h5:
            h5.remove_node(h5.get_node("/astromon_cat_src"))

        with pytest.raises(utils.MissingTableException, match="astromon_cat_src"):
            save_with_lock(
                dbfile, {"astromon_cat_src": _cat_src_row(catalog="RFC", obsid=7009)}
            )


def test_save_with_lock_cat_src_write_drops_stale_xcorr():
    """Rewriting cat_src for an obsid clears that obsid's xcorr rows first.

    astromon_cat_src has no detect_method column, so db.save keys it on obsid alone
    and a rewrite renumbers every c_id for that obsid. xcorr rows from a
    detect_method this run did not redo would keep their old c_id and then join to
    the wrong catalog source in get_cross_matches. Other obsids are untouched.
    """
    from astromon.scripts.maintenance.process_one_obsid import save_with_lock

    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "pipeline.h5"
        db.create_empty_tables(dbfile)
        db.save(
            "astromon_xcorr",
            vstack(
                [
                    _xcorr_row(select_name="astromon_24", detect_method="celldetect"),
                    _xcorr_row(
                        select_name="astromon_25", detect_method="gaussian_detect"
                    ),
                    _xcorr_row(
                        select_name="astromon_24",
                        detect_method="celldetect",
                        obsid=7002,
                    ),
                ],
                metadata_conflicts="silent",
            ),
            dbfile,
        )

        # Rerun obsid 7001 with celldetect only.
        save_with_lock(
            dbfile,
            {
                "astromon_cat_src": _cat_src_row(catalog="RFC", obsid=7001),
                "astromon_xcorr": _xcorr_row(
                    select_name="astromon_24", detect_method="celldetect", c_id=99
                ),
            },
            obsid=7001,
        )

        xcorr = db.get_table("astromon_xcorr", dbfile)
        this_obsid = xcorr[xcorr["obsid"] == 7001]
        assert list(this_obsid["detect_method"]) == ["celldetect"], (
            "stale gaussian_detect rows must not survive a cat_src rewrite"
        )
        assert list(this_obsid["c_id"]) == [99]
        assert len(xcorr[xcorr["obsid"] == 7002]) == 1, "other obsids untouched"


def test_save_with_lock_drops_stale_xcorr_even_with_no_new_matches():
    """A rerun that produces no matches still clears the obsid's stale xcorr rows."""
    from astromon.scripts.maintenance.process_one_obsid import save_with_lock

    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "pipeline.h5"
        db.create_empty_tables(dbfile)
        db.save(
            "astromon_xcorr",
            _xcorr_row(select_name="astromon_25", detect_method="gaussian_detect"),
            dbfile,
        )

        save_with_lock(
            dbfile,
            {"astromon_cat_src": _cat_src_row(catalog="RFC", obsid=7001)},
            obsid=7001,
        )

        xcorr = db.get_table("astromon_xcorr", dbfile)
        assert len(xcorr[xcorr["obsid"] == 7001]) == 0


def test_save_with_lock_without_cat_src_keeps_xcorr():
    """A write that does not touch cat_src leaves c_id numbering, and xcorr, alone."""
    from astromon.scripts.maintenance.process_one_obsid import save_with_lock

    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "pipeline.h5"
        db.create_empty_tables(dbfile)
        db.save(
            "astromon_xcorr",
            _xcorr_row(select_name="astromon_25", detect_method="gaussian_detect"),
            dbfile,
        )

        save_with_lock(
            dbfile,
            {"astromon_obs": db.get_table("astromon_obs", dbfile)},
            obsid=7001,
        )

        assert len(db.get_table("astromon_xcorr", dbfile)) == 1
