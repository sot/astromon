"""Tests for filling astromon_obs_gti from local work directories."""

from pathlib import Path

import numpy as np
from astropy.table import Table

from astromon import db
from astromon.scripts.maintenance import backfill_obs_gti as bf
from astromon.tests.test_obs_gti import (
    ACIS_GTIS,
    FIRST_START,
    LAST_STOP,
    create_pre_gti_db,
    gti_row,
    write_evt2,
)

FILTERED_GTIS = [([758670000.0], [758678000.0])]


def _primary(workdir: Path, obsid: int) -> Path:
    return workdir / f"obs{obsid // 1000:02d}" / str(obsid) / "primary"


def _seed(tmp_path: Path) -> tuple[Path, Path]:
    """A pre-gti database with five obsids in astromon_obs, and their work dirs.

    1001: archive evt2, plus a filtered copy with a different span.
    1002: only the filtered copy (a cleaned-up work dir).
    1003: no work dir.
    1004: a work dir with no primary/.
    2005: an archive evt2 with no GTI extension.
    3006: a work dir with an archive evt2, but not in astromon_obs.
    """
    dbfile = tmp_path / "astromon.h5"
    workdir = tmp_path / "work"
    create_pre_gti_db(dbfile)
    obs = Table(np.zeros(5, dtype=db.ASTROMON_OBS_DTYPE))
    obs["obsid"] = [1001, 1002, 1003, 1004, 2005]
    db.save("astromon_obs", obs, dbfile)

    write_evt2(_primary(workdir, 1001) / "acisf01001N002_evt2.fits.gz", ACIS_GTIS)
    write_evt2(_primary(workdir, 1001) / "1001_evt2_filtered.fits.gz", FILTERED_GTIS)
    write_evt2(_primary(workdir, 1002) / "1002_evt2_filtered.fits.gz", FILTERED_GTIS)
    (workdir / "obs01" / "1004" / "sources").mkdir(parents=True)
    write_evt2(_primary(workdir, 2005) / "acisf02005N001_evt2.fits.gz", [])
    write_evt2(_primary(workdir, 3006) / "acisf03006N001_evt2.fits.gz", ACIS_GTIS)
    return dbfile, workdir


def _spans(dbfile: Path) -> dict[int, tuple[float, float]]:
    obs_gti = db.get_table("astromon_obs_gti", dbfile)
    return {int(row["obsid"]): (row["gti_start"], row["gti_stop"]) for row in obs_gti}


def test_backfill_fills_each_astromon_obs_obsid_from_its_archive_evt2(tmp_path):
    dbfile, workdir = _seed(tmp_path)

    summary = bf.backfill(dbfile, workdir)

    assert _spans(dbfile) == {1001: (FIRST_START, LAST_STOP)}
    assert summary["filled"] == [1001]
    assert summary["no_archive_evt2"] == [1002, 1004]
    assert summary["no_work_dir"] == [1003]
    assert list(summary["unreadable"]) == [2005]
    assert "no GTI" in summary["unreadable"][2005]


def test_backfill_dry_run_writes_nothing(tmp_path):
    dbfile, workdir = _seed(tmp_path)

    summary = bf.backfill(dbfile, workdir, dry_run=True)

    assert summary["filled"] == [1001]
    with db.connect(dbfile) as h5:
        assert "astromon_obs_gti" not in h5.root


def test_backfill_keeps_rows_for_obsids_it_does_not_fill(tmp_path):
    """A row the pipeline wrote for an obsid whose work dir is gone stays."""
    dbfile, workdir = _seed(tmp_path)
    db.save("astromon_obs_gti", gti_row(1003, 10.0, 20.0), dbfile)

    bf.backfill(dbfile, workdir)

    assert _spans(dbfile) == {1001: (FIRST_START, LAST_STOP), 1003: (10.0, 20.0)}


def test_backfill_in_parallel_matches_the_serial_result(tmp_path):
    dbfile, workdir = _seed(tmp_path)

    summary = bf.backfill(dbfile, workdir, processes=2)

    assert _spans(dbfile) == {1001: (FIRST_START, LAST_STOP)}
    assert summary["no_archive_evt2"] == [1002, 1004]
    assert list(summary["unreadable"]) == [2005]


def test_main_runs_the_backfill(tmp_path):
    dbfile, workdir = _seed(tmp_path)

    bf.main(["--db", str(dbfile), "--workdir", str(workdir)])

    assert _spans(dbfile) == {1001: (FIRST_START, LAST_STOP)}
