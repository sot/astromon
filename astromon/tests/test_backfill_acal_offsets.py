"""Tests for filling astromon_xray_src's applied-alignment columns from acal1 files."""

import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack
from Quaternion import Quat

from astromon import db
from astromon.scripts.maintenance import backfill_acal_offsets as bf


def _write_acal(path, obsid, caldb_version, offsets):
    """A minimal acal1 file: OBS_ID, CALDBVER and one row of applied matrices."""
    dy, dz = offsets
    matrices = {
        "aca_align": np.eye(3),
        "aca_misalign": Quat(equatorial=[dy / 3600, dz / 3600, 0]).transform,
        "fts_misalign": np.eye(3),
    }
    hdu = fits.BinTableHDU.from_columns(
        [
            fits.Column(name=name, format="9D", dim="(3,3)", array=matrix[np.newaxis])
            for name, matrix in matrices.items()
        ]
    )
    hdu.header["OBS_ID"] = str(obsid)
    hdu.header["CALDBVER"] = caldb_version
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path)


def _xray(obsid, source_id, method, caldb_version="4.9.5"):
    row = Table(np.zeros(1, dtype=db.ASTROMON_XRAY_SRC_DTYPE))
    row["obsid"] = obsid
    row["id"] = source_id
    row["detect_method"] = method
    row["caldb_version"] = caldb_version
    row["acal_dy"] = np.nan
    row["acal_dz"] = np.nan
    return row


def _seed(tmp_path):
    """Four obsids, their rows recorded with an evt2-style version 4.9.5.

    1: one acal1 file (4.9.2), in both detect methods' rows.
    2: two acal1 files sharing 4.9.4 but not the matrix (OBIs across a change).
    3: two acal1 files that disagree on CALDBVER.
    4: no acal1 file.
    5: an acal1 file from a later reprocessing (4.12.6) than its rows' (4.9.3),
       as for obsid 62649: that aspect run did not make these rows' positions.
    6: rows with no version recorded ("0.0"), so no way to tell.
    """
    acal_dir = tmp_path / "acal"
    _write_acal(acal_dir / "pcadf1N001_acal1.fits.gz", 1, "4.9.2", (84.15, 51.25))
    _write_acal(acal_dir / "pcadf2N005_acal1.fits.gz", 2, "4.9.4", (84.15, 51.25))
    _write_acal(acal_dir / "pcadf3N004_acal1.fits.gz", 2, "4.9.4", (84.40, 51.10))
    _write_acal(acal_dir / "pcadf4N006_acal1.fits.gz", 3, "4.9.4", (84.15, 51.25))
    _write_acal(acal_dir / "pcadf5N006_acal1.fits.gz", 3, "4.9.5", (84.15, 51.25))
    _write_acal(acal_dir / "pcadf6N002_acal1.fits.gz", 5, "4.12.6", (86.77, 56.60))
    _write_acal(acal_dir / "pcadf7N001_acal1.fits.gz", 6, "4.9.3", (86.68, 55.98))
    dbfile = tmp_path / "backfill.h5"
    rows = [
        _xray(1, 1, "celldetect"),
        _xray(1, 2, "celldetect"),
        _xray(1, 1, "gaussian_detect"),
        _xray(2, 1, "celldetect"),
        _xray(3, 1, "celldetect"),
        _xray(4, 1, "celldetect"),
        _xray(5, 1, "celldetect", caldb_version="4.9.3"),
        _xray(6, 1, "celldetect", caldb_version="0.0"),
    ]
    db.save(
        "astromon_xray_src",
        vstack(rows, metadata_conflicts="silent"),
        dbfile,
        ignore_obsid=True,
    )
    return dbfile, acal_dir


def _rows(dbfile, obsid):
    xray = db.get_table("astromon_xray_src", dbfile)
    return xray[np.asarray(xray["obsid"]) == obsid]


def test_an_obsid_with_one_acal_file_gets_its_matrix_and_version(tmp_path):
    dbfile, acal_dir = _seed(tmp_path)

    bf.backfill(dbfile, acal_dir)

    rows = _rows(dbfile, 1)
    assert len(rows) == 3  # every detect method's rows
    assert list(rows["caldb_version"]) == ["4.9.2"] * 3
    assert list(rows["caldb_version_source"]) == ["acal1"] * 3
    np.testing.assert_allclose(rows["acal_dy"], 84.15)
    np.testing.assert_allclose(rows["acal_dz"], 51.25)


def test_acal_files_that_disagree_on_the_matrix_record_only_the_version(tmp_path):
    dbfile, acal_dir = _seed(tmp_path)

    bf.backfill(dbfile, acal_dir)

    rows = _rows(dbfile, 2)
    assert list(rows["caldb_version"]) == ["4.9.4"]
    assert list(rows["caldb_version_source"]) == ["acal1"]
    assert np.isnan(rows["acal_dy"][0])


def test_acal_files_that_disagree_on_the_version_leave_the_rows_alone(tmp_path):
    dbfile, acal_dir = _seed(tmp_path)

    bf.backfill(dbfile, acal_dir)

    for obsid in (3, 4):  # 4 has no acal1 file at all
        rows = _rows(dbfile, obsid)
        assert list(rows["caldb_version"]) == ["4.9.5"]
        assert list(rows["caldb_version_source"]) == [""]
        assert np.isnan(rows["acal_dy"][0])


def test_an_acal_file_newer_than_the_rows_processing_is_not_used(tmp_path):
    """A later reprocessing's acal1 does not describe the stored positions.

    Obsid 62649 was stored from CALDB 4.9.3 processing; the acal1 at hand is
    from its 4.12.6 reprocessing, whose matrix differs by 0.55" in dy. Rows
    with no recorded version cannot be checked, so they are left alone too.
    """
    dbfile, acal_dir = _seed(tmp_path)

    bf.backfill(dbfile, acal_dir)

    for obsid, version in ((5, "4.9.3"), (6, "0.0")):
        rows = _rows(dbfile, obsid)
        assert list(rows["caldb_version"]) == [version]
        assert list(rows["caldb_version_source"]) == [""]
        assert np.isnan(rows["acal_dy"][0])


def test_the_summary_counts_each_outcome_and_version_change(tmp_path):
    dbfile, acal_dir = _seed(tmp_path)

    result = bf.backfill(dbfile, acal_dir)

    assert result["matrix recorded"] == 1
    assert result["acal1 files disagree on the matrix"] == 1
    assert result["acal1 newer than the recorded processing"] == 2
    assert result["acal1 files disagree on CALDBVER"] == 1
    assert result["no acal1 file"] == 1
    assert result["version changes"] == {("4.9.5", "4.9.2"): 1, ("4.9.5", "4.9.4"): 1}


def test_a_dry_run_writes_nothing(tmp_path):
    dbfile, acal_dir = _seed(tmp_path)
    before = db.get_table("astromon_xray_src", dbfile)

    result = bf.backfill(dbfile, acal_dir, dry_run=True)

    after = db.get_table("astromon_xray_src", dbfile)
    assert result["matrix recorded"] == 1
    assert list(after["caldb_version"]) == list(before["caldb_version"])
    assert np.all(np.isnan(after["acal_dy"]))


def test_recorded_versions_are_normalized(tmp_path):
    """A stored "4.9.6." (which get_calalign_offsets refuses) becomes "4.9.6"."""
    dbfile, acal_dir = _seed(tmp_path)
    xray = db.get_table("astromon_xray_src", dbfile)
    xray["caldb_version"][np.asarray(xray["obsid"]) == 4] = "4.9.6."
    db.save("astromon_xray_src", xray, dbfile, ignore_obsid=True)

    result = bf.backfill(dbfile, acal_dir)

    assert list(_rows(dbfile, 4)["caldb_version"]) == ["4.9.6"]
    assert result["versions normalized"] == 1
