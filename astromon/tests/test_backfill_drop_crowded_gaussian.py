"""Tests for removing gaussian_detect rows whose celldetect seed was crowded."""

import numpy as np
from astropy.table import Table, vstack

from astromon import db
from astromon.scripts.maintenance import backfill_drop_crowded_gaussian as bf


def _xray(obsid, source_id, method, *, y, z, nnd=None):
    row = Table(np.zeros(1, dtype=db.ASTROMON_XRAY_SRC_DTYPE))
    row["obsid"] = obsid
    row["id"] = source_id
    row["detect_method"] = method
    row["y_angle"] = y
    row["z_angle"] = z
    if nnd is not None:
        row["near_neighbor_dist"] = nnd
    return row


def _seed_db(tmp_path, rows):
    dbfile = tmp_path / "drop.h5"
    db.save(
        "astromon_xray_src",
        vstack(rows, metadata_conflicts="silent"),
        dbfile,
        ignore_obsid=True,
    )
    return dbfile


def test_a_gaussian_row_seeded_from_a_crowded_celldetect_source_is_dropped(tmp_path):
    """obsid 7263 in miniature: a gaussian row whose seed's own nnd is <= 6"."""
    dbfile = _seed_db(
        tmp_path,
        [
            _xray(7263, 5, "celldetect", y=-11.52, z=5.77, nnd=3.71),
            _xray(7263, 7, "celldetect", y=-15.22, z=5.55, nnd=3.71),
            _xray(7263, 5, "gaussian_detect", y=-12.56, z=5.72, nnd=8.24),
        ],
    )

    result = bf.backfill(dbfile)

    assert result["dropped"] == 1
    xray = db.get_table("astromon_xray_src", dbfile)
    assert (
        len(xray[np.asarray(xray["detect_method"]).astype(str) == "gaussian_detect"])
        == 0
    )


def test_a_gaussian_row_from_an_isolated_seed_is_kept(tmp_path):
    dbfile = _seed_db(
        tmp_path,
        [
            _xray(9000, 1, "celldetect", y=0.0, z=0.0, nnd=np.inf),
            _xray(9000, 1, "gaussian_detect", y=0.05, z=0.0, nnd=np.inf),
        ],
    )

    result = bf.backfill(dbfile)

    assert result["dropped"] == 0
    xray = db.get_table("astromon_xray_src", dbfile)
    assert len(xray) == 2


def test_celldetect_rows_are_never_dropped(tmp_path):
    """The filter targets gaussian_detect rows; celldetect rows are untouched."""
    dbfile = _seed_db(
        tmp_path,
        [
            _xray(7263, 5, "celldetect", y=-11.52, z=5.77, nnd=3.71),
            _xray(7263, 7, "celldetect", y=-15.22, z=5.55, nnd=3.71),
        ],
    )

    result = bf.backfill(dbfile)

    assert result["dropped"] == 0
    assert len(db.get_table("astromon_xray_src", dbfile)) == 2


def test_boundary_matches_the_shared_threshold(tmp_path):
    """A seed with near_neighbor_dist exactly at the threshold is crowded, the
    same as simple_cross_match's own `> NEAR_NEIGHBOR_DIST_ARCSEC` convention."""
    dbfile = _seed_db(
        tmp_path,
        [
            _xray(1, 1, "celldetect", y=0.0, z=0.0, nnd=6.0),
            _xray(1, 1, "gaussian_detect", y=0.0, z=0.0, nnd=6.0),
        ],
    )

    result = bf.backfill(dbfile)

    assert result["dropped"] == 1


def test_a_gaussian_row_with_no_matching_seed_is_left_alone(tmp_path):
    """Should not happen, but must not raise: nothing to compare against."""
    dbfile = _seed_db(
        tmp_path,
        [
            _xray(1, 1, "gaussian_detect", y=0.0, z=0.0, nnd=np.inf),
        ],
    )

    result = bf.backfill(dbfile)

    assert result["dropped"] == 0


def test_dry_run_writes_nothing(tmp_path):
    dbfile = _seed_db(
        tmp_path,
        [
            _xray(7263, 5, "celldetect", y=-11.52, z=5.77, nnd=3.71),
            _xray(7263, 7, "celldetect", y=-15.22, z=5.55, nnd=3.71),
            _xray(7263, 5, "gaussian_detect", y=-12.56, z=5.72, nnd=8.24),
        ],
    )
    before = dbfile.read_bytes()

    result = bf.backfill(dbfile, dry_run=True)

    assert result["dropped"] == 1, "the report is still real"
    assert dbfile.read_bytes() == before
