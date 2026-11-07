"""Tests for removing gaussian_detect rows whose celldetect seed is excluded."""

import numpy as np
from astropy.table import Table, vstack

from astromon import db
from astromon.scripts.maintenance import backfill_drop_excluded_gaussian_seeds as bf


def _xray(
    obsid,
    source_id,
    method,
    *,
    y,
    z,
    nnd=None,
    grating_arm=False,
    acis_streak=False,
    brightest=False,
    snr=0.0,
):
    row = Table(np.zeros(1, dtype=db.ASTROMON_XRAY_SRC_DTYPE))
    row["snr"] = snr
    row["obsid"] = obsid
    row["id"] = source_id
    row["detect_method"] = method
    row["y_angle"] = y
    row["z_angle"] = z
    if nnd is not None:
        row["near_neighbor_dist"] = nnd
    row["grating_arm"] = int(grating_arm)
    row["acis_streak"] = int(acis_streak)
    row["brightest"] = int(brightest)
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
    assert result["dropped_crowded"] == 1
    assert result["dropped_grating_arm"] == 0
    assert result["dropped_acis_streak"] == 0
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


def test_a_gaussian_row_seeded_from_a_grating_arm_source_is_dropped(tmp_path):
    """A seed flagged grating_arm is dropped, no exemption."""
    dbfile = _seed_db(
        tmp_path,
        [
            _xray(2462, 1, "celldetect", y=0.0, z=0.0, nnd=np.inf, grating_arm=True),
            _xray(2462, 1, "gaussian_detect", y=0.0, z=0.0, nnd=np.inf),
        ],
    )

    result = bf.backfill(dbfile)

    assert result["dropped"] == 1
    assert result["dropped_grating_arm"] == 1
    assert result["dropped_crowded"] == 0
    assert result["dropped_acis_streak"] == 0


def test_a_gaussian_row_seeded_from_an_acis_streak_source_is_dropped(tmp_path):
    dbfile = _seed_db(
        tmp_path,
        [
            _xray(15, 1, "celldetect", y=0.0, z=0.0, nnd=np.inf, acis_streak=True),
            _xray(15, 1, "gaussian_detect", y=0.0, z=0.0, nnd=np.inf),
        ],
    )

    result = bf.backfill(dbfile)

    assert result["dropped"] == 1
    assert result["dropped_acis_streak"] == 1


def test_the_brightest_source_is_exempted_from_the_streak_drop(tmp_path):
    """Mirrors _drop_acis_streak_seeds's own brightest exemption at fit time."""
    dbfile = _seed_db(
        tmp_path,
        [
            _xray(
                15,
                1,
                "celldetect",
                y=0.0,
                z=0.0,
                nnd=np.inf,
                acis_streak=True,
                brightest=True,
            ),
            _xray(15, 1, "gaussian_detect", y=0.0, z=0.0, nnd=np.inf),
        ],
    )

    result = bf.backfill(dbfile)

    assert result["dropped"] == 0


def test_brightest_does_not_exempt_crowding_or_grating_arm(tmp_path):
    """The brightest exemption is specific to the streak drop at fit time --
    it must not also suppress the other two, independent filters."""
    dbfile = _seed_db(
        tmp_path,
        [
            _xray(
                1,
                1,
                "celldetect",
                y=0.0,
                z=0.0,
                nnd=3.0,
                grating_arm=True,
                brightest=True,
            ),
            _xray(1, 1, "gaussian_detect", y=0.0, z=0.0, nnd=3.0),
        ],
    )

    result = bf.backfill(dbfile)

    assert result["dropped"] == 1
    assert result["dropped_crowded"] == 1
    assert result["dropped_grating_arm"] == 1


def test_a_row_matching_multiple_reasons_is_only_counted_once_in_dropped(tmp_path):
    """The union, not the sum, is what actually gets removed from the table."""
    dbfile = _seed_db(
        tmp_path,
        [
            _xray(
                1,
                1,
                "celldetect",
                y=0.0,
                z=0.0,
                nnd=3.0,
                grating_arm=True,
                acis_streak=True,
            ),
            _xray(1, 1, "gaussian_detect", y=0.0, z=0.0, nnd=3.0),
        ],
    )

    result = bf.backfill(dbfile)

    assert result["dropped"] == 1
    assert result["dropped_crowded"] == 1
    assert result["dropped_grating_arm"] == 1
    assert result["dropped_acis_streak"] == 1
    xray = db.get_table("astromon_xray_src", dbfile)
    assert (
        len(xray[np.asarray(xray["detect_method"]).astype(str) == "gaussian_detect"])
        == 0
    )


def _three_sources_one_on_an_arm(gaussian_ids: list[int]) -> list[Table]:
    """Seeds 1 (on a grating arm, the brightest gaussian fit), 2 and 3, 10" and 30"
    apart; gaussian rows for `gaussian_ids`, flagged as fit over all three."""
    celldetect = [
        _xray(2462, 1, "celldetect", y=0.0, z=0.0, nnd=10.0, grating_arm=True),
        _xray(2462, 2, "celldetect", y=10.0, z=0.0, nnd=10.0),
        _xray(2462, 3, "celldetect", y=40.0, z=0.0, nnd=30.0),
    ]
    gaussian = {
        1: _xray(
            2462, 1, "gaussian_detect", y=0.0, z=0.0, nnd=10.0, snr=90.0, brightest=True
        ),
        2: _xray(2462, 2, "gaussian_detect", y=10.0, z=0.0, nnd=10.0, snr=20.0),
        3: _xray(2462, 3, "gaussian_detect", y=40.0, z=0.0, nnd=30.0, snr=5.0),
    }
    return celldetect + [gaussian[source_id] for source_id in gaussian_ids]


def _gaussian_rows(dbfile) -> Table:
    xray = db.get_table("astromon_xray_src", dbfile)
    return xray[np.asarray(xray["detect_method"]).astype(str) == "gaussian_detect"]


def test_dropping_the_brightest_row_flags_the_brightest_one_left(tmp_path):
    """The brightest gaussian source was on a grating arm: a fresh fit never makes
    it, so the brightest of the rest is flagged, and it is nobody's neighbour."""
    dbfile = _seed_db(tmp_path, _three_sources_one_on_an_arm([1, 2, 3]))

    result = bf.backfill(dbfile)

    assert result["dropped"] == 1
    gaussian = _gaussian_rows(dbfile)
    assert list(gaussian["id"]) == [2, 3]
    assert list(gaussian["brightest"]) == [1, 0]
    np.testing.assert_allclose(gaussian["near_neighbor_dist"], [30.0, 30.0])


def test_rows_dropped_before_the_recompute_have_their_flags_repaired(tmp_path):
    """Rows dropped by an earlier run left the rest's flags as they were: a rerun has
    nothing to drop but still repairs them, and leaves celldetect's alone."""
    dbfile = _seed_db(tmp_path, _three_sources_one_on_an_arm([2, 3]))

    result = bf.backfill(dbfile)

    assert result["dropped"] == 0
    assert result["brightest_changed"] == 1
    assert result["near_neighbor_dist_changed"] == 1
    gaussian = _gaussian_rows(dbfile)
    assert list(gaussian["brightest"]) == [1, 0]
    np.testing.assert_allclose(gaussian["near_neighbor_dist"], [30.0, 30.0])
    xray = db.get_table("astromon_xray_src", dbfile)
    celldetect = xray[np.asarray(xray["detect_method"]).astype(str) == "celldetect"]
    np.testing.assert_allclose(celldetect["near_neighbor_dist"], [10.0, 10.0, 30.0])


def test_flags_already_right_are_not_counted_or_rewritten(tmp_path):
    dbfile = _seed_db(tmp_path, _three_sources_one_on_an_arm([2, 3]))
    bf.backfill(dbfile)
    before = dbfile.read_bytes()

    result = bf.backfill(dbfile)

    assert result["brightest_changed"] == 0
    assert result["near_neighbor_dist_changed"] == 0
    assert dbfile.read_bytes() == before
