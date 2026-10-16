"""Tests for recomputing astromon_cat_src.y_angle/z_angle from evt2 pointing.

Context: requery_cat_src.py and backfill_icrf2.py used to derive their pointing
quaternion from mica's local obspar archive / astromon_obs's nominal target
values -- both documented elsewhere as unreliable, disagreeing with the evt2
header's measured pointing by tens of arcsec in real cases. That has since been
fixed to read the evt2 header instead (see evt2_pointing), matching
astromon_xray_src's own y_angle/z_angle. This module repairs the rows already
written with the wrong pointing: y_angle/z_angle recomputed in place from each
row's own (already-correct) ra/dec and the now-correct evt2 quaternion, leaving
id/obsid/catalog/ra/dec/name/mag/separation/celldetect_x_id untouched so
astromon_xcorr's c_id references are unaffected by this write.
"""

from unittest.mock import patch

import numpy as np
import pytest
from astropy.table import Table
from cxotime import CxoTime
from Quaternion import Quat

from astromon import db
from astromon.scripts.maintenance import fix_cat_src_pointing as fix_pointing


def _cat_src_row(obsid, source_id, catalog, ra, dec, y_angle, z_angle):
    row = Table(np.zeros(1, dtype=db.ASTROMON_CAT_SRC_DTYPE))
    row["obsid"] = obsid
    row["id"] = source_id
    row["catalog"] = catalog
    row["ra"] = ra
    row["dec"] = dec
    row["y_angle"] = y_angle
    row["z_angle"] = z_angle
    return row


def _fake_pointing(ra, dec):
    """A quaternion pointed exactly at (ra, dec): y_angle/z_angle both 0 there."""
    return Quat(equatorial=(ra, dec, 0.0)), CxoTime("2020-01-01T00:00:00")


def test_recompute_updates_y_angle_z_angle_from_the_evt2_pointing():
    cat_src = Table(
        np.concatenate(
            [
                _cat_src_row(
                    7001, 0, "RFC", ra=150.0, dec=2.0, y_angle=999, z_angle=999
                ),
            ]
        )
    )

    with patch.object(
        fix_pointing, "evt2_pointing", return_value=_fake_pointing(150.0, 2.0)
    ):
        updated, skipped, n_updated = fix_pointing.recompute_cat_src_pointing(
            cat_src, [7001], workdir="/wd", archive_dir="/ad"
        )

    assert n_updated == 1
    assert skipped == []
    assert updated["y_angle"][0] == pytest.approx(0.0, abs=1e-6)
    assert updated["z_angle"][0] == pytest.approx(0.0, abs=1e-6)


def test_recompute_leaves_id_obsid_catalog_ra_dec_untouched():
    """The whole point: only y_angle/z_angle change, so c_id references survive."""
    cat_src = Table(
        np.concatenate(
            [
                _cat_src_row(
                    7001, 5, "RFC", ra=150.0, dec=2.0, y_angle=999, z_angle=999
                ),
            ]
        )
    )

    with patch.object(
        fix_pointing, "evt2_pointing", return_value=_fake_pointing(150.0, 2.0)
    ):
        updated, _, _ = fix_pointing.recompute_cat_src_pointing(
            cat_src, [7001], workdir="/wd", archive_dir="/ad"
        )

    assert updated["id"][0] == 5
    assert updated["obsid"][0] == 7001
    assert str(updated["catalog"][0]) == "RFC"
    assert updated["ra"][0] == pytest.approx(150.0)
    assert updated["dec"][0] == pytest.approx(2.0)


def test_recompute_covers_every_catalog_for_a_touched_obsid():
    """Every row for a touched obsid is recomputed, not just the catalogs the
    bug happened to touch -- the same evt2 quaternion applies to all of them,
    and this keeps the obsid's rows uniformly consistent."""
    cat_src = Table(
        np.concatenate(
            [
                _cat_src_row(
                    7001, 0, "RFC", ra=150.0, dec=2.0, y_angle=999, z_angle=999
                ),
                _cat_src_row(
                    7001, 1, "Tycho2", ra=150.001, dec=2.001, y_angle=999, z_angle=999
                ),
            ]
        )
    )

    with patch.object(
        fix_pointing, "evt2_pointing", return_value=_fake_pointing(150.0, 2.0)
    ):
        updated, _, n_updated = fix_pointing.recompute_cat_src_pointing(
            cat_src, [7001], workdir="/wd", archive_dir="/ad"
        )

    assert n_updated == 2
    assert not np.any(updated["y_angle"] == 999)


def test_recompute_skips_an_obsid_with_no_local_evt2_file():
    """No evt2 file is a skip, not a fall back to a worse pointing source --
    and the row's existing (wrong) values are left alone rather than guessed at."""
    cat_src = Table(
        np.concatenate(
            [
                _cat_src_row(
                    7001, 0, "RFC", ra=150.0, dec=2.0, y_angle=999, z_angle=999
                ),
            ]
        )
    )

    with patch.object(fix_pointing, "evt2_pointing", return_value=None):
        updated, skipped, n_updated = fix_pointing.recompute_cat_src_pointing(
            cat_src, [7001], workdir="/wd", archive_dir="/ad"
        )

    assert skipped == [7001]
    assert n_updated == 0
    assert updated["y_angle"][0] == 999


def test_recompute_only_returns_rows_for_the_given_obsids():
    """The subset returned is meant to go straight to db.save's default
    merge-by-obsid-key write -- it must not include other obsids' rows, or a
    caller could accidentally believe untouched obsids need rewriting too."""
    cat_src = Table(
        np.concatenate(
            [
                _cat_src_row(
                    7001, 0, "RFC", ra=150.0, dec=2.0, y_angle=999, z_angle=999
                ),
                _cat_src_row(8001, 0, "RFC", ra=10.0, dec=-5.0, y_angle=1, z_angle=1),
            ]
        )
    )

    with patch.object(
        fix_pointing, "evt2_pointing", return_value=_fake_pointing(150.0, 2.0)
    ):
        updated, _, _ = fix_pointing.recompute_cat_src_pointing(
            cat_src, [7001], workdir="/wd", archive_dir="/ad"
        )

    assert list(updated["obsid"]) == [7001]
