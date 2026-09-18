"""Tests for astromon.utils helpers that are not tied to a single pipeline stage."""

from unittest.mock import patch

import numpy as np
import pytest
from astropy.table import Table
from cxotime import CxoTime

from astromon import utils


def _fake_calalign_table():
    """A single-version calalign table, as calalign_from_files would return."""
    aca_misalign = np.tile(np.eye(3), (1, 1, 1))
    fts_misalign = np.tile(np.eye(3), (1, 1, 1))
    dy, dz = utils.get_offsets(aca_misalign)
    return Table(
        {
            "start": CxoTime(["1999:001:00:00:00"]),
            "stop": CxoTime(["2050:001:00:00:00"]),
            "detector": ["ACIS-S"],
            "caldb_version": ["4.10.0"],
            "since": CxoTime(["1999:001:00:00:00"]),
            "aca_misalign": aca_misalign,
            "fts_misalign": fts_misalign,
            "dy": dy,
            "dz": dz,
        }
    )


def test_get_calalign_offsets_row_order():
    """A shuffled input order must not silently misalign the output rows.

    get_calalign_offsets joins per-source rows against the CALALIGN table and then
    checks that the join preserved row order relative to the input. That check used
    ``np.all(a != b)``, which only fires when *every* row disagrees -- it misses a
    partial reorder where at least one row coincidentally lands back in its original
    position. With obsids given out of order ([2, 1, 3]), astropy's join() re-sorts
    by the join keys (producing [1, 2, 3]); the last position happens to match (3 ==
    3), so the old check's ``np.all(!=)`` was False and the misalignment slipped
    through silently. The fix uses ``not np.all(a == b)`` instead, which raises
    whenever any row -- not just every row -- is out of place.
    """
    all_matches = Table(
        {
            "obsid": [2, 1, 3],
            "x_id": [1, 1, 1],
            "detector": ["ACIS-S", "ACIS-S", "ACIS-S"],
            "time": CxoTime(["2015:001:00:00:00"] * 3),
            # one row uses a differently-shaped version string so this column stays
            # ragged (dtype=object) rather than collapsing into a uniform 2D array,
            # matching how get_calalign_offsets expects mixed-length caldb versions.
            "caldb_version": ["4.10.0", "4.10.0.0", "4.10.0"],
        }
    )

    with patch.object(
        utils, "calalign_from_files", return_value=_fake_calalign_table()
    ):
        with pytest.raises(RuntimeError, match="all_matches.obsid != result.obsid"):
            utils.get_calalign_offsets(all_matches)
