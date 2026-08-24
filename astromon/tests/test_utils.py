"""Tests for astromon.utils helpers that are not tied to a single pipeline stage."""

from unittest.mock import patch

import numpy as np
import pytest
from astropy.table import Table
from cxotime import CxoTime

from astromon import utils


def _fake_ciao_env(prefix):
    """A CIAO environment delta as Ska.Shell.getenv would return it."""
    return {"ASCDS_INSTALL": str(prefix), "PATH": "/usr/bin:/bin"}


@pytest.fixture
def clean_ciao_env_cache():
    """Isolate the module-level CIAO_ENV cache from other tests."""
    saved = dict(utils.CIAO_ENV)
    utils.CIAO_ENV.clear()
    yield utils.CIAO_ENV
    utils.CIAO_ENV.clear()
    utils.CIAO_ENV.update(saved)


def test_ciao_env_cache_not_polluted_by_workdir(tmp_path, clean_ciao_env_cache):
    """The cached CIAO environment holds no instance-specific parameter paths.

    Ciao caches the expensive `source ciao.sh` result per prefix. Storing the live
    instance dict let the per-observation ASCDS_WORK_PATH and PFILES leak into the
    cache, so a later Ciao(prefix) with no workdir inherited a param directory that
    Observation.get_ciao had already removed.
    """
    prefix = tmp_path / "ciao"
    (prefix / "param").mkdir(parents=True)
    workdir_a = tmp_path / "obs_a" / "param"

    with patch.object(utils.Ska.Shell, "getenv", return_value=_fake_ciao_env(prefix)):
        ciao_a = utils.Ciao(prefix=prefix, workdir=workdir_a, logger="astromon")

        assert ciao_a.env["ASCDS_WORK_PATH"] == str(workdir_a)
        cached = clean_ciao_env_cache[prefix]
        assert "ASCDS_WORK_PATH" not in cached
        assert "PFILES" not in cached

        # A later instance with no workdir must not inherit obs_a's param path.
        ciao_b = utils.Ciao(prefix=prefix, logger="astromon")

    assert "ASCDS_WORK_PATH" not in ciao_b.env
    assert "PFILES" not in ciao_b.env


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
