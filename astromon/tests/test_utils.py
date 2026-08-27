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


def test_ciao_env_cache_avoids_repeat_getenv_call(tmp_path, clean_ciao_env_cache):
    """A cached prefix must not re-invoke the expensive `source ciao.sh` call.

    CIAO_ENV.get(prefix, Ska.Shell.getenv(...)) evaluates the default argument
    eagerly, so the subprocess call ran on every Ciao() construction regardless
    of whether prefix was already cached.
    """
    prefix = tmp_path / "ciao"
    (prefix / "param").mkdir(parents=True)

    with patch.object(
        utils.Ska.Shell, "getenv", return_value=_fake_ciao_env(prefix)
    ) as mock_getenv:
        utils.Ciao(prefix=prefix, logger="astromon")
        assert mock_getenv.call_count == 1

        utils.Ciao(prefix=prefix, logger="astromon")
        assert mock_getenv.call_count == 1, (
            "a cached prefix must not re-invoke Ska.Shell.getenv"
        )


def _fake_calalign_table():
    """A minimal two-version calalign table, as calalign_from_files would return."""
    aca_misalign = np.tile(np.eye(3), (2, 1, 1))
    fts_misalign = np.tile(np.eye(3), (2, 1, 1))
    dy, dz = utils.get_offsets(aca_misalign)
    return Table(
        {
            "start": CxoTime(["1999:001:00:00:00", "2010:001:00:00:00"]),
            "stop": CxoTime(["2010:001:00:00:00", "2050:001:00:00:00"]),
            "detector": ["ACIS-S", "ACIS-S"],
            "caldb_version": ["4.0.0", "4.10.0"],
            "since": CxoTime(["1999:001:00:00:00", "2010:001:00:00:00"]),
            "aca_misalign": aca_misalign,
            "fts_misalign": fts_misalign,
            "dy": dy,
            "dz": dz,
        }
    )


def test_get_calalign_offsets_row_order():
    """A shuffled input row order must come back out in that same order.

    join() sorts by the join keys internally and gives no guarantee that the
    input row order survives. get_calalign_offsets used to raise whenever the
    join's output order didn't already match all_matches' input order -- an
    unnecessary restriction on a valid, arbitrarily-ordered input table, since
    match_id_keys (obsid, x_id[, detect_method]) already uniquely identify
    each row. It must instead restore all_matches' own order rather than
    reject it.
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
        result = utils.get_calalign_offsets(all_matches)

    assert list(result["obsid"]) == [2, 1, 3]


def test_get_calalign_offsets_raises_on_malformed_version_string():
    """A non-numeric version segment must raise rather than being silently dropped.

    get_calalign_offsets used to build the caldb_version/calalign_version
    comparison tuples with ``if f.isdigit()``, which drops any non-numeric
    dot-separated segment instead of raising. A truncated tuple then compares
    incorrectly (lexicographically shorter-vs-longer) against a full-length
    one, silently picking the wrong CalDB row as "actual" or "reference"
    instead of failing loudly on the unexpected input.
    """
    all_matches = Table(
        {
            "obsid": [1],
            "x_id": [1],
            "detector": ["ACIS-S"],
            "time": CxoTime(["2015:001:00:00:00"]),
            "caldb_version": ["4.N0.0"],
        }
    )

    with patch.object(
        utils, "calalign_from_files", return_value=_fake_calalign_table()
    ):
        with pytest.raises(ValueError):
            utils.get_calalign_offsets(all_matches)


def test_get_calalign_offsets_disambiguates_x_id_by_detect_method():
    """x_id numbering restarts per detect_method, so obsid+x_id alone can collide.

    celldetect and gaussian_detect each number their sources from 1 for a given
    obsid, so a matches table spanning both methods can have two physically
    different sources sharing (obsid, x_id). Without 'detect_method' in the
    grouping key, get_calalign_offsets collapses/mismatches those rows and used
    to raise RuntimeError("len(all_matches) != len(actual)").
    """
    all_matches = Table(
        {
            "obsid": [1, 1],
            "x_id": [1, 1],
            "detect_method": ["celldetect", "gaussian_detect"],
            "detector": ["ACIS-S", "ACIS-S"],
            "time": CxoTime(["2015:001:00:00:00", "2015:001:00:00:00"]),
            "caldb_version": ["4.10.0", "4.10.0"],
        }
    )

    with patch.object(
        utils, "calalign_from_files", return_value=_fake_calalign_table()
    ):
        result = utils.get_calalign_offsets(all_matches)

    assert len(result) == 2
    assert list(result["detect_method"]) == ["celldetect", "gaussian_detect"]


def test_get_calalign_offsets_without_detect_method_column():
    """Tables without a 'detect_method' column (e.g. from older DBs) still work."""
    all_matches = Table(
        {
            "obsid": [1, 2],
            "x_id": [1, 1],
            "detector": ["ACIS-S", "ACIS-S"],
            "time": CxoTime(["2015:001:00:00:00", "2015:001:00:00:00"]),
            "caldb_version": ["4.10.0", "4.10.0"],
        }
    )

    with patch.object(
        utils, "calalign_from_files", return_value=_fake_calalign_table()
    ):
        result = utils.get_calalign_offsets(all_matches)

    assert list(result["obsid"]) == [1, 2]
    assert "detect_method" not in result.colnames
