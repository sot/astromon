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
