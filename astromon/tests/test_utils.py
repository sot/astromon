"""Tests for astromon.utils helpers that are not tied to a single pipeline stage."""

from unittest.mock import patch

import pytest

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
