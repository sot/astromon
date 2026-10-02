"""Tests for the production-cache guard in conftest.

A guard that silently stops working is worse than no guard, so the wiring is
checked end-to-end in a subprocess as well as unit-tested: an in-process check
cannot retarget the guard, because monkeypatch is torn down before the guard's
teardown runs.
"""

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from astropy import units as u
from astropy.coordinates import SkyCoord

from astromon import cross_match
from astromon.tests import conftest as guard

_REPO_ROOT = Path(__file__).resolve().parents[2]


def test_snapshot_of_a_missing_directory_is_empty(tmp_path):
    assert guard.snapshot(tmp_path / "nope") == {}


def test_snapshot_records_every_entry(tmp_path):
    (tmp_path / "a.fits").write_text("aa")
    (tmp_path / "b.txt").write_text("bbbb")

    result = guard.snapshot(tmp_path)

    assert set(result) == {"a.fits", "b.txt"}
    assert result["a.fits"][0] == 2
    assert result["b.txt"][0] == 4


def test_describe_changes_is_none_when_nothing_moved(tmp_path):
    (tmp_path / "a.fits").write_text("aa")
    before = guard.snapshot(tmp_path)

    assert guard.describe_changes(before, guard.snapshot(tmp_path)) is None


def test_describe_changes_reports_an_added_file(tmp_path):
    before = guard.snapshot(tmp_path)
    (tmp_path / "new_catalog.fits").write_text("x")

    assert guard.describe_changes(before, guard.snapshot(tmp_path)) == (
        "added ['new_catalog.fits']"
    )


def test_describe_changes_reports_a_modified_file(tmp_path):
    victim = tmp_path / "cat.fits"
    victim.write_text("original")
    before = guard.snapshot(tmp_path)
    victim.write_text("a longer body, so the size differs")

    assert guard.describe_changes(before, guard.snapshot(tmp_path)) == (
        "modified ['cat.fits']"
    )


def test_describe_changes_reports_a_removed_file(tmp_path):
    victim = tmp_path / "cat.fits"
    victim.write_text("x")
    before = guard.snapshot(tmp_path)
    victim.unlink()

    assert guard.describe_changes(before, guard.snapshot(tmp_path)) == (
        "removed ['cat.fits']"
    )


def _run_pytest_with_guard(tmp_path, body, guarded_dir):
    """Run pytest on a one-test file in tmp_path, with the guard retargeted."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    shutil.copy(Path(guard.__file__), tmp_path / "conftest.py")
    (tmp_path / "test_subject.py").write_text(body)
    return subprocess.run(
        [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", str(tmp_path)],
        capture_output=True,
        text=True,
        check=False,
        cwd=tmp_path,
        # The worktree, not the installed astromon: cwd is a scratch directory,
        # so without this the subprocess imports whichever version is on the env.
        env={
            **os.environ,
            guard.GUARDED_DIR_ENV: str(guarded_dir),
            "PYTHONPATH": str(_REPO_ROOT),
        },
    )


def test_guard_fails_a_test_that_writes_to_the_guarded_directory(tmp_path):
    guarded = tmp_path / "guarded"
    guarded.mkdir()
    run = _run_pytest_with_guard(
        tmp_path / "run",
        "from pathlib import Path\n"
        "def test_writes():\n"
        f"    Path({str(guarded)!r}).joinpath('leaked.fits').write_text('x')\n",
        guarded,
    )

    assert run.returncode != 0, run.stdout
    assert "touched the production catalog cache" in run.stdout
    assert "added ['leaked.fits']" in run.stdout


def test_guard_passes_a_test_that_writes_nowhere(tmp_path):
    guarded = tmp_path / "guarded"
    guarded.mkdir()
    run = _run_pytest_with_guard(
        tmp_path / "run", "def test_quiet():\n    assert True\n", guarded
    )

    assert run.returncode == 0, run.stdout


def test_local_catalog_dispatch_does_not_call_the_refreshing_getter():
    """RFC dispatch must not go through get_rfc, whose refresh can write."""

    def explode(*args, **kwargs):
        raise AssertionError("get_rfc must not be called from a dispatched lookup")

    original = cross_match.get_rfc
    cross_match.get_rfc = explode
    try:
        position = SkyCoord([187.2779], [2.0524], unit="deg")
        cross_match.LOCAL_CATALOG_GETTERS["RFC"](
            pos=position, radius=5 * u.arcsec, logging_tag=""
        )
    finally:
        cross_match.get_rfc = original


def test_local_catalog_dispatch_still_finds_a_real_source():
    """Bypassing the refresh must not break the lookup: 3C 273 is in RFC."""
    position = SkyCoord([187.2779], [2.0524], unit="deg")

    result = cross_match.LOCAL_CATALOG_GETTERS["RFC"](
        pos=position, radius=5 * u.arcsec, logging_tag=""
    )

    assert len(result) >= 1


def test_max_age_days_none_is_not_enough_on_its_own(tmp_path):
    """Why the fixture bypasses the refresh instead of just passing None.

    get_rfc re-checks upstream when the cache is stale *or* the release marker is
    missing. Production has no release marker, so max_age_days=None leaves the
    check due and a test can still reach the network and write.
    """
    cache = tmp_path / "rfc_catalog.txt"
    cache.write_text("body")
    old = time.time() - 30 * 86400
    os.utime(cache, (old, old))

    assert not cross_match._cache_is_stale(cache, None), (
        "max_age_days=None must make staleness alone stop mattering"
    )
    # Named by convention rather than through the helper, which is
    # _rfc_release_marker on some branches and _release_marker on others.
    release_marker = cache.with_name(cache.name + ".release")
    assert not release_marker.exists(), (
        "the missing release marker is the other half of the condition"
    )
