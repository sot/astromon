"""Test-suite guards.

The whole-catalog caches (RFC, ICRF3, Milliquas, Quaia, GaiaAGN, GaiaQSO,
MilliquasGaia) default to ``$SKA/data/astromon``, resolved at import time in
:data:`astromon.cross_match._ASTROMON_DATA_DIR`. That is a real production
directory holding hundreds of megabytes of catalogs -- and astromon.h5 itself --
so a test that calls a catalog getter without an explicit ``cache_path`` writes
there rather than into a temporary directory.

That is not hypothetical. A refactor once moved the network call in
``get_milliquas_gaia`` behind a new helper; two tests still patched the old
internals, so their patches silently became no-ops -- ``unittest.mock.patch`` on
a name nothing calls does not fail -- and the suite performed a live 130 MB CDS
X-Match and wrote a 32 MB cache into the production directory. Nothing in the
test output said so; the only signal was the suite taking three times as long.
"""

import os
from pathlib import Path

import pytest
from astropy.table import Table

from astromon import cross_match

#: Overrides the guarded directory. Exists so the guard itself can be tested
#: against a scratch directory; nothing in the package reads it.
GUARDED_DIR_ENV = "ASTROMON_TEST_GUARDED_CACHE_DIR"


def guarded_directory() -> Path:
    """The directory the guard watches, resolved fresh on every call."""
    override = os.environ.get(GUARDED_DIR_ENV)
    if override:
        return Path(override)
    return Path(cross_match._ASTROMON_DATA_DIR)


def snapshot(directory: Path) -> dict[str, tuple[int, int]]:
    """Map each entry in *directory* to its (size, mtime_ns), or {} if absent."""
    if not directory.is_dir():
        return {}
    result = {}
    for entry in os.scandir(directory):
        try:
            stat = entry.stat()
        except OSError:
            # Vanished between scandir and stat; treated as absent either way.
            continue
        result[entry.name] = (stat.st_size, stat.st_mtime_ns)
    return result


def describe_changes(
    before: dict[str, tuple[int, int]], after: dict[str, tuple[int, int]]
) -> str | None:
    """Summarise how two snapshots differ, or None when they do not."""
    added = sorted(set(after) - set(before))
    removed = sorted(set(before) - set(after))
    modified = sorted(n for n in set(before) & set(after) if before[n] != after[n])
    parts = [
        text
        for text in (
            f"added {added}" if added else "",
            f"modified {modified}" if modified else "",
            f"removed {removed}" if removed else "",
        )
        if text
    ]
    return "; ".join(parts) if parts else None


@pytest.fixture(autouse=True)
def guard_production_catalog_cache():
    """Fail any test that adds to, modifies, or removes a production cache file.

    Detects rather than prevents. Redirecting the cache paths to a temporary
    directory would be stronger, but the cache-gated tests compute flags like
    ``HAS_GAIA_AGN_CACHE`` from the real paths at module import -- before any
    fixture runs -- so they would pass their own gate and then read an empty
    directory, turning a clean skip into a download.

    A test that legitimately needs a real cache should read it and not write it.
    A test that needs to exercise a write should pass ``cache_path=tmp_path/...``
    and patch the network boundary (``cross_match.requests``), not an internal
    helper that a later refactor can route around.

    For a run that must not touch the shared caches at all, set
    ``ASTROMON_DATA_DIR`` to a scratch directory. Every cache moves before
    import, so this guard then watches the scratch directory instead. Note what
    that actually does on a *cold* scratch directory: the tests which reach a
    catalogue getter download into it, and this guard fires on each of them --
    correctly, since a test did write a cache. Seed the directory first (or run
    once and accept those reports) and the suite is then clean against it. It is
    opt-in for that reason, not because it does not work.

    Note for anyone extending this: the directory is resolved on each snapshot,
    but a test that redirects it with ``monkeypatch`` still will not be caught --
    monkeypatch is torn down before this fixture resumes. Use
    ``GUARDED_DIR_ENV`` if you need to retarget the guard itself.
    """
    directory = guarded_directory()
    before = snapshot(directory)
    yield
    changes = describe_changes(before, snapshot(guarded_directory()))
    if changes is not None:
        raise AssertionError(
            f"test touched the production catalog cache at {directory}: {changes}. "
            "Pass an explicit cache_path into a temporary directory, and patch "
            "cross_match.requests rather than an internal helper."
        )


#: For each local whole-catalog getter: the module attribute naming its cache,
#: and the reader the real getter uses on a good cache.
_LOCAL_CATALOG_CACHES = {
    "RFC": ("_RFC_CACHE_PATH", lambda path: cross_match._parse_rfc(path.read_text())),
    "ICRF3": (
        "_ICRF3_CACHE_PATH",
        lambda path: Table.read(path, format="ascii.ecsv"),
    ),
}


@pytest.fixture(autouse=True)
def local_catalogs_never_refresh(monkeypatch):
    """Read the real local catalog caches; never let a test refresh one.

    ``_get()`` routes RFC and ICRF3 to their whole-catalog getters with no
    ``cache_path``, so any test reaching them through ``rough_match`` or ``_get``
    -- ``KNOWN_VIZIER_SOURCES`` includes RFC, so the network-gated tests do --
    works against the production cache and can rewrite it.

    Suppressing the refresh with ``max_age_days=None`` is not enough: ``get_rfc``
    re-checks when the cache is stale **or** the release marker is missing, and
    the production cache has no release marker. So this bypasses the refresh
    entirely and reads the cache file, which is all these tests need.

    When a cache is absent the real getter is left in place: there is nothing to
    protect, and downloading is the intended way to populate it. Tests that
    exercise refresh logic call ``get_rfc(cache_path=...)`` directly rather than
    through ``LOCAL_CATALOG_GETTERS``, so they are unaffected either way.
    """
    for name, (path_attr, reader) in _LOCAL_CATALOG_CACHES.items():
        if name not in cross_match.LOCAL_CATALOG_GETTERS:
            continue
        cache_path = getattr(cross_match, path_attr, None)
        if cache_path is None or not Path(cache_path).exists():
            continue

        def replacement(
            pos, radius, logging_tag="", _name=name, _path=cache_path, _reader=reader
        ):
            catalog = cross_match._read_catalog_cached(Path(_path), _reader)
            return cross_match._local_catalog_near(
                catalog, _name, pos, radius, logging_tag
            )

        monkeypatch.setitem(cross_match.LOCAL_CATALOG_GETTERS, name, replacement)
