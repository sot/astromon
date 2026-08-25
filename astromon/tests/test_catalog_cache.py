"""Tests for the locally cached catalogs' refresh policy.

get_rfc checks astrogeo.org for a newer release at most once per `max_age_days`,
so that a pipeline run does not re-check once per obsid.  These tests cover the
location, refresh, and fallback behavior that would otherwise be spread across
the catalogue-query tests. Keeping them together makes it easier to review the
cache contract without paging through matching logic.
"""

import inspect
import json
import os
import subprocess
import sys
import time
import urllib.error
from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np
import pytest
from astropy.table import Table

from astromon import cross_match

# Real records from the RFC fixed-width catalog, enough for _parse_rfc to return
# a usable table so these tests exercise the real read path.
RFC_SAMPLE = """\
# Radio Fundamental Catalogue
# The Radio Fundamental Catalogue, data release: rfc_2026b
#
RFC J0000-0022  2357-006  00 00 01.659837 -00 22 10.02335    1.26   2.68    0.230
RFC J1229+0203  3C273B    12 29 06.699755 +02 03 08.59833    0.09   0.17   -0.013
"""


def _resolved_cache_paths(env_overrides):
    """Import cross_match in a subprocess and report where its caches resolve."""
    script = (
        "import json\n"
        "from astromon import cross_match\n"
        "from astromon import observation, utils\n"
        "print(json.dumps({\n"
        "    'dir': str(utils.ASTROMON_DATA_DIR),\n"
        "    'archive': str(observation.ARCHIVE_DIR),\n"
        "    'rfc': str(cross_match._RFC_CACHE_PATH),\n"
        "    'quaia': str(cross_match._QUAIA_CACHE_PATH),\n"
        "    'all': {n: str(getattr(cross_match, n))\n"
        "            for n in dir(cross_match) if n.endswith('_CACHE_PATH')},\n"
        "}))\n"
    )
    repo_root = Path(__file__).resolve().parents[2]
    child_env = {k: v for k, v in os.environ.items() if k != "ASTROMON_DATA_DIR"}
    child_env["PYTHONPATH"] = str(repo_root)
    child_env.update(env_overrides)
    completed = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        check=True,
        env=child_env,
    )
    return json.loads(completed.stdout)


def _stale_cache(tmp_path, age_days=10):
    """An RFC cache old enough that the release check is due."""
    cache = tmp_path / "rfc_catalog.txt"
    cache.write_text(RFC_SAMPLE)
    old = time.time() - age_days * 86400
    os.utime(cache, (old, old))
    return cache


def _failing_urlretrieve(calls):
    def urlretrieve(url, filename):
        calls.append((url, filename))
        raise urllib.error.HTTPError(url, 404, "Not Found", None, None)

    return urlretrieve


def test_cache_is_stale_absent_file(tmp_path):
    """A cache_path that does not exist is always stale."""
    assert cross_match._cache_is_stale(tmp_path / "missing.txt", max_age_days=30)


def test_cache_is_stale_respects_max_age(tmp_path):
    """An existing file is stale only once it exceeds max_age_days."""
    cache_path = tmp_path / "cached.txt"
    cache_path.write_text("data")

    assert not cross_match._cache_is_stale(cache_path, max_age_days=30)

    old_time = time.time() - 31 * 86400
    os.utime(cache_path, (old_time, old_time))
    assert cross_match._cache_is_stale(cache_path, max_age_days=30)


def test_cache_is_stale_none_max_age_never_refreshes(tmp_path):
    """max_age_days=None means only absence makes the cache stale."""
    cache_path = tmp_path / "cached.txt"
    cache_path.write_text("data")

    old_time = time.time() - 10_000 * 86400
    os.utime(cache_path, (old_time, old_time))
    assert not cross_match._cache_is_stale(cache_path, max_age_days=None)


def test_cache_dir_defaults_under_ska():
    """Without an override the caches stay where they have always been."""
    resolved = _resolved_cache_paths({"SKA": "/tmp/fake-ska"})

    assert resolved["dir"] == "/tmp/fake-ska/data/astromon"
    assert resolved["rfc"] == "/tmp/fake-ska/data/astromon/rfc_catalog.txt"


def test_astromon_data_dir_overrides_the_cache_location():
    """ASTROMON_DATA_DIR isolates the catalog caches the way ASTROMON_FILE does the DB."""
    resolved = _resolved_cache_paths(
        {"SKA": "/tmp/fake-ska", "ASTROMON_DATA_DIR": "/tmp/isolated-catalogs"}
    )

    assert resolved["dir"] == "/tmp/isolated-catalogs"
    assert resolved["rfc"] == "/tmp/isolated-catalogs/rfc_catalog.txt"
    assert resolved["quaia"] == "/tmp/isolated-catalogs/quaia_catalog.fits"


def test_archive_dir_defaults_under_the_shared_data_dir():
    """The observation archive lives under the same root as the catalogs."""
    resolved = _resolved_cache_paths({"SKA": "/tmp/fake-ska"})

    assert resolved["archive"] == "/tmp/fake-ska/data/astromon/xray_observations"


def test_astromon_data_dir_moves_the_archive_too():
    """ARCHIVE_DIR follows ASTROMON_DATA_DIR, so one variable isolates a run."""
    resolved = _resolved_cache_paths(
        {"SKA": "/tmp/fake-ska", "ASTROMON_DATA_DIR": "/tmp/isolated-catalogs"}
    )

    assert resolved["archive"] == "/tmp/isolated-catalogs/xray_observations"


def test_astromon_data_dir_moves_every_catalog_cache():
    """No cache may be left behind pointing at the default location."""
    resolved = _resolved_cache_paths(
        {"SKA": "/tmp/fake-ska", "ASTROMON_DATA_DIR": "/tmp/isolated-catalogs"}
    )

    assert resolved["all"], "no cache constants found to check"
    stragglers = {
        name: path
        for name, path in resolved["all"].items()
        if not path.startswith("/tmp/isolated-catalogs/")
    }
    assert not stragglers, f"still under the default location: {stragglers}"


_RFC_LANDING_HTML = """
<html><body>
<B><A HREF="/sol/rfc/rfc_2026b"> rfc_2026b catalogue of compact radio sources</A></B>.
Please cite data release rfc_2026b and include DOI 10.25966/dhrk-zh08.
<A HREF="/rfc/rfc_2026a_map.pdf"><IMG SRC="/rfc/rfc_2026a_map.png"></A>
</body></html>
"""


def test_discover_latest_rfc_release_reads_the_announced_release():
    """Discovery reports the release astrogeo announces, not the newest directory."""
    fake_response = Mock(text=_RFC_LANDING_HTML)
    fake_response.raise_for_status = Mock()
    with patch.object(cross_match.requests, "get", return_value=fake_response):
        release = cross_match._discover_latest_rfc_release()
    assert release == "rfc_2026b"


def test_discover_latest_rfc_release_ignores_releases_named_only_as_assets():
    """A map PDF for another release must not be mistaken for the announcement."""
    html = (
        '<A HREF="/sol/rfc/rfc_2025d">rfc_2025d</A>'
        '<A HREF="/rfc/rfc_2026b_map.pdf">map</A>'
    )
    fake_response = Mock(text=html)
    fake_response.raise_for_status = Mock()
    with patch.object(cross_match.requests, "get", return_value=fake_response):
        assert cross_match._discover_latest_rfc_release() == "rfc_2025d"


def test_discover_latest_rfc_release_no_releases_found_raises():
    """A page announcing no release raises rather than returning junk."""
    fake_response = Mock(text="<html><body>nothing here</body></html>")
    fake_response.raise_for_status = Mock()
    with patch.object(cross_match.requests, "get", return_value=fake_response):
        with pytest.raises(RuntimeError, match="No RFC release"):
            cross_match._discover_latest_rfc_release()


def _write_fake_rfc_file(dest, *_args, **_kwargs):
    """Write a minimal valid RFC catalog to *dest* directly."""
    Path(dest).write_text(
        "RFC J1229+0203  3C273B    12 29 06.699755 +02 03 08.59833    0.09   0.17   -0.013\n"
    )


def _urlretrieve_writes_fake_rfc(_url, dest, *_args, **_kwargs):
    """Stand-in for urllib.request.urlretrieve(url, dest, ...)."""
    _write_fake_rfc_file(dest)


def test_get_rfc_no_cache_downloads_and_writes_release_marker(tmp_path):
    """First-ever call with no cache discovers the release, downloads, and records it."""
    cache = tmp_path / "rfc_catalog.txt"
    with (
        patch.object(
            cross_match, "_discover_latest_rfc_release", return_value="rfc_2026b"
        ),
        patch(
            "astromon.cross_match.urllib.request.urlretrieve",
            side_effect=_urlretrieve_writes_fake_rfc,
        ) as mock_urlretrieve,
    ):
        rfc = cross_match.get_rfc(cache_path=cache)

    assert len(rfc) == 1
    assert mock_urlretrieve.call_count == 1
    assert mock_urlretrieve.call_args[0][0] == (
        "http://astrogeo.org/sol/rfc/rfc_2026b/rfc_2026b_cat.txt"
    )
    marker = cache.with_name(cache.name + ".release")
    assert marker.read_text() == "rfc_2026b"


def test_get_rfc_fresh_cache_skips_network_entirely(tmp_path):
    """A cache newer than max_age_days never touches the network."""
    cache = tmp_path / "rfc_catalog.txt"
    _write_fake_rfc_file(cache)
    cache.with_name(cache.name + ".release").write_text("rfc_2026b")

    with (
        patch.object(cross_match, "_discover_latest_rfc_release") as mock_discover,
        patch("astromon.cross_match.urllib.request.urlretrieve") as mock_urlretrieve,
    ):
        rfc = cross_match.get_rfc(cache_path=cache, max_age_days=1)

    assert len(rfc) == 1
    mock_discover.assert_not_called()
    mock_urlretrieve.assert_not_called()


def test_get_rfc_stale_but_same_release_skips_download(tmp_path):
    """When the check interval elapses but the release has not changed, no download happens."""
    cache = tmp_path / "rfc_catalog.txt"
    _write_fake_rfc_file(cache)
    cache.with_name(cache.name + ".release").write_text("rfc_2026b")
    old_time = time.time() - 2 * 86400
    os.utime(cache, (old_time, old_time))

    with (
        patch.object(
            cross_match, "_discover_latest_rfc_release", return_value="rfc_2026b"
        ),
        patch("astromon.cross_match.urllib.request.urlretrieve") as mock_urlretrieve,
    ):
        rfc = cross_match.get_rfc(cache_path=cache, max_age_days=1)

    assert len(rfc) == 1
    mock_urlretrieve.assert_not_called()
    assert cache.stat().st_mtime > old_time


def test_get_rfc_new_release_downloads(tmp_path):
    """A genuinely newer release name triggers a real download and marker update."""
    cache = tmp_path / "rfc_catalog.txt"
    _write_fake_rfc_file(cache)
    cache.with_name(cache.name + ".release").write_text("rfc_2026a")
    old_time = time.time() - 2 * 86400
    os.utime(cache, (old_time, old_time))

    with (
        patch.object(
            cross_match, "_discover_latest_rfc_release", return_value="rfc_2026b"
        ),
        patch(
            "astromon.cross_match.urllib.request.urlretrieve",
            side_effect=_urlretrieve_writes_fake_rfc,
        ) as mock_urlretrieve,
    ):
        cross_match.get_rfc(cache_path=cache, max_age_days=1)

    assert mock_urlretrieve.call_count == 1
    marker = cache.with_name(cache.name + ".release")
    assert marker.read_text() == "rfc_2026b"


def test_get_rfc_download_failure_falls_back_to_stale_cache(tmp_path):
    """A newly discovered but broken release falls back to the existing cache."""
    cache = tmp_path / "rfc_catalog.txt"
    _write_fake_rfc_file(cache)
    cache.with_name(cache.name + ".release").write_text("rfc_2026b")
    old_time = time.time() - 2 * 86400
    os.utime(cache, (old_time, old_time))
    size_before = cache.stat().st_size

    with (
        patch.object(
            cross_match, "_discover_latest_rfc_release", return_value="rfc_2026c"
        ),
        patch(
            "astromon.cross_match.urllib.request.urlretrieve",
            side_effect=OSError("HTTP Error 404: Not Found"),
        ),
        patch.object(cross_match.logger, "warning") as mock_warning,
    ):
        rfc = cross_match.get_rfc(cache_path=cache, max_age_days=1)

    assert len(rfc) == 1
    assert cache.stat().st_size == size_before
    assert cache.with_name(cache.name + ".release").read_text() == "rfc_2026b"
    assert "RFC download failed" in mock_warning.call_args[0][0]


def test_get_rfc_discovery_failure_no_cache_raises(tmp_path):
    """With no cache to fall back to, a discovery failure raises."""
    cache = tmp_path / "rfc_catalog.txt"
    with patch.object(
        cross_match, "_discover_latest_rfc_release", side_effect=OSError("network down")
    ):
        with pytest.raises(OSError, match="network down"):
            cross_match.get_rfc(cache_path=cache)


def test_get_rfc_discovery_failure_with_cache_falls_back(tmp_path):
    """A transient failure checking for a new release falls back to the existing cache."""
    cache = tmp_path / "rfc_catalog.txt"
    _write_fake_rfc_file(cache)
    cache.with_name(cache.name + ".release").write_text("rfc_2026b")
    old_time = time.time() - 2 * 86400
    os.utime(cache, (old_time, old_time))

    with (
        patch.object(
            cross_match,
            "_discover_latest_rfc_release",
            side_effect=OSError("network down"),
        ),
        patch.object(cross_match.logger, "warning") as mock_warning,
    ):
        rfc = cross_match.get_rfc(cache_path=cache, max_age_days=1)

    assert len(rfc) == 1
    assert "RFC release check failed" in mock_warning.call_args[0][0]


def test_get_rfc_download_failure_applies_the_check_interval(tmp_path):
    """A release whose catalog file 404s must not be retried on every call.

    astrogeo.org publishes the release directory before the catalog file, so this
    is a real state, not a hypothetical one.  Retrying per call turns one obsid's
    worth of work into one failed HTTP request per obsid.
    """
    cache = _stale_cache(tmp_path)
    calls = []

    with (
        patch.object(
            cross_match, "_discover_latest_rfc_release", return_value="rfc_2026c"
        ),
        patch.object(
            cross_match.urllib.request, "urlretrieve", _failing_urlretrieve(calls)
        ),
    ):
        first = cross_match.get_rfc(cache_path=cache, max_age_days=1)
        second = cross_match.get_rfc(cache_path=cache, max_age_days=1)

    assert len(first) == 2, "the stale cache is still returned"
    assert len(second) == 2
    assert len(calls) == 1, "a broken endpoint must be attempted once, not per call"


def test_get_rfc_release_check_failure_applies_the_check_interval(tmp_path):
    """Same for the directory listing that discovers the release name."""
    cache = _stale_cache(tmp_path)
    checks = []

    def failing_discover():
        checks.append(1)
        raise urllib.error.URLError("no route to host")

    with patch.object(cross_match, "_discover_latest_rfc_release", failing_discover):
        cross_match.get_rfc(cache_path=cache, max_age_days=1)
        cross_match.get_rfc(cache_path=cache, max_age_days=1)

    assert len(checks) == 1


def test_get_rfc_retries_once_the_interval_has_passed(tmp_path):
    """The interval suppresses retries; it does not disable them."""
    cache = _stale_cache(tmp_path)
    calls = []

    with (
        patch.object(
            cross_match, "_discover_latest_rfc_release", return_value="rfc_2026c"
        ),
        patch.object(
            cross_match.urllib.request, "urlretrieve", _failing_urlretrieve(calls)
        ),
    ):
        cross_match.get_rfc(cache_path=cache, max_age_days=1)
        # Age whatever the first attempt recorded, as a day passing would.
        for marker in tmp_path.iterdir():
            old = time.time() - 10 * 86400
            os.utime(marker, (old, old))
        cross_match.get_rfc(cache_path=cache, max_age_days=1)

    assert len(calls) == 2


def test_get_rfc_keeps_the_cache_when_a_download_is_cut_short(tmp_path):
    """A download that dies midway must not leave a truncated catalog in place.

    urlretrieve writes to its destination as it goes, so downloading straight
    onto the cache means a dropped connection replaces a good catalog with a
    partial one -- and the fallback path then happily parses it.
    """
    cache = _stale_cache(tmp_path)

    def truncating_urlretrieve(url, filename):
        with open(filename, "w") as fh:
            fh.write("# Radio Fundamental Catalogue\nRFC J0000-00")
        raise ConnectionResetError("connection reset by peer")

    with (
        patch.object(
            cross_match, "_discover_latest_rfc_release", return_value="rfc_2026c"
        ),
        patch.object(cross_match.urllib.request, "urlretrieve", truncating_urlretrieve),
    ):
        rfc = cross_match.get_rfc(cache_path=cache, max_age_days=1)

    assert len(rfc) == 2, "the previous catalog survives a failed download"
    assert cache.read_text() == RFC_SAMPLE


def test_get_rfc_raises_when_there_is_no_cache_to_fall_back_on(tmp_path):
    """With nothing cached there is no catalog to return, so the failure surfaces."""
    calls = []

    with (
        patch.object(
            cross_match, "_discover_latest_rfc_release", return_value="rfc_2026c"
        ),
        patch.object(
            cross_match.urllib.request, "urlretrieve", _failing_urlretrieve(calls)
        ),
        pytest.raises(urllib.error.HTTPError),
    ):
        cross_match.get_rfc(cache_path=tmp_path / "absent.txt", max_age_days=1)


# ─── recorded catalog versions ───────────────────────────────────────────────


ICRF3_SAMPLE_ROWS = 2


def _icrf3_cache(tmp_path, version=None, age_days=0):
    """An ICRF3 cache, optionally with a recorded version and an age."""
    cache = tmp_path / "icrf3_catalog.ecsv"
    Table(
        {
            "name": np.array(["ICRF J000020.3-322101", "ICRF J122906.6+020308"]),
            "ra": np.array([0.0848, 187.2779]),
            "dec": np.array([-32.3504, 2.0524]),
        }
    ).write(cache, format="ascii.ecsv", overwrite=True)
    if version is not None:
        cross_match._release_marker(cache).write_text(version)
    if age_days:
        old = time.time() - age_days * 86400
        os.utime(cache, (old, old))
    return cache


def _exploding_vizier(*args, **kwargs):
    raise AssertionError("no catalog download should have been attempted")


def _stub_vizier(calls):
    """A Vizier stand-in returning two ICRF3 rows in VizieR's column names."""

    class StubVizier:
        def __init__(self, *args, **kwargs):
            pass

        def get_catalogs(self, catalog_id):
            calls.append(catalog_id)
            return [
                Table(
                    {
                        "ICRF": np.array(["J000020.3-322101", "J122906.6+020308"]),
                        "_RAJ2000": np.array([0.0848, 187.2779]),
                        "_DEJ2000": np.array([-32.3504, 2.0524]),
                    }
                )
            ]

    return StubVizier


def test_declared_catalog_versions_match_the_identifiers_used_to_fetch_them():
    """The version and the thing it describes must not drift apart.

    Each of these catalogs is pinned by an identifier that already names its
    release -- a VizieR designation, a CDS catalogue number, a Zenodo record, a
    Gaia data-release schema. The recorded version is that identifier, so assert
    it is still a substring of what the code actually fetches rather than a
    string someone updated in one place only.
    """
    assert cross_match._ICRF3_VERSION in cross_match._ICRF3_VIZIER_ID
    assert cross_match._MILLIQUAS_VERSION in cross_match._MILLIQUAS_URL
    assert cross_match._QUAIA_VERSION.split("/")[-1] in cross_match._QUAIA_URL

    for getter in (cross_match.get_gaia_agn_catalog, cross_match.get_gaia_qso_catalog):
        assert f"{cross_match._GAIA_DR3_VERSION}." in inspect.getsource(getter)

    # The crossmatch is pinned to a release on each side, so its version names both
    # and each half has to still match the identifier that side is fetched by.
    milliquas_pin, gaia_pin = cross_match._MILLIQUAS_GAIA_VERSION.split("+")
    assert milliquas_pin in cross_match._MILLIQUAS_VIZIER_ID
    assert gaia_pin in cross_match._GAIA_DR3_VIZIER_ID


def test_icrf3_is_not_refetched_because_the_cache_got_old(tmp_path):
    """ICRF3 is a published realisation; age says nothing about it being current.

    A new realisation means a new VizieR designation, which means a code change,
    which the recorded version detects. A timer only forced a pointless refetch.
    """
    cache = _icrf3_cache(tmp_path, version=cross_match._ICRF3_VERSION, age_days=400)

    with patch.object(cross_match, "Vizier", _exploding_vizier):
        icrf3 = cross_match.get_icrf3(cache_path=cache)

    assert len(icrf3) == ICRF3_SAMPLE_ROWS


def test_cache_is_refetched_when_the_recorded_version_is_not_the_declared_one(tmp_path):
    """Bumping the pinned identifier is what invalidates a cache."""
    cache = _icrf3_cache(tmp_path, version="J/A+A/000/A000")
    calls = []

    with patch.object(cross_match, "Vizier", _stub_vizier(calls)):
        icrf3 = cross_match.get_icrf3(cache_path=cache)

    assert calls == [cross_match._ICRF3_VIZIER_ID]
    assert len(icrf3) == ICRF3_SAMPLE_ROWS
    recorded = cross_match._release_marker(cache).read_text().strip()
    assert recorded == cross_match._ICRF3_VERSION


def test_cache_predating_version_recording_is_adopted_not_refetched(tmp_path):
    """There are already hundreds of MB of these caches on disk with no marker.

    Treating an unrecorded version as a mismatch would re-download every one of
    them. Adopt what is there instead -- and do not write a marker claiming a
    version that was never verified.
    """
    cache = _icrf3_cache(tmp_path, version=None, age_days=400)

    with patch.object(cross_match, "Vizier", _exploding_vizier):
        cross_match.get_icrf3(cache_path=cache)

    assert not cross_match._release_marker(cache).exists()


def test_catalog_versions_reports_recorded_unrecorded_and_absent(tmp_path):
    """Provenance of what is actually on disk, for a run to log."""
    recorded = tmp_path / "recorded.fits"
    recorded.write_text("not really a catalog")
    cross_match._release_marker(recorded).write_text("VII/294")
    unrecorded = tmp_path / "unrecorded.fits"
    unrecorded.write_text("not really a catalog")

    versions = cross_match.catalog_versions(
        {
            "Recorded": recorded,
            "Unrecorded": unrecorded,
            "Absent": tmp_path / "absent.fits",
        }
    )

    assert versions == {"Recorded": "VII/294", "Unrecorded": None, "Absent": None}


def test_catalog_versions_covers_every_locally_cached_catalog():
    """A new local cache should show up here without anyone remembering to add it."""
    assert set(cross_match.CATALOG_CACHE_PATHS) == {
        "RFC",
        "ICRF3",
        "Milliquas",
        "Quaia",
        "GaiaAGN",
        "GaiaQSO",
        "MilliquasGaia",
    }
