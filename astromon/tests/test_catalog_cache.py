"""Tests for the locally cached catalogs' refresh policy.

get_rfc checks astrogeo.org for a newer release at most once per `max_age_days`,
so that a pipeline run does not re-check once per obsid.  These tests cover the
case the interval originally missed -- a check or download that *fails* -- which
matters as soon as anything loops over thousands of obsids in one process.
"""

import inspect
import os
import time
import urllib.error
from unittest.mock import patch

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
