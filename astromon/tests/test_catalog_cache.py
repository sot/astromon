"""Tests for the locally cached catalogs' refresh policy.

get_rfc checks astrogeo.org for a newer release at most once per `max_age_days`,
so that a pipeline run does not re-check once per obsid.  These tests cover the
case the interval originally missed -- a check or download that *fails* -- which
matters as soon as anything loops over thousands of obsids in one process.
"""

import os
import time
import urllib.error
from unittest.mock import patch

import pytest

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
