"""Tests for cross-match selection semantics and anchor invariants.

These cases exercise the matching rules that are easy to accidentally change
when the catalogue source table grows new columns or new detection methods.
Keeping them in their own file makes it easier to review the intended contract:
matching is positional, the celldetect anchor is metadata, and selection still
groups by method as well as source.
"""

from pathlib import Path
from unittest.mock import patch

import numpy as np
from astropy import units as u
from astropy.table import Table, vstack
from cxotime import CxoTime

from astromon import cross_match
from astromon.tests.utils import (
    minimal_cat_src as _minimal_cat_src,
    minimal_obs as _minimal_obs,
    minimal_xray_src as _minimal_xray_src,
    two_method_xray_src as _two_method_xray_src,
)

DATA_DIR = Path(__file__).parent / "data"
VIZIER_QUERY_TIME = CxoTime("2013-02-11T07:39:09")


def test_the_anchor_always_names_a_celldetect_source():
    """Every anchor resolves to a celldetect source in the same observation."""
    cat_src = Table.read(DATA_DIR / "astromon_cat_src.ecsv")
    xray_src = Table.read(DATA_DIR / "astromon_xray_src.ecsv")

    if "detect_method" in xray_src.colnames:
        xray_src = xray_src[xray_src["detect_method"] == "celldetect"]
    celldetect = set(
        zip(
            np.asarray(xray_src["obsid"]).tolist(),
            np.asarray(xray_src["id"]).tolist(),
            strict=True,
        )
    )

    unresolved = [
        (int(obsid), int(anchor))
        for obsid, anchor in zip(
            cat_src["obsid"], cat_src["celldetect_x_id"], strict=True
        )
        if (int(obsid), int(anchor)) not in celldetect
    ]

    assert not unresolved, f"anchors naming no celldetect source: {unresolved[:5]}"


def test_matches_are_unchanged_when_the_anchor_is_scrambled():
    """cat_src's anchor must not decide which X-ray source a catalogue row matches."""
    obs = _minimal_obs(99900)
    xray = _minimal_xray_src(99900)
    good = _minimal_cat_src(99900)

    scrambled = good.copy()
    scrambled["celldetect_x_id"] = 4242

    from_good = cross_match.compute_cross_matches(
        "tycho2", astromon_obs=obs, astromon_xray_src=xray, astromon_cat_src=good
    )
    from_scrambled = cross_match.compute_cross_matches(
        "tycho2", astromon_obs=obs, astromon_xray_src=xray, astromon_cat_src=scrambled
    )

    assert len(from_good) == 1
    for column in ("obsid", "x_id", "c_id", "dr", "detect_method"):
        assert list(from_good[column]) == list(from_scrambled[column]), column


def test_a_method_missing_the_source_yields_no_match_for_that_method():
    """gaussian dropping a source means no gaussian match, not a substituted one."""
    xray = _two_method_xray_src(99900)
    xray["id"] = [1, 2]
    xray["y_angle"] = [0.0, 36.0]

    matches = cross_match.compute_cross_matches(
        "tycho2",
        astromon_obs=_minimal_obs(99900),
        astromon_xray_src=xray,
        astromon_cat_src=_minimal_cat_src(99900),
    )

    assert list(matches["detect_method"]) == ["celldetect"]
    assert list(matches["x_id"]) == [1]


def test_empty_inputs_still_return_the_joined_dtype():
    """The empty path must not reference a column the schema no longer has."""
    empty_cat = _minimal_cat_src(99900)[:0]

    matches = cross_match.compute_cross_matches(
        "tycho2",
        astromon_obs=_minimal_obs(99900),
        astromon_xray_src=_minimal_xray_src(99900),
        astromon_cat_src=empty_cat,
    )

    assert len(matches) == 0
    for column in ("obsid", "x_id", "c_id", "celldetect_x_id"):
        assert column in matches.colnames, column


def test_selection_keeps_one_row_per_method():
    """A selection without detect_method_filter must not collapse the two methods."""
    matches = cross_match.compute_cross_matches(
        "tycho2",
        astromon_obs=_minimal_obs(99900),
        astromon_xray_src=_two_method_xray_src(99900),
        astromon_cat_src=_minimal_cat_src(99900),
    )

    assert sorted(matches["detect_method"]) == ["celldetect", "gaussian_detect"]


def test_selection_still_picks_one_match_per_source_within_a_method():
    """The grouping itself is unchanged: still one best catalogue match per source."""
    cat = vstack(
        [
            _minimal_cat_src(99900, name="near", ra=83.82028),
            _minimal_cat_src(99900, name="far", ra=83.8206),
        ],
        metadata_conflicts="silent",
    )
    cat["id"] = [0, 1]

    matches = cross_match.compute_cross_matches(
        "tycho2",
        astromon_obs=_minimal_obs(99900),
        astromon_xray_src=_minimal_xray_src(99900),
        astromon_cat_src=cat,
    )

    assert len(matches) == 1
    assert matches["name"][0] == "near"


def test_rough_match_deduplicates_repeated_catalog_rows():
    """One row per (catalog source, X-ray source) pair, not per returned copy."""
    sources = Table({"id": [1, 2], "ra": [10.0, 10.0005], "dec": [20.0, 20.0]})
    returned_twice = Table(
        {
            "catalog": ["Tycho2", "Tycho2"],
            "name": ["dup-source", "dup-source"],
            "ra": [10.00025, 10.00025],
            "dec": [20.0, 20.0],
            "mag": [11.0, 11.0],
        }
    )

    with patch.object(cross_match, "_get", return_value=returned_twice):
        result = cross_match.rough_match(
            sources, VIZIER_QUERY_TIME, radius=5 * u.arcsec, catalogs=["Tycho2"]
        )

    pairs = sorted(zip(result["name"], result["x_id"], strict=True))
    assert pairs == [("dup-source", 1), ("dup-source", 2)]


def test_rough_match_keeps_same_named_sources_at_different_positions():
    """Deduplication keys on position too, because SDSS names are not unique."""
    sources = Table({"id": [1], "ra": [10.0], "dec": [20.0]})
    two_objects_one_name = Table(
        {
            "catalog": ["SDSS", "SDSS"],
            "name": ["shared-name", "shared-name"],
            "ra": [10.0002, 10.0004],
            "dec": [20.0, 20.0],
            "mag": [19.0, 20.0],
        }
    )

    with patch.object(cross_match, "_get", return_value=two_objects_one_name):
        result = cross_match.rough_match(
            sources, VIZIER_QUERY_TIME, radius=5 * u.arcsec, catalogs=["SDSS"]
        )

    assert len(result) == 2
    assert sorted(round(float(r), 4) for r in result["ra"]) == [10.0002, 10.0004]
