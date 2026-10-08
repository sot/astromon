"""Tests for discovering new, completed observations via the local ocat cache."""

from astropy.table import Table
from cxotime import CxoTime

from astromon.scripts.maintenance import find_new_observations as fno


def _ocat(rows):
    return Table(rows=rows, names=["obsid", "status", "start_date"])


def test_an_obsid_observed_after_since_is_included():
    ocat = _ocat([(100, "archived", "2026-09-20 00:00:00")])

    result = fno.find_new_candidate_obsids(
        ocat, known_obsids=set(), since=CxoTime("2026-09-19")
    )

    assert result == [100]


def test_an_obsid_observed_before_since_is_excluded():
    ocat = _ocat([(100, "archived", "2026-09-18 00:00:00")])

    result = fno.find_new_candidate_obsids(
        ocat, known_obsids=set(), since=CxoTime("2026-09-19")
    )

    assert result == []


def test_an_obsid_exactly_at_since_is_excluded():
    """`since` is normally the latest date_obs already known -- that exact
    timestamp belongs to an already-known obsid, so the boundary itself must
    not be re-included."""
    ocat = _ocat([(100, "archived", "2026-09-19 12:00:00")])

    result = fno.find_new_candidate_obsids(
        ocat, known_obsids=set(), since=CxoTime("2026-09-19 12:00:00")
    )

    assert result == []


def test_scheduled_or_unobserved_status_is_excluded():
    """A real, not-yet-happened ocat row often carries the literal string
    '0' as its start_date (no real date assigned yet), not a parseable
    date -- these statuses must be dropped before any date parsing is
    attempted, not filtered out only afterward."""
    ocat = _ocat(
        [
            (100, "scheduled", "0"),
            (101, "unobserved", "0"),
            (102, "untriggered", "0"),
        ]
    )

    result = fno.find_new_candidate_obsids(
        ocat, known_obsids=set(), since=CxoTime("2026-09-19")
    )

    assert result == []


def test_canceled_or_discarded_status_is_excluded():
    ocat = _ocat(
        [
            (100, "canceled", "2026-09-20 00:00:00"),
            (101, "discarded", "2026-09-20 00:00:00"),
        ]
    )

    result = fno.find_new_candidate_obsids(
        ocat, known_obsids=set(), since=CxoTime("2026-09-19")
    )

    assert result == []


def test_an_observed_row_with_a_placeholder_start_date_is_excluded():
    """In practice this is ocat's engineering/replay ('ER' type, blank
    target_name) rows: ~95% of 'observed'-status rows carry the literal
    string '0' instead of a real start_date, unlike 'archived' rows, which
    never do. Without a real date there is no way to tell whether this is
    new, so it is dropped rather than guessed at -- obsid >= 38000 rejects
    most of these anyway via Observation.is_selected, but ocat assigns some
    of this block below that threshold too, so the date check must not
    assume the obsid-range filter alone will catch every one of them."""
    ocat = _ocat([(100, "observed", "0")])

    result = fno.find_new_candidate_obsids(
        ocat, known_obsids=set(), since=CxoTime("2026-09-19")
    )

    assert result == []


def test_observed_status_is_included_same_as_archived():
    """'observed' means the exposure happened but products may not be fully
    archived/calibrated yet -- process_one_obsid can still attempt it, same
    as any other obsid; a real download failure is handled downstream, not
    pre-filtered here."""
    ocat = _ocat([(100, "observed", "2026-09-20 00:00:00")])

    result = fno.find_new_candidate_obsids(
        ocat, known_obsids=set(), since=CxoTime("2026-09-19")
    )

    assert result == [100]


def test_an_already_known_obsid_is_excluded_even_if_the_date_matches():
    ocat = _ocat([(100, "archived", "2026-09-20 00:00:00")])

    result = fno.find_new_candidate_obsids(
        ocat, known_obsids={100}, since=CxoTime("2026-09-19")
    )

    assert result == []


def test_results_are_sorted_plain_ints():
    ocat = _ocat(
        [
            (300, "archived", "2026-09-21 00:00:00"),
            (100, "archived", "2026-09-20 00:00:00"),
            (200, "observed", "2026-09-22 00:00:00"),
        ]
    )

    result = fno.find_new_candidate_obsids(
        ocat, known_obsids=set(), since=CxoTime("2026-09-19")
    )

    assert result == [100, 200, 300]
    assert all(type(o) is int for o in result)


def test_no_candidates_returns_an_empty_list():
    ocat = _ocat([(100, "archived", "2026-09-10 00:00:00")])

    result = fno.find_new_candidate_obsids(
        ocat, known_obsids=set(), since=CxoTime("2026-09-19")
    )

    assert result == []
