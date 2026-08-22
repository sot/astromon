"""Tests for telling a VizieR server error apart from a field with no counterparts.

astroquery 0.4.11 reads neither of the two signals VizieR uses to report failure,
so an overloaded server arrives looking exactly like an empty sky. For a
cross-match that is the worst available confusion: the result is written down as
"no counterpart here" and the run reports success. Upstream fixed it in
astropy/astroquery#3632, merged 2026-07-24 for 0.4.12, which is not released --
and not in the newest published 0.4.12.dev snapshot either.
"""

import numpy as np
import pytest
from astropy import coordinates as coords
from astropy import units as u
from astroquery.vizier import Vizier

from astromon import cross_match

# A VizieR error response: HTTP 200, a VOTable carrying QUERY_STATUS=ERROR and no
# table at all. This is what an overloaded server sends.
ERROR_VOTABLE = b"""<?xml version="1.0" encoding="UTF-8"?>
<VOTABLE version="1.4" xmlns="http://www.ivoa.net/xml/VOTable/v1.3">
<INFO ID="Error" name="QUERY_STATUS" value="ERROR">****Service unavailable: overloaded</INFO>
</VOTABLE>
"""

# A successful query over a field with no counterparts: same empty result, but the
# server says so. The whole point is that these two must not be conflated.
OK_EMPTY_VOTABLE = b"""<?xml version="1.0" encoding="UTF-8"?>
<VOTABLE version="1.4" xmlns="http://www.ivoa.net/xml/VOTable/v1.3">
<INFO name="QUERY_STATUS" value="OK"/>
<RESOURCE name="yCat">
<TABLE name="II/246/out"><FIELD name="RAJ2000" datatype="double" unit="deg"/>
<DATA><TABLEDATA></TABLEDATA></DATA></TABLE>
</RESOURCE>
</VOTABLE>
"""

OK_ONE_ROW_VOTABLE = b"""<?xml version="1.0" encoding="UTF-8"?>
<VOTABLE version="1.4" xmlns="http://www.ivoa.net/xml/VOTable/v1.3">
<INFO name="QUERY_STATUS" value="OK"/>
<RESOURCE name="yCat">
<TABLE name="II/246/out"><FIELD name="RAJ2000" datatype="double" unit="deg"/>
<DATA><TABLEDATA><TR><TD>187.2779</TD></TR></TABLEDATA></DATA></TABLE>
</RESOURCE>
</VOTABLE>
"""

# What astroquery 0.4.11 does with anything that is not a VOTable, TSV or FITS
# payload: _parse_result has no else branch, so it falls through to None.
HTML_ERROR_PAGE = b"<html><head><title>503 Service Unavailable</title></head></html>"

# A row-limited query: VizieR still returns rows (up to the limit), but the
# response also carries an OVERFLOW status saying the result was truncated.
OVERFLOW_VOTABLE = b"""<?xml version="1.0" encoding="UTF-8"?>
<VOTABLE version="1.4" xmlns="http://www.ivoa.net/xml/VOTable/v1.3">
<INFO name="QUERY_STATUS" value="OVERFLOW">truncated at 50 rows</INFO>
<RESOURCE name="yCat">
<TABLE name="II/246/out"><FIELD name="RAJ2000" datatype="double" unit="deg"/>
<DATA><TABLEDATA><TR><TD>187.2779</TD></TR></TABLEDATA></DATA></TABLE>
</RESOURCE>
</VOTABLE>
"""

# What a real overflowed response looks like: an "OK" INFO first, then a
# separate "OVERFLOW" INFO -- reading only the first QUERY_STATUS INFO element
# (as astroquery does, and as _vizier_query_info used to) reads "OK" and never
# sees the OVERFLOW that says the result is incomplete.
OK_THEN_OVERFLOW_VOTABLE = b"""<?xml version="1.0" encoding="UTF-8"?>
<VOTABLE version="1.4" xmlns="http://www.ivoa.net/xml/VOTable/v1.3">
<INFO name="QUERY_STATUS" value="OK"/>
<INFO name="QUERY_STATUS" value="OVERFLOW">truncated at 50 rows</INFO>
<RESOURCE name="yCat">
<TABLE name="II/246/out"><FIELD name="RAJ2000" datatype="double" unit="deg"/>
<DATA><TABLEDATA><TR><TD>187.2779</TD></TR></TABLEDATA></DATA></TABLE>
</RESOURCE>
</VOTABLE>
"""


class _StubResponse:
    def __init__(self, content):
        self.content = content
        self.status_code = 200

    def raise_for_status(self):
        return None


def _vizier_returning(content):
    """A real Vizier, with only the HTTP round trip replaced."""
    vizier = Vizier()
    vizier.query_region_async = lambda *args, **kwargs: _StubResponse(content)
    return vizier


POSITION = coords.SkyCoord([187.2779], [2.0524], unit="deg")


def _query(content):
    return cross_match._query_vizier_region(
        _vizier_returning(content), POSITION, 3 * u.arcsec, "II/246"
    )


def test_server_error_raises_instead_of_looking_like_an_empty_field():
    with pytest.raises(cross_match.VizierServerError, match="II/246"):
        _query(ERROR_VOTABLE)


def test_server_error_message_carries_what_the_server_said():
    """The operator needs to know it was overload, not a coding mistake."""
    with pytest.raises(cross_match.VizierServerError, match="overloaded"):
        _query(ERROR_VOTABLE)


def test_a_genuinely_empty_field_is_still_empty():
    """The distinction only means something if the honest empty case survives."""
    result = _query(OK_EMPTY_VOTABLE)

    assert len(result) == 0


def test_a_successful_query_is_parsed_as_before():
    """The error check must not disturb the normal path it wraps."""
    result = _query(OK_ONE_ROW_VOTABLE)

    tables = list(result)
    assert len(tables) == 1
    assert np.isclose(tables[0]["RAJ2000"][0], 187.2779)


def test_an_unparseable_response_raises_rather_than_returning_none():
    """astroquery's _parse_result falls through to None on unrecognised content.

    That is still true on astroquery main, so it is not covered by the upstream
    fix: an HTML error page becomes None, and None then propagates into code that
    expects a TableList.
    """
    with pytest.raises(cross_match.VizierServerError, match="could not be parsed"):
        _query(HTML_ERROR_PAGE)


def test_query_status_is_read_without_parsing_rows():
    """The status check is a cheap read of the INFO elements, not a second parse."""
    assert cross_match._votable_query_status(ERROR_VOTABLE) == "ERROR"
    assert cross_match._votable_query_status(OK_EMPTY_VOTABLE) == "OK"
    assert cross_match._votable_query_status(HTML_ERROR_PAGE) is None


def test_overflow_raises_instead_of_returning_a_truncated_result():
    """A row-limited result must not be recorded as a complete one."""
    with pytest.raises(cross_match.VizierServerError, match="truncated"):
        _query(OVERFLOW_VOTABLE)


def test_overflow_after_ok_is_still_caught():
    """VizieR can emit an OK QUERY_STATUS INFO followed by a separate OVERFLOW
    one for the same response. Reading only the first INFO element reads "OK"
    and never sees the OVERFLOW -- the truncated result then looks complete.
    """
    assert cross_match._votable_query_status(OK_THEN_OVERFLOW_VOTABLE) == "OVERFLOW"
    with pytest.raises(cross_match.VizierServerError, match="truncated"):
        _query(OK_THEN_OVERFLOW_VOTABLE)


def test_get_vizier_requests_the_unlimited_row_count():
    """get_vizier's astroquery path must not rely on Vizier's default row_limit
    of 50 -- a field with more real matches than that would otherwise be
    silently truncated and reported as a complete result.

    ``Vizier`` (imported from astroquery) is a pre-built singleton instance,
    not a class -- calling it goes through ``VizierClass.__call__``, so the
    module-level name in cross_match is stubbed directly rather than trying to
    intercept a constructor that never runs.
    """
    pos = coords.SkyCoord([187.2779], [2.0524], unit="deg", obstime="2020-01-01")

    captured = {}

    def fake_vizier(*args, **kwargs):
        captured.update(kwargs)
        return object()

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(cross_match, "Vizier", fake_vizier)
        mp.setattr(cross_match, "_query_vizier_region", lambda vizier, *a, **k: [])
        cross_match.get_vizier(pos, "2MASS", "II/246/out", ["_2MASS"], {}, raw=True)

    assert captured.get("row_limit") == -1
