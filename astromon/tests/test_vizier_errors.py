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
    assert cross_match._vizier_query_status(ERROR_VOTABLE) == "ERROR"
    assert cross_match._vizier_query_status(OK_EMPTY_VOTABLE) == "OK"
    assert cross_match._vizier_query_status(HTML_ERROR_PAGE) is None
