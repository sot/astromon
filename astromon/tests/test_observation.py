import re
from unittest.mock import Mock

import pytest

from astromon import observation

# The real HRC plate scale, per the Chandra Proposers' Observatory Guide. filter_events and
# filter_sources both filter events/sources in a circle around the optical axis, and need to
# convert the filter radius from arcsec to detector pixels. Using the wrong scale here does not
# raise an exception: it silently changes the filtered region size (see astromon/observation.py).
HRC_ARCSEC_PER_PIXEL = 0.13180
ACIS_ARCSEC_PER_PIXEL = 0.5

CIRCLE_RE = re.compile(r"circle\([^,]+,[^,]+,([^)]+)\)")


def _fake_obs(is_hrc):
    obs = Mock()
    obs.is_hrc = is_hrc
    obs.ciao = Mock()
    obs.ciao.pget = Mock(side_effect=["10.0", "-20.0", "4096.5", "4096.5"])
    return obs


def _dmcopy_circle_radius(ciao_mock):
    (dmcopy_call,) = [c for c in ciao_mock.call_args_list if c.args[0] == "dmcopy"]
    match = CIRCLE_RE.search(dmcopy_call.args[1])
    return float(match.group(1))


@pytest.mark.parametrize(
    "is_hrc,arcsec_per_pixel",
    [(True, HRC_ARCSEC_PER_PIXEL), (False, ACIS_ARCSEC_PER_PIXEL)],
)
def test_filter_events_pixel_scale(is_hrc, arcsec_per_pixel):
    obs = _fake_obs(is_hrc)
    inputs = {"events": ["evt2.fits"], "fov": ["fov1.fits"]}
    outputs = {"events": "evt2_filtered.fits.gz"}

    observation.filter_events.func(obs, inputs, outputs)

    radius = 180  # matches the fixed radius in filter_events
    assert _dmcopy_circle_radius(obs.ciao) == pytest.approx(radius / arcsec_per_pixel)


@pytest.mark.parametrize(
    "is_hrc,arcsec_per_pixel",
    [(True, HRC_ARCSEC_PER_PIXEL), (False, ACIS_ARCSEC_PER_PIXEL)],
)
def test_filter_sources_pixel_scale(is_hrc, arcsec_per_pixel):
    obs = _fake_obs(is_hrc)
    inputs = {"events": "evt2_filtered.fits.gz", "src": "baseline.src"}
    outputs = {"src": "filtered.src"}

    observation.filter_sources.func(obs, inputs, outputs)

    radius = 180  # matches the fixed radius in filter_sources
    assert _dmcopy_circle_radius(obs.ciao) == pytest.approx(radius / arcsec_per_pixel)
