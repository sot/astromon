"""Tests for astromon.scripts.analysis.solo_catalog_matches."""

import numpy as np
from astropy.table import Table

from astromon import db


def _cat_src_row(obsid, celldetect_x_id, catalog, y_angle=0.0, z_angle=0.0):
    row = np.zeros(1, dtype=db.ASTROMON_CAT_SRC_DTYPE)
    row["obsid"] = obsid
    row["id"] = 1
    row["celldetect_x_id"] = celldetect_x_id
    row["catalog"] = catalog
    row["name"] = "test-source"
    row["separation"] = 1.0
    row["y_angle"] = y_angle
    row["z_angle"] = z_angle
    table = Table(row)
    table.convert_bytestring_to_unicode()
    return table


def _xray_row(obsid, id_, snr=10.0, y_angle=0.0, z_angle=0.0):
    row = np.zeros(1, dtype=db.ASTROMON_XRAY_SRC_DTYPE)
    row["obsid"] = obsid
    row["id"] = id_
    row["snr"] = snr
    row["y_angle"] = y_angle
    row["z_angle"] = z_angle
    row["detect_method"] = "celldetect"
    table = Table(row)
    table.convert_bytestring_to_unicode()
    return table


def _obs_row(obsid, target="Test Target", tstart=500000000.0):
    row = np.zeros(1, dtype=db.ASTROMON_OBS_DTYPE)
    row["obsid"] = obsid
    row["target"] = target
    row["tstart"] = tstart
    table = Table(row)
    table.convert_bytestring_to_unicode()
    return table


def test_build_solo_matches_reads_celldetect_x_id_not_x_id():
    """astromon_cat_src has celldetect_x_id, not the retired x_id column.

    build_solo_matches used to read cat_src["x_id"], which raised KeyError on
    a table from the current schema. This confirms it reads the real column
    and produces a match for a celldetect source with an RFC counterpart at
    the same y/z angle (zero separation).
    """
    from astromon.scripts.analysis.solo_catalog_matches import build_solo_matches

    obsid = 100
    cat_src = _cat_src_row(obsid, celldetect_x_id=1, catalog="RFC")
    xray_src = _xray_row(obsid, id_=1)
    obs = _obs_row(obsid)

    matches = build_solo_matches(cat_src, xray_src, obs)

    assert len(matches) == 1
    assert matches["x_id"][0] == 1
    assert matches["catalog"][0].strip() == "RFC"
    assert matches["dr"][0] == 0.0


def test_build_solo_matches_ignores_non_celldetect_xray_sources():
    """cat_src's celldetect_x_id only ever matches a celldetect xray row."""
    from astromon.scripts.analysis.solo_catalog_matches import build_solo_matches

    obsid = 100
    cat_src = _cat_src_row(obsid, celldetect_x_id=1, catalog="RFC")
    xray_src = _xray_row(obsid, id_=1)
    xray_src["detect_method"] = "gaussian_detect"
    obs = _obs_row(obsid)

    try:
        build_solo_matches(cat_src, xray_src, obs)
    except RuntimeError as exc:
        assert "No matches found" in str(exc)
    else:
        raise AssertionError("expected no celldetect match to raise RuntimeError")
