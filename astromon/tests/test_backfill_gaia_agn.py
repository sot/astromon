"""Tests for astromon.scripts.maintenance.backfill_gaia_agn.

Regression coverage for the bug caught in review: run_xcorr_for_obsids used to
hand-set candidates["x_id"] per detect_method (wrong column name -- the schema is
celldetect_x_id) inside a loop wrapped in a broad except, so a real obsid with more
than one detect_method silently produced xcorr rows for at most one of them.
"""

import numpy as np
from astropy.table import Table, vstack

from astromon import db
from astromon.scripts.maintenance import backfill_gaia_agn as bf

OBSID = 88101


def _obs_row(obsid=OBSID, ra=150.0, dec=2.0):
    return Table(
        {
            "obsid": [obsid],
            "detector": ["ACIS-S"],
            "target": ["test"],
            "grating": ["NONE"],
            "sim_z": [-190.143],
            "date_obs": ["2020-01-01T00:00:00"],
            "tstart": [0.0],
            "ascdsver": ["10.0"],
            "ra": [ra],
            "dec": [dec],
            "roll": [0.0],
            "category_id": [50],
            "version": [10.0],
        }
    )


def _xray_row(obsid, source_id, detect_method, ra, dec):
    return Table(
        {
            "obsid": [obsid],
            "id": [source_id],
            "ra": [ra],
            "dec": [dec],
            "net_counts": [500.0],
            "y_angle": [0.0],
            "z_angle": [0.0],
            "r_angle": [10.0],
            "snr": [10.0],
            "near_neighbor_dist": [60.0],
            "acis_streak": [False],
            "grating_arm": [False],
            "brightest": [False],
            "caldb_version": ["4.10.0"],
            "detect_method": [detect_method],
        }
    )


def _old_gaia_row(obsid, ra, dec, name="GaiaAGN-1"):
    return Table(
        {
            "obsid": [obsid],
            "id": [0],
            "catalog": ["GaiaAGN"],
            "name": [name],
            "ra": [ra],
            "dec": [dec],
            "mag": [15.0],
        }
    )


def test_build_gaia_cat_src_anchors_to_celldetect_not_x_id():
    """The rebuilt cat_src must anchor via celldetect_x_id, not the removed x_id."""
    ra, dec = 150.0, 2.0
    new_xray = vstack(
        [
            _xray_row(OBSID, 3, "celldetect", ra, dec),
            _xray_row(OBSID, 9, "gaussian_detect", ra, dec),
        ],
        metadata_conflicts="silent",
    )
    new_obs = _obs_row(ra=ra, dec=dec)
    old_gaia = _old_gaia_row(OBSID, ra + 1.0 / 3600, dec)
    existing_cat = Table(dtype=db.ASTROMON_CAT_SRC_DTYPE)

    new_cat = bf.build_gaia_cat_src(old_gaia, new_xray, new_obs, existing_cat)

    assert "celldetect_x_id" in new_cat.colnames
    assert "x_id" not in new_cat.colnames
    assert list(new_cat["celldetect_x_id"]) == [3]


def test_run_xcorr_for_obsids_populates_both_detect_methods():
    """The exact failure mode from review: xcorr populated for only one method."""
    ra, dec = 150.0, 2.0
    new_xray = vstack(
        [
            _xray_row(OBSID, 3, "celldetect", ra, dec),
            _xray_row(OBSID, 9, "gaussian_detect", ra, dec),
        ],
        metadata_conflicts="silent",
    )
    new_obs = _obs_row(ra=ra, dec=dec)
    old_gaia = _old_gaia_row(OBSID, ra + 1.0 / 3600, dec)
    existing_cat = Table(dtype=db.ASTROMON_CAT_SRC_DTYPE)

    gaia_cat = bf.build_gaia_cat_src(old_gaia, new_xray, new_obs, existing_cat)
    xcorr = bf.run_xcorr_for_obsids(
        np.unique(gaia_cat["obsid"]), gaia_cat, new_xray, new_obs
    )

    methods = sorted(str(m) for m in xcorr["detect_method"])
    assert methods == ["celldetect", "gaussian_detect"], (
        f"expected a gaia_agn xcorr row for both detect_methods, got {methods}"
    )
    assert np.all(np.array(xcorr["select_name"]).astype(str) == "gaia_agn")
