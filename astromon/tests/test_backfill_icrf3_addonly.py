"""Tests for astromon.scripts.maintenance.backfill_icrf3_addonly.

Regression coverage for the bug caught in review: run_icrf3_xcorr used to hand-set
candidates["x_id"] per detect_method (wrong column name -- the schema is
celldetect_x_id) inside a loop wrapped in a broad except, so a real obsid with more
than one detect_method silently produced icrf3 xcorr rows for at most one of them.
"""

import numpy as np
from astropy.table import Table, vstack

from astromon import db
from astromon.scripts.maintenance import backfill_icrf3_addonly as bf

OBSID = 88201


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


def _icrf3_catalog(ra, dec, name="ICRF3-test"):
    """Shaped like get_icrf3()'s return: name, ra, dec."""
    return Table({"name": [name], "ra": [ra], "dec": [dec]})


def test_build_icrf3_cat_src_anchors_to_celldetect_not_x_id():
    """The rebuilt cat_src must anchor via celldetect_x_id, not the removed x_id."""
    ra, dec = 150.0, 2.0
    obspar = _obs_row(ra=ra, dec=dec)
    celldetect_xray = _xray_row(OBSID, 3, "celldetect", ra, dec)
    icrf3 = _icrf3_catalog(ra + 1.0 / 3600, dec)
    existing_cat = Table(dtype=db.ASTROMON_CAT_SRC_DTYPE)

    new_cat = bf.build_icrf3_cat_src(
        icrf3, celldetect_xray, obspar, np.array([OBSID]), existing_cat
    )

    assert "celldetect_x_id" in new_cat.colnames
    assert "x_id" not in new_cat.colnames
    assert list(new_cat["celldetect_x_id"]) == [3]


def test_run_icrf3_xcorr_populates_both_detect_methods():
    """The exact failure mode from review: xcorr populated for only one method."""
    ra, dec = 150.0, 2.0
    obspar = _obs_row(ra=ra, dec=dec)
    celldetect_xray = _xray_row(OBSID, 3, "celldetect", ra, dec)
    full_xray = vstack(
        [celldetect_xray, _xray_row(OBSID, 9, "gaussian_detect", ra, dec)],
        metadata_conflicts="silent",
    )
    icrf3 = _icrf3_catalog(ra + 1.0 / 3600, dec)
    existing_cat = Table(dtype=db.ASTROMON_CAT_SRC_DTYPE)

    icrf3_cat = bf.build_icrf3_cat_src(
        icrf3, celldetect_xray, obspar, np.array([OBSID]), existing_cat
    )
    xcorr = bf.run_icrf3_xcorr(np.array([OBSID]), icrf3_cat, full_xray, obspar)

    methods = sorted(str(m) for m in xcorr["detect_method"])
    assert methods == ["celldetect", "gaussian_detect"], (
        f"expected an icrf3 xcorr row for both detect_methods, got {methods}"
    )
    assert np.all(np.array(xcorr["select_name"]).astype(str) == "icrf3")
