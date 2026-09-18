"""Tests for astromon.scripts.analysis.review_xray_images."""

import tempfile
from pathlib import Path

import numpy as np
from astropy.table import Table

from astromon import db


def test_load_cat_src_rfc_reads_celldetect_x_id_not_x_id():
    """astromon_cat_src's physical HDF5 column is celldetect_x_id, not x_id.

    load_cat_src_rfc used to read row["x_id"] directly off the raw pytables
    row, which no longer exists once the column was renamed. This confirms it
    reads the real column and only picks up RFC/ICRS rows.
    """
    from astromon.scripts.analysis.review_xray_images import load_cat_src_rfc

    obsid = 100
    rows = np.zeros(2, dtype=db.ASTROMON_CAT_SRC_DTYPE)
    rows["obsid"] = obsid
    rows["celldetect_x_id"] = [1, 2]
    rows["catalog"] = ["RFC", "Tycho2"]

    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "test.h5"
        db.save("astromon_cat_src", Table(rows), dbfile)

        rfc_ids = load_cat_src_rfc(obsid, dbfile)

    assert rfc_ids == {1}
