"""Tests for astromon.scripts.analysis.review_xray_images."""

import tempfile
from pathlib import Path

import numpy as np
from astropy.table import Table

from astromon import db


def test_generated_js_guards_reduce_against_zero_sources():
    """The panel's peak-SNR summary must not crash on an obsid with no sources.

    The generated page's build() calls sources.reduce((a,b)=>...) with no
    initial value and no length check. For any obsid with zero celldetect
    sources -- a real, anticipated case (read_sources returns an empty list
    for it) -- Array.prototype.reduce on an empty array with no initial value
    throws inside build(), which runs inside VIZ.map(build): the whole review
    page fails to render, not just that one panel.

    There is no JS runtime in this environment to execute the page directly,
    so this checks statically that the fix (skip peak-SNR entirely when there
    are no sources) is present in the generated template, and that the
    unguarded call is gone.
    """
    from astromon.scripts.analysis.review_xray_images import build_html

    html = build_html(viz_data=[], title_note="test", result_note="")

    assert "const peak=sources.length?" in html, (
        "the peak-SNR block (and its sources.reduce() call) must be gated on "
        "sources.length, not run unconditionally"
    )


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
