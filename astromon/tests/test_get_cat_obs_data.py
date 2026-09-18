"""Tests for the per-obsid pipeline entry point in astromon.scripts.get_cat_obs_data."""

import tempfile
from pathlib import Path

import numpy as np
from astropy.table import Table

from astromon import db


def _obs_row(obsid: int) -> Table:
    row = np.zeros(1, dtype=db.ASTROMON_OBS_DTYPE)
    row["obsid"] = obsid
    return Table(row)


def _xray_row(obsid: int, id_: int) -> Table:
    row = np.zeros(1, dtype=db.ASTROMON_XRAY_SRC_DTYPE)
    row["obsid"] = obsid
    row["id"] = id_
    return Table(row)


def _cat_src_row(obsid: int, id_: int, celldetect_x_id: int) -> Table:
    row = np.zeros(1, dtype=db.ASTROMON_CAT_SRC_DTYPE)
    row["obsid"] = obsid
    row["id"] = id_
    row["celldetect_x_id"] = celldetect_x_id
    row["catalog"] = "RFC"
    return Table(row)


def _xcorr_row(obsid: int, c_id: int, x_id: int) -> Table:
    row = np.zeros(1, dtype=db.ASTROMON_XCORR_DTYPE)
    row["obsid"] = obsid
    row["c_id"] = c_id
    row["x_id"] = x_id
    row["select_name"] = "astromon_21"
    return Table(row)


def test_save_drops_stale_cat_src_and_xcorr_on_an_authoritative_empty_rerun():
    """An empty cat_src/xcorr from a genuine full success removes old rows.

    save() has no --skip-catalog-match: reaching this point with ok=True means
    a real catalog-match pass ran, so an empty astromon_cat_src is a real "found
    nothing", not "didn't check". Without this, the obsid never appears in the
    vstacked data handed to db.save (there is nothing to vstack for an empty
    table), and its old matches survive untouched -- while astromon_status (not
    exercised by this test) would still record success.
    """
    from astromon.scripts.get_cat_obs_data import save

    obsid = 555
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "test.h5"
        db.save("astromon_obs", _obs_row(obsid), dbfile)
        db.save("astromon_xray_src", _xray_row(obsid, 1), dbfile)
        db.save("astromon_cat_src", _cat_src_row(obsid, 1, 1), dbfile)
        db.save("astromon_xcorr", _xcorr_row(obsid, 1, 1), dbfile)

        save(
            [
                {
                    "ok": True,
                    "msg": "",
                    "obsid": obsid,
                    "astromon_obs": _obs_row(obsid),
                    "astromon_xray_src": _xray_row(obsid, 1),
                    "astromon_cat_src": db.create_table("astromon_cat_src"),
                    "astromon_xcorr": db.create_table("astromon_xcorr"),
                }
            ],
            dbfile,
        )

        cat_src = db.get_table("astromon_cat_src", dbfile)
        xcorr = db.get_table("astromon_xcorr", dbfile)
        assert obsid not in np.asarray(cat_src["obsid"])
        assert obsid not in np.asarray(xcorr["obsid"])
        # The genuinely new data for this obsid was still saved.
        assert obsid in np.asarray(db.get_table("astromon_obs", dbfile)["obsid"])


def test_save_leaves_other_obsids_cat_src_and_xcorr_alone():
    """Only the obsid(s) in this save() call are affected, not the whole table."""
    from astromon.scripts.get_cat_obs_data import save

    rerun_obsid = 555
    untouched_obsid = 777
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "test.h5"
        for obsid in (rerun_obsid, untouched_obsid):
            db.save("astromon_obs", _obs_row(obsid), dbfile)
            db.save("astromon_xray_src", _xray_row(obsid, 1), dbfile)
            db.save("astromon_cat_src", _cat_src_row(obsid, 1, 1), dbfile)
            db.save("astromon_xcorr", _xcorr_row(obsid, 1, 1), dbfile)

        save(
            [
                {
                    "ok": True,
                    "msg": "",
                    "obsid": rerun_obsid,
                    "astromon_obs": _obs_row(rerun_obsid),
                    "astromon_xray_src": _xray_row(rerun_obsid, 1),
                    "astromon_cat_src": db.create_table("astromon_cat_src"),
                    "astromon_xcorr": db.create_table("astromon_xcorr"),
                }
            ],
            dbfile,
        )

        cat_src = db.get_table("astromon_cat_src", dbfile)
        xcorr = db.get_table("astromon_xcorr", dbfile)
        assert untouched_obsid in np.asarray(cat_src["obsid"])
        assert untouched_obsid in np.asarray(xcorr["obsid"])


def test_save_keeps_genuine_new_cat_src_and_xcorr_rows():
    """A rerun with real new matches saves them normally, not just drops old ones."""
    from astromon.scripts.get_cat_obs_data import save

    obsid = 555
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "test.h5"
        db.save("astromon_obs", _obs_row(obsid), dbfile)
        db.save("astromon_xray_src", _xray_row(obsid, 1), dbfile)
        db.save("astromon_cat_src", _cat_src_row(obsid, 1, 1), dbfile)
        db.save("astromon_xcorr", _xcorr_row(obsid, 1, 1), dbfile)

        new_cat_src = _cat_src_row(obsid, 2, 1)
        new_xcorr = _xcorr_row(obsid, 2, 1)
        save(
            [
                {
                    "ok": True,
                    "msg": "",
                    "obsid": obsid,
                    "astromon_obs": _obs_row(obsid),
                    "astromon_xray_src": _xray_row(obsid, 1),
                    "astromon_cat_src": new_cat_src,
                    "astromon_xcorr": new_xcorr,
                }
            ],
            dbfile,
        )

        cat_src = db.get_table("astromon_cat_src", dbfile)
        assert list(cat_src["id"]) == [2]


def test_split_versions_separates_peak_gaussian():
    """peak_gaussian_detect is run for its .src output, not saved as a method.

    Its fitted position never won on match count; what it carries is the
    peak_offset diagnostic, which get_sources reads off the .src file. So it is
    run but not persisted.
    """
    from astromon.scripts.get_cat_obs_data import _split_versions

    saved, diagnostic = _split_versions(
        ("celldetect", "gaussian_detect", "peak_gaussian_detect")
    )
    assert saved == ["celldetect", "gaussian_detect"]
    assert diagnostic == ["peak_gaussian_detect"]


def test_split_versions_preserves_order_and_handles_absence():
    from astromon.scripts.get_cat_obs_data import _split_versions

    assert _split_versions(("gaussian_detect", "celldetect")) == (
        ["gaussian_detect", "celldetect"],
        [],
    )
    assert _split_versions(("peak_gaussian_detect",)) == ([], ["peak_gaussian_detect"])
    assert _split_versions(()) == ([], [])
