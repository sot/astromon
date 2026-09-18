"""Tests for astromon.scripts.backfill_icrf2.

Regression coverage for two bugs caught in review:

- db.save(..., replace_keys=("catalog",)) deletes every existing ICRF2 row
  globally, but the xcorr rebuild only covers obsids backfill() actually
  produced a new row for. An obsid that had an ICRF2 candidate before but
  gets none this time loses its cat_src row with no matching xcorr cleanup,
  leaving a dangling icrf2 xcorr row.
- The per-detect_method loop called the since-deleted
  cross_match.remap_x_id_to_sources, which raised AttributeError on every
  call -- silently swallowed by a broad except, so cat_src got ICRF2 rows
  but xcorr stayed empty for every detect_method. A test that only checks
  "doesn't raise" would not have caught that; this checks actual xcorr row
  counts per detect_method.
"""

import sys
from unittest.mock import patch

import numpy as np
from astropy.table import Table

from astromon import cross_match, db
from astromon.scripts import backfill_icrf2 as bf

OBSID_KEEP = 3001
OBSID_DROP = 3002


def _cat_src_row(obsid, source_id, catalog="ICRF2"):
    row = np.zeros(1, dtype=db.ASTROMON_CAT_SRC_DTYPE)
    row["obsid"] = obsid
    row["id"] = source_id
    row["catalog"] = catalog
    return Table(row)


def _xcorr_row(obsid, c_id, x_id, select_name="icrf2"):
    row = np.zeros(1, dtype=db.ASTROMON_XCORR_DTYPE)
    row["obsid"] = obsid
    row["c_id"] = c_id
    row["x_id"] = x_id
    row["select_name"] = select_name
    return Table(row)


def test_main_drops_stale_xcorr_for_an_obsid_that_lost_its_icrf2_candidate(
    tmp_path, monkeypatch
):
    dbfile = tmp_path / "astromon.h5"
    db.create_empty_tables(dbfile)

    db.save("astromon_cat_src", _cat_src_row(OBSID_KEEP, 0), dbfile)
    db.save("astromon_cat_src", _cat_src_row(OBSID_DROP, 0), dbfile)
    db.save("astromon_xcorr", _xcorr_row(OBSID_KEEP, 0, 1), dbfile)
    db.save("astromon_xcorr", _xcorr_row(OBSID_DROP, 0, 2), dbfile)

    # backfill() only produces a fresh row for OBSID_KEEP this run.
    new_cat_src = _cat_src_row(OBSID_KEEP, 0)
    new_xcorr = _xcorr_row(OBSID_KEEP, 0, 1)

    monkeypatch.setattr(sys, "argv", ["backfill_icrf2", "--dbfile", str(dbfile)])
    with patch.object(bf, "backfill", return_value=(new_cat_src, new_xcorr)):
        bf.main()

    xcorr = db.get_table("astromon_xcorr", dbfile)
    cat_src = db.get_table("astromon_cat_src", dbfile)

    assert OBSID_KEEP in np.asarray(xcorr["obsid"]), (
        "the obsid that got a real ICRF2 match must keep its xcorr row"
    )
    assert OBSID_DROP not in np.asarray(xcorr["obsid"]), (
        "the obsid that lost its ICRF2 candidate must not keep a dangling "
        "icrf2 xcorr row"
    )
    assert OBSID_DROP not in np.asarray(cat_src["obsid"])


OBSID = 88001


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
    )[0]


def _xray_row(obsid, source_id, detect_method, ra, dec):
    row = Table(np.zeros(1, dtype=db.ASTROMON_XRAY_SRC_DTYPE))
    row["obsid"] = obsid
    row["id"] = source_id
    row["ra"] = ra
    row["dec"] = dec
    row["net_counts"] = 500.0
    row["snr"] = 10.0
    row["r_angle"] = 10.0
    row["near_neighbor_dist"] = 60.0
    row["acis_streak"] = False
    row["detect_method"] = detect_method
    return row


def _seed_db(tmp_path, obs_row, xray_rows):
    dbfile = tmp_path / "icrf2.h5"
    db.save("astromon_obs", Table([obs_row]), dbfile, ignore_obsid=True)
    db.save(
        "astromon_xray_src",
        Table(np.concatenate(xray_rows)),
        dbfile,
        ignore_obsid=True,
    )
    return dbfile


def _icrf2_candidate_near(ra, dec):
    """One ICRF2 source ~1 arcsec from (ra, dec), shaped like a get_vizier result."""
    return Table(
        {
            "catalog": ["ICRF2"],
            "name": ["ICRF2-test"],
            "ra": [ra + 1.0 / 3600 / np.cos(np.radians(dec))],
            "dec": [dec],
            "mag": [np.nan],
        }
    )


def test_backfill_populates_xcorr_for_every_detect_method(tmp_path):
    """The exact failure mode Codex's review caught: cat_src populated, xcorr empty.

    One obsid with celldetect and gaussian_detect sources at the same position.
    ICRF2 rough-matches once against celldetect, then compute_cross_matches must
    run over both detect_methods -- not just the one used to build the rough match.
    """
    ra, dec = 150.0, 2.0
    dbfile = _seed_db(
        tmp_path,
        _obs_row(ra=ra, dec=dec),
        [
            _xray_row(OBSID, 1, "celldetect", ra, dec),
            _xray_row(OBSID, 1, "gaussian_detect", ra, dec),
        ],
    )

    with patch.object(cross_match, "_get", return_value=_icrf2_candidate_near(ra, dec)):
        cat_src, xcorr = bf.backfill(dbfile)

    assert len(cat_src) == 1
    assert cat_src["catalog"][0] == "ICRF2"
    assert "celldetect_x_id" in cat_src.colnames
    assert cat_src["celldetect_x_id"][0] == 1

    methods = sorted(str(m) for m in xcorr["detect_method"])
    assert methods == ["celldetect", "gaussian_detect"], (
        f"expected an icrf2 xcorr row for both detect_methods, got {methods}"
    )
    assert np.all(np.array(xcorr["select_name"]).astype(str) == "icrf2")


def test_backfill_does_not_raise_and_saves_cleanly(tmp_path):
    """backfill() must not raise AttributeError from a stale remap call, and its
    output must be directly usable by db.save (correct column names)."""
    ra, dec = 30.0, -10.0
    dbfile = _seed_db(
        tmp_path,
        _obs_row(obsid=OBSID + 1, ra=ra, dec=dec),
        [_xray_row(OBSID + 1, 5, "celldetect", ra, dec)],
    )

    with patch.object(cross_match, "_get", return_value=_icrf2_candidate_near(ra, dec)):
        cat_src, xcorr = bf.backfill(dbfile)

    assert len(cat_src) == 1
    db.save("astromon_cat_src", cat_src, dbfile=dbfile, replace_keys=("catalog",))
    db.save("astromon_xcorr", xcorr, dbfile=dbfile, select_name_key=True)

    saved_cat = db.get_table("astromon_cat_src", dbfile)
    assert list(saved_cat["celldetect_x_id"]) == [5]
