"""Tests for repairing y_angle/z_angle/r_angle that disagree with the evt2 pointing.

Context: a ``get_sources`` cache written while sources' angles still came from the
quick-look obspar pointing was later reused and saved back to the database, so
those obsids' rows hold RA/Dec that are right and y_angle/z_angle that are not.
The repair recomputes the angles from each row's own RA/Dec and the pointing in
the obsid's local evt2 header, for exactly the (obsid, detect_method) groups whose
stored angles disagree with that recomputation.
"""

import gzip
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table
from chandra_aca.transform import radec_to_yagzag, yagzag_to_radec
from Quaternion import Quat

from astromon import db
from astromon.scripts.maintenance import fix_pointing_angles as fix

#: Two pointings ~25" apart, roughly what the quick-look obspar and the evt2
#: header disagree by in the affected obsids.
EVT2_POINTING = (150.0, 2.0, 30.0)
OBSPAR_POINTING = (150.0060, 2.0035, 30.004)


def write_evt2(path: Path, pointing: tuple[float, float, float]) -> Path:
    """A gzipped evt2 with an empty primary HDU and an EVENTS header carrying `pointing`."""
    events = fits.BinTableHDU.from_columns(
        [fits.Column(name="TIME", format="D", array=np.arange(1000.0))], name="EVENTS"
    )
    events.header["RA_PNT"], events.header["DEC_PNT"], events.header["ROLL_PNT"] = (
        pointing
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wb") as fh:
        fits.HDUList([fits.PrimaryHDU(), events]).writeto(fh)
    return path


def primary_dir(workdir: Path, obsid: int) -> Path:
    return workdir / f"obs{obsid // 1000:02d}" / str(obsid) / "primary"


def angles(ra, dec, pointing):
    return radec_to_yagzag(np.asarray(ra), np.asarray(dec), Quat(equatorial=pointing))


def xray_rows(obsid, method, ra, dec, pointing, first_id=1):
    rows = Table(np.zeros(len(ra), dtype=db.ASTROMON_XRAY_SRC_DTYPE))
    rows["obsid"] = obsid
    rows["id"] = np.arange(first_id, first_id + len(ra))
    rows["detect_method"] = method
    rows["ra"], rows["dec"] = ra, dec
    rows["snr"] = 7.5
    rows["y_angle"], rows["z_angle"] = angles(ra, dec, pointing)
    rows["r_angle"] = np.hypot(rows["y_angle"], rows["z_angle"])
    return rows


def cat_rows(obsid, ra, dec, pointing):
    rows = Table(np.zeros(len(ra), dtype=db.ASTROMON_CAT_SRC_DTYPE))
    rows["obsid"] = obsid
    rows["id"] = np.arange(len(ra))
    rows["celldetect_x_id"] = 1
    rows["catalog"] = "Tycho2"
    rows["ra"], rows["dec"] = ra, dec
    rows["y_angle"], rows["z_angle"] = angles(ra, dec, pointing)
    return rows


RA = [150.001, 150.010, 149.995]
DEC = [2.002, 1.990, 2.008]


# --- local_evt2_pointing ---------------------------------------------------------


def test_local_evt2_pointing_reads_the_events_header(tmp_path):
    write_evt2(
        primary_dir(tmp_path, 7001) / "acisf07001N002_evt2.fits.gz", EVT2_POINTING
    )

    assert fix.local_evt2_pointing(7001, [tmp_path]) == pytest.approx(EVT2_POINTING)


def test_local_evt2_pointing_prefers_the_archive_evt2_to_the_filtered_copy(tmp_path):
    primary = primary_dir(tmp_path, 7001)
    write_evt2(primary / "7001_evt2_filtered.fits.gz", OBSPAR_POINTING)
    write_evt2(primary / "acisf07001N002_evt2.fits.gz", EVT2_POINTING)

    assert fix.local_evt2_pointing(7001, [tmp_path]) == pytest.approx(EVT2_POINTING)


def test_local_evt2_pointing_falls_back_to_the_filtered_copy(tmp_path):
    """filter_events copies the evt2 header, so a cleaned-up work dir still has it."""
    write_evt2(
        primary_dir(tmp_path, 62877) / "62877_evt2_filtered.fits.gz", EVT2_POINTING
    )

    assert fix.local_evt2_pointing(62877, [tmp_path]) == pytest.approx(EVT2_POINTING)


def test_local_evt2_pointing_searches_roots_in_order(tmp_path):
    first, second = tmp_path / "work", tmp_path / "other"
    write_evt2(primary_dir(second, 7001) / "acisf07001N002_evt2.fits.gz", EVT2_POINTING)

    assert fix.local_evt2_pointing(7001, [first, second]) == pytest.approx(
        EVT2_POINTING
    )


def test_local_evt2_pointing_is_none_without_a_local_evt2(tmp_path):
    (primary_dir(tmp_path, 7001)).mkdir(parents=True)

    assert fix.local_evt2_pointing(7001, [tmp_path]) is None
    assert fix.local_evt2_pointing(7002, [tmp_path]) is None


# --- angle_errors / recompute_angles ---------------------------------------------


def test_angle_errors_is_zero_for_rows_consistent_with_the_pointing():
    rows = xray_rows(7001, "celldetect", RA, DEC, EVT2_POINTING)

    errors = fix.angle_errors(rows, {7001: EVT2_POINTING})

    assert np.all(errors < 1e-3)


def test_angle_errors_measures_rows_computed_with_another_pointing():
    rows = xray_rows(7001, "celldetect", RA, DEC, OBSPAR_POINTING)

    errors = fix.angle_errors(rows, {7001: EVT2_POINTING})

    assert np.all(errors > 10)


def test_angle_errors_is_nan_for_an_obsid_without_a_pointing():
    rows = xray_rows(7001, "celldetect", RA, DEC, OBSPAR_POINTING)

    assert np.all(np.isnan(fix.angle_errors(rows, {})))


def test_recompute_angles_uses_the_given_pointing_and_updates_r_angle():
    rows = xray_rows(7001, "celldetect", RA, DEC, OBSPAR_POINTING)

    fixed = fix.recompute_angles(rows, {7001: EVT2_POINTING})

    y_expected, z_expected = angles(RA, DEC, EVT2_POINTING)
    assert np.allclose(fixed["y_angle"], y_expected, atol=1e-3)
    assert np.allclose(fixed["z_angle"], z_expected, atol=1e-3)
    assert np.allclose(fixed["r_angle"], np.hypot(y_expected, z_expected), atol=1e-3)
    assert np.all(fix.angle_errors(fixed, {7001: EVT2_POINTING}) < 1e-3)


def test_recompute_angles_changes_nothing_but_the_angles():
    rows = xray_rows(7001, "celldetect", RA, DEC, OBSPAR_POINTING)

    fixed = fix.recompute_angles(rows, {7001: EVT2_POINTING})

    for name in rows.colnames:
        if name not in ("y_angle", "z_angle", "r_angle"):
            assert np.all(fixed[name] == rows[name]), name
    # and the input is left alone
    assert np.all(fix.angle_errors(rows, {7001: EVT2_POINTING}) > 10)


def test_recompute_angles_handles_cat_src_which_has_no_r_angle():
    rows = cat_rows(7001, RA, DEC, OBSPAR_POINTING)

    fixed = fix.recompute_angles(rows, {7001: EVT2_POINTING})

    assert "r_angle" not in fixed.colnames
    assert np.all(fix.angle_errors(fixed, {7001: EVT2_POINTING}) < 1e-3)


def test_recompute_angles_leaves_obsids_without_a_pointing_alone():
    rows = xray_rows(7001, "celldetect", RA, DEC, OBSPAR_POINTING)

    fixed = fix.recompute_angles(rows, {})

    assert np.all(fixed["y_angle"] == rows["y_angle"])
    assert np.all(fixed["z_angle"] == rows["z_angle"])


# --- stale_groups ----------------------------------------------------------------


def test_stale_groups_names_each_obsid_and_method_with_a_bad_row():
    rows = Table(
        np.concatenate(
            [
                xray_rows(7001, "celldetect", RA, DEC, OBSPAR_POINTING).as_array(),
                xray_rows(7001, "gaussian_detect", RA, DEC, EVT2_POINTING).as_array(),
                xray_rows(7002, "celldetect", RA, DEC, EVT2_POINTING).as_array(),
                xray_rows(7002, "gaussian_detect", RA, DEC, OBSPAR_POINTING).as_array(),
            ]
        )
    )
    pointings = {7001: EVT2_POINTING, 7002: EVT2_POINTING}

    groups = fix.stale_groups(rows, pointings, ["obsid", "detect_method"])

    assert groups == [(7001, "celldetect"), (7002, "gaussian_detect")]


def test_stale_groups_ignores_differences_below_the_tolerance():
    rows = xray_rows(7001, "celldetect", RA, DEC, EVT2_POINTING)
    rows["y_angle"] += 0.05

    assert fix.stale_groups(rows, {7001: EVT2_POINTING}, ["obsid"]) == []
    assert fix.stale_groups(rows, {7001: EVT2_POINTING}, ["obsid"], tolerance=0.01) == [
        (7001,)
    ]


def test_stale_groups_skips_obsids_without_a_pointing():
    rows = xray_rows(7001, "celldetect", RA, DEC, OBSPAR_POINTING)

    assert fix.stale_groups(rows, {}, ["obsid"]) == []


# --- repair (end to end on a small database) -------------------------------------


def _seed(tmp_path: Path) -> tuple[Path, Path]:
    """A database and a work dir.

    7001: celldetect stale (obspar angles), gaussian_detect fine, cat_src fine.
    7002: everything fine.
    7003: cat_src stale, and the only local evt2 is the filtered copy.
    7004: celldetect stale, but no local evt2 at all -- left alone.
    """
    dbfile, workdir = tmp_path / "astromon.h5", tmp_path / "work"
    db.create_empty_tables(dbfile)
    xray = [
        xray_rows(7001, "celldetect", RA, DEC, OBSPAR_POINTING),
        xray_rows(7001, "gaussian_detect", RA, DEC, EVT2_POINTING),
        xray_rows(7002, "celldetect", RA, DEC, EVT2_POINTING),
        xray_rows(7003, "celldetect", RA, DEC, EVT2_POINTING),
        xray_rows(7004, "celldetect", RA, DEC, OBSPAR_POINTING),
    ]
    db.save(
        "astromon_xray_src", Table(np.concatenate([t.as_array() for t in xray])), dbfile
    )
    cat = [
        cat_rows(7001, RA, DEC, EVT2_POINTING),
        cat_rows(7002, RA, DEC, EVT2_POINTING),
        cat_rows(7003, RA, DEC, OBSPAR_POINTING),
    ]
    db.save(
        "astromon_cat_src", Table(np.concatenate([t.as_array() for t in cat])), dbfile
    )
    for obsid in (7001, 7002):
        write_evt2(
            primary_dir(workdir, obsid) / f"acisf0{obsid}N002_evt2.fits.gz",
            EVT2_POINTING,
        )
    write_evt2(primary_dir(workdir, 7003) / "7003_evt2_filtered.fits.gz", EVT2_POINTING)
    return dbfile, workdir


def _sorted(table: Table, keys: list[str]) -> Table:
    table = table.copy()
    table.sort(keys)
    return table


def test_repair_fixes_exactly_the_stale_groups(tmp_path):
    dbfile, workdir = _seed(tmp_path)
    xray_before = _sorted(
        db.get_table("astromon_xray_src", dbfile), ["obsid", "detect_method", "id"]
    )
    cat_before = _sorted(db.get_table("astromon_cat_src", dbfile), ["obsid", "id"])

    summary = fix.repair(dbfile, [workdir])

    assert summary["astromon_xray_src"] == [(7001, "celldetect")]
    assert summary["astromon_cat_src"] == [(7003,)]
    assert summary["no_pointing"] == [7004]
    assert summary["touched_obsids"] == [7001, 7003]

    xray = _sorted(
        db.get_table("astromon_xray_src", dbfile), ["obsid", "detect_method", "id"]
    )
    cat = _sorted(db.get_table("astromon_cat_src", dbfile), ["obsid", "id"])
    assert len(xray) == len(xray_before) and len(cat) == len(cat_before)
    pointings = {7001: EVT2_POINTING, 7002: EVT2_POINTING, 7003: EVT2_POINTING}
    has_pointing = np.asarray(xray["obsid"]) != 7004
    assert np.all(fix.angle_errors(xray[has_pointing], pointings) < 1e-3)
    assert np.all(fix.angle_errors(cat, pointings) < 1e-3)

    # untouched groups are bit-for-bit what they were, 7004 included
    unchanged = ~(
        (np.asarray(xray["obsid"]) == 7001)
        & (np.asarray(xray["detect_method"]).astype(str) == "celldetect")
    )
    for name in xray.colnames:
        assert np.all(xray[name][unchanged] == xray_before[name][unchanged]), name
    for name in ("obsid", "id", "celldetect_x_id", "catalog", "ra", "dec"):
        assert np.all(cat[name] == cat_before[name]), name
    other_cat = np.asarray(cat["obsid"]) != 7003
    assert np.all(cat["y_angle"][other_cat] == cat_before["y_angle"][other_cat])


def test_repair_dry_run_writes_nothing(tmp_path):
    dbfile, workdir = _seed(tmp_path)
    before = dbfile.read_bytes()

    summary = fix.repair(dbfile, [workdir], dry_run=True)

    assert summary["touched_obsids"] == [7001, 7003]
    assert dbfile.read_bytes() == before


def test_repair_can_be_restricted_to_given_obsids(tmp_path):
    dbfile, workdir = _seed(tmp_path)

    summary = fix.repair(dbfile, [workdir], obsids=[7003, 7004])

    assert summary["astromon_xray_src"] == []
    assert summary["astromon_cat_src"] == [(7003,)]
    assert summary["no_pointing"] == [7004]
    xray = db.get_table("astromon_xray_src", dbfile)
    still_stale = (np.asarray(xray["obsid"]) == 7001) & (
        np.asarray(xray["detect_method"]).astype(str) == "celldetect"
    )
    assert np.all(fix.angle_errors(xray[still_stale], {7001: EVT2_POINTING}) > 10)


def test_repair_reports_rows_it_moves_past_180_arcsec(tmp_path):
    """get_sources drops sources at r_angle >= 180"; the repair reports, not drops, them."""
    dbfile, workdir = _seed(tmp_path)
    ra, dec = yagzag_to_radec(181.0, 0.0, Quat(equatorial=EVT2_POINTING))
    edge = xray_rows(7002, "gaussian_detect", [ra], [dec], EVT2_POINTING, first_id=9)
    edge["y_angle"], edge["z_angle"], edge["r_angle"] = 179.0, 0.0, 179.0
    db.save("astromon_xray_src", edge, dbfile)

    summary = fix.repair(dbfile, [workdir])

    assert summary["past_180"] == [(7002, "gaussian_detect", 9)]
    xray = db.get_table("astromon_xray_src", dbfile)
    kept = xray[(xray["obsid"] == 7002) & (xray["id"] == 9)]
    assert kept["r_angle"][0] == pytest.approx(181.0, abs=1e-3)


def test_repair_is_idempotent(tmp_path):
    dbfile, workdir = _seed(tmp_path)
    fix.repair(dbfile, [workdir])

    summary = fix.repair(dbfile, [workdir])

    assert summary["touched_obsids"] == []
