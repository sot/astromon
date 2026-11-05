"""Tests for the per-obsid audit of X-ray y/z against catalog y/z and the sky.

Each astromon_cat_src row names its nearest celldetect source (celldetect_x_id).
If both sides' y/z come from the same pointing, the y/z distance between them is
the sky separation. When one side's y/z were computed with another pointing, the
y/z distance picks up that pointing difference, so median(dr_yz - sky_sep) over an
obsid's rows measures it.
"""

import numpy as np
import pytest
from astropy.table import Table, vstack
from chandra_aca.transform import radec_to_yagzag, yagzag_to_radec
from Quaternion import Quat

from astromon import db
from astromon.scripts.maintenance import audit_yz_consistency as audit

EVT2_POINTING = (150.0, 2.0, 30.0)
#: ~25" from EVT2_POINTING, like the quick-look obspar in the affected obsids
OBSPAR_POINTING = (150.0060, 2.0035, 30.004)

#: (yag, zag) of three X-ray sources, arcsec
XRAY_YZ = [(10.0, -20.0), (-60.0, 35.0), (90.0, 80.0)]
#: each catalog counterpart's offset from its X-ray source, arcsec
COUNTERPART_OFFSET = (0.4, -0.3)


def xray_rows(
    obsid: int,
    pointing_for_angles: tuple[float, float, float],
    method: str = "celldetect",
) -> Table:
    """Celldetect sources at XRAY_YZ (under EVT2_POINTING), angles stored with `pointing_for_angles`."""
    quat = Quat(equatorial=EVT2_POINTING)
    yag, zag = np.array(XRAY_YZ).T
    ra, dec = yagzag_to_radec(yag, zag, quat)
    rows = Table(np.zeros(len(yag), dtype=db.ASTROMON_XRAY_SRC_DTYPE))
    rows["obsid"] = obsid
    rows["id"] = np.arange(1, len(yag) + 1)
    rows["detect_method"] = method
    rows["ra"], rows["dec"] = ra, dec
    rows["y_angle"], rows["z_angle"] = radec_to_yagzag(
        ra, dec, Quat(equatorial=pointing_for_angles)
    )
    rows["r_angle"] = np.hypot(rows["y_angle"], rows["z_angle"])
    return rows


def cat_rows(obsid: int, anchors: tuple[int, int, int] = (1, 2, 3)) -> Table:
    """One catalog counterpart per X-ray source, angles from EVT2_POINTING."""
    quat = Quat(equatorial=EVT2_POINTING)
    yag, zag = np.array(XRAY_YZ).T
    ra, dec = yagzag_to_radec(
        yag + COUNTERPART_OFFSET[0], zag + COUNTERPART_OFFSET[1], quat
    )
    rows = Table(np.zeros(len(yag), dtype=db.ASTROMON_CAT_SRC_DTYPE))
    rows["obsid"] = obsid
    rows["id"] = np.arange(len(yag))
    rows["celldetect_x_id"] = anchors
    rows["catalog"] = "Tycho2"
    rows["ra"], rows["dec"] = ra, dec
    rows["y_angle"], rows["z_angle"] = radec_to_yagzag(ra, dec, quat)
    return rows


def test_anchor_pairs_measure_y_z_distance_and_sky_separation():
    pairs = audit.anchor_pairs(xray_rows(7001, EVT2_POINTING), cat_rows(7001))

    assert len(pairs) == 3
    offset = np.hypot(*COUNTERPART_OFFSET)
    assert np.allclose(pairs["dr_yz"], offset, atol=1e-3)
    assert np.allclose(pairs["sky_sep"], offset, atol=1e-3)


def test_excess_is_zero_when_both_sides_use_one_pointing():
    excess = audit.yz_excess_by_obsid(xray_rows(7001, EVT2_POINTING), cat_rows(7001))

    assert list(excess["obsid"]) == [7001]
    assert excess["n_pairs"][0] == 3
    assert excess["excess"][0] == pytest.approx(0.0, abs=1e-3)


def test_excess_measures_x_ray_angles_from_another_pointing():
    excess = audit.yz_excess_by_obsid(xray_rows(7001, OBSPAR_POINTING), cat_rows(7001))

    assert excess["excess"][0] > 15


def test_only_celldetect_sources_are_anchors():
    """celldetect_x_id is a celldetect id; a gaussian source with the same id is not its anchor."""
    xray = vstack(
        [
            xray_rows(7001, EVT2_POINTING),
            xray_rows(7001, OBSPAR_POINTING, method="gaussian_detect"),
        ]
    )

    excess = audit.yz_excess_by_obsid(xray, cat_rows(7001))

    assert excess["n_pairs"][0] == 3
    assert excess["excess"][0] == pytest.approx(0.0, abs=1e-3)


def test_rows_without_an_anchor_are_left_out():
    """celldetect ids start at 1, so celldetect_x_id 0 (or an id that is gone) anchors nothing."""
    excess = audit.yz_excess_by_obsid(
        xray_rows(7001, OBSPAR_POINTING), cat_rows(7001, anchors=(0, 2, 9))
    )

    assert excess["n_pairs"][0] == 1


def test_obsids_with_no_anchored_rows_are_not_listed():
    xray = vstack([xray_rows(7001, EVT2_POINTING), xray_rows(7002, OBSPAR_POINTING)])

    excess = audit.yz_excess_by_obsid(xray, cat_rows(7001))

    assert list(excess["obsid"]) == [7001]


def test_inconsistent_obsids_uses_the_threshold():
    xray = vstack([xray_rows(7001, EVT2_POINTING), xray_rows(7002, OBSPAR_POINTING)])
    cat = vstack([cat_rows(7001), cat_rows(7002)])
    excess = audit.yz_excess_by_obsid(xray, cat)

    assert audit.inconsistent_obsids(excess) == [7002]
    assert audit.inconsistent_obsids(excess, threshold=100.0) == []


def test_audit_reads_a_database_and_writes_nothing_to_it(tmp_path):
    dbfile = tmp_path / "astromon.h5"
    db.create_empty_tables(dbfile)
    db.save(
        "astromon_xray_src",
        vstack([xray_rows(7001, EVT2_POINTING), xray_rows(7002, OBSPAR_POINTING)]),
        dbfile,
    )
    db.save("astromon_cat_src", vstack([cat_rows(7001), cat_rows(7002)]), dbfile)
    before = dbfile.read_bytes()

    excess = audit.audit(dbfile)

    assert audit.inconsistent_obsids(excess) == [7002]
    assert dbfile.read_bytes() == before


def test_no_anchored_rows_gives_an_empty_audit():
    """A database with nothing anchored (e.g. a new one) has nothing to flag."""
    no_xray = Table(np.zeros(0, dtype=db.ASTROMON_XRAY_SRC_DTYPE))
    no_cat = Table(np.zeros(0, dtype=db.ASTROMON_CAT_SRC_DTYPE))

    excess = audit.yz_excess_by_obsid(no_xray, no_cat)

    assert len(excess) == 0
    assert audit.inconsistent_obsids(excess) == []
