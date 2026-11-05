"""Tests for recomputing celldetect grating_arm flags from the .src and the arm mask.

Context: for 13 grating observations the arm mask was never made (NaN event
positions crashed the zero-order finder), so _on_grating_arm flagged nothing and
their celldetect rows were stored with grating_arm False everywhere. Once the mask
exists, the flags can be recomputed from each stored source's sky position in the
.src, without re-deriving anything else.
"""

from pathlib import Path

import numpy as np
from astropy.table import Table

from astromon import db
from astromon.scripts.maintenance import backfill_grating_arm_flags as backfill

#: A zero order at (100, 100) and one dispersed arm running up from it.
ZERO_ORDER = (100.0, 100.0, 5.0)
ARM = (100.0, 300.0, 3.0, 150.0, 0.0)  # centre x, y, half width, half height, angle


def obsid_dir(root: Path, obsid: int) -> Path:
    return root / f"obs{obsid // 1000:02d}" / str(obsid)


def write_mask(root: Path, obsid: int) -> None:
    rows = [
        {
            "TG_PART": 0,
            "X": ZERO_ORDER[0],
            "Y": ZERO_ORDER[1],
            "R": [ZERO_ORDER[2], ZERO_ORDER[2]],
            "ROTANG": 0.0,
        },
        {
            "TG_PART": 1,
            "X": ARM[0],
            "Y": ARM[1],
            "R": [ARM[2], ARM[3]],
            "ROTANG": ARM[4],
        },
    ]
    path = obsid_dir(root, obsid) / "images" / f"{obsid}_grating_arm_mask.fits"
    path.parent.mkdir(parents=True, exist_ok=True)
    Table(rows=rows).write(path, format="fits")


def write_celldetect_src(root: Path, obsid: int, xy: list[tuple[float, float]]) -> None:
    """COMPONENT 1..n at the given sky pixel positions."""
    path = obsid_dir(root, obsid) / "sources" / f"{obsid}_celldetect.src"
    path.parent.mkdir(parents=True, exist_ok=True)
    Table(
        {
            "COMPONENT": np.arange(1, len(xy) + 1, dtype=np.int32),
            "X": np.array([x for x, _ in xy]),
            "Y": np.array([y for _, y in xy]),
        }
    ).write(path, format="fits")


def xray_rows(
    obsid: int, method: str, source_ids: list[int], on_arm: list[bool]
) -> Table:
    rows = Table(np.zeros(len(source_ids), dtype=db.ASTROMON_XRAY_SRC_DTYPE))
    rows["obsid"] = obsid
    rows["id"] = source_ids
    rows["detect_method"] = method
    rows["ra"] = 150.0 + np.arange(len(source_ids)) / 3600
    rows["snr"] = 9.0
    rows["grating_arm"] = on_arm
    return rows


#: obsid 7001: id 1 at the zero order, 2 on the arm, 3 off it, 4 on the arm but
#: beyond 180" (not stored). obsid 7002: a .src but no mask.
SOURCES_7001 = [(100.0, 100.0), (100.0, 300.0), (400.0, 400.0), (100.0, 420.0)]


def _seed(tmp_path: Path) -> tuple[Path, Path]:
    dbfile, work = tmp_path / "astromon.h5", tmp_path / "work"
    db.create_empty_tables(dbfile)
    rows = [
        xray_rows(7001, "celldetect", [1, 2, 3], [False, False, False]),
        xray_rows(7001, "gaussian_detect", [1, 2], [False, False]),
        xray_rows(7002, "celldetect", [1], [False]),
    ]
    db.save(
        "astromon_xray_src", Table(np.concatenate([r.as_array() for r in rows])), dbfile
    )
    write_mask(work, 7001)
    write_celldetect_src(work, 7001, SOURCES_7001)
    write_celldetect_src(work, 7002, [(100.0, 300.0)])
    return dbfile, work


def _flags(dbfile: Path, obsid: int, method: str = "celldetect") -> dict[int, bool]:
    xray = db.get_table("astromon_xray_src", dbfile)
    rows = xray[
        (xray["obsid"] == obsid)
        & (np.asarray(xray["detect_method"]).astype(str) == method)
    ]
    return {
        int(i): bool(g) for i, g in zip(rows["id"], rows["grating_arm"], strict=True)
    }


def test_arm_source_ids_reads_the_src_and_the_mask(tmp_path):
    _, work = _seed(tmp_path)

    assert backfill.arm_source_ids(7001, [work]) == {2, 4}


def test_arm_source_ids_is_none_without_a_mask(tmp_path):
    _, work = _seed(tmp_path)

    assert backfill.arm_source_ids(7002, [work]) is None
    assert backfill.arm_source_ids(7003, [work]) is None


def test_backfill_sets_the_flags_of_stored_celldetect_rows(tmp_path):
    dbfile, work = _seed(tmp_path)

    summary = backfill.backfill(dbfile, [7001, 7002], [work])

    assert _flags(dbfile, 7001) == {1: False, 2: True, 3: False}
    assert summary["updated"] == {7001: 1}
    assert summary["no_mask"] == [7002]


def test_backfill_changes_nothing_else(tmp_path):
    dbfile, work = _seed(tmp_path)
    before = db.get_table("astromon_xray_src", dbfile)

    backfill.backfill(dbfile, [7001, 7002], [work])

    after = db.get_table("astromon_xray_src", dbfile)
    for table in (before, after):
        table.sort(["obsid", "detect_method", "id"])
    changed = np.asarray(before["grating_arm"]) != np.asarray(after["grating_arm"])
    assert int(changed.sum()) == 1
    for name in before.colnames:
        if name != "grating_arm":
            assert np.all(before[name] == after[name]), name
    assert _flags(dbfile, 7001, "gaussian_detect") == {1: False, 2: False}


def test_backfill_dry_run_writes_nothing(tmp_path):
    dbfile, work = _seed(tmp_path)
    before = dbfile.read_bytes()

    summary = backfill.backfill(dbfile, [7001], [work], dry_run=True)

    assert summary["updated"] == {7001: 1}
    assert dbfile.read_bytes() == before


def test_backfill_is_idempotent(tmp_path):
    dbfile, work = _seed(tmp_path)
    backfill.backfill(dbfile, [7001], [work])

    summary = backfill.backfill(dbfile, [7001], [work])

    assert summary["updated"] == {}


def test_backfill_can_invalidate_the_cached_sources_it_corrected(tmp_path):
    """The sources cache still holds the old flags; a later read must recompute them."""
    dbfile, work = _seed(tmp_path)

    backfill.backfill(dbfile, [7001], [work], invalidate_caches=True)

    cache = obsid_dir(work, 7001) / "cache"
    markers = sorted(p.name for p in cache.glob(".sources__*.invalidated"))
    assert markers == [
        ".sources__62080401290e26adc129bc6513a2fdad.fits.gz.invalidated",
        ".sources__d9c8dab1a259c6355c91f8e9dd79041f.fits.gz.invalidated",
    ]
