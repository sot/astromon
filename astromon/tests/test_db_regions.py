"""Tests for astromon region-table behaviors.

These checks are all about region identity and synchronization: allocating new
ids, keeping them unique, and merging copies of the same region table without
reusing or duplicating ids. Keeping them separate from the generic DB tests
makes the region-specific invariants much easier to review.
"""

import tempfile
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

from astromon import db, utils

DATA_DIR = Path(__file__).parent / "data"


@pytest.mark.filterwarnings("error")
def test_regions():
    dbfile_ext = "h5"
    tables = {}
    names = ["astromon_xcorr", "astromon_cat_src", "astromon_xray_src", "astromon_obs"]

    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / f"test_regions.{dbfile_ext}"
        for name in names:
            tables[name] = Table.read(DATA_DIR / f"{name}.ecsv")
            db.save(name, tables[name], dbfile)

        with pytest.raises(utils.MissingTableException):
            db.get_table("astromon_regions", dbfile)
        db.save("astromon_regions", db.create_table("astromon_regions"), dbfile)

        regions_1 = Table(
            [
                {
                    "ra": 0.0,
                    "dec": 0.0,
                    "radius": 5,
                    "obsid": 0,
                    "user": "me",
                    "comments": "",
                }
            ]
        )
        db.add_regions(regions_1, dbfile=dbfile)
        regions_2 = Table(
            [
                {
                    "ra": 1.0,
                    "dec": 0.0,
                    "radius": 5,
                    "obsid": 0,
                    "user": "them",
                    "comments": "",
                }
            ]
        )
        db.add_regions(regions_2, dbfile=str(dbfile))
        regions = db.get_table("astromon_regions", dbfile=dbfile)
        regions_ref = Table(
            np.array(
                [
                    (
                        1,
                        "97dffe105c",
                        "2025:317:14:50:54.677",
                        0.0,
                        0.0,
                        5.0,
                        0,
                        "me",
                        "",
                    ),
                    (
                        2,
                        "6f4940d9ee",
                        "2025:317:14:50:54.691",
                        1.0,
                        0.0,
                        5.0,
                        0,
                        "them",
                        "",
                    ),
                ],
                dtype=db.ASTROMON_REGION_DTYPE,
            )
        )
        assert regions.colnames == regions_ref.colnames
        dtypes_differ = [
            (n, regions.dtype[n], db.DTYPES["astromon_regions"][n])
            for n in regions.dtype.names
            if regions.dtype[n] != db.DTYPES["astromon_regions"][n]
            and (dbfile_ext != "h5" or db.DTYPES["astromon_regions"][n].char != "S")
        ]
        assert not dtypes_differ, f" dtypes differ: {dtypes_differ}"
        for name in regions.colnames:
            if db.DTYPES["astromon_regions"][name].char != "S":
                assert np.all(regions[name] == regions_ref[name])

        db.remove_regions([1], dbfile=dbfile)
        regions = db.get_table("astromon_regions", dbfile=dbfile)
        regions_ref = Table(
            np.array(
                [
                    (
                        2,
                        "6f4940d9ee",
                        "2025:317:14:50:54.691",
                        1.0,
                        0.0,
                        5.0,
                        0,
                        "them",
                        "",
                    )
                ],
                dtype=db.ASTROMON_REGION_DTYPE,
            )
        )
        assert regions.colnames == regions_ref.colnames
        dtypes_differ = [
            (n, regions.dtype[n], db.DTYPES["astromon_regions"][n])
            for n in regions.dtype.names
            if regions.dtype[n] != db.DTYPES["astromon_regions"][n]
            and (dbfile_ext != "h5" or db.DTYPES["astromon_regions"][n].char != "S")
        ]
        assert not dtypes_differ, f" dtypes differ: {dtypes_differ}"
        for name in regions.colnames:
            if db.DTYPES["astromon_regions"][name].char != "S":
                assert np.all(regions[name] == regions_ref[name])

        db.remove_regions([2], dbfile=dbfile)
        regions = db.get_table("astromon_regions", dbfile=dbfile)
        assert regions.colnames == regions_ref.colnames
        dtypes_differ = [
            (n, regions.dtype[n], db.DTYPES["astromon_regions"][n])
            for n in regions.dtype.names
            if regions.dtype[n] != db.DTYPES["astromon_regions"][n]
            and (dbfile_ext != "h5" or db.DTYPES["astromon_regions"][n].char != "S")
        ]
        assert not dtypes_differ, f" dtypes differ: {dtypes_differ}"
        assert len(regions) == 0

        db.remove_regions([3], dbfile=dbfile)

        regions_1 = Table(
            [
                {
                    "ra": 0.0,
                    "dec": 0.0,
                    "radius": 5,
                    "obsid": 0,
                    "user": "me",
                    "comments": "",
                },
                {
                    "ra": 1.0,
                    "dec": 0.0,
                    "radius": 5,
                    "obsid": 0,
                    "user": "them",
                    "comments": "",
                },
            ]
        )
        db.add_regions(regions_1, dbfile=dbfile)
        regions = db.get_table("astromon_regions", dbfile=dbfile)
        regions_ref = Table(
            np.array(
                [
                    (
                        3,
                        "ef17c1872d",
                        "2025:317:14:59:18.139",
                        0.0,
                        0.0,
                        5.0,
                        0,
                        "me",
                        "",
                    ),
                    (
                        4,
                        "3c688d7cda",
                        "2025:317:14:59:18.139",
                        1.0,
                        0.0,
                        5.0,
                        0,
                        "them",
                        "",
                    ),
                ],
                dtype=db.ASTROMON_REGION_DTYPE,
            )
        )
        assert regions.colnames == regions_ref.colnames
        for name in regions.colnames:
            if db.DTYPES["astromon_regions"][name].char != "S":
                assert np.all(regions[name] == regions_ref[name])


def _seed_regions(dbfile: Path, ids=(1, 2, 3), last_region_id=None) -> None:
    rows = Table(np.zeros(len(ids), dtype=db.ASTROMON_REGION_DTYPE))
    rows["region_id"] = list(ids)
    rows["region_id_str"] = [f"seed{i:06d}" for i in ids]
    rows["ra"] = [10.0 * i for i in ids]
    rows["dec"] = [0.0] * len(ids)
    rows["radius"] = [10.0] * len(ids)
    db.save("astromon_regions", rows, dbfile, ignore_obsid=True)
    if last_region_id is not None:
        meta = Table(np.zeros(1, dtype=db.DTYPES["astromon_meta"]))
        meta["last_region_id"] = last_region_id
        db.save("astromon_meta", meta, dbfile, ignore_obsid=True)


def _one_region() -> Table:
    return Table(
        [
            {
                "ra": 99.0,
                "dec": 1.0,
                "radius": 10.0,
                "obsid": 12345,
                "user": "test",
                "comments": "new region",
            }
        ]
    )


def _add_and_report(dbfile: Path):
    new = db.add_regions(_one_region(), dbfile=dbfile)
    ids = [int(v) for v in db.get_table("astromon_regions", dbfile)["region_id"]]
    counter = int(db.get_table("astromon_meta", dbfile)["last_region_id"][0])
    return int(new["region_id"][0]), ids, counter


def test_add_regions_with_correct_counter_continues_the_sequence():
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "regions.h5"
        _seed_regions(dbfile, ids=(1, 2, 3), last_region_id=3)
        new_id, ids, counter = _add_and_report(dbfile)
        assert new_id == 4
        assert ids == [1, 2, 3, 4]
        assert counter == 4


def test_add_regions_with_absent_meta_does_not_duplicate_ids():
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "regions.h5"
        _seed_regions(dbfile, ids=(1, 2, 3), last_region_id=None)
        new_id, ids, counter = _add_and_report(dbfile)
        assert len(ids) == len(set(ids)), f"duplicate region_id allocated: {ids}"
        assert new_id == 4
        assert counter == 4


def test_add_regions_with_empty_meta_table_does_not_duplicate_ids():
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "regions.h5"
        _seed_regions(dbfile, ids=(1, 2, 3), last_region_id=None)
        db.save(
            "astromon_meta",
            db.create_table("astromon_meta"),
            dbfile,
            ignore_obsid=True,
        )
        assert len(db.get_table("astromon_meta", dbfile)) == 0
        new_id, ids, counter = _add_and_report(dbfile)
        assert len(ids) == len(set(ids)), f"duplicate region_id allocated: {ids}"
        assert new_id == 4


def test_add_regions_with_stale_counter_does_not_duplicate_ids():
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "regions.h5"
        _seed_regions(dbfile, ids=(1, 2, 3), last_region_id=1)
        new_id, ids, counter = _add_and_report(dbfile)
        assert len(ids) == len(set(ids)), f"duplicate region_id allocated: {ids}"
        assert new_id == 4
        assert counter == 4


def test_add_regions_counter_never_goes_backwards():
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "regions.h5"
        _seed_regions(dbfile, ids=(1, 2), last_region_id=9)
        new_id, ids, counter = _add_and_report(dbfile)
        assert new_id == 10
        assert len(ids) == len(set(ids))
        assert counter == 10


def test_duplicate_region_id_makes_removal_delete_two_regions():
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "regions.h5"
        _seed_regions(dbfile, ids=(1, 2, 3, 1), last_region_id=3)
        db.remove_regions([1], dbfile=dbfile)
        remaining = [
            int(v) for v in db.get_table("astromon_regions", dbfile)["region_id"]
        ]
        assert remaining == [2, 3]


def _regions_file(path: Path, ids, strs, last_region_id=0) -> None:
    rows = Table(np.zeros(len(ids), dtype=db.ASTROMON_REGION_DTYPE))
    if len(ids):
        rows["region_id"] = list(ids)
        rows["region_id_str"] = list(strs)
        rows["ra"] = [10.0 * i for i in range(len(ids))]
        rows["radius"] = [10.0] * len(ids)
        rows["last_modified"] = ["2025:001:00:00:00.000"] * len(ids)
    db.save("astromon_regions", rows, path, ignore_obsid=True)
    meta = Table(np.zeros(1, dtype=db.DTYPES["astromon_meta"]))
    meta["last_region_id"] = last_region_id
    db.save("astromon_meta", meta, path, ignore_obsid=True)


def _synced(primary_ids, primary_strs, dev_ids, dev_strs, dev_last=0, primary_last=0):
    with tempfile.TemporaryDirectory() as tmpdir:
        primary = Path(tmpdir) / "primary.h5"
        dev = Path(tmpdir) / "dev.h5"
        _regions_file(primary, primary_ids, primary_strs, primary_last)
        _regions_file(dev, dev_ids, dev_strs, dev_last)
        db.sync_regions(primary, dev)
        regions = db.get_table("astromon_regions", dev)
        ids = [int(v) for v in regions["region_id"]]
        strs = [str(v) for v in regions["region_id_str"]]
        last = int(db.get_table("astromon_meta", dev)["last_region_id"][0])
        return ids, strs, last


def test_sync_regions_brings_over_new_rows():
    ids, strs, last = _synced([1, 2, 3], ["aa", "bb", "cc"], [], [], primary_last=3)
    assert sorted(strs) == ["aa", "bb", "cc"]
    assert ids == sorted(set(ids))
    assert last == max(ids)


def test_sync_regions_repairs_a_duplicate_carried_from_primary():
    ids, strs, last = _synced(
        [1, 2, 3, 1], ["aa", "bb", "cc", "dd"], [], [], primary_last=1
    )
    assert len(strs) == 4
    assert len(ids) == len(set(ids)), f"duplicate region_id after sync: {ids}"
    assert last == max(ids)


def test_sync_regions_does_not_invent_new_duplicates():
    ids, strs, last = _synced(
        [1, 2, 3, 1], ["aa", "bb", "cc", "dd"], [1], ["zz"], dev_last=1, primary_last=1
    )
    assert len(strs) == 5
    assert len(ids) == len(set(ids)), f"duplicate region_id after sync: {ids}"
    assert last == max(ids)


def test_sync_regions_keeps_destination_ids_stable():
    ids, strs, last = _synced(
        [7, 8], ["pp", "qq"], [7], ["zz"], dev_last=7, primary_last=8
    )
    by_str = dict(zip(strs, ids, strict=True))
    assert by_str["zz"] == 7
    assert len(ids) == len(set(ids)), f"duplicate region_id after sync: {ids}"
    assert last == max(ids)


def test_sync_regions_is_idempotent():
    with tempfile.TemporaryDirectory() as tmpdir:
        primary = Path(tmpdir) / "primary.h5"
        dev = Path(tmpdir) / "dev.h5"
        _regions_file(primary, [1, 2, 3, 1], ["aa", "bb", "cc", "dd"], 1)
        _regions_file(dev, [], [], 0)
        db.sync_regions(primary, dev)
        first = db.get_table("astromon_regions", dev)
        db.sync_regions(primary, dev)
        second = db.get_table("astromon_regions", dev)
        assert len(first) == len(second)
        assert sorted(str(v) for v in first["region_id_str"]) == sorted(
            str(v) for v in second["region_id_str"]
        )
        ids = [int(v) for v in second["region_id"]]
        assert len(ids) == len(set(ids))
