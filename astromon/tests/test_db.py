import tempfile
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table, vstack

from astromon import db, utils

DATA_DIR = Path(__file__).parent / "data"


def test_get():
    table_names = [
        "astromon_cat_src",
        "astromon_meta",
        "astromon_obs",
        "astromon_regions",
        "astromon_xcorr",
        "astromon_xray_src",
    ]
    for name in table_names:
        db.get_table(name)


@pytest.mark.filterwarnings("error")
def test_dtypes():
    dbfile_ext = "h5"
    # this checks that the dtypes in db.DTYPES
    # and those of the returned tables match
    names = [
        "astromon_xcorr",
        "astromon_cat_src",
        "astromon_xray_src",
        "astromon_obs",
        "astromon_regions",
    ]
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / f"test_dtypes.{dbfile_ext}"
        for name in names:
            # t = db.get_table(name, dbfile)  # fails, file does not exist
            db.save(
                name, Table.read(DATA_DIR / f"{name}.ecsv"), dbfile
            )  # populate the table
            t = db.get_table(name, dbfile)
            assert t.dtype.names == db.DTYPES[name].names
            dtypes_differ = [
                (n, t.dtype[n], db.DTYPES[name][n], db.DTYPES[name][n].char)
                for n in t.dtype.names
                if t.dtype[n] != db.DTYPES[name][n]
                # skip comparison on HDF5 because of string/unicode differences
                and (dbfile_ext != "h5" or db.DTYPES[name][n].char != "S")
            ]
            assert not dtypes_differ, f"{name} dtypes differ: {dtypes_differ}"


@pytest.mark.filterwarnings("error")
def test_save_and_get():
    tables = {}
    tables_h5 = {}
    names = ["astromon_xcorr", "astromon_cat_src", "astromon_xray_src", "astromon_obs"]
    for name in names:
        tables[name] = Table.read(DATA_DIR / f"{name}.ecsv")

    with tempfile.TemporaryDirectory() as tmpdir:
        print("HDF5 format")
        dbfile = Path(tmpdir) / "test_save_and_read.h5"
        for name in names:
            print(name)
            db.save(name, tables[name], dbfile)
            tables_h5[name] = db.get_table(name, dbfile)

        for name in names:
            assert len(tables[name]) == len(tables_h5[name]), (
                f"{name} table has different lengths"
            )
            # string types are different
            for colname in tables[name].colnames:
                if tables[name][colname].dtype.char not in ["U", "S"]:
                    assert np.all(
                        (tables[name][colname] == tables_h5[name][colname])
                        | (
                            np.isnan(tables[name][colname])
                            & np.isnan(tables_h5[name][colname])
                        )
                    ), f"{name} {colname} column differs"


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

        # warnings are filtered as errors here
        with pytest.raises(utils.MissingTableException):
            db.get_table("astromon_regions", dbfile)
        # adding an empty table to prevent exception
        db.save("astromon_regions", db.create_table("astromon_regions"), dbfile)

        # adding one at a time
        db.get_table("astromon_regions", dbfile)
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
        assert regions.colnames == regions_ref.colnames, "col names"
        dtypes_differ = [
            (n, regions.dtype[n], db.DTYPES["astromon_regions"][n])
            for n in regions.dtype.names
            if regions.dtype[n] != db.DTYPES["astromon_regions"][n]
            # skip comparison on HDF5 because of string/unicode differences
            and (dbfile_ext != "h5" or db.DTYPES["astromon_regions"][n].char != "S")
        ]
        assert not dtypes_differ, f" dtypes differ: {dtypes_differ}"
        for name in regions.colnames:
            if db.DTYPES["astromon_regions"][name].char != "S":
                assert np.all(regions[name] == regions_ref[name])

        # removing
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
        assert regions.colnames == regions_ref.colnames, "col names"
        dtypes_differ = [
            (n, regions.dtype[n], db.DTYPES["astromon_regions"][n])
            for n in regions.dtype.names
            if regions.dtype[n] != db.DTYPES["astromon_regions"][n]
            # skip comparison on HDF5 because of string/unicode differences
            and (dbfile_ext != "h5" or db.DTYPES["astromon_regions"][n].char != "S")
        ]
        assert not dtypes_differ, f" dtypes differ: {dtypes_differ}"
        for name in regions.colnames:
            if db.DTYPES["astromon_regions"][name].char != "S":
                assert np.all(regions[name] == regions_ref[name])

        # removing so it is empty
        db.remove_regions([2], dbfile=dbfile)
        regions = db.get_table("astromon_regions", dbfile=dbfile)
        assert regions.colnames == regions_ref.colnames, "col names"
        # assert regions.dtype == regions_ref.dtype, 'dtypes'
        dtypes_differ = [
            (n, regions.dtype[n], db.DTYPES["astromon_regions"][n])
            for n in regions.dtype.names
            if regions.dtype[n] != db.DTYPES["astromon_regions"][n]
            # skip comparison on HDF5 because of string/unicode differences
            and (dbfile_ext != "h5" or db.DTYPES["astromon_regions"][n].char != "S")
        ]
        assert not dtypes_differ, f" dtypes differ: {dtypes_differ}"
        assert len(regions) == 0

        # remove non-existent
        db.remove_regions([3], dbfile=dbfile)  # silently removes nothing

        # adding a few and with autoincrementing region_id
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
        assert regions.colnames == regions_ref.colnames, "col names"
        for name in regions.colnames:
            if db.DTYPES["astromon_regions"][name].char != "S":
                assert np.all(regions[name] == regions_ref[name])


def setup_files_for_region_sync(tmpdir):
    # these are the columns we will check in this test
    cols = [
        "region_id",
        "region_id_str",
        "ra",
        "dec",
        "radius",
        "obsid",
        "user",
        "comments",
    ]

    # we start with this table, make some independent changes in two separate db files
    # and then sync them in different ways to verify the results are as expected
    date = "2025:317:14:47:27.097"
    basis = [
        (9, "8f53da3773", date, 128.836, -45.176, 200.0, 0, "j", "0"),
        (11, "40fc4e84fb", date, 187.702, 12.392, 5.0, 0, "jgonzalez", "0"),
        (12, "63b66ed331", date, 187.701, 12.392, 5.0, 0, "jgonzalez", "0"),
        (13, "a1c54bdab6", date, 187.704, 12.391, 5.0, 0, "jgonzalez", "0"),
    ]

    meta = Table(np.array([(13,)], dtype=db.ASTROMON_META_DTYPE))
    basis = Table(np.array(basis, dtype=db.ASTROMON_REGION_DTYPE))
    db.set_formats(basis)

    dbfile_1 = Path(tmpdir) / "test_region_sync_1.h5"
    db.save("astromon_meta", meta, dbfile_1)
    db.save("astromon_regions", basis, dbfile_1)

    dbfile_2 = Path(tmpdir) / "test_region_sync_2.h5"
    db.save("astromon_meta", meta, dbfile_2)
    db.save("astromon_regions", basis, dbfile_2)

    # make changes in dbfile_1
    add_1 = Table(
        [
            {
                "ra": 0.0,
                "dec": 0.0,
                "radius": 5,
                "obsid": 0,
                "user": "me",
                "comments": "1",
            },
        ]
    )
    db.add_regions(add_1, dbfile=dbfile_1)
    db.remove_regions([13], dbfile=dbfile_1)
    regions_1 = db.get_table("astromon_regions", dbfile=dbfile_1)
    regions_1["ra"][:2] = [128.8, 187.6]
    db.update_regions(regions_1, dbfile=dbfile_1)

    print("\n# sanity check after modifications dbfile_1")
    regions_1 = db.get_table("astromon_regions", dbfile=dbfile_1)
    ref_1 = Table(
        np.array(
            [
                (9, "8f53da3773", 128.8, -45.176, 200.0, 0, "j", "0"),
                (11, "40fc4e84fb", 187.6, 12.392, 5.0, 0, "jgonzalez", "0"),
                (12, "63b66ed331", 187.701, 12.392, 5.0, 0, "jgonzalez", "0"),
                # this id string will be different
                (14, "", 0.0, 0.0, 5.0, 0, "me", "1"),
            ],
            dtype=db.ASTROMON_REGION_DTYPE[cols],
        )
    )
    db.set_formats(ref_1)
    # all id strings should be the same except for the region we just added
    assert np.all(regions_1["region_id_str"][:-1] == ref_1["region_id_str"][:-1]), (
        "ID after modifications dbfile_1"
    )
    ref_1["region_id_str"] = regions_1["region_id_str"]
    print("ref")
    ref_1.pprint(max_width=-1, max_lines=-1)
    print("test")
    regions_1[cols].pprint(max_width=-1, max_lines=-1)
    for col in cols:
        assert np.all(regions_1[col] == ref_1[col]), (
            f"{col} after modifications dbfile_1"
        )

    # make changes in dbfile_2
    add_2 = Table(
        [
            {
                "ra": 1.0,
                "dec": 0.0,
                "radius": 5,
                "obsid": 0,
                "user": "them",
                "comments": "2",
            }
        ]
    )
    db.add_regions(add_2, dbfile=dbfile_2)
    db.remove_regions([12], dbfile=dbfile_2)
    regions_2 = db.get_table("astromon_regions", dbfile=dbfile_2)
    regions_2["ra"][1] = 187.7
    db.update_regions(regions_2, dbfile=dbfile_2)

    # print("\n# sanity check after modifications dbfile_2")
    regions_2 = db.get_table("astromon_regions", dbfile=dbfile_2)
    ref_2 = Table(
        np.array(
            [
                (9, "8f53da3773", 128.836, -45.176, 200.0, 0, "j", "0"),
                (11, "40fc4e84fb", 187.7, 12.392, 5.0, 0, "jgonzalez", "0"),
                (13, "a1c54bdab6", 187.704, 12.391, 5.0, 0, "jgonzalez", "0"),
                # this id string will be different
                (14, "", 1.0, 0.0, 5.0, 0, "them", "2"),
            ],
            dtype=db.ASTROMON_REGION_DTYPE[cols],
        )
    )
    db.set_formats(ref_2)
    # all id strings should be the same except for the region we just added
    assert np.all(regions_2["region_id_str"][:-1] == ref_2["region_id_str"][:-1]), (
        "ID after modifications dbfile_2"
    )
    ref_2["region_id_str"] = regions_2["region_id_str"]
    # print("ref")
    # ref_2.pprint(max_width=-1, max_lines=-1)
    # print("test")
    # regions_2[cols].pprint(max_width=-1, max_lines=-1)
    for col in cols:
        assert np.all(regions_2[col] == ref_2[col]), (
            f"{col} after modifications dbfile_2"
        )

    return dbfile_1, dbfile_2


def test_region_sync_1_2():
    # these are the columns we will check in this test
    cols = [
        "region_id",
        "region_id_str",
        "ra",
        "dec",
        "radius",
        "obsid",
        "user",
        "comments",
    ]

    # in this test:
    # - region 9 is synced (because it has a more recent last_modified date in dbfile_1)
    # - region 11 is not synced (because it has a more recent last_modified date in dbfile_2)
    # - region 14 from dbfile_1 is added and assigned region_id 15 to avoid repeated region_id

    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile_1, dbfile_2 = setup_files_for_region_sync(tmpdir)

        print("\n# syncing 1 -> 2")

        regions_1 = db.get_table("astromon_regions", dbfile=dbfile_1)
        regions_2 = db.get_table("astromon_regions", dbfile=dbfile_2)

        db.sync_regions(dbfile_1, dbfile_2)

        result = db.get_table("astromon_regions", dbfile=dbfile_2)

        id1 = regions_2["region_id_str"][-1]
        id2 = regions_1["region_id_str"][-1]
        ref_3 = Table(
            np.array(
                [
                    (9, "8f53da3773", 128.8, -45.176, 200.0, 0, "j", "0"),
                    (11, "40fc4e84fb", 187.7, 12.392, 5.0, 0, "jgonzalez", "0"),
                    (12, "63b66ed331", 187.701, 12.392, 5.0, 0, "jgonzalez", "0"),
                    (13, "a1c54bdab6", 187.704, 12.391, 5.0, 0, "jgonzalez", "0"),
                    (14, id1, 1.0, 0.0, 5.0, 0, "them", "2"),
                    (15, id2, 0.0, 0.0, 5.0, 0, "me", "1"),
                ],
                dtype=db.ASTROMON_REGION_DTYPE[cols],
            )
        )
        db.set_formats(ref_3)

        print("regions_1")
        regions_1.pprint()
        print("regions_2")
        regions_2.pprint()
        print("ref")
        ref_3.pprint(max_width=-1, max_lines=-1)
        print("test")
        result[cols].pprint(max_width=-1, max_lines=-1)
        for col in cols:
            assert np.all(result[col] == ref_3[col]), f"{col} after rsync 1 -> 3"


def test_region_sync_1_2_rm():
    # these are the columns we will check in this test
    cols = [
        "region_id",
        "region_id_str",
        "ra",
        "dec",
        "radius",
        "obsid",
        "user",
        "comments",
    ]

    # in this test:
    # - region 9 is synced (because it has a more recent last_modified date in dbfile_1)
    # - region 11 is not synced (because it has a more recent last_modified date in dbfile_2)
    # - region 13 is removed from dbfile_2
    # - region 14 is removed from dbfile_2
    # - region 14 from dbfile_1 is added

    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile_1, dbfile_2 = setup_files_for_region_sync(tmpdir)

        print("\n# syncing 1 -> 2 with removal")

        regions_1 = db.get_table("astromon_regions", dbfile=dbfile_1)
        regions_2 = db.get_table("astromon_regions", dbfile=dbfile_2)

        db.sync_regions(dbfile_1, dbfile_2, remove=True)

        result = db.get_table("astromon_regions", dbfile=dbfile_2)

        id2 = regions_1["region_id_str"][-1]
        ref_3 = Table(
            np.array(
                [
                    (9, "8f53da3773", 128.8, -45.176, 200.0, 0, "j", "0"),
                    (11, "40fc4e84fb", 187.7, 12.392, 5.0, 0, "jgonzalez", "0"),
                    (12, "63b66ed331", 187.701, 12.392, 5.0, 0, "jgonzalez", "0"),
                    (14, id2, 0.0, 0.0, 5.0, 0, "me", "1"),
                ],
                dtype=db.ASTROMON_REGION_DTYPE[cols],
            )
        )
        db.set_formats(ref_3)

        print("regions_1")
        regions_1.pprint()
        print("regions_2")
        regions_2.pprint()
        print("ref")
        ref_3.pprint(max_width=-1, max_lines=-1)
        print("test")
        result[cols].pprint(max_width=-1, max_lines=-1)
        for col in cols:
            assert np.all(result[col] == ref_3[col]), f"{col} after rsync 1 -> 3"


def test_region_sync_2_1():
    # these are the columns we will check in this test
    cols = [
        "region_id",
        "region_id_str",
        "ra",
        "dec",
        "radius",
        "obsid",
        "user",
        "comments",
    ]

    # in this test:
    # - region 9 is not synced (because it has a more recent last_modified date in dbfile_1)
    # - region 11 is synced (because it has a more recent last_modified date in dbfile_2)
    # - region 14 from dbfile_2 is added and assigned region_id 15 to avoid repeated region_id

    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile_1, dbfile_2 = setup_files_for_region_sync(tmpdir)

        print("\n# syncing 2 -> 1")

        regions_1 = db.get_table("astromon_regions", dbfile=dbfile_1)
        regions_2 = db.get_table("astromon_regions", dbfile=dbfile_2)

        db.sync_regions(dbfile_2, dbfile_1)

        result = db.get_table("astromon_regions", dbfile=dbfile_1)

        id1 = regions_2["region_id_str"][-1]
        id2 = regions_1["region_id_str"][-1]
        ref_4 = Table(
            np.array(
                [
                    (9, "8f53da3773", 128.8, -45.176, 200.0, 0, "j", "0"),
                    (11, "40fc4e84fb", 187.7, 12.392, 5.0, 0, "jgonzalez", "0"),
                    (12, "63b66ed331", 187.701, 12.392, 5.0, 0, "jgonzalez", "0"),
                    (13, "a1c54bdab6", 187.704, 12.391, 5.0, 0, "jgonzalez", "0"),
                    (14, id2, 0.0, 0.0, 5.0, 0, "me", "1"),
                    (15, id1, 1.0, 0.0, 5.0, 0, "them", "2"),
                ],
                dtype=db.ASTROMON_REGION_DTYPE[cols],
            )
        )
        db.set_formats(ref_4)

        print("regions_1")
        regions_1.pprint()
        print("regions_2")
        regions_2.pprint()
        print("ref")
        ref_4.pprint(max_width=-1, max_lines=-1)
        print("test")
        result[cols].pprint(max_width=-1, max_lines=-1)
        for col in cols:
            assert np.all(result[col] == ref_4[col]), f"{col} after rsync 1 -> 3"


def test_region_sync_2_1_rm():
    # these are the columns we will check in this test
    cols = [
        "region_id",
        "region_id_str",
        "ra",
        "dec",
        "radius",
        "obsid",
        "user",
        "comments",
    ]

    # in this test:
    # - region 9 is not synced (because it has a more recent last_modified date in dbfile_1)
    # - region 11 is synced (because it has a more recent last_modified date in dbfile_2)
    # - region 14 from dbfile_2 is added to dbfile_1
    # - region 12 is removed from dbfile_1
    # - region 14 is removed from dbfile_1

    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile_1, dbfile_2 = setup_files_for_region_sync(tmpdir)

        print("\n# syncing 2 -> 1 with removal")

        regions_1 = db.get_table("astromon_regions", dbfile=dbfile_1)
        regions_2 = db.get_table("astromon_regions", dbfile=dbfile_2)

        db.sync_regions(dbfile_2, dbfile_1, remove=True)

        result = db.get_table("astromon_regions", dbfile=dbfile_1)

        id1 = regions_2["region_id_str"][-1]
        ref_4 = Table(
            np.array(
                [
                    (9, "8f53da3773", 128.8, -45.176, 200.0, 0, "j", "0"),
                    (11, "40fc4e84fb", 187.7, 12.392, 5.0, 0, "jgonzalez", "0"),
                    (13, "a1c54bdab6", 187.704, 12.391, 5.0, 0, "jgonzalez", "0"),
                    (14, id1, 1.0, 0.0, 5.0, 0, "them", "2"),
                ],
                dtype=db.ASTROMON_REGION_DTYPE[cols],
            )
        )

        db.set_formats(ref_4)

        print("regions_1")
        regions_1.pprint()
        print("regions_2")
        regions_2.pprint()
        print("ref")
        ref_4.pprint(max_width=-1, max_lines=-1)
        print("test")
        result[cols].pprint(max_width=-1, max_lines=-1)
        for col in cols:
            assert np.all(result[col] == ref_4[col]), f"{col} after rsync 1 -> 3"


# ---------------------------------------------------------------------------
# _cast_to_dtype tests
# ---------------------------------------------------------------------------


def test_cast_to_dtype_identity_when_dtypes_match():
    """Returns the input array unchanged when the dtype already matches."""
    dtype = np.dtype([("obsid", np.int32), ("snr", np.float32)])
    arr = np.zeros(3, dtype=dtype)
    arr["obsid"] = [10, 20, 30]
    arr["snr"] = [5.0, 8.0, 3.0]

    result = db._cast_to_dtype(arr, dtype)

    assert result is arr  # same object, not a copy


def test_cast_to_dtype_zero_fills_missing_field():
    """New fields not present in the source array are zero-filled."""
    old_dtype = np.dtype([("obsid", np.int32), ("snr", np.float32)])
    new_dtype = np.dtype(
        [("obsid", np.int32), ("snr", np.float32), ("brightest", np.int32)]
    )

    arr = np.zeros(4, dtype=old_dtype)
    arr["obsid"] = [1, 2, 3, 4]
    arr["snr"] = [5.0, 10.0, 3.0, 7.0]

    result = db._cast_to_dtype(arr, new_dtype)

    assert result.dtype == new_dtype
    np.testing.assert_array_equal(result["obsid"], [1, 2, 3, 4])
    np.testing.assert_array_almost_equal(result["snr"], [5.0, 10.0, 3.0, 7.0])
    np.testing.assert_array_equal(result["brightest"], [0, 0, 0, 0])


def test_cast_to_dtype_preserves_existing_values():
    """Columns present in both old and new dtype keep their values."""
    old_dtype = np.dtype(
        [("obsid", np.int32), ("snr", np.float32), ("acis_streak", np.int32)]
    )
    new_dtype = np.dtype(
        [
            ("obsid", np.int32),
            ("snr", np.float32),
            ("acis_streak", np.int32),
            ("brightest", np.int32),
        ]
    )

    arr = np.zeros(3, dtype=old_dtype)
    arr["obsid"] = [10, 20, 30]
    arr["snr"] = [2.0, 9.0, 4.0]
    arr["acis_streak"] = [0, 1, 0]

    result = db._cast_to_dtype(arr, new_dtype)

    np.testing.assert_array_equal(result["obsid"], [10, 20, 30])
    np.testing.assert_array_almost_equal(result["snr"], [2.0, 9.0, 4.0])
    np.testing.assert_array_equal(result["acis_streak"], [0, 1, 0])
    np.testing.assert_array_equal(result["brightest"], [0, 0, 0])


def test_save_schema_migration_zero_fills_new_column():
    """save() correctly migrates an existing HDF5 table that lacks a new column.

    Simulates the real scenario: DB written before 'brightest' was added,
    then save() is called with data conforming to the new dtype. Existing rows
    must have the new column zero-filled, while the newly saved rows carry
    the real values.
    """
    import tables as tb

    # Old dtype: astromon_xray_src minus 'brightest'
    old_dtype = np.dtype(
        [
            (n, db.ASTROMON_XRAY_SRC_DTYPE[n])
            for n in db.ASTROMON_XRAY_SRC_DTYPE.names
            if n != "brightest"
        ]
    )
    new_dtype = db.ASTROMON_XRAY_SRC_DTYPE

    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "migration_test.h5"

        # Write a row for obsid=1111 directly with the OLD schema.
        old_row = np.zeros(1, dtype=old_dtype)
        old_row["obsid"] = 1111
        old_row["id"] = 1
        old_row["snr"] = 12.0
        old_row["acis_streak"] = 1
        with tb.open_file(str(dbfile), "a") as h5:
            h5.create_table("/", "astromon_xray_src", old_dtype, "xray src")
            node = h5.root.astromon_xray_src
            node.append(old_row)

        # Now save new data for obsid=2222 using the current (new) dtype.
        new_row = Table(np.zeros(1, dtype=new_dtype))
        new_row["obsid"] = 2222
        new_row["id"] = 1
        new_row["snr"] = 8.0
        new_row["brightest"] = 1
        db.save("astromon_xray_src", new_row, dbfile)

        result = db.get_table("astromon_xray_src", dbfile)

    # Both obsids present.
    assert set(result["obsid"]) == {1111, 2222}

    # Non-string column names must match exactly (HDF5 returns Unicode strings
    # for byte-string columns, so skip those — same caveat as test_dtypes above).
    non_string_names = [n for n in new_dtype.names if new_dtype[n].char != "S"]
    assert result.dtype.names == new_dtype.names, "column names and order must match"
    mismatched = [
        (n, result.dtype[n], new_dtype[n])
        for n in non_string_names
        if result.dtype[n] != new_dtype[n]
    ]
    assert not mismatched, f"non-string dtypes differ: {mismatched}"

    old_row_result = result[result["obsid"] == 1111]
    new_row_result = result[result["obsid"] == 2222]

    # Old row: 'brightest' was not in the original table, must be zero-filled.
    assert int(old_row_result["brightest"][0]) == 0
    assert int(old_row_result["acis_streak"][0]) == 1

    # New row: 'brightest' was explicitly set to 1.
    assert int(new_row_result["brightest"][0]) == 1


def _xcorr_row(
    *,
    select_name: str,
    obsid: int = 7001,
    detect_method: str = "celldetect",
    c_id: int = 10,
    x_id: int = 20,
    dy: float = 0.1,
    dz: float = 0.2,
    dr: float = 0.3,
) -> Table:
    """Build a one-row astromon_xcorr table using the current DB dtype."""
    row = Table(np.zeros(1, dtype=db.ASTROMON_XCORR_DTYPE))
    row["select_name"] = select_name
    row["obsid"] = obsid
    row["c_id"] = c_id
    row["x_id"] = x_id
    row["dy"] = dy
    row["dz"] = dz
    row["dr"] = dr
    row["detect_method"] = detect_method
    return row


def test_save_xcorr_select_name_key_preserves_other_select_names():
    """select_name_key=True only replaces the matching select_name rows.

    When saving a partial astromon_xcorr backfill (e.g. gaia_agn only), rows for
    other select_names with the same obsid+detect_method must remain untouched.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "select_name_key_test.h5"

        original = vstack(
            [
                _xcorr_row(select_name="astromon_21", dy=1.1, dz=1.2, dr=1.3),
                _xcorr_row(select_name="gaia_agn", dy=2.1, dz=2.2, dr=2.3),
            ],
            metadata_conflicts="silent",
        )
        db.save("astromon_xcorr", original, dbfile)

        updated = _xcorr_row(select_name="gaia_agn", dy=9.1, dz=9.2, dr=9.3)
        db.save("astromon_xcorr", updated, dbfile, select_name_key=True)

        result = db.get_table("astromon_xcorr", dbfile)

    assert len(result) == 2

    astromon_21 = result[result["select_name"] == "astromon_21"]
    gaia_agn = result[result["select_name"] == "gaia_agn"]

    assert len(astromon_21) == 1
    assert len(gaia_agn) == 1
    np.testing.assert_allclose(astromon_21["dy"], [1.1])
    np.testing.assert_allclose(astromon_21["dz"], [1.2])
    np.testing.assert_allclose(astromon_21["dr"], [1.3])
    np.testing.assert_allclose(gaia_agn["dy"], [9.1])
    np.testing.assert_allclose(gaia_agn["dz"], [9.2])
    np.testing.assert_allclose(gaia_agn["dr"], [9.3])


def test_save_xcorr_without_select_name_key_replaces_all_select_names_for_method():
    """Default save() behavior still replaces all rows for obsid+detect_method."""
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "detect_method_key_test.h5"

        original = vstack(
            [
                _xcorr_row(select_name="astromon_21", dy=1.1),
                _xcorr_row(select_name="gaia_agn", dy=2.2),
            ],
            metadata_conflicts="silent",
        )
        db.save("astromon_xcorr", original, dbfile)

        replacement = _xcorr_row(select_name="gaia_agn", dy=7.7)
        db.save("astromon_xcorr", replacement, dbfile)

        result = db.get_table("astromon_xcorr", dbfile)

    assert len(result) == 1
    assert result["select_name"].tolist() == ["gaia_agn"]
    np.testing.assert_allclose(result["dy"], [7.7])


def _stored_dtype(dtype: np.dtype, extra_fields=()) -> np.dtype:
    """On-disk (pytables) form of a DB dtype: unicode columns become bytes."""
    fields = [
        (
            name,
            f"S{dtype[name].itemsize // 4}" if dtype[name].kind == "U" else dtype[name],
        )
        for name in dtype.names
    ]
    return np.dtype(fields + list(extra_fields))


def _write_raw_db(dbfile, xray_extra_fields=(), xray_extra_values=None):
    """Create an astromon HDF5 file directly, bypassing db.save's dtype cast.

    This is how a database written by an older code version looks: it can carry
    columns that the current schema no longer has.
    """
    import tables

    xray_dtype = _stored_dtype(db.ASTROMON_XRAY_SRC_DTYPE, xray_extra_fields)
    xray = np.zeros(3, dtype=xray_dtype)
    xray["obsid"] = [1, 2, 3]
    xray["ra"] = [10.0, 20.0, 30.0]
    for name, values in (xray_extra_values or {}).items():
        xray[name] = values

    xcorr = np.zeros(1, dtype=_stored_dtype(db.ASTROMON_XCORR_DTYPE))

    with tables.open_file(str(dbfile), "w") as h5:
        h5.create_table("/", "astromon_xray_src", xray)
        h5.create_table("/", "astromon_xcorr", xcorr)


def test_migrate_db_drops_retired_sb_radius():
    """migrate() removes the retired sb_radius column and preserves other data."""
    from astromon.scripts import migrate_db

    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "old_schema.h5"
        _write_raw_db(
            dbfile,
            xray_extra_fields=[("sb_radius", "<f4")],
            xray_extra_values={"sb_radius": [7.0, 8.0, 9.0]},
        )

        migrate_db.migrate(dbfile)

        result = db.get_table("astromon_xray_src", dbfile)
        assert "sb_radius" not in result.colnames
        np.testing.assert_allclose(sorted(result["ra"]), [10.0, 20.0, 30.0])
        # A backup of the pre-migration file is left next to the original.
        assert dbfile.with_suffix(".h5.pre_migrate").exists()


def test_migrate_db_dry_run_leaves_file_unchanged():
    """migrate(dry_run=True) reports but does not modify the database."""
    from astromon.scripts import migrate_db

    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "old_schema.h5"
        _write_raw_db(
            dbfile,
            xray_extra_fields=[("sb_radius", "<f4")],
            xray_extra_values={"sb_radius": [7.0, 8.0, 9.0]},
        )
        before = dbfile.read_bytes()

        migrate_db.migrate(dbfile, dry_run=True)

        assert dbfile.read_bytes() == before
        assert not dbfile.with_suffix(".h5.pre_migrate").exists()


def test_migrate_db_refuses_unknown_extra_column():
    """A stored column not in the schema or the retired list raises instead of dropping."""
    from astromon.scripts import migrate_db

    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "unknown_extra.h5"
        _write_raw_db(dbfile, xray_extra_fields=[("mystery_column", "<f8")])

        with pytest.raises(ValueError, match="mystery_column"):
            migrate_db.migrate(dbfile)


def test_migrate_db_up_to_date_is_noop():
    """A database already at the current schema is left untouched."""
    from astromon.scripts import migrate_db

    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "current_schema.h5"
        _write_raw_db(dbfile)
        before = dbfile.read_bytes()

        migrate_db.migrate(dbfile)

        assert dbfile.read_bytes() == before


def _cat_src_row(
    *,
    catalog: str,
    obsid: int = 7001,
    x_id: int = 1,
    name: str = "src",
) -> Table:
    """Build a one-row astromon_cat_src table using the current DB dtype."""
    row = Table(np.zeros(1, dtype=db.ASTROMON_CAT_SRC_DTYPE))
    row["obsid"] = obsid
    row["x_id"] = x_id
    row["catalog"] = catalog
    row["name"] = name
    return row


def _seeded_cat_src_db(dbfile: Path) -> None:
    """astromon_cat_src with Tycho2 and RFC rows across two obsids."""
    original = vstack(
        [
            _cat_src_row(catalog="Tycho2", obsid=7001, name="tycho-a"),
            _cat_src_row(catalog="RFC", obsid=7001, name="rfc-a"),
            _cat_src_row(catalog="Tycho2", obsid=7002, name="tycho-b"),
            _cat_src_row(catalog="RFC", obsid=7002, name="rfc-b"),
        ],
        metadata_conflicts="silent",
    )
    db.save("astromon_cat_src", original, dbfile)


def test_save_replace_keys_by_catalog():
    """replace_keys=("catalog",) replaces all rows of the incoming catalogs, all obsids.

    This is the backfill pattern: rebuild one catalog's cat_src rows from scratch
    without touching other catalogs (previously done by hand with a full-table
    strip + vstack + ignore_obsid=True rewrite).
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "replace_by_catalog.h5"
        _seeded_cat_src_db(dbfile)

        new_rfc = _cat_src_row(catalog="RFC", obsid=7001, name="rfc-new")
        db.save("astromon_cat_src", new_rfc, dbfile, replace_keys=("catalog",))

        result = db.get_table("astromon_cat_src", dbfile)

    tycho = result[result["catalog"] == "Tycho2"]
    rfc = result[result["catalog"] == "RFC"]
    assert sorted(tycho["name"]) == ["tycho-a", "tycho-b"]
    # Both old RFC rows (including obsid 7002's) are gone; only the new row remains.
    assert list(rfc["name"]) == ["rfc-new"]


def test_save_replace_keys_by_obsid_and_catalog():
    """replace_keys=("obsid", "catalog") only replaces the pairs present in the data."""
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "replace_by_obsid_catalog.h5"
        _seeded_cat_src_db(dbfile)

        new_rfc = _cat_src_row(catalog="RFC", obsid=7001, name="rfc-new")
        db.save("astromon_cat_src", new_rfc, dbfile, replace_keys=("obsid", "catalog"))

        result = db.get_table("astromon_cat_src", dbfile)

    rfc = result[result["catalog"] == "RFC"]
    # obsid 7001's RFC row is replaced; obsid 7002's RFC row survives.
    assert sorted(rfc["name"]) == ["rfc-b", "rfc-new"]
    assert len(result[result["catalog"] == "Tycho2"]) == 2


def test_save_replace_keys_matches_select_name_key():
    """replace_keys spelling gives the same result as the legacy select_name_key flag."""
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile_flag = Path(tmpdir) / "with_flag.h5"
        dbfile_keys = Path(tmpdir) / "with_keys.h5"

        original = vstack(
            [
                _xcorr_row(select_name="astromon_21", dy=1.1),
                _xcorr_row(select_name="gaia_agn", dy=2.2),
            ],
            metadata_conflicts="silent",
        )
        updated = _xcorr_row(select_name="gaia_agn", dy=9.9)

        db.save("astromon_xcorr", original, dbfile_flag)
        db.save("astromon_xcorr", updated, dbfile_flag, select_name_key=True)

        db.save("astromon_xcorr", original, dbfile_keys)
        db.save(
            "astromon_xcorr",
            updated,
            dbfile_keys,
            replace_keys=("obsid", "detect_method", "select_name"),
        )

        result_flag = db.get_table("astromon_xcorr", dbfile_flag)
        result_keys = db.get_table("astromon_xcorr", dbfile_keys)

    result_flag.sort("select_name")
    result_keys.sort("select_name")
    assert result_flag.colnames == result_keys.colnames
    for col in result_flag.colnames:
        assert np.array_equal(result_flag[col], result_keys[col]), col


def test_save_replace_keys_unknown_column_raises():
    """A replace_keys column absent from the table schema raises ValueError."""
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "bad_key.h5"
        _seeded_cat_src_db(dbfile)

        with pytest.raises(ValueError, match="not_a_column"):
            db.save(
                "astromon_cat_src",
                _cat_src_row(catalog="RFC"),
                dbfile,
                replace_keys=("not_a_column",),
            )


def test_save_replace_keys_conflicting_flags_raise():
    """replace_keys cannot be combined with ignore_obsid or select_name_key."""
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "conflict.h5"
        _seeded_cat_src_db(dbfile)
        row = _cat_src_row(catalog="RFC")

        with pytest.raises(ValueError, match="ignore_obsid"):
            db.save(
                "astromon_cat_src",
                row,
                dbfile,
                ignore_obsid=True,
                replace_keys=("catalog",),
            )
        with pytest.raises(ValueError, match="select_name_key"):
            db.save(
                "astromon_cat_src",
                row,
                dbfile,
                select_name_key=True,
                replace_keys=("catalog",),
            )


def test_create_empty_tables_creates_all_known_tables():
    """create_empty_tables() populates a fresh file with every schema, zero rows."""
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "fresh.h5"
        db.create_empty_tables(dbfile)

        for name, dtype in db.DTYPES.items():
            result = db.get_table(name, dbfile)
            assert len(result) == 0, name
            assert result.dtype.names == dtype.names, name


def test_create_empty_tables_preserves_existing_rows():
    """create_empty_tables() is idempotent: it never clobbers a populated table."""
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "populated.h5"
        _seeded_cat_src_db(dbfile)

        db.create_empty_tables(dbfile)

        result = db.get_table("astromon_cat_src", dbfile)
        assert len(result) == 4
        assert sorted(result["name"]) == ["rfc-a", "rfc-b", "tycho-a", "tycho-b"]


def test_save_expect_existing_raises_when_table_absent():
    """expect_existing=True refuses to create a table that should already be there."""
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "no_table.h5"
        _seeded_cat_src_db(dbfile)  # creates astromon_cat_src only

        with pytest.raises(utils.MissingTableException, match="astromon_xcorr"):
            db.save(
                "astromon_xcorr",
                _xcorr_row(select_name="astromon_21"),
                dbfile,
                expect_existing=True,
            )


def test_save_expect_existing_allows_normal_append():
    """expect_existing=True is a no-op when the table is present."""
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "append.h5"
        _seeded_cat_src_db(dbfile)

        db.save(
            "astromon_cat_src",
            _cat_src_row(catalog="Quaia", obsid=7003, name="quaia-c"),
            dbfile,
            expect_existing=True,
        )

        result = db.get_table("astromon_cat_src", dbfile)
        assert len(result) == 5
        assert "quaia-c" in list(result["name"])


def test_save_expect_existing_catches_destroyed_table():
    """A table destroyed mid-write is reported loudly instead of silently recreated.

    db.save writes by remove_node + create_table, which is not atomic: a process
    killed between the two (run_all SIGKILLs workers on timeout) leaves the table
    missing. Without the guard the next save recreates it holding only its own
    rows, silently discarding every other obsid.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "destroyed.h5"
        _seeded_cat_src_db(dbfile)
        assert len(db.get_table("astromon_cat_src", dbfile)) == 4

        # Simulate the crash window: node removed, replacement never written.
        with db.connect(dbfile, mode="r+") as h5:
            h5.remove_node(h5.get_node("/astromon_cat_src"))

        with pytest.raises(utils.MissingTableException, match="astromon_cat_src"):
            db.save(
                "astromon_cat_src",
                _cat_src_row(catalog="RFC", obsid=7001, name="rfc-new"),
                dbfile,
                expect_existing=True,
            )


def test_cast_to_dtype_fills_missing_float_with_nan():
    """Float columns absent from the stored data read back as NaN, not 0.0.

    0.0 is a meaningful value for these columns (psfratio, concentration_ratio),
    so zero-filling an unmeasured column makes it indistinguishable from a real
    measurement and silently admits those rows into threshold filters. NaN is
    the honest "not measured". Integer columns keep 0 -- there is no integer NaN,
    and 0 is the correct conservative default for the boolean-ish flags.
    """
    old_dtype = np.dtype([("obsid", np.int32), ("snr", np.float32)])
    new_dtype = np.dtype(
        [
            ("obsid", np.int32),
            ("snr", np.float32),
            ("psfratio", np.float32),
            ("brightest", np.int32),
        ]
    )

    arr = np.zeros(2, dtype=old_dtype)
    arr["obsid"] = [1, 2]
    arr["snr"] = [5.0, 10.0]

    result = db._cast_to_dtype(arr, new_dtype)

    assert np.all(np.isnan(result["psfratio"])), "unmeasured float must be NaN"
    np.testing.assert_array_equal(result["brightest"], [0, 0])
    np.testing.assert_array_almost_equal(result["snr"], [5.0, 10.0])


def test_get_table_does_not_invent_zero_for_absent_float_column():
    """Reading a pre-schema-change file must not fabricate 0.0 measurements."""
    import tables as tb

    old_dtype = np.dtype(
        [
            (n, db.ASTROMON_XRAY_SRC_DTYPE[n])
            for n in db.ASTROMON_XRAY_SRC_DTYPE.names
            if n != "psfratio"
        ]
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "legacy_schema.h5"
        rows = np.zeros(3, dtype=old_dtype)
        rows["obsid"] = [1, 2, 3]
        with tb.open_file(str(dbfile), "a") as h5:
            h5.create_table("/", "astromon_xray_src", old_dtype, "xray src")
            h5.root.astromon_xray_src.append(rows)

        result = db.get_table("astromon_xray_src", dbfile)

    assert np.all(np.isnan(np.asarray(result["psfratio"], dtype=float)))


# --- region_id allocation ---------------------------------------------------
#
# add_regions takes its next id from astromon_meta.last_region_id and never
# consults the regions table it has already loaded. When the counter is absent
# or stale, numbering silently restarts and collides with existing rows --
# which matters because remove_regions() deletes by region_id, so one duplicate
# id makes a removal take out two regions.


def _seed_regions(dbfile: Path, ids=(1, 2, 3), last_region_id=None) -> None:
    """astromon_regions with the given ids, and astromon_meta only if asked.

    last_region_id=None leaves astromon_meta out of the file entirely, which is
    the state a database is in after create_empty_tables() has never run and
    also after a hand-built file.
    """
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
    """The baseline: a counter in step with the data allocates the next id."""
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "regions.h5"
        _seed_regions(dbfile, ids=(1, 2, 3), last_region_id=3)
        new_id, ids, counter = _add_and_report(dbfile)
        assert new_id == 4
        assert ids == [1, 2, 3, 4]
        assert counter == 4


def test_add_regions_with_absent_meta_does_not_duplicate_ids():
    """A missing counter means "unknown", not "no regions were ever allocated".

    Reproduces how $SKA/data/astromon/astromon.h5 ended up with ids 1,2,3,1.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "regions.h5"
        _seed_regions(dbfile, ids=(1, 2, 3), last_region_id=None)
        new_id, ids, counter = _add_and_report(dbfile)
        assert len(ids) == len(set(ids)), f"duplicate region_id allocated: {ids}"
        assert new_id == 4
        assert counter == 4


def test_add_regions_with_empty_meta_table_does_not_duplicate_ids():
    """create_empty_tables() makes astromon_meta with zero rows, not one.

    That hits the same len(...) > 0 false branch as an absent table, so a
    freshly initialized database is exposed too.
    """
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
    """A counter behind the data must not re-issue ids that are already in use."""
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "regions.h5"
        _seed_regions(dbfile, ids=(1, 2, 3), last_region_id=1)
        new_id, ids, counter = _add_and_report(dbfile)
        assert len(ids) == len(set(ids)), f"duplicate region_id allocated: {ids}"
        assert new_id == 4
        assert counter == 4, "the counter should be repaired, not left behind"


def test_add_regions_counter_never_goes_backwards():
    """A counter ahead of the data is authoritative: ids must not be reused.

    sync_regions can leave last_region_id above max(region_id) after rows are
    removed, and those retired ids must not come back.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "regions.h5"
        _seed_regions(dbfile, ids=(1, 2), last_region_id=9)
        new_id, ids, counter = _add_and_report(dbfile)
        assert new_id == 10, f"expected to continue from the counter, got {new_id}"
        assert len(ids) == len(set(ids))
        assert counter == 10


def test_duplicate_region_id_makes_removal_delete_two_regions():
    """Why the above matters: remove_regions deletes by id, so a duplicate
    id means one removal silently takes out an unrelated exclusion region."""
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / "regions.h5"
        _seed_regions(dbfile, ids=(1, 2, 3, 1), last_region_id=3)
        db.remove_regions([1], dbfile=dbfile)
        remaining = [
            int(v) for v in db.get_table("astromon_regions", dbfile)["region_id"]
        ]
        assert remaining == [2, 3], (
            "both rows sharing region_id=1 are removed together -- unique ids are"
            " what keeps a removal scoped to one region"
        )


# --- sync_regions id uniqueness --------------------------------------------
#
# sync_regions matches rows across files on region_id_str (a stable hash) but
# still carries an integer region_id, and remove_regions() deletes by that
# integer. So the merged table has to come out with unique ids no matter what
# the two inputs looked like.


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
    assert ids == sorted(set(ids)), f"ids must be unique and sorted: {ids}"
    assert last == max(ids)


def test_sync_regions_repairs_a_duplicate_carried_from_primary():
    """A primary file that already holds a duplicate id must not propagate it.

    $SKA/data/astromon/astromon.h5 holds region_ids 1,2,3,1, so this is the
    actual promotion case.
    """
    ids, strs, last = _synced(
        [1, 2, 3, 1], ["aa", "bb", "cc", "dd"], [], [], primary_last=1
    )
    assert len(strs) == 4, "all four regions should arrive"
    assert len(ids) == len(set(ids)), f"duplicate region_id after sync: {ids}"
    assert last == max(ids)


def test_sync_regions_does_not_invent_new_duplicates():
    """Renumbering must consider the whole merged set, not just the destination.

    Reassigning only the rows that collide with the destination, from a counter
    that ignores the other incoming ids, lands on ids the incoming rows already
    use -- turning one duplicate into several.
    """
    ids, strs, last = _synced(
        [1, 2, 3, 1], ["aa", "bb", "cc", "dd"], [1], ["zz"], dev_last=1, primary_last=1
    )
    assert len(strs) == 5, "one existing plus four new"
    assert len(ids) == len(set(ids)), f"duplicate region_id after sync: {ids}"
    assert last == max(ids)


def test_sync_regions_keeps_destination_ids_stable():
    """Rows already in dev must keep their ids; only incoming collisions move."""
    ids, strs, last = _synced(
        [7, 8], ["pp", "qq"], [7], ["zz"], dev_last=7, primary_last=8
    )
    by_str = dict(zip(strs, ids, strict=True))
    assert by_str["zz"] == 7, "the pre-existing dev row should not be renumbered"
    assert len(ids) == len(set(ids)), f"duplicate region_id after sync: {ids}"


def test_sync_regions_is_idempotent():
    """Running it twice changes nothing the second time."""
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
