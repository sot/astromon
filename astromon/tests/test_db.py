import tempfile
from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

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
