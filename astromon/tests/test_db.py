from pathlib import Path
import tempfile
import numpy as np
import pytest

from astropy.table import Table

from astromon import db
from astromon import utils


DATA_DIR = Path(__file__).parent / 'data'


@pytest.mark.filterwarnings("error")
@pytest.mark.parametrize('dbfile_ext', ['h5', 'db3'])
def test_dtypes(dbfile_ext):
    # this checks that the dtypes in db.DTYPES, the ones from astromon/sql/tables/*sql,
    # and those of the returned tables match
    names = [
        'astromon_xcorr', 'astromon_cat_src', 'astromon_xray_src', 'astromon_obs',
        'astromon_regions'
    ]
    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / f'test_dtypes.{dbfile_ext}'
        for name in names:
            # t = db.get_table(name, dbfile)  # fails, file does not exist
            db.save(name, Table.read(DATA_DIR / f'{name}.ecsv'), dbfile)  # populate the table
            t = db.get_table(name, dbfile)
            assert t.dtype.names == db.DTYPES[name].names
            dtypes_differ = [
                (n, t.dtype[n], db.DTYPES[name][n], db.DTYPES[name][n].char)
                for n in t.dtype.names if t.dtype[n] != db.DTYPES[name][n]
                # skip comparison on HDF5 because of string/unicode differences
                and (dbfile_ext != 'h5' or db.DTYPES[name][n].char != 'S')
            ]
            assert not dtypes_differ, f'{name} dtypes differ: {dtypes_differ}'


@pytest.mark.filterwarnings("error")
def test_save_and_get():
    tables = {}
    tables_h5 = {}
    tables_sql = {}
    names = ['astromon_xcorr', 'astromon_cat_src', 'astromon_xray_src', 'astromon_obs']
    for name in names:
        tables[name] = Table.read(DATA_DIR / f'{name}.ecsv')

    with tempfile.TemporaryDirectory() as tmpdir:
        print('HDF5 format')
        dbfile = Path(tmpdir) / 'test_save_and_read.h5'
        for name in names:
            print(name)
            db.save(name, tables[name], dbfile)
            tables_h5[name] = db.get_table(name, dbfile)

        print('SQLite format')
        dbfile = Path(tmpdir) / 'test_save_and_read.db3'
        for name in names:
            db.save(name, tables[name], dbfile)
            tables_sql[name] = db.get_table(name, dbfile)

        for name in names:
            assert len(tables_sql[name]) == len(tables_h5[name]), \
                f'{name} table has different lengths'
            # string types are different
            for colname in tables_sql[name].colnames:
                if tables_sql[name][colname].dtype.char != 'S':
                    assert (np.all(
                        (tables_sql[name][colname] == tables_h5[name][colname])
                        | (np.isnan(tables_sql[name][colname]) & np.isnan(tables_h5[name][colname]))
                    )), f'{name} {colname} column differs'


@pytest.mark.filterwarnings("error")
@pytest.mark.parametrize('dbfile_ext', ['h5', 'db3'])
def test_regions(dbfile_ext):
    tables = {}
    names = ['astromon_xcorr', 'astromon_cat_src', 'astromon_xray_src', 'astromon_obs']

    with tempfile.TemporaryDirectory() as tmpdir:
        dbfile = Path(tmpdir) / f'test_regions.{dbfile_ext}'
        for name in names:
            tables[name] = Table.read(DATA_DIR / f'{name}.ecsv')
            db.save(name, tables[name], dbfile)

        # warnings are filtered as errors here
        with pytest.raises(utils.MissingTableException):
            db.get_table('astromon_regions', dbfile)
        # adding an empty table to prevent exception
        db.save(
            'astromon_regions',
            db.create_table('astromon_regions'),
            dbfile
        )

        # adding one at a time
        db.get_table('astromon_regions', dbfile)
        regions_1 = Table([
            {'ra': 0., 'dec': 0., 'radius': 5, 'obsid': 0, 'user': 'me', 'comments': ''}
        ])
        db.add_regions(regions_1, dbfile=dbfile)
        regions_2 = Table(
            [{'ra': 1., 'dec': 0., 'radius': 5, 'obsid': 0, 'user': 'them', 'comments': ''}])
        db.add_regions(regions_2, dbfile=str(dbfile))
        regions = db.get_table('astromon_regions', dbfile=dbfile)
        regions_ref = Table(
            [[1, 2], [0.0, 1.0], [0.0, 0.0], [5, 5], [0, 0], ['me', 'them'], ['', '']],
            dtype=db.ASTROMON_REGION_DTYPE,
            names=[name for name in db.ASTROMON_REGION_DTYPE.names]
        )
        assert regions.colnames == regions_ref.colnames, 'col names'
        dtypes_differ = [
            (n, regions.dtype[n], db.DTYPES['astromon_regions'][n])
            for n in regions.dtype.names if regions.dtype[n] != db.DTYPES['astromon_regions'][n]
            # skip comparison on HDF5 because of string/unicode differences
            and (dbfile_ext != 'h5' or db.DTYPES['astromon_regions'][n].char != 'S')
        ]
        assert not dtypes_differ, f' dtypes differ: {dtypes_differ}'
        for name in regions.colnames:
            if db.DTYPES['astromon_regions'][name].char != 'S':
                assert np.all(regions[name] == regions_ref[name])

        # removing
        db.remove_regions([1], dbfile=dbfile)
        regions = db.get_table('astromon_regions', dbfile=dbfile)
        regions_ref = Table(
            [[2], [1.0], [0.0], [5], [0], ['them'], ['']],
            dtype=db.ASTROMON_REGION_DTYPE,
            names=[name for name in db.ASTROMON_REGION_DTYPE.names]
        )
        assert regions.colnames == regions_ref.colnames, 'col names'
        dtypes_differ = [
            (n, regions.dtype[n], db.DTYPES['astromon_regions'][n])
            for n in regions.dtype.names if regions.dtype[n] != db.DTYPES['astromon_regions'][n]
            # skip comparison on HDF5 because of string/unicode differences
            and (dbfile_ext != 'h5' or db.DTYPES['astromon_regions'][n].char != 'S')
        ]
        assert not dtypes_differ, f' dtypes differ: {dtypes_differ}'
        for name in regions.colnames:
            if db.DTYPES['astromon_regions'][name].char != 'S':
                assert np.all(regions[name] == regions_ref[name])

        # removing so it is empty
        db.remove_regions([2], dbfile=dbfile)
        regions = db.get_table('astromon_regions', dbfile=dbfile)
        assert regions.colnames == regions_ref.colnames, 'col names'
        # assert regions.dtype == regions_ref.dtype, 'dtypes'
        dtypes_differ = [
            (n, regions.dtype[n], db.DTYPES['astromon_regions'][n])
            for n in regions.dtype.names if regions.dtype[n] != db.DTYPES['astromon_regions'][n]
            # skip comparison on HDF5 because of string/unicode differences
            and (dbfile_ext != 'h5' or db.DTYPES['astromon_regions'][n].char != 'S')
        ]
        assert not dtypes_differ, f' dtypes differ: {dtypes_differ}'
        assert len(regions) == 0

        # remove non-existent
        db.remove_regions([3], dbfile=dbfile)  # silently removes nothing

        # adding a few and with autoincrementing region_id
        regions_1 = Table([
            {'ra': 0., 'dec': 0., 'radius': 5, 'obsid': 0, 'user': 'me', 'comments': ''},
            {'ra': 1., 'dec': 0., 'radius': 5, 'obsid': 0, 'user': 'them', 'comments': ''}
        ])
        db.add_regions(regions_1, dbfile=dbfile)
        regions = db.get_table('astromon_regions', dbfile=dbfile)
        regions_ref = Table(
            [[3, 4], [0.0, 1.0], [0.0, 0.0], [5, 5], [0, 0], ['me', 'them'], ['', '']],
            dtype=db.ASTROMON_REGION_DTYPE,
            names=[name for name in db.ASTROMON_REGION_DTYPE.names]
        )
        assert regions.colnames == regions_ref.colnames, 'col names'
        for name in regions.colnames:
            if db.DTYPES['astromon_regions'][name].char != 'S':
                assert np.all(regions[name] == regions_ref[name])
