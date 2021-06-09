import os
import sqlite3
import numpy as np
from pathlib import Path
from Ska.DBI import DBI

from astropy import table

if 'ASTROMON_FILE' in os.environ:
    FILE = Path(os.environ['ASTROMON_FILE'])
else:
    FILE = Path(os.environ['SKA']) / 'data' / 'astromon' / 'astromon.db3'


ASTROMON_XCORR_DTYPE = np.dtype([
    ('select_name', '<U14'),
    ('obsid', np.int8),
    ('c_id', np.int8),
    ('x_id', np.int8),
    ('dy', np.float),
    ('dz', np.float),
    ('dr', np.float)
])


ASTROMON_XRAY_SRC_DTYPE = np.dtype([
    ('obsid', np.int8),
    ('id', np.int8),
    ('name', '<U20'),
    ('ra', np.float64),
    ('dec', np.float64),
    ('net_counts', np.float),
    ('y_angle', np.float),
    ('z_angle', np.float),
    ('r_angle', np.float),
    ('snr', np.float),
    ('near_neighbor_dist', np.float),
    ('double_id', np.int8),
    ('status_id', np.int8)
])


ASTROMONCAT_SRC_DTYPE = np.dtype([
    ('obsid', np.int8),
    ('id', np.int8),
    ('catalog', '<U11'),
    ('name', '<U24'),
    ('ra', np.float64),
    ('dec', np.float64),
    ('mag', np.float),
    ('y_angle', np.float),
    ('z_angle', np.float)
])


ASTROMON_OBS_DTYPE = np.dtype([
    ('obsid', np.int8),
    ('version', np.float),
    ('detector', '<U6'),
    ('target', '<U28'),
    ('grating', '<U4'),
    ('sim_z', np.float),
    ('date_obs', '<U20'),
    ('tstart', ),
    ('fids', '<U20'),
    ('ascdsver', '<U32'),
    ('ra', np.float64),
    ('dec', np.float64),
    ('roll', np.float64),
    ('process_status', '<U229'),
    ('category_id', np.int)
])


def get(table_name, dbfile=None):
    """
    Get an ENTIRE table.
    """
    if dbfile is None:
        dbfile = FILE
    dbi = DBI('sqlite', dbfile)
    return table.Table(dbi.fetchall(f'select * from {table_name}'))


def _save(con, table_name, data):
    cur = con.cursor()

    if type(table_name) is list:
        all_table_names = table_name
        all_data = data
    else:
        all_table_names = [table_name]
        all_data = [data]

    for table_name, data in zip(all_table_names, all_data):
        tables = [
            row[0] for row in
            cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
        ]
        if table_name not in tables:
            with open(Path(__file__).parent / 'sql' / 'tables' / f'{table_name}.sql') as fh:
                cur.execute(fh.read())
        columns = cur.execute(f"PRAGMA table_info('{table_name}')").fetchall()
        table_column_names = [row[1] for row in columns]
        data_column_names = [name for name in data.colnames if name in table_column_names]
        if not data_column_names:
            raise Exception('Input data has no columns in common with table in DB')

        values = [str(tuple(row)) for row in data[data_column_names]]
        insert_query = f"insert into {table_name} ({', '.join(data_column_names)}) values"
        insert_query += ", ".join(values)
        insert_query += ";"
        insert_query = insert_query.replace('inf', '1e30')  # HORRIBLE HACK
        insert_query = insert_query.replace('masked', 'NULL')  # HORRIBLE HACK

        if 'obsid' in data.colnames:
            obsids = ', '.join(np.unique(data['obsid']).astype(str))
            cur.execute(f"DELETE FROM {table_name} WHERE OBSID IN ({obsids})")

        cur.execute(insert_query)


def save(db, table_name, data):
    """
    Insert data into a table, deleting previous entries for the same OBSID.

    If the table does not exist, it is created using pre-existing table definitions.
    """
    if type(db) is sqlite3.Connection:
        return _save(db, table_name, data)

    retries = 3
    for retry in range(retries):
        try:
            with sqlite3.connect(db) as con:
                _save(con, table_name, data)
        except sqlite3.OperationalError as e:
            # only OperationalError with msg='database is locked' are ignored up to {retries} times
            if str(e) != 'database is locked':
                raise
            if retry == retries - 1:
                raise Exception(f'Failed to save to DB after {retries} retries')
