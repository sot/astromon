"""
"""

import logging
import requests
from pathlib import Path

import numpy as np

from astropy import table
from astropy import coordinates as coords, units as u
from astropy import io

from astroquery.vizier import Vizier

from Ska.DBI import DBI
from cxotime import CxoTime
import astromon
from astromon import db

logger = logging.getLogger('astromon')


def _get_vizier(source, ra, dec, time, radius):
    """
    This fetches the vizier url, but it doesn't parse the result.
    """
    time = CxoTime(time)
    url = (
        'http://vizier.cfa.harvard.edu/viz-bin/asu-tsv?'
        '-source={source}'
        '&-sort=_r'
        '&-c={ra:.5f}+{dec:.4f}'
        '&-out.add=_RA(J2000,J{year})'
        '&-c.rs={radius}'
    )
    opts = {'source': source, 'ra': ra, 'dec': dec, 'year': time.frac_year, 'radius': radius}
    rc = requests.get(url.format(**opts))
    return io.ascii.read(rc.content.decode())


CROSS_MATCH_DTYPE = np.dtype([
    # ('obsid', 'int'),
    # ('id', 'int'),
    ('catalog', '<U24'),
    ('name', '<U16'),
    ('ra', float),
    ('dec', float),
    ('mag', float)
])


def get_vizier(pos, catalog, cat_identifier, name_cols, columns):
    vizier_result = Vizier.query_region(pos, radius=3 * u.arcsec, catalog=cat_identifier)
    vizier_result = [r for r in vizier_result]
    if len(vizier_result) == 0:
        logger.debug(f'{catalog} at {pos} has no results')
        return table.Table(dtype=CROSS_MATCH_DTYPE)
    vizier_result = table.vstack(vizier_result, metadata_conflicts='silent')
    vizier_result['catalog'] = catalog
    vizier_result['name'] = ['-'.join([str(row[n]) for n in name_cols]) for row in vizier_result]

    result = table.Table(data=np.zeros(len(vizier_result), dtype=CROSS_MATCH_DTYPE))
    for col in CROSS_MATCH_DTYPE.names:
        src_col = columns.get(col, col)
        if src_col in vizier_result.colnames:
            result[col] = vizier_result[src_col]
        else:
            result[col] = table.MaskedColumn(
                dtype=CROSS_MATCH_DTYPE[col],
                length=len(vizier_result),
                mask=np.ones(len(vizier_result))
            )
    logger.debug(f'{catalog} at {pos} has {len(vizier_result)} results')
    return result


def rough_match(sources, time):
    pos = coords.SkyCoord(
        ra=sources['ra'], dec=sources['dec'],
        unit='deg',
        frame='icrs',
        obstime=time
    )

    res = [
        get_vizier(
            pos,
            catalog='Tycho2',
            cat_identifier='I/259/tyc2',
            name_cols=['TYC1', 'TYC2', 'TYC3'],
            columns={'ra': 'RA_ICRS_', 'dec': 'DE_ICRS_', 'mag': 'VTmag'}
        ),
        get_vizier(
            pos,
            catalog='ICRS',
            cat_identifier='I/251',
            name_cols=['ICRS'],
            columns={'ra': 'RAJ2000', 'dec': 'DEJ2000'}
        ),
        get_vizier(
            pos,
            catalog='USNO-B1.0',
            cat_identifier='USNO-B1.0',
            name_cols=['USNO-B1.0'],
            columns={'ra': 'RAJ2000', 'dec': 'DEJ2000', 'mag': 'R1mag'}
        ),
        get_vizier(
            pos,
            catalog='2MASS',
            cat_identifier='II/246/out',
            name_cols=['_2MASS'],
            columns={'ra': 'RAJ2000', 'dec': 'DEJ2000', 'mag': 'Kmag'}
        ),
        get_vizier(
            pos,
            catalog='SDSS',
            cat_identifier='II/294',
            name_cols=['SDSS'],
            columns={'ra': 'RAJ2000', 'dec': 'DEJ2000', 'mag': 'rmag'}
        ),
    ]

    return table.vstack([r for r in res], metadata_conflicts='silent')


def do_sql_cross_match(selection_name):
    sql_script = Path(astromon.__file__).parent / 'sql' / 'x-corr' / f'{selection_name}.sql'
    if not sql_script.exists():
        logging.error(f'File {sql_script} does not exist')
        sql_script = Path(selection_name)
    if not sql_script.exists():
        logging.error(f'File {sql_script} does not exist')
        msg = f'{selection_name} is not a known selection'
        logging.error(msg)
        raise Exception(msg)

    with open(sql_script) as fh:
        sql_query = fh.read()
        dbi = DBI('sqlite', db.FILE, numpy=False)
        return table.Table(dbi.fetchall(sql_query))
