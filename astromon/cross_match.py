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
from astromon import db, utils

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
        logger.debug(f'{catalog:>24s} has no results')
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
    logger.debug(f'{catalog:>24s} has {len(vizier_result)} results')
    return result


def rough_match(sources, time):
    logger.debug('rough_match started')
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


@utils.logging_call_decorator
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


@utils.logging_call_decorator
def cross_match(name):
    if name == 'standard_xcorr':
        return _standard_cross_match(
            name=name,
            catalogs=['Tycho2', 'SIMBAD_high', 'CELMON', 'ICRS', 'ASTROMON'],
            snr=5.0,
            r_angle=120,
            start=CxoTime() - 5 * 365 * u.day,
        )
    elif name == 'oaa4_snr4':
        return _standard_cross_match(
            name=name,
            catalogs=['Tycho2', 'SIMBAD_high', 'CELMON', 'ICRS', 'ASTROMON',
                      'SDSS', '2MASS', 'USNO-B1.0'],
            snr=4.0,
            r_angle=240,
        )
    else:
        raise Exception(f'Unknown x-matching name: "{name}"')


def _standard_cross_match(name, catalogs, snr, r_angle, start=None, update=False):
    """
    Standard cross X-ray -- Catalog cross correlation query

    - X-ray sources within 3 arcmin off-axis (NONE) or 0.4 arcmin (grating)
    - X-ray snr >  5.0
    - X-ray extr_rad_grating,r,a,0.4,,,"Extr. rad. around grating source (arcmin)"
    - X-ray - Catalog position match within 3 arcsec
    - Catalog is high precision 'Tycho2', 'SIMBAD_high', 'CELMON', 'ICRS', 'ASTROMON'
    """
    astromon_xcorr = db.get('astromon_xcorr')
    astromon_xray_src = db.get('astromon_xray_src')
    astromon_cat_src = db.get('astromon_cat_src')
    astromon_obs = db.get('astromon_obs')

    astromon_obs = astromon_obs[(astromon_obs['process_status'] == 'OK')]
    if update:
        astromon_obs = astromon_obs[~np.in1d(astromon_obs['obsid'], astromon_xcorr['obsid'])]

    date_obs = CxoTime(astromon_obs['date_obs'].astype(str))
    if start is not None:
        astromon_obs = astromon_obs[date_obs > start]

    # only X-ray sources with SNR > 5 and no near-neighbors
    astromon_xray_src = astromon_xray_src[
        (astromon_xray_src['snr'] > snr)
        & (astromon_xray_src['near_neighbor_dist'] > 6.0)
        # & (astromon_xray_src['status_id'] = None or astromon_xray_src['status_id'] = 0)
    ]

    # only some catalogs
    astromon_cat_src = astromon_cat_src[np.in1d(astromon_cat_src['catalog'], catalogs)]

    # renaming columns before join to not have to deal with name collisions
    astromon_cat_src.rename_columns(
        ['name', 'id', 'ra', 'dec', 'y_angle', 'z_angle'],
        ['c_name', 'c_id', 'c_ra', 'c_dec', 'c_y_angle', 'c_z_angle']
    )
    astromon_xray_src.rename_columns(
        ['name', 'id', 'ra', 'dec', 'y_angle', 'z_angle'],
        ['x_name', 'x_id', 'x_ra', 'x_dec', 'x_y_angle', 'x_z_angle']
    )

    result = table.join(
        astromon_obs,
        astromon_xray_src,
        keys=['obsid'],
    )
    result = table.join(
        result,
        astromon_cat_src,
        keys=['obsid'],
    )

    result['dy'] = result['c_y_angle'] - result['x_y_angle']
    result['dz'] = result['c_z_angle'] - result['x_z_angle']
    result['dr2'] = result['dy']**2 + result['dz']**2
    result['dr'] = np.sqrt(result['dr2'])
    result['select_name'] = name

    result = result[
        ((result['r_angle'] < 24) | ((result['grating'] == 'NONE') & (result['r_angle'] < r_angle)))
        & (result['dr2'] < 9.0)
    ]

    return result[['select_name', 'obsid', 'c_id', 'x_id', 'dy', 'dz', 'dr']]
