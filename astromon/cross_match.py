"""
"""

import re
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
    opts = {
        'source': source,
        'ra': float(ra),
        'dec': float(dec),
        'year': time.frac_year,
        'radius': float(radius)
    }
    rc = requests.get(url.format(**opts))
    output = rc.content.decode()
    n_rows = len([line for line in output.split('\n') if not re.match('#', line) and line.strip()])
    if n_rows:
        return io.ascii.read(output, data_start=3)


CROSS_MATCH_DTYPE = np.dtype([
    # ('obsid', 'int'),
    # ('id', 'int'),
    ('catalog', '<U24'),
    ('name', '<U16'),
    ('ra', float),
    ('dec', float),
    ('mag', float)
])


def get_vizier(pos, catalog, cat_identifier, name_cols, columns,
               radius=3 * u.arcsec, logging_tag='', use_astroquery=True, raw=False):
    if use_astroquery:
        vizier = Vizier(
            columns=[
                f'_RA(J2000,{pos.obstime.frac_year:8.3f})',
                f'_DE(J2000,{pos.obstime.frac_year:8.3f})'
            ]
        )
        vizier_result = vizier.query_region(pos, radius=radius, catalog=cat_identifier)
        vizier_result = [r for r in vizier_result]
    else:
        vizier_result = [
            _get_vizier(cat_identifier, p.ra / u.deg, p.dec / u.deg, pos.obstime, radius / u.arcsec)
            for p in pos
        ]
        vizier_result = [r for r in vizier_result if r is not None]

    if raw:
        return vizier_result

    if len(vizier_result) == 0:
        logger.debug(f'{logging_tag} {catalog:>24s} has no results')
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
    logger.debug(f'{logging_tag} {catalog:>24s} has {len(vizier_result)} results')
    return result


VIZIER_CATALOGS = {
    'Tycho2': dict(
        catalog='Tycho2',
        cat_identifier='I/259/tyc2',
        name_cols=['TYC1', 'TYC2', 'TYC3'],
        columns={
            'ra': '_RAJ2000_{time.frac_year:.3f}',
            'dec': '_DEJ2000_{time.frac_year:.3f}',
            'mag': 'VTmag'
        },
    ),
    # https://vizier.u-strasbg.fr/viz-bin/VizieR?-source=I/323
    'ICRS': dict(
        catalog='ICRS',
        cat_identifier='I/323',
        name_cols=['ICRF'],
        columns={
            'ra': '_RAJ2000',
            'dec': '_DEJ2000'
        },
    ),
    # https://vizier.u-strasbg.fr/viz-bin/VizieR?-source=I/284
    'USNO-B1.0': dict(
        catalog='USNO-B1.0',
        cat_identifier='USNO-B1.0',
        name_cols=['USNO-B1.0'],
        columns={
            'ra': '_RAJ2000',
            'dec': '_DEJ2000',
            'mag': 'R1mag'
        },
    ),
    # http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=I/239
    'HIP': dict(
        catalog='HIP',
        cat_identifier='I/239/hip_main',
        name_cols=['HIP'],
        columns={
            'ra': '_RAJ2000',
            'dec': '_DEJ2000',
            'mag': 'Vmag'
        },
    ),
    # http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=I/311
    # https://heasarc.gsfc.nasa.gov/W3Browse/all/hipnewcat.html
    'HIP2': dict(
        catalog='HIP',
        cat_identifier='I/311/hip2',
        name_cols=['HIP'],
        columns={
            'ra': '_RAJ2000',
            'dec': '_DEJ2000',
            'mag': 'Vmag',
        }
    ),
    # http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=I/322
    'UCAC4': dict(
        catalog='UCAC4',
        cat_identifier='I/322',
        name_cols=['UCAC4'],
        columns={
            'ra': '_RAJ2000',
            'dec': '_DEJ2000',
            'mag': 'f.mag'
        },
    ),
    '2MASS': dict(
        catalog='2MASS',
        cat_identifier='II/246/out',
        name_cols=['_2MASS'],
        columns={
            'ra': '_RAJ2000',
            'dec': '_DEJ2000',
            'mag': 'Kmag'},
    ),
    'SDSS': dict(
        catalog='SDSS',
        cat_identifier='II/294',
        name_cols=['SDSS'],
        columns={
            'ra': '_RAJ2000_{time.frac_year:.3f}',
            'dec': '_DEJ2000_{time.frac_year:.3f}',
            'mag': 'rmag'
        },
    ),
    # https://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=I/345/gaia2
    'Gaia2': dict(
        catalog='Gaia2',
        cat_identifier='I/345/gaia2',
        name_cols=['Source'],
        columns={
            'ra': '_RA_ICRS',
            'dec': '_DE_ICRS',
            'mag': 'Gmag'
        },
    ),
}


def _get(
    catalog, time, pos=None, ra=None, dec=None,
    radius=3 * u.arcsec, logging_tag='', raw=False
):
    assert (ra is not None and dec is not None) or pos is not None, 'pos or ra/dec required'
    if pos is None:
        pos = coords.SkyCoord(
            ra=ra, dec=dec,
            unit='deg',
            frame='icrs',
            obstime=time
        )
    else:
        assert ra is None and dec is None, 'ra/dec not required if giving pos'
    params = VIZIER_CATALOGS[catalog].copy()
    columns = {}
    for name, col in params['columns'].items():
        columns[name] = col.format(time=time)
    params['columns'] = columns
    return get_vizier(
        pos,
        radius=radius,
        logging_tag=logging_tag,
        **params,
        raw=raw
    )


def rough_match(
    sources, time, radius=3 * u.arcsec,
    catalogs=('Tycho2', 'ICRS', 'USNO-B1.0', '2MASS', 'SDSS')
):
    """
    Find sources in a set of standard catalogs around the x-ray sources given in `sources`,
    within an angular separation of at most `radius` (an astropy quantity).

    This function uses the 'ra', 'dec', 'id' and 'obsid' columns of the `sources` table.
    The obsid column is optional and is used only to check that the query corresponds to a single
    OBSID. Trying to query more than one OBSID at a time causes unexpected results after the
    astropy table join
    (see warning on Custom Join Functions: https://docs.astropy.org/en/stable/table/operations.html)

    WARNING: the current implementation of this function does create a temporary table with the
    cartesian product of the x-ray sources and the catalog sources tables
    (length ~ n_x_ray_sources^2). This is intended to be used with few x-ray sources at a time.
    """
    if len(sources) == 0:
        return []

    if 'obsid' in sources.dtype.names:
        assert len(np.unique(sources['obsid'])) <= 1, 'rough_match only handles one OBSID at a time'
        logging_tag = f'OBSID={sources["obsid"][0]} '
    else:
        logging_tag = ''

    logger.debug(f'{logging_tag} rough_match started')
    pos = coords.SkyCoord(
        ra=sources['ra'], dec=sources['dec'],
        unit='deg',
        frame='icrs',
        obstime=time
    )

    res = [
        _get(
            pos=pos,
            time=time,
            radius=radius,
            catalog=name,
            logging_tag=logging_tag
        ) for name in catalogs]
    res = table.vstack([r for r in res], metadata_conflicts='silent')

    if len(sources) and len(res):
        sources['coord_xray'] = coords.SkyCoord(sources['ra'], sources['dec'], unit='deg')
        res['coord_cat'] = coords.SkyCoord(res['ra'], res['dec'], unit='deg')
        sources['x_id'] = sources['id']
        res = table.join(res, sources[['coord_xray', 'x_id']], join_type='cartesian')
        sep = res['coord_xray'].separation(res['coord_cat'])
        res['separation'] = sep.to_value(unit='arcsec')
        res.remove_columns(['coord_xray', 'coord_cat'])
        res = res[sep < radius]
    return res


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


def cross_match(
    name,
    astromon_obs=None,
    astromon_xray_src=None,
    astromon_cat_src=None,
    dbfile=None,
    logging_tag='',
):
    """
    Cross-match x-ray sources with catalog counterparts.
    """
    if astromon_xray_src is None:
        astromon_xray_src = db.get_table('astromon_xray_src', dbfile)
    if astromon_cat_src is None:
        astromon_cat_src = db.get_table('astromon_cat_src', dbfile)
    if astromon_obs is None:
        astromon_obs = db.get_table('astromon_obs', dbfile)

    if name == 'standard_xcorr':
        return _standard_cross_match(
            astromon_obs,
            astromon_xray_src,
            astromon_cat_src,
            name=name,
            catalogs=['Tycho2', 'SIMBAD_high', 'CELMON', 'ICRS', 'ASTROMON'],
            snr=5.0,
            r_angle=120,
            start=CxoTime() - 5 * 365 * u.day,
        )
    elif name == 'oaa4_snr4':
        return _standard_cross_match(
            astromon_obs,
            astromon_xray_src,
            astromon_cat_src,
            name=name,
            catalogs=['Tycho2', 'SIMBAD_high', 'CELMON', 'ICRS', 'ASTROMON',
                      'SDSS', '2MASS', 'USNO-B1.0'],
            snr=4.0,
            r_angle=240,
        )
    elif name == 'astromon_21':
        return _simple_cross_match(
            astromon_obs,
            astromon_xray_src,
            astromon_cat_src,
            name=name,
            catalogs=['ICRS', 'Tycho2'],
            snr=3,
            r_angle=120.,
            logging_tag=logging_tag,
        )
    else:
        raise Exception(f'Unknown x-matching name: "{name}"')


def _standard_cross_match(
    astromon_obs,
    astromon_xray_src,
    astromon_cat_src,
    name,
    catalogs,
    snr,
    r_angle,
    start=None
):
    """
    Standard cross X-ray -- Catalog cross correlation query

    - X-ray sources within 3 arcmin off-axis (NONE) or 0.4 arcmin (grating)
    - X-ray snr >  5.0
    - X-ray extr_rad_grating,r,a,0.4,,,"Extr. rad. around grating source (arcmin)"
    - X-ray - Catalog position match within 3 arcsec
    - Catalog is high precision 'Tycho2', 'SIMBAD_high', 'CELMON', 'ICRS', 'ASTROMON'
    """

    # astromon_obs = astromon_obs[(astromon_obs['process_status'] == 'OK')]

    if start is not None:
        date_obs = CxoTime(astromon_obs['date_obs'].astype(str))
        astromon_obs = astromon_obs[date_obs > start]

    # only X-ray sources with SNR > 5 and no near-neighbors
    astromon_xray_src = astromon_xray_src[
        (astromon_xray_src['snr'] > snr)
        & (astromon_xray_src['near_neighbor_dist'] > 6.0)
        # astromon_xray_src['status_id'] = 0)
    ]

    # only some catalogs
    astromon_cat_src = astromon_cat_src[np.in1d(astromon_cat_src['catalog'], catalogs)]

    result = table.join(
        astromon_obs,
        astromon_xray_src,
        keys=['obsid'],
    )
    result = table.join(
        result,
        astromon_cat_src,
        keys=['obsid', 'x_id'],
    )

    astromon_xray_src.rename_columns(
        ['name_1', 'id_1', 'ra_1', 'dec_1', 'y_angle_1', 'z_angle_1'],
        ['x_name', 'x_id', 'x_ra', 'x_dec', 'x_y_angle', 'x_z_angle']
    )
    astromon_cat_src.rename_columns(
        ['name_2', 'id_2', 'ra_2', 'dec_2', 'y_angle_2', 'z_angle_2'],
        ['c_name', 'c_id', 'c_ra', 'c_dec', 'c_y_angle', 'c_z_angle']
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


def _simple_cross_match(
    astromon_obs,
    astromon_xray_src,
    astromon_cat_src,
    name,
    catalogs,
    snr,
    r_angle=120.,
    start=None,
    logging_tag='',
):
    astromon_xray_src = astromon_xray_src[astromon_xray_src['snr'] > snr]
    astromon_cat_src = astromon_cat_src[np.in1d(astromon_cat_src['catalog'], catalogs)]

    if start is not None:
        date_obs = CxoTime(astromon_obs['date_obs'].astype(str))
        astromon_obs = astromon_obs[date_obs > start]

    if len(astromon_xray_src) == 0 or len(astromon_cat_src) == 0:
        logger.debug(f'{logging_tag} No xray or cat sources')
        return []
    matches = table.join(
        astromon_obs,
        astromon_xray_src,
        keys=['obsid'],
    )
    matches = table.join(
        matches,
        astromon_cat_src,
        keys=['obsid', 'x_id'],
        table_names=['xray', 'cat']
    )
    assert np.all(matches['id_xray'] == matches['x_id'])
    matches.rename_column('id_cat', 'c_id')

    matches['dz'] = matches['z_angle_xray'] - matches['z_angle_cat']
    matches['dy'] = matches['y_angle_xray'] - matches['y_angle_cat']
    matches['dr'] = np.sqrt(matches['dy']**2 + matches['dz']**2)
    matches['cat_order'] = 200
    for i, k in enumerate(catalogs):
        matches[matches['catalog'] == k]['cat_order'] = i

    matches = matches[matches['r_angle'] < r_angle]

    if len(matches) == 0:
        return []

    matches = matches.group_by('obsid')
    for g in matches.groups:
        g.sort(['cat_order', 'dr'])
    result = table.vstack([g[0] for g in matches.groups])

    result['select_name'] = name

    return result[['select_name', 'obsid', 'c_id', 'x_id', 'dy', 'dz', 'dr']]
