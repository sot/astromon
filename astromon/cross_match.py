"""
"""

import logging
from numpy.core.shape_base import vstack
import requests
from pathlib import Path

from astropy import table
from astropy import coordinates as coords, units as u
from astropy import io

from astroquery.vizier import Vizier

from Ska.DBI import DBI
from cxotime import CxoTime
import astromon
from astromon import db

logger = logging.getLogger('astromon')


def get_vizier(source, ra, dec, time, radius):
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


def _get_tycho(sources, time, radius):
    matches = []
    for source in sources:
        match = get_vizier('I/259/TYC2', source['RA'], source['DEC'], time, radius)
        matches.append(match)
    return vstack(matches)


def rough_match(sources, time):
    pos = coords.SkyCoord(
        ra=sources['ra'], dec=sources['dec'],
        unit='deg',
        frame='icrs',
        obstime=time
    )
    # res = Vizier.query_region(pos, radius=3 * u.arcsec, catalog='Tycho2')
    res = Vizier.query_region(pos, radius=3 * u.arcsec, catalog='I/259/tyc2')
    res = [r for r in res]
    if not res:
        return []
    res = table.vstack([r for r in res], metadata_conflicts='silent')
    return res


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
