#!/usr/bin/env python
"""
Script to find x-ray sources in observations, and a tentative set of optical/radio counterparts.
"""

import os
import re
import logging
import argparse
import tempfile
from pathlib import Path
import sqlite3

from multiprocessing import Pool

from cxotime import CxoTime, units as u

from astropy.table import Table

import stk

from Ska.arc5gl import Arc5gl
from astromon.observation import Observation
from astromon.cross_match import rough_match
from astromon import db, utils

logger = logging.getLogger('astromon')


def get_obsids(tstart, tstop):
    tstart = CxoTime(tstart)
    tstop = CxoTime(tstop)
    with tempfile.TemporaryDirectory() as td, utils.chdir(td):
        path = Path(td)
        arc5gl = Arc5gl()
        arc5gl.sendline(f'tstart={tstart.date}')
        arc5gl.sendline(f'tstop={tstop.date}')
        arc5gl.sendline('get obspar')
        names = [str(p) for p in path.glob('*')]
        return [re.search('axaff([0-9]+)_', name).group(1) for name in names]


def process(obsid, workdir, db_file):
    """
    This is where the actual work is done.
    """
    try:
        logger.info(f'Processing OBSID {obsid}')
        logger.info('-----------------------')
        observation = Observation(obsid, workdir)
        observation.process()
        obspar = Table([observation.get_obspar()])
        sources = observation.get_sources()
        matches = rough_match(sources, CxoTime(observation.get_obsid_info()['date']))

        logger.debug(f'About to update {db_file}')
        with sqlite3.connect(db_file) as con:
            db.save(con, 'astromon_obs', obspar)
            if len(sources):
                db.save(con, 'astromon_xray_src', sources)
            if len(matches):
                db.save(con, 'astromon_cat_src', matches)
    except Exception as e:
        logger.error(f'OBSID {obsid} failed: {e}')


def get_parser():
    """
    Get the argument parser.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--db-file',
        help='SQLite file where data is saved',
        default='ASTROMON_table.rdb',
        type=Path
    )
    parser.add_argument(
        '--workdir',
        help='Working directory. Default is to create a tmp directory',
        default=None
    )
    parser.add_argument(
        '--start',
        help='Start of the time range to process',
    )
    parser.add_argument(
        '--stop',
        help='End of the time range to process',
    )
    parser.add_argument(
        '--obsid',
        help='An OBSID to process or a file with a list of OBSIDs',
    )
    parser.add_argument(
        '--log-level',
        help='Logging severity level',
        choices=['debug', 'info', 'error', 'fatal'],
        default='debug'
    )
    return parser


def main():
    """
    Main routine that deals with processing arguments, output and such.
    """

    parser = get_parser()
    args = parser.parse_args()

    # not using pyyaks because CIAO overrides the default logger class using logging.setLoggerClass
    logger = logging.getLogger('astromon')
    hdlr = logging.StreamHandler()
    hdlr.setLevel(args.log_level.upper())
    logger.addHandler(hdlr)
    logger.setLevel(args.log_level.upper())
    logging.getLogger('astromon.utils').setLevel('WARNING')

    if not args.db_file.exists():
        logger.info(f'File does not exist: {args.db_file}. Will create')
    args.db_file.parent.mkdir(exist_ok=True, parents=True)

    if args.workdir is None:
        tmp = tempfile.TemporaryDirectory()
        args.workdir = tmp.name
    Path(args.workdir).mkdir(exist_ok=True, parents=True)

    args.stop = CxoTime(args.stop) if args.stop is not None else CxoTime()
    args.start = CxoTime(args.start) if args.start is not None else args.stop - 30 * u.day
    obsids = get_obsids(args.start, args.stop)

    if args.obsid is not None:
        if os.path.exists(args.obsid) and os.path.isfile(args.obsid):
            # '@-' builds stack but does not include path name
            obsids_2 = stk.build("@-" + args.obsid)
        else:
            obsids_2 = stk.build(args.obsid)
        obsids = [obsid for obsid in obsids_2 if obsid in obsids]

    logger.info(f'will process the following obsids: {", ".join(obsids)}')
    task_args = [(obsid, args.workdir, args.db_file) for obsid in obsids]

    # '9' because the archive limits number of concurrent connections
    # from same IP address (well, it did when it was 'ftp').
    with Pool(processes=9) as pool:
        pool.starmap(process, task_args)


if __name__ == '__main__':
    main()
