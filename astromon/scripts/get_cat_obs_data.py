#!/usr/bin/env python
"""
Script to find x-ray sources in observations, and a tentative set of optical/radio counterparts.
"""

import os
# import sys
import logging
import argparse
import tempfile
from pathlib import Path
import sqlite3

from cxotime import CxoTime

from astropy.table import Table

import stk

from astromon.observation import Observation
from astromon.cross_match import rough_match
from astromon import db

logger = logging.getLogger('astromon')


def process(obsid, workdir, db_file):
    """
    This is where the actual work is done.
    """
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
        '--obsid',
        help='An OBSID to process or a file with a list of OBSIDs',
        required=True
    )
    parser.add_argument(
        '--log-level',
        help='Logging severity level',
        choices=['debug', 'info', 'error', 'fatal']
    )
    return parser


def main():
    """
    Main routine that deals wit processing arguments, output and such.
    """

    from ciao_contrib._tools.taskrunner import TaskRunner

    args = get_parser().parse_args()

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

    if os.path.exists(args.obsid) and os.path.isfile(args.obsid):
        # '@-' builds stack but does not include path name
        obsids = stk.build("@-" + args.obsid)
    else:
        obsids = stk.build(args.obsid)

    taskrunner = TaskRunner()
    for obsid in obsids:
        taskrunner.add_task(f"OBS_ID={obsid}", "", process, obsid, args.workdir, args.db_file)

    # '9' because the archive limits number of concurrent connections
    # from same IP address (well, it did when it was 'ftp').
    taskrunner.run_tasks(processes=9)


if __name__ == '__main__':
    main()
