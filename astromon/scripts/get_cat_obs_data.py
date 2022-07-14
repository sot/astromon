#!/usr/bin/env python
"""
Script to find x-ray sources in observations, and a tentative set of optical/radio counterparts.
"""

import sys
import os
import re
import argparse
import tempfile
import logging
import traceback
from pathlib import Path
import shutil

from multiprocessing import Pool, set_start_method
import numpy as np

from cxotime import CxoTime, units as u

from astropy.table import Table, Column, vstack

import stk

from chandra_aca.transform import radec_to_yagzag
from Quaternion import Quat
from Ska.arc5gl import Arc5gl
from astromon.observation import Observation, Skipped, SkippedWithWarning
from astromon.cross_match import rough_match, compute_cross_matches
from astromon import db, utils
import pyyaks.logger


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


def save(data, db_file):
    logger = logging.getLogger(name='astromon')

    errors = {}
    skipped = {}
    n_skip = 0
    for d in [d for d in data if d['msg']]:
        if d['ok']:
            obs = skipped.get(d['msg'], [])
            skipped[d['msg']] = obs + [str(d['obsid'])]
        else:
            obs = errors.get(d['msg'], [])
            errors[d['msg']] = obs + [str(d['obsid'])]
        n_skip += 1
    n = len(data)
    if n_skip:
        logger.warning(f'WARNING: {n} observations , {n_skip} skipped.')
    else:
        logger.info(f'{n} observations , {n_skip} skipped.')
    for k in skipped:
        logger.warning(f'WARNING: {len(skipped[k])} {k} ({", ".join(skipped[k])})')
    for k in errors:
        logger.warning(f'ERROR: {len(errors[k])} {k} ({", ".join(errors[k])})')

    data = [d for d in data if 'astromon_obs' in d and len(d['astromon_obs']) > 0]
    if len(data) == 0:
        logger.info('Nothing to save')
        return

    names = ['astromon_obs', 'astromon_xray_src', 'astromon_cat_src', 'astromon_xcorr']
    # these following baroque lines are here because there are some columns we cannot vstack
    # so I decided to only vstack the columns that are in the dtype,
    # but it turns out that not all columns in the dtype are actually in the data...
    # and the data itself could be an empty array...
    data = {
        name: [d[name] for d in data if len(d[name])] for name in names
    }
    for name in names:
        if len(data[name]):
            cols = [col for col in db.DTYPES[name].names if col in data[name][0].colnames]
            data[name] = [d[cols] for d in data[name]]
            data[name] = vstack(data[name], metadata_conflicts='silent')

    logger.debug(f'About to write {len(data["astromon_obs"])} observations to {db_file}')
    with db.connect(db_file, mode='r+') as con:
        for name in names:
            if len(data[name]):
                db.save(name, data[name], con)


def process(obsid, workdir, log_level, archive_dir):
    """
    This is where the actual work is done.
    """
    logger = pyyaks.logger.get_logger(name='astromon', level=log_level)
    try:
        logger.info(f'OBSID={obsid} *** Processing OBSID {obsid} ***')
        observation = Observation(obsid, workdir, archive_dir=archive_dir)
        observation.process()
        if archive_dir:
            observation.archive()

        obspar = Table([observation.get_info()])
        sources = observation.get_sources()
        match_candidates = rough_match(sources, CxoTime(observation.get_obspar()['date_obs']))

        if len(match_candidates):
            q = Quat(equatorial=(obspar['ra_pnt'][0], obspar['dec_pnt'][0], obspar['roll_pnt'][0]))
            match_candidates['obsid'] = obsid
            match_candidates['id'] = np.arange(len(match_candidates))
            match_candidates['y_angle'], match_candidates['z_angle'] = radec_to_yagzag(
                match_candidates['ra'], match_candidates['dec'], q
            )
        else:
            match_candidates['obsid'] = Column(dtype=int)
            match_candidates['id'] = Column(dtype=int)
            match_candidates['y_angle'] = Column(dtype=np.float32)
            match_candidates['z_angle'] = Column(dtype=np.float32)

        logger.debug(f'OBSID={obsid} About to cross-match')
        matches = compute_cross_matches(
            'astromon_21',
            astromon_obs=obspar,
            astromon_xray_src=sources,
            astromon_cat_src=match_candidates,
            logging_tag=f'OBSID={obsid}',
        )
        if matches:
            matches = matches[['select_name', 'obsid', 'c_id', 'x_id', 'dy', 'dz', 'dr']]

        result = {
            'ok': True,
            'msg': '',
            'obsid': obsid,
            'astromon_obs': obspar,
            'astromon_xray_src': sources,
            'astromon_cat_src': match_candidates,
            'astromon_xcorr': matches,
        }

        (observation.workdir / 'results').mkdir(exist_ok=True)
        for name in ['astromon_obs', 'astromon_xray_src', 'astromon_cat_src', 'astromon_xcorr']:
            fn = observation.workdir / 'results' / f'{name}.fits'
            fn.unlink(missing_ok=True)
            if len(result[name]) > 0:
                result[name].write(fn)

        if archive_dir:
            observation.archive(
                'astromon_obs.fits',
                'astromon_xray_src.fits',
                'astromon_cat_src.fits',
                'astromon_xcorr.fits',
            )
        return result

    except Skipped as e:
        ok = True
        msg = f'skipped: {e}'
        logger.info(f'OBSID={obsid} skipped')
    except SkippedWithWarning as e:
        ok = True
        msg = f'skipped: {e}'
        logger.warning(f'OBSID={obsid} WARNING - skipped: {e}')
    except Exception as e:
        ok = False
        msg = f'error: {e}'
        logger.error(f'OBSID={obsid} failed: {e}')
        exc_type, exc_value, exc_traceback = sys.exc_info()
        trace = traceback.extract_tb(exc_traceback)
        for step in trace:
            logger.error(f'OBSID={obsid}            in {step.filename}:{step.lineno}/{step.name}:')
            logger.error(f'OBSID={obsid}            {step.line}')

    return {
        'ok': ok,
        'msg': msg,
        'obsid': obsid,
        'astromon_obs': [],
        'astromon_xray_src': [],
        'astromon_cat_src': [],
        'astromon_xcorr': [],
    }


def get_parser():
    """
    Get the argument parser.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--db-file',
        help='SQLite file where data is saved',
        default='astromon.h5',
        type=Path
    )
    parser.add_argument(
        '--workdir',
        help='Working directory. A temp directory is created for each observation.',
        default=None
    )
    parser.add_argument(
        '--archive-dir',
        help='Archive directory.',
        default=None
    )
    parser.add_argument(
        '--start',
        help='Start of the time range to process. Default: stop - 30 days.',
    )
    parser.add_argument(
        '--stop',
        help='End of the time range to process. Default: now.',
    )
    parser.add_argument(
        '--obsid',
        help='An OBSID to process or a file with a list of OBSIDs',
    )
    parser.add_argument(
        '--no-exception',
        help='Do not skip observations already in file',
        action='store_true',
        default=False,
    )
    parser.add_argument(
        '--threads',
        help='The number of processes to use.',
        # '9' because the archive limits number of concurrent connections
        # from same IP address (well, it did when it was 'ftp').
        default=9,
        type=int,
    )
    parser.add_argument(
        '--log-level',
        help='Logging severity level',
        choices=['debug', 'info', 'warning', 'error', 'fatal'],
        default='debug'
    )
    return parser


def main():
    """
    Main routine that deals with processing arguments, output and such.
    """
    set_start_method('spawn')

    parser = get_parser()
    args = parser.parse_args()

    logger = pyyaks.logger.get_logger(name='astromon', level=args.log_level.upper())

    if args.workdir:
        workdir = Path(args.workdir)
        workdir.mkdir(exist_ok=True, parents=True)
    else:
        # if args.workdir is not given, create one for the result file
        # do not pass this downstream, because it can fill the tmp directory
        tmpdir = tempfile.TemporaryDirectory()  # this variable should live until the end
        workdir = Path(str(tmpdir))

    # all changes go in this file which is copied back to its final destination at the end.
    db_file = workdir / args.db_file.name
    db_file.parent.mkdir(exist_ok=True, parents=True)
    if args.db_file.exists():
        shutil.copy(args.db_file, db_file)
    else:
        logger.info(f'File does not exist: {args.db_file}. Will create')
        for name in db.DTYPES:
            tab = db.create_table(name)
            db.save(name, tab, dbfile=db_file)

    if args.obsid is not None:
        if os.path.exists(args.obsid) and os.path.isfile(args.obsid):
            # '@-' builds stack but does not include path name
            obsids = stk.build("@-" + args.obsid)
        else:
            obsids = stk.build(args.obsid)
    else:
        args.stop = CxoTime(args.stop) if args.stop is not None else CxoTime()
        args.start = CxoTime(args.start) if args.start is not None else args.stop - 30 * u.day
        obsids = get_obsids(args.start, args.stop)

    if len(obsids) == 0:
        logger.info('No OBSIDs found as specified')
        return

    if db_file.exists() and not args.no_exception:
        # "no_exception" means all OBSIDs are processes,
        # not no_exception means we need to look for existing OBSIDs and skip them
        obsids = np.array(obsids, dtype=int)
        try:
            astromon_obs = db.get_table('astromon_obs', dbfile=db_file)
            exceptions = np.in1d(obsids, astromon_obs['obsid'])
            n_exceptions = np.sum(exceptions)
            if n_exceptions:
                exceptions_str = str(obsids[exceptions])[1:-1]
                logger.info(
                    f'skipping {n_exceptions} OBSID{"s" if n_exceptions > 1 else ""} '
                    f'that {"are" if n_exceptions > 1 else "is"} on file already: {exceptions_str}'
                )
            obsids = obsids[~exceptions].astype(str)
        except utils.MissingTableException:
            pass

    if len(obsids) == 0:
        logger.info('All OBSIDs already processed')
        return

    logger.info(f'will process the following obsids: {", ".join(obsids)}')
    task_args = [
        (int(obsid), args.workdir, args.log_level.upper(), args.archive_dir)
        for obsid in obsids
    ]

    with Pool(processes=args.threads) as pool:
        results = pool.starmap(process, task_args)

    save(results, db_file)

    args.db_file.parent.mkdir(exist_ok=True, parents=True)
    shutil.copy(db_file, args.db_file)


if __name__ == '__main__':
    main()
