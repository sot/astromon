#!/usr/bin/env python

"""
Compute the Celestial location radius RMS corresponding to the PRD requirement
of 1.0 arcsec.
"""

import numpy as np
import json
import argparse
from pathlib import Path

from astropy import units as u

from cxotime import CxoTime
from astromon import db, cross_match


def get_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--offsets-file', type=Path)
    return parser


def main():
    args = get_parser().parse_args()

    dat = db.get_cross_matches()
    dat = dat[~cross_match.get_bad_target_mask(dat)]

    select_name = np.unique(dat['select_name'])
    assert len(select_name) <= 1, f'More than one selection'
    assert len(dat), 'no data'
    select_name = select_name[0]

    if args.offsets_file:
        print(f'Subtracting offsets from {args.offsets_file}')
        with open(args.offsets_file) as fh:
            offsets = json.load(fh)
        times = CxoTime([of['tstart'] for of in offsets])
        dy = np.array([0.] + list([of['dy'] for of in offsets]))
        dz = np.array([0.] + list([of['dz'] for of in offsets]))
        time_bin = np.digitize(dat['time'].cxcsec, times.cxcsec)
        dat['dy'] -= dy[time_bin]
        dat['dz'] -= dz[time_bin]
        dat['dr'] = np.sqrt(dat['dy']**2 + dat['dz']**2)

    # For all data
    print(f"Using all '{select_name}' data (last 5 years)")
    print("N srcs: {}".format(len(dat)))
    print('RMS radius', np.sqrt(np.mean(dat['dr'] ** 2)))
    # 0.410

    percentile = np.percentile(dat['dr'], 90)
    print(f"90 percentile radius = {percentile:.2f} arcsec")

    for detector in ['ACIS-S', 'ACIS-I', 'HRC-S', 'HRC-I']:
        det = dat['detector'] == detector
        percentile = np.percentile(dat['dr'][det], 90) if np.count_nonzero(det) else np.nan
        print(f"90 percentile radius for {detector} is {percentile:.2f} arcsec")

    print("{:.1f} percent outside a 1 arcsec radius".format(
        100.0 * np.count_nonzero(dat['dr'] > 1.0) / len(dat['dr'])))

    now = CxoTime()

    for start, stop in [(now - 4 * u.year, now - 2 * u.year), (now - 2 * u.year, now)]:
        print("---")
        ok = ((CxoTime(dat['date_obs']) > start)
              & (CxoTime(dat['date_obs']) < stop))
        print("{} to {}".format(start.date, stop.date))

        print("N srcs: {}".format(len(dat[ok])))
        print('RMS radius', np.sqrt(np.mean(dat[ok]['dr'] ** 2)))
        percentile = np.percentile(dat[ok]['dr'], 90)
        print(f"90 percentile radius = {percentile:.2f} arcsec")

        for detector in ['ACIS-S', 'ACIS-I', 'HRC-S', 'HRC-I']:
            det = dat[ok]['detector'] == detector
            percentile = np.percentile(dat[ok]['dr'][det], 90)
            print(f"90 percentile radius for {detector} is {percentile:.2f} arcsec")

        print("{:.1f} percent outside a 1 arcsec radius".format(
            100.0 * np.count_nonzero(dat[ok]['dr'] > 1.0) / len(dat[ok]['dr'])))


if __name__ == '__main__':
    main()
