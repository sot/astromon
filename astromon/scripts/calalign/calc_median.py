#!/usr/bin/env python

"""
Calculate median offset in regular time bins.
"""

import json
import argparse
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.style
from astropy.table import Table

from cxotime import CxoTime

from astromon import db


matplotlib.style.use('bmh')

assert db.FILE.name[-3:] == '.h5', 'Using old astromon file, define ASTROMON_FILE env var.'


def get_bins(years, bins_per_year=2):
    years = np.atleast_1d(years)
    ymin = np.floor(np.min(years * bins_per_year)) / bins_per_year
    ymax = np.ceil(np.max(years * bins_per_year)) / bins_per_year
    nbins = int((ymax - ymin) * bins_per_year)
    bins = np.linspace(ymin, ymax, nbins + 1)
    return bins


def binned_median(years, vals, bins_per_year=2):
    """
    Compute the median value in each bin of the input data.
    """
    years = np.atleast_1d(years)
    bins = get_bins(years, bins_per_year=bins_per_year)
    years_bin = np.digitize(years, bins)
    t = Table([years_bin, years, vals], names=['years_bin', 'years', 'vals'])
    tg = t.group_by('years_bin')
    tga = tg.groups.aggregate(np.median)
    return tga['years'], tga['vals'], bins


def filter_xcorr(xcorr, *,
                 snr=None, r_angle=None,
                 start=None, stop=None,
                 **kwargs):
    """Filter the matched X-ray and catalog sources table.

    Any additional keyword arguments are used to filter the table like:
    - xcorr[key] == val  # scalar values
    - np.isin(xcorr[key], val)  # list values

    Returns the filtered Table.
    """

    ok = np.ones(len(xcorr), dtype=bool)

    if snr is not None:
        ok &= xcorr['snr'] >= snr

    if r_angle is not None:
        ok &= xcorr['r_angle'] <= r_angle

    if start is not None:
        ok &= xcorr['tstart'] >= CxoTime(start).secs

    if stop is not None:
        ok &= xcorr['tstart'] <= CxoTime(stop).secs

    for key, val in kwargs.items():
        if isinstance(val, list):
            ok &= np.isin(xcorr[key], val)
        else:
            ok &= xcorr[key] == val

    return xcorr[ok]


def plot_corr_srcs(xcorr, medians, start, stop, ms=1, **kwargs):
    """Plot astrometry offsets and return a new table with offsets included.
    """
    x = np.zeros((2, len(medians)))
    x[0] = medians['start'].frac_year
    x[1] = medians['stop'].frac_year
    sel = (medians['stop'] <= stop) & (medians['start'] > start)

    fig, (ax0, ax1) = plt.subplots(nrows=2, ncols=1, figsize=(8, 6))
    ax0.plot(xcorr['year'], xcorr['dy'], '.', ms=ms)
    ax0.plot(x, np.tile(medians['dy'], (2, 1)), '-', color='C0')
    ax0.plot(x[:, sel], np.tile(medians['dy'][sel], (2, 1)), '-', color='C1')
    ax1.plot(xcorr['year'], xcorr['dz'], '.', ms=ms)
    ax1.plot(x, np.tile(medians['dz'], (2, 1)), '-', color='C0')
    ax1.plot(x[:, sel], np.tile(medians['dz'][sel], (2, 1)), '-', color='C1')
    ax0.set_ylim(-2, 2)
    ax1.set_ylim(-2, 2)
    ax1.set_xlabel('Year')
    ax0.set_ylabel(r'$\Delta$y (arcsec)')
    ax1.set_ylabel(r'$\Delta$z (arcsec)')
    ax0.set_title(f'X-ray - Catalog astrometry')


def get_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--start',
        type=CxoTime,
        default=CxoTime('1999:001'),
        help="Only consider time bins starting after this time.")
    parser.add_argument(
        '--stop',
        type=CxoTime,
        default=CxoTime(),
        help=(
            "Only consider time bins ending before this time. "
            "The actual stop time is the minimum of this and the time of the last observation."))
    parser.add_argument(
        '--out',
        type=Path,
        help="JSON filename where to write offsets.")
    parser.add_argument(
        '--bins-per-year',
        default=2,
        help="Number of bins per year. Default=2")
    parser.add_argument(
        '--snr',
        default=5,
        help="Minimum SNR to consider an x-ray source")
    parser.add_argument(
        '--plot', action='store_true',
        default=False,
        help="Show plots")
    parser.add_argument(
        '--stdout', action='store_true',
        default=False,
        help="Write offsets to stdout.")
    return parser


def main():
    args = get_parser().parse_args()

    matches = db.get_cross_matches()
    matches = filter_xcorr(matches, snr=args.snr)

    args.stop = min(args.stop, np.max(matches['time']))

    start = CxoTime('1999:001')

    matches['year'] = CxoTime(matches['time']).decimalyear

    ok = matches['year'] > CxoTime(start).decimalyear
    bad_targets = ['RW Aur', 'Tau Boo', '70 OPH', '16 Cyg', 'M87', 'Orion', 'HD 97950' 'HD4915']
    bad_targets = [x.replace(' ', '').lower() for x in bad_targets]
    for ii, target in enumerate(matches['target']):
        target = target.replace(' ', '').lower()
        for bad_target in bad_targets:
            if target.startswith(bad_target):
                ok[ii] = False
    matches['ok'] = ok

    years_dy_binned, dy_binned, dy_bins = \
        binned_median(matches['year'], matches['dy'], bins_per_year=args.bins_per_year)
    years_dz_binned, dz_binned, dz_bins = \
        binned_median(matches['year'], matches['dz'], bins_per_year=args.bins_per_year)

    assert np.all(dy_bins == dy_bins)

    medians = Table()
    medians['start'] = CxoTime(dy_bins[:-1], format='frac_year')
    medians['stop'] = CxoTime(dy_bins[1:], format='frac_year')
    medians['years_dy'] = years_dy_binned
    medians['dy'] = dy_binned
    medians['years_dz'] = years_dz_binned
    medians['dz'] = dz_binned

    matches['dy_sub'] = (
        matches['dy'] - np.interp(x=matches['year'], xp=medians['years_dy'], fp=medians['dy'])
    )
    matches['dz_sub'] = (
        matches['dz'] - np.interp(x=matches['year'], xp=medians['years_dz'], fp=medians['dz'])
    )

    shifts = [
        {
            'tstart': row['start'].isot,
            'tstop': row['stop'].isot,
            'dy': row['dy'],
            'dz': row['dz'],
        } for row in medians[(medians['stop'] <= args.stop) & (medians['start'] > args.start)]
    ]

    if args.out:
        with open(args.out, 'w') as fh:
            json.dump(shifts, fh, indent=2)
    if args.stdout:
        print(json.dumps(shifts, indent=2))

    if args.plot:
        plot_corr_srcs(matches, medians, args.start, args.stop)
        plt.show()


if __name__ == '__main__':
    main()
