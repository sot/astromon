#!/usr/bin/env python

import json
import matplotlib.pyplot as plt

from Ska.Matplotlib import plot_cxctime
from astropy.table import Table


def plot_offset_change():
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(12, 8), num='offset_change')
    plt.sca(ax0)
    plot_cxctime(
        [OFFSETS['tstart'], OFFSETS['tstop']],
        [OFFSETS['dy'], OFFSETS['dy']],
        '-',
        color='tab:blue'
    )
    plot_cxctime(
        SOURCE['date_obs'],
        -(SOURCE['dy_new'] - SOURCE['dy_old']),
        '.',
        color='tab:blue'
    )
    ax0.set_ylabel(r'$-(dy_{new} - dy_{old})$')
    plt.grid()

    plt.sca(ax1)
    plot_cxctime(
        [OFFSETS['tstart'], OFFSETS['tstop']],
        [OFFSETS['dz'], OFFSETS['dz']],
        '-',
        color='tab:orange'
    )
    plot_cxctime(
        SOURCE['date_obs'],
        -(SOURCE['dz_new'] - SOURCE['dz_old']),
        '.',
        color='tab:orange'
    )
    ax1.set_ylabel(r'$-(dz_{new} - dz_{old})$')
    ax0.set_title('Change in X-ray Source Offsets')
    plt.grid()
    fig.savefig(f'{fig.get_label()}.png')


def plot_offsets():
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(12, 8), num='offsets_old')
    plt.sca(ax0)
    plot_cxctime(
        [OFFSETS['tstart'], OFFSETS['tstop']],
        [OFFSETS['dy'], OFFSETS['dy']],
        '-',
        color='k'
    )
    plot_cxctime(
        SOURCE['date_obs'], SOURCE['dy_old'], '.',
        markersize=10, color='k', label='old'
    )
    ax0.set_ylabel(r'$\Delta Y$')
    # plt.legend(loc='upper left')
    plt.ylim((-1.5, 1.5))
    plt.grid()

    plt.sca(ax1)
    plot_cxctime(
        [OFFSETS['tstart'], OFFSETS['tstop']],
        [OFFSETS['dz'], OFFSETS['dz']],
        '-',
        color='k'
    )
    plot_cxctime(
        SOURCE['date_obs'], SOURCE['dz_old'], '.',
        markersize=10, color='k', label='old'
    )
    ax1.set_ylabel(r'$\Delta Z$')
    # plt.legend(loc='upper left')
    ax0.set_title(f'X-ray Source Offsets (CALALIGN Version {OLD_VERSION})')
    plt.ylim((-1.5, 1.5))
    plt.grid()
    fig.savefig(f'{fig.get_label()}.png')

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(12, 8), num='offsets_new')
    plt.sca(ax0)
    plot_cxctime(
        [OFFSETS['tstart'][0], OFFSETS['tstop'][-1]],
        [0, 0],
        '-',
        color='tab:orange'
    )
    plot_cxctime(
        SOURCE['date_obs'], SOURCE['dy_new'], '.',
        markersize=10, color='tab:orange', label='new'
    )
    ax0.set_ylabel(r'$\Delta Y$')
    # plt.legend(loc='upper left')
    plt.ylim((-1.5, 1.5))
    plt.grid()

    plt.sca(ax1)
    plot_cxctime(
        [OFFSETS['tstart'][0], OFFSETS['tstop'][-1]],
        [0, 0],
        '-',
        color='tab:orange'
    )
    plot_cxctime(
        SOURCE['date_obs'], SOURCE['dz_new'], '.',
        markersize=10, color='tab:orange', label='new'
    )
    ax1.set_ylabel(r'$\Delta Z$')
    # plt.legend(loc='upper left')
    ax0.set_title(f'X-ray Source Offsets (CALALIGN Version {NEW_VERSION})')
    plt.grid()
    plt.ylim((-1.5, 1.5))
    fig.savefig(f'{fig.get_label()}.png')

    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(12, 8), num='offsets')
    plt.sca(ax0)
    plot_cxctime(
        [OFFSETS['tstart'], OFFSETS['tstop']],
        [OFFSETS['dy'], OFFSETS['dy']],
        '-',
        color='k'
    )
    plot_cxctime(
        SOURCE['date_obs'], SOURCE['dy_old'], '.',
        markersize=10, color='k', label=f'Version {OLD_VERSION}'
    )
    plot_cxctime(
        SOURCE['date_obs'], SOURCE['dy_new'], '.',
        markersize=10, color='tab:blue', label=f'Version {NEW_VERSION}'
    )
    ax0.set_ylabel(r'$\Delta Y$')
    plt.legend(loc='upper left')
    plt.grid()
    plt.ylim((-1.5, 1.5))

    plt.sca(ax1)
    plot_cxctime(
        [OFFSETS['tstart'], OFFSETS['tstop']],
        [OFFSETS['dz'], OFFSETS['dz']],
        '-',
        color='k'
    )
    plot_cxctime(
        SOURCE['date_obs'], SOURCE['dz_old'], '.',
        markersize=10, color='k', label=f'Version {OLD_VERSION}'
    )
    plot_cxctime(
        SOURCE['date_obs'], SOURCE['dz_new'], '.',
        markersize=10, color='tab:orange', label=f'Version {NEW_VERSION}'
    )
    ax1.set_ylabel(r'$\Delta Z$')
    plt.legend(loc='upper left')
    ax0.set_title('X-ray Source Offsets')
    plt.grid()
    plt.ylim((-1.5, 1.5))
    fig.savefig(f'{fig.get_label()}.png')


def plot_attitude_change():
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(12, 8), num='attitude_change')
    plt.sca(ax0)
    plot_cxctime(
        [OFFSETS['tstart'], OFFSETS['tstop']],
        [OFFSETS['dy'], OFFSETS['dy']],
        '-',
        color='tab:blue'
    )
    plot_cxctime(OBS['date-obs'], -OBS['sol_y_mean'], '.')
    ax0.set_ylabel(r'$-(Y_{sol\,new} - Y_{sol\,orig})$')
    plt.grid()
    plt.sca(ax1)
    plot_cxctime(
        [OFFSETS['tstart'], OFFSETS['tstop']],
        [OFFSETS['dz'], OFFSETS['dz']],
        '-',
        color='tab:blue'
    )
    plot_cxctime(OBS['date-obs'], -OBS['sol_z_mean'], '.')
    ax1.set_ylabel(r'$-(Z_{sol\,new} - Z_{sol\,orig})$')
    ax0.set_title('Change in Aspect Solution')
    plt.grid()
    fig.savefig(f'{fig.get_label()}.png')


def plot_event_change():
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(12, 8), num='event_change')
    plt.sca(ax0)
    plot_cxctime(
        [OFFSETS['tstart'], OFFSETS['tstop']],
        [OFFSETS['dy'], OFFSETS['dy']],
        '-',
        color='tab:blue'
    )
    plot_cxctime(OBS['date-obs'], -OBS['evt_y_mean'], '.')
    ax0.set_ylabel(r'$-(Y_{evt\,new} - Y_{evt\,orig})$')
    plt.grid()
    plt.sca(ax1)
    plot_cxctime(
        [OFFSETS['tstart'], OFFSETS['tstop']],
        [OFFSETS['dz'], OFFSETS['dz']],
        '-',
        color='tab:blue'
    )
    plot_cxctime(OBS['date-obs'], -OBS['evt_z_mean'], '.')
    ax1.set_ylabel(r'$-(Z_{evt\,new} - Z_{evt\,orig})$')
    ax0.set_title('Change in Events Y/Z')
    plt.grid()
    fig.savefig(f'{fig.get_label()}.png')


if __name__ == '__main__':
    with open('offsets.json') as fh:
        OFFSETS = Table(json.load(fh))

    SOURCE = Table.read('source_data.fits')
    OBS = Table.read('obs_data.fits')

    OLD_VERSION = 'N0009'
    NEW_VERSION = 'N0010'

    plot_offsets()
    plot_offset_change()
    plot_attitude_change()
    plot_event_change()
