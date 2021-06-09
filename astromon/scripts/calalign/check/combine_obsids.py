#!/usr/bin/env python


from pathlib import Path
import numpy as np
from tqdm import tqdm


from astropy.table import Table, vstack, join, join_skycoord
from astropy.io import fits, ascii
from astropy.coordinates import SkyCoord
from astropy import units as u

from chandra_aca.transform import radec_to_yagzag
from Quaternion import Quat

from astromon import db


from astromon.utils import Ciao


DBFILE = '/proj/sot/ska/jgonzalez/aca_cal_align/update_2022-feb/check/ASTROMON_table.h5'

IMAGE_CENTER = {
    'acis-i': (4096.5, 4096.5),
    'acis-s': (4096.5, 4096.5),
    'hrc-i': (16384.5, 16484.5),
    'hrc-s': (32768.5, 32768.5),
}
IMG_CENTER_LABEL = {
    'acis-i': ('TCRPX11', 'TCRPX12'),
    'acis-s': ('TCRPX11', 'TCRPX12'),
    'hrc-i': ('TCRPX18', 'TCRPX19'),
    'hrc-s': ('TCRPX21', 'TCRPX22'),
}


def do_events():
    print('doing events')
    reprodir = Path('reprodata')
    obsdirs = sorted(reprodir.glob('*/*'))
    result = []
    for obsdir in tqdm(obsdirs[::-1]):
        CIAO = Ciao(workdir=obsdir)
        obsid = int(obsdir.name)

        obspar_files = list(obsdir.glob('archive/obspar/*obs0a.par'))
        if not obspar_files:
            print(f'{obsid} skipped (no obspar)')
            continue
        obsid_info = {
            r[0]: r[3] for r in ascii.read(obspar_files[0])
        }
        det = obsid_info['detector'].lower()
        # print(obsid, det)

        orig_asol_files = list((obsdir / 'archive' / 'uncorr' / 'out1').glob('pcad*asol*'))
        if len(orig_asol_files) < 1:
            print(f'{obsid} skipped (no asol)')
            continue
        orig_asol_file = orig_asol_files[0]
        corr_asol_file = list(obsdir.glob('archive/ASP*/out1/pcad*asol*'))[0]

        orig_hdu = fits.open(orig_asol_file)
        q = Quat([orig_hdu[1].header['RA_PNT'],
                  orig_hdu[1].header['DEC_PNT'],
                  orig_hdu[1].header['ROLL_PNT']])
        corr_hdu = fits.open(corr_asol_file)
        # Find the original and new yag and zag positions for
        # all of the aspect solution samples.  Deltas used later.
        orig_ys, orig_zs = radec_to_yagzag(
            orig_hdu[1].data['ra'],
            orig_hdu[1].data['dec'],
            q
        )
        new_ys, new_zs = radec_to_yagzag(
            corr_hdu[1].data['ra'],
            corr_hdu[1].data['dec'],
            q
        )
        orig_files = list(obsdir.glob('reproject/*uncorr_evt2*'))
        if len(orig_files) < 1:
            print(f'{obsid} skipped (no sources)')
            continue
        orig_file = orig_files[0]
        corr_file = list(obsdir.glob('reproject/*_corr_evt2*'))[0]
        orig_evt_hdu = fits.open(orig_file)
        corr_evt_hdu = fits.open(corr_file)
        # Find the delta sky positions for the events
        d_sky_x = corr_evt_hdu[1].data['x'] - orig_evt_hdu[1].data['x']
        d_sky_y = corr_evt_hdu[1].data['y'] - orig_evt_hdu[1].data['y']
        # Also find the radial change to find the point with the largest offset
        d_sky_r = np.sqrt(d_sky_x**2 + d_sky_y**2)
        sky_0 = {
            v: orig_evt_hdu[1].header[f'TCRPX{int(k[-2:])}']
            for k, v in orig_evt_hdu[1].header.items() if k.startswith('TTYPE') and v in ('x', 'y')
        }
        sky_x_0 = sky_0['x']
        sky_y_0 = sky_0['y']

        # Transform the sky stddev and means to ra/dec by applying them as
        # as deltas from the center position and converting with dmcoords.
        # this is clunky
        CIAO('punlearn', 'dmcoords')
        CIAO(
            'dmcoords',
            orig_file, orig_asol_file,
            x=sky_x_0 + np.std(d_sky_x),
            y=sky_y_0,
            option='sky',
            celfmt='deg',
        )
        stdpos_x = {k: float(CIAO.pget('dmcoords', k)) for k in ['ra', 'dec']}
        CIAO(
            'dmcoords',
            orig_file, orig_asol_file,
            x=sky_x_0,
            y=sky_y_0 + np.std(d_sky_y),
            option='sky',
            celfmt='deg',
        )
        stdpos_y = {k: float(CIAO.pget('dmcoords', k)) for k in ['ra', 'dec']}
        # Find the new and old positions of the events with the largest and smallest
        # radial offsets
        max_dpos_n = CIAO(
            'dmcoords',
            orig_file, orig_asol_file,
            x=corr_evt_hdu[1].data['x'][d_sky_r == np.max(d_sky_r)][0],
            y=corr_evt_hdu[1].data['y'][d_sky_r == np.max(d_sky_r)][0],
            option='sky',
            celfmt='deg',
        )
        max_dpos_n = {k: float(CIAO.pget('dmcoords', k)) for k in ['ra', 'dec']}
        CIAO(
            'dmcoords',
            orig_file, orig_asol_file,
            x=orig_evt_hdu[1].data['x'][d_sky_r == np.max(d_sky_r)][0],
            y=orig_evt_hdu[1].data['y'][d_sky_r == np.max(d_sky_r)][0],
            option='sky',
            celfmt='deg',
        )
        max_dpos_o = {k: float(CIAO.pget('dmcoords', k)) for k in ['ra', 'dec']}
        CIAO(
            'dmcoords',
            orig_file, orig_asol_file,
            x=corr_evt_hdu[1].data['x'][d_sky_r == np.min(d_sky_r)][0],
            y=corr_evt_hdu[1].data['y'][d_sky_r == np.min(d_sky_r)][0],
            option='sky',
            celfmt='deg',
        )
        min_dpos_n = {k: float(CIAO.pget('dmcoords', k)) for k in ['ra', 'dec']}
        CIAO(
            'dmcoords',
            orig_file, orig_asol_file,
            x=orig_evt_hdu[1].data['x'][d_sky_r == np.min(d_sky_r)][0],
            y=orig_evt_hdu[1].data['y'][d_sky_r == np.min(d_sky_r)][0],
            option='sky',
            celfmt='deg',
        )
        min_dpos_o = {k: float(CIAO.pget('dmcoords', k)) for k in ['ra', 'dec']}
        # Transform the "mean" offset into ra / dec using the center point as a reference
        CIAO(
            'dmcoords',
            orig_file, orig_asol_file,
            x=sky_x_0 + np.mean(d_sky_x),
            y=sky_y_0 + np.mean(d_sky_y),
            option='sky',
            celfmt='deg',
        )
        mean_dpos = {k: float(CIAO.pget('dmcoords', k)) for k in ['ra', 'dec']}
        CIAO(
            'dmcoords',
            orig_file, orig_asol_file,
            x=sky_x_0,
            y=sky_y_0,
            option='sky',
            celfmt='deg',
        )
        center_pos = {k: float(CIAO.pget('dmcoords', k)) for k in ['ra', 'dec']}
        stdxy, stdxz = radec_to_yagzag(stdpos_x['ra'], stdpos_x['dec'], q)
        stdyy, stdyz = radec_to_yagzag(stdpos_y['ra'], stdpos_y['dec'], q)
        orig_y, orig_z = radec_to_yagzag(center_pos['ra'], center_pos['dec'], q)
        new_y, new_z = radec_to_yagzag(mean_dpos['ra'], mean_dpos['dec'], q)
        max_new_y, max_new_z = radec_to_yagzag(max_dpos_n['ra'], max_dpos_n['dec'], q)
        max_orig_y, max_orig_z = radec_to_yagzag(max_dpos_o['ra'], max_dpos_o['dec'], q)
        min_new_y, min_new_z = radec_to_yagzag(min_dpos_n['ra'], min_dpos_n['dec'], q)
        min_orig_y, min_orig_z = radec_to_yagzag(min_dpos_o['ra'], min_dpos_o['dec'], q)
        obsinfo = {
            'det': det,
            'date-obs': orig_evt_hdu[1].header['DATE-OBS'],
            'obsid': obsid,
            'ra': float(orig_hdu[1].header['RA_PNT']),
            'dec': float(orig_hdu[1].header['DEC_PNT']),
            'roll': float(orig_hdu[1].header['ROLL_PNT']),
            'sol_ra_old_mean': np.mean(orig_hdu[1].data['ra']),
            'sol_dec_old_mean': np.mean(orig_hdu[1].data['dec']),
            'sol_ra_new_mean': np.mean(corr_hdu[1].data['ra']),
            'sol_dec_new_mean': np.mean(corr_hdu[1].data['dec']),
            'sol_ra_change_mean': np.mean(corr_hdu[1].data['ra'] - orig_hdu[1].data['ra']),
            'sol_dec_change_mean': np.mean(corr_hdu[1].data['dec'] - orig_hdu[1].data['dec']),
            'sol_y_mean': np.mean(new_ys - orig_ys),
            'sol_y_max': np.max(new_ys - orig_ys),
            'sol_y_min': np.min(new_ys - orig_ys),
            'sol_y_std': np.std(new_ys - orig_ys),
            'sol_z_mean': np.mean(new_zs - orig_zs),
            'sol_z_max': np.max(new_zs - orig_zs),
            'sol_z_min': np.min(new_zs - orig_zs),
            'sol_z_std': np.std(new_zs - orig_zs),
            'evt_y_mean': (new_y - orig_y),
            'evt_z_mean': (new_z - orig_z),
            'evt_y_max': (max_new_y - max_orig_y),
            'evt_z_max': (max_new_z - max_orig_z),
            'evt_y_min': (min_new_y - min_orig_y),
            'evt_z_min': (min_new_z - min_orig_z),
            'evt_stdx_y': (stdxy - orig_y),
            'evt_stdx_z': (stdxz - orig_z),
            'evt_stdy_y': (stdyy - orig_y),
            'evt_stdy_z': (stdyz - orig_z)
        }
        result.append(obsinfo)
    t = Table(result)
    t.write('obs_data.fits', overwrite=True)


def do_sources():
    print('doing sources')
    matches = db.get_cross_matches(dbfile=DBFILE)
    reprodir = Path('reprodata')
    obsdirs = sorted(reprodir.glob('*/*'))
    result = []
    for obsdir in tqdm(obsdirs[::-1]):
        obsid = int(obsdir.name)
        # print(obsid)
        new_file = obsdir / 'reproject' / 'new_src2.fits'
        old_file = obsdir / 'reproject' / 'old_src2.fits'
        obspar_files = list(obsdir.glob('archive/obspar/*obs0a.par'))
        if not new_file.exists() or not old_file.exists():
            print(f'{obsid} skipped (no sources)')
            continue
        if not obspar_files:
            print(f'{obsid} skipped (no obspar)')
            continue
        new = Table(fits.open(new_file)[1].data)
        old = Table(fits.open(old_file)[1].data)

        obsid_info = {
            r[0]: r[3] for r in ascii.read(obspar_files[0])
        }
        q = Quat(equatorial=(obsid_info['ra_pnt'], obsid_info['dec_pnt'], obsid_info['roll_pnt']))
        new['Y'], new['Z'] = radec_to_yagzag(new['RA'], new['DEC'], q)
        old['Y'], old['Z'] = radec_to_yagzag(old['RA'], old['DEC'], q)
        matches['coord'] = SkyCoord(matches['x_ra'], matches['x_dec'], unit='deg')
        new['coord'] = SkyCoord(new['RA'], new['DEC'], unit='deg')
        old['coord'] = SkyCoord(old['RA'], old['DEC'], unit='deg')
        t = join(
            matches[matches['obsid'] == obsid], new,
            join_funcs={'coord': join_skycoord(2 * u.arcsec)}
        )
        t.rename_column('coord_1', 'coord')
        t = join(
            t, old,
            keys=['coord'],
            join_funcs={'coord': join_skycoord(2 * u.arcsec)},
            table_names=['new', 'old']
        )
        # c_y_angle/c_z_angle are the yag/zag of the counterpart's ra/dec at the reference attitude
        t['dy_new'] = t['Y_new'] - t['c_y_angle']  # dy with new calalign
        t['dz_new'] = t['Z_new'] - t['c_z_angle']  # dz with new calalign
        t['dy_old'] = t['Y_old'] - t['c_y_angle']  # dy with previous calalign
        t['dz_old'] = t['Z_old'] - t['c_z_angle']  # dz with previous calalign
        t['ra_pnt'] = float(obsid_info['ra_pnt'])
        t['dec_pnt'] = float(obsid_info['dec_pnt'])
        t['roll_pnt'] = float(obsid_info['roll_pnt'])
        result.append(t[[
            'obsid', 'date_obs', 'tstart', 'ra_pnt', 'dec_pnt', 'roll_pnt', 'dy', 'dz', 'snr',
            'c_id', 'c_ra', 'c_dec', 'c_y_angle', 'c_z_angle',
            'x_id', 'x_ra', 'x_dec', 'x_y_angle', 'x_z_angle',
            'RA_old', 'DEC_old', 'SNR_old', 'Y_old', 'Z_old', 'dy_old', 'dz_old',
            'RA_new', 'DEC_new', 'SNR_new', 'Y_new', 'Z_new', 'dy_new', 'dz_new',
        ]])
    result = vstack(result)
    result.write('source_data.fits', overwrite=True)


def main():
    do_sources()
    do_events()


if __name__ == '__main__':
    main()
