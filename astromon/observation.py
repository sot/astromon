#!/usr/bin/env python


import os
# import sys
import logging
import argparse
import tempfile
import subprocess
from pathlib import Path

try:
    from ciao_contrib import runtool
    from ciao_contrib.runtool import make_tool, dmstat
    import paramio
    with_ciao = True
except ModuleNotFoundError:
    with_ciao = False

from astromon.utils import ciao_context_function, logging_call_decorator, chdir

from astropy import table
from astropy.io import fits, ascii

import numpy as np

import Ska.arc5gl

logger = logging.getLogger('astromon')


_multi_obi_obsids = [
    82, 108, 279, 380, 400, 433, 800, 861, 897, 906, 943, 1411, 1431,
    1456, 1561, 1578, 2010, 2042, 2077, 2365, 2783, 3057, 3182, 3764,
    4175, 60879, 60880, 62249, 62264, 62796
]


def download_chandra_obsids(obsids, filetypes):
    r = subprocess.run(['download_chandra_obsid', '-t', '-q'], stdout=subprocess.PIPE)
    available_types = r.stdout.decode().strip().split(':')[-1].split()
    exclude = [t for t in available_types if t not in filetypes]
    r = subprocess.run(
        ['download_chandra_obsid', ','.join(obsids), '--exclude', ','.join(exclude)],
        stdout=subprocess.PIPE
    )


class Observation:
    def __init__(self, obsid, workdir=None, source='arc5gl'):
        self.tmp = tempfile.TemporaryDirectory() if workdir is None else None
        self.workdir = Path(self.tmp.name if workdir is None else workdir).expanduser()
        self.obsid = obsid
        self._source = source
        logger.info(f'Starting {self}. Context will be {self.workdir / self.obsid}')
        self._rebin = False

    def __str__(self):
        return f'OBSID {self.obsid}'

    @ciao_context_function
    def get_obsid_info(self):
        obsid_info = {'obsid': self.obsid}
        obsid_info.update({
            k: None for k in
            ['multi-obi', 'obs_mode', 'grating', 'instrument', 'read_mode', 'dtycycle', 'error']
        })

        obsid_info['multi-obi'] = self.obsid in _multi_obi_obsids

        msk = list((self.workdir / self.obsid / 'secondary').glob("*msk1.fits*"))
        if len(msk) != 1:
            obsid_info['error'] = f'n_mask == {len(msk)}'
            return obsid_info

        try:
            def get_key_value(hdul, k):
                for hdu in hdul:
                    if k in hdu.header.keys():
                        return hdu.header[k]
            hdul = fits.open(str(msk[0]))
            obsid_info.update({
                'obs_mode': get_key_value(hdul, "OBS_MODE"),
                'grating': get_key_value(hdul, "GRATING"),
                'instrument': get_key_value(hdul, "INSTRUME"),
                'date': get_key_value(hdul, "DATE-OBS"),
                'tstart': get_key_value(hdul, "TSTART"),
                'tstop': get_key_value(hdul, "TSTOP"),
                'sim_x': get_key_value(hdul, "SIM_X"),
                'sim_y': get_key_value(hdul, "SIM_Y"),
                'sim_Z': get_key_value(hdul, "SIM_Z"),
            })

            if obsid_info['instrument'] == 'ACIS':
                obsid_info.update({
                    'read_mode': get_key_value(hdul, "READMODE"),
                    'dtycycle': get_key_value(hdul, "DTYCYCLE"),
                })
        except Exception as e:
            obsid_info['error'] = f"error reading mask file: {e}"

        return obsid_info

    def get_obspar(self):
        (self.workdir / self.obsid).mkdir(exist_ok=True, parents=True)
        obspar_file = list((self.workdir / self.obsid).glob('*obs0*'))
        if not obspar_file:
            self.download(['obspar'])
        obspar_file = list((self.workdir / self.obsid).glob('*obs0*'))
        assert len(obspar_file) > 0
        obspar_file = str(obspar_file[0])
        t = ascii.read(obspar_file)
        self._obsid_info = {
            r[0]: r[3] for r in t
        }
        self._obsid_info['instrument'] = self._obsid_info['instrume'].lower()
        self._obsid_info['obsid'] = int(self.obsid)
        return(self._obsid_info)

    @logging_call_decorator
    @ciao_context_function
    def _download_archive(self, ftypes):
        """
        Download data from chandra public archive
        """
        if not self.workdir.exists():
            self.workdir.mkdir(exist_ok=True, parents=True)

        secondary = self.workdir / self.obsid / 'secondary'
        if secondary.exists():
            logger.debug(f'{secondary} exists, skipping download')
            return

        repro = self.workdir / self.obsid / 'repro'
        if repro.exists():
            logger.debug(f'{repro} exists, skipping download')
            return

        with chdir(self.workdir):
            good = download_chandra_obsids([self.obsid], filetypes=ftypes)

        if good[0] is not True:
            raise RuntimeError(f"Can't download {self.obsid}")

    def _download_arc5gl(self, ftypes):
        """
        Download data from chandra public archive
        """
        for ftype in ftypes:
            if ftype == 'obspar':
                src, dest = 'obspar', '.'
            else:
                locs = self._archive_file_locations()
                if ftype not in locs:
                    logger.error(f'{ftype=} skipped because it is not in known locations')
                    continue
                src, dest = locs[ftype]
            logger.info(f'  {ftype=}')
            dest_files = list((self.workdir / self.obsid / dest).glob(f'*{ftype}*'))
            if len(dest_files):
                logger.info(f'    skipping download of *{ftype}*')
                continue
            logger.info(f'    {src} -> {dest}')
            (self.workdir / self.obsid / dest).mkdir(exist_ok=True, parents=True)
            with chdir(self.workdir / self.obsid / dest):
                arc5gl = Ska.arc5gl.Arc5gl()
                arc5gl.sendline(f'obsid={self.obsid}')
                arc5gl.sendline(f'get {src}')
                del arc5gl

    @logging_call_decorator
    def download(self, ftypes=None):
        if ftypes is None:
            ftypes = ["evt2", "asol", "msk"]
        if self._source and self._source == 'archive':
            return self._download_archive(ftypes)
        return self._download_arc5gl(ftypes)

    @logging_call_decorator
    @ciao_context_function
    def repro(self):
        """
        Reprocess data
        """
        images = self.workdir / self.obsid / 'images'
        if images.exists():
            # Skip repro, already done
            return

        chandra_repro = make_tool("chandra_repro")
        verb = chandra_repro(str(self.workdir / self.obsid), outdir="", cleanup=True, clobber=True)
        if verb:
            logger.info(f'repro {verb}')

    @logging_call_decorator
    @ciao_context_function
    def make_images(self, evt=None):
        'Create image. Also creates the exposure map and psfmap'

        is_hrc = self.get_obsid_info()['instrument'].lower() == 'hrc'

        if evt is None:
            evtfiles = list((self.workdir / self.obsid / 'primary').glob('*_evt2_filtered.fits*'))
            assert len(evtfiles) == 1, f"Expected 1 evt file, there are {len(evtfiles)}"
            evt = evtfiles[0]

        outdir = self.workdir / self.obsid / 'images'
        if outdir.exists():
            logger.info(f'  directory {outdir} exists, skipping')
            return

        outdir.mkdir(exist_ok=True)

        fimg = make_tool("fluximage")
        fimg.infile = str(evt)
        fimg.outroot = str(outdir / self.obsid)
        if is_hrc is True:
            fimg.bands = "wide"
            fimg.binsize = 4 if self._rebin else 1
        else:
            fimg.bands = "broad"
            fimg.binsize = 1
        fimg.psfecf = 0.9
        fimg.background = "none"
        verb = fimg(clobber=True, parallel=False)
        if verb:
            logger.info(f'images {verb}')

    @logging_call_decorator
    @ciao_context_function
    def run_wavdetect(self, edition, skip_exist=False, scales="1.4 2 4 8 16 32"):
        'Run wavdetect'

        is_hrc = self.get_obsid_info()['instrument'].lower() == 'hrc'

        imgdir = self.workdir / self.obsid / 'images'
        detdir = self.workdir / self.obsid / 'sources'
        detdir.mkdir(parents=True, exist_ok=True)

        band = "wide" if is_hrc else "broad"

        # inputs
        img = self.obsid + "_" + band + "_thresh.img"
        exp = self.obsid + "_" + band + "_thresh.expmap"
        psf = self.obsid + "_" + band + "_thresh.psfmap"

        # outputs
        root = f"{self.obsid}_{edition}"
        src = root + ".src"
        cel = root + ".cell"
        nbk = root + ".nbkg"
        rec = root + ".recon"

        if (detdir / src).exists() and skip_exist:
            return

        wavdetect = make_tool("wavdetect")
        wavdetect.infile = str(imgdir / img)
        wavdetect.psffile = str(imgdir / psf)  # PSF
        wavdetect.expfile = str(imgdir / exp)  # exposure map

        wavdetect.outfile = str(detdir / src)
        wavdetect.scellfile = str(detdir / cel)
        wavdetect.imagefile = str(detdir / nbk)
        wavdetect.defnbkgfile = str(detdir / rec)

        wavdetect.scales = scales

        verb = wavdetect(clobber=True)
        if verb:
            logger.info(f'wavdetect {verb}')

        # ~ import subprocess as subprocess
        # ~ subprocess.run("gzip -f {}".format(wavdetect.scellfile).split(" "))
        # ~ subprocess.run("gzip -f {}".format(wavdetect.imagefile).split(" "))
        # ~ subprocess.run("gzip -f {}".format(wavdetect.defnbkgfile).split(" "))

    def run_celldetect(self, snr=3):
        # Find sources in the small field
        imgdir = self.workdir / self.obsid / 'images'
        is_hrc = self.get_obsid_info()['instrument'].lower() == 'hrc'
        band = "wide" if is_hrc else "broad"
        celldetect = make_tool("celldetect")
        celldetect(
            imgdir / f"{self.obsid}_{band}_thresh.img",
            self.workdir / self.obsid / 'sources' / f"{self.obsid}_celldetect.src",
            psffile=imgdir / f"{self.obsid}_{band}_thresh.psfmap",  # either this or set fixedcell=
            thresh=snr,
            maxlogicalwindow=2048,
            clobber=True
        )

    @logging_call_decorator
    @ciao_context_function
    def filter_events(self, radius=180, psfratio=1):  # radius in arcsec
        is_hrc = self.get_obsid_info()['instrument'].lower() == 'hrc'
        # I'm using a fixed pixel size of 0.5 arcsec, but this might need fixing
        pixel = 1 if is_hrc else 0.5
        if self._rebin and is_hrc:
            pixel *= 2
        try:
            evt = list((self.workdir / self.obsid / 'primary').glob('*evt2.fits*'))[0]
        except Exception:
            raise Exception(f'evt2 file not found   ') from None

        evt2 = str(evt).replace('evt2', 'evt2_filtered')

        if Path(evt2).exists():
            logger.info(f'  file {evt2} exists, skipping')
            return

        runtool.dmkeypar(evt, 'RA_PNT')
        ra = runtool.dmkeypar.value
        runtool.dmkeypar(evt, 'DEC_PNT')
        dec = runtool.dmkeypar.value
        runtool.dmcoords(evt, op='cel', celfmt='deg', ra=ra, dec=dec)
        filters = [
            f'(x,y)=circle({runtool.dmcoords.x},{runtool.dmcoords.y},{radius/pixel})'
        ]
        filters = ','.join(filters)
        runtool.dmcopy(f'{evt}[{filters}]', evt2)
        return evt2

    @logging_call_decorator
    @ciao_context_function
    def filter_sources(self, radius=180, psfratio=1):  # radius in arcsec
        # I'm using a fixed pixel size of 0.5 arcsec, but this might need fixing
        is_hrc = self.get_obsid_info()['instrument'].lower() == 'hrc'
        # I'm using a fixed pixel size of 0.5 arcsec, but this might need fixing
        pixel = 1 if is_hrc else 0.5
        if self._rebin and is_hrc:
            pixel *= 2
        try:
            evt = list((self.workdir / self.obsid / 'primary').glob('*evt2.fits*'))[0]
        except Exception:
            raise Exception(f'evt2 file not found   ') from None

        src = self.workdir / self.obsid / 'sources' / f'{self.obsid}_baseline.src'
        if not src.exists():
            raise Exception(f'src file not found   ')
        src2 = str(src).replace('baseline', 'filtered')

        runtool.dmkeypar(evt, 'RA_PNT')
        ra = runtool.dmkeypar.value
        runtool.dmkeypar(evt, 'DEC_PNT')
        dec = runtool.dmkeypar.value
        runtool.dmcoords(evt, op='cel', celfmt='deg', ra=ra, dec=dec)
        filters = [
            f'psfratio=:{psfratio}',
            f'(x,y)=circle({runtool.dmcoords.x},{runtool.dmcoords.y},{radius/pixel})'
        ]
        filters = ','.join(filters)
        runtool.dmcopy(f'{src}[{filters}]', src2)
        return src2

    @logging_call_decorator
    @ciao_context_function
    def calculate_centroids(self):
        src_hdus = fits.open(self.workdir / self.obsid / 'sources' / f'{self.obsid}_baseline.src')
        img = self.workdir / self.obsid / 'images' / f'{self.obsid}_wide_flux.img'

        assert img.exists()
        src = table.Table(src_hdus[1].data)
        result = []
        for row in src:
            x, y = row[['X', 'Y']]
            r = row['R'].max()
            dmstat.punlearn()
            stdout = dmstat(
                f'{img}[sky=circle({x},{y},{r})]',
                centroid='yes',
                # clip='yes'
            )
            logger.info(stdout)
            dmstat.write_params()
            x, y = np.array(paramio.pget('dmstat', 'out_cntrd_phys').split(',')).astype(float)
            result.append([x, y])
        return np.array(result, dtype=[('x', float), ('y', float)])

    @logging_call_decorator
    def process(self):
        'Main routine to process a single obsid'

        self.download()
        obsid_info = self.get_obsid_info()
        if obsid_info['error']:
            raise Exception(obsid_info['error'])

        ok = (
            not obsid_info['error']
            and not obsid_info['multi-obi']
            and obsid_info['obs_mode'] == 'POINTING'
            # and obsid_info['grating'] == 'NONE'
            and (
                obsid_info['instrument'] == 'HRC'
                or (
                    obsid_info['instrument'] == 'ACIS'
                    and obsid_info['read_mode'] == "TIMED"
                    and obsid_info['dtycycle'] == 0
                )
            )
        )
        if not ok:
            raise Exception('Skipping')

        # Repro
        # repro(self.obsid)

        # Analysis
        self.filter_events()
        self.make_images()
        self.run_wavdetect("baseline", skip_exist=True)
        self.run_celldetect()

    @logging_call_decorator
    def get_sources(self, version='celldetect'):
        """
        Returns a table of sources formatted for the astromon_xray_source SQL table
        """
        from Quaternion import Quat
        from chandra_aca.transform import radec_to_yagzag
        obspar = self.get_obspar()
        q = Quat(equatorial=(obspar['ra_pnt'], obspar['dec_pnt'], obspar['roll_pnt']))

        hdu_list = fits.open(
            self.workdir / self.obsid / 'sources' / f'{self.obsid}_{version}.src')
        sources = table.Table(hdu_list[1].data)

        if len(sources) == 0:
            raise Exception('No x-ray sources found')

        # SNR is there only in the celldetect case, and this is how it should be calculated
        # (https://cxc.harvard.edu/ciao/download/doc/detect_manual/cell_theory.html)
        # but this gives a slightly different result from what celldetect actually gives
        # wavdetect does not calculate this.
        C = sources['NET_COUNTS']
        B = sources['BKG_COUNTS']
        sigma_c = sources['NET_COUNTS_ERR']
        sigma_b = sources['BKG_COUNTS_ERR']
        sources['SNR'] = (C - B) / np.sqrt(sigma_c**2 + sigma_b**2)

        columns = [
            c for c in zip(
                ['RA', 'DEC', 'COMPONENT', 'NET_COUNTS', 'SNR', 'PSFRATIO'],
                ['ra', 'dec', 'id', 'net_counts', 'snr', 'psfratio']
            ) if c[0] in sources.colnames
        ]
        sources.rename_columns(*list(zip(*columns)))

        sources['obsid'] = self.obsid
        sources['y_angle'], sources['z_angle'] = radec_to_yagzag(sources['ra'], sources['dec'], q)
        sources['r_angle'] = np.sqrt(sources['y_angle']**2 + sources['z_angle']**2)

        # calculate the distance to the closest source
        src1 = sources.as_array()[None]
        src2 = sources.as_array()[:, None]
        i = np.diag([np.inf] * len(sources))
        distance = np.sqrt((src1['y_angle'] - src2['y_angle'])**2
                           + (src1['z_angle'] - src2['z_angle'])**2) + i
        sources['near_neighbor_dist'] = np.min(distance, axis=0)

        cols = [
            'obsid',
            'id',
            # 'name',  # not setting it for now. Is it used?
            'ra',
            'dec',
            'net_counts',
            'y_angle',
            'z_angle',
            'r_angle',
            'snr',
            'near_neighbor_dist',
            # 'double_id',  # does not seem to be set in current astromon
            # 'status_id',  # does not seem to be set in current astromon
            'psfratio',
        ]
        return sources[cols]

    def _archive_file_locations(self):
        instrument = self.get_obspar()['instrument']
        return {
            'obspar': ('obspar', '.'),
            'evt2': (f'{instrument}2{{evt2}}', 'primary'),
            'evt1': (f'{instrument}1{{evt1}}', 'secondary'),
            'fov': (f'{instrument}1{{fov}}', 'primary'),
            'msk': (f'{instrument}1{{msk}}', 'secondary'),
            'mtl': (f'{instrument}1{{mtl}}', 'secondary'),
            'bpix': (f'{instrument}1[*bpix*]', 'primary'),
            'flt': (f'{instrument}1[*flt*]', 'secondary'),
            'stat': (f'{instrument}1[*stat*]', 'secondary'),
            'asol': (f'asp1[*{self.obsid}*asol*]', 'primary'),
            # 'dtf': f'',
            # 'pbk': f'{instrument}0[*pbk0*]',
            # 'bias': f'{instrument}0[*bias*]',
        }


def get_parser():
    """
    Get the argument parser
    """
    ska = Path(os.environ['SKA'])
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'obsid',
        help='OBSID(s) to process or a file with a list of OBSIDs',
        nargs='+'
    )
    parser.add_argument(
        '--workdir',
        help='Working directory. Default is to create a tmp directory',
        default=ska / 'data' / 'absolute_astrometry',
        type=Path
    )
    parser.add_argument(
        '--log-level',
        help='logging level',
        choices=['debug', 'info', 'warning'],
        default='info'
    )
    return parser


def main():
    """
    Main routine
    """

    from ciao_contrib._tools.taskrunner import TaskRunner

    args = get_parser().parse_args()

    # not using pyyaks because CIAO overrides the default logger class using logging.setLoggerClass
    logging.basicConfig(level=args.log_level.upper())

    logger.info(f'workdir: {args.workdir}')
    if args.workdir is None:
        tmp = tempfile.TemporaryDirectory()
        args.workdir = tmp.name
    Path(args.workdir).mkdir(exist_ok=True, parents=True)

    obsids = []
    for obsid in args.obsid:
        if os.path.exists(obsid) and os.path.isfile(obsid):
            with open(obsid) as fh:
                obsids += [l.strip() for l in fh.readlines()]
        else:
            obsids += obsid.split(',')

    taskrunner = TaskRunner()
    for obsid in obsids:
        observation = Observation(obsid, args.workdir)
        taskrunner.add_task(f"OBS_ID={obsid}", "", observation.process)

    # '9' because the archive limits number of concurrent connections
    # from same IP address (well, it did when it was 'ftp').
    taskrunner.run_tasks(processes=9)


if __name__ == '__main__':
    main()
