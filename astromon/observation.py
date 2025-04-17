#!/usr/bin/env python


import argparse
import collections
import functools

# import sys
import logging
import os
import re
import shutil
import subprocess
import tempfile
import warnings
from pathlib import Path

import bs4
import chardet
import numpy as np
import regions
import requests
from astropy import table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii, fits
from astropy.wcs import WCS, FITSFixedWarning

from astromon.cache import Cache
from astromon.utils import Ciao, FlowException, chdir, logging_call_decorator

__all__ = ["Observation"]


logger = logging.getLogger("astromon")


_multi_obi_obsids = [
    82,
    108,
    279,
    380,
    400,
    433,
    800,
    861,
    897,
    906,
    943,
    1411,
    1431,
    1456,
    1561,
    1578,
    2010,
    2042,
    2077,
    2365,
    2783,
    3057,
    3182,
    3764,
    4175,
    60879,
    60880,
    62249,
    62264,
    62796,
]


ID_CATEGORY_MAP = {
    10: "Normal Stars and WD",
    20: "WD Binaries and CVs",
    30: "BH and NS Binaries",
    40: "Normal Galaxies",
    50: "Active Galaxies and Quasars",
    60: "Extragalactic Diffuse Emission & Surveys",
    70: "Galactic Diffuse Emission & Surveys",
    100: "Solar System and Misc",
    110: "SN, SNR, and Isolated NS",
    120: "Clusters of Galaxies",
    200: "Unknown",
}
"""
Mapping between observation ID and category names.
"""


CATEGORY_ID_MAP = collections.defaultdict(
    lambda: 200, {k.lower(): v for v, k in ID_CATEGORY_MAP.items()}
)
"""
Mapping between observation category names and numerical values.
"""


class Skipped(FlowException):
    """
    Exception class used to abort and silently skip processing an observation.
    """


class SkippedWithWarning(FlowException):
    """
    Exception class used to abort, issue a warning, and skip processing an observation.
    """


@functools.cache
def ciao_fails(ciao_prefix, workdir):
    try:
        Ciao(
            prefix=ciao_prefix,
            workdir=workdir,
            logger="astromon",
        )
        return ""
    except Exception as e:
        return f"CIAO could not be initialized: {e}"


class Observation:
    """
    Class to encapsulate calls to CIAO and arc5gl.

    Parameters
    ----------
    obsid: int
    workdir : pathlib.Path or str.
        Top-level working directory. Used to set PFILES and ASCDS_WORK_PATH and to download files.
        The actual working directory will be {workdir}/obs{int(obsid)//1000:02d}/{obsid}.
        If not given, the working directory will be a temporary directory.
    source: str
        One of 'arc5gl' (to get data using arc5gl) or 'archive' (to get data from chandra public
        archive using CIAO's download_chandra_obsid).
    archive_dir: pathlib.Path or str.
        Top-level archive directory. Used to permanently store the result files.
        The actual archive directory will be {archive_dir}/obs{int(obsid)//1000:02d}/{obsid}.
    ciao_prefix : str
        The location of CIAO.
    logger: logging.Logger.
        If not provided, the root logger is used.
    archive_regex: list of str.
        Files matching any regex in the list are archived.
    use_ciao: bool
        If this is False, there will be no calls to CIAO tools (and some methods will just fail).
    """

    def __init__(
        self,
        obsid,
        workdir=None,
        source="arc5gl",
        ciao_prefix=None,
        archive_dir=None,
        archive_regex=None,
        use_ciao=True,
    ):
        self.use_ciao = use_ciao or ciao_prefix
        self.ciao_prefix = ciao_prefix
        self._clear = workdir is None
        self.tmp = tempfile.TemporaryDirectory() if workdir is None else None
        self.obsid = str(obsid)
        subdir = f"obs{int(obsid) // 1000:02d}"
        self.workdir = (
            Path(self.tmp.name if workdir is None else workdir).expanduser()
            / subdir
            / self.obsid
        )
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.archive_dir = (
            (Path(archive_dir).expanduser() / subdir / self.obsid)
            if archive_dir
            else None
        )

        self.cache_dir = self.workdir / "cache"

        if archive_regex is None:
            self.archive_regex = [
                "*.par.gz",
                "cache",
                # '*evt2_filtered.fits.gz',
                "*_wide_flux.img",
                "*_broad_flux.img",
                "*_acis_streaks.txt",
                "*.src",
            ]
        else:
            self.archive_regex = archive_regex
        self._source = source
        logger.info(f"{self} starting. Context: {self.workdir}")
        self._rebin = False
        self.ciao_ = None
        # ciao is initialized only if needed, but we make this check at the very beginning
        # so the user gets a warning at the top if calling CIAO will likely fail
        if self.use_ciao and (msg := ciao_fails(ciao_prefix, workdir)):
            self.use_ciao = False
            logger.warning(msg)

    is_multi_obi = property(lambda self: int(self.get_obspar()["obsid"]) in _multi_obi_obsids)

    is_hrc = property(lambda self: self.get_obspar()["instrume"].lower() == "hrc")

    is_acis = property(lambda self: self.get_obspar()["instrume"].lower() == "acis")

    def create_archive_symlinks(self):
        """
        Create symlinks from working directory to the archive.

        Normally this should not be needed, but it might be convenient so one works only in the
        working directory. It is a bit slow (about 0.07 seconds, which translates to 5 minutes when
        creating instances for 5000 observations).
        """
        if self.archive_dir is not None:
            for file in self.archive_dir.glob("**/*"):
                if not file.is_dir():
                    dest = self.workdir / file.relative_to(self.archive_dir)
                    if not dest.exists():
                        dest.parent.mkdir(parents=True, exist_ok=True)
                        os.symlink(file, dest)

    def get_ciao(self):
        if self.use_ciao and self.ciao_ is None:
            try:
                self.ciao_ = Ciao(
                    prefix=self.ciao_prefix,
                    workdir=self.workdir / "param",
                    logger="astromon",
                )
            except Exception as e:
                self.use_ciao = False
                logger.warning(f"CIAO could not be initialized: {e}")
        return self.ciao_

    ciao = property(get_ciao)

    def __del__(self):
        if hasattr(self, "_clear") and self._clear:
            logger.info(f"Clearing {self.workdir}")
            shutil.rmtree(self.workdir, ignore_errors=True)

    def __str__(self):
        return f"OBSID={self.obsid}"

    def get_info(self):
        """
        Get observation info.

        Observation info is a combination of OBSPAR (from the obs0 file) and the sequence summary
        from https://icxc.harvard.edu/cgi-bin/mp/target.cgi.
        """
        obsid_info = self.get_obspar()
        obsid_info.update(self._get_sequence_summary())
        return obsid_info

    def file_glob(self, pattern):
        files = {}
        if self.archive_dir is not None:
            files.update(
                {
                    pth.relative_to(self.archive_dir): pth
                    for pth in self.archive_dir.glob(pattern)
                }
            )
        files.update(
            {pth.relative_to(self.workdir): pth for pth in self.workdir.glob(pattern)}
        )
        return list(files.values())

    def file_path(self, path):
        if Path(path).is_absolute():
            raise Exception(
                f"argument to Observation.file_path must be a relative path ({path})"
            )
        wp = self.workdir / path
        if self.archive_dir is not None:
            ap = self.archive_dir / path
            if not wp.exists() and ap.exists():
                return ap
        return wp

    @Cache("seq_summary", fmt="json", dir_attribute="cache_dir")
    def _get_sequence_summary(self):
        def parse_name_value(child):
            names = ["Title", "PI", "Observer", "Subject Category", "Cycle"]
            if m := re.match("(.+):(.+)", child):
                _name, _value = m.groups()
                _value = _value.strip()
                if _name in names and _value:
                    return {_name: _value}
            return {}

        obspar = self.get_obspar()
        url = "https://icxc.harvard.edu/cgi-bin/mp/target.cgi?{seq_num}"
        r = requests.get(url.format(**obspar))
        soup = bs4.BeautifulSoup(
            r.content.decode(chardet.detect(r.content)["encoding"]), features="lxml"
        )
        info = {}
        for h4 in soup.find_all("h4"):
            for child in h4.contents:
                if not isinstance(child, bs4.element.Tag):
                    info.update(parse_name_value(child))

        if "Subject Category" not in info:
            info["Subject Category"] = "NO MATCH"
        info["category_id"] = CATEGORY_ID_MAP[info["Subject Category"].lower()]

        return info

    @Cache("obspar", fmt="json", dir_attribute="cache_dir")
    def get_obspar(self):
        """
        Get the contents of the obs0 file as a dictionary.
        """
        obspar_file = self.file_glob("*obs0*")
        if not obspar_file:
            logger.debug(f"{self} No obspar file for OBSID {self.obsid}. Downloading")
            self.download(["obspar"])
        obspar_file = self.file_glob("*obs0*")
        if len(obspar_file) == 0:
            raise Exception(f"{self} No obspar file for OBSID {self.obsid}.")
        obspar_file = str(obspar_file[0])
        t = ascii.read(obspar_file)
        self._obsid_info = dict(t[["col1", "col4"]].as_array().tolist())
        self._obsid_info["instrument"] = self._obsid_info["instrume"].lower()
        self._obsid_info["obsid"] = int(self.obsid)
        self._obsid_info["target"] = self._obsid_info["object"]
        self._obsid_info["date_obs"] = self._obsid_info["date-obs"]
        self._obsid_info["dec"] = float(self._obsid_info["dec_nom"])
        self._obsid_info["ra"] = float(self._obsid_info["ra_nom"])
        self._obsid_info["roll"] = float(self._obsid_info["roll_nom"])

        import re

        m = re.match(r"(\d+).(\d+).(\d+)", self._obsid_info["ascdsver"])
        if m:
            version = [int(v) for v in m.groups()]
            self._obsid_info["version"] = (
                version[0] + version[1] / 100 + version[2] / 10000
            )
        else:
            self._obsid_info["version"] = 0.0

        return self._obsid_info

    @logging_call_decorator
    def _download_archive(self, ftypes):
        """
        Download data from chandra public archive using CIAO's download_chandra_obsid
        """
        if not self.workdir.parent.exists():
            self.workdir.parent.mkdir(exist_ok=True, parents=True)

        secondary = self.workdir / "secondary"
        if secondary.exists():
            logger.debug(f"{self} {secondary} exists, skipping download")
            return

        repro = self.file_path("repro")
        if repro.exists():
            logger.debug(f"{self} {repro} exists, skipping download")
            return

        with chdir(self.workdir.parent):
            r = subprocess.run(
                ["download_chandra_obsid", "-t", "-q"],
                stdout=subprocess.PIPE,
                env=self.ciao.env,
                check=True,
            )
            available_types = r.stdout.decode().strip().split(":")[-1].split()
            exclude = [t for t in available_types if t not in ftypes]
            process = subprocess.Popen(
                ["download_chandra_obsid", self.obsid, "--exclude", ",".join(exclude)],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                env=self.ciao.env,
            )
            output, _ = process.communicate().decode()
            output = "\n".join([f"{self} {line}" for line in output.split("\n")])
            logger.info(output)
            if process.returncode:
                raise Exception(f"{self.obsid} failed to download")

    def _download_arc5gl(self, ftypes):
        """
        Download data from chandra public archive
        """
        import Ska.arc5gl

        for ftype in ftypes:
            if ftype == "obspar":
                src, dest = "obspar", "."
            else:
                locs = self.archive_file_locations
                if ftype not in locs:
                    logger.error(
                        f"{self} {ftype=} skipped because it is not in known locations"
                    )
                    continue
                src, dest = locs[ftype]
            if src[:4] == "none":
                # if instrument is NONE, there will be no files, no point in trying
                continue
            logger.info(f"{self}   {ftype=}")
            dest_files = self.file_glob(f"{dest}/*{ftype}*")
            if dest_files:
                logger.info(f"{self}     skipping download of *{ftype}*")
                continue
            logger.info(f"{self}     {src} -> {dest}")
            (self.workdir / dest).mkdir(exist_ok=True, parents=True)
            with chdir(self.workdir / dest):
                arc5gl = Ska.arc5gl.Arc5gl()
                arc5gl.sendline(f"obsid={self.obsid}")
                arc5gl.sendline(f"get {src}")
                del arc5gl

    @logging_call_decorator
    def download(self, ftypes=None):
        """
        Download observation files using the Chandra archive or arc5gl.

        Parameters
        ----------
        ftypes: list of str
            Each must be a key in self.archive_file_locations.
        """
        if ftypes is None:
            ftypes = ["evt2", "asol", "msk", "fov", "bpix", "dtf"]
        if self._source == "archive":
            return self._download_archive(ftypes)
        elif self._source == "arc5gl":
            return self._download_arc5gl(ftypes)
        if self._source is None:
            raise Exception("No data source has been specified as fallback")
        raise Exception(f'Unknown data source: "{self._source}"')

    @logging_call_decorator
    def archive(self, *regex, destination=None):
        """
        Move observation files to an archive location.

        Parameters
        ----------
        regex: list of str
            Optional. If not given, self.archive_regex is used.
            Files matching any of the strings are arcived in a long-term location.
        destination: pathlib.Path or str.
            Optional. If not given, self.archive_dir is used.
            Long-term location where to store files.
        """
        if not regex:
            regex = self.archive_regex
        if destination is None:
            destination = self.archive_dir
        else:
            destination = Path(destination) / self.obsid
        if destination is None:
            raise Exception("archive destination was not specified")
        logger.debug(f"{self} Archiving to {destination}:")
        for pattern in regex:
            for src in self.workdir.rglob(f"**/{pattern}"):
                logger.debug(f"{self}   - {src}")
                dest = destination / src.relative_to(self.workdir)
                dest.parent.mkdir(exist_ok=True, parents=True)
                if src.is_dir():
                    dest.mkdir(exist_ok=True, parents=True)
                    for src_2 in src.glob("**/*"):
                        dest_2 = destination / src_2.relative_to(self.workdir)
                        try:
                            shutil.copy(src_2, dest_2)
                        except shutil.SameFileError:
                            pass
                else:
                    try:
                        shutil.copy(src, dest)
                    except shutil.SameFileError:
                        pass

    @logging_call_decorator
    def repro(self):
        """
        Reprocess data.
        """
        if self.file_path("images").exists():
            # Skip repro, already done
            return

        self.ciao(
            "chandra_repro",
            self.workdir,
            # outdir="",
            cleanup="yes",
            clobber="yes",
            logging_tag=str(self),
        )

    @logging_call_decorator
    def make_images(self, evt=None):
        """
        Create image. Also creates the exposure map and psfmap.
        """
        band = "wide" if self.is_hrc is True else "broad"

        outdir = self.file_path("images")
        output_files = [
            f"images/{self.obsid}_{band}_thresh.img",
            f"images/{self.obsid}.fov",
            f"images/{self.obsid}_{band}_thresh.expmap",
            f"images/{self.obsid}_{band}_thresh.psfmap",
            f"images/{self.obsid}_{band}_flux.img",
        ]
        missing_output = [fn for fn in output_files if not self.file_path(fn).exists()]
        if not missing_output:
            logger.info(f"{self}   directory {outdir} exists, skipping")
            return

        logger.info(f"Making images ({','.join([Path(m).name for m in missing_output])})")

        self.download(["evt2", "fov", "asol", "msk", "bpix", "dtf"])

        if evt is None:
            evtfiles = self.file_glob("primary/*_evt2_filtered.fits*")
            if len(evtfiles) != 1:
                raise Exception(f"Expected 1 evt file, there are {len(evtfiles)}")
            evt = evtfiles[0]

        process = subprocess.Popen(
            ["dmlist", f"{evt}[2]", "count"],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            env=self.ciao.env,
        )
        output, _ = process.communicate()
        n_events = int(output)
        if n_events <= 0:
            raise SkippedWithWarning(
                f"No events in {evt.absolute().relative_to(evt.parent.parent.absolute())}"
            )

        # are there more? is it always level1?
        fov_files = self.file_glob("primary/*_fov1.fits*")
        if len(fov_files) != 1:
            raise Exception(f"Expected 1 FOV file, there are {len(fov_files)}")
        fov_file = fov_files[0]

        outdir.mkdir(exist_ok=True)

        self.ciao(
            "fluximage",
            infile=evt,
            outroot=outdir / self.obsid,
            bands=band,
            binsize=(4 if self._rebin else 1) if self.is_hrc is True else 1,
            psfecf=0.9,
            background="none",
            logging_tag=str(self),
            clobber="yes",  # clobber will overwrite partially archived results
        )
        if self.is_acis:
            pileup_file = outdir / (self.obsid + "_pileup.img")
            pileup_max_file = outdir / (self.obsid + "_pileup_max.img")
            self.ciao(
                "pileup_map",
                infile=outdir / (self.obsid + "_" + band + "_thresh.img"),
                outfile=pileup_file,
                clobber="yes",
                logging_tag=str(self),
            )
            self.ciao(
                "dmimgfilt",
                infile=pileup_file,
                outfile=pileup_max_file,
                fun="max",
                mask="circle(0,0,3)",
                clobber="yes",
                logging_tag=str(self),
            )
            try:
                self.ciao(
                    "acis_streak_map",
                    infile=str(evt).replace("_filtered", ""),
                    fovfile=fov_file,
                    bkgroot=outdir / (self.obsid + "_acis_streaks_bkg.fits"),
                    regfile=outdir / (self.obsid + "_acis_streaks.txt"),
                    msigma="4",
                    clobber="yes",
                    logging_tag=str(self),
                )
            except Exception:
                logger.warning(f"{self}   acis_streak_map failed")

    @logging_call_decorator
    def run_wavdetect(self, edition, skip_exist=False, scales="1.4 2 4 8 16 32"):
        """
        Run wavdetect.
        """
        band = "wide" if self.is_hrc else "broad"
        root = f"{self.obsid}_{edition}"

        # this only checks that the output file sources/{self.obsid}_{edition}.src exists
        # but this function also creates:
        #    sources/{self.obsid}_{edition}.cell
        #    sources/{self.obsid}_{edition}.img
        #    sources/{self.obsid}_{edition}.nbkg
        outfile = self.file_path(f"sources/{root}.src")
        if outfile.exists() and skip_exist:
            logger.debug(f"{self} Wavdetect {edition} file already there. Skipping.")
            return
        else:
            if not skip_exist:
                logger.debug("Not skipping if wavdetect file exists")
            if not outfile.exists():
                logger.debug(
                    f"{self} Wavdetect {edition} file not there ({outfile}). Will process."
                )

        required_files = {
            "image_file": self.file_path(f"images/{self.obsid}_{band}_thresh.img"),
            "psf_file": self.file_path(f"images/{self.obsid}_{band}_thresh.psfmap"),
            "exposure_file": self.file_path(
                f"images/{self.obsid}_{band}_thresh.expmap"
            ),
        }
        if missing := [
            f"{k} ({v})" for k, v in required_files.items() if not v.exists()
        ]:
            filenames = ", ".join(missing)
            raise RuntimeError(f"Missing file(s) required for wavdetect: {filenames}")

        detdir = self.file_path("sources")
        detdir.mkdir(parents=True, exist_ok=True)

        scales = scales.split()
        # if wavdetect fails, it tries again removing the largest two scales
        for _ in range(2):
            try:
                self.ciao(
                    "wavdetect",
                    infile=required_files["image_file"],
                    expfile=required_files["exposure_file"],  # exposure map
                    psffile=required_files["psf_file"],  # PSF
                    outfile=outfile,
                    scellfile=detdir / (root + ".cell"),
                    imagefile=detdir / (root + ".img"),
                    defnbkgfile=detdir / (root + ".nbkg"),
                    scales=" ".join(scales),
                    clobber="yes",
                    logging_tag=str(self),
                )
            except Exception:
                scales = scales[:-1]
                if len(scales) < 3:
                    raise

    @logging_call_decorator
    def run_celldetect(self, snr=3):
        """
        Run celldetect.
        """
        # checks that the file sources/{self.obsid}_celldetect.src does not exist
        if self.file_path(f"sources/{self.obsid}_celldetect.src").exists():
            logger.debug(f"{self} Celldetect file already there. Skipping.")
            return

        # output is in the workdir
        out = self.workdir / "sources" / f"{self.obsid}_celldetect.src"
        band = "wide" if self.is_hrc else "broad"

        required_files = {
            "image_file": self.file_path(f"images/{self.obsid}_{band}_thresh.img"),
            "psf_file": self.file_path(f"images/{self.obsid}_{band}_thresh.psfmap"),
        }
        if missing := [
            f"{k} ({v})" for k, v in required_files.items() if not v.exists()
        ]:
            filenames = ", ".join(missing)
            raise RuntimeError(f"Missing file(s) required for celldetect: {filenames}")

        out.parent.mkdir(parents=True, exist_ok=True)

        self.ciao(
            "celldetect",
            required_files["image_file"],
            out,
            psffile=required_files["psf_file"],  # either this or set fixedcell=
            thresh=snr,
            maxlogicalwindow=2048,
            clobber="yes",
            logging_tag=str(self),
        )

    @logging_call_decorator
    def filter_events(self, radius=180, psfratio=1):  # radius in arcsec
        """
        Filter x-ray events outside a radius around the optical axis.
        """
        self.download(["evt2"])

        # I'm using a fixed pixel size of 0.5 arcsec, but this might need fixing
        pixel = 1 if self.is_hrc else 0.5
        if self._rebin and self.is_hrc:
            pixel *= 2
        try:
            evt = self.file_glob("primary/*evt2.fits*")[0]
        except Exception:
            raise SkippedWithWarning("evt2 file not found") from None

        evt2 = self.file_path(
            str(evt.relative_to(evt.parent.parent)).replace("evt2", "evt2_filtered")
        )

        if evt2.exists():
            logger.info(f"{self}   file {evt2} exists, skipping")
            return

        self.ciao("dmkeypar", evt, "RA_PNT", logging_tag=str(self))
        ra = self.ciao.pget("dmkeypar", logging_tag=str(self))
        self.ciao("dmkeypar", evt, "DEC_PNT", logging_tag=str(self))
        dec = self.ciao.pget("dmkeypar", logging_tag=str(self))

        if not ra:
            raise Exception("RA is not set")
        if not dec:
            raise Exception("DEC is not set")
        self.ciao(
            "dmcoords",
            evt,
            op="cel",
            celfmt="deg",
            ra=ra,
            dec=dec,
            logging_tag=str(self),
        )
        x = self.ciao.pget("dmcoords", "x", logging_tag=str(self))
        y = self.ciao.pget("dmcoords", "y", logging_tag=str(self))
        logger.info(f"{self} filtering circle({x},{y},{radius / pixel}).")
        self.ciao(
            "dmcopy",
            f"{evt}[(x,y)=circle({x},{y},{radius / pixel})]",
            evt2,
            logging_tag=str(self),
        )

        return evt2

    @logging_call_decorator
    def filter_sources(self, radius=180, psfratio=1):  # radius in arcsec
        """
        Filter detected sources outside a radius around the optical axis.
        """
        self.download(["evt2"])

        # I'm using a fixed pixel size of 0.5 arcsec, but this might need fixing
        pixel = 1 if self.is_hrc else 0.5
        if self._rebin and self.is_hrc:
            pixel *= 2
        try:
            evt = self.file_glob("primary/*evt2.fits*")[0]
        except Exception:
            raise Exception("evt2 file not found   ") from None

        src = self.file_path(f"sources/{self.obsid}_baseline.src")
        if not src.exists():
            raise Exception("src file not found   ")
        out = self.file_path(f"sources/{self.obsid}_filtered.src")

        self.ciao("dmkeypar", evt, "RA_PNT", logging_tag=str(self))
        ra = self.ciao.pget("dmkeypar", logging_tag=str(self))
        self.ciao("dmkeypar", evt, "DEC_PNT", logging_tag=str(self))
        dec = self.ciao.pget("dmkeypar", logging_tag=str(self))

        self.ciao("dmcoords", evt, op="cel", celfmt="deg", ra=ra, dec=dec)
        x = self.ciao.pget("dmcoords", "x", logging_tag=str(self))
        y = self.ciao.pget("dmcoords", "y", logging_tag=str(self))
        filters = [f"psfratio=:{psfratio}", f"(x,y)=circle({x},{y},{radius / pixel})"]
        filters = ",".join(filters)
        self.ciao("dmcopy", f"{src}[{filters}]", out, logging_tag=str(self))

        return out

    @logging_call_decorator
    def calculate_centroids(self):
        """
        Re-compute centroids.
        """
        src_hdus = fits.open(self.file_path(f"sources/{self.obsid}_baseline.src"))
        band = "wide" if self.is_hrc else "broad"
        img = self.file_path(f"images/{self.obsid}_{band}_flux.img")

        if not img.exists():
            raise Exception(f"Image file not found {img}")
        src = table.Table(src_hdus[1].data)
        result = []
        for row in src:
            x, y = row[["X", "Y"]]
            r = row["R"].max()
            self.ciao(
                "dmstat",
                f"{img}[sky=circle({x},{y},{r})]",
                centroid="yes",
                # clip='yes',
                logging_tag=str(self),
            )
            x, y = np.array(
                self.ciao.pget("dmstat", "out_cntrd_phys", logging_tag=str(self)).split(
                    ","
                )
            ).astype(float)
            result.append([x, y])
        return np.array(result, dtype=[("x", float), ("y", float)])

    @logging_call_decorator
    def process(self):
        """
        Main routine to process a single obsid with "standard" steps.

        This function skips:
            - observations with obs_mode other than "POINTING".
            - ACIS observations with readmode other than "TIMED" and dtycycle different than 0.
        """

        obsid_info = self.get_info()

        ok = (
            int(obsid_info["obsid"]) not in _multi_obi_obsids
            # and obsid_info['category_id'] not in [110]
            and obsid_info["obs_mode"] == "POINTING"
            # and obsid_info['grating'] == 'NONE'
            and (
                obsid_info["instrume"] == "HRC"
                or (
                    obsid_info["instrume"] == "ACIS"
                    and obsid_info["readmode"] == "TIMED"
                    and int(obsid_info["dtycycle"]) == 0
                )
            )
        )
        if not ok:
            raise Skipped("does not fulfill observation requirements")

        # Analysis
        self.filter_events()
        self.make_images()
        self.run_wavdetect("baseline", skip_exist=True)
        self.run_celldetect()

    @Cache("sources", fmt="table", dir_attribute="cache_dir")
    @logging_call_decorator
    def get_sources(self, version="celldetect"):
        """
        Returns a table of sources formatted for the astromon_xray_source SQL table
        """
        from chandra_aca.transform import radec_to_yagzag
        from Quaternion import Quat

        obspar = self.get_obspar()
        q = Quat(equatorial=(obspar["ra_pnt"], obspar["dec_pnt"], obspar["roll_pnt"]))

        hdu_list = fits.open(self.file_path(f"sources/{self.obsid}_{version}.src"))
        sources = table.Table(hdu_list[1].data)

        if len(sources) == 0:
            raise SkippedWithWarning("No x-ray sources found")

        sources["pileup"] = self._pileup_value(sources)
        sources["acis_streak"] = self._on_acis_streak(sources)

        if "SNR" not in sources.colnames:
            logger.debug(f"{self} adding masked column for SNR.")
            sources["SNR"] = table.MaskedColumn(
                length=len(sources), mask=np.ones(len(sources))
            )

        columns = [
            c
            for c in zip(
                ["RA", "DEC", "COMPONENT", "NET_COUNTS", "SNR", "PSFRATIO"],
                ["ra", "dec", "id", "net_counts", "snr", "psfratio"],
                strict=True,
            )
            if c[0] in sources.colnames
        ]
        sources.rename_columns(*list(zip(*columns, strict=True)))

        sources["obsid"] = int(self.obsid)
        sources["y_angle"], sources["z_angle"] = radec_to_yagzag(
            sources["ra"], sources["dec"], q
        )
        sources["r_angle"] = np.sqrt(sources["y_angle"] ** 2 + sources["z_angle"] ** 2)

        # calculate the distance to the closest source
        src1 = sources.as_array()[None]
        src2 = sources.as_array()[:, None]
        i = np.diag([np.inf] * len(sources))
        distance = (
            np.sqrt(
                (src1["y_angle"] - src2["y_angle"]) ** 2
                + (src1["z_angle"] - src2["z_angle"]) ** 2
            )
            + i
        )
        sources["near_neighbor_dist"] = np.min(distance, axis=0)
        sources["caldb_version"] = self.get_calalign()["caldb_version"]

        cols = [
            "obsid",
            "id",
            "ra",
            "dec",
            "net_counts",
            "y_angle",
            "z_angle",
            "r_angle",
            "snr",
            "near_neighbor_dist",
            "psfratio",
            "pileup",
            "acis_streak",
            "caldb_version",
        ]
        return sources[cols]

    def _pileup_value(self, src):
        pileup_file = self.file_path(f"images/{self.obsid}_pileup_smeared.img")
        if not pileup_file.exists():
            return np.zeros(len(src))
        hdus = fits.open(pileup_file)
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message="'datfix' made the change 'Set DATEREF",
                category=FITSFixedWarning,
            )
            wcs = WCS(hdus[0].header)
        loc = SkyCoord(src["RA"] * u.deg, src["DEC"] * u.deg)
        pix = np.round(wcs.world_to_pixel(loc)).astype(int)
        return hdus[0].data[(pix[1], pix[0])]

    def _on_acis_streak(self, src):
        acis_streaks_file = self.file_path(f"images/{self.obsid}_acis_streaks.txt")
        result = np.zeros(len(src), dtype=bool)
        if acis_streaks_file.exists():
            reg = regions.Regions.read(acis_streaks_file)
            pos = regions.PixCoord(x=src["X"], y=src["Y"])
            for pol in reg.regions:
                result += pol.contains(pos)
        return result

    @property
    def archive_file_locations(self):
        instrument = self.get_obspar()["instrument"]
        return {
            "obspar": ("obspar", "."),
            "evt2": (f"{instrument}2{{evt2}}", "primary"),
            "evt1": (f"{instrument}1{{evt1}}", "secondary"),
            "fov": (f"{instrument}1{{fov}}", "primary"),
            "msk": (f"{instrument}1{{msk}}", "secondary"),
            "mtl": (f"{instrument}1{{mtl}}", "secondary"),
            "bpix": (f"{instrument}1[*bpix*]", "primary"),
            "flt": (f"{instrument}1[*flt*]", "secondary"),
            "stat": (f"{instrument}1[*stat*]", "secondary"),
            "asol": ("asp1[*asol*]", "secondary"),
            "acal": ("asp1[*acal*]", "secondary"),
            "dtf": (f"{instrument}1[*dtf1*]", "primary"),
            # 'pbk': f'{instrument}0[*pbk0*]',
            # 'bias': f'{instrument}0[*bias*]',
        }

    @Cache("calalign", fmt="json", dir_attribute="cache_dir")
    def get_calalign(self):
        self.download(["acal"])
        cal_file = self.file_glob("secondary/*acal*fits*")
        if cal_file:
            hdus = fits.open(cal_file[0])
            calalign = {
                k: hdus[1].data[k]
                for k in ["aca_align", "aca_misalign", "fts_misalign"]
            }
            calalign = {k: v.tolist() for k, v in calalign.items()}
            calalign["obsid"] = int(self.obsid)
            calalign["caldb_version"] = hdus[1].header["CALDBVER"]
            return calalign


def get_parser():
    """
    Get the argument parser to process a set of observations.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "obsid", help="OBSID(s) to process or a file with a list of OBSIDs", nargs="+"
    )
    parser.add_argument(
        "--workdir",
        help="Working directory. Default is to create a tmp directory",
        default=".",
        type=Path,
    )
    parser.add_argument(
        "--ciao-prefix", help="CIAO prefix", default="/soft/ciao", type=Path
    )
    parser.add_argument(
        "--data-source",
        help="Where to get the data",
        default="archive",
        choices=["archive", "arc5gl"],
    )
    parser.add_argument(
        "--log-level",
        help="logging level",
        choices=["debug", "info", "warning"],
        default="info",
    )
    return parser


def _process(o, **kwargs):
    observation = Observation(o, **kwargs)
    observation.process()


def main():
    """
    Main routine to process a set of observations.
    """
    from functools import partial
    from multiprocessing import Pool

    import pyyaks.logger

    args = get_parser().parse_args()

    pyyaks.logger.get_logger(name="astromon", level=args.log_level.upper())

    logger.info(f"workdir: {args.workdir}")
    if args.workdir is None:
        tmp = tempfile.TemporaryDirectory()
        args.workdir = tmp.name
    Path(args.workdir).mkdir(exist_ok=True, parents=True)

    obsids = []
    for obsid in args.obsid:
        if os.path.exists(obsid) and os.path.isfile(obsid):
            with open(obsid) as fh:
                obsids += [line.strip() for line in fh.readlines()]
        else:
            obsids += obsid.split(",")

    process = partial(
        _process,
        workdir=args.workdir,
        ciao_prefix=args.ciao_prefix,
        source=args.data_source,
    )
    # '9' because the archive limits number of concurrent connections
    # from same IP address (well, it did when it was 'ftp').
    with Pool(processes=9) as pool:
        pool.starmap(process, [(o,) for o in obsids])


if __name__ == "__main__":
    main()
