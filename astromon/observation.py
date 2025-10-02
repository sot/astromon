#!/usr/bin/env python


import argparse
import collections
import functools

# import sys
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
import Ska.arc5gl
from astropy import table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii, fits
from astropy.wcs import WCS, FITSFixedWarning
from chandra_aca.transform import radec_to_yagzag
from Quaternion import Quat
from ska_helpers.logging import basic_logger

from astromon.stored_result import Storage, stored_result
from astromon.task import TASKS, ReturnCode, ReturnValue, dependencies, run_tasks, task
from astromon.utils import Ciao, chdir, logging_call_decorator

__all__ = ["Observation"]


logger = basic_logger(__name__, level="WARNING")

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


ARCHIVE_DIR = Path(os.environ["SKA"]) / "data" / "astromon" / "archive"


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
        # checking if the workdir exists before creating is faster than passing the
        # exist_ok=True argument to mkdir, and there are thousands of observations
        # so this is a significant speedup
        if not self.workdir.exists():
            self.workdir.mkdir(parents=True)
        self.archive_dir = (
            (Path(archive_dir).expanduser() / subdir / self.obsid)
            if archive_dir
            else ARCHIVE_DIR / subdir / self.obsid
        )
        self.storage = Storage(
            workdir=self.workdir,
            archive_dir=self.archive_dir,
        )

        self.cache_dir = self.workdir / "cache"

        if archive_regex is None:
            self.archive_regex = [
                "*.par.gz",
                "cache",
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

    is_multi_obi = property(
        lambda self: int(self.get_obspar()["obsid"]) in _multi_obi_obsids
    )

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

    @stored_result("seq_summary", fmt="json", subdir="cache")
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

    @stored_result("obspar", fmt="json", subdir="cache")
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

        types = {
            "r": float,
            "s": str,
            "i": int,
        }
        self._obsid_info = {row["col1"]: types[row["col2"]](row["col4"]) for row in t}
        self._obsid_info["instrument"] = self._obsid_info["instrume"].lower()
        self._obsid_info["obsid"] = int(self.obsid)
        self._obsid_info["target"] = self._obsid_info["object"]
        self._obsid_info["date_obs"] = self._obsid_info["date-obs"]
        self._obsid_info["dec"] = float(self._obsid_info["dec_nom"])
        self._obsid_info["ra"] = float(self._obsid_info["ra_nom"])
        self._obsid_info["roll"] = float(self._obsid_info["roll_nom"])

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

    def _download_arc5gl(self, ftypes, revision=None, force=False):
        """
        Download data from chandra public archive
        """

        if ftypes != ["obspar"] and revision is None:
            # the first thing we do is to get the obspar file, to know the revision number.
            # self.get_obspar calls self.download with ftypes=['obspar'], and that is why we just
            # checked that ftypes is not ['obspar'], to avoid infinite recursion.
            # The obspar is cached, so this happens only once.
            obspar = self.get_obspar()
            revision = obspar["revision"]

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
            if dest_files and not force:
                logger.info(f"{self}     skipping download of *{ftype}*")
                continue
            logger.info(f"{self}     {src} -> {dest}")
            dest = self.workdir / dest
            if not dest.exists():
                dest.mkdir(exist_ok=True, parents=True)
            with chdir(dest):
                arc5gl = Ska.arc5gl.Arc5gl()
                arc5gl.sendline(f"obsid={self.obsid}")
                if revision is not None:
                    arc5gl.sendline(f"version={revision}")
                arc5gl.sendline(f"get {src}")
                del arc5gl

    @logging_call_decorator
    def download(self, ftypes=None, revision=None, force=False):
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
            return self._download_arc5gl(ftypes, revision=revision, force=force)
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
        regex = regex if regex else self.archive_regex
        destination = (
            Path(destination) / self.obsid if destination else self.archive_dir
        )
        if destination is None:
            raise Exception("archive destination was not specified")
        logger.debug(f"{self} Archiving to {destination}:")
        for pattern in regex:
            for src in self.workdir.rglob(f"**/{pattern}"):
                logger.debug(f"{self}   - {src}")
                dest = destination / src.relative_to(self.workdir)
                if not dest.parent.exists():
                    dest.parent.mkdir(exist_ok=True, parents=True)
                if src.is_dir():
                    if not dest.exists():
                        dest.mkdir(exist_ok=True, parents=True)
                    for src_2 in src.glob("**/*"):
                        if src_2.is_dir():
                            # Only files are copied. Parent directories are created automatically.
                            # empty directories are not copied
                            continue
                        dest_2 = destination / src_2.relative_to(self.workdir)
                        if not dest_2.parent.exists():
                            dest_2.parent.mkdir(exist_ok=True, parents=True)
                        try:
                            shutil.copy(src_2, dest_2)
                        except shutil.SameFileError:
                            # links to the same file are not copied
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
    def make_images(self):
        """
        Create image. Also creates the exposure map and psfmap.
        """
        return make_images.run(self, run_dependencies=True)

    @logging_call_decorator
    def run_wavdetect(self):
        """
        Run wavdetect.
        """
        return wavdetect.run(self, run_dependencies=True)

    @logging_call_decorator
    def run_celldetect(self):
        """
        Run celldetect.
        """
        return celldetect.run(self, run_dependencies=True)

    @logging_call_decorator
    def filter_events(self):
        """
        Filter x-ray events outside a radius around the optical axis.
        """
        return filter_events.run(self, run_dependencies=True)

    @logging_call_decorator
    def filter_sources(self):
        """
        Filter detected sources outside a radius around the optical axis.
        """
        return filter_sources.run(self, run_dependencies=True)

    @dependencies(
        optional_files={
            "sources": "sources/{obsid}_baseline.src",
            "images": "images/{obsid}_{band}_flux.img",
        },
        variables={
            "band": lambda obs: "wide" if obs.is_hrc else "broad",
        },
    )
    @logging_call_decorator
    def calculate_centroids(self):
        """
        Re-compute centroids.
        """
        dtype = [("x", float), ("y", float)]
        if not (
            wavdetect.get_result(self).return_code == ReturnCode.OK
            and make_images.get_result(self).return_code == ReturnCode.OK
        ):
            return np.array([], dtype=dtype)

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
        return np.array(result, dtype=dtype)

    @property
    def is_selected(self):
        """
        Check if the observation fulfills the requirements for astromon processing.

        This function skips:
            - observations with obs_mode other than "POINTING".
            - ACIS observations with readmode other than "TIMED" and dtycycle different than 0.

        Returns
        -------
        bool
            True if the observation is suitable for astromon processing, False otherwise.
        """
        obsid_info = self.get_info()

        return (
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

    @logging_call_decorator
    def process(self):
        """
        Main routine to process a single obsid with "standard" steps.
        """

        if not self.is_selected:
            return ReturnValue(
                return_code=ReturnCode.SKIP,
                msg="observation does not fulfill requirements",
            )

        band = "wide" if self.is_hrc else "broad"
        rv = run_tasks(
            self,
            requested_files=[
                f"images/{self.obsid}_{band}_flux.img",
                f"sources/{self.obsid}_celldetect.src",
                f"sources/{self.obsid}_baseline.src",
            ],
        )
        if any(v.return_code.value >= ReturnCode.ERROR.value for v in rv.values()):
            return [
                v for v in rv.values() if v.return_code.value >= ReturnCode.ERROR.value
            ][0]
        elif any(v.return_code.value >= ReturnCode.SKIP.value for v in rv.values()):
            return [
                v for v in rv.values() if v.return_code.value >= ReturnCode.SKIP.value
            ][0]
        else:
            return ReturnValue(return_code=ReturnCode.OK)

    @stored_result("sources", fmt="table", subdir="cache")
    @dependencies(
        optional_files={
            "sources": "sources/{obsid}_{version}.src",
            "pileup": "images/{obsid}_pileup_smeared.img",
            "acis_streaks": "images/{obsid}_acis_streaks.txt",
        },
    )
    @logging_call_decorator
    def get_sources(self, *, version="celldetect"):
        """
        Returns a table of sources formatted for the astromon_xray_source SQL table
        """
        # default dtype to be used when returning an empty list
        dtype = [
            ("obsid", ">i8"),
            ("id", ">i4"),
            ("ra", ">f8"),
            ("dec", ">f8"),
            ("net_counts", ">f4"),
            ("y_angle", ">f8"),
            ("z_angle", ">f8"),
            ("r_angle", ">f8"),
            ("snr", ">f4"),
            ("near_neighbor_dist", ">f8"),
            ("psfratio", ">f4"),
            ("pileup", ">f4"),
            ("acis_streak", "?"),
            ("caldb_version", "<U6"),
        ]

        if TASKS.tasks[version].get_result(self).return_code != ReturnCode.OK:
            return table.Table(dtype=dtype)

        obspar = self.get_obspar()
        q = Quat(equatorial=(obspar["ra_pnt"], obspar["dec_pnt"], obspar["roll_pnt"]))

        hdu_list = fits.open(self.file_path(f"sources/{self.obsid}_{version}.src"))
        sources = table.Table(hdu_list[1].data)

        if len(sources) == 0:
            return table.Table(dtype=dtype)

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

    @stored_result("calalign", fmt="json", subdir="cache")
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

    def get_asol_files(self):
        return [self.file_path(f) for f in self._get_asol_files_cached()]

    @stored_result("asol", subdir="cache")
    def _get_asol_files_cached(self):
        # this function is cached because it downloads the events file if needed
        # and that can be slow.
        # asol file is specified in events file header
        # trying filtered event file first, because they are usually cached
        evt = self.file_path(f"primary/{self.obsid}_evt2_filtered.fits.gz")
        if not evt.exists():
            # if the filtered event file is not found, fall back to the unfiltered version
            self.download(["evt2"])
            evt_files = self.file_glob("primary/*_evt2.fits*")
            if len(evt_files) == 0:
                raise Exception("Expected one evt file, there are none")
            evt = evt_files[0]

        self.download(["asol"])

        hdu = fits.open(evt)
        if "ASOLFILE" in hdu[1].header:
            filenames = hdu[1].header["ASOLFILE"].split(",")
            filenames = [
                (fi + ".gz" if fi.endswith(".fits") else fi) for fi in filenames
            ]
            asol = [self.file_path(f"secondary/{filename}") for filename in filenames]

            # the file does not exist, check if the file revision is different
            # from the current revision
            if not asol[0].exists() and (
                mr := re.search(r"N(?P<revision>\d\d\d)_asol1.fits", asol[0].name)
            ):
                revision = int(mr.group("revision"))
                cur_revision = int(self.get_obspar()["revision"])
                if revision != cur_revision:
                    self.download(["asol"], revision=revision, force=True)
        else:
            # asol file is not given in the event file,
            # see if there is something and hope for the best
            asol = self.file_glob("secondary/*asol*fits*")

        if any(not f.exists() for f in asol):
            missing = ", ".join([str(f) for f in asol if not f.exists()])
            raise Exception(f"Aspect solution files not found: {missing}")

        # the result needs to be a relative path so it can be cached for later use
        # because the files might currently be in a temporary directory
        asol = [f.relative_to(f.parent.parent) for f in asol]
        return asol

    def dmcoords(self, name, **kwargs):
        """
        Call dmcoords with the given arguments.

        Examples
        --------

        Get the off-axis angle, given (x, y) "sky" coordinates:

            obs.dmcoords("theta", option="sky", x=4069.94266994267, y=4076.716625716626)

        Get the off-axis angle, given (ra, dec) in celestial ("cel") coordinates:

            obs.dmcoords("theta", option="cel", celfmt="deg", ra=20.46451186, dec=-28.34952557)

        Parameters
        ----------
        name: str
            Name of the output parameter.
        kwargs: dict
            Arguments to pass to dmcoords.
        """

        # trying filtered event file first, because they are usually cached
        evt = self.file_path(f"primary/{self.obsid}_evt2_filtered.fits.gz")
        if not evt.exists():
            # if the filtered event file is not found, fall back to the unfiltered version
            self.download(["evt2"])
            evt_files = self.file_glob("primary/*_evt2.fits*")
            if len(evt_files) == 0:
                raise Exception("Expected one evt file, there are none")
            evt = evt_files[0]

        # asol file is specified in events file header
        asol = self.get_asol_files()[0]

        args = [evt, f"asolfile={asol}"]
        # if I do not unlearn, the following two calls in succession will hang
        # obs.dmcoords("theta", option="sky", x=4069.94266994267, y=4076.716625716626)
        # obs.dmcoords("theta", option="cel", celfmt="deg", ra=20.46451186, dec=-28.34952557)
        self.ciao("punlearn", "dmcoords", logging_tag=str(self))
        self.ciao("dmcoords", *args, logging_tag=str(self), **kwargs)
        value = self.ciao.pget("dmcoords", name, logging_tag=str(self))
        return np.array(value).astype(float)


@task(
    name="make_images",
    inputs={
        "events": "primary/*_evt2.fits*",
        "filtered_events": "primary/{obsid}_evt2_filtered.fits.gz",
        "fov": "primary/*_fov1.fits*",
        "asol_file": "secondary/*asol*fits*",
        "badpixfile": "primary/*_bpix1.fits*",
        "maskfile": "secondary/*_msk1.fits*",
    },
    optional_inputs={
        "dtffile": "primary/*_dtf1.fits*",
    },
    outputs={
        "image_file": "images/{obsid}_{band}_thresh.img",
        "fov": "images/{obsid}.fov",
        "exposure_file": "images/{obsid}_{band}_thresh.expmap",
        "psf_file": "images/{obsid}_{band}_thresh.psfmap",
        "flux": "images/{obsid}_{band}_flux.img",
        "pileup": "images/{obsid}_pileup.img",
        "pileup_max": "images/{obsid}_pileup_max.img",
        "acis_streaks_bkg": "images/{obsid}_acis_streaks_bkg.fits",
        "acis_streaks": "images/{obsid}_acis_streaks.txt",
    },
    download=(["evt2", "fov", "asol", "msk", "bpix", "dtf"]),
    variables={
        "band": lambda obs: "wide" if obs.is_hrc else "broad",
    },
)
def make_images(obs, inputs, outputs):
    """
    Create image. Also creates the exposure map and psfmap.
    """
    bin_size = (4 if obs._rebin else 1) if obs.is_hrc is True else 1

    evt_filt = inputs["filtered_events"]
    fov_file = inputs["fov"][0]

    logging_tag = str(obs)
    band = "wide" if obs.is_hrc is True else "broad"
    ciao = obs.ciao
    obsid = obs.obsid

    # this makes sure the asol file is downloaded, even in rare cases where it is not downloaded
    # when "asol" is downloaded
    asol = ",".join([str(a) for a in obs.get_asol_files()])

    process = subprocess.Popen(
        ["dmlist", f"{evt_filt}[2]", "count"],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        env=ciao.env,
    )
    output, _ = process.communicate()
    n_events = int(output)
    if n_events <= 0:
        msg = f"{obsid}   No events in {evt_filt.relative_to(evt_filt.parent.parent)}"
        return ReturnCode.SKIP, msg

    kwargs = {}
    if "dtffile" in inputs:
        kwargs["dtffile"] = inputs["dtffile"][0]

    ciao(
        "fluximage",
        infile=evt_filt,
        outroot=outputs["image_file"].parent / f"{obsid}",
        asolfile=asol,
        badpixfile=inputs["badpixfile"][0],
        maskfile=inputs["maskfile"][0],
        bands=band,
        binsize=bin_size,
        psfecf=0.9,
        background="none",
        logging_tag=logging_tag,
        clobber="yes",  # clobber will overwrite partially archived results
        **kwargs,
    )
    if obs.is_acis:
        pileup_file = outputs["pileup"]
        pileup_max_file = outputs["pileup_max"]
        ciao(
            "pileup_map",
            infile=outputs["image_file"],
            outfile=pileup_file,
            clobber="yes",
            logging_tag=logging_tag,
        )
        ciao(
            "dmimgfilt",
            infile=pileup_file,
            outfile=pileup_max_file,
            fun="max",
            mask="circle(0,0,3)",
            clobber="yes",
            logging_tag=logging_tag,
        )
        try:
            ciao(
                "acis_streak_map",
                infile=inputs["events"][0],
                fovfile=fov_file,
                bkgroot=outputs["acis_streaks_bkg"],
                regfile=outputs["acis_streaks"],
                msigma="4",
                clobber="yes",
                logging_tag=logging_tag,
            )
        except Exception:
            logger.warning(f"{obsid}   acis_streak_map failed")


@task(
    name="wavdetect",
    inputs={
        "image_file": "images/{obsid}_{band}_thresh.img",
        "psf_file": "images/{obsid}_{band}_thresh.psfmap",
        "exposure_file": "images/{obsid}_{band}_thresh.expmap",
    },
    outputs={
        "src": "sources/{obsid}_baseline.src",
        "cell": "sources/{obsid}_baseline.cell",
        "img": "sources/{obsid}_baseline.img",
        "nbkg": "sources/{obsid}_baseline.nbkg",
    },
    variables={
        "band": lambda obs: "wide" if obs.is_hrc else "broad",
    },
)
def wavdetect(obs, inputs, outputs):
    """
    Run wavdetect.
    """
    # possible parameter:
    scales = ["1.4", "2", "4", "8", "16", "32"]

    # if wavdetect fails, it tries again removing the largest two scales
    for _ in range(2):
        try:
            obs.ciao(
                "wavdetect",
                infile=inputs["image_file"],
                expfile=inputs["exposure_file"],  # exposure map
                psffile=inputs["psf_file"],  # PSF
                outfile=outputs["src"],
                scellfile=outputs["cell"],
                imagefile=outputs["img"],
                defnbkgfile=outputs["nbkg"],
                scales=" ".join(scales),
                clobber="yes",
                logging_tag=str(obs),
            )
        except Exception:
            scales = scales[:-1]
            if len(scales) < 3:
                raise
    return ReturnCode.OK


@task(
    name="celldetect",
    inputs={
        "image_file": "images/{obsid}_{band}_thresh.img",
        "psf_file": "images/{obsid}_{band}_thresh.psfmap",
    },
    outputs={
        "src": "sources/{obsid}_celldetect.src",
    },
    variables={
        "band": lambda obs: "wide" if obs.is_hrc else "broad",
    },
)
def celldetect(obs, inputs, outputs):
    """
    Run celldetect.
    """
    # possible parameter:
    snr = 3

    obs.ciao(
        "celldetect",
        inputs["image_file"],
        outputs["src"],
        psffile=inputs["psf_file"],  # either this or set fixedcell=
        thresh=snr,
        maxlogicalwindow=2048,
        clobber="yes",
        logging_tag=str(obs),
    )

    return ReturnCode.OK


@task(
    name="filter_events",
    inputs={
        "events": "primary/*_evt2.fits*",
        "fov": "primary/*_fov1.fits*",
    },
    outputs={
        "events": "primary/{obsid}_evt2_filtered.fits.gz",
    },
    download=(["evt2", "fov"]),
)
def filter_events(obs, inputs, outputs):
    """
    Filter x-ray events outside a radius around the optical axis.
    """
    # possible parameter:
    radius = 180  # radius in arcsec

    # I'm using a fixed pixel size of 0.5 arcsec, but this might need fixing
    pixel = 1 if obs.is_hrc else 0.5

    evt = inputs["events"][0]
    evt2 = outputs["events"]

    obs.ciao("dmkeypar", evt, "RA_PNT", logging_tag=str(obs))
    ra = obs.ciao.pget("dmkeypar", logging_tag=str(obs))
    obs.ciao("dmkeypar", evt, "DEC_PNT", logging_tag=str(obs))
    dec = obs.ciao.pget("dmkeypar", logging_tag=str(obs))

    if not ra:
        raise Exception("RA is not set")
    if not dec:
        raise Exception("DEC is not set")
    obs.ciao("punlearn", "dmcoords", logging_tag=str(obs))
    obs.ciao(
        "dmcoords",
        evt,
        op="cel",
        celfmt="deg",
        ra=ra,
        dec=dec,
        logging_tag=str(obs),
    )
    x = obs.ciao.pget("dmcoords", "x", logging_tag=str(obs))
    y = obs.ciao.pget("dmcoords", "y", logging_tag=str(obs))
    logger.info(f"{obs} filtering circle({x},{y},{radius / pixel}).")
    obs.ciao(
        "dmcopy",
        f"{evt}[(x,y)=circle({x},{y},{radius / pixel})]",
        evt2,
        logging_tag=str(obs),
        clobber="yes",
    )


@task(
    name="filter_sources",
    inputs={
        "events": "primary/{obsid}_evt2_filtered.fits.gz",
        "src": "sources/{obsid}_baseline.src",
    },
    outputs={
        "src": "sources/{obsid}_filtered.src",
    },
    download=(["evt2"]),
)
def filter_sources(obs, inputs, outputs):
    """
    Filter detected sources outside a radius around the optical axis.
    """
    # possible parameters:
    radius = 180  # radius in arcsec
    psfratio = 1

    # I don't think the following should be a parameter
    # I'm using a fixed pixel size of 0.5 arcsec, but this might need fixing
    pixel = 1 if obs.is_hrc else 0.5

    evt = inputs["events"]

    obs.ciao("dmkeypar", evt, "RA_PNT", logging_tag=str(obs))
    ra = obs.ciao.pget("dmkeypar", logging_tag=str(obs))
    obs.ciao("dmkeypar", evt, "DEC_PNT", logging_tag=str(obs))
    dec = obs.ciao.pget("dmkeypar", logging_tag=str(obs))

    obs.ciao("punlearn", "dmcoords", logging_tag=str(obs))
    obs.ciao("dmcoords", evt, op="cel", celfmt="deg", ra=ra, dec=dec)
    x = obs.ciao.pget("dmcoords", "x", logging_tag=str(obs))
    y = obs.ciao.pget("dmcoords", "y", logging_tag=str(obs))
    filters = [f"psfratio=:{psfratio}", f"(x,y)=circle({x},{y},{radius / pixel})"]
    filters = ",".join(filters)

    obs.ciao(
        "dmcopy", f"{inputs['src']}[{filters}]", outputs["src"], logging_tag=str(obs)
    )


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
    from functools import partial  # noqa: PLC0415
    from multiprocessing import Pool  # noqa: PLC0415

    args = get_parser().parse_args()

    logger.setLevel(args.log_level.upper())

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
