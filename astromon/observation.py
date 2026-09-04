#!/usr/bin/env python


import argparse
import collections
import functools
import gzip

# import sys
import os
import re
import shutil
import subprocess
import tempfile
import warnings
from datetime import datetime
from pathlib import Path

import numpy as np
import regions
import Ska.arc5gl
from astropy import table
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import ascii, fits
from astropy.wcs import WCS, FITSFixedWarning
from chandra_aca.transform import radec_to_yagzag, yagzag_to_radec
from mica.archive import cda
from Quaternion import Quat
from ska_helpers.logging import basic_logger

from astromon import source_detection, utils
from astromon.stored_result import Storage, stored_result
from astromon.task import TASKS, ReturnCode, ReturnValue, dependencies, run_tasks, task
from astromon.utils import Ciao, chdir, logging_call_decorator

__all__ = ["Observation", "ObsidNotPubliclyAvailable", "ObsparUnavailable"]


class ObsparUnavailable(RuntimeError):
    """Raised when no obspar can be obtained for an obsid.

    The obspar has exactly two sources: the local mica archive, or an arc5gl
    download. The public archive is not one of them -- CDA does not publish the
    observation parameter file among download_chandra_obsid's filetypes.
    """


class ObsidNotPubliclyAvailable(RuntimeError):
    """Raised when an obsid is not available via the public CDA archive.

    Treated as a "skip" rather than a pipeline failure — the obsid is
    intentionally unavailable (embargoed, permanently proprietary, or
    simply not yet released), and retrying will not help.
    """


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

# CDA ocat `category` strings differ from ID_CATEGORY_MAP labels (e.g. "STARS AND WD"
# vs "Normal Stars and WD"), so map them directly to numeric IDs.
CDA_CATEGORY_ID_MAP = collections.defaultdict(
    lambda: 200,
    {
        "STARS AND WD": 10,
        "WD BINARIES AND CVS": 20,
        "BH AND NS BINARIES": 30,
        "NORMAL GALAXIES": 40,
        "ACTIVE GALAXIES AND QUASARS": 50,
        "EXTRAGALACTIC DIFFUSE EMISSION AND SURVEYS": 60,
        "GALACTIC DIFFUSE EMISSION AND SURVEYS": 70,
        "SOLAR SYSTEM AND MISC": 100,
        "SN, SNR AND ISOLATED NS": 110,
        "CLUSTERS OF GALAXIES": 120,
    },
)
"""Mapping from CDA ocat category strings to numeric category IDs."""


ARCHIVE_DIR = Path(os.environ["SKA"]) / "data" / "astromon" / "xray_observations"


def _flag_brightest_source(snr_arr: np.ndarray) -> np.ndarray:
    """Return a bool array with True for the single highest-SNR source.

    NaN values are ignored. If all values are NaN, returns all False.

    Parameters
    ----------
    snr_arr : np.ndarray
        1-D array of SNR values (any float dtype, may contain NaN).

    Returns
    -------
    np.ndarray
        Boolean array of the same length; exactly one element is True
        (the argmax over finite values), or all False if no finite values exist.

    Examples
    --------
    >>> import numpy as np
    >>> _flag_brightest_source(np.array([3.0, 10.0, 5.0]))
    array([False,  True, False])
    >>> _flag_brightest_source(np.array([float('nan'), float('nan')]))
    array([False, False])
    """
    snr = np.asarray(snr_arr, dtype=float)
    result = np.zeros(len(snr), dtype=bool)
    if np.any(np.isfinite(snr)):
        result[int(np.argmax(np.where(np.isfinite(snr), snr, -np.inf)))] = True
    return result


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
                "*_acis_streaks.fits",
                "*.src",
                "*_psf_size_*.fits",
            ]
        else:
            self.archive_regex = archive_regex
        self._source = source
        logger.info(f"{self} starting. Context: {self.workdir}")
        self._rebin = True
        self.ciao_ = None
        # ciao is initialized only if needed, but we make this check at the very beginning
        # so the user gets a warning at the top if calling CIAO will likely fail
        if self.use_ciao and (msg := ciao_fails(ciao_prefix, workdir)):
            self.use_ciao = False
            logger.warning(msg)

    is_multi_obi = property(lambda self: int(self.obsid) in _multi_obi_obsids)

    is_hrc = property(lambda self: self.get_evt2_info()["instrument"] == "hrc")

    is_acis = property(lambda self: self.get_evt2_info()["instrument"] == "acis")

    def cleanup_downloads(self):
        """Remove intermediate files that are only needed during processing.

        Call this after all detection versions and cross-matches are complete.
        Files that the pipeline never reads (level-1 events, VV report, preview
        JPEGs, etc.) are excluded at download time by ``_download_archive`` and
        therefore never hit the disk.

        This method handles the processing intermediates:

        Deleted
        -------
        primary/
            Raw level-2 event file (``*_evt2.fits``, excluding the filtered copy).
            The filtered copy (created by filter_events) is sufficient to re-run
            any source detection step; re-downloading is all-or-nothing anyway.
        sources/
            wavdetect auxiliary outputs (``*.cell``, ``*.img``, ``*.nbkg``).  The
            actual source lists (``*.src``) and PSF-size tables (``*.fits``) are
            kept so gaussian_detect can be re-run without re-running wavdetect.
        param/
            CIAO parameter files written during the run.
        results/
            Per-obsid FITS copies of the DB tables (redundant with astromon.h5).
        images/
            fluximage outputs (flux map, PSF map, exposure map), pileup maps,
            and acis_streak images.  Not needed for aggregate analysis or
            cross-matching.

            Regenerating them means re-downloading: ``make_images`` declares
            ``primary/*_evt2.fits*`` -- the raw level-2 events -- among its
            inputs, and this method deletes exactly that, keeping only the
            filtered copy.  So a cleaned workdir cannot re-image in place.
            ``clear_invalid_make_images.py`` exists because that bit in
            production: ~437 obsids were marked for re-imaging and failed with
            "Missing input files for task make_images: events".

        Kept
        ----
        primary/
            Filtered evt2, asol, and bpix — sufficient to rerun source
            detection on a single outlier without re-downloading from CDA.
        sources/
            Source lists (``*.src``) and PSF-size tables (``*.fits``).
        secondary/
            Mask file.
        cache/
            Task-state JSON files so reruns skip completed steps.
        """
        deleted_bytes = 0

        def _remove(path):
            nonlocal deleted_bytes
            if path.is_file() and not path.is_symlink():
                deleted_bytes += path.stat().st_size
                path.unlink()

        # Raw evt2 — keep the filtered copy, delete the full original.
        # Both share "*_evt2.fits*" but the filtered file has "_filtered" in its name.
        primary = self.workdir / "primary"
        for path in primary.glob("*_evt2.fits*"):
            if "_filtered" not in path.name:
                _remove(path)

        # wavdetect cell/image/background auxiliaries — not needed for any re-run;
        # gaussian_detect works from the filtered evt2 + *.src source list directly.
        sources_dir = self.workdir / "sources"
        for pattern in ("*.cell", "*.img", "*.nbkg"):
            for path in sources_dir.glob(pattern):
                _remove(path)

        # CIAO parameter files, DB-redundant results, and flux/pileup/streak images.
        # images/ is entirely regenerable from the filtered evt2 + asol.
        for dirname in ("param", "results", "images"):
            dirpath = self.workdir / dirname
            if dirpath.exists():
                for path in dirpath.rglob("*"):
                    _remove(path)
                shutil.rmtree(dirpath, ignore_errors=True)

        # images/ is gone, so make_images is no longer done. Its stored result
        # lives in cache/, which this method deliberately keeps, so without this it
        # would go on reporting success while its outputs were deleted -- and
        # pileup, acis_streak and grating_arm would come back as measured zeros on
        # any later recompute. This does not make the images regenerable in place;
        # the raw evt2 make_images needs is deleted above. What it buys is that a
        # recompute fails loudly on the missing input rather than quietly writing
        # zeros. follow_dependents=False because the detection results stay valid:
        # keeping sources/*.src is the point of this cleanup.
        make_images.invalidate_result(self, follow_dependents=False)
        logger.info(f"{self} cleanup_downloads removed {deleted_bytes / 1e6:.1f} MB")

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
                # Clear any stale .par files from previous interrupted runs.
                # A partial or interrupted fluximage (or other CIAO tool) can
                # leave a malformed param file that causes subsequent runs to
                # fail with "punlearn() Could not open parameter file".
                param_dir = self.workdir / "param"
                if param_dir.exists():
                    shutil.rmtree(param_dir)
                self.ciao_ = Ciao(
                    prefix=self.ciao_prefix,
                    workdir=param_dir,
                    logger="astromon",
                )
            except Exception as e:
                self.use_ciao = False
                logger.warning(f"CIAO could not be initialized: {e}")
        if self.ciao_ is None:
            # Fail here with an actionable message instead of letting the caller use None as
            # if it were a Ciao instance, which surfaces as a bare TypeError/AttributeError
            # far from the actual cause.
            raise RuntimeError(
                f"{self} requires CIAO, but CIAO is not available "
                f"(ciao_prefix={self.ciao_prefix}). See earlier log messages for why."
            )
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

        Observation info is a combination of pointing/configuration metadata read
        directly from the evt2 file's own header (:any:`get_evt2_info`) and
        proposal-level metadata -- title, PI, category, and the ocat's detector label
        -- from the Chandra ocat (:any:`_get_sequence_summary`).  The obspar is not
        used: it is a separately-versioned archive product that is not guaranteed to
        be in sync with the evt2 file actually being processed (see
        :any:`get_evt2_info`).
        """
        obsid_info = self.get_evt2_info()
        seq_summary = self._get_sequence_summary()
        obsid_info.update(seq_summary)
        obsid_info["detector"] = seq_summary.get("instr", "")
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
        """Get observation metadata (title, PI, category) from the Chandra ocat.

        Tries the local HDF5 ocat first (fast, available at CXC), then falls
        back to the public CDA web service (works anywhere).  If both fail,
        returns a default with category_id=200 (Unknown).

        Returns a dict with keys Title, PI, Observer, Subject Category, Cycle,
        category_id, instr (the ocat detector label, e.g. "ACIS-S"), mode (readout
        mode, "TE" or "CC"), and d_cyc ("Y"/"N", duty cycling in use).  instr/mode/
        d_cyc are also used by :any:`is_selected` to gate observations, in place of
        obspar's obs_mode/instrume/readmode/dtycycle: non-science secondary/duplicate
        segments come back from the ocat with a blank or "NONE" instr, and mode/d_cyc
        line up with readmode/dtycycle (TE<->TIMED, CC<->CONTINUOUS) to better than
        99.99% across the archive.
        """
        obsid_int = int(self.obsid)
        ocat_row = None
        for fetch in (
            lambda: cda.get_ocat_local(obsid_int),
            lambda: cda.get_ocat_web(obsid_int),
        ):
            try:
                ocat_row = fetch()
                break
            except Exception as exc:
                logger.debug(f"{self} ocat fetch failed: {exc}")

        if ocat_row is None:
            logger.warning(f"{self} could not fetch ocat data; using Unknown category")
            return {
                "Title": "",
                "PI": "",
                "Observer": "",
                "Subject Category": "NO MATCH",
                "Cycle": "",
                "category_id": 200,
                "instr": "",
                "mode": "",
                "d_cyc": "",
            }

        cda_category = str(ocat_row.get("category", "")).upper()
        category_id = CDA_CATEGORY_ID_MAP[cda_category]
        return {
            "Title": str(ocat_row.get("prop_title", "")),
            "PI": str(ocat_row.get("pi_name", "")),
            "Observer": str(ocat_row.get("observer", "")),
            "Subject Category": cda_category.title(),
            "Cycle": str(ocat_row.get("obs_cycle", "")),
            "category_id": category_id,
            "instr": str(ocat_row.get("instr", "")),
            "mode": str(ocat_row.get("mode", "")),
            "d_cyc": str(ocat_row.get("d_cyc", "")),
        }

    def _get_mica_obspar(self) -> dict | None:
        """Obspar from the local mica archive, or None if it has no entry.

        mica.archive.obspar reads $SKA/data/mica/archive/obspar and returns None
        when the obsid is absent.  Any other failure (no $SKA, archive not synced)
        is also treated as a miss so the caller can fall back to downloading.
        """
        from mica.archive import obspar as mica_obspar  # noqa: PLC0415

        try:
            return mica_obspar.get_obspar(int(self.obsid))
        except Exception as exc:
            logger.debug(f"{self} mica obspar lookup failed: {exc}")
            return None

    @stored_result("obspar", fmt="json", subdir="cache")
    def get_obspar(self):
        """
        Get the contents of the obs0 file as a dictionary.

        Sources are tried in order: an obspar already in the working directory,
        then the local mica obspar archive, then an arc5gl download.

        mica reads $SKA/data/mica/archive/obspar, so the common case needs neither
        the network nor arc5gl.  That ordering matters beyond convenience: every
        task parameter substitution resolves through here, so downloading before
        consulting mica made the whole task framework unreachable off the HEAD
        network -- including the dependency and cache-invalidation tests.
        """

        def _find_obspar():
            # arc5gl downloads the obspar to workdir root; download_chandra_obsid
            # places it in secondary/.  Check both.
            return self.file_glob("*obs0*") or self.file_glob("secondary/*obs0*")

        obspar_file = _find_obspar()
        mica_info = None
        if not obspar_file:
            mica_info = self._get_mica_obspar()
            if mica_info is None:
                logger.debug(
                    f"{self} No obspar file or mica entry for OBSID {self.obsid}."
                    " Downloading"
                )
                self.download(["obspar"])
                obspar_file = _find_obspar()

        if mica_info is not None:
            logger.debug(f"{self} obspar from mica.archive.obspar")
            # mica's parser already normalises 'date-obs' to 'date_obs'.
            self._obsid_info = dict(mica_info)
        elif len(obspar_file) > 0:
            obspar_file = str(obspar_file[0])
            t = ascii.read(obspar_file)
            types = {"r": float, "s": str, "i": int}
            self._obsid_info = {
                row["col1"]: types[row["col2"]](row["col4"]) for row in t
            }
            # par files use 'date-obs' (hyphen); normalise to 'date_obs'
            self._obsid_info["date_obs"] = self._obsid_info.pop("date-obs", "")
        else:
            raise ObsparUnavailable(
                f"{self} no obspar available for OBSID {self.obsid}: not in the"
                " working directory and not in the local mica archive"
                " ($SKA/data/mica/archive/obspar). Its only other source is an"
                ' arc5gl download, which needs source="arc5gl" and the CXC'
                " network."
            )

        self._obsid_info["instrument"] = self._obsid_info["instrume"].lower()
        self._obsid_info["obsid"] = int(self.obsid)
        self._obsid_info["target"] = self._obsid_info["object"]
        self._obsid_info["dec"] = float(self._obsid_info["dec_nom"])
        self._obsid_info["ra"] = float(self._obsid_info["ra_nom"])
        self._obsid_info["roll"] = float(self._obsid_info["roll_nom"])

        m = re.match(r"(\d+).(\d+).(\d+)", str(self._obsid_info.get("ascdsver", "")))
        if m:
            version = [int(v) for v in m.groups()]
            self._obsid_info["version"] = (
                version[0] + version[1] / 100 + version[2] / 10000
            )
        else:
            self._obsid_info["version"] = 0.0

        return self._obsid_info

    @stored_result("evt2_info", fmt="json", subdir="cache")
    @dependencies(download=["evt2"])
    def get_evt2_info(self):
        """
        Get pointing and observation metadata directly from the evt2 file's own header.

        This is the source of truth for anything that must be self-consistent with the
        events actually being processed -- most importantly the pointing
        (``ra_pnt``/``dec_pnt``/``roll_pnt``), which is what ``filter_events``'s 180"
        circular cut is centered on.  ``get_obspar`` is a separately-versioned archive
        product and can genuinely disagree with the evt2 file's own header for the same
        obsid: confirmed empirically for a real observation, where obspar's cached
        RA_PNT differed from the evt2 header's RA_PNT (and from the actual measured
        aspect trajectory) by 38.6 arcsec, because the two products are revisioned
        independently by the archive and are not reprocessed together.

        Returns
        -------
        dict
            ra_pnt/dec_pnt/roll_pnt (and ra/dec/roll aliases), ra_targ/dec_targ,
            revision, obs_mode, readmode, dtycycle, grating, instrume, instrument
            (lowercase), detnam, target, sim_x/sim_y/sim_z, datamode, ascdsver, obsid,
            and a derived ``version`` float -- all read straight from the evt2 file's
            primary event-list header, the same file ``filter_events`` reads via
            ``dmkeypar``.
        """
        event_files = self.file_glob("primary/*_evt2.fits*")
        if len(event_files) == 0:
            raise Exception(f"{self} No evt2 file for OBSID {self.obsid}.")

        header = fits.getheader(str(event_files[0]), 1)

        info = {
            "ra_pnt": float(header["RA_PNT"]),
            "dec_pnt": float(header["DEC_PNT"]),
            "roll_pnt": float(header["ROLL_PNT"]),
            "ra": float(header["RA_PNT"]),
            "dec": float(header["DEC_PNT"]),
            "roll": float(header["ROLL_PNT"]),
            "ra_targ": float(header["RA_TARG"]),
            "dec_targ": float(header["DEC_TARG"]),
            "obs_mode": str(header.get("OBS_MODE", "")),
            "readmode": str(header.get("READMODE", "")),
            "dtycycle": int(header["DTYCYCLE"]) if "DTYCYCLE" in header else None,
            "grating": str(header.get("GRATING", "")),
            "instrume": str(header.get("INSTRUME", "")),
            "instrument": str(header.get("INSTRUME", "")).lower(),
            "detnam": str(header.get("DETNAM", "")),
            "target": str(header.get("OBJECT", "")),
            "sim_x": float(header["SIM_X"]) if "SIM_X" in header else None,
            "sim_y": float(header["SIM_Y"]) if "SIM_Y" in header else None,
            "sim_z": float(header["SIM_Z"]) if "SIM_Z" in header else None,
            "datamode": str(header.get("DATAMODE", "")),
            "ascdsver": str(header.get("ASCDSVER", "")),
            "date_obs": str(header.get("DATE-OBS", "")),
            "obsid": int(self.obsid),
        }

        m = re.search(r"N(?P<revision>\d{3})_evt2", str(event_files[0]))
        info["revision"] = int(m.group("revision")) if m else 0

        vm = re.match(r"(\d+).(\d+).(\d+)", info["ascdsver"])
        if vm:
            v = [int(x) for x in vm.groups()]
            info["version"] = v[0] + v[1] / 100 + v[2] / 10000
        else:
            info["version"] = 0.0

        return info

    @stored_result("ra_dec_wcs", fmt="pickle", subdir="cache")
    @dependencies(download=["evt2"])
    def get_ra_dec_wcs(self):
        """
        Get WCS to convert between sky coordinates (x, y) and (ra, dec).
        """
        event_files = self.file_glob("primary/*_evt2.fits*")
        if len(event_files) == 0:
            raise Exception(f"{self} No evt2 file for OBSID {self.obsid}.")

        wcs = utils.get_wcs_from_fits_header(event_files[0], hdu=1)

        return wcs

    def get_off_axis_angle(self, ra=None, dec=None, x=None, y=None):
        """
        Get APPROXIMATE off-axis angle in arcmin for given (ra, dec) or sky (x, y) coordinates.

        This is not the most accurate calculation of the off-axis angle. It uses the coordinates in
        the sky coordinate system, which is a tangent plane at the nominal pointing. It is still
        within 1 arcsec of the true off-axis angle for all cases we care about.

        Parameters
        ----------
        ra: np.array
            Right ascension in degrees.
        dec: np.array
            Declination in degrees.
        x: np.array
            Sky X coordinate in pixels.
        y: np.array
            Sky Y coordinate in pixels.

        Returns
        -------
        np.array
            Off-axis angle in arcmin. Same shape as input arrays.
        """
        errors = []
        if (ra is None) != (dec is None):
            errors.append("Both ra and dec must be provided together.")
        if (x is None) != (y is None):
            errors.append("Both x and y must be provided together.")
        radec = ra is not None and dec is not None
        skyxy = x is not None and y is not None
        if radec and skyxy:
            errors.append("Either ra/dec or x/y must be provided, not both.")
        if not (radec or skyxy):
            errors.append("Either ra/dec or x/y must be provided.")
        if errors:
            raise Exception(" ".join(errors))

        wcs = self.get_ra_dec_wcs()
        if skyxy:
            xy = np.array([x, y]).T
        else:
            xy = np.array(wcs.world_to_pixel_values(ra, dec)).T
        uv = (wcs.wcs.cdelt * (xy - wcs.wcs.crpix)).T
        return np.sqrt(np.sum(uv**2, axis=0)) * 60

    # Filetypes needed by the pipeline.  download_chandra_obsid accepts a
    # comma-separated list so we request only what we use, rather than
    # downloading the full set (evt1/evt1a, VV reports, spectra, JPEGs, …)
    # and deleting the excess afterward.
    _ARCHIVE_FILETYPES = "asol,bpix,dtf,evt2,fov,msk"

    # Requested filetypes the public archive cannot supply at all, as opposed to
    # the ones _download_archive covers by always fetching _ARCHIVE_FILETYPES.
    _ARCHIVE_UNAVAILABLE = ("obspar",)

    def _archive_download_marker(self) -> Path:
        """File written once an archive download has fully succeeded.

        The presence of ``secondary/`` cannot serve as this signal:
        download_chandra_obsid creates it early, and neither the timeout nor the
        error path removes it, so an interrupted download left behind a directory
        that every later attempt read as "already done" -- turning one transient
        failure into a workdir that needed manual cleanup before it would ever
        download again. It lives inside ``secondary/`` so that removing that
        directory resets the state deliberately.
        """
        return self.workdir / "secondary" / ".download_complete"

    def _fix_archive_bpix(self) -> None:
        """Gunzip the bad pixel file from secondary/ into primary/ as a real file.

        CIAO fluximage cannot find the bad pixel file when it is a symlink (even
        if the symlink target is valid).  ``download_chandra_obsid`` places the
        bpix file in ``secondary/`` and the old code created a symlink in
        ``primary/``; this replaces that with a real decompressed copy.

        Also removes any stale ``.fits.gz`` symlink left by older runs so that
        the ``primary/*_bpix1.fits*`` glob in ``make_images`` returns exactly
        one candidate.

        Idempotent: safe to call on every run regardless of whether the download
        was cached.
        """
        secondary = self.workdir / "secondary"
        primary = self.workdir / "primary"
        for bpix_file in secondary.glob("*bpix*fits*"):
            is_gz = bpix_file.suffix == ".gz"
            if is_gz:
                dest = primary / bpix_file.stem  # strip .gz → plain .fits
                stale = primary / bpix_file.name  # old .gz symlink, if any
            else:
                dest = primary / bpix_file.name
                stale = primary / (bpix_file.name + ".gz")
            if stale.exists() or stale.is_symlink():
                stale.unlink()
            if not dest.exists():
                # Only the .gz case needs decompressing; opening a plain FITS
                # file with gzip.open raises BadGzipFile.
                opener = gzip.open if is_gz else open
                with opener(bpix_file, "rb") as f_in, open(dest, "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)

    @logging_call_decorator
    def _download_archive(self, ftypes):  # noqa: PLR0912
        """
        Download pipeline-required files for this obsid from the Chandra public archive.

        Uses CIAO's ``download_chandra_obsid`` to fetch only the filetypes the
        pipeline actually reads (``_ARCHIVE_FILETYPES``).  ``ftypes`` is accepted
        for API compatibility with ``_download_arc5gl`` and is otherwise ignored --
        except for the filetypes CDA does not publish at all, listed in
        ``_ARCHIVE_UNAVAILABLE``, which raise :class:`ObsparUnavailable`. Ignoring
        those silently meant a request for the obspar fetched a full evt2 and five
        other filetypes and then failed claiming the download produced nothing.
        A marker written after the download and its layout fix-ups have all
        succeeded is the "already downloaded" signal -- see
        :meth:`_archive_download_marker` for why the ``secondary/`` directory
        cannot be. Subsequent calls are then no-ops.
        """
        unavailable = sorted(set(ftypes or ()) & set(self._ARCHIVE_UNAVAILABLE))
        if unavailable:
            raise ObsparUnavailable(
                f"obsid {self.obsid}: {', '.join(unavailable)} cannot be downloaded"
                " from the public archive -- CDA does not publish it among"
                " download_chandra_obsid's filetypes. It has to come from the local"
                " mica archive ($SKA/data/mica/archive/obspar) or from arc5gl with"
                ' source="arc5gl".'
            )

        secondary = self.workdir / "secondary"
        marker = self._archive_download_marker()
        if marker.exists():
            logger.debug(f"{self} {marker} present, skipping download")
            self._fix_archive_bpix()
            return

        # Check public release date before attempting download, so we get a
        # clear failure message rather than a silent "not found on archive site".
        # Use get_ocat_local (mica's on-disk cache) to avoid a network round-trip;
        # fall back to get_ocat_web only if the local cache misses.
        try:
            try:
                info = cda.get_ocat_local(int(self.obsid))
            except Exception:
                info = cda.get_ocat_web(self.obsid, summary=True)
            pa_raw = info.get("public_avail", "")
            # get_ocat_local returns '0' for embargoed/no-release-date observations.
            # get_ocat_web returns np.ma.masked for the same case (np.ma.is_masked
            # needed because bool(masked) raises TypeError in newer numpy).
            pa_str = str(pa_raw).strip() if not np.ma.is_masked(pa_raw) else ""
            if not pa_str or pa_str == "0":
                # No release date set — embargoed or permanently proprietary.
                raise ObsidNotPubliclyAvailable(
                    f"obsid {self.obsid} has no public release date "
                    f"(public_avail={pa_raw!r}); not available via public archive"
                )
            release_date = datetime.strptime(pa_str, "%Y-%m-%d %H:%M:%S")
            if release_date > datetime.now():
                raise ObsidNotPubliclyAvailable(
                    f"obsid {self.obsid} not yet publicly available "
                    f"(release date: {pa_str})"
                )
        except ObsidNotPubliclyAvailable:
            raise
        except Exception as exc:
            logger.warning(f"{self} could not check public release date: {exc}")

        if not self.workdir.parent.exists():
            self.workdir.parent.mkdir(exist_ok=True, parents=True)

        # Maximum time to wait for download_chandra_obsid.  A large observation
        # at a slow network rate might take ~2 hours; beyond that the connection
        # is almost certainly hung (broken TCP with no data flowing).
        _DOWNLOAD_TIMEOUT_SEC = 7200

        with chdir(self.workdir.parent):
            process = subprocess.Popen(
                ["download_chandra_obsid", self.obsid, self._ARCHIVE_FILETYPES],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                env=self.ciao.env,
            )
            try:
                stdout_bytes, _ = process.communicate(timeout=_DOWNLOAD_TIMEOUT_SEC)
            except subprocess.TimeoutExpired:
                process.kill()
                process.communicate()
                raise RuntimeError(
                    f"obsid {self.obsid} download_chandra_obsid timed out "
                    f"after {_DOWNLOAD_TIMEOUT_SEC}s (hung connection)"
                ) from None
            stdout_text = stdout_bytes.decode(errors="replace")
            for line in stdout_text.splitlines():
                logger.info(f"{self} {line}")
            # download_chandra_obsid exits 0 even when the obsid is not found on
            # the archive site — detect that case from its output.  Raise a plain
            # RuntimeError (retryable failure) rather than ObsidNotPubliclyAvailable
            # because the same message appears on transient network failures; only the
            # ocat check above has enough context to declare an obsid permanently
            # unavailable.
            if "not found on the archive site" in stdout_text:
                raise RuntimeError(
                    f"obsid {self.obsid} not found on the CDA archive site "
                    f"(download_chandra_obsid: {stdout_text.strip()!r})"
                )
            if process.returncode:
                raise Exception(
                    f"{self.obsid} download_chandra_obsid failed "
                    f"(exit code {process.returncode})"
                )

        # download_chandra_obsid follows the official CDA layout, which differs
        # from the arc5gl layout that the rest of the pipeline expects:
        #   asol: CDA puts it in primary/, pipeline expects secondary/
        #   bpix: CDA puts it in secondary/, pipeline expects primary/
        # Create symlinks so downstream code finds files in the expected places.
        primary = self.workdir / "primary"
        for asol_file in primary.glob("*asol*fits*"):
            link = secondary / asol_file.name
            if not link.exists():
                link.symlink_to(asol_file)
        self._fix_archive_bpix()

        # Last, so that anything raising above leaves the download retryable.
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.write_text("download_chandra_obsid completed\n")

    def _download_arc5gl(self, ftypes, revision=None, force=False):
        """
        Download data from chandra public archive.

        With ``revision=None`` (the default), each ftype is fetched at its own current
        archive default, independently.  This used to default to obspar's own revision
        number and pass that as ``version=`` for every ftype, on the assumption that one
        revision number meant the same processing epoch across product types.  It
        doesn't: evt2, asol, and obspar each carry their own independently-incremented
        revision counter for the same obsid (confirmed empirically -- N003/N001/N002 for
        the same real obsid, all "current" at the same time), so pinning one product's
        counter onto a request for another was never actually achieving self-consistency,
        and at worst could silently fetch a superseded file for a product type whose
        revision history happens to include one at that same number.  Pass ``revision``
        explicitly only when there's a specific, known reason to request an exact one.
        """

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
            - obsids >= 38000, the reserved block for engineering/calibration-replay
              sequences, which are never real science pointings.
            - obsids with no usable Chandra ocat instrument value.  Non-science
              secondary/duplicate segments come back from the ocat with a blank or
              "NONE" instr even when a real ocat row exists for the obsid.
            - ACIS observations not in TE (timed exposure) readout mode, or with
              duty cycling in use.  Continuous Clocking (CC) mode collapses one
              spatial dimension for high-time-resolution timing, so 2-D source
              detection cannot produce a valid position from it.

        All of this comes from the Chandra ocat rather than the obspar: the obspar is
        a separately-versioned archive product that is not guaranteed to be in sync
        with the obsid actually being processed (see :any:`get_evt2_info`), while the
        ocat fields used here (instr, mode, d_cyc) line up with the obspar's
        obs_mode/instrume/readmode/dtycycle to better than 99.99% across the archive.

        Returns
        -------
        bool
            True if the observation is suitable for astromon processing, False otherwise.
        """
        obsid_int = int(self.obsid)
        if obsid_int >= 38000 or obsid_int in _multi_obi_obsids:
            return False

        seq = self._get_sequence_summary()
        instr = str(seq.get("instr", "")).strip().upper()
        if instr in ("", "NONE"):
            return False

        if instr.startswith("ACIS"):
            return seq.get("mode", "") == "TE" and seq.get("d_cyc", "") == "N"
        return True

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

    @dependencies(
        optional_files={
            "sources": "sources/{obsid}_{version}.src",
            "pileup": "images/{obsid}_pileup_max.img",
            "acis_streaks": "images/{obsid}_acis_streaks.fits",
            "grating_arm_mask": "images/{obsid}_grating_arm_mask.fits",
            "psf_size": "sources/{obsid}_psf_size_{version}.fits",
        },
    )
    @logging_call_decorator
    def _get_sources(self, *, version="celldetect"):
        """
        Read the output from CIAO source detection, add some columns, and return an astropy Table.
        """

        hdu_list = fits.open(self.file_path(f"sources/{self.obsid}_{version}.src"))
        sources = table.Table(hdu_list[1].data)
        # this sets the native endianness. Some packages do not support big endian.
        sources = table.Table(sources.as_array())

        if len(sources) == 0:
            return table.Table()

        if "y_angle" not in sources.colnames or "z_angle" not in sources.colnames:
            evt2_info = self.get_evt2_info()
            q = Quat(
                equatorial=(
                    evt2_info["ra_pnt"],
                    evt2_info["dec_pnt"],
                    evt2_info["roll_pnt"],
                )
            )
            sources["y_angle"], sources["z_angle"] = radec_to_yagzag(
                sources["RA"], sources["DEC"], q
            )
        sources["r_angle"] = np.sqrt(sources["y_angle"] ** 2 + sources["z_angle"] ** 2)

        # filter_events cuts input events to a 180" radius (see that task), but
        # celldetect (and gaussian_detect, which seeds from it) fits cells/Gaussians
        # on the full rectangular flux image built from those events, so a detection's
        # position is not itself constrained to stay inside that circle.  Drop the ones
        # that land beyond it.
        sources = sources[sources["r_angle"] < 180]

        if "near_neighbor_dist" not in sources.colnames:
            # this uses y_angle and z_angle
            sources["near_neighbor_dist"] = utils.get_near_neighbor_dist(
                sources, sources
            )

        # checking upper and lower case because it gets renamed below
        if "SNR" not in sources.colnames and "snr" not in sources.colnames:
            logger.debug(f"{self} adding masked column for SNR.")
            sources["SNR"] = table.MaskedColumn(
                length=len(sources), mask=np.ones(len(sources))
            )

        if "ecf_radius" not in sources.colnames:
            ecf = table.Table.read(
                self.file_path(f"sources/{self.obsid}_psf_size_{version}.fits")
            )
            pixel_size = 0.4920 if self.is_acis else 0.13175
            sources["ecf_radius"] = ecf["R"] * pixel_size

        sources["obsid"] = int(self.obsid)
        sources["caldb_version"] = self.get_calalign()["caldb_version"]
        sources["pileup"] = self._pileup_value(sources)
        sources["acis_streak"] = self._on_acis_streak(sources)
        sources["grating_arm"] = self._on_grating_arm(sources)

        # Flag the single brightest source (highest SNR) in the observation.
        # This is typically the calibration target; the flag lets downstream
        # cross-matching include it even if it is also flagged acis_streak=True
        # (e.g. when the grating dispersed spectrum triggers a streak detection).
        snr_col = "SNR" if "SNR" in sources.colnames else "snr"
        sources["brightest"] = _flag_brightest_source(sources[snr_col])

        # Computed before the COMPONENT -> id rename below.
        sources["peak_offset"] = self._peak_offset(sources)

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

        return sources

    @stored_result("sources", fmt="table", subdir="cache")
    def get_sources(self, *, version="celldetect", astromon_format=True):
        """
        Returns a table of sources formatted for the astromon_xray_source SQL table.

        If astromon_format is False, returns all the columns in the source detection output.
        """

        # default dtype to be used when returning an empty list
        dtype = np.dtype(
            [
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
                ("concentration_ratio", ">f4"),
                ("pileup", ">f4"),
                ("acis_streak", "?"),
                ("grating_arm", "?"),
                ("brightest", "?"),
                ("caldb_version", "<U20"),
                ("detect_method", "<U24"),
                ("peak_offset", ">f4"),
            ]
        )

        if (
            TASKS.tasks[version].get_result(self) is not None
            and TASKS.tasks[version].get_result(self).return_code != ReturnCode.OK
        ):
            return table.Table(dtype=dtype)

        sources = self._get_sources(version=version)

        if len(sources) == 0:
            return table.Table(dtype=dtype)

        sources["detect_method"] = version

        # Fill columns present in dtype but not produced by this detection version
        # (e.g. concentration_ratio is only computed by gaussian_detect and
        # peak_gaussian_detect, not celldetect). Imported here because db imports
        # this module at load time.
        from astromon import db  # noqa: PLC0415

        db.conform_to_dtype(sources, "astromon_xray_src")

        return (
            table.Table(sources[dtype.names], dtype=dtype)
            if astromon_format
            else sources
        )

    def _pileup_value(self, src):
        pileup_file = self.file_path(f"images/{self.obsid}_pileup_max.img")
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
        acis_streaks_file = self.file_path(f"images/{self.obsid}_acis_streaks.fits")
        result = np.zeros(len(src), dtype=bool)
        if acis_streaks_file.exists():
            acis_streaks = table.Table.read(acis_streaks_file)
            polygons = []
            for row in acis_streaks:
                vertices = regions.PixCoord(x=row["X"], y=row["Y"])
                polygons.append(regions.PolygonPixelRegion(vertices=vertices))
            pos = regions.PixCoord(x=src["X"], y=src["Y"])
            for pol in polygons:
                result |= pol.contains(pos)
        return result

    def _peak_offset(self, src):
        """Angular distance from each source to the peak-seeded Gaussian fit, arcsec.

        The peak-seeded fit re-seeds from the local emission maximum instead of
        the detection centroid, so how far it moves is a measure of how stable
        the centroid is: a large offset means extended, confused, or faint
        emission. Storing that distance is what lets peak_gaussian_detect stop
        being persisted as its own detect_method.

        Sources are paired on ``COMPONENT``, which every ``.src`` file in the
        chain carries with the same values -- it is the celldetect component
        number, threaded through both Gaussian fits.

        Returns NaN where the offset is not measurable: no peak fit on disk, no
        COMPONENT column, or a component the peak fit has no row for (the fit
        can fail on a source celldetect found). NaN rather than 0 because 0 is a
        meaningful value here -- perfect agreement -- and must not be
        indistinguishable from "never measured".
        """
        result = np.full(len(src), np.nan)
        peak_file = self.file_path(f"sources/{self.obsid}_peak_gaussian_detect.src")
        if "COMPONENT" not in src.colnames or not peak_file.exists():
            return result

        peak = table.Table.read(peak_file)
        if len(peak) == 0:
            return result
        peak_pos = {
            int(component): (float(ra), float(dec))
            for component, ra, dec in zip(
                peak["COMPONENT"], peak["RA"], peak["DEC"], strict=True
            )
        }

        rows = [
            (i, peak_pos[int(component)])
            for i, component in enumerate(src["COMPONENT"])
            if int(component) in peak_pos
        ]
        if not rows:
            return result

        idx = [i for i, _ in rows]
        here = SkyCoord(
            np.asarray(src["RA"], dtype=float)[idx] * u.deg,
            np.asarray(src["DEC"], dtype=float)[idx] * u.deg,
        )
        there = SkyCoord(
            [p[0] for _, p in rows] * u.deg, [p[1] for _, p in rows] * u.deg
        )
        result[idx] = here.separation(there).arcsec
        return result

    def _on_grating_arm(self, src):
        """Return bool array: True if source falls inside a grating arm corridor.

        Sources within the zero-order circle (TG_PART=0) are not flagged; only
        sources inside a dispersed-arm rotbox (TG_PART=1 for HEG, 2 for MEG,
        3 for LETG) but outside the zero-order circle are marked.
        """
        mask_file = self.file_path(f"images/{self.obsid}_grating_arm_mask.fits")
        result = np.zeros(len(src), dtype=bool)
        if not mask_file.exists():
            return result

        mask_regions = table.Table.read(mask_file)
        pos = regions.PixCoord(x=src["X"], y=src["Y"])

        # Identify the zero-order circle row (TG_PART=0) to build the exclusion zone.
        zo_rows = mask_regions[mask_regions["TG_PART"] == 0]
        on_zero_order = np.zeros(len(src), dtype=bool)
        for row in zo_rows:
            zo_center = regions.PixCoord(x=row["X"], y=row["Y"])
            zo_circle = regions.CirclePixelRegion(
                center=zo_center, radius=float(row["R"][0])
            )
            on_zero_order |= zo_circle.contains(pos)

        # Flag sources inside any arm rotbox that are not in the zero-order circle.
        arm_rows = mask_regions[mask_regions["TG_PART"] != 0]
        for row in arm_rows:
            arm_rect = regions.RectanglePixelRegion(
                center=regions.PixCoord(x=row["X"], y=row["Y"]),
                width=2.0 * float(row["R"][0]),
                height=2.0 * float(row["R"][1]),
                angle=float(row["ROTANG"]) * u.deg,
            )
            result |= arm_rect.contains(pos) & ~on_zero_order

        return result

    @property
    def archive_file_locations(self):
        # Needed to decide how to download evt2 itself, so it must come from something
        # cheap and pre-download -- the ocat's instr label (e.g. "ACIS-S"), not
        # get_evt2_info (which requires evt2 to already exist -- that would be circular).
        instrument = str(self._get_sequence_summary().get("instr", "")).split("-")[0].lower()
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
        # acal file not available (e.g. when data was downloaded via
        # download_chandra_obsid which does not include CXC-internal products).
        # Fall back to reading CALDBVER from the L2 event file, which records
        # the CALDB version used when CXC processed the observation.
        caldb_version = "0.0"
        evt_files = self.file_glob("primary/*_evt2.fits*")
        if evt_files:
            try:
                with fits.open(evt_files[0]) as hdus:
                    caldb_version = hdus[1].header.get("CALDBVER", "0.0")
            except Exception as exc:
                logger.debug(f"{self} could not read CALDBVER from event file: {exc}")
        logger.debug(
            f"{self} no acal file; using CALDBVER={caldb_version!r} from event file"
        )
        return {"obsid": int(self.obsid), "caldb_version": caldb_version}

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

            # If the expected asol file (per evt2's own ASOLFILE header) is missing,
            # download it at the exact revision that filename names.
            if not asol[0].exists() and (
                mr := re.search(r"N(?P<revision>\d\d\d)_asol1.fits", asol[0].name)
            ):
                revision = int(mr.group("revision"))
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
        "acis_streaks": "images/{obsid}_acis_streaks.fits",
        "grating_arm_mask": "images/{obsid}_grating_arm_mask.fits",
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
            # the acis_streak_map command only creates ascii files
            # (even though there is one fits example in the CIAO docs)
            # and the regions package does not read that text file, so I call acis_streak_map
            # and then convert the file to fits using dmmakereg
            acis_streaks_file_ascii = outputs["acis_streaks"].parent / outputs[
                "acis_streaks"
            ].name.replace(".fits", ".reg")
            ciao(
                "acis_streak_map",
                infile=inputs["events"][0],
                fovfile=fov_file,
                bkgroot=outputs["acis_streaks_bkg"],
                regfile=acis_streaks_file_ascii,
                msigma="4",
                ssigma="3",
                clobber="yes",
                logging_tag=logging_tag,
            )
            # to make matters worse, dmmakereg chokes on txt files with no regions
            with open(acis_streaks_file_ascii, "r") as fh:
                n_regions = len(
                    [line for line in fh.readlines() if re.search("Polygon", line)]
                )
            if n_regions > 0:
                obs.ciao.logger.info(
                    f"{logging_tag} acis_streak_map found {n_regions} streaks"
                )
                ciao(
                    "dmmakereg",
                    f"region({acis_streaks_file_ascii})",
                    outputs["acis_streaks"],
                    "kernel=fits",
                    "clobber=yes",
                )
            else:
                # we want feedback to go to the same place as the output from CIAO
                obs.ciao.logger.info(f"{logging_tag} acis_streak_map found no streaks")
        except Exception:
            logger.warning(f"{obsid}   acis_streak_map failed")

    obsid_info = obs.get_info()
    if obsid_info["grating"] != "NONE":
        try:
            _make_grating_arm_mask(obs, inputs, outputs, logging_tag)
        except Exception:
            logger.warning(f"{obsid}   grating arm mask failed")


def _find_zeroth_order(evt_file: Path) -> tuple[float, float]:
    """Find the zero-order position by centroiding the brightest event cluster.

    Bins the event file coarsely to locate the dominant peak, then computes
    an event-weighted centroid within a window around that peak.  For grating
    observations the zero-order is the dominant source in the filtered event
    file; dispersed-arm events contribute negligible bias.

    This avoids tg_findzo (ACIS-only) and chandra_repro, and is more accurate
    than using the nominal pointing coordinates (RA_TARG/DEC_TARG) from the
    obspar, which can differ from the true zero-order by a few arcseconds.

    Parameters
    ----------
    evt_file
        Path to the filtered level-2 event file.

    Returns
    -------
    x, y
        Sky pixel coordinates of the zero-order.
    """
    with fits.open(evt_file) as hdus:
        data = hdus[1].data

    x = data["X"].astype(float)
    y = data["Y"].astype(float)

    # Coarse bin to locate the peak region.
    bin_size = 32
    x_min = int(x.min())
    y_min = int(y.min())
    nx = (int(x.max()) - x_min) // bin_size + 1
    ny = (int(y.max()) - y_min) // bin_size + 1
    xi = ((x - x_min) / bin_size).astype(int)
    yi = ((y - y_min) / bin_size).astype(int)
    coarse = np.zeros((ny, nx), dtype=np.intp)
    np.add.at(coarse, (yi, xi), 1)

    peak_yi, peak_xi = np.unravel_index(coarse.argmax(), coarse.shape)
    cx = x_min + (peak_xi + 0.5) * bin_size
    cy = y_min + (peak_yi + 0.5) * bin_size

    # Fine centroid: all events within 4 coarse bins of the peak.
    window = bin_size * 4
    in_window = (np.abs(x - cx) <= window) & (np.abs(y - cy) <= window)
    n_in = int(in_window.sum())
    if n_in < 5:
        raise ValueError(
            f"Too few events ({n_in}) within {window} px of peak "
            f"for zero-order centroid"
        )

    return float(x[in_window].mean()), float(y[in_window].mean())


def _make_grating_arm_mask(obs, inputs, outputs, logging_tag):
    """Create the grating arm mask via tg_create_mask.

    The zero-order position is found by centroiding the brightest event
    cluster in the filtered event file (see _find_zeroth_order), which
    gives the actual measured zero-order rather than the nominal pointing.
    """
    evt_file = inputs["events"][0]
    evt2_info = obs.get_evt2_info()
    ciao = obs.ciao

    zo_x, zo_y = _find_zeroth_order(evt_file)
    logger.debug(f"{obs.obsid}   zero-order from events: x={zo_x:.2f}, y={zo_y:.2f}")

    # Write a minimal position table that tg_create_mask accepts as input_pos_tab.
    # RA/DEC are included as metadata; tg_create_mask uses X/Y for positioning.
    zo_pos_table = table.Table(
        {
            "X": np.array([zo_x]),
            "Y": np.array([zo_y]),
            "RA": np.array([evt2_info["ra_targ"]]),
            "DEC": np.array([evt2_info["dec_targ"]]),
        }
    )
    zo_pos_file = (
        outputs["grating_arm_mask"].parent / f"{obs.obsid}_grating_zo_pos.fits"
    )
    zo_pos_table.write(zo_pos_file, overwrite=True)
    # tg_create_mask looks for the SRCLIST extension by name
    with fits.open(zo_pos_file, mode="update") as hdul:
        hdul[1].name = "SRCLIST"
        hdul.flush()

    ciao(
        "tg_create_mask",
        infile=str(evt_file),
        outfile=str(outputs["grating_arm_mask"]),
        input_pos_tab=f"{zo_pos_file}[SRCLIST]",
        grating_obs="header_value",
        clobber="yes",
        logging_tag=logging_tag,
    )


@task(
    name="wavdetect",
    inputs={
        "events": "primary/{obsid}_evt2_filtered.fits.gz",
        "image_file": "images/{obsid}_{band}_thresh.img",
        "psf_file": "images/{obsid}_{band}_thresh.psfmap",
        "exposure_file": "images/{obsid}_{band}_thresh.expmap",
    },
    outputs={
        "src": "sources/{obsid}_baseline.src",
        "cell": "sources/{obsid}_baseline.cell",
        "img": "sources/{obsid}_baseline.img",
        "nbkg": "sources/{obsid}_baseline.nbkg",
        "psf_size": "sources/{obsid}_psf_size_baseline.fits",
    },
    variables={
        "band": lambda obs: "wide" if obs.is_hrc else "broad",
    },
    download=(["evt2"]),
)
def wavdetect(obs, inputs, outputs):
    """
    Run wavdetect.
    """
    # possible parameter:
    scales = ["1.4", "2", "4", "8", "16", "32"]
    band = "wide" if obs.is_hrc else "broad"
    ecf = 0.9

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
            break
        except Exception:
            # if wavdetect fails, try again with fewer scales (but at least 3)
            scales = scales[:-1]
            if len(scales) < 3:
                raise

    sources = table.Table.read(outputs["src"])
    # CIAO psfsize_srcs chokes on an empty source list (gives misleading error:
    # "position value must be a filename or have a +, -, or comma separator")
    if len(sources) > 0:
        obs.ciao(
            "psfsize_srcs",
            inputs["events"],
            outputs["src"],
            outputs["psf_size"],
            f"energy={band}",
            f"ecf={ecf}",
            "clobber=yes",
        )

    return ReturnCode.OK


@task(
    name="gaussian_detect",
    inputs={
        "events": "primary/{obsid}_evt2_filtered.fits.gz",
        "src": "sources/{obsid}_celldetect.src",
    },
    optional_inputs={
        "psf_size": "sources/{obsid}_psf_size_celldetect.fits",
    },
    outputs={
        "src": "sources/{obsid}_gaussian_detect.src",
        "psf_size": "sources/{obsid}_psf_size_gaussian_detect.fits",
    },
    variables={
        "band": lambda obs: "wide" if obs.is_hrc else "broad",
    },
    # No download: the inputs above are derived products. filter_events is the
    # task that reads the raw event file, and declares that download itself.
)
def gaussian_detect(obs, inputs, outputs):
    return _fit_gaussian_sources(obs, inputs, outputs, seed_from_peak=False)


def _drop_crowded_seeds(input_sources):
    """Remove celldetect sources too close to another to trust a fit on.

    A gaussian fit is not attempted for any source whose nearest neighbour is
    within the crowding radius.

    A fit seeded inside another source's footprint risks blending with it: at
    obsid 7263, two celldetect sources 3.71" apart produced one fit that
    absorbed both sources' counts and one fit that failed to converge at all.
    Both outcomes are useless for astrometry, and every current selection
    already excludes near_neighbor_dist <= 6" downstream -- so this is work the
    pipeline was always going to discard, done before it is spent rather than
    after.

    Parameters
    ----------
    input_sources : astropy.table.Table
        Celldetect sources with ``COMPONENT``, ``y_angle`` and ``z_angle``
        columns, as read from the ``.src`` file.
    """
    if len(input_sources) == 0:
        return input_sources
    nnd = utils.get_near_neighbor_dist(input_sources, input_sources, id_col="COMPONENT")
    return input_sources[nnd > utils.NEAR_NEIGHBOR_DIST_ARCSEC]


def _drop_grating_arm_seeds(input_sources, obs):
    """Remove celldetect sources that fall on a grating observation's dispersed arm.

    Each arm knot is celldetect sampling one point of a continuous dispersed
    structure, not an independent source: at obsid 23308 (MKN421, LETG), all 87
    celldetect sources in the field land on the same two arms. A gaussian fit
    seeded there fits one knot as if it were isolated, which is exactly the
    kind of situation a large ``peak_offset`` already flags after the fact --
    this drops the ones celldetect itself can already tell are on the arm,
    before a fit is spent on them.

    A no-op on observations with no ``grating_arm_mask.fits`` (not a grating
    observation, or the mask step was skipped) -- see
    :meth:`Observation._on_grating_arm`, which this reuses directly.

    Parameters
    ----------
    input_sources : astropy.table.Table
        Celldetect sources with ``X`` and ``Y`` columns, as read from the
        ``.src`` file.
    obs : Observation
        Used to locate this obsid's ``grating_arm_mask.fits``.
    """
    if len(input_sources) == 0:
        return input_sources
    return input_sources[~obs._on_grating_arm(input_sources)]


def _drop_acis_streak_seeds(input_sources, obs):
    """Remove celldetect sources that fall on an ACIS readout streak.

    A readout streak is instrumental, not astrophysical: a fit seeded on one
    is fitting a CCD artifact, not a source. A no-op on observations with no
    ``acis_streaks.fits`` (not ACIS, or no source bright enough to produce
    one) -- see :meth:`Observation._on_acis_streak`, which this reuses
    directly.

    Parameters
    ----------
    input_sources : astropy.table.Table
        Celldetect sources with ``X`` and ``Y`` columns, as read from the
        ``.src`` file.
    obs : Observation
        Used to locate this obsid's ``acis_streaks.fits``.
    """
    if len(input_sources) == 0:
        return input_sources
    return input_sources[~obs._on_acis_streak(input_sources)]


def _seed_and_select_events(events_yag, events_zag, source, box_size, seed_from_peak):
    """Return ``(fit_source, event_mask)`` for one source.

    With `seed_from_peak`, the fit is re-seeded from the local emission peak instead
    of the celldetect position -- this avoids being pulled to an off-peak catalog or
    ICM centroid when the true emission maximum is offset -- and the event box is
    centred on that same peak so the fit window stays symmetric about its own seed.
    Anchoring the box to the celldetect position instead truncates the events
    asymmetrically, dropping the far side of the peak and biasing the fitted centroid
    back toward the position the re-seed was meant to escape.

    Parameters
    ----------
    events_yag, events_zag : np.ndarray
        Y- and Z-angle of all events in arcsec.
    source : dict or :any:`astropy.table.Row`
        Source to fit, with 'y_angle' and 'z_angle' in arcsec.
    box_size : float
        Half-width of the fit box in arcsec.
    seed_from_peak : bool
        Whether to re-seed from the local peak.

    Returns
    -------
    tuple
        ``(fit_source, event_mask)``. `fit_source` is `source` itself when
        `seed_from_peak` is False, otherwise a copy with the peak position.
    """
    fit_source = source
    center_yag = float(source["y_angle"])
    center_zag = float(source["z_angle"])

    if seed_from_peak:
        center_yag, center_zag = source_detection.find_local_peak(
            events_yag,
            events_zag,
            center_yag,
            center_zag,
            box_size=box_size,
        )
        fit_source = dict(source)
        fit_source["y_angle"] = center_yag
        fit_source["z_angle"] = center_zag

    event_mask = (np.abs(events_yag - center_yag) < box_size) & (
        np.abs(events_zag - center_zag) < box_size
    )
    return fit_source, event_mask


def _fit_gaussian_sources(  # noqa: PLR0915
    obs, inputs, outputs, seed_from_peak
):
    """Shared implementation of gaussian_detect and peak_gaussian_detect.

    The two tasks are identical except for the fit seed: gaussian_detect seeds
    from the celldetect position, peak_gaussian_detect re-seeds from the local
    image peak (``seed_from_peak=True``).
    """
    box_size = 4

    dtype = np.dtype(
        [
            ("COMPONENT", ">i4"),
            ("RA", ">f8"),
            ("DEC", ">f8"),
            ("X", ">f8"),
            ("Y", ">f8"),
            ("y_angle", ">f8"),
            ("z_angle", ">f8"),
            ("RA_ERR", ">f8"),
            ("DEC_ERR", ">f8"),
            ("X_ERR", ">f8"),
            ("Y_ERR", ">f8"),
            ("params", ">f8", (6,)),
            ("hess_inv", ">f8", (6, 6)),
            ("ndof", ">i8"),
            ("fit_ok", "?"),
            ("p_signal", ">f8"),
            ("sigma", ">f8", (2,)),
            ("rot_angle", ">f8"),
            ("sigma_y_angle", ">f8"),
            ("sigma_z_angle", ">f8"),
            ("corr_y_angle_z_angle", ">f8"),
            ("source_area", ">f8"),
            ("NET_COUNTS", ">i8"),
            ("signal", ">f8"),
            ("background", ">f8"),
            ("snr", ">f8"),
            ("ks_y_angle", ">f8"),
            ("ks_z_angle", ">f8"),
            ("ks_p_value_y_angle", ">f8"),
            ("ks_p_value_z_angle", ">f8"),
            ("ks_sign_y_angle", ">i8"),
            ("ks_sign_z_angle", ">i8"),
            ("ecf_radius", ">f8"),
            ("PSFRATIO", ">f8"),
            ("concentration_ratio", ">f8"),
            ("CORR_RA_DEC", ">f8"),
            ("CORR_X_Y", ">f8"),
        ]
    )

    input_sources = table.Table.read(inputs["src"])
    if len(input_sources) == 0:
        results = table.Table(dtype=dtype)
        results.write(outputs["src"], format="fits", overwrite=True)
        return ReturnCode.OK

    events = table.Table.read(inputs["events"], hdu=1)
    wcs = utils.get_wcs_from_fits_header(inputs["events"], hdu=1)
    events["RA"], events["DEC"] = wcs.pixel_to_world_values(events["x"], events["y"])

    obs_info = obs.get_info()
    att = Quat([obs_info["ra_nom"], obs_info["dec_nom"], obs_info["roll_nom"]])
    events["y_angle"], events["z_angle"] = radec_to_yagzag(
        events["RA"], events["DEC"], att
    )
    input_sources["y_angle"], input_sources["z_angle"] = radec_to_yagzag(
        input_sources["RA"], input_sources["DEC"], att
    )

    if not seed_from_peak:
        # Only for gaussian_detect. peak_gaussian_detect is never persisted as its
        # own detect_method (it exists to feed the peak_offset diagnostic) and its
        # rows never enter matching, so filtering its input would only lose
        # coverage -- exactly the coverage that made it possible to check, for a
        # crowded pair, whether re-seeding from the peak recovers a usable fit.
        input_sources = _drop_crowded_seeds(input_sources)
        input_sources = _drop_grating_arm_seeds(input_sources, obs)
        input_sources = _drop_acis_streak_seeds(input_sources, obs)
        if len(input_sources) == 0:
            results = table.Table(dtype=dtype)
            results.write(outputs["src"], format="fits", overwrite=True)
            return ReturnCode.OK

    events_yag = np.asarray(events["y_angle"])
    events_zag = np.asarray(events["z_angle"])

    results = []
    for source in input_sources:
        fit_source, sel = _seed_and_select_events(
            events_yag, events_zag, source, box_size, seed_from_peak
        )
        res = source_detection.fit_gaussian_2d(
            events[sel],
            fit_source,
            box_size=box_size,
        )
        results.append(res)

    results = table.Table(results)

    ecf = table.Table.read(inputs["psf_size"])
    pixel_size = 0.4920 if obs.is_acis else 0.13175
    results["ecf_radius"] = ecf["R"] * pixel_size
    results["PSFRATIO"] = (
        np.sqrt(results["sigma"][:, 0] * results["sigma"][:, 1]) / results["ecf_radius"]
    )

    results["concentration_ratio"] = [
        source_detection.concentration_ratio(
            events_yag,
            events_zag,
            float(row["y_angle"]),
            float(row["z_angle"]),
        )
        for row in results
    ]

    results.rename_column("n", "NET_COUNTS")

    results["RA"], results["DEC"] = yagzag_to_radec(
        results["y_angle"], results["z_angle"], att
    )
    results["X"], results["Y"] = wcs.world_to_pixel_values(
        results["RA"], results["DEC"]
    )

    # the uncertainties in yag/zag are given by the marginalized inverse Hessian matrix
    yag_zag_cov = np.asarray(results["hess_inv"][:, :2, :2])

    # and to propagate the uncertainties from yag/zag to ra/dec and x/y
    # linearize the transformations around the nominal attitude
    # I do this the dumb way: by finite differences

    d_arc = 1e-3
    ra_nom, dec_nom = yagzag_to_radec(0, 0, att)

    # RA/dec
    yagzag_to_radec_matrix = (
        np.vstack(
            [
                np.array(yagzag_to_radec(d_arc, 0, att)) - np.array([ra_nom, dec_nom]),
                np.array(yagzag_to_radec(0, d_arc, att)) - np.array([ra_nom, dec_nom]),
            ]
        ).T
        / d_arc
    )

    ra_dec_cov = yagzag_to_radec_matrix @ yag_zag_cov @ yagzag_to_radec_matrix.T

    results["RA_ERR"] = np.sqrt(ra_dec_cov[:, 0, 0])
    results["DEC_ERR"] = np.sqrt(ra_dec_cov[:, 1, 1])
    results["CORR_RA_DEC"] = ra_dec_cov[:, 0, 1] / (
        results["RA_ERR"] * results["DEC_ERR"]
    )

    # X/Y
    yagzag_to_pixel_matrix = np.vstack(
        [
            wcs.world_to_pixel_values(
                ra_nom + yagzag_to_radec_matrix.T[0, 0],
                dec_nom + yagzag_to_radec_matrix.T[0, 1],
            ),
            wcs.world_to_pixel_values(
                ra_nom + yagzag_to_radec_matrix.T[1, 0],
                dec_nom + yagzag_to_radec_matrix.T[1, 1],
            ),
        ]
    ).T

    pixel_cov = yagzag_to_pixel_matrix @ yag_zag_cov @ yagzag_to_pixel_matrix.T

    results["X_ERR"] = np.sqrt(pixel_cov[:, 0, 0])
    results["Y_ERR"] = np.sqrt(pixel_cov[:, 1, 1])
    results["CORR_X_Y"] = pixel_cov[:, 0, 1] / (results["X_ERR"] * results["Y_ERR"])

    units = {
        "RA": u.deg,
        "DEC": u.deg,
        "X": u.pixel,
        "Y": u.pixel,
        "y_angle": u.arcsec,
        "z_angle": u.arcsec,
        "RA_ERR": u.deg,
        "DEC_ERR": u.deg,
        "X_ERR": u.pixel,
        "Y_ERR": u.pixel,
    }
    for col, unit in units.items():
        results[col].unit = unit

    cols = [
        "COMPONENT",
        "RA",
        "DEC",
        "X",
        "Y",
        "y_angle",
        "z_angle",
        "RA_ERR",
        "DEC_ERR",
        "X_ERR",
        "Y_ERR",
    ]
    cols += [col for col in results.colnames if col not in cols]
    results = results[cols]

    results = results[results["fit_ok"]]

    if not outputs["src"].parent.exists():
        outputs["src"].parent.mkdir(exist_ok=True, parents=True)

    results.write(outputs["src"], format="fits", overwrite=True)

    shutil.copyfile(inputs["psf_size"], outputs["psf_size"])


@task(
    name="peak_gaussian_detect",
    inputs={
        "events": "primary/{obsid}_evt2_filtered.fits.gz",
        "src": "sources/{obsid}_celldetect.src",
    },
    optional_inputs={
        "psf_size": "sources/{obsid}_psf_size_celldetect.fits",
    },
    outputs={
        "src": "sources/{obsid}_peak_gaussian_detect.src",
        "psf_size": "sources/{obsid}_psf_size_peak_gaussian_detect.fits",
    },
    variables={
        "band": lambda obs: "wide" if obs.is_hrc else "broad",
    },
    # No download: the inputs above are derived products. filter_events is the
    # task that reads the raw event file, and declares that download itself.
)
def peak_gaussian_detect(obs, inputs, outputs):
    """
    Gaussian centroiding seeded from the local image peak rather than the celldetect position.

    For each celldetect source, the local brightness maximum in a smoothed event-density
    image is found within the fit box.  The Gaussian fit is then seeded from that peak
    instead of from the celldetect position.  This avoids being pulled to an off-peak
    catalog or ICM centroid when the true emission maximum is offset.

    Output columns are identical to gaussian_detect, so get_sources() works unchanged.
    """
    return _fit_gaussian_sources(obs, inputs, outputs, seed_from_peak=True)


@task(
    name="celldetect",
    inputs={
        "events": "primary/{obsid}_evt2_filtered.fits.gz",
        "image_file": "images/{obsid}_{band}_thresh.img",
        "psf_file": "images/{obsid}_{band}_thresh.psfmap",
    },
    outputs={
        "src": "sources/{obsid}_celldetect.src",
        # "psf_size": "sources/{obsid}_psf_size_celldetect.fits",
    },
    variables={
        "band": lambda obs: "wide" if obs.is_hrc else "broad",
    },
    download=(["evt2"]),
)
def celldetect(obs, inputs, outputs):
    """
    Run celldetect.
    """
    # possible parameter:
    snr = 3

    # Note that the maxlogicalwindow needs to be bigger than whatever the
    # binned filtered image is.
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
    name="celldetect_psf_size",
    inputs={
        "events": "primary/{obsid}_evt2_filtered.fits.gz",
        "src": "sources/{obsid}_celldetect.src",
    },
    outputs={
        "psf_size": "sources/{obsid}_psf_size_celldetect.fits",
    },
)
def celldetect_psf_size(obs, inputs, outputs):
    """
    Run celldetect.
    """
    band = "wide" if obs.is_hrc else "broad"
    ecf = 0.9

    sources = table.Table.read(inputs["src"])
    # it seems CIAO chokes on an empty source list
    if len(sources) > 0:
        obs.ciao(
            "psfsize_srcs",
            inputs["events"],
            inputs["src"],
            outputs["psf_size"],
            f"energy={band}",
            f"ecf={ecf}",
            "clobber=yes",
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

    pixel = 0.13180 if obs.is_hrc else 0.5

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
