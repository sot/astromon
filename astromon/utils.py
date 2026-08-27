import functools
import logging
import os
import re
import subprocess
import tempfile
import warnings
from contextlib import AbstractContextManager
from pathlib import Path

import numpy as np
import Ska.Shell
from astropy.io import fits
from astropy.table import Table, hstack, join, vstack
from astropy.wcs import WCS, FITSFixedWarning
from cxotime import CxoTime
from Quaternion import Quat

__all__ = [
    "communicate",
    "Ciao",
    "chdir",
    "ciao_context",
    "logging_call_decorator",
    "calalign_from_files",
    "get_calalign_offsets",
    "get_latest_calalign_matrix",
    "get_rebased_offsets",
    "get_wcs_from_fits_header",
    "get_near_neighbor_dist",
    "NEAR_NEIGHBOR_DIST_ARCSEC",
    "ASTROMON_DATA_DIR",
]


# The one near_neighbor_dist threshold used everywhere it matters: excluding a
# source from matching in cross_match.simple_cross_match's default args, and
# excluding a celldetect source from even being handed to a gaussian fit in
# observation._drop_crowded_seeds. One constant so the two cannot drift apart --
# they are answering the same question ("is this position crowded?") at two
# different points in the pipeline, and disagreeing between them defeats the
# point of filtering early.
NEAR_NEIGHBOR_DIST_ARCSEC = 6.0


# Root of astromon's own data: the locally cached whole catalogs and the archive
# of downloaded observations. Defined here, in the module both cross_match and
# observation import, so there is one answer rather than one per consumer.
#
# ASTROMON_DATA_DIR isolates it the way ASTROMON_FILE isolates the database.
# Without it the only way to move this tree was to move SKA itself, which takes
# mica, CALDB and everything else along -- so a development run could point at
# its own database and still read, refresh and overwrite the shared catalogs and
# archive into the shared tree.
if "ASTROMON_DATA_DIR" in os.environ:
    ASTROMON_DATA_DIR = Path(os.environ["ASTROMON_DATA_DIR"])
else:
    ASTROMON_DATA_DIR = (
        Path(os.environ.get("SKA", Path.home() / "ska")) / "data" / "astromon"
    )


CIAO_ENV = {}


class MissingTableException(Exception):
    """
    Exception class in case a table is missing in the DB file.
    """


class CiaoProcessFailure(RuntimeError):
    def __init__(self, *args, lines=(), **kwargs):
        super().__init__(*args, **kwargs)
        self.lines = tuple(lines)


def communicate(process, logger=None, level="WARNING", text=False, logging_tag=""):
    """
    Real-time reading of a subprocess stdout.

    Parameters
    ----------
    process:
        process returned by :any:`python:subprocess.Popen`
    logger: :any:`python:logging.Logger`
    level: str or int
        the logging level
    text: bool
        whether the process output is text or binary
    """
    if type(level) is str:
        level = getattr(logging, level)
    if logger is None:
        logger = logging.getLogger()

    lines = []
    while True:
        if process.poll() is not None:
            break
        line = process.stdout.readline()
        line = line if text else line.decode()
        line = line.rstrip("\n")
        logger.log(level, f"{logging_tag} {line}")
        if line:
            lines.append(line)

    # in case the buffer is still not empty after the process ended
    for line in process.stdout.readlines():
        line = line if text else line.decode()  # noqa: PLW2901
        line = line.rstrip("\n")  # noqa: PLW2901
        logger.log(level, f"{logging_tag} {line}")
        if line:
            lines.append(line)

    return lines


class Ciao:
    """
    Encapsulate calls to CIAO tools.

    This class does two things:

    * Set the CIAO context (PFILES and ASCDS_WORK_PATH).
    * Handle calls to subprocesses.

    Parameters
    ----------
    prefix : str
        The location of CIAO.
    workdir : :any:`python:pathlib.Path` or str
        Working directory. Used to set PFILES and ASCDS_WORK_PATH.
    logger: :any:`python:logging.Logger`
        If not provided, the root logger is used.
    """

    def __init__(self, prefix=None, workdir=None, logger=None):
        if prefix is None:
            prefix = "/soft/ciao"
        if workdir is None:
            self.workdir = tempfile.TemporaryDirectory()
        else:
            self.workdir = Path(workdir).absolute
        if logger is None:
            self.logger = logging
        elif type(logger) is str:
            self.logger = logging.getLogger(logger)
        else:
            self.logger = logger
        prefix = Path(prefix)
        if not prefix.exists():
            raise FileNotFoundError(
                f"CIAO prefix {prefix} does not exist. Install CIAO there, or pass a "
                "different ciao_prefix (--ciao-prefix on the command line)."
            )

        self.env = CIAO_ENV.get(
            prefix, Ska.Shell.getenv(f"source {prefix}/bin/ciao.sh")
        ).copy()

        if "ASCDS_INSTALL" not in self.env:
            # conda-packaged CIAO: ciao.sh is a no-op so getenv returns an
            # empty delta.  Replicate what the conda activation scripts set
            # (etc/conda/activate.d/*.sh) so CIAO tools can find CALDB etc.
            self.env = dict(os.environ)
            self.env["ASCDS_INSTALL"] = str(prefix)
            self.env["ASCDS_CONTRIB"] = str(prefix)
            self.env["ASCDS_CALIB"] = str(prefix / "data")
            caldb = str(prefix / "CALDB")
            self.env["CALDB"] = caldb
            self.env["CALDBCONFIG"] = caldb + "/software/tools/caldb.config"
            self.env["CALDBALIAS"] = caldb + "/software/tools/alias_config.fits"
            self.env["PATH"] = (
                str(prefix / "bin") + ":" + self.env.get("PATH", "/usr/bin:/bin")
            )

        # Cache a copy, not the live instance dict: the workdir block below writes
        # ASCDS_WORK_PATH and PFILES into self.env, and those are specific to this
        # observation. Caching the live dict let them leak into every later Ciao for
        # the same prefix -- including one constructed with no workdir, which would
        # then inherit a param directory that has since been removed.
        CIAO_ENV[prefix] = dict(self.env)

        if workdir is not None:
            workdir = Path(workdir)
            workdir.mkdir(parents=True, exist_ok=True)
            self.env["ASCDS_WORK_PATH"] = str(workdir)
            if ":" in str(workdir.absolute()):
                raise RuntimeError("CIAO workdir cannot contain colon")
            ascds = self.env["ASCDS_INSTALL"]
            pf_dirs = [ascds + "/param"]
            contrib_param = ascds + "/contrib/param"
            if Path(contrib_param).exists():
                pf_dirs.append(contrib_param)
            self.env["PFILES"] = "{};{}".format(
                str(workdir.absolute()), ":".join(pf_dirs)
            )

    def __call__(
        self,
        name,
        *args,
        logging_tag="",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        **kwargs,
    ):
        """
        Run a CIAO command.
        """
        args = [str(a) for a in args]
        kwargs = {k: str(v) for k, v in kwargs.items()}
        cmd = [name] + args
        for k in kwargs:
            cmd.append(f"{k}={kwargs[k]}")
        self.logger.info(logging_tag + " " + " ".join(cmd))
        process = subprocess.Popen(cmd, stdout=stdout, stderr=stderr, env=self.env)
        lines = communicate(
            process, self.logger, level=20, logging_tag=logging_tag
        )  # 20 is logging.INFO
        if process.returncode:
            raise CiaoProcessFailure(f"{logging_tag} {name} failed", lines=lines)

    def pget(self, command, param="value", logging_tag=""):
        """
        Call CIAO's pget and return the standard output.
        """
        self.logger.info(f"{logging_tag} pget {command} {param}")
        out = (
            subprocess.check_output(["pget", command, param], env=self.env)
            .decode()
            .strip()
        )
        self.logger.info(f"{logging_tag} pget {out}")
        return (
            subprocess.check_output(["pget", command, param], env=self.env)
            .decode()
            .strip()
        )


def set_ciao_context(directory):
    """
    Set the environment variables used by CIAO (PFILES and ASCDS_WORK_PATH).
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    os.environ["ASCDS_WORK_PATH"] = str(directory)
    if ":" in str(directory.absolute()):
        raise RuntimeError("CIAO context directory cannot contain colon")
    pf = "{};{}:{}".format(
        str(directory.absolute()),
        os.environ["ASCDS_INSTALL"] + "/param",
        os.environ["ASCDS_INSTALL"] + "/contrib/param",
    )
    os.environ["PFILES"] = pf


class chdir(AbstractContextManager):
    """
    Context manager to temporarily change to a directory.
    """

    def __init__(self, directory):
        self.cwd = os.getcwd()
        self.directory = directory
        os.chdir(self.directory)

    def __exit__(self, exc_type, exc_value, exc_traceback):
        os.chdir(self.cwd)
        return exc_type is None


class CIAOContext(AbstractContextManager):
    """
    Context manager to set CIAO context (PFILES and ASCDS_WORK_PATH)

    Parameters
    ----------
    directory : :any:`python:pathlib.Path` or str
        Working directory. Used to set PFILES and ASCDS_WORK_PATH.
    logger: :any:`python:logging.Logger`
        If not provided, the root logger is used.
    """

    def __init__(self, directory, logger=None):
        self.directory = directory
        self.logger = logging if logger is None else logger
        self.logger.debug(f"entering ciao context {self.directory}")
        self.env = {k: os.environ.get(k, None) for k in ["ASCDS_WORK_PATH", "PFILES"]}
        set_ciao_context(self.directory)

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.logger.debug(f"leaving ciao context {self.directory}")
        os.environ.update(
            {key: val for key, val in self.env.items() if val is not None}
        )
        return exc_type is None


def ciao_context_function(func):
    """
    Decorator to set the CIAO parameters directory within the context of a member function.

    By default, CIAO stores function parameters in $HOME/cxcds_param4. Within a function decorated
    with this decorator, this changes to ${self.workdir}/${self.obsid}/params by setting the PFILES
    variable. This also sets ASCDS_WORK_PATH.

    The decorated class must have a 'workdir' and an 'obsid' data members.
    """

    @functools.wraps(func)
    def ciao_context_wrapper(self, *args, **kwargs):
        with CIAOContext(
            self.workdir / self.obsid / "params",
            logger=logging.getLogger("astromon.utils"),
        ):
            return func(self, *args, **kwargs)

    return ciao_context_wrapper


ciao_context = CIAOContext
"""
Context manager to set CIAO context (PFILES and ASCDS_WORK_PATH)


Parameters
----------
directory : pathlib.Path or str
    Working directory. Used to set PFILES and ASCDS_WORK_PATH.
logger: :any:`python:logging.Logger`
    If not provided, the root logger is used.
"""


class LoggingCallDecorator:
    """
    Decorator to add logging messages at the start/end of the decorated function.
    """

    def __init__(self, logger=None, log_level=logging.DEBUG, ignore_exceptions=False):
        self.logger = logger if logger is not None else logging
        self.log_level = log_level
        self.ignore_exceptions = ignore_exceptions

    def __call__(self, func):
        @functools.wraps(func)
        def _logging_call_decorator(*args, **kwargs):
            args_str = []
            instance = ""
            if args:
                instance = f"{args[0]} "
                args_str += [", ".join([str(arg) for arg in args[1:]])]
            if kwargs:
                args_str += [", ".join([f"{k}={v}" for k, v in kwargs.items()])]
            args_str = ", ".join(args_str)
            self.logger.log(
                self.log_level, f"{instance}{func.__name__}({args_str}) Started"
            )
            try:
                result = func(*args, **kwargs)
                self.logger.log(
                    self.log_level, f"{instance}{func.__name__}({args_str}) Finished"
                )
                return result
            except Exception as e:
                self.logger.log(
                    self.log_level, f"{instance}{func.__name__}({args_str}) Except: {e}"
                )
                if not self.ignore_exceptions:
                    raise

        return _logging_call_decorator


logging_call_decorator = LoggingCallDecorator(
    log_level=logging.INFO, logger=logging.getLogger("astromon")
)
"""
Decorator to add logging messages at the start/end of the decorated function.
"""


def calalign_from_files(calalign_dir=None):
    """
    Get data from calalign files in `calalign_dir`.

    Parameters
    ----------
    calalign_dir: :any:`pathlib.Path`
        Directory where to find the calalign files.

    Returns
    -------
    :any:`astropy.table.Table`
    """

    if calalign_dir is None:
        calalign_dir = "/data/caldb/data/chandra/pcad/align"

    calalign_files = sorted(Path(calalign_dir).glob("*.fits"))
    if not calalign_files:
        raise Exception(f"CALALIGN files does not exist at {calalign_dir}")

    caldb_info = {
        "N0008": {
            "CalDB": "4.4.4",
            "since": "2011-06-28T20:45:00",
        },
        "N0009": {
            "CalDB": "4.6.2",
            "since": "2014-07-09T21:00:00",
        },
        "N0010": {
            "CalDB": "4.10.0",
            "since": "2022-06-28T14:00:00",
        },
    }

    transforms = []
    for filename in sorted(calalign_files):
        # pcadD2001-01-01alignN0008.fits
        file_re = re.search(
            "(?P<date>[0-9]{4}-[0-9]{2}-[0-9]{2})align(?P<version>N[0-9]{4}).fits",
            filename.name,
        )
        if not file_re:
            raise Exception(f"filename {filename} does not conform to expected format")

        info = file_re.groupdict()
        hdus = fits.open(filename)
        cals = hdus[1].data
        cvsd = (
            hdus[1].header["CVSD0001"]
            if "CVSD0001" in hdus[1].header
            else "2050:001:00:00:00.00"
        )
        cved = (
            hdus[1].header["CVED0001"]
            if "CVED0001" in hdus[1].header
            else "2050:001:00:00:00.00"
        )
        transforms.extend(
            [
                (
                    {
                        "start": CxoTime(cvsd),
                        "stop": CxoTime(cved),
                        "detector": cal["INSTR_ID"].strip(),
                        "date": CxoTime(info["date"]),
                        "version": info["version"],
                        "caldb_version": caldb_info[info["version"]]["CalDB"],
                        "since": caldb_info[info["version"]]["since"],
                        "aca_misalign": cal["ACA_MISALIGN"],
                        "fts_misalign": cal["FTS_MISALIGN"],
                    }
                )
                for cal in cals
            ]
        )

    calalign = Table(transforms)
    calalign["dy"], calalign["dz"] = get_offsets(calalign["aca_misalign"])
    return calalign


def get_offsets(aca_misalign):
    """
    Get the yag/zag offsets from the aca_misalign matrix.

    Parameters
    ----------
    aca_misalign: np.array
        The ACA misalignment matrices stored in calalign files as "aca_misalign"

    Returns
    -------
    (np.array, np.array)
        Two arrays with shape (aca_misalign.shape[1:])
    """

    aca_misalign = np.atleast_1d(aca_misalign)
    dys = []
    dzs = []
    for xform in aca_misalign:
        q = Quat(transform=xform)
        dzs.append(-q.pitch * 3600)
        dys.append(q.yaw * 3600)
    return np.array(dys), np.array(dzs)


def get_calalign_offsets(all_matches, ref_calalign=None, calalign_dir=None):
    """
    Get a table with the yag/zag offsets subtracted by the aca_misalign matrix.

    The table includes two sets of offsets: the offsets used when the observation was processed,
    and the offsets one would get when using a reference CALALIGN.

    Parameters
    ----------
    all_matches: :any:`astropy.table.Table`
        The table of x-ray sources. The required columns are 'obsid' and 'x_id'.
        If a 'detect_method' column is present, it is included in the match
        identifier: 'x_id' is only unique within a given (obsid, detect_method)
        pair (celldetect and gaussian_detect each number their own sources
        starting from scratch for a given obsid), so a table that mixes
        detect methods needs 'detect_method' to tell those sources apart.

    Returns
    -------
    :any:`astropy.table.Table`
        A table with keys:

        - calalign_dy
        - calalign_dz
        - aca_misalign
        - fts_misalign
        - calalign_version
        - ref_calalign_dy
        - ref_calalign_dz
        - ref_aca_misalign
        - ref_fts_misalign
        - ref_calalign_version
    """
    # obsid and x_id identify a match, unless detect_method is present, in which
    # case x_id is only unique within a given (obsid, detect_method) pair.
    match_id_keys = ["obsid", "x_id"]
    if "detect_method" in all_matches.colnames:
        match_id_keys = match_id_keys + ["detect_method"]

    # this is a table of all calalign matrices
    calalign = calalign_from_files(calalign_dir=calalign_dir)
    calalign.rename_columns(
        ["dy", "dz", "caldb_version"],
        ["calalign_dy", "calalign_dz", "calalign_version"],
    )

    # steps:
    # - join on detector
    # - filter out calaligns that do not correspond to the time of the match
    # - filter out calalign versions after the desired one
    # - group by match (match_id_keys identify a match uniquely)
    # - take the row corresponding to the latest calibration
    match_w_cal = join(all_matches, calalign, keys=["detector"])
    match_w_cal = match_w_cal[
        (match_w_cal["start"] <= match_w_cal["time"])
        & (match_w_cal["time"] < match_w_cal["stop"])
    ]
    match_w_cal["cdbv"] = [
        [int(f) for f in s.split(".") if f.isdigit()]
        for s in match_w_cal["caldb_version"]
    ]
    match_w_cal["cav"] = [
        [int(f) for f in s.split(".") if f.isdigit()]
        for s in match_w_cal["calalign_version"]
    ]

    # I actually would prefer to not sort the whole table
    # what I want is to sort the groups, which are small
    match_w_cal.sort(keys=match_id_keys + ["cav", "start"], reverse=True)

    if ref_calalign is None:
        ref_calalign = max(
            [
                [int(f) for f in v.split(".")]
                for v in np.unique(calalign["calalign_version"])
            ]
        )
    else:
        ref_calalign = [int(f) for f in ref_calalign.split(".")]

    # the actual calalign used must be a CalDB version no later than one used to process the OBSID
    # Use tuple() for both sides so the comparison is lexicographic (returns a scalar bool),
    # not element-wise (which would return a numpy bool array and break the list comprehension).
    actual = match_w_cal[[tuple(r["cdbv"]) >= tuple(r["cav"]) for r in match_w_cal]]
    actual = actual.group_by(match_id_keys)
    actual = actual[actual.groups.indices[:-1]]

    # the reference calalign must be a CalDB version no later than the given CalDB reference
    reference = match_w_cal[
        [tuple(ref_calalign) >= tuple(r["cav"]) for r in match_w_cal]
    ]
    reference = reference.group_by(match_id_keys)
    reference = reference[reference.groups.indices[:-1]]

    cols = [
        "calalign_dy",
        "calalign_dz",
        "aca_misalign",
        "fts_misalign",
        "calalign_version",
    ]
    ref_cols = ["ref_" + c for c in cols]
    reference.rename_columns(cols, ref_cols)

    # these should be true here:
    if len(all_matches) != len(actual):
        raise RuntimeError("len(all_matches) != len(actual)")
    if len(all_matches) != len(reference):
        raise RuntimeError("len(all_matches) != len(reference)")
    for key in match_id_keys:
        if np.all(reference[key] != actual[key]):
            raise RuntimeError(f"reference.{key} != actual.{key}")

    result = join(
        all_matches[match_id_keys],
        hstack([actual[match_id_keys + cols], reference[ref_cols + ["since"]]]),
        keys=match_id_keys,
    )

    for key in match_id_keys:
        if np.all(all_matches[key] != result[key]):
            raise RuntimeError(f"all_matches.{key} != result.{key}")

    # these are observations that happened after the reference calalign was added
    result["after_caldb"] = result["since"] < all_matches["time"]

    return result


def get_latest_calalign_matrix(calalign_dir=None):
    """
    Get the single most-recently-dated CALALIGN alignment matrix for each detector.

    CALDB packages many distinct alignment solutions -- issued periodically to track
    slow periscope drift -- under the same version label (e.g. every file from
    2013-01-19 through 2021-07-02 in a typical CALDB checkout is labeled "N0010"), so
    "the latest CalDB version" and "the latest matrix" are not the same thing: see
    :any:`get_calalign_offsets`'s ``ref_calalign_dy``/``ref_calalign_dz``, which pick
    a reference by version label and can land on a much older matrix than the most
    recent one as a result. This instead picks the single most-recently-dated file per
    detector, by its CVSD0001 start time, ignoring the version label entirely.

    Parameters
    ----------
    calalign_dir: :any:`pathlib.Path`
        Directory where to find the calalign files. Passed straight through to
        :any:`calalign_from_files`.

    Returns
    -------
    dict
        {detector: (dy, dz)} for the most-recently-dated matrix per detector.
    """
    calalign_files = calalign_from_files(calalign_dir=calalign_dir)
    latest = {}
    for detector in np.unique(calalign_files["detector"]):
        rows = calalign_files[calalign_files["detector"] == detector]
        row = rows[np.argmax(rows["start"])]
        latest[detector] = (row["dy"], row["dz"])
    return latest


def get_rebased_offsets(all_matches, ref_matrix=None, calalign_dir=None):
    """
    Rebase dy/dz to a single fixed CALALIGN alignment matrix.

    Each observation's dy/dz has whatever CALALIGN alignment was in effect when it
    was processed baked in, so a dy/dz time series mixes real signal with step
    changes purely from CALALIGN updates. This moves every match onto one fixed
    reference matrix instead, so the whole mission reads as one continuous signal:
    ``dy_rebased = dy - (calalign_dy - ref_dy)``, per detector.

    This is the rebase method used throughout the ``absolute_astrometry`` project
    (e.g. ``celmon-story-1-source-rebase.ipynb``, independently checked there against
    real ASCDS pipeline reprocessing). It is deliberately *not* built on
    :any:`get_calalign_offsets`'s ``ref_calalign_dy``/``ref_calalign_dz`` columns
    (used by ``scripts/calc_astrometry_rms.py``'s ``--use-latest-calalign`` and
    ``web/celmon.py``'s ``use_reference_calalign``): those pick a reference by CalDB
    *version* label, which can be a much older matrix than the most recent one --
    see :any:`get_latest_calalign_matrix`, this function's default reference when
    `ref_matrix` is not given. Pass your own `ref_matrix` to pin a specific epoch
    instead of always tracking the newest one.

    Parameters
    ----------
    all_matches: :any:`astropy.table.Table`
        Table of x-ray/catalog matches. Required columns: 'obsid', 'x_id', 'detector',
        'time', 'dy', 'dz', 'caldb_version'. Include 'detect_method' if `all_matches`
        spans more than one detection method (celldetect and gaussian_detect each
        number 'x_id' independently per obsid) -- see :any:`get_calalign_offsets`.
    ref_matrix: dict
        {detector: (dy, dz)}, e.g. from :any:`get_latest_calalign_matrix`. Computed
        automatically from `calalign_dir` if not given.
    calalign_dir: :any:`pathlib.Path`
        Directory where to find the calalign files.

    Returns
    -------
    :any:`astropy.table.Table`
        A copy of `all_matches` with two extra columns, 'dy_rebased'/'dz_rebased'.
        Rows with `caldb_version == "0.0"` (no real CalDB version recorded, so there
        is no calalign entry to rebase from) get NaN instead of being dropped.
    """
    if ref_matrix is None:
        ref_matrix = get_latest_calalign_matrix(calalign_dir=calalign_dir)

    # caldb_version "0.0" means no real CalDB version was recorded for that
    # observation -- get_calalign_offsets has no calalign entry to match it against.
    has_caldb = all_matches["caldb_version"] != "0.0"

    no_caldb = all_matches[~has_caldb].copy()
    no_caldb["dy_rebased"] = np.nan
    no_caldb["dz_rebased"] = np.nan

    sub = all_matches[has_caldb].copy()
    if len(sub) == 0:
        return no_caldb

    cal = get_calalign_offsets(sub, calalign_dir=calalign_dir)
    ref_dy = np.array([ref_matrix[d][0] for d in sub["detector"]])
    ref_dz = np.array([ref_matrix[d][1] for d in sub["detector"]])
    sub["dy_rebased"] = sub["dy"] - (cal["calalign_dy"] - ref_dy)
    sub["dz_rebased"] = sub["dz"] - (cal["calalign_dz"] - ref_dz)

    return vstack([sub, no_caldb]) if len(no_caldb) > 0 else sub


def get_wcs_from_fits_header(filename=None, hdu=0, header=None):
    if filename is None and header is None:
        raise ValueError("Either filename or header must be provided")
    if filename is not None and header is not None:
        raise ValueError("Either filename or header must be provided, not both")

    if header is None:
        hdus = fits.open(filename)
        header = hdus[hdu].header

    ctypes = {val: key for key, val in header["TCTYP*"].items()}

    ra_col = ctypes.pop("RA---TAN", None)
    dec_col = ctypes.pop("DEC--TAN", None)

    c1 = (
        m.group(1)
        if ra_col is not None and (m := re.search(r"(\d+)$", ra_col))
        else None
    )
    c2 = (
        m.group(1)
        if dec_col is not None and (m := re.search(r"(\d+)$", dec_col))
        else None
    )

    if c1 is None or c2 is None:
        raise ValueError("Could not determine WCS columns")

    params = {
        "CTYPE1": header[f"TCTYP{c1}"],
        "CRVAL1": header[f"TCRVL{c1}"],
        "CRPIX1": header[f"TCRPX{c1}"] + 1,
        "CDELT1": header[f"TCDLT{c1}"],
        "CUNIT1": header[f"TCUNI{c1}"],
        "CTYPE2": header[f"TCTYP{c2}"],
        "CRVAL2": header[f"TCRVL{c2}"],
        "CRPIX2": header[f"TCRPX{c2}"] + 1,
        "CDELT2": header[f"TCDLT{c2}"],
        "CUNIT2": header[f"TCUNI{c2}"],
    }
    with warnings.catch_warnings():
        # For some reason FITS/WCS seems to think many of the CXC header
        # keywords are non-standard and need to be fixed.
        warnings.simplefilter("ignore", category=FITSFixedWarning)
        wcs = WCS(params)
    return wcs


def get_near_neighbor_dist(sources, all_sources=None, id_col="COMPONENT"):
    """
    Calculate the distance from each source in `sources` to the nearest source in `all_sources`.

    The distance is the euclidean distance in the y_angle/z_angle plane.

    If `all_sources` is None, `sources` is used for both.

    Parameters
    ----------
    sources : astropy.table.Table
        The sources to calculate the nearest neighbor distance for.
    all_sources : astropy.table.Table
        The sources to use as potential neighbors.
    id_col : str
        The column to use as the unique identifier for each source.
    """
    all_sources = sources if all_sources is None else all_sources

    src1 = sources.as_array()[None]
    src2 = all_sources.as_array()[:, None]
    distance = np.sqrt(
        (src1["y_angle"] - src2["y_angle"]) ** 2
        + (src1["z_angle"] - src2["z_angle"]) ** 2
    )
    # add a large value to the distance where the sources are the same
    # (a source should not be its own closest neighbor)
    distance += np.where(src1[id_col] == src2[id_col], np.inf, 0).reshape(
        distance.shape
    )
    return np.min(distance, axis=0)
