import functools
import logging
import os
import subprocess
import tempfile
from contextlib import AbstractContextManager
from pathlib import Path

import numpy as np
import Ska.Shell
from astropy.table import join, hstack

__all__ = [
    "communicate",
    "Ciao",
    "chdir",
    "ciao_context",
    "logging_call_decorator",
    "FlowException",
    "calalign_from_files",
    "get_calalign_offsets",
]


CIAO_ENV = {}


class MissingTableException(Exception):
    """
    Exception class in case a table is missing in the DB file.
    """

    pass


class FlowException(Exception):
    """
    Exception class to interrupt the execution flow.

    This exception class is used by :any:`logging_call_decorator` and is silently ignored unless
    instructed otherwise.
    """

    pass


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

    while True:
        if process.poll() is not None:
            break
        line = process.stdout.readline()
        line = line if text else line.decode()
        line = line.rstrip("\n")
        logger.log(level, f"{logging_tag} {line}")

    # in case the buffer is still not empty after the process ended
    for line in process.stdout.readlines():
        line = line if text else line.decode()
        line = line.rstrip("\n")
        logger.log(level, f"{logging_tag} {line}")


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
            raise Exception(f"CIAO prefix {prefix} does not exist")

        self.env = CIAO_ENV.get(
            prefix, Ska.Shell.getenv(f"source {prefix}/bin/ciao.sh")
        ).copy()

        if workdir is not None:
            workdir = Path(workdir)
            workdir.mkdir(parents=True, exist_ok=True)
            self.env["ASCDS_WORK_PATH"] = str(workdir)
            assert ":" not in str(
                workdir.absolute()
            ), "CIAO workdir cannot contain colon"
            pf = "{};{}:{}".format(
                str(workdir.absolute()),
                self.env["ASCDS_INSTALL"] + "/param",
                self.env["ASCDS_INSTALL"] + "/contrib/param",
            )
            self.env["PFILES"] = pf

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
        communicate(
            process, self.logger, level=20, logging_tag=logging_tag
        )  # 20 is logging.INFO
        if process.returncode:
            raise Exception(f"{logging_tag} {name} failed")

    def pget(self, command, param="value", logging_tag=""):
        """
        Call CIAO's pget and return the standard output.
        """
        self.logger.info(f"{logging_tag} pget {command} {param}")
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
    assert ":" not in str(
        directory.absolute()
    ), "CIAO context directory can not contain colon"
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
            except FlowException:
                if not self.ignore_exceptions:
                    raise
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
    import re
    from astropy.io import fits
    from cxotime import CxoTime
    from astropy.table import Table

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
        for cal in cals:
            transforms.append(
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
    from Quaternion import Quat
    import numpy as np

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
        The table of x-ray sources. The only required columns are 'obsid' and 'x_id'

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
    # - group by match (obsid and x_id are unique identifiers)
    # - take the row corresponding to the latest calibration
    match_w_cal = join(all_matches, calalign, keys=["detector"])
    match_w_cal = match_w_cal[
        (match_w_cal["start"] <= match_w_cal["time"])
        & (match_w_cal["time"] < match_w_cal["stop"])
    ]
    match_w_cal["cdbv"] = [
        [int(f) for f in s.split(".")] for s in match_w_cal["caldb_version"]
    ]
    match_w_cal["cav"] = [
        [int(f) for f in s.split(".")] for s in match_w_cal["calalign_version"]
    ]

    # I actually would prefer to not sort the whole table
    # what I want is to sort the groups, which are small
    match_w_cal.sort(keys=["obsid", "x_id", "cav", "start"], reverse=True)

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
    # NOTE: cdbv's dtype is "object" because versions have different lengthd
    # but cav's dtype is array, so need to convert to list
    actual = match_w_cal[[r["cdbv"] >= r["cav"].tolist() for r in match_w_cal]]
    actual = actual.group_by(["obsid", "x_id"])
    actual = actual[actual.groups.indices[:-1]]

    # the reference calalign must be a CalDB version no later than the given CalDB reference
    reference = match_w_cal[[ref_calalign >= r["cav"].tolist() for r in match_w_cal]]
    reference = reference.group_by(["obsid", "x_id"])
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
    assert len(all_matches) == len(actual)
    assert len(all_matches) == len(reference)
    assert np.all(reference["obsid"] == actual["obsid"])
    assert np.all(reference["x_id"] == actual["x_id"])

    result = join(
        all_matches[["obsid", "x_id"]],
        hstack([actual[["obsid", "x_id"] + cols], reference[ref_cols + ["since"]]]),
        keys=["obsid", "x_id"],
    )

    assert np.all(all_matches["obsid"] == result["obsid"])
    assert np.all(all_matches["x_id"] == result["x_id"])

    # these are observations that happened after the reference calalign was added
    result["after_caldb"] = result["since"] < all_matches["time"]

    return result
