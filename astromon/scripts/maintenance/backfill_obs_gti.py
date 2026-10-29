"""Fill astromon_obs_gti from the archive evt2 files in local work directories.

Why this is needed
------------------
astromon_obs_gti holds the time span each observation's event data cover: the
first GTI START and last GTI STOP over every GTI extension of the archive evt2.
The pipeline writes it from the version that added the table on; obsids
processed before that have no row. Their work directories usually still hold
the archive evt2, which is all this needs.

For every obsid in astromon_obs, this reads the archive evt2 in
``<workdir>/obsNN/<obsid>/primary/`` -- never ``<obsid>_evt2_filtered.fits.gz``
next to it -- headers and GTI extensions only, and saves one row per obsid.
Obsids with no work directory, no archive evt2 (e.g. removed by
Observation.cleanup_downloads) or an unreadable one get no row and are listed in
the log. Rows for obsids it does not fill are left as they are.

Usage
-----
::

    python -m astromon.scripts.maintenance.backfill_obs_gti \\
        --db astromon_copy.h5 --workdir /Volumes/Black/data/astromon/work \\
        --processes 4
"""

import argparse
import logging
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import NamedTuple

import numpy as np
from astropy.table import Table

from astromon import db, utils

logger = logging.getLogger("backfill_obs_gti")

#: What happened to each obsid, in the order they are reported.
OUTCOMES = ("filled", "no_archive_evt2", "no_work_dir", "unreadable")


class SpanReading(NamedTuple):
    """One obsid's outcome: its span when `outcome` is "filled", else why not."""

    obsid: int
    outcome: str
    gti_start: float = np.nan
    gti_stop: float = np.nan
    reason: str = ""


def obsid_work_dir(workdir: Path, obsid: int) -> Path:
    """`obsid`'s directory under `workdir`, laid out the way Observation does it."""
    return Path(workdir) / f"obs{obsid // 1000:02d}" / str(obsid)


def read_span(workdir: Path, obsid: int) -> SpanReading:
    """Read `obsid`'s GTI span from the archive evt2 in its work directory."""
    obsid_dir = obsid_work_dir(workdir, obsid)
    if not obsid_dir.is_dir():
        return SpanReading(obsid, "no_work_dir")
    try:
        event_file = utils.find_archive_evt2(
            list((obsid_dir / "primary").glob("*_evt2.fits*"))
        )
        if event_file is None:
            return SpanReading(obsid, "no_archive_evt2")
        gti_start, gti_stop = utils.read_gti_span(event_file)
    except Exception as exc:
        # One damaged file (truncated gzip, two revisions, no GTI) is reported,
        # not allowed to stop a pass over thousands of them.
        return SpanReading(obsid, "unreadable", reason=f"{type(exc).__name__}: {exc}")
    return SpanReading(obsid, "filled", gti_start, gti_stop)


def _read_spans(workdir: Path, obsids: list[int], processes: int) -> list[SpanReading]:
    """read_span for every obsid, in order, logging progress."""
    progress_every = 1000

    def logged(readings):
        for count, reading in enumerate(readings, start=1):
            if count % progress_every == 0:
                logger.info(f"read {count:,} of {len(obsids):,}")
            yield reading

    if processes <= 1:
        return list(logged(read_span(workdir, obsid) for obsid in obsids))
    with ProcessPoolExecutor(max_workers=processes) as pool:
        readings = pool.map(read_span, [workdir] * len(obsids), obsids, chunksize=16)
        return list(logged(readings))


def backfill(
    db_file: Path, workdir: Path, dry_run: bool = False, processes: int = 1
) -> dict:
    """Fill astromon_obs_gti in `db_file` for the obsids in its astromon_obs.

    Returns a summary: the obsids "filled", those with "no_archive_evt2" and
    "no_work_dir", and "unreadable", a dict of obsid to the reason.
    """
    db_file = Path(db_file)
    obsids = sorted(
        {int(obsid) for obsid in db.get_table("astromon_obs", db_file)["obsid"]}
    )
    logger.info(f"astromon_obs obsids: {len(obsids):,}")
    logger.info(f"reading archive evt2 GTIs under {workdir}")

    readings = _read_spans(Path(workdir), obsids, processes)
    by_outcome = {
        outcome: [reading for reading in readings if reading.outcome == outcome]
        for outcome in OUTCOMES
    }
    for outcome in OUTCOMES:
        logger.info(f"{outcome}: {len(by_outcome[outcome]):,} obsids")
    for outcome in ("no_archive_evt2", "no_work_dir"):
        if by_outcome[outcome]:
            obsid_list = " ".join(str(reading.obsid) for reading in by_outcome[outcome])
            logger.info(f"{outcome} obsids: {obsid_list}")
    for reading in by_outcome["unreadable"]:
        logger.warning(f"unreadable: {reading.obsid}: {reading.reason}")

    filled = by_outcome["filled"]
    if dry_run:
        logger.info("Dry run -- nothing written.")
    elif filled:
        obs_gti = Table(
            {
                "obsid": [reading.obsid for reading in filled],
                "gti_start": [reading.gti_start for reading in filled],
                "gti_stop": [reading.gti_stop for reading in filled],
            },
            dtype=[
                db.ASTROMON_OBS_GTI_DTYPE[name]
                for name in db.ASTROMON_OBS_GTI_DTYPE.names
            ],
        )
        db.save("astromon_obs_gti", obs_gti, db_file, expect_existing=True)
        logger.info(f"wrote {len(obs_gti):,} astromon_obs_gti rows to {db_file}")

    return {
        "filled": [reading.obsid for reading in filled],
        "no_archive_evt2": [reading.obsid for reading in by_outcome["no_archive_evt2"]],
        "no_work_dir": [reading.obsid for reading in by_outcome["no_work_dir"]],
        "unreadable": {
            reading.obsid: reading.reason for reading in by_outcome["unreadable"]
        },
    }


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", required=True, type=Path, dest="db_file")
    parser.add_argument(
        "--workdir",
        required=True,
        type=Path,
        help="Top-level work directory, holding obsNN/<obsid>/primary/.",
    )
    parser.add_argument(
        "--processes",
        type=int,
        default=1,
        help="Worker processes reading evt2 files. Default: 1.",
    )
    parser.add_argument("--dry-run", action="store_true", dest="dry_run")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    backfill(args.db_file, args.workdir, dry_run=args.dry_run, processes=args.processes)


if __name__ == "__main__":
    main()
