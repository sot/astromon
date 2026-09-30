"""Record the applied alignment of existing rows from a directory of acal1 files.

Why this is needed
------------------
Rows written before astromon_xray_src had acal_dy/acal_dz and
caldb_version_source read them as NaN and "", so every rebase reconstructs
their alignment matrix from caldb_version -- and that version came from the
event file for anything downloaded from the public archive, which can have been
reprocessed after the aspect run that applied the matrix. For the obsids whose
acal1 files are at hand (e.g. pulled once from the CXC archive), this records
what Observation.get_calalign now records at ingest. Per obsid, from its acal1
files (utils.read_acal_files):

- files that agree on CALDBVER: caldb_version becomes that version, and
  caldb_version_source "acal1";
- files that also agree on the applied matrix: acal_dy/acal_dz become its
  offsets. OBIs observed across an alignment change have no single matrix, so
  those rows keep NaN and go on being reconstructed;
- files that disagree on CALDBVER, and obsids with no acal1 file, are left as
  they are.

Only rows the acal1 files can describe are touched: an acal1 file newer than a
row's recorded caldb_version comes from a later reprocessing than the one that
made the stored positions (obsid 62649: stored from CALDB 4.9.3, acal1 from
4.12.6, 0.55" apart in dy), and a row with no version ("0.0") cannot be
checked. Those rows are left as they are. An acal1 file as old as, or older
than, the recorded version is the aspect run the events used -- the events can
have been reprocessed since without redoing it.

Every recorded caldb_version is also normalized (utils.normalized_caldb_version:
no trailing "."), as get_calalign now records them, since get_calalign_offsets
deliberately refuses a string like "4.9.6.".

The acal1 files are matched to obsids by their OBS_ID header, so the directory
can be flat. The whole table is read and written back: run this while nothing
else writes to the database.

Usage
-----
::

    python -m astromon.scripts.maintenance.backfill_acal_offsets \\
        --db /Volumes/Black/data/astromon/astromon.h5 --acal-dir ~/acal --dry-run
"""

import argparse
import logging
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits

from astromon import db, utils

logger = logging.getLogger("backfill_acal_offsets")


def index_acal_dir(acal_dir):
    """Each obsid's acal1 files in `acal_dir`, by the OBS_ID in their headers."""
    files = defaultdict(list)
    for path in sorted(Path(acal_dir).glob("*acal1.fits*")):
        files[int(fits.getval(path, "OBS_ID", ext=1))].append(path)
    return dict(files)


def _row_groups(xray):
    """Row indices of `xray` per (obsid, recorded caldb_version), in one pass."""
    obsids = np.asarray(xray["obsid"])
    versions = np.array([str(v) for v in xray["caldb_version"]])
    order = np.lexsort((versions, obsids))
    keys = list(zip(obsids[order].tolist(), versions[order].tolist(), strict=True))
    groups = defaultdict(list)
    for key, row in zip(keys, order.tolist(), strict=True):
        groups[key].append(row)
    return {key: np.array(rows) for key, rows in groups.items()}


def backfill(db_file, acal_dir, dry_run=False):
    """Fill astromon_xray_src's acal columns in `db_file` from `acal_dir`.

    Returns a summary counting (obsid, recorded caldb_version) row groups: how
    many got a matrix ("matrix recorded"), only a version ("acal1 files
    disagree on the matrix"), or nothing ("acal1 newer than the recorded
    processing", "acal1 files disagree on CALDBVER", "no acal1 file"); and
    "version changes", the number of groups per (old caldb_version, new
    caldb_version); and "versions normalized", the number of rows whose
    recorded version was malformed.
    """
    db_file = Path(db_file)
    xray = db.get_table("astromon_xray_src", db_file)
    acal_files = index_acal_dir(acal_dir)
    logger.info(f"astromon_xray_src rows: {len(xray):,}")
    logger.info(f"acal1 files for {len(acal_files):,} obsids in {acal_dir}")

    summary = Counter()
    recorded_versions = [str(v) for v in xray["caldb_version"]]
    normalized = [utils.normalized_caldb_version(v) for v in recorded_versions]
    summary["versions normalized"] = sum(
        old != new for old, new in zip(recorded_versions, normalized, strict=True)
    )
    xray["caldb_version"] = normalized
    groups = _row_groups(xray)
    version_changes = Counter()
    for (obsid, recorded), rows in groups.items():
        paths = acal_files.get(obsid)
        if not paths:
            summary["no acal1 file"] += 1
            continue
        acal = utils.read_acal_files(paths)
        if acal["caldb_version"] is None:
            summary["acal1 files disagree on CALDBVER"] += 1
            continue
        if recorded == "0.0" or utils.caldb_version_order(
            acal["caldb_version"]
        ) > utils.caldb_version_order(recorded):
            summary["acal1 newer than the recorded processing"] += 1
            continue
        if recorded != acal["caldb_version"]:
            version_changes[(recorded, acal["caldb_version"])] += 1
        xray["caldb_version"][rows] = acal["caldb_version"]
        xray["caldb_version_source"][rows] = "acal1"
        if acal["aca_misalign"] is None:
            summary["acal1 files disagree on the matrix"] += 1
            continue
        dys, dzs = utils.get_offsets(np.array([acal["aca_misalign"]]))
        xray["acal_dy"][rows] = dys[0]
        xray["acal_dz"][rows] = dzs[0]
        summary["matrix recorded"] += 1

    for outcome in (
        "matrix recorded",
        "acal1 files disagree on the matrix",
        "acal1 newer than the recorded processing",
        "acal1 files disagree on CALDBVER",
        "no acal1 file",
    ):
        logger.info(f"{outcome}: {summary[outcome]:,} (obsid, version) groups")
    for (old_version, new_version), count in sorted(version_changes.items()):
        logger.info(f"caldb_version {old_version} -> {new_version}: {count:,} groups")
    logger.info(
        f"malformed caldb_version normalized: {summary['versions normalized']:,} rows"
    )

    if dry_run:
        logger.info("Dry run -- nothing written.")
    else:
        db.save(
            "astromon_xray_src",
            xray,
            db_file,
            ignore_obsid=True,
            expect_existing=True,
        )
        logger.info("written")

    return {"db_file": db_file, **summary, "version changes": dict(version_changes)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", required=True, type=Path, dest="db_file")
    parser.add_argument("--acal-dir", required=True, type=Path, dest="acal_dir")
    parser.add_argument("--dry-run", action="store_true", dest="dry_run")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    backfill(args.db_file, args.acal_dir, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
