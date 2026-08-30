"""Backfill astromon_status success rows for obsids already in astromon_obs.

astromon_status did not exist before this branch (see db.save_status). Every
obsid currently in astromon_obs got there by actually succeeding -- an obsid
with zero celldetect sources never gets its astromon_obs row saved at all (see
process_obsid() in get_cat_obs_data.py), so presence in astromon_obs already
means "processed successfully". This script makes that fact explicit in
astromon_status instead of leaving the table empty for every obsid that
predates this feature, using data already in the DB:

  - status: "success" for every obsid (see above for why that is safe to
    assume).
  - ascdsver: read straight from astromon_obs.ascdsver.
  - note: actual source/xcorr counts for that obsid, across all detect
    methods, so it looks like what a live run would have recorded (see
    save_with_lock's note in process_one_obsid.py) rather than an empty
    placeholder.
  - timestamp: when this backfill ran, not the original processing time,
    which was never recorded and cannot be recovered.

Purely additive: obsids that already have an astromon_status row (from a live
run since this feature landed) are left alone, so this is safe to run more
than once and cannot downgrade a real recorded outcome to a backfilled one.

Usage
-----
    python -m astromon.scripts.maintenance.backfill_status_from_obs --db PATH [--dry-run]
"""

import argparse
import logging
from collections import Counter

import numpy as np
from astropy.table import Table
from cxotime import CxoTime

from astromon import db, utils

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("backfill_status_from_obs")


def build_backfill_rows(dbfile) -> Table:
    """Return the astromon_status rows to add: one per obsid not already present."""
    obs = db.get_table("astromon_obs", dbfile)
    xray_src = db.get_table("astromon_xray_src", dbfile)
    xcorr = db.get_table("astromon_xcorr", dbfile)

    n_sources = Counter(np.asarray(xray_src["obsid"]).tolist())
    n_xcorr = Counter(np.asarray(xcorr["obsid"]).tolist())

    try:
        existing = db.get_table("astromon_status", dbfile)
        already_tracked = set(np.asarray(existing["obsid"]).tolist())
    except utils.MissingTableException:
        already_tracked = set()

    obsids = np.asarray(obs["obsid"])
    ascdsvers = np.asarray(obs["ascdsver"]).astype(str)
    to_add = [
        {
            "obsid": int(obsid),
            "status": "success",
            "note": f"{n_sources[int(obsid)]} sources, {n_xcorr[int(obsid)]} xcorr "
            "(backfilled from existing astromon_obs)",
            "ascdsver": ascdsver,
        }
        for obsid, ascdsver in zip(obsids, ascdsvers, strict=True)
        if int(obsid) not in already_tracked
    ]

    logger.info(
        f"{len(obs)} obsids in astromon_obs, {len(already_tracked)} already tracked "
        f"in astromon_status, {len(to_add)} to backfill"
    )
    if not to_add:
        return Table(dtype=db.DTYPES["astromon_status"])

    timestamp = CxoTime.now().isot
    return Table(
        {
            "obsid": [row["obsid"] for row in to_add],
            "status": [row["status"] for row in to_add],
            "note": [row["note"] for row in to_add],
            "ascdsver": [row["ascdsver"] for row in to_add],
            "timestamp": [timestamp] * len(to_add),
        },
        dtype=[np.int32, "S24", "S200", "S32", "S32"],
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--db", required=True, type=str, help="Path to the astromon HDF5 file."
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        default=False,
        dest="dry_run",
        help="Compute and show what would be added, but do not write to the DB.",
    )
    args = parser.parse_args()

    rows = build_backfill_rows(args.db)

    if len(rows) == 0:
        logger.info("Nothing to backfill.")
        return

    logger.info("Sample of rows to add:")
    for row in rows[:5]:
        logger.info(
            f"  obsid={row['obsid']} status={row['status']} "
            f"ascdsver={row['ascdsver']!r} note={row['note']!r}"
        )

    if args.dry_run:
        logger.info(f"Dry run -- would add {len(rows)} rows. Not writing.")
        return

    db.create_empty_tables(args.db, table_names=["astromon_status"])
    db.save("astromon_status", rows, args.db, expect_existing=True)
    logger.info(f"Added {len(rows)} astromon_status rows.")


if __name__ == "__main__":
    main()
