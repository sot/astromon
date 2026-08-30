"""Backfill versions_done/catalog_matched for astromon_status rows that predate them.

versions_done and catalog_matched did not exist when backfill_status_from_obs.py
first populated astromon_status, so every row it wrote got the conservative
save_status() defaults (versions_done="", catalog_matched=False) -- which is
wrong for them specifically: they all predate the --skip-catalog-match feature,
so catalog matching was always attempted for them, and versions_done can be
derived exactly the way the live pipeline derives it now, from
astromon_xray_src.detect_method.

This is a one-time schema-migration follow-up, not a general-purpose backfill:
unlike backfill_status_from_obs.py it targets a specific tell (status="success"
with an empty versions_done) rather than "no row at all", so it is safe to run
more than once -- a row it has already updated has a non-empty versions_done
and is skipped the second time -- but it is not the tool for anything broader
than this one migration.

Usage
-----
    python -m astromon.scripts.maintenance.backfill_status_versions_done --db PATH [--dry-run]
"""

import argparse
import logging
from collections import defaultdict

import numpy as np
from astropy.table import Table

from astromon import db

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("backfill_status_versions_done")


def build_updated_rows(dbfile) -> Table:
    """Return the astromon_status rows needing versions_done/catalog_matched filled in.

    Only status="success" rows with an empty versions_done qualify -- that is
    the tell that a row predates this migration, per the module docstring.
    Everything else about the row (status, note, ascdsver, timestamp) is
    carried over unchanged.
    """
    status = db.get_table("astromon_status", dbfile)
    xray_src = db.get_table("astromon_xray_src", dbfile)

    is_success = np.asarray(status["status"]).astype(str) == "success"
    is_unmigrated = np.asarray(status["versions_done"]).astype(str) == ""
    to_update = status[is_success & is_unmigrated]

    logger.info(
        f"{len(status)} astromon_status rows, {is_success.sum()} status=success, "
        f"{len(to_update)} still need versions_done backfilled"
    )
    if len(to_update) == 0:
        return Table(dtype=db.DTYPES["astromon_status"])

    methods_by_obsid = defaultdict(set)
    for obsid, method in zip(
        np.asarray(xray_src["obsid"]),
        np.asarray(xray_src["detect_method"]).astype(str),
        strict=True,
    ):
        methods_by_obsid[int(obsid)].add(method)

    to_update["versions_done"] = [
        ",".join(sorted(methods_by_obsid.get(int(obsid), ())))
        for obsid in to_update["obsid"]
    ]
    to_update["catalog_matched"] = 1

    n_no_xray = sum(
        1 for obsid in to_update["obsid"] if not methods_by_obsid.get(int(obsid))
    )
    if n_no_xray:
        logger.warning(
            f"{n_no_xray} obsid(s) marked status=success have no astromon_xray_src "
            "rows at all -- versions_done left empty for them, which will look "
            "like they still need this migration on a rerun. Investigate before "
            "relying on this script being idempotent for those obsids."
        )

    return to_update


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
        help="Compute and show what would be updated, but do not write to the DB.",
    )
    args = parser.parse_args()

    rows = build_updated_rows(args.db)

    if len(rows) == 0:
        logger.info("Nothing to backfill.")
        return

    logger.info("Sample of rows to update:")
    for row in rows[:5]:
        logger.info(
            f"  obsid={row['obsid']} versions_done={row['versions_done']!r} "
            f"catalog_matched={row['catalog_matched']}"
        )

    if args.dry_run:
        logger.info(f"Dry run -- would update {len(rows)} rows. Not writing.")
        return

    db.save("astromon_status", rows, args.db, expect_existing=True)
    logger.info(f"Updated {len(rows)} astromon_status rows.")


if __name__ == "__main__":
    main()
