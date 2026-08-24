"""Remove gaussian_detect rows whose celldetect seed was already crowded.

Why this is needed
------------------
Every gaussian_detect row in the database today was fit unconditionally, even
when its celldetect seed sat within near_neighbor_dist of another source --
the crowding that every current selection already excludes. Rows already
written have to be removed to bring the database into line with the fix in
observation._drop_crowded_seeds, which now stops those fits from being
attempted at all.

This is pure removal, not detection: celldetect's own near_neighbor_dist is
already stored, untouched by anything here, so no fit has to be re-run and no
CIAO or event data is needed -- only the same comparison the fit-time filter
makes, applied to rows that predate it.

Follow this with rebuild_xcorr, since the removed sources may carry matches.

Usage
-----
::

    python -m astromon.scripts.maintenance.backfill_drop_crowded_gaussian \\
        --db /Volumes/Black/data/astromon/astromon.h5 --dry-run
"""

import argparse
import logging
from pathlib import Path

import numpy as np

from astromon import db, utils

logger = logging.getLogger("backfill_drop_crowded_gaussian")


def crowded_gaussian_mask(xray):
    """True for every gaussian_detect row whose celldetect seed is crowded.

    A row with no matching celldetect seed (should not happen; the seed that
    produced it is the same row that would carry its near_neighbor_dist) is
    left alone rather than guessed at.
    """
    method = np.asarray(xray["detect_method"]).astype(str)
    obsid = np.asarray(xray["obsid"])
    source_id = np.asarray(xray["id"])
    nnd = np.asarray(xray["near_neighbor_dist"]).astype(float)

    cd = method == "celldetect"
    cd_nnd = {
        (int(o), int(i)): n
        for o, i, n in zip(obsid[cd], source_id[cd], nnd[cd], strict=True)
    }

    mask = np.zeros(len(xray), dtype=bool)
    for k in np.where(method == "gaussian_detect")[0]:
        seed_nnd = cd_nnd.get((int(obsid[k]), int(source_id[k])))
        if seed_nnd is not None and seed_nnd <= utils.NEAR_NEIGHBOR_DIST_ARCSEC:
            mask[k] = True
    return mask


def backfill(db_file: Path, dry_run: bool = False) -> dict:
    """Drop crowded-seed gaussian_detect rows from `db_file`."""
    db_file = Path(db_file)
    xray = db.get_table("astromon_xray_src", db_file)
    drop = crowded_gaussian_mask(xray)
    n_dropped = int(drop.sum())

    logger.info(f"astromon_xray_src rows: {len(xray):,}")
    logger.info(
        f"gaussian_detect rows seeded from a crowded celldetect source: {n_dropped:,}"
    )

    if dry_run:
        logger.info("Dry run -- nothing written.")
    else:
        db.save(
            "astromon_xray_src",
            xray[~drop],
            db_file,
            ignore_obsid=True,
            expect_existing=True,
        )
        logger.info(
            "written; run rebuild_xcorr next, since the removed rows may carry matches"
        )

    return {"db_file": db_file, "dropped": n_dropped}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", required=True, type=Path, dest="db_file")
    parser.add_argument("--dry-run", action="store_true", dest="dry_run")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    backfill(args.db_file, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
