"""Recompute astromon_cat_src.y_angle/z_angle from the evt2-header pointing.

Why this is needed
-------------------
requery_cat_src.py and backfill_icrf2.py used to derive their pointing
quaternion from mica's local obspar archive / astromon_obs's nominal target
ra/dec/roll. Both are documented elsewhere in this codebase as unreliable --
obspar can disagree with the evt2 header's measured pointing by tens of
arcsec for the same obsid (see Observation.get_evt2_info's docstring), and
astromon_obs's ra/dec/roll are nominal target values, not the pointing. That
has since been fixed (both scripts now use evt2_pointing, the same source
astromon_xray_src's own y_angle/z_angle are computed from) -- but the fix only
stops it happening again. Rows already written with the old, mismatched
pointing are still wrong: dy = x_y_angle - c_y_angle stops canceling a shared
pointing error, and picks up the mismatch as if it were a real astrometric
offset, the moment x_y_angle (already evt2-based) and c_y_angle (still
obspar/nominal-based) disagree.

This repairs those rows directly: y_angle/z_angle recomputed from each row's
own (already-correct) ra/dec and the now-correct evt2 quaternion. Nothing else
changes -- id/obsid/catalog/ra/dec/name/mag/separation/celldetect_x_id are
untouched, so astromon_xcorr's c_id references are unaffected by this write
(no rebuild_xcorr needed for this step alone). An obsid with no evt2 file
available locally is skipped and reported, not guessed at with a worse
pointing source.

Usage
-----
::

    python -m astromon.scripts.maintenance.fix_cat_src_pointing \\
        --db astromon.h5 --workdir WORKDIR --archive-dir ARCHIVE_DIR \\
        --obsid-file touched_obsids.txt [--dry-run]
"""

import argparse
import logging
from pathlib import Path
from typing import Sequence

import numpy as np
from astropy.table import Table
from chandra_aca.transform import radec_to_yagzag

from astromon import db
from astromon.scripts.maintenance.requery_cat_src import evt2_pointing

logger = logging.getLogger("fix_cat_src_pointing")


def recompute_cat_src_pointing(
    cat_src: Table, obsids: Sequence[int], workdir, archive_dir
) -> tuple[Table, list[int], int]:
    """Recompute y_angle/z_angle for every `cat_src` row belonging to `obsids`.

    Returns (subset of `cat_src` for just these obsids with y_angle/z_angle
    corrected, obsids skipped for lack of a local evt2 file, count of rows
    actually updated). The returned subset is meant to go straight to
    db.save's default merge-by-obsid-key write -- it deliberately excludes
    every other obsid's rows.
    """
    obsid_col = np.asarray(cat_src["obsid"])
    subset = cat_src[np.isin(obsid_col, list(obsids))].copy()
    subset_obsid = np.asarray(subset["obsid"])

    skipped = []
    n_updated = 0
    for obsid in obsids:
        pointing = evt2_pointing(obsid, workdir, archive_dir)
        if pointing is None:
            skipped.append(int(obsid))
            continue
        quat, _ = pointing
        mask = subset_obsid == obsid
        if not mask.any():
            continue
        new_y, new_z = radec_to_yagzag(
            np.asarray(subset["ra"])[mask], np.asarray(subset["dec"])[mask], quat
        )
        subset["y_angle"][mask] = new_y
        subset["z_angle"][mask] = new_z
        n_updated += int(mask.sum())
    return subset, skipped, n_updated


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", required=True, type=Path, dest="db_file")
    parser.add_argument("--workdir", required=True, type=Path)
    parser.add_argument("--archive-dir", required=True, type=Path, dest="archive_dir")
    parser.add_argument(
        "--obsid-file",
        required=True,
        type=Path,
        help="file with one obsid per line",
    )
    parser.add_argument("--dry-run", action="store_true", dest="dry_run")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    obsids = [int(line) for line in args.obsid_file.read_text().split() if line.strip()]
    logger.info(f"recomputing astromon_cat_src pointing for {len(obsids):,} obsid(s)")

    cat_src = db.get_table("astromon_cat_src", args.db_file)
    updated, skipped, n_updated = recompute_cat_src_pointing(
        cat_src, obsids, args.workdir, args.archive_dir
    )

    logger.info(f"rows recomputed: {n_updated:,}")
    logger.info(f"obsids skipped (no local evt2 file): {len(skipped):,}")
    if skipped:
        logger.info(f"skipped obsids: {skipped}")

    if args.dry_run:
        logger.info("dry run -- nothing written")
        return

    db.save("astromon_cat_src", updated, dbfile=args.db_file)
    logger.info("written")


if __name__ == "__main__":
    main()
