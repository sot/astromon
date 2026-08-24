"""Delete task_make_images.json for ACIS obsids that have stale streak files.

Background
----------
acis_streak_map runs as a sub-step inside the make_images task.  Commit 9f0e441
added ssigma=3 to that call; observations processed before that commit have streak
detection results computed with the old (looser) default.

The 1,785 ACIS obsids in preserve-workdir that have an acis_streaks.fits file
are the only ones affected — obsids with no streak file either have no streaks
regardless of ssigma, or are HRC (no streaks at all).  Deleting
cache/task_make_images.json for those obsids causes process_one_obsid.py to
re-run make_images (and therefore acis_streak_map with ssigma=3) during the
subsequent full pipeline rerun.

Usage::

    python -m astromon.scripts.maintenance.invalidate_streak_caches \\
        [--preserve-workdir DIR] [--dry-run]
"""

import argparse
import json
import logging
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
logger = logging.getLogger("invalidate_streak_caches")

_DEFAULT_PRESERVE_WORKDIR = Path("/Volumes/Black/data/astromon/work")
_TASK_CACHE_NAME = "task_make_images.json"


def find_streak_obsids(preserve_workdir: Path) -> list[tuple[int, Path]]:
    """Return (obsid, streak_fits_path) for all obsids with a streak file.

    Parameters
    ----------
    preserve_workdir
        Root of the preserve-workdir tree, e.g. /Volumes/Black/data/astromon/work.

    Returns
    -------
    List of (obsid, streak_file_path) sorted by obsid.
    """
    streak_files = sorted(preserve_workdir.glob("obs*/*/**/*_acis_streaks.fits"))
    results = []
    for path in streak_files:
        # Path pattern: .../obs{NN}/{obsid}/images/{obsid}_acis_streaks.fits
        try:
            obsid = int(path.parent.parent.name)
        except ValueError:
            logger.warning(f"Could not parse obsid from {path}, skipping")
            continue
        results.append((obsid, path))
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--preserve-workdir",
        type=Path,
        default=_DEFAULT_PRESERVE_WORKDIR,
        help=f"Root of preserve-workdir (default: {_DEFAULT_PRESERVE_WORKDIR})",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be deleted without actually deleting anything.",
    )
    args = parser.parse_args()

    preserve_workdir = args.preserve_workdir
    if not preserve_workdir.exists():
        logger.error(f"Preserve-workdir {preserve_workdir} does not exist")
        return

    logger.info(f"Scanning {preserve_workdir} for ACIS streak files …")
    streak_obsids = find_streak_obsids(preserve_workdir)
    logger.info(f"Found {len(streak_obsids)} ACIS obsids with streak files")

    n_deleted = 0
    n_missing = 0
    n_already_invalid = 0

    for obsid, streak_path in streak_obsids:
        obsid_dir = streak_path.parent.parent  # .../obs{NN}/{obsid}/
        cache_file = obsid_dir / "cache" / _TASK_CACHE_NAME

        if not cache_file.exists():
            n_missing += 1
            logger.debug(f"OBSID={obsid}: no cache file (already absent)")
            continue

        # Check if it's already invalid so we don't double-count.
        try:
            with open(cache_file) as fh:
                data = json.load(fh)
            if data.get("return_code") == "INVALID":
                n_already_invalid += 1
                logger.debug(f"OBSID={obsid}: cache already INVALID, skipping")
                continue
        except (json.JSONDecodeError, KeyError):
            pass  # Corrupted — delete it anyway.

        if args.dry_run:
            logger.info(f"  (dry-run) would delete {cache_file}")
        else:
            cache_file.unlink()
            logger.debug(f"OBSID={obsid}: deleted {cache_file.name}")
        n_deleted += 1

    logger.info(
        f"{'(dry-run) would delete' if args.dry_run else 'Deleted'} "
        f"{n_deleted} make_images cache files; "
        f"{n_missing} already absent; "
        f"{n_already_invalid} already INVALID"
    )
    logger.info(
        "Next step: run a full pipeline rerun over all obsids — those 1,785 "
        "obsids will re-run make_images with ssigma=3; all others only re-run "
        "the catalog/cross-match stage (fast)."
    )


if __name__ == "__main__":
    main()
