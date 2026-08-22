"""Clear the download sentinel for obsids with an INVALID make_images cache.

Background
----------
``invalidate_streak_caches.py`` set ``task_make_images.json`` to
``return_code: INVALID`` for ~437 ACIS obsids so the pipeline would re-run
``make_images`` with the corrected ``ssigma=3`` parameter.

The re-run needs to download a fresh raw evt2 file.  But ``observation.py``'s
download method skips the download entirely when the ``secondary/`` subdirectory
already exists (the "already downloaded" sentinel).  The result: the pipeline
tries to run ``make_images``, can't find ``primary/*_evt2.fits*`` (only the
filtered version is there), and fails with::

    Missing input files for task make_images: events: primary/*_evt2.fits*

Fix: delete the ``secondary/`` sentinel so the next pipeline run re-downloads
the raw data and re-runs ``make_images`` properly.

Usage::

    # Dry run — show what would be cleared, write obsid list:
    python -m astromon.scripts.maintenance.clear_invalid_make_images --workdir /Volumes/Black/data/astromon/work

    # Actually clear:
    python -m astromon.scripts.maintenance.clear_invalid_make_images --workdir /Volumes/Black/data/astromon/work --apply

    # Write cleared obsids to a file for use as --obsid-list:
    python -m astromon.scripts.maintenance.clear_invalid_make_images --workdir /Volumes/Black/data/astromon/work --apply \\
        --obsid-list /Volumes/Black/data/astromon/obsids_reimage.txt
"""

import argparse
import json
import logging
import shutil
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def find_invalid_make_images(workdir: Path) -> list[tuple[int, Path]]:
    """Return (obsid, obsid_workdir) pairs where task_make_images.json is INVALID.

    Parameters
    ----------
    workdir
        Root of the pipeline working directory (e.g. ``/Volumes/Black/data/astromon/work``).
    """
    results = []
    for cache_json in sorted(workdir.rglob("task_make_images.json")):
        try:
            data = json.loads(cache_json.read_text())
        except (json.JSONDecodeError, OSError) as exc:
            logger.warning(f"Could not read {cache_json}: {exc}")
            continue
        if data.get("return_code") != "INVALID":
            continue
        # Expected layout: {workdir}/obs{NN}/{obsid}/cache/task_make_images.json
        obsid_dir = cache_json.parent.parent
        try:
            obsid = int(obsid_dir.name)
        except ValueError:
            logger.warning(f"Unexpected path layout at {cache_json}, skipping")
            continue
        results.append((obsid, obsid_dir))
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--workdir",
        required=True,
        type=Path,
        help="Root of the pipeline working directory.",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        default=False,
        help="Actually delete secondary/ directories.  Without this flag the "
        "script only reports what it would do (dry run).",
    )
    parser.add_argument(
        "--obsid-list",
        default=None,
        type=Path,
        dest="obsid_list",
        help="Write the list of affected obsids to this file (one per line, "
        "descending order) for use as --obsid-list in run_all.py.",
    )
    args = parser.parse_args()

    mode = "APPLY" if args.apply else "DRY RUN"
    logger.info(f"Mode: {mode}")
    logger.info(f"Workdir: {args.workdir}")

    affected = find_invalid_make_images(args.workdir)
    if not affected:
        logger.info("No obsids with INVALID make_images cache found.")
        return

    logger.info(f"Found {len(affected)} obsid(s) with INVALID make_images cache.")

    cleared = []
    skipped_no_secondary = []

    for obsid, obsid_dir in affected:
        secondary = obsid_dir / "secondary"
        if not secondary.exists():
            logger.debug(f"  {obsid}: no secondary/ to clear, skipping")
            skipped_no_secondary.append(obsid)
            continue

        if args.apply:
            logger.info(f"  {obsid}: deleting {secondary}")
            shutil.rmtree(secondary)
            cleared.append(obsid)
        else:
            logger.info(f"  {obsid}: would delete {secondary}")
            cleared.append(obsid)

    verb = "Cleared" if args.apply else "Would clear"
    logger.info(
        f"\n{verb} {len(cleared)} secondary/ directories. "
        f"{len(skipped_no_secondary)} obsid(s) had no secondary/ (nothing to do)."
    )

    if args.obsid_list is not None and cleared:
        # Sort descending (recent first) to match the backfill convention.
        obsids_sorted = sorted(cleared, reverse=True)
        args.obsid_list.write_text("\n".join(str(o) for o in obsids_sorted) + "\n")
        logger.info(f"Wrote {len(obsids_sorted)} obsids to {args.obsid_list}")
        if not args.apply:
            logger.info(
                "(This list reflects what would be cleared; re-run with --apply "
                "to actually clear secondary/ and then use the list.)"
            )


if __name__ == "__main__":
    main()
