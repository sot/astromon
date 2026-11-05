"""Recompute celldetect grating_arm flags from the .src positions and the arm mask.

Why this is needed
------------------
``Observation._on_grating_arm`` flags sources inside a grating observation's
dispersed arms using ``images/<obsid>_grating_arm_mask.fits``, and flags nothing
when that file is missing. For 13 grating observations it was never made: the
zero-order finder crashed on events with NaN sky positions and make_images swallowed
the error. Their celldetect rows were stored with grating_arm False everywhere, so
arm knots entered cross-matching and gaussian_detect fit seeds it should drop.

What it does
------------
For each obsid with a mask and a celldetect ``.src`` (under ``--root``, work first,
then archive), recomputes which sources are on an arm from their sky X/Y in the
``.src`` and rewrites ``grating_arm`` for that obsid's stored celldetect rows. No
other column changes and nothing is re-detected, so neither CIAO nor a download is
needed. ``--invalidate-caches`` also marks the obsid's cached celldetect and
gaussian_detect source tables invalid, since they still hold the old flags.

What it does not do
-------------------
- Make masks: put them in place first (make_images makes them).
- Drop the gaussian_detect rows now seeded from an arm source: run
  ``backfill_drop_excluded_gaussian_seeds.py`` next, then ``rebuild_xcorr.py`` for
  any obsid with stored matches.

Usage
-----
::

    python -m astromon.scripts.maintenance.backfill_grating_arm_flags \\
        --db astromon.h5 --obsid-file obsids.txt \\
        --root /Volumes/Black/data/astromon/work \\
        --root /Volumes/Black/data/astromon/archive [--invalidate-caches] [--dry-run]
"""

import argparse
import logging
from collections.abc import Sequence
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack

from astromon import db
from astromon.observation import Observation

logger = logging.getLogger("backfill_grating_arm_flags")

#: The source tables whose cached copies carry grating_arm (directly, or through
#: the seed rule applied to gaussian_detect).
CACHED_VERSIONS = ("celldetect", "gaussian_detect")


def _observation(obsid: int, roots: Sequence[Path]) -> Observation:
    """An Observation reading from roots[0], falling back to roots[1] (no CIAO)."""
    roots = [Path(root) for root in roots]
    archive = roots[1] if len(roots) > 1 else roots[0]
    return Observation(obsid, workdir=roots[0], archive_dir=archive, use_ciao=False)


def arm_source_ids(obsid: int, roots: Sequence[Path]) -> set[int] | None:
    """COMPONENT ids of the celldetect sources on a grating arm, or None if no mask.

    Example::

        work = Path("/Volumes/Black/data/astromon/work")
        on_arm = arm_source_ids(19864, [work])
    """
    obs = _observation(obsid, roots)
    mask = obs.file_path(f"images/{obsid}_grating_arm_mask.fits")
    src_file = obs.file_path(f"sources/{obsid}_celldetect.src")
    if not (mask.exists() and src_file.exists()):
        return None
    sources = Table(fits.getdata(src_file, 1))
    on_arm = obs._on_grating_arm(sources)
    return {int(component) for component in sources["COMPONENT"][on_arm]}


def backfill(
    dbfile: Path,
    obsids: Sequence[int],
    roots: Sequence[Path],
    invalidate_caches: bool = False,
    dry_run: bool = False,
) -> dict:
    """Rewrite grating_arm for the celldetect rows of `obsids` in `dbfile`.

    Returns ``{"updated": {obsid: rows whose flag changed}, "no_mask": [obsids]}``.
    """
    xray = db.get_table("astromon_xray_src", dbfile)
    is_celldetect = np.asarray(xray["detect_method"]).astype(str) == "celldetect"
    updated = {}
    no_mask = []
    corrected = []
    for obsid in obsids:
        arm_ids = arm_source_ids(obsid, roots)
        if arm_ids is None:
            no_mask.append(int(obsid))
            continue
        rows = xray[is_celldetect & (np.asarray(xray["obsid"]) == obsid)]
        flags = np.isin(np.asarray(rows["id"]), list(arm_ids))
        n_changed = int(np.sum(flags != np.asarray(rows["grating_arm"]).astype(bool)))
        if n_changed == 0:
            continue
        rows["grating_arm"] = flags
        corrected.append(rows)
        updated[int(obsid)] = n_changed
        logger.info(
            f"OBSID={obsid}: {int(flags.sum())} of {len(rows)} on an arm,"
            f" {n_changed} changed"
        )

    if corrected and not dry_run:
        db.save("astromon_xray_src", vstack(corrected), dbfile, expect_existing=True)
        if invalidate_caches:
            for obsid in updated:
                obs = _observation(obsid, roots)
                for version in CACHED_VERSIONS:
                    obs._cached_sources.invalidate_result(version=version)
    return {"updated": updated, "no_mask": no_mask}


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--db", required=True, type=Path, dest="db_file")
    parser.add_argument("--obsid-file", required=True, type=Path)
    parser.add_argument(
        "--root",
        required=True,
        type=Path,
        action="append",
        dest="roots",
        help="work dir tree (obsNN/<obsid>); a second --root is the archive fallback",
    )
    parser.add_argument("--invalidate-caches", action="store_true")
    parser.add_argument("--dry-run", action="store_true", dest="dry_run")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    if not args.db_file.exists():
        parser.error(f"db does not exist: {args.db_file}")

    obsids = [int(x) for x in args.obsid_file.read_text().split()]
    summary = backfill(
        args.db_file, obsids, args.roots, args.invalidate_caches, args.dry_run
    )
    logger.info(f"rows changed per obsid: {summary['updated']}")
    logger.info(f"obsids with no mask or no celldetect .src: {summary['no_mask']}")
    if args.dry_run:
        logger.info("dry run -- nothing written")
    else:
        logger.info(
            "next: backfill_drop_excluded_gaussian_seeds, then rebuild_xcorr for"
            " obsids with stored matches"
        )


if __name__ == "__main__":
    main()
