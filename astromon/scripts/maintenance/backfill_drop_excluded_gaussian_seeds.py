"""Remove gaussian_detect rows whose celldetect seed is now excluded pre-fit.

Why this is needed
------------------
Every gaussian_detect row in the database today was fit unconditionally, even
for celldetect seeds that current selection logic now excludes before a fit
is ever attempted: a seed within near_neighbor_dist of another source
(crowded), on a grating dispersion arm, or on an ACIS readout streak (unless
it is the observation's single brightest source, which is exempted from the
streak drop the same way the fit-time filter exempts it). Rows already
written have to be removed to bring the database into line with
observation._drop_crowded_seeds, _drop_grating_arm_seeds and
_drop_acis_streak_seeds, which now stop those fits from being attempted at
all. The archived .src files these rows came from reflect whatever version
of that selection logic was live when each obsid was originally processed,
not necessarily today's.

This is pure removal, not detection: celldetect's own near_neighbor_dist,
grating_arm, acis_streak and brightest columns are already stored, untouched
by anything here, so no fit has to be re-run and no CIAO or event data is
needed -- only the same comparisons the fit-time filters make, applied to
rows that predate them. This is deliberately cheaper and safer than
reapplying the fit: a direct refit of all ~10,296 gaussian_detect obsids was
verified (read-only, against a scratch copy of the outputs) to produce
exactly the same row-level drop decisions as this column-based approach,
because nothing else about the fit changed between the archived rows'
vintage and now -- only this pre-fit seed selection.

Follow this with rebuild_xcorr, since the removed sources may carry matches.

Usage
-----
::

    python -m astromon.scripts.maintenance.backfill_drop_excluded_gaussian_seeds \\
        --db /Volumes/Black/data/astromon/astromon.h5 --dry-run
"""

import argparse
import logging
from pathlib import Path

import numpy as np

from astromon import db, utils

logger = logging.getLogger("backfill_drop_excluded_gaussian_seeds")


def excluded_seed_gaussian_mask(xray) -> dict[str, np.ndarray]:
    """Per-reason boolean masks for gaussian_detect rows whose seed is now excluded.

    Returns a dict with keys "crowded", "grating_arm", "acis_streak", each a
    boolean array the length of `xray`, true where that specific reason
    applies to a gaussian_detect row (based on its originating celldetect
    seed's own stored columns). A row matching more than one reason appears
    true in more than one array -- the caller decides how to report and drop
    the union.

    A gaussian_detect row with no matching celldetect seed (should not
    happen; the seed that produced it is the same row that would carry these
    columns) is left alone in every mask rather than guessed at.
    """
    method = np.asarray(xray["detect_method"]).astype(str)
    obsid = np.asarray(xray["obsid"])
    source_id = np.asarray(xray["id"])
    nnd = np.asarray(xray["near_neighbor_dist"]).astype(float)
    grating_arm = np.asarray(xray["grating_arm"]).astype(bool)
    acis_streak = np.asarray(xray["acis_streak"]).astype(bool)
    brightest = np.asarray(xray["brightest"]).astype(bool)

    cd = method == "celldetect"
    cd_seed = {
        (int(o), int(i)): (float(n), bool(g), bool(s), bool(b))
        for o, i, n, g, s, b in zip(
            obsid[cd],
            source_id[cd],
            nnd[cd],
            grating_arm[cd],
            acis_streak[cd],
            brightest[cd],
            strict=True,
        )
    }

    masks = {
        reason: np.zeros(len(xray), dtype=bool)
        for reason in ("crowded", "grating_arm", "acis_streak")
    }
    for k in np.where(method == "gaussian_detect")[0]:
        seed = cd_seed.get((int(obsid[k]), int(source_id[k])))
        if seed is None:
            continue
        seed_nnd, seed_grating_arm, seed_acis_streak, seed_brightest = seed
        if seed_nnd <= utils.NEAR_NEIGHBOR_DIST_ARCSEC:
            masks["crowded"][k] = True
        if seed_grating_arm:
            masks["grating_arm"][k] = True
        if seed_acis_streak and not seed_brightest:
            masks["acis_streak"][k] = True
    return masks


def backfill(db_file: Path, dry_run: bool = False) -> dict:
    """Drop excluded-seed gaussian_detect rows from `db_file`."""
    db_file = Path(db_file)
    xray = db.get_table("astromon_xray_src", db_file)
    masks = excluded_seed_gaussian_mask(xray)
    drop = masks["crowded"] | masks["grating_arm"] | masks["acis_streak"]
    n_dropped = int(drop.sum())

    logger.info(f"astromon_xray_src rows: {len(xray):,}")
    logger.info(
        f"gaussian_detect rows seeded from an excluded celldetect source: "
        f"{n_dropped:,} (crowded: {int(masks['crowded'].sum()):,}, "
        f"grating_arm: {int(masks['grating_arm'].sum()):,}, "
        f"acis_streak: {int(masks['acis_streak'].sum()):,})"
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

    return {
        "db_file": db_file,
        "dropped": n_dropped,
        "dropped_crowded": int(masks["crowded"].sum()),
        "dropped_grating_arm": int(masks["grating_arm"].sum()),
        "dropped_acis_streak": int(masks["acis_streak"].sum()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", required=True, type=Path, dest="db_file")
    parser.add_argument("--dry-run", action="store_true", dest="dry_run")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(message)s")
    backfill(args.db_file, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
