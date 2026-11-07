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

brightest and near_neighbor_dist describe a source against the rest of its
observation's gaussian_detect sources, so they are recomputed over the rows left,
as a fit without the dropped seeds has them. Earlier runs dropped rows without
that, so the recompute covers every obsid, whether or not this run drops anything.

Then run drop_retired_rows: the removed rows may carry astromon_21 matches, stored
before astromon_21 took celldetect sources only.

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

from astromon import db
from astromon.observation import excluded_seed_mask, with_recomputed_set_columns

logger = logging.getLogger("backfill_drop_excluded_gaussian_seeds")

# A recomputed near_neighbor_dist this close to the stored one is the same value: it
# is stored as float32, so recomputing it from the stored angles moves the last bits.
_SAME_DISTANCE_ARCSEC = 1e-3


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

    # The rule itself is shared with Observation.get_sources, which applies it
    # when reading gaussian_detect sources (observation.excluded_seed_mask).
    cd = method == "celldetect"
    cd_reasons = excluded_seed_mask(xray[cd])
    reason_names = ("crowded", "grating_arm", "acis_streak")
    cd_seed = {
        (int(o), int(i)): tuple(bool(cd_reasons[name][k]) for name in reason_names)
        for k, (o, i) in enumerate(zip(obsid[cd], source_id[cd], strict=True))
    }

    masks = {reason: np.zeros(len(xray), dtype=bool) for reason in reason_names}
    for k in np.where(method == "gaussian_detect")[0]:
        seed = cd_seed.get((int(obsid[k]), int(source_id[k])))
        if seed is None:
            continue
        for name, applies in zip(reason_names, seed, strict=True):
            masks[name][k] = applies
    return masks


def recomputed_gaussian_set_columns(xray) -> tuple[np.ndarray, np.ndarray]:
    """brightest and near_neighbor_dist for each row of `xray`.

    They are recomputed over each obsid's gaussian_detect rows; other methods' rows
    keep their stored values.
    """
    brightest = np.asarray(xray["brightest"]).copy()
    near_neighbor_dist = np.asarray(xray["near_neighbor_dist"]).copy()
    gaussian = np.flatnonzero(
        np.asarray(xray["detect_method"]).astype(str) == "gaussian_detect"
    )
    if len(gaussian) == 0:
        return brightest, near_neighbor_dist
    columns = xray["id", "snr", "y_angle", "z_angle", "brightest", "near_neighbor_dist"]
    obsids = np.asarray(xray["obsid"])[gaussian]
    by_obsid = gaussian[np.argsort(obsids, kind="stable")]
    boundaries = np.flatnonzero(np.diff(np.asarray(xray["obsid"])[by_obsid])) + 1
    for rows in np.split(by_obsid, boundaries):
        sources = with_recomputed_set_columns(columns[rows])
        brightest[rows] = sources["brightest"]
        near_neighbor_dist[rows] = sources["near_neighbor_dist"]
    return brightest, near_neighbor_dist


def backfill(db_file: Path, dry_run: bool = False) -> dict:
    """Drop excluded-seed gaussian_detect rows from `db_file`.

    brightest and near_neighbor_dist are then recomputed over the gaussian_detect rows
    left, in every obsid.
    """
    db_file = Path(db_file)
    xray = db.get_table("astromon_xray_src", db_file)
    masks = excluded_seed_gaussian_mask(xray)
    drop = masks["crowded"] | masks["grating_arm"] | masks["acis_streak"]
    n_dropped = int(drop.sum())

    kept = xray[~drop]
    brightest, near_neighbor_dist = recomputed_gaussian_set_columns(kept)
    brightest_changed = brightest != np.asarray(kept["brightest"])
    distance_changed = ~np.isclose(
        near_neighbor_dist,
        np.asarray(kept["near_neighbor_dist"]),
        rtol=0,
        atol=_SAME_DISTANCE_ARCSEC,
        equal_nan=True,
    )
    kept["brightest"][brightest_changed] = brightest[brightest_changed]
    kept["near_neighbor_dist"][distance_changed] = near_neighbor_dist[distance_changed]
    changed_obsids = np.unique(
        np.asarray(kept["obsid"])[brightest_changed | distance_changed]
    )

    logger.info(f"astromon_xray_src rows: {len(xray):,}")
    logger.info(
        f"gaussian_detect rows seeded from an excluded celldetect source: "
        f"{n_dropped:,} (crowded: {int(masks['crowded'].sum()):,}, "
        f"grating_arm: {int(masks['grating_arm'].sum()):,}, "
        f"acis_streak: {int(masks['acis_streak'].sum()):,})"
    )

    logger.info(
        f"recomputed over the gaussian_detect rows left: brightest changes in"
        f" {int(brightest_changed.sum()):,} rows, near_neighbor_dist in"
        f" {int(distance_changed.sum()):,} rows, {len(changed_obsids):,} obsids"
    )

    if dry_run:
        logger.info("Dry run -- nothing written.")
    elif n_dropped or len(changed_obsids):
        db.save(
            "astromon_xray_src",
            kept,
            db_file,
            ignore_obsid=True,
            expect_existing=True,
        )
        logger.info(
            "written; run drop_retired_rows next, since the removed rows may carry"
            " astromon_21 matches"
        )
    else:
        logger.info("Nothing to change.")

    return {
        "db_file": db_file,
        "dropped": n_dropped,
        "dropped_crowded": int(masks["crowded"].sum()),
        "dropped_grating_arm": int(masks["grating_arm"].sum()),
        "dropped_acis_streak": int(masks["acis_streak"].sum()),
        "brightest_changed": int(brightest_changed.sum()),
        "near_neighbor_dist_changed": int(distance_changed.sum()),
        "changed_obsids": changed_obsids.tolist(),
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
