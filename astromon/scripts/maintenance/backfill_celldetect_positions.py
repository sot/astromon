"""Re-persist celldetect positions computed with the current evt2-header attitude.

Background
----------
gaussian_detect and celldetect both derive y_angle/z_angle from the same
evt2-header quaternion (`Observation.get_evt2_info()`'s ra_pnt/dec_pnt/roll_pnt
-- see `Observation._get_sources`). Neither one has ever used a different
attitude source. What went stale is `astromon_xray_src`'s celldetect rows
themselves: a full gaussian_detect-only reprocessing campaign wrote
`(obsid, "gaussian_detect")`-keyed rows via `db.save()`, which replaces rows
matching that key only -- the existing `(obsid, "celldetect")` rows for the
same obsids were never touched, so they still reflect whatever attitude was
current when celldetect last ran for that obsid (in some cases, before an
earlier pointing correction).

For any obsid whose `sources/{obsid}_celldetect.src`, its psf_size sidecar,
and a cached evt2_info are already on disk, `Observation.get_sources(version=
"celldetect")` recomputes the correct positions from those files alone --
no CIAO call, no download, no subprocess. This script forces that
recomputation (bypassing get_sources' own stored-result cache, which does not
know the attitude fix ever happened) for every obsid in a list, and re-persists
the results to astromon_xray_src in batches, exactly like `run_all.py`'s
--batch-size does for the same reason: db.save()'s remove-then-recreate
rewrite must not run once per obsid.

Any obsid this fast path cannot handle (missing .src, missing psf_size, or
any other exception) is left for a full `run_all.py --versions celldetect`
pass instead -- reported at the end, not silently dropped.

Usage::

    python backfill_celldetect_positions.py \\
        --obsid-list obsids_gaussian_redo.txt \\
        --db-file astromon_gaussian_redo.h5 \\
        --workdir /Volumes/Black/data/astromon/work \\
        --archive-dir /Volumes/Black/data/astromon/archive \\
        --batch-size 500
"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path

from astropy.table import vstack

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
logger = logging.getLogger("backfill_celldetect_positions")


def compact_db(db_file: Path) -> tuple[int, int]:
    """Repack `db_file` in place with ptrepack, reclaiming dead node space.

    See run_all.py's compact_db -- same rationale: db.save()'s
    remove-then-recreate rewrite leaves the previous copy of whatever it just
    replaced as unreachable but still-allocated space, and only a full
    repack actually reclaims it.
    """
    before = db_file.stat().st_size
    ptrepack = Path(sys.executable).parent / "ptrepack"
    tmp_path = db_file.with_suffix(".compacting.h5")
    tmp_path.unlink(missing_ok=True)
    subprocess.run(
        [
            str(ptrepack),
            "--chunkshape=auto",
            "--propindexes",
            "--complevel=6",
            str(db_file),
            str(tmp_path),
        ],
        check=True,
        capture_output=True,
    )
    tmp_path.replace(db_file)
    after = db_file.stat().st_size
    return before, after


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--obsid-list", required=True, type=Path)
    parser.add_argument("--db-file", required=True, type=Path)
    parser.add_argument("--workdir", required=True, type=Path)
    parser.add_argument("--archive-dir", required=True, type=Path)
    parser.add_argument("--batch-size", default=500, type=int)
    args = parser.parse_args()

    from astromon import db  # noqa: PLC0415
    from astromon.observation import Observation  # noqa: PLC0415

    obsids = [int(line.strip()) for line in open(args.obsid_list) if line.strip()]
    logger.info(f"{len(obsids)} obsids to backfill")

    batch: list = []
    needs_full_processing: list[tuple[int, str]] = []
    n_done = 0
    n_empty = 0

    def flush(batch: list) -> None:
        if not batch:
            return
        combined = vstack(batch, metadata_conflicts="silent")
        db.save("astromon_xray_src", combined, args.db_file, expect_existing=True)
        before, after = compact_db(args.db_file)
        logger.info(
            f"merged {len(combined)} row(s) from {len(batch)} obsid(s); "
            f"compacted {before / 1e6:.1f} MB -> {after / 1e6:.1f} MB"
        )

    for i, obsid in enumerate(obsids, start=1):
        try:
            obs = Observation(
                obsid,
                args.workdir,
                archive_dir=args.archive_dir,
                source="archive",
            )
            obs.get_sources.invalidate_result(version="celldetect")
            sources = obs.get_sources(version="celldetect")
        except Exception as exc:
            needs_full_processing.append((obsid, f"{type(exc).__name__}: {exc}"))
            continue

        if len(sources) == 0:
            n_empty += 1
            continue

        batch.append(sources)
        n_done += 1

        if len(batch) >= args.batch_size:
            flush(batch)
            batch = []

        if i % 500 == 0:
            logger.info(
                f"[{i}/{len(obsids)}] {n_done} backfilled, {n_empty} empty, "
                f"{len(needs_full_processing)} need full processing"
            )

    flush(batch)

    logger.info(
        f"done: {n_done} backfilled, {n_empty} empty (no celldetect sources), "
        f"{len(needs_full_processing)} need full processing"
    )
    if needs_full_processing:
        fallback_path = args.obsid_list.with_name(
            args.obsid_list.stem + "_needs_full_celldetect.txt"
        )
        fallback_path.write_text(
            "\n".join(str(obsid) for obsid, _ in needs_full_processing) + "\n"
        )
        logger.info(
            f"obsids needing a full run_all.py --versions celldetect pass "
            f"written to {fallback_path}"
        )
        for obsid, reason in needs_full_processing[:20]:
            logger.info(f"  {obsid}: {reason}")


if __name__ == "__main__":
    main()
