"""Re-persist gaussian_detect y_angle/z_angle computed with the current evt2 pointing.

Background
----------
``Observation._get_sources`` used to reuse whatever y_angle/z_angle columns
were already present in a detect method's ``.src`` file, only computing them
from RA/DEC when missing. celldetect's ``.src`` never has those columns, so
it always got a fresh computation. gaussian_detect's ``.src`` always has them
(written at fit time from the attitude current when the fit ran), so if an
obsid was later reprocessed with a newer evt2/attitude but its existing
gaussian ``.src`` was reused rather than refit, the stored y_angle/z_angle
went stale while RA/DEC (read off the fit's own WCS) stayed correct. See the
fix in ``Observation._get_sources``, which now always recomputes
y_angle/z_angle/r_angle from RA/DEC with the current evt2_info.

This script forces that recomputation to be re-persisted for a given list of
obsids, bypassing ``get_sources``' own stored-result cache (which holds the
now-stale values from before the fix). It is the gaussian_detect counterpart
of ``backfill_celldetect_positions.py`` and reuses the same batching approach
for the same reason: ``db.save()``'s remove-then-recreate rewrite must not
run once per obsid.

Any obsid this fast path cannot handle (missing .src, no evt2 file on disk
and a download is needed, or any other exception) is left for separate
attention -- reported at the end, not silently dropped.

Before calling get_sources, each obsid is checked with
``Observation.sources_would_recompute()``: if the archived .src (or any other
input the version's task needs) is missing, reading it back would silently
fall through to a real detection rerun instead of a cheap re-persist -- see
that method's docstring. Obsids that would recompute are routed to
needs_attention rather than run, so this script never triggers a real refit.

Usage::

    python backfill_gaussian_yz.py \\
        --obsid-list gaussian_yz_inconsistent_obsids.txt \\
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
logger = logging.getLogger("backfill_gaussian_yz")


def compact_db(db_file: Path) -> tuple[int, int]:
    """Repack `db_file` in place with ptrepack, reclaiming dead node space."""
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
    parser.add_argument("--version", default="gaussian_detect")
    parser.add_argument("--batch-size", default=500, type=int)
    args = parser.parse_args()

    from astromon import db  # noqa: PLC0415
    from astromon.observation import Observation  # noqa: PLC0415

    obsids = [
        int(line.strip())
        for line in open(args.obsid_list)
        if line.strip() and not line.startswith("#")
    ]
    logger.info(f"{len(obsids)} obsids to backfill (version={args.version})")

    batch: list = []
    needs_attention: list[tuple[int, str]] = []
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
            if obs.sources_would_recompute(version=args.version):
                needs_attention.append(
                    (obsid, "would trigger a real detection rerun, not a re-persist")
                )
                continue
            obs.get_sources.invalidate_result(version=args.version)
            sources = obs.get_sources(version=args.version)
        except Exception as exc:
            needs_attention.append((obsid, f"{type(exc).__name__}: {exc}"))
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
                f"{len(needs_attention)} need attention"
            )

    flush(batch)

    logger.info(
        f"done: {n_done} backfilled, {n_empty} empty (no {args.version} sources), "
        f"{len(needs_attention)} need attention"
    )
    if needs_attention:
        fallback_path = args.obsid_list.with_name(
            args.obsid_list.stem + "_needs_attention.txt"
        )
        fallback_path.write_text(
            "\n".join(str(obsid) for obsid, _ in needs_attention) + "\n"
        )
        logger.info(f"obsids needing attention written to {fallback_path}")
        for obsid, reason in needs_attention[:20]:
            logger.info(f"  {obsid}: {reason}")


if __name__ == "__main__":
    main()
