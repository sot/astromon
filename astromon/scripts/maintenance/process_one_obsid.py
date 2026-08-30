"""
Process and save a single obsid as a subprocess.

The actual per-obsid work (detection, catalog queries, cross-matches,
archiving) is done by process_obsid() in get_cat_obs_data.py -- the
single source of truth. This script adds:

  - Subprocess isolation: run by run_all.py one obsid at a time, so a
    hard crash in CIAO only kills this process, not the whole run.
  - RESULT line: always prints "RESULT: {json}" last, giving run_all.py
    a reliable outcome signal independent of the return code.
  - File-locked DB writes: safe for --parallel runs.
  - Per-obsid logging: all output captured in logs/<obsid>.log by
    run_all.py.
  - Persisted per-obsid status: every outcome (success/skipped/
    skipped_not_public/failure) is recorded in the astromon_status table
    (see db.save_status), not just in run_all.py's --tracking-csv or the
    log files. This is what lets a query against the DB itself tell "never
    attempted" (no astromon_status row) apart from "attempted, skipped, no
    x-ray sources found" (a row with status="skipped") -- both currently
    look identical from astromon_obs alone, since a skipped obsid never
    gets one of those rows either. Where known, the archive's processing
    version (ascdsver) at the time of the attempt is recorded alongside it,
    so a later rerun can tell a genuinely-reprocessed obsid (worth
    retrying) from one that has not changed since it was last skipped.

Usage::

    python -m astromon.scripts.maintenance.process_one_obsid \\
        <obsid> <workdir> <db_file> [<archive_dir>]
        [--source arc5gl|archive] [--ciao-prefix /path/to/ciao]
        [--versions celldetect gaussian_detect ...]
"""

import fcntl
import json
import logging
import os
import signal
import threading
import time
import traceback
from pathlib import Path

import numpy as np
from ska_helpers.logging import basic_logger

from astromon import db
from astromon.observation import ObsidNotPubliclyAvailable
from astromon.scripts.get_cat_obs_data import (
    Skipped,
    SkippedWithWarning,
    process_obsid,
)


def _start_parent_watchdog() -> None:
    """Kill this process group when the orchestrating parent process exits.

    When run_all.py is killed externally (e.g. disk eject during sleep,
    SIGKILL from the OS) it cannot clean up its children before dying.
    This daemon thread detects that by watching os.getppid(): once the
    parent PID becomes 1 (launchd/init) the parent is gone, and we
    terminate ourselves and every CIAO subprocess we've spawned.

    Called once at the top of main(); runs as a daemon thread so it does
    not prevent normal exit.
    """
    expected_ppid = os.getppid()

    def _watch():
        while True:
            time.sleep(5)
            if os.getppid() != expected_ppid:
                # Parent is gone — kill our entire process group.
                for sig in (signal.SIGTERM, signal.SIGKILL):
                    try:
                        os.killpg(os.getpgrp(), sig)
                    except Exception:
                        pass
                    time.sleep(3)

    threading.Thread(target=_watch, daemon=True).start()


def emit_result(**fields):
    print("RESULT: " + json.dumps(fields), flush=True)


def _drop_xcorr_for_obsid(con, obsid: int) -> int:
    """Remove all astromon_xcorr rows for `obsid`, whatever their detect_method.

    Returns the number of rows removed.
    """
    xcorr = db.get_table("astromon_xcorr", con)
    keep = xcorr["obsid"] != obsid
    n_dropped = int(len(xcorr) - keep.sum())
    if n_dropped:
        db.save(
            "astromon_xcorr", xcorr[keep], con, ignore_obsid=True, expect_existing=True
        )
    return n_dropped


def _drop_cat_src_for_obsid(con, obsid: int) -> int:
    """Remove all astromon_cat_src rows for `obsid`. Returns how many went.

    Only called when the caller has said its (empty) table is authoritative for
    the obsid. db.save cannot express this on its own: an empty table is a no-op
    there, and one with the obsid's rows removed leaves them in place, which is
    documented behaviour rather than an oversight.
    """
    cat_src = db.get_table("astromon_cat_src", con)
    keep = np.asarray(cat_src["obsid"]) != obsid
    n_removed = int((~keep).sum())
    if n_removed:
        db.save(
            "astromon_cat_src",
            cat_src[keep],
            con,
            ignore_obsid=True,
            expect_existing=True,
        )
    return n_removed


def save_with_lock(
    db_file: Path,
    tables: dict,
    obsid: int | None = None,
    replace_cat_src: bool = False,
    status: str | None = None,
    note: str = "",
    ascdsver: str = "",
):
    """Serialize DB writes across concurrently-running worker processes.

    HDF5/PyTables doesn't support safe concurrent multi-process writes.
    The slow detection/catalog work above this point is unaffected --
    only the brief final write is serialized.

    Also fills any columns required by the schema but absent from the
    data (columns produced by only some detect methods) with NaN for
    floats or 0 for integers.

    Saves use expect_existing=True: db.save replaces a table non-atomically
    (remove_node then create_table), so a worker SIGKILLed in that window --
    run_all does this on its per-obsid timeout -- leaves the table missing,
    and an unguarded save would recreate it holding only this obsid's rows.

    `replace_cat_src` says the caller's astromon_cat_src is the authoritative
    content for this obsid, so an empty table means "this obsid has none" and the
    obsid's stored cat_src rows are removed along with its xcorr rows. Without it
    an empty table is skipped and the stored rows are left alone, which is the
    safe default: a rerun against a smaller catalog set legitimately produces
    fewer -- possibly zero -- candidates, and inferring deletion from a row count
    would silently discard real rows.

    Both tables have to go together in that case. Dropping only xcorr would leave
    catalog rows with no matches, which is indistinguishable from a genuine
    no-match result and so a worse state than leaving both.

    When `obsid` is given and astromon_cat_src is among the tables being written,
    every astromon_xcorr row for that obsid is dropped first. astromon_cat_src has
    no detect_method column, so db.save keys it on obsid alone and rewriting it
    renumbers every c_id for the obsid; an xcorr row from a detect_method this run
    did not redo would keep its old c_id and then join to the wrong catalog source
    in db.get_cross_matches. The consequence is deliberate: a rerun with a subset
    of --versions leaves xcorr only for the versions it actually ran.

    `status`/`note`/`ascdsver`, when given, additionally record this obsid's
    outcome in astromon_status under the same lock -- see db.save_status. Passing
    `tables={}` with only `status` set (no obsid-keyed tables to write) is
    expected for the skip/failure paths, which have nothing else to save.
    """
    lock_path = Path(str(db_file) + ".lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with open(lock_path, "w") as lockfile:
        fcntl.flock(lockfile, fcntl.LOCK_EX)
        try:
            # Initialize a brand-new database with all tables present, so that from
            # here on a missing table means corruption rather than a first write.
            # Done under the lock so concurrent workers cannot race on creation, and
            # only when the file itself is absent -- a table missing from a file that
            # already exists is exactly the damage expect_existing=True must catch.
            if not Path(db_file).exists():
                db.create_empty_tables(db_file)

            with db.connect(db_file, mode="r+") as con:
                cat_src = tables.get("astromon_cat_src")
                rewriting_cat_src = obsid is not None and cat_src is not None
                if rewriting_cat_src and (len(cat_src) or replace_cat_src):
                    n_dropped = _drop_xcorr_for_obsid(con, obsid)
                    if n_dropped:
                        logging.getLogger("astromon").debug(
                            f"OBSID={obsid} dropped {n_dropped} stale astromon_xcorr "
                            "row(s) before rewriting astromon_cat_src"
                        )
                if rewriting_cat_src and replace_cat_src and not len(cat_src):
                    n_removed = _drop_cat_src_for_obsid(con, obsid)
                    if n_removed:
                        logging.getLogger("astromon").info(
                            f"OBSID={obsid} removed {n_removed} astromon_cat_src "
                            "row(s): this run found no candidates for it"
                        )

                for name, data in tables.items():
                    if data is None or len(data) == 0:
                        continue
                    if name in db.DTYPES:
                        db.conform_to_dtype(data, name)
                    db.save(name, data, con, expect_existing=True)

                if status is not None and obsid is not None:
                    db.save_status(con, obsid, status, note=note, ascdsver=ascdsver)
        finally:
            fcntl.flock(lockfile, fcntl.LOCK_UN)


def main():
    import argparse  # noqa: PLC0415

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("obsid", type=int)
    parser.add_argument("workdir", type=Path)
    parser.add_argument("db_file", type=Path)
    # archive_dir is positional but optional so that the run_all.py call
    # signature (which passes it as the fourth positional) stays compatible.
    parser.add_argument("archive_dir", nargs="?", default=None, type=Path)
    parser.add_argument(
        "--source",
        default="arc5gl",
        choices=["arc5gl", "archive"],
        help="Data download method (default: arc5gl)",
    )
    parser.add_argument(
        "--ciao-prefix",
        default=None,
        dest="ciao_prefix",
        help="Path to CIAO installation (default: /soft/ciao)",
    )
    parser.add_argument(
        "--versions",
        nargs="+",
        default=["celldetect"],
        help="Detection versions to run (default: celldetect)",
    )
    parser.add_argument(
        "--cleanup",
        action="store_true",
        default=False,
        help="Delete downloaded files not needed after source detection",
    )
    parser.add_argument(
        "--skip-catalog-match",
        action="store_true",
        default=False,
        dest="skip_catalog_match",
        help="Detect sources but skip the catalog rough_match/cross-match step "
        "entirely (astromon_cat_src/astromon_xcorr come back empty). For a "
        "detection-only bulk pass whose catalog candidates get filled in "
        "afterward, in batches, by requery_cat_src.py.",
    )
    args = parser.parse_args()

    _start_parent_watchdog()

    obsid = args.obsid
    workdir = str(args.workdir)
    db_file = args.db_file
    # Guard against "" → Path(".") via argparse type=Path: treat the CWD sentinel as None.
    _raw_archive_dir = args.archive_dir
    archive_dir = (
        str(_raw_archive_dir)
        if _raw_archive_dir and _raw_archive_dir != Path(".")
        else None
    )
    versions = tuple(args.versions)

    basic_logger(
        "astromon",
        level="DEBUG",
        format="%(asctime)s %(funcName)-25s: %(message)s",
    )
    logger = logging.getLogger("astromon")

    try:
        result = process_obsid(
            obsid,
            workdir,
            archive_dir=archive_dir,
            versions=versions,
            source=args.source,
            ciao_prefix=args.ciao_prefix,
            cleanup=args.cleanup,
            skip_catalog_match=args.skip_catalog_match,
        )
        save_with_lock(
            db_file,
            {
                "astromon_obs": result["astromon_obs"],
                "astromon_xray_src": result["astromon_xray_src"],
                "astromon_cat_src": result["astromon_cat_src"],
                "astromon_xcorr": result["astromon_xcorr"],
            },
            obsid=obsid,
            # The pipeline recomputes this obsid's catalog candidates from
            # scratch, so its result is authoritative -- and already treated that
            # way whenever it is non-empty, since db.save replaces every row for
            # the obsids present in the data. Passing this makes the empty case
            # behave the same instead of being the one exception that silently
            # keeps the previous run's rows.
            replace_cat_src=True,
            status="success",
            note=f"{len(result['astromon_xray_src'])} sources, "
            f"{len(result['astromon_xcorr'])} xcorr"
            + (" (catalog match skipped)" if args.skip_catalog_match else ""),
            ascdsver=result["astromon_obs"]["ascdsver"][0]
            if len(result["astromon_obs"])
            else "",
        )
        emit_result(
            obsid=obsid,
            status="success",
            n_sources=len(result["astromon_xray_src"]),
            n_xcorr=len(result["astromon_xcorr"]),
        )
    except ObsidNotPubliclyAvailable as e:
        logger.info(f"OBSID={obsid} not publicly available, skipping: {e}")
        save_with_lock(
            db_file, {}, obsid=obsid, status="skipped_not_public", note=str(e)
        )
        emit_result(obsid=obsid, status="skipped_not_public", note=str(e))
    except (Skipped, SkippedWithWarning) as e:
        logger.info(f"OBSID={obsid} skipped: {e}")
        save_with_lock(
            db_file, {}, obsid=obsid, status="skipped", note=str(e), ascdsver=e.ascdsver
        )
        emit_result(obsid=obsid, status="skipped", note=str(e))
    except Exception as e:
        logger.error(f"OBSID={obsid} failed: {e}")
        logger.error(traceback.format_exc())
        save_with_lock(db_file, {}, obsid=obsid, status="failure", note=str(e))
        emit_result(obsid=obsid, status="failure", stage="exception", error=str(e))


if __name__ == "__main__":
    main()
