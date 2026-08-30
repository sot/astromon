"""Process obsids one per subprocess and record per-obsid outcomes.

Runs process_one_obsid.py per obsid (not the astromon-cross-match CLI --
we confirmed that swallows per-obsid failures silently: no traceback, no
non-zero exit code, no visible log output, and astromon.h5 left
untouched on a real failure). process_one_obsid.py does the same actual
work (Observation.process() -> get_sources() -> rough_match() ->
compute_cross_matches() -> db.save()) but with its own exception
handling, and always prints one "RESULT: {json}" line at the end -- that
line, not the subprocess return code, is what this script reads to
determine the real outcome. If that line is missing (e.g. a hard crash
in a CIAO subprocess call that takes the whole worker down before it can
print anything), the obsid is recorded as a failure with that noted.

Subprocess-per-obsid is kept deliberately (rather than importing and
calling process_one_obsid in-process): a hard crash in one obsid's CIAO
calls only kills that one subprocess, not this orchestrator or the rest
of the run.

Usage:
    python -m astromon.scripts.maintenance.run_all \\
        --db-file /proj/sot/ska3/test-astromon/data/astromon/tmp/astromon.h5 \\
        --workdir /proj/sot/ska3/test-astromon/data/astromon/work \\
        --archive-dir /proj/sot/ska3/test-astromon/data/astromon/archive \\
        --obsid-list obsids.txt \\
        --log-dir logs \\
        --tracking-csv tracking.csv \\
        --parallel 1 \\
        --versions celldetect gaussian_detect peak_gaussian_detect

Resumable: obsids already recorded as "success" in --tracking-csv are
skipped on rerun. To force reprocessing everything, use a fresh
--tracking-csv path or delete the old one.
"""

import argparse
import csv
import json
import os
import shutil
import signal
import subprocess
import sys
import threading
import time
import traceback
from datetime import datetime
from pathlib import Path

WORKER_MODULE = "astromon.scripts.maintenance.process_one_obsid"

# Default per-obsid wall-clock timeout.  acis_streak_map has been observed to
# hang indefinitely on certain HRC observations; this caps the damage at one
# missed obsid rather than stalling the whole run for hours.
_WORKER_TIMEOUT_SEC = 1800  # 30 minutes


def _kill_process_group(pgid: int) -> None:
    """Send SIGTERM then SIGKILL to an entire process group.

    Used to tear down a worker and all its CIAO subprocesses when the
    per-obsid timeout expires.  With start_new_session=True the worker PID
    equals its PGID, so os.killpg(proc.pid, ...) reaches every descendant.

    Both ProcessLookupError (ESRCH: the group is already gone) and
    PermissionError (EPERM) are treated as "nothing left to kill" -- macOS
    has been observed to raise EPERM rather than ESRCH for a pgid whose
    leader already exited (e.g. a concurrent external kill of the same
    worker), and there is nothing more this function can do about a
    process group it cannot signal either way. Letting either propagate
    would crash the whole orchestrator (see process_and_record's own
    try/except for why that must never happen for a single obsid).
    """
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(pgid, sig)
        except (ProcessLookupError, PermissionError):
            return  # already gone (or unsignalable) — nothing left to kill
        time.sleep(3)


def load_done_obsids(tracking_csv: Path) -> set[int]:
    done = set()
    if tracking_csv.exists():
        with open(tracking_csv) as f:
            for row in csv.DictReader(f):
                if row["status"] in (
                    "success",
                    "skipped",
                    "skipped_no_sources",
                    "skipped_not_public",
                ):
                    done.add(int(row["obsid"]))
    return done


def parse_result_line(stdout: str) -> dict | None:
    for line in reversed(stdout.splitlines()):
        if line.startswith("RESULT: "):
            return json.loads(line[len("RESULT: ") :])
    return None


def run_one(  # noqa: PLR0917
    obsid: int,
    db_file: Path,
    workdir: Path,
    log_dir: Path,
    versions: tuple = ("celldetect",),
    archive_dir: Path | None = None,
    source: str = "arc5gl",
    ciao_prefix: str | None = None,
    cleanup: bool = False,
    skip_catalog_match: bool = False,
    mirror: str | None = None,
    worker_timeout: int = _WORKER_TIMEOUT_SEC,
) -> dict:
    cmd = [
        sys.executable,
        "-m",
        WORKER_MODULE,
        str(obsid),
        str(workdir),
        str(db_file),
        # Only pass archive_dir when actually set; omitting it lets process_one_obsid.py
        # default to None rather than converting "" → Path(".") → truthy.
        *([str(archive_dir)] if archive_dir else []),
        "--source",
        source,
        "--versions",
        *versions,
    ]
    if ciao_prefix:
        cmd += ["--ciao-prefix", str(ciao_prefix)]
    if cleanup:
        cmd += ["--cleanup"]
    if skip_catalog_match:
        cmd += ["--skip-catalog-match"]
    if mirror:
        cmd += ["--mirror", mirror]

    start = time.time()
    stdout = ""
    stderr = ""
    returncode = -1
    timed_out = False

    # start_new_session=True puts the worker in its own process group so
    # that on timeout we can kill it *and* every CIAO grandchild it spawned
    # (acis_streak_map, celldetect, etc.) with a single os.killpg() call.
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        start_new_session=True,
    )
    try:
        raw_stdout, raw_stderr = proc.communicate(timeout=worker_timeout)
        stdout = raw_stdout.decode(errors="replace")
        stderr = raw_stderr.decode(errors="replace")
        returncode = proc.returncode
    except subprocess.TimeoutExpired:
        timed_out = True
        # Kill the entire process group; proc.pid == pgid with start_new_session.
        _kill_process_group(proc.pid)
        raw_stdout, raw_stderr = proc.communicate()
        stdout = raw_stdout.decode(errors="replace")
        stderr = raw_stderr.decode(errors="replace")

    elapsed = time.time() - start

    log_path = log_dir / f"{obsid}.log"
    with open(log_path, "w") as f:
        f.write(f"$ {' '.join(cmd)}\n\n")
        f.write(f"--- stdout ---\n{stdout}\n")
        f.write(f"--- stderr ---\n{stderr}\n")
        if timed_out:
            f.write(f"--- TIMED OUT after {worker_timeout}s ---\n")
        else:
            f.write(f"--- returncode: {returncode} ---\n")

    if timed_out:
        status = "failure"
        note = (
            f"worker timed out after {worker_timeout}s (likely stuck CIAO subprocess)"
        )
    else:
        result = parse_result_line(stdout)
        if result is None:
            status = "failure"
            note = "no RESULT line -- worker process likely crashed hard (segfault/OOM/killed)"
        else:
            status = result["status"]
            note = result.get("error", "")

    return {
        "obsid": obsid,
        "status": status,
        "note": note,
        "returncode": returncode,
        "elapsed_sec": round(elapsed, 1),
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "log_file": str(log_path),
    }


def main():  # noqa: PLR0915
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db-file", required=True, type=Path)
    parser.add_argument("--workdir", required=True, type=Path)
    parser.add_argument(
        "--obsid-list", required=True, type=Path, help="text file, one obsid per line"
    )
    parser.add_argument(
        "--log-dir",
        default=Path("logs"),
        type=Path,
        help="per-obsid full stdout/stderr logs, always written",
    )
    parser.add_argument(
        "--tracking-csv",
        default=Path("tracking.csv"),
        type=Path,
        help="structured per-obsid outcome log, appended incrementally",
    )
    parser.add_argument(
        "--parallel",
        default=1,
        type=int,
        help="number of obsids to run concurrently as separate worker "
        "subprocesses; database writes are serialized across them "
        "via a file lock regardless of this setting",
    )
    parser.add_argument(
        "--archive-dir",
        default=None,
        type=Path,
        help="optional pre-populated Chandra archive directory; "
        "if given, obsid data is symlinked from here rather "
        "than downloaded",
    )
    parser.add_argument(
        "--source",
        default="arc5gl",
        choices=["arc5gl", "archive"],
        help="data download method: 'arc5gl' (default, CXC-internal) "
        "or 'archive' (public CDA via download_chandra_obsid)",
    )
    parser.add_argument(
        "--ciao-prefix",
        default=None,
        dest="ciao_prefix",
        type=Path,
        help="path to CIAO installation (default: /soft/ciao); "
        "required when running outside the CXC network",
    )
    parser.add_argument(
        "--versions",
        nargs="+",
        default=["celldetect"],
        help="detection versions to run per obsid "
        "(default: celldetect). Pass multiple to run "
        "gaussian_detect and/or peak_gaussian_detect as well.",
    )
    parser.add_argument(
        "--cleanup",
        action="store_true",
        default=False,
        help="after processing each obsid, delete downloaded files "
        "that are not needed after source detection (level-1 "
        "events, VV report, preview JPEGs, etc.)",
    )
    parser.add_argument(
        "--skip-catalog-match",
        action="store_true",
        default=False,
        dest="skip_catalog_match",
        help="detect sources but skip the catalog rough_match/cross-match step "
        "entirely (astromon_cat_src/astromon_xcorr come back empty). For a "
        "detection-only bulk pass whose catalog candidates get filled in "
        "afterward, in batches, by requery_cat_src.py.",
    )
    parser.add_argument(
        "--mirror",
        default=None,
        help="CDA mirror site to pass to download_chandra_obsid --mirror (only "
        "used with --source archive). download_chandra_obsid does not fall "
        "back to the primary CDA site if the mirror does not have the obsid "
        '-- it is skipped, same as a genuine "not found". Default: use the '
        "primary CDA site directly (no mirror).",
    )
    parser.add_argument(
        "--timeout",
        default=_WORKER_TIMEOUT_SEC,
        type=int,
        dest="worker_timeout",
        help=f"per-obsid wall-clock timeout in seconds (default: {_WORKER_TIMEOUT_SEC}). "
        "Obsids exceeding this are killed and recorded as failure.",
    )
    parser.add_argument(
        "--preserve-workdir",
        default=None,
        type=Path,
        dest="preserve_workdir",
        help="after each obsid finishes, move its workdir subtree "
        "here for later reprocessing without re-downloading. "
        "Uses the same obs{NN}/{obsid} layout as the working "
        "directory. Useful when workdir is on fast local storage "
        "and preserve-workdir is on larger external storage.",
    )
    args = parser.parse_args()

    args.log_dir.mkdir(parents=True, exist_ok=True)
    args.workdir.mkdir(parents=True, exist_ok=True)

    with open(args.obsid_list) as f:
        obsids = [int(line.strip()) for line in f if line.strip()]

    done = load_done_obsids(args.tracking_csv)
    todo = [o for o in obsids if o not in done]
    print(
        f"{len(obsids)} obsids total, {len(done)} already recorded done, "
        f"{len(todo)} to process"
    )

    write_lock = threading.Lock()
    csv_is_new = not args.tracking_csv.exists()
    csv_file = open(args.tracking_csv, "a", newline="")
    csv_writer = csv.DictWriter(
        csv_file,
        fieldnames=[
            "obsid",
            "status",
            "note",
            "returncode",
            "elapsed_sec",
            "timestamp",
            "log_file",
        ],
    )
    if csv_is_new:
        csv_writer.writeheader()
        csv_file.flush()

    counts: dict[str, int] = {}

    versions = tuple(args.versions)

    def process_and_record(obsid):
        start = time.time()
        try:
            result = run_one(
                obsid,
                args.db_file,
                args.workdir,
                args.log_dir,
                versions,
                archive_dir=args.archive_dir,
                source=args.source,
                ciao_prefix=args.ciao_prefix,
                cleanup=args.cleanup,
                skip_catalog_match=args.skip_catalog_match,
                mirror=args.mirror,
                worker_timeout=args.worker_timeout,
            )

            if args.preserve_workdir is not None:
                # Mirror the Observation workdir layout: {base}/obs{obsid//1000:02d}/{obsid}
                slot = f"obs{obsid // 1000:02d}"
                src = args.workdir / slot / str(obsid)
                dst = args.preserve_workdir / slot / str(obsid)
                if src.exists():
                    dst.parent.mkdir(parents=True, exist_ok=True)
                    if dst.exists():
                        # An older copy is already there (e.g. from a prior rsync or an
                        # earlier run). The tree we just produced is the newer one, so
                        # it wins -- keeping the stale copy would silently discard this
                        # reprocessing and let a later rerun reuse the old cache/images.
                        print(
                            f"  obsid {obsid}: replacing stale {dst} with the "
                            "workdir from this run"
                        )
                        shutil.rmtree(str(dst))
                    shutil.move(str(src), str(dst))
        except Exception as exc:
            # This must never escape: a single obsid's own cleanup/bookkeeping
            # going wrong (e.g. os.killpg raising something other than
            # ProcessLookupError -- see _kill_process_group) is not the same
            # kind of failure as that obsid's actual processing, but letting
            # ANY exception propagate out of a ThreadPoolExecutor worker
            # re-raises in the main thread on the next pool.map() iteration
            # and kills the whole orchestrator -- every other in-flight and
            # not-yet-started obsid in this run along with it. The module
            # docstring's "a hard crash in one obsid's CIAO calls only kills
            # that one subprocess, not this orchestrator or the rest of the
            # run" is exactly the invariant this restores.
            print(
                f"  obsid {obsid}: process_and_record raised {exc!r}, recording "
                "as failure and continuing"
            )
            traceback.print_exc()
            result = {
                "obsid": obsid,
                "status": "failure",
                "note": f"orchestrator-side error: {exc}",
                "returncode": -1,
                "elapsed_sec": round(time.time() - start, 1),
                "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
                "log_file": "",
            }

        with write_lock:
            csv_writer.writerow(result)
            csv_file.flush()
            counts[result["status"]] = counts.get(result["status"], 0) + 1
            n_done = sum(counts.values())
            tally = ", ".join(f"{v} {k}" for k, v in counts.items())
            print(
                f"{datetime.now().strftime('%Y-%m-%dT%H:%M:%S')} "
                f"[{n_done}/{len(todo)}] obsid {obsid}: {result['status']} "
                f"({result['elapsed_sec']}s)  -- running total: {tally}"
            )

        return result

    if args.parallel <= 1:
        for obsid in todo:
            process_and_record(obsid)
    else:
        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(max_workers=args.parallel) as pool:
            list(pool.map(process_and_record, todo))

    csv_file.close()

    print()
    print("done:", ", ".join(f"{v} {k}" for k, v in counts.items()))
    if counts.get("failure"):
        print(
            f"see {args.tracking_csv} for which obsids failed, "
            f"and {args.log_dir}/<obsid>.log for each one's full output (a real "
            f"traceback, not empty stdout/stderr)"
        )


if __name__ == "__main__":
    main()
