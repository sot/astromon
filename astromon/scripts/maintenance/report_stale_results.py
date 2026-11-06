"""Report the task and cached results that the current code versions make stale.

Why this is needed
------------------
Every task (``@task(version=N)``) and every cached function
(``@stored_result(version=N)``) stores its version with each result, and a result
from another version is not reused: the task runs again, or the cached value is
recomputed, the next time that obsid is processed -- and a stale task reruns
everything downstream of it. Bumping a cheap derivation (the ``sources`` table,
rebuilt from the ``.src``) costs little; bumping make_images, celldetect or
gaussian_detect means CIAO reruns for every obsid processed afterwards. This says
how many, and which of them lack the filtered evt2 those reruns read.

What it does
------------
Reads the stored versions -- ``task_<name>.json``, and the ``<name>__*`` files in
each obsid's ``cache/`` (the work copy if there is one, else the archive copy, as
the pipeline would) -- and lists every result whose version differs from the code's,
plus the tasks downstream of a stale task. ``--bump NAME`` adds one to that task's or
cached function's version first, to see what a planned bump would cost. Read only.

Usage
-----
::

    python -m astromon.scripts.maintenance.report_stale_results \\
        --root /Volumes/Black/data/astromon/work \\
        --root /Volumes/Black/data/astromon/archive \\
        [--bump gaussian_detect] [--out stale.ecsv]
"""

import argparse
import json
import logging
from collections import Counter
from collections.abc import Iterable, Sequence
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path

from astropy.io import fits
from astropy.table import Table

from astromon.scripts.maintenance.backfill_obs_gti import obsid_work_dir
from astromon.stored_result import FORMATS, TABLE_VERSION_KEY

logger = logging.getLogger("report_stale_results")

Versions = dict[str, int]


def stored_task_version(path: Path) -> int:
    """The version stored in a ``task_<name>.json`` (0 if it predates versions)."""
    return int(json.loads(Path(path).read_text()).get("version", 0))


def stored_cache_version(path: Path) -> int:
    """The version stored in a cached result file (0 if it predates versions)."""
    path = Path(path)
    if path.name.endswith((".fits", ".fits.gz")):
        return int(fits.getheader(path, 1).get(TABLE_VERSION_KEY, 0))
    fmt = "json" if path.name.endswith(".json") else "pickle"
    return FORMATS[fmt].load_versioned(path)[1]


def bumped(
    task_versions: Versions, cache_versions: Versions, bump: Iterable[str]
) -> tuple[Versions, Versions]:
    """Copies of the two version maps with one added to each name in `bump`."""
    task_versions, cache_versions = dict(task_versions), dict(cache_versions)
    for name in bump:
        if name in task_versions:
            task_versions[name] += 1
        elif name in cache_versions:
            cache_versions[name] += 1
        else:
            raise ValueError(
                f"nothing named {name!r} among the tasks or cached results"
            )
    return task_versions, cache_versions


def _effective_cache_files(obsid: int, roots: Sequence[Path]) -> dict[str, Path]:
    """File name -> path in ``cache/``, the first root's copy winning (work over archive)."""
    files = {}
    for root in reversed(list(roots)):
        cache = obsid_work_dir(root, obsid) / "cache"
        if cache.is_dir():
            files.update(
                {f.name: f for f in cache.iterdir() if not f.name.startswith(".")}
            )
    return files


def obsid_stale_rows(
    obsid: int,
    roots: Sequence[Path],
    task_versions: Versions,
    cache_versions: Versions,
    dependents: dict[str, set[str]],
) -> list[dict]:
    """The stale results of one obsid, one row each."""
    files = _effective_cache_files(obsid, roots)
    filtered_evt2 = any(
        (
            obsid_work_dir(root, obsid) / "primary" / f"{obsid}_evt2_filtered.fits.gz"
        ).exists()
        for root in roots
    )
    rows = []

    def row(kind: str, name: str, reason: str, filename: str) -> dict:
        return {
            "obsid": obsid,
            "kind": kind,
            "name": name,
            "reason": reason,
            "file": filename,
            "filtered_evt2": filtered_evt2,
        }

    task_files = {
        name: files[f"task_{name}.json"]
        for name in task_versions
        if f"task_{name}.json" in files
    }
    own = {}
    for name, path in task_files.items():
        stored = stored_task_version(path)
        if stored != task_versions[name]:
            own[name] = f"version {stored}, code {task_versions[name]}"
    for name, path in task_files.items():
        if name in own:
            reason = own[name]
        else:
            upstream = sorted(
                stale for stale in own if name in dependents.get(stale, set())
            )
            if not upstream:
                continue
            reason = "downstream of " + ", ".join(upstream)
        rows.append(row("task", name, reason, path.name))

    for filename, path in sorted(files.items()):
        name = filename.split("__")[0]
        if "__" not in filename or name not in cache_versions:
            continue
        stored = stored_cache_version(path)
        if stored != cache_versions[name]:
            rows.append(
                row(
                    "cache",
                    name,
                    f"version {stored}, code {cache_versions[name]}",
                    filename,
                )
            )
    return rows


def _obsids(roots: Sequence[Path]) -> list[int]:
    return sorted(
        {
            int(directory.name)
            for root in roots
            for directory in Path(root).glob("obs*/*")
            if directory.is_dir() and directory.name.isdigit()
        }
    )


def stale_results(
    roots: Sequence[Path],
    task_versions: Versions,
    cache_versions: Versions,
    dependents: dict[str, set[str]],
    processes: int = 1,
) -> list[dict]:
    """Every stale result under `roots` (see :func:`obsid_stale_rows`)."""
    one = partial(
        obsid_stale_rows,
        roots=[Path(root) for root in roots],
        task_versions=task_versions,
        cache_versions=cache_versions,
        dependents=dependents,
    )
    obsids = _obsids(roots)
    if processes > 1:
        with ProcessPoolExecutor(processes) as pool:
            per_obsid = list(pool.map(one, obsids, chunksize=100))
    else:
        per_obsid = [one(obsid) for obsid in obsids]
    return [row for rows in per_obsid for row in rows]


def current_versions() -> tuple[Versions, Versions, dict[str, set[str]]]:
    """The code's task and cached-result versions, and each task's dependents."""
    from astromon.observation import Observation  # noqa: PLC0415  (registers the tasks)
    from astromon.stored_result import StoredResultWrapper  # noqa: PLC0415
    from astromon.task import TASKS  # noqa: PLC0415

    task_versions = {name: task.version for name, task in TASKS.tasks.items()}
    cache_versions = {
        wrapper.name: wrapper.version
        for wrapper in vars(Observation).values()
        if isinstance(wrapper, StoredResultWrapper)
    }
    dependents = {name: set(TASKS.get_dependents(name)) for name in TASKS.tasks}
    return task_versions, cache_versions, dependents


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--root",
        required=True,
        type=Path,
        action="append",
        dest="roots",
        help="obsNN/<obsid> tree; give the work tree first, then the archive",
    )
    parser.add_argument(
        "--bump",
        action="append",
        default=[],
        help="add one to this task's or cached result's version",
    )
    parser.add_argument("--processes", type=int, default=8)
    parser.add_argument(
        "--out", type=Path, help="write one row per stale result (ECSV)"
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")

    task_versions, cache_versions, dependents = current_versions()
    task_versions, cache_versions = bumped(task_versions, cache_versions, args.bump)
    rows = stale_results(
        args.roots, task_versions, cache_versions, dependents, args.processes
    )

    counts = Counter((row["kind"], row["name"]) for row in rows)
    for (kind, name), count in sorted(counts.items()):
        logger.info(f"{kind:5s} {name:24s} {count:7,d} stale")
    task_obsids = {row["obsid"] for row in rows if row["kind"] == "task"}
    no_input = {
        row["obsid"]
        for row in rows
        if row["kind"] == "task" and not row["filtered_evt2"]
    }
    logger.info(
        f"obsids with a task to rerun: {len(task_obsids):,}"
        f" ({len(no_input):,} without a filtered evt2 in any root)"
    )
    if args.out:
        if rows:
            Table(rows=rows).write(args.out, format="ascii.ecsv", overwrite=True)
        logger.info(f"wrote {len(rows):,} row(s) to {args.out}")


if __name__ == "__main__":
    main()
