"""Re-query catalog candidates for obsids that are already in the database.

``astromon_cat_src`` is missing catalog candidates for an unknown set of obsids.
The gap cannot be measured from inside the database: a catalog that was never
queried for a field looks exactly like a catalog that was queried and had nothing
there, and a wholesale table rewrite renumbers ``id`` contiguously so it leaves no
gap behind either. The only way to find out what is missing is to ask the catalogs
again and compare.

Doing that does not require reprocessing the observations. Every input the catalog
queries need is already stored or available locally:

* the x-ray source positions, from ``astromon_xray_src`` (``celldetect``, the
  method the pipeline itself queries from);
* the pointing quaternion and observation date, from the local mica obspar
  archive.

So this runs without CIAO, without event files and without the archive -- which
also means it does not touch ``astromon_xray_src`` and cannot disturb the
detection results or the ``peak_offset`` diagnostic.

Two properties keep the write safe:

* Rows are written keyed on ``(obsid, catalog)``, so re-querying RFC cannot
  disturb the obsid's 2MASS rows. The earlier one-off backfills rewrote the whole
  table from a merged copy, which made every run decide the fate of every catalog.
* New ids are allocated above the obsid's highest existing id, so adding a catalog
  the obsid never had leaves every ``astromon_xcorr.c_id`` still pointing at the
  same row.

Replacing a catalog the obsid *did* have does renumber those rows, so the obsids
where that happened are reported: their ``astromon_xcorr`` rows have to be
recomputed with ``rebuild_xcorr`` afterwards.

Usage::

    # cheap first pass: the catalogs served from a local cache, no network
    python -m astromon.scripts.maintenance.requery_cat_src --db astromon.h5 --all

    # everything, including the VizieR and Gaia TAP queries
    python -m astromon.scripts.maintenance.requery_cat_src --db astromon.h5 --all \
        --catalogs all --workers 4

    python -m astromon.scripts.maintenance.requery_cat_src --db astromon.h5 \
        --obsids 8541 8571 --dry-run
"""

import argparse
import json
import logging
from collections.abc import Callable, Sequence
from multiprocessing.pool import ThreadPool
from pathlib import Path

import numpy as np
from astropy import coordinates as coords
from astropy import units as u
from astropy.table import Table, vstack
from cxotime import CxoTime
from Quaternion import Quat
from ska_helpers.logging import basic_logger

from astromon import db
from astromon.cross_match import (
    assign_cat_src_ids,
    get_desi_v161_candidates,
    get_gaia_agn,
    get_gaia_qso_candidates,
    get_gaia_var_stars,
    get_milliquas_gaia,
    get_quaia_candidates,
    rough_match,
)

logger = basic_logger("astromon", level="INFO")

SEARCH_RADIUS = 3 * u.arcsec

# Served from a local cache, so a full pass over every obsid costs no network at all.
LOCAL_CATALOGS = ("RFC", "ICRF3", "GaiaAGN", "GaiaQSO", "Quaia")

# One VizieR cone search or Gaia TAP round trip per obsid each. Orders of magnitude
# slower than the local ones, which is why --catalogs defaults to the local set.
REMOTE_CATALOGS = (
    "Tycho2",
    "USNO-B1.0",
    "2MASS",
    "SDSS",
    "MilliquasGaia",
    "DESIV161",
    "GaiaVarStar",
)

ALL_CATALOGS = LOCAL_CATALOGS + REMOTE_CATALOGS


class CatalogQueryError(Exception):
    """A catalog query failed, so this obsid's result would be incomplete."""


def _rough_match_getter(catalog: str) -> Callable:
    """Wrap :func:`rough_match` so one catalog can be queried on its own.

    rough_match emits one row per (source, candidate) pair within the radius,
    carrying the per-pair anchor as ``x_id`` until assign_cat_src_ids renames it
    to ``celldetect_x_id``.
    """

    # logging_tag is part of the shared getter signature but unused here: rough_match
    # builds its own tag from the obsid column of the table it is handed.
    def getter(
        sources: Table,
        obs_time: CxoTime,
        logging_tag: str = "",  # noqa: ARG001
    ) -> Table:
        # Pass a column subset: rough_match adds a working column to what it is
        # given, and the caller reuses `sources` for every catalog.
        return Table(
            rough_match(
                sources[["obsid", "id", "ra", "dec"]], obs_time, catalogs=(catalog,)
            )
        )

    return getter


def _position_getter(query: Callable) -> Callable:
    """Wrap a ``get_*(pos, radius, logging_tag)`` catalog function."""

    def getter(sources: Table, obs_time: CxoTime, logging_tag: str = "") -> Table:
        # obstime matters: get_gaia_var_stars corrects for proper motion with it.
        pos = coords.SkyCoord(
            sources["ra"], sources["dec"], unit="deg", obstime=obs_time
        )
        return Table(query(pos, radius=SEARCH_RADIUS, logging_tag=logging_tag))

    return getter


CATALOG_GETTERS: dict[str, Callable] = {
    "RFC": _rough_match_getter("RFC"),
    "ICRF3": _rough_match_getter("ICRF3"),
    "Tycho2": _rough_match_getter("Tycho2"),
    "USNO-B1.0": _rough_match_getter("USNO-B1.0"),
    "2MASS": _rough_match_getter("2MASS"),
    "SDSS": _rough_match_getter("SDSS"),
    "GaiaAGN": _position_getter(get_gaia_agn),
    "GaiaQSO": _position_getter(get_gaia_qso_candidates),
    "Quaia": _position_getter(get_quaia_candidates),
    "MilliquasGaia": _position_getter(get_milliquas_gaia),
    "DESIV161": _position_getter(get_desi_v161_candidates),
    "GaiaVarStar": _position_getter(get_gaia_var_stars),
}


def load_completed(progress_file: Path | None) -> dict[int, set[str]]:
    """What a previous run finished, as ``{obsid: {catalog, ...}}``.

    Records are ``"<obsid> <catalog>"``, one per line, appended as the run goes.
    Keying on the pair rather than the obsid alone is what makes resuming across
    different ``--catalogs`` sets safe: an obsid that got ten catalogs in one pass
    is not complete for a run that also wants the eleventh.

    Two kinds of line are ignored. A torn final record, because a run killed
    mid-append leaves one, and only newline-terminated records count -- "70" is
    both a plausible obsid and the first half of "7001". And a bare obsid with no
    catalog, which is the older format: it says nothing about which catalogs ran,
    and assuming it covered everything is precisely the mistake this keying
    removes, so those obsids are re-queried.
    """
    if progress_file is None or not Path(progress_file).exists():
        return {}
    records = Path(progress_file).read_text().split("\n")
    if records and records[-1] != "":
        logger.debug(f"ignoring unterminated final record {records[-1]!r}")
    completed: dict[int, set[str]] = {}
    legacy = 0
    for record in records[:-1]:
        fields = record.split()
        if not fields:
            continue
        if len(fields) < 2:
            legacy += 1
            continue
        try:
            obsid = int(fields[0])
        except ValueError:
            logger.debug(f"ignoring unparsable progress record {record!r}")
            continue
        completed.setdefault(obsid, set()).add(fields[1])
    if legacy:
        logger.warning(
            f"{progress_file}: {legacy} record(s) name an obsid but no catalog."
            " They predate per-catalog progress and cannot say what ran, so those"
            " obsids will be re-queried."
        )
    return completed


def record_completed(progress_file: Path, done: dict) -> None:
    """Append the catalogs finished for each obsid, as ``{obsid: [catalog, ...]}``.

    Appended rather than rewritten, so a killed run loses at most its last record,
    and accumulating: a later pass adds catalogs to an obsid instead of replacing
    what it already had.
    """
    with open(progress_file, "a") as handle:
        handle.write(
            "".join(
                f"{int(obsid)} {catalog}\n"
                for obsid, catalogs in done.items()
                for catalog in catalogs
            )
        )
        handle.flush()


def is_complete(completed: dict, obsid: int, catalogs: Sequence[str]) -> bool:
    """Whether `obsid` has already been re-queried for every one of `catalogs`."""
    return set(catalogs) <= completed.get(int(obsid), set())


def empty_cat_src() -> Table:
    """A zero-row astromon_cat_src table with every column the schema requires."""
    return Table(np.zeros(0, dtype=db.ASTROMON_CAT_SRC_DTYPE))


def celldetect_sources(xray: Table, obsid: int) -> Table:
    """The celldetect rows of `xray` for `obsid`.

    celldetect is what the pipeline queries catalogs from: it is the one method
    that always runs, and the positions barely move between methods, so the
    candidate set it produces is the one every method is matched against.
    """
    selected = xray[
        (np.asarray(xray["obsid"]) == obsid)
        & (np.asarray(xray["detect_method"]).astype(str) == "celldetect")
    ]
    if len(selected) == 0:
        raise ValueError(f"OBSID={obsid} no celldetect sources in astromon_xray_src")
    return selected


def next_cat_src_id(cat: Table, obsid: int) -> int:
    """The first id that is free for `obsid`, i.e. above every id it already has.

    ``astromon_cat_src.id`` only has to be unique within an obsid, and existing
    ``astromon_xcorr.c_id`` values have to keep pointing at the rows they point at
    now, so new rows go above the current maximum rather than filling any gaps.
    """
    existing = np.asarray(cat["id"])[np.asarray(cat["obsid"]) == obsid]
    return int(existing.max()) + 1 if len(existing) else 0


def obspar_pointing(obsid: int) -> tuple[Quat, CxoTime] | None:
    """The pointing quaternion and observation start from the local mica archive.

    Returns None when mica has no entry, which is the caller's cue to skip the
    obsid rather than fall back on ``astromon_obs``: the ra/dec/roll stored there
    are the nominal target values, not the pointing, so the y/z angles derived
    from them would not match the ones the pipeline writes.
    """
    from mica.archive import obspar as mica_obspar  # noqa: PLC0415

    try:
        obspar = mica_obspar.get_obspar(int(obsid))
    except Exception as exc:
        logger.debug(f"OBSID={obsid} mica obspar lookup failed: {exc}")
        return None
    if not obspar:
        return None
    quat = Quat(equatorial=(obspar["ra_pnt"], obspar["dec_pnt"], obspar["roll_pnt"]))
    return quat, CxoTime(obspar["date_obs"])


def requery_obsid(  # noqa: PLR0917
    obsid: int,
    sources: Table,
    quat: Quat,
    obs_time: CxoTime,
    catalogs: Sequence[str],
    id_offset: int = 0,
    getters: dict[str, Callable] | None = None,
) -> Table:
    """Query each of `catalogs` around `sources` and return one cat_src table.

    Ids run consecutively from `id_offset` across all the catalogs queried, which
    is how the pipeline numbers them within an obsid.

    Raises :class:`CatalogQueryError` if any catalog fails: a partial result would
    be written as though those catalogs had legitimately returned nothing, which
    is the exact failure mode this script exists to undo.
    """
    getters = CATALOG_GETTERS if getters is None else getters
    unknown = [name for name in catalogs if name not in getters]
    if unknown:
        raise KeyError(f"no catalog getter for: {', '.join(unknown)}")

    found = []
    next_id = id_offset
    for catalog in catalogs:
        try:
            candidates = getters[catalog](
                sources, obs_time, logging_tag=f"OBSID={obsid}"
            )
        except Exception as exc:
            raise CatalogQueryError(
                f"OBSID={obsid} {catalog} query failed: {exc}"
            ) from exc
        if len(candidates) == 0:
            logger.debug(f"OBSID={obsid} {catalog}: no candidates")
            continue
        assign_cat_src_ids(candidates, obsid, sources, quat, id_offset=next_id)
        next_id += len(candidates)
        found.append(db.conform_to_dtype(candidates, "astromon_cat_src"))
        logger.debug(f"OBSID={obsid} {catalog}: {len(candidates)} candidate(s)")

    if not found:
        return empty_cat_src()
    columns = list(db.ASTROMON_CAT_SRC_DTYPE.names)
    return vstack(
        [found_table[columns] for found_table in found], metadata_conflicts="silent"
    )


def write_candidates(dbfile: Path, candidates: Table, dry_run: bool = False) -> dict:
    """Write `candidates` into astromon_cat_src, one (obsid, catalog) key at a time.

    Returns ``{"added_rows": n, "replaced_obsids": [...], "replaced_pairs": n}``.
    ``replaced_obsids`` are the obsids where a catalog that already had rows was
    rewritten: those rows got new ids, so the obsid's ``astromon_xcorr`` entries
    are stale and need ``rebuild_xcorr``. Obsids that only gained a catalog they
    never had are not listed -- their new ids sit above anything an existing
    ``c_id`` refers to.

    Note the one thing keying on ``(obsid, catalog)`` cannot do: a catalog that
    returns nothing writes no rows, so it also removes none. This adds and
    refreshes candidates; it never concludes that stored rows should not exist.
    """
    if len(candidates) == 0:
        return {"added_rows": 0, "replaced_obsids": [], "replaced_pairs": 0}

    stored = db.get_table("astromon_cat_src", dbfile) if Path(dbfile).exists() else None
    replaced_obsids: set[int] = set()
    replaced_pairs = 0
    if stored is not None and len(stored):
        stored_pairs = set(
            zip(
                np.asarray(stored["obsid"]).tolist(),
                np.asarray(stored["catalog"]).astype(str).tolist(),
                strict=True,
            )
        )
        incoming_pairs = set(
            zip(
                np.asarray(candidates["obsid"]).tolist(),
                np.asarray(candidates["catalog"]).astype(str).tolist(),
                strict=True,
            )
        )
        for obsid, catalog in incoming_pairs:
            if (obsid, catalog) in stored_pairs:
                replaced_obsids.add(int(obsid))
                replaced_pairs += 1

    result = {
        "added_rows": len(candidates),
        "replaced_obsids": sorted(replaced_obsids),
        "replaced_pairs": replaced_pairs,
    }
    if dry_run:
        logger.info(
            f"dry run: would write {len(candidates):,} row(s), "
            f"replacing {replaced_pairs} existing (obsid, catalog) pair(s)"
        )
        return result

    db.save(
        "astromon_cat_src",
        candidates,
        dbfile,
        replace_keys=("obsid", "catalog"),
        expect_existing=True,
    )
    return result


def _query_one(args: tuple) -> tuple[int, Table | None, str | None]:
    """Query one obsid; return (obsid, candidates, error). Runs in a worker thread."""
    obsid, sources, id_offset, catalogs = args
    pointing = obspar_pointing(obsid)
    if pointing is None:
        return obsid, None, "no mica obspar entry"
    quat, obs_time = pointing
    try:
        candidates = requery_obsid(
            obsid, sources, quat, obs_time, catalogs, id_offset=id_offset
        )
    except CatalogQueryError as exc:
        return obsid, None, str(exc)
    return obsid, candidates, None


def requery(  # noqa: PLR0917
    dbfile: Path,
    obsids: Sequence[int],
    catalogs: Sequence[str],
    dry_run: bool = False,
    workers: int = 1,
    chunk_size: int = 200,
    progress_file: Path | None = None,
    resume: bool = False,
) -> dict:
    """Re-query `catalogs` for `obsids`, writing a chunk at a time.

    Writing per chunk rather than once at the end bounds both memory and how much
    work a crash throws away. Each chunk re-reads the stored ids so a catalog
    replaced in an earlier chunk does not hand out ids twice.

    With `resume`, obsids already listed in `progress_file` are skipped -- keyed
    by (obsid, catalog), so a run with a narrower catalog set cannot mark an
    obsid complete for catalogs it never queried. Finished records are appended
    only after each chunk is safely written.
    """
    xray = db.get_table("astromon_xray_src", dbfile)
    requested = list(obsids)
    already_done = load_completed(progress_file) if resume else {}
    obsids = [
        obsid for obsid in requested if not is_complete(already_done, obsid, catalogs)
    ]
    resumed = len(requested) - len(obsids)
    if resumed:
        logger.info(f"resuming: {resumed:,} obsid(s) already recorded as done")
    logger.info(
        f"re-querying {len(catalogs)} catalog(s) for {len(obsids):,} obsid(s): "
        f"{', '.join(catalogs)}"
    )

    summary = {
        "resumed": resumed,
        "obsids": len(obsids),
        "catalogs": list(catalogs),
        "added_rows": 0,
        "replaced_obsids": [],
        "queried": 0,
        "skipped": {},
    }

    for start in range(0, len(obsids), chunk_size):
        chunk = list(obsids[start : start + chunk_size])
        cat = db.get_table("astromon_cat_src", dbfile)

        work = []
        for obsid in chunk:
            try:
                sources = celldetect_sources(xray, obsid)
            except ValueError as exc:
                summary["skipped"][str(obsid)] = str(exc)
                continue
            work.append((obsid, sources, next_cat_src_id(cat, obsid), catalogs))

        if workers > 1:
            with ThreadPool(workers) as pool:
                results = pool.map(_query_one, work)
        else:
            results = [_query_one(item) for item in work]

        found, completed = [], []
        for obsid, candidates, error in results:
            if error is not None:
                logger.warning(f"OBSID={obsid} skipped: {error}")
                summary["skipped"][str(obsid)] = error
                continue
            summary["queried"] += 1
            completed.append(obsid)
            if len(candidates):
                found.append(candidates)

        if found:
            written = write_candidates(
                dbfile, vstack(found, metadata_conflicts="silent"), dry_run=dry_run
            )
            summary["added_rows"] += written["added_rows"]
            summary["replaced_obsids"].extend(written["replaced_obsids"])

        # Only after the rows are safely written, so a resume never skips an obsid
        # whose candidates never landed.
        if completed and progress_file is not None and not dry_run:
            record_completed(progress_file, dict.fromkeys(completed, catalogs))
        logger.info(
            f"{min(start + chunk_size, len(obsids)):,}/{len(obsids):,} obsids -- "
            f"{summary['added_rows']:,} rows written, "
            f"{len(summary['skipped'])} skipped"
        )

    summary["replaced_obsids"] = sorted(set(summary["replaced_obsids"]))
    return summary


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--db", required=True, type=Path, help="astromon h5 file")
    obsid_group = parser.add_mutually_exclusive_group(required=True)
    obsid_group.add_argument(
        "--all", action="store_true", help="every obsid in astromon_obs"
    )
    obsid_group.add_argument("--obsids", nargs="+", type=int, help="obsids to re-query")
    obsid_group.add_argument(
        "--obsid-file", type=Path, help="file with one obsid per line, or a JSON list"
    )
    parser.add_argument(
        "--catalogs",
        nargs="+",
        default=list(LOCAL_CATALOGS),
        help=(
            "catalogs to re-query, or 'local' / 'remote' / 'all'. "
            f"Default: {' '.join(LOCAL_CATALOGS)}"
        ),
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="threads for the catalog queries; only worth raising for remote catalogs",
    )
    parser.add_argument("--chunk-size", type=int, default=200)
    parser.add_argument(
        "--progress-file",
        type=Path,
        help=(
            "append each finished (obsid, catalog) pair here, so --resume can skip"
            " what is already done -- per catalog, so a run with a different"
            " --catalogs set cannot mark an obsid complete for catalogs it skipped"
        ),
    )
    parser.add_argument(
        "--resume",
        action="store_true",
        help=(
            "skip obsids whose every requested catalog is already recorded in"
            " --progress-file"
        ),
    )
    parser.add_argument("--dry-run", action="store_true", dest="dry_run")
    parser.add_argument(
        "--summary-file", type=Path, help="write the JSON summary here as well"
    )
    parser.add_argument("--log-level", default="INFO")
    return parser


def _resolve_catalogs(names: Sequence[str]) -> list[str]:
    groups = {"local": LOCAL_CATALOGS, "remote": REMOTE_CATALOGS, "all": ALL_CATALOGS}
    if len(names) == 1 and names[0] in groups:
        return list(groups[names[0]])
    unknown = [name for name in names if name not in CATALOG_GETTERS]
    if unknown:
        raise SystemExit(
            f"unknown catalog(s): {', '.join(unknown)}. "
            f"Known: {', '.join(ALL_CATALOGS)}"
        )
    return list(names)


def _resolve_obsids(args: argparse.Namespace) -> list[int]:
    if args.all:
        return sorted(
            int(obsid) for obsid in db.get_table("astromon_obs", args.db)["obsid"]
        )
    if args.obsid_file:
        content = args.obsid_file.read_text().strip()
        if content.startswith("["):
            return [int(obsid) for obsid in json.loads(content)]
        return [int(line) for line in content.split() if line]
    return list(args.obsids)


def main() -> None:
    args = get_parser().parse_args()
    logging.getLogger("astromon").setLevel(args.log_level.upper())

    catalogs = _resolve_catalogs(args.catalogs)
    obsids = _resolve_obsids(args)

    summary = requery(
        args.db,
        obsids,
        catalogs,
        dry_run=args.dry_run,
        workers=args.workers,
        chunk_size=args.chunk_size,
    )

    logger.info(
        f"done: {summary['queried']:,} obsid(s) queried, "
        f"{summary['added_rows']:,} row(s) written, "
        f"{len(summary['skipped'])} skipped"
    )
    if summary["replaced_obsids"]:
        logger.info(
            f"{len(summary['replaced_obsids'])} obsid(s) had a catalog replaced and so "
            "need `rebuild_xcorr --obsids ...` to repair their c_id references"
        )
    if args.summary_file:
        args.summary_file.write_text(json.dumps(summary, indent=2))
        logger.info(f"summary written to {args.summary_file}")


if __name__ == "__main__":
    main()
