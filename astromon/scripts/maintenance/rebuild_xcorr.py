"""Recompute astromon_xcorr for given obsids from data already in the database.

Why this is needed
------------------
astromon_cat_src has no detect_method column, so db.save keys it on obsid alone
and rewriting it renumbers every c_id for that obsid. An xcorr row written by an
earlier run then still holds the old c_id, which either resolves to nothing or --
worse -- resolves to a different catalog source.
``find_obsids_with_dangling_c_id`` catches the first case.

The second is not detectable from these two tables. It would need cat_src to say
which X-ray source each row was matched to, and it does not: celldetect_x_id is a
celldetect-scoped anchor, deliberately not a per-method match key, so a c_id that
resolves to the wrong row still looks entirely well-formed. A renumbered c_id
pointing at a real but different catalog source is exactly what
save_with_lock's xcorr drop exists to prevent rather than detect.

``process_one_obsid.save_with_lock`` now drops an obsid's xcorr whenever it
rewrites that obsid's cat_src, so new occurrences should not arise. This repairs
the ones already in a database.

What it does
------------
For the target obsids, deletes their xcorr rows and recomputes every select_name
from the astromon_obs, astromon_xray_src and astromon_cat_src rows already
stored. No detection is re-run, nothing is downloaded, and CIAO is not needed --
the inputs are all in the file.

Usage::

    python -m astromon.scripts.maintenance.rebuild_xcorr \\
        --db /Volumes/Black/data/astromon/astromon.h5 --dry-run

By default it repairs exactly the obsids that fail the consistency check; pass
``--obsids`` to target a specific list instead (for example after re-querying
catalogs for fields near RA=0).
"""

import argparse
import logging
from collections import Counter
from collections.abc import Sequence
from pathlib import Path

import numpy as np
from astropy.table import Table, vstack

from astromon import db
from astromon.cross_match import (
    CROSS_MATCHES_ARGS,
    compute_cross_matches,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger("rebuild_xcorr")

_XCORR_COLS = list(db.ASTROMON_XCORR_DTYPE.names)


def find_obsids_with_dangling_c_id(xcorr: Table, cat: Table) -> dict:
    """Obsids with xcorr rows whose c_id names no cat_src row.

    Returns ``{obsid: n_dangling}``.

    This used also to flag rows where ``xcorr.x_id`` disagreed with the cat_src
    row's ``x_id``, on the theory that a resolvable-but-wrong c_id could be caught
    that way. It cannot: cat_src's ``celldetect_x_id`` is a celldetect-scoped
    anchor, not the X-ray source the row was matched to, so the two differ for
    perfectly good rows -- any gaussian match, and any row whose nearest celldetect
    source is not the one the dr cut paired it with. Comparing them reported
    healthy data as broken.
    """
    known = {(int(o), int(i)) for o, i in zip(cat["obsid"], cat["id"], strict=True)}
    per_obsid: dict[int, int] = {}
    for obsid, c_id in zip(
        np.asarray(xcorr["obsid"]), np.asarray(xcorr["c_id"]), strict=True
    ):
        if (int(obsid), int(c_id)) not in known:
            per_obsid[int(obsid)] = per_obsid.get(int(obsid), 0) + 1
    return per_obsid


def find_orphaned_selections(xcorr: Table, cat: Table) -> dict:
    """xcorr rows whose select_name needs a catalog the obsid no longer has.

    A rebuild recomputes from the stored cat_src, so if an obsid's candidates for
    some catalog have gone missing the recomputed result silently drops those
    matches. That is data loss dressed up as a repair, and it is not what this
    script is for: it means cat_src itself is incomplete and the catalogs have to
    be re-queried first.

    Returns ``{select_name: count}``.
    """
    catalogs_by_obsid: dict[int, set[str]] = {}
    obsid_col = np.asarray(cat["obsid"])
    catalog_col = np.asarray(cat["catalog"]).astype(str)
    for obsid in set(obsid_col.tolist()):
        catalogs_by_obsid[int(obsid)] = set(catalog_col[obsid_col == obsid].tolist())

    orphaned: Counter = Counter()
    for obsid, select_name in zip(
        np.asarray(xcorr["obsid"]),
        np.asarray(xcorr["select_name"]).astype(str),
        strict=True,
    ):
        args = CROSS_MATCHES_ARGS.get(select_name)
        if args is None:
            continue
        present = catalogs_by_obsid.get(int(obsid), set())
        if not (set(args["catalogs"]) & present):
            orphaned[select_name] += 1
    return dict(orphaned)


def select_names_in(xcorr: Table) -> list[str]:
    """The select_names actually present, so a rebuild does not invent new ones."""
    return sorted(set(np.asarray(xcorr["select_name"]).astype(str).tolist()))


def _check_cat_src_is_complete(xcorr: Table, cat: Table, force: bool = False) -> None:
    """Refuse to rebuild from a cat_src that is missing catalogs the matches used."""
    orphaned = find_orphaned_selections(xcorr, cat)
    if not orphaned:
        return
    detail = ", ".join(f"{k}={v}" for k, v in sorted(orphaned.items()))
    message = (
        f"{sum(orphaned.values())} existing xcorr row(s) reference a select_name"
        f" whose catalogs are absent from the stored cat_src ({detail})."
        " Recomputing would drop those matches rather than repair them:"
        " astromon_cat_src is itself incomplete for these obsids, so the catalogs"
        " need re-querying before a rebuild means anything."
    )
    if not force:
        raise RuntimeError(message + " Pass force=True to rebuild anyway.")
    logger.warning(message + " Proceeding because force was requested.")


def recompute(
    select_names: Sequence[str], obs: Table, xray: Table, cat: Table
) -> Table | None:
    """Recompute matches for `select_names` from the stored tables.

    One pass per select_name over every detect method at once. That used to need a
    pass per method, against `cat` with its ``x_id`` re-pointed at that method's
    source ids, because the join keyed on ``(obsid, x_id)`` and a single stored
    anchor could only ever join to one method -- which halved every selection
    without a ``detect_method_filter``.

    Neither is needed now. simple_cross_match joins on obsid and lets the dr cut
    pair each catalogue source, and it groups per detect_method, so one pass yields
    one row per (catalogue source, X-ray source, method) exactly as the pipeline
    writes them.
    """
    methods = sorted(set(np.asarray(xray["detect_method"]).astype(str)))
    logger.info(f"detect methods present: {', '.join(methods)}")

    rebuilt = []
    for name in select_names:
        if name not in CROSS_MATCHES_ARGS:
            logger.warning(f"  {name}: unknown select_name, skipping")
            continue
        result = compute_cross_matches(
            name, astromon_obs=obs, astromon_xray_src=xray, astromon_cat_src=cat
        )
        if len(result):
            rebuilt.append(result[_XCORR_COLS])
        logger.info(f"  {name:16s} {len(result):>6,} matches")

    return vstack(rebuilt, metadata_conflicts="silent") if rebuilt else None


def rebuild(
    dbfile: Path,
    obsids: list[int] | None = None,
    select_names: list[str] | None = None,
    dry_run: bool = False,
    force: bool = False,
) -> dict:
    """Recompute xcorr for `obsids` (default: the inconsistent ones)."""
    xcorr = db.get_table("astromon_xcorr", dbfile)
    cat = db.get_table("astromon_cat_src", dbfile)
    xray = db.get_table("astromon_xray_src", dbfile)
    obs = db.get_table("astromon_obs", dbfile)

    if obsids is None:
        bad = find_obsids_with_dangling_c_id(xcorr, cat)
        obsids = sorted(bad)
        logger.info(
            f"{len(obsids)} obsid(s) with dangling references: "
            f"{sum(bad.values())} xcorr row(s) whose c_id names no cat_src row"
        )
    else:
        logger.info(f"{len(obsids)} obsid(s) requested explicitly")
    if not obsids:
        logger.info("Nothing to do.")
        return {"obsids": [], "removed": 0, "rebuilt": 0}

    if select_names is None:
        select_names = select_names_in(xcorr)
    logger.info(f"rebuilding select_names: {', '.join(select_names)}")

    target = np.isin(np.asarray(xcorr["obsid"]), obsids)
    logger.info(f"existing xcorr rows for these obsids: {int(target.sum()):,}")

    _check_cat_src_is_complete(
        xcorr[target], cat[np.isin(np.asarray(cat["obsid"]), obsids)], force=force
    )

    # Restrict the inputs once, then run each selection over the subset.
    obs_sub = obs[np.isin(np.asarray(obs["obsid"]), obsids)]
    xray_sub = xray[np.isin(np.asarray(xray["obsid"]), obsids)]
    cat_sub = cat[np.isin(np.asarray(cat["obsid"]), obsids)]
    if len(obs_sub) == 0 or len(xray_sub) == 0 or len(cat_sub) == 0:
        raise RuntimeError(
            "no obs/xray/cat rows for the requested obsids -- nothing to rebuild from"
        )

    new = recompute(select_names, obs_sub, xray_sub, cat_sub)
    n_new = len(new) if new is not None else 0
    logger.info(f"rebuilt {n_new:,} rows, replacing {int(target.sum()):,}")

    kept = xcorr[~target]
    combined = vstack([kept, new], metadata_conflicts="silent") if n_new else kept

    # Verify before writing, not after. Checking the file once it was already
    # saved meant a failure raised with the bad rows in place, leaving the
    # database worse than when the rebuild started.
    remaining = find_obsids_with_dangling_c_id(combined, cat)
    still = [o for o in obsids if o in remaining]
    if still:
        raise RuntimeError(
            f"{len(still)} obsid(s) would still have dangling c_id after rebuild: "
            f"{still[:8]}. Nothing was written."
        )

    if dry_run:
        logger.info("Dry run -- nothing written.")
    else:
        db.save(
            "astromon_xcorr", combined, dbfile, ignore_obsid=True, expect_existing=True
        )
        logger.info(f"astromon_xcorr: {len(xcorr):,} -> {len(combined):,} rows")
        logger.info("verified: no inconsistent references remain for these obsids")

    return {"obsids": obsids, "removed": int(target.sum()), "rebuilt": n_new}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", required=True, type=Path)
    parser.add_argument(
        "--obsids",
        type=Path,
        default=None,
        help="text file of obsids; default is the inconsistent ones",
    )
    parser.add_argument(
        "--select-names",
        nargs="+",
        default=None,
        help="select_names to rebuild; default is those present in the DB",
    )
    parser.add_argument("--dry-run", action="store_true", dest="dry_run")
    parser.add_argument(
        "--force",
        action="store_true",
        help="rebuild even where the stored cat_src is missing catalogs the"
        " existing matches used -- this DROPS those matches",
    )
    args = parser.parse_args()
    if not args.db.exists():
        parser.error(f"db does not exist: {args.db}")
    obsids = [int(x) for x in args.obsids.read_text().split()] if args.obsids else None
    rebuild(
        args.db, obsids=obsids, select_names=args.select_names, dry_run=args.dry_run
    )
    logger.info("Done.")


if __name__ == "__main__":
    main()
