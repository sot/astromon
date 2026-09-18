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

If an obsid's astromon_cat_src is itself missing candidates a stored match used,
recomputing from it would silently drop that match rather than repair it -- this
refuses that obsid (see ``_check_rebuild_does_not_lose_matches``) instead of
writing a smaller result. The fix is to re-query that obsid's catalogs with
``requery_cat_src.py`` first, then rerun this script; ``--force`` overrides the
refusal but does not fix the missing candidates, it just accepts losing them.
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


def find_selections_with_fewer_matches(before: Table, after: Table | None) -> dict:
    """(obsid, select_name) pairs where a rebuild would store fewer xcorr rows.

    A rebuild recomputes from the stored cat_src, so if an obsid's candidates for
    some catalog have gone missing the recomputed result silently drops those
    matches. That is data loss dressed up as a repair, and it is not what this
    script is for: it means cat_src itself is incomplete and the catalogs have to
    be re-queried first.

    This used to be judged by checking which catalogs a select_name's hierarchy
    still had candidates for in cat_src, but that reasoning cannot actually tell
    "this catalog was queried and genuinely found nothing" apart from "this
    catalog's candidates existed and were later dropped" -- both look identical
    as catalog-absent-for-this-obsid, and the former is the common case (most
    obsids will not have a GaiaVarStar or DESIV161 candidate, and that is normal,
    not loss). Any input-side heuristic runs into the same wall, because the
    completeness information it would need does not exist in the stored tables.

    Comparing outcomes sidesteps the question entirely: `before` and `after` are
    the existing and recomputed xcorr rows for the obsids being rebuilt, and this
    just counts rows per (obsid, select_name) in each. A catalog that was queried
    and found nothing never contributed a row either way, so it cannot produce a
    deficit here. A catalog whose candidates existed and vanished can only show up
    as fewer recomputed rows than before -- which is exactly, and only, the harm
    this check exists to catch, regardless of which catalog or how many were
    involved.

    What this does not catch: a source whose match silently changes to a worse
    (but still valid) catalog keeps the same row count, so a quality regression
    like that passes clean. That is a separate, already-documented limitation --
    see the module docstring's note on `save_with_lock`'s xcorr drop -- not
    something this check was ever able to detect.

    Returns ``{(obsid, select_name): deficit}`` for pairs where `after` has fewer
    rows than `before` (a `select_name`/`obsid` present in `after` but not
    `before` scores 0, since that is a new selection, not a loss).
    """
    before_counts = Counter(
        zip(
            np.asarray(before["obsid"]).tolist(),
            np.asarray(before["select_name"]).astype(str).tolist(),
            strict=True,
        )
    )
    after_counts: Counter = Counter()
    if after is not None and len(after):
        after_counts = Counter(
            zip(
                np.asarray(after["obsid"]).tolist(),
                np.asarray(after["select_name"]).astype(str).tolist(),
                strict=True,
            )
        )
    return {
        key: n_before - after_counts.get(key, 0)
        for key, n_before in before_counts.items()
        if after_counts.get(key, 0) < n_before
    }


def select_names_in(xcorr: Table) -> list[str]:
    """The select_names actually present, so a rebuild does not invent new ones."""
    return sorted(set(np.asarray(xcorr["select_name"]).astype(str).tolist()))


def _check_rebuild_does_not_lose_matches(
    before: Table, after: Table | None, force: bool = False
) -> None:
    """Refuse a rebuild that would store fewer xcorr rows than exist today."""
    fewer = find_selections_with_fewer_matches(before, after)
    if not fewer:
        return
    detail = ", ".join(
        f"obsid {obsid} {select_name}: -{deficit}"
        for (obsid, select_name), deficit in sorted(fewer.items())
    )
    message = (
        f"{sum(fewer.values())} existing xcorr row(s) would be lost, not repaired,"
        f" by this rebuild ({detail})."
        " Something removed candidates from the stored astromon_cat_src that these"
        " matches used, so recomputing would silently drop them. Re-query the"
        " affected obsids with requery_cat_src.py first, then retry this rebuild --"
        " that is the actual fix. force=True is not a substitute for that: it just"
        " accepts the loss and writes the smaller result anyway."
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

    _check_rebuild_does_not_lose_matches(xcorr[target], new, force=force)

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
        help="rebuild even where recomputing would store fewer xcorr rows than"
        " exist today -- this DROPS the matches that would otherwise be lost."
        " Usually the actual fix is to re-query the affected obsids with"
        " requery_cat_src.py first and retry without --force",
    )
    args = parser.parse_args()
    if not args.db.exists():
        parser.error(f"db does not exist: {args.db}")
    obsids = [int(x) for x in args.obsids.read_text().split()] if args.obsids else None
    rebuild(
        args.db,
        obsids=obsids,
        select_names=args.select_names,
        dry_run=args.dry_run,
        force=args.force,
    )
    logger.info("Done.")


if __name__ == "__main__":
    main()
