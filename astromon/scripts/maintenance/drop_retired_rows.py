"""Drop what astromon no longer queries or computes from a database.

* astromon_cat_src rows from the retired catalogs, :data:`RETIRED_CATALOGS`;
* astromon_xcorr rows of any select_name the pipeline no longer stores (anything but
  :data:`astromon.cross_match.STORED_SELECTIONS`): astromon_22-27 and the
  single-catalog selections. The default selection is computed on demand instead.

Every other row keeps its id, so the stored astromon_21 rows keep pointing at the same
catalog sources. The run refuses, writing nothing, if a stored row would be left
pointing at a dropped catalog source.

Usage::

    python -m astromon.scripts.maintenance.drop_retired_rows --db astromon.h5 --dry-run
    python -m astromon.scripts.maintenance.drop_retired_rows --db astromon.h5
"""

import argparse
import logging
from collections import Counter
from pathlib import Path

import numpy as np

from astromon import db
from astromon.cross_match import STORED_SELECTIONS

logger = logging.getLogger("astromon")

RETIRED_CATALOGS = ("Quaia", "GaiaVarStar")


def drop_retired(dbfile: Path, dry_run: bool = False) -> dict:
    """Drop retired catalogs' candidates and retired selections' matches from `dbfile`.

    Returns ``{"cat_src_dropped": {catalog: n}, "xcorr_dropped": {select_name: n}}``.

    Raises
    ------
    RuntimeError
        If a stored match would be left pointing at a dropped catalog source.
    """
    cat_src = db.get_table("astromon_cat_src", dbfile)
    xcorr = db.get_table("astromon_xcorr", dbfile)

    catalogs = np.asarray(cat_src["catalog"]).astype(str)
    drop_cat = np.isin(catalogs, RETIRED_CATALOGS)
    select_names = np.asarray(xcorr["select_name"]).astype(str)
    keep_xcorr = np.isin(select_names, STORED_SELECTIONS)

    dropped_sources = set(
        zip(
            np.asarray(cat_src["obsid"])[drop_cat].tolist(),
            np.asarray(cat_src["id"])[drop_cat].tolist(),
            strict=True,
        )
    )
    orphaned = [
        f"obsid {obsid} {name} c_id {c_id}"
        for obsid, c_id, name in zip(
            np.asarray(xcorr["obsid"])[keep_xcorr].tolist(),
            np.asarray(xcorr["c_id"])[keep_xcorr].tolist(),
            select_names[keep_xcorr].tolist(),
            strict=True,
        )
        if (obsid, c_id) in dropped_sources
    ]
    if orphaned:
        raise RuntimeError(
            f"{len(orphaned)} stored match(es) point at a retired catalog's source, "
            f"so dropping it would leave them dangling: {orphaned[:5]}. Nothing written."
        )

    summary = {
        "cat_src_dropped": dict(sorted(Counter(catalogs[drop_cat].tolist()).items())),
        "xcorr_dropped": dict(
            sorted(Counter(select_names[~keep_xcorr].tolist()).items())
        ),
    }
    logger.info(
        f"astromon_cat_src: dropping {int(drop_cat.sum()):,} of {len(cat_src):,} rows"
        f" {summary['cat_src_dropped']}"
    )
    logger.info(
        f"astromon_xcorr: dropping {int((~keep_xcorr).sum()):,} of {len(xcorr):,} rows"
        f" {summary['xcorr_dropped']}"
    )
    if not dry_run:
        db.save(
            "astromon_cat_src",
            cat_src[~drop_cat],
            dbfile,
            ignore_obsid=True,
            expect_existing=True,
        )
        db.save(
            "astromon_xcorr",
            xcorr[keep_xcorr],
            dbfile,
            ignore_obsid=True,
            expect_existing=True,
        )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true", dest="dry_run")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    if not args.db.exists():
        parser.error(f"db does not exist: {args.db}")
    drop_retired(args.db, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
