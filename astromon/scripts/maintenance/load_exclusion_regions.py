"""Load the absolute-astrometry calibration's exclusions into astromon_regions.

The list lives in ``astromon/data/excluded_sources.ecsv`` so it is versioned with the
code. Each entry becomes an exclusion region, which
``db.get_cross_matches(exclude_regions=True)`` applies:

* a catalog source, in every observation: a circle of
  :data:`CATALOG_SOURCE_RADIUS_ARCSEC` around the source's position in
  astromon_cat_src. That is the cross-match's own dr limit, so it removes every
  match to the source, whatever the observation's target is called;
* a whole observation: a circle of :data:`WHOLE_OBSERVATION_RADIUS_ARCSEC` around
  its pointing, for that obsid only.

Loading is idempotent: an entry whose region is already in the table (recognised by
the ``exclusion <key>:`` prefix of its comments) is not added again. An entry the
database cannot place, because it holds no such catalog source or observation, is
reported and skipped rather than guessed.

Usage::

    python -m astromon.scripts.maintenance.load_exclusion_regions --db astromon.h5 --dry-run
    python -m astromon.scripts.maintenance.load_exclusion_regions --db astromon.h5
"""

import argparse
import logging
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table

from astromon import db

logger = logging.getLogger("load_exclusion_regions")

EXCLUSIONS_FILE = Path(__file__).parents[2] / "data" / "excluded_sources.ecsv"

#: The cross-match's dr limit: every match to a catalog source lies within it.
CATALOG_SOURCE_RADIUS_ARCSEC = 3.0
#: Larger than any detector's field around its pointing (HRC-S is ~99 arcmin long).
WHOLE_OBSERVATION_RADIUS_ARCSEC = 3600.0
#: astromon_regions.comments is S200.
_COMMENTS_LENGTH = 200


def read_exclusions(exclusions_file: Path | None = None) -> Table:
    """The exclusion list, with empty catalog/name for whole-observation entries.

    >>> exclusions = read_exclusions()
    >>> int(np.sum(exclusions["obsid"] > 0))
    3
    """
    exclusions = Table.read(exclusions_file or EXCLUSIONS_FILE, format="ascii.ecsv")
    for name in ("catalog", "name", "reason"):
        column = exclusions[name]
        exclusions[name] = column.filled("") if hasattr(column, "filled") else column
    return exclusions


def _key(catalog: str, name: str, obsid: int) -> str:
    return f"obsid {obsid}" if obsid > 0 else f"{catalog} {name}"


def _comments(key: str, reason: str) -> str:
    return f"{_comments_prefix(key)}{reason}"[:_COMMENTS_LENGTH]


def _comments_prefix(key: str) -> str:
    """What every region loaded for `key` starts its comments with."""
    return f"exclusion {key}: "


def _mean_position(ra: np.ndarray, dec: np.ndarray) -> tuple[float, float]:
    """Mean of positions on the sphere, so a source straddling RA=0 stays there."""
    cartesian = SkyCoord(ra * u.deg, dec * u.deg).cartesian
    mean = SkyCoord(
        x=np.mean(cartesian.x),
        y=np.mean(cartesian.y),
        z=np.mean(cartesian.z),
        representation_type="cartesian",
    ).spherical
    return float(mean.lon.to_value(u.deg)) % 360.0, float(mean.lat.to_value(u.deg))


def exclusion_regions(
    exclusions: Table, cat_src: Table, obs: Table
) -> tuple[Table, list[str]]:
    """Place each exclusion on the sky.

    Parameters
    ----------
    exclusions : astropy.table.Table
        As returned by :func:`read_exclusions`.
    cat_src : astropy.table.Table
        astromon_cat_src (``catalog``, ``name``, ``ra``, ``dec``).
    obs : astropy.table.Table
        astromon_obs (``obsid``, ``ra``, ``dec``).

    Returns
    -------
    tuple
        ``(regions, missing)``: the astromon_regions rows (``ra``, ``dec``, ``radius``,
        ``obsid``, ``comments``, plus the entry's ``key``, which astromon_regions does
        not store) and the keys of the entries that could not be placed.
    """
    catalogs = np.asarray(cat_src["catalog"]).astype(str)
    names = np.asarray(cat_src["name"]).astype(str)
    obsids = np.asarray(obs["obsid"]).astype(int)
    rows = []
    missing = []
    for catalog, name, obsid, reason in zip(
        np.asarray(exclusions["catalog"]).astype(str),
        np.asarray(exclusions["name"]).astype(str),
        np.asarray(exclusions["obsid"]).astype(int),
        np.asarray(exclusions["reason"]).astype(str),
        strict=True,
    ):
        key = _key(catalog, name, obsid)
        if obsid > 0:
            pointing = obsids == obsid
            if not pointing.any():
                missing.append(key)
                continue
            index = int(np.flatnonzero(pointing)[0])
            ra, dec = float(obs["ra"][index]), float(obs["dec"][index])
            radius, region_obsid = WHOLE_OBSERVATION_RADIUS_ARCSEC, obsid
        else:
            chosen = (catalogs == catalog) & (names == name)
            if not chosen.any():
                missing.append(key)
                continue
            ra, dec = _mean_position(
                np.asarray(cat_src["ra"], dtype=float)[chosen],
                np.asarray(cat_src["dec"], dtype=float)[chosen],
            )
            radius, region_obsid = CATALOG_SOURCE_RADIUS_ARCSEC, 0
        rows.append((ra, dec, radius, region_obsid, _comments(key, reason), key))
    regions = Table(
        rows=rows or None,
        names=["ra", "dec", "radius", "obsid", "comments", "key"],
        dtype=[float, float, float, int, str, str],
    )
    return regions, missing


def load(
    dbfile: Path,
    exclusions_file: Path | None = None,
    user: str = "astromon",
    dry_run: bool = False,
) -> dict:
    """Add the exclusions not yet in `dbfile`'s astromon_regions.

    Returns ``{"added": n, "already_present": n, "missing": [key, ...]}``.
    """
    exclusions = read_exclusions(exclusions_file)
    regions, missing = exclusion_regions(
        exclusions,
        db.get_table("astromon_cat_src", dbfile),
        db.get_table("astromon_obs", dbfile),
    )
    existing = np.asarray(db.get_table("astromon_regions", dbfile)["comments"])
    existing = existing.astype(str)
    present = np.array(
        [
            any(comment.startswith(_comments_prefix(key)) for comment in existing)
            for key in regions["key"]
        ],
        dtype=bool,
    )
    new = regions[~present]
    new["user"] = user
    new.remove_column("key")
    for key in missing:
        logger.warning(f"exclusion {key}: not in the database, so not placed")
    logger.info(
        f"exclusion regions: {len(new)} to add, {int(present.sum())} already present,"
        f" {len(missing)} not placeable"
    )
    if len(new) and not dry_run:
        db.add_regions(new, dbfile)
    return {
        "added": len(new),
        "already_present": int(present.sum()),
        "missing": missing,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--db", required=True, type=Path)
    parser.add_argument(
        "--exclusions-file", type=Path, default=None, help="default: the bundled list"
    )
    parser.add_argument("--user", default="astromon")
    parser.add_argument("--dry-run", action="store_true", dest="dry_run")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    if not args.db.exists():
        parser.error(f"db does not exist: {args.db}")
    load(args.db, args.exclusions_file, user=args.user, dry_run=args.dry_run)


if __name__ == "__main__":
    main()
