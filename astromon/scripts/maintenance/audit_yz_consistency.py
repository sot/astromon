"""Flag obsids whose X-ray and catalog y/z disagree with the sky.

What it checks
--------------
Every astromon_cat_src row names its nearest celldetect source
(``celldetect_x_id``). When both sides' y_angle/z_angle come from the same
pointing, the y/z distance between the two is their sky separation. When one
side was computed with another pointing -- e.g. a ``get_sources`` cache written
with the quick-look obspar and reused after the switch to the evt2 header -- the
y/z distance also carries the pointing difference, which is the same for the
whole obsid. So, per obsid,

    excess = median over anchored cat_src rows of (dr_yz - sky_sep)

is about 0 when the obsid is consistent and about the pointing difference (6-40"
in practice) when it is not. Obsids with ``excess > --threshold`` (1") are listed.

This is the check that found the 593 obsids repaired on 2026-11-04. It needs
only the database, so rerun it after any reprocessing or repair. It reads and
writes nothing else.

What it does not see
--------------------
- An obsid with no anchored catalog row (no candidates, or none whose
  ``celldetect_x_id`` names a stored celldetect source). On 2026-11-04 that hid
  42 of the 634 affected obsids.
- Which side is wrong: it flags a mismatch between the tables (62877's was in
  cat_src), and gaussian_detect rows are not anchors, so their angles are not
  tested.

``fix_pointing_angles.py --dry-run`` tests every row of both tables against the
evt2 pointing instead, but needs the evt2 files on disk.

Usage
-----
::

    python -m astromon.scripts.maintenance.audit_yz_consistency \\
        --db astromon.h5 --obsids-out yz_inconsistent.txt [--out yz_excess.ecsv]
"""

import argparse
import logging
from pathlib import Path

import numpy as np
from astropy.coordinates import angular_separation
from astropy.table import Table, join

from astromon import db

logger = logging.getLogger("audit_yz_consistency")

DEFAULT_THRESHOLD = 1.0  # arcsec


def anchor_pairs(xray_src: Table, cat_src: Table) -> Table:
    """Each cat_src row joined to its celldetect anchor, with dr_yz, sky_sep and excess (arcsec).

    Rows whose ``celldetect_x_id`` names no stored celldetect source (0 means none)
    are left out.
    """
    celldetect = xray_src[
        np.asarray(xray_src["detect_method"]).astype(str) == "celldetect"
    ]
    anchors = Table(
        {
            "obsid": np.asarray(celldetect["obsid"], dtype=np.int64),
            "celldetect_x_id": np.asarray(celldetect["id"], dtype=np.int64),
            "x_ra": np.asarray(celldetect["ra"], dtype=float),
            "x_dec": np.asarray(celldetect["dec"], dtype=float),
            "x_y_angle": np.asarray(celldetect["y_angle"], dtype=float),
            "x_z_angle": np.asarray(celldetect["z_angle"], dtype=float),
        }
    )
    candidates = Table(
        {
            "obsid": np.asarray(cat_src["obsid"], dtype=np.int64),
            "c_id": np.asarray(cat_src["id"], dtype=np.int64),
            "celldetect_x_id": np.asarray(cat_src["celldetect_x_id"], dtype=np.int64),
            "c_ra": np.asarray(cat_src["ra"], dtype=float),
            "c_dec": np.asarray(cat_src["dec"], dtype=float),
            "c_y_angle": np.asarray(cat_src["y_angle"], dtype=float),
            "c_z_angle": np.asarray(cat_src["z_angle"], dtype=float),
        }
    )
    pairs = join(candidates, anchors, keys=["obsid", "celldetect_x_id"])
    pairs["dr_yz"] = np.hypot(
        pairs["c_y_angle"] - pairs["x_y_angle"],
        pairs["c_z_angle"] - pairs["x_z_angle"],
    )
    pairs["sky_sep"] = (
        np.degrees(
            angular_separation(
                np.radians(pairs["c_ra"]),
                np.radians(pairs["c_dec"]),
                np.radians(pairs["x_ra"]),
                np.radians(pairs["x_dec"]),
            )
        )
        * 3600
    )
    pairs["excess"] = pairs["dr_yz"] - pairs["sky_sep"]
    return pairs


def yz_excess_by_obsid(xray_src: Table, cat_src: Table) -> Table:
    """Per obsid with at least one anchored cat_src row: n_pairs and median excess (arcsec).

    Example::

        excess = yz_excess_by_obsid(xray_src, cat_src)
        flagged = inconsistent_obsids(excess)
    """
    pairs = anchor_pairs(xray_src, cat_src)
    obsids = np.asarray(pairs["obsid"])
    order = np.argsort(obsids, kind="stable")
    unique, starts = np.unique(obsids[order], return_index=True)
    groups = np.split(np.asarray(pairs["excess"])[order], starts[1:])
    return Table(
        {
            "obsid": unique,
            "n_pairs": [len(group) for group in groups],
            "excess": [float(np.median(group)) for group in groups],
        }
    )


def inconsistent_obsids(
    excess: Table, threshold: float = DEFAULT_THRESHOLD
) -> list[int]:
    """Obsids whose median excess is above `threshold` arcsec."""
    return [int(obsid) for obsid in excess["obsid"][excess["excess"] > threshold]]


def audit(dbfile: Path) -> Table:
    """:func:`yz_excess_by_obsid` for the tables in `dbfile`."""
    return yz_excess_by_obsid(
        db.get_table("astromon_xray_src", dbfile),
        db.get_table("astromon_cat_src", dbfile),
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--db", required=True, type=Path, dest="db_file")
    parser.add_argument("--threshold", type=float, default=DEFAULT_THRESHOLD)
    parser.add_argument(
        "--obsids-out", type=Path, help="write the flagged obsids here, one per line"
    )
    parser.add_argument(
        "--out", type=Path, help="write every obsid's n_pairs and excess here (ECSV)"
    )
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
    if not args.db_file.exists():
        parser.error(f"db does not exist: {args.db_file}")

    excess = audit(args.db_file)
    flagged = inconsistent_obsids(excess, args.threshold)
    logger.info(
        f"{len(excess):,} obsid(s) with anchored catalog rows;"
        f' {len(flagged):,} with median(dr_yz - sky_sep) > {args.threshold}"'
    )
    if args.out:
        excess.write(args.out, format="ascii.ecsv", overwrite=True)
        logger.info(f"wrote {args.out}")
    if args.obsids_out:
        args.obsids_out.write_text("".join(f"{obsid}\n" for obsid in flagged))
        logger.info(f"wrote {args.obsids_out}")


if __name__ == "__main__":
    main()
